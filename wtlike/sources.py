# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/03_sources.ipynb (unless otherwise specified).

__all__ = ['get_wtzip_index', 'WeightMan', 'weight_radius_plots', 'findsource', 'SourceLookup', 'PointSource']

# Cell
import os, sys, pickle, healpy, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

from scipy.integrate import quad
from astropy.coordinates import SkyCoord, Angle
from .config import *

# Cell
def get_wtzip_index(config, update=False):

    wtzipfile = config.datapath/'weight_files.zip'
    assert wtzipfile.is_file(), 'Expect a zip file'

    with  zipfile.ZipFile(wtzipfile) as wtzip:
        if 'index.pkl' in wtzip.namelist() and not update:
            return pickle.load(wtzip.open('index.pkl'))

        if config.verbose>0:
            print(f'Extracting info from {wtzipfile}')
        name=[]; glat=[]; glon=[]
        for filename in wtzip.namelist():
            if filename=='index.pkl': continue
            with wtzip.open(filename) as file:
                wtd = pickle.load(file, encoding='latin1')
                l,b = wtd['source_lb']
                name.append(Path(filename).name.split('_weights.pkl')[0].replace('_',' ').replace('p','+') )
                glon.append(l)
                glat.append(b)
        zip_index = dict(name=name,
                    coord=SkyCoord(glon, glat, unit='deg', frame='galactic').fk5
               )
        ### write to temp file, insert back into the zip
        ### SHould be a way to just stream
        pickle.dump(zip_index, open('/tmp/wtfile_index.pkl', 'wb'))
        with zipfile.ZipFile(wtzipfile, mode='a') as wtzip:
            wtzip.write('/tmp/wtfile_index.pkl', 'index.pkl')
    return zip_index

# Cell
class WeightMan(dict):
    """ Weight Management

    * Load weight tables
    * Assign weights to photons
    """

    def __init__(self, config, source):
        """
        """
        self.source = source
        nickname = source.nickname
        datapath =Path(config.datapath)

        filename = 'weight_files/'+nickname.replace(' ','_').replace('+','p')+'_weights.pkl'

        if (datapath/filename).exists():
#             print('found in directory')
            with open(datapath/filename, 'rb') as inp:
                wtd =  pickle.load(inp, encoding='latin1')

        elif (datapath/'weight_files.zip').exists():
            # check the zip file
#             print('load from zip')
            with  zipfile.ZipFile(datapath/'weight_files.zip') as wtzip:
                wtd =  pickle.load(wtzip.open(filename), encoding='latin1')

        else:
            raise Exception(f'No weight info found for {nickname}')

        self.update(wtd)
        self.__dict__.update(wtd)
        self.filename=filename
        self.config = config
#         pos = self['source_lb']
#         print(f'\tSource is {self["source_name"]} at ({pos[0]:.2f}, {pos[1]:.2f})')

        # check format--old has pixels, weights at tome
        srcfile = f'file "{self.filename}"' if self.source is None else f'file from source "{source.filename}"_weights.pkl'

        if hasattr(self, 'nside'):
            self.format=0
            if config.verbose>0:
                print(f'WeightMan: {srcfile} old format, nside={self.nside}')

            test_elements = 'energy_bins pixels weights nside model_name radius order roi_name'.split()
            assert np.all([x in wtd.keys() for x in test_elements]),f'Dict missing one of the keys {test_elements}'
            if config.verbose>0:
                print(f'Load weights from file {os.path.realpath(filename)}')
                pos = self['source_lb']
                print(f'\tFound: {self["source_name"]} at ({pos[0]:.2f}, {pos[1]:.2f})')
            # extract pixel ids and nside used
            self.wt_pix   = self['pixels']
            self.nside_wt = self['nside']

            # merge the weights into a table, with default nans
            # indexing is band id rows by weight pixel columns
            # append one empty column for photons not in a weight pixel
            # calculated weights are in a dict with band id keys
            self.wts = np.full((32, len(self.wt_pix)+1), np.nan, dtype=np.float32)
            weight_dict = self['weights']
            for k in weight_dict.keys():
                t = weight_dict[k]
                if len(t.shape)==2:
                    t = t.T[0] #???
                self.wts[k,:-1] = t

        else:
            self.format=1
            wtdict = self.wt_dict
            nsides = [v['nside'] for v in wtdict.values() ];

            if config.verbose>1:
                print(f'WeightMan: {srcfile} : {len(nsides)} bamds'\
                      f' with nsides {nsides[0]} to {nsides[-1]}')
            if self.source is not None:
                self.source.fit_info = self.fitinfo
                if config.verbose>2:
                    print(f'\tAdded fit info {self.fitinfo} to source')

    def _new_format(self, photons):

        wt_tables =self.wt_dict
        data_nside=1024
        photons.loc[:,'weight'] = np.nan

        if self.config.verbose>1:
            print(f'WeightMan: processing {len(photons):,} photons')

        def load_data( band_id):
            """ fetch pixels and weights for the band;
                adjust pixels to the band nside
                generate mask for pixels, weights
            """
            band = photons[photons.band==band_id] #.query('band== @band_id')
            wt_table = wt_tables[band_id]
            nside =  wt_table['nside']
            new_weights = wt_table['wts'].astype(np.float16)
            to_shift = int(2*np.log2(data_nside//nside))
            data_pixels = np.right_shift(band.nest_index, to_shift)
            wt_pixels=wt_table['pixels']
            good = np.isin( data_pixels, wt_pixels)
            if self.config.verbose>2:
                print(f'\t {band_id:2}: {len(band):8,} -> {sum(good ):8,}')
            return data_pixels, new_weights, good

        def set_weights(band_id):
            if band_id not in wt_tables.keys(): return

            data_pixels, new_weights, good = load_data(band_id)
            wt_pixels = wt_tables[band_id]['pixels']
            indicies = np.searchsorted( wt_pixels, data_pixels[good])
            new_wts = new_weights[indicies]
            # get subset of photons in this band, with new weights
            these_photons = photons[photons.band==band_id][good]
            these_photons.loc[:,'weight']=new_wts
            photons.loc[photons.band==band_id,'weight'] = (these_photons.weight).astype(np.float16)
    #         if self.config.verbose>1:
    #             print(f' -> {len(new_wts):8,}')

        for band_id in range(16):
            set_weights(band_id)

        return photons

    def add_weights(self, photons):
        """
        get the photon pixel ids, convert to NEST (if not already) and right shift them
        add column 'weight', remove `nest_index'
        remove rows with nan weight

        """
        assert photons is not None
        photons = self._new_format(photons)
        assert photons is not None

        # don't need these columns now (add flag to config to control??)
        photons.drop(['nest_index'], axis=1, inplace=True)
        noweight = np.isnan(photons.weight.values)
        if self.config.verbose>1:
            print(f'\tremove {sum(noweight):,} events without weight')

        ret = photons[np.logical_not(noweight)]
        assert ret is not None
        return ret

def weight_radius_plots(photons):
    """
    """
    import matplotlib.pyplot as plt

    fig, axx = plt.subplots(2,8, figsize=(16,5), sharex=True, sharey=True)
    plt.subplots_adjust(hspace=0.02, wspace=0)
    for id,ax in enumerate(axx.flatten()):
        subset = photons.query('band==@id & weight>0')
        ax.semilogy(subset.radius, subset.weight, '.', label=f'{id}');
        ax.legend(loc='upper right', fontsize=10)
        ax.grid(alpha=0.5)
    ax.set(ylim=(8e-4, 1.2), xlim=(0,4.9))
    plt.suptitle('Weights vs. radius per band')

# Cell
def findsource(*pars, gal=False):
    """
    Return a SkyCoord, looking up a name, or interpreting args

    Optional inputs: len(pars) is 1 for a source name or Jxxx, or 2 for coordinate pair
    - name -- look up the name, return None if not found
    - Jxxx -- intrepret interpret to get ra, dec, then convert
    - ra,dec -- assume frame=fk5
    - l,b, gal=True -- assume degrees, frame=galactic
    """

    import astropy.units as u
    if len(pars)==1:
        name = pars[0]
        if name.startswith('J'):
            # parse the name for (ra,dec)
            ra=(name[1:3]+'h'+name[3:7]+'m')
            dec = (name[7:10]+'d'+name[10:12]+'m')
            (ra,dec) = map(lambda a: float(Angle(a, unit=u.deg).to_string(decimal=True)),(ra,dec))
            skycoord = SkyCoord(ra, dec, unit='deg', frame='fk5')
        else:
            try:
                skycoord = SkyCoord.from_name(name)
            except Exception as e:
                # not found
                return None
    elif len(pars)==2:
        name = f'({pars[0]},{pars[1]})'
        #gal = kwargs.get('gal', False)
        skycoord=SkyCoord(*pars, unit='deg', frame='galactic' if gal else 'fk5')
    else:
        raise TypeError('require name or ra,dec or l,b,gal=True')
    return skycoord

# Cell
class SourceLookup():
    """ Use lists of the pointlike and catalog sources to check for correspondence of a name or position
    """

    max_sep = 0.1

    def __init__(self, config):
        from astropy.io import fits
        import pandas as pd

        zip_index = get_wtzip_index(config)
        self.pt_dirs=zip_index['coord']
        self.pt_names = zip_index['name']

        catalog_file = Path(config.catalog_fits).expanduser()

        # make this optional
        assert(catalog_file).is_file(), f'4FGL file not found: {catalog_file}'
        hdus = fits.open(catalog_file)
        data = hdus[1].data
        self.cat_names = list(map(lambda x: x.strip() , data['Source_Name']))
        self.cat_dirs = SkyCoord(data['RAJ2000'], data['DEJ2000'], unit='deg', frame='fk5')

    def __call__(self, *pars, **kwargs):
        """
        Search the catalog lists. Options are:

        * name of a pointlike source
        * name of a source found by astropy, or a coordinate, which is close to a source in the pointlike list
        * a coordinate pair (ra,dec), or (l,b, gal=True)

        Returns the pointlike name. If  with the latter `None` if not in 4FGL.
        """
        self.psep=self.csep=99 # flag not found

        # first check pointlike list
        if len(pars)==1 and pars[0] in self.pt_names:
            idx = list(self.pt_names).index(pars[0])
            skycoord = self.pt_dirs[idx]
        else:
            # get coord either by known catalog name, or explict coordinate pair
            try:
                skycoord = findsource(*pars, **kwargs)
            except TypeError as err:
                print(err)
                return None
            if skycoord is None:
                error = f'*** Name "{pars}" not found by astropy, and is not in the pointlike list'
                print(error, file=sys.stderr)
                return None

        self.psep=0

        idx, sep2d, _= skycoord.match_to_catalog_sky(self.pt_dirs)
        self.psep = sep = sep2d.deg[0]
        pt_name =  self.pt_names[idx]
        if sep > self.max_sep:
            error = f'*** Name "{pars}" is {sep:.2f} deg > {self.max_sep} from pointlike source {pt_name}'

            print(error, file=sys.stderr)

            return None

        # check for 4FGL correspondence
        idx, sep2d, _= skycoord.match_to_catalog_sky(self.cat_dirs)
        self.csep = sep2d.deg[0]
        self.cat_name = self.cat_names[idx] if self.csep < self.max_sep else None
        self.skycoord = skycoord

        return pt_name

# Cell

class PointSource():
    """Manage the position and name of a point source
    """
    def __init__(self, *pars, **kwargs): # name,  position=None, nickname=None, config=None,):
        """

        """
        config = self.config = kwargs.pop('config', Config())
        lookup = SourceLookup(config)
        gal = kwargs.get('gal', False)
        self.nickname = pt_name = lookup(*pars, ** kwargs )
        if pt_name is None:
            raise Exception('Source not found')
        self.skycoord = lookup.skycoord
        #print(pars)
        if len(pars)==1:
            name = pars[0]
            if name==pt_name and lookup.cat_name is not None:
                name = lookup.cat_name
        else:
            gal = kwargs.get('gal', False)
            name=f'{"fk5" if gal else "gal"} ({pars[0]},{pars[1]}) '
        self.name = name
        gal = self.skycoord.galactic
        self.l, self.b = gal.l.deg, gal.b.deg
        self.cat_name = lookup.cat_name

        try:
            self.wtman = WeightMan(self.config, self)
            # add wtman attribute references
            self.__dict__.update(self.wtman.__dict__)
        except Exception as e:
            print(f'Unexpected Weigthman failure: {e}', file=sys.stderr)
            raise


    def __str__(self):
        return f'Source "{self.name}" at: (l,b)=({self.l:.3f},{self.b:.3f}), nickname {self.nickname}'
    def __repr__(self): return str(self)

    @property
    def ra(self):
        sk = self.skycoord.transform_to('fk5')
        return sk.ra.value
    @property
    def dec(self):
        sk = self.skycoord.transform_to('fk5')
        return sk.dec.value

    @property
    def filename(self):
        """Modified name for file system"""
        return self.name.replace(' ', '_').replace('+','p') if getattr(self,'nickname',None) is None else self.nickname

    @classmethod
    def fk5(cls, name, position):
        """position: (ra,dec) tuple """
        ra,dec = position
        sk = SkyCoord(ra, dec, unit='deg',  frame='fk5').transform_to('galactic')
        return cls(name, (sk.l.value, sk.b.value))

    @property
    def spectral_model(self):
        if not hasattr(self, 'fit_info'): return None
        modelname = self.fit_info['modelname']
        pars = self.fit_info['pars']
        if modelname=='LogParabola':
            return self.LogParabola(pars)
        elif modelname=='PLSuperExpCutoff':
            return self.PLSuperExpCutoff(pars)
        else:
            raise Exception(f'PointSource: Unrecognized spectral model name {fi["modelname"]}')

    def __call__(self, energy):
        """if wtman set, return photon flux at energy"""
        return self.spectral_model(energy) if self.spectral_model else None

    def sed_plot(self, ax=None, figsize=(5,4), **kwargs):
        """Make an SED for the source

        - kwargs -- for the Axes object (xlim, ylim, etc.)
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=figsize) if ax is None else (ax.figure, ax)
        x =np.logspace(2,5,61)
        y = self(x)
        ax.loglog(x/1e3, y*x**2 * 1e6, '-')
        ax.grid(alpha=0.5)
        kw = dict(xlabel='Energy (GeV)',
                  ylabel=r'$\mathrm{Energy\ Flux\ (eV\ cm^{-2}\ s^{-1})}$',
                  title=f'{self.name}',
                  xlim=(x.min(),x.max()),
                 )
        kw.update(kwargs)
        ax.set(**kw)

    class FluxModel():
        emin, emax = 1e2, 1e5
        def __init__(self, pars, e0=1000):
            self.pars=pars
            self.e0=e0

        def photon_flux(self):
            return quad(self, self.emin, self.emax)[0]

        def energy_flux(self):
            func = lambda e: self(e) * e**2
            return quad(func, self.emin, self.emax)[0]

    class LogParabola(FluxModel):

        def __call__(self, e):
            n0,alpha,beta,e_break=self.pars
            x = np.log(e_break/e)
            y = (alpha - beta*x)*x
            return n0*np.exp(y)

    class PLSuperExpCutoff(FluxModel):

        def __call__(self,e):
            print('WARNING: check e0!')
            n0,gamma,cutoff,b=self.pars
            return n0*(self.e0/e)**gamma*np.exp(-(e/cutoff)**b)