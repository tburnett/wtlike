# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/03_sources.ipynb.

# %% auto 0
__all__ = ['get_wtzip_index', 'WeightMan', 'weight_radius_plots', 'findsource', 'FermiCatalog', 'SourceLookup', 'PointSource']

# %% ../nbs/03_sources.ipynb 7
import os, sys, pickle, healpy, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

from scipy.integrate import quad
from astropy.coordinates import SkyCoord, Angle
from .config import *

# %% ../nbs/03_sources.ipynb 8
def get_wtzip_index(config, update=False):

    wtzipfile = config.datapath/'weight_files.zip'
    if not  wtzipfile.is_file():
        print( f'Did not find the zip file {wtzipfile}', file=sys.stderr)
        return None

    with  zipfile.ZipFile(wtzipfile) as wtzip:
        if 'index.pkl' in wtzip.namelist() and not update:
            zi =  pickle.load(wtzip.open('index.pkl'))
#             if 'coord' in zi:
#                 print('Updating coord in index.pkl')

            zi['coord'] = SkyCoord(zi['glon'], zi['glat'], unit='deg', frame='galactic').fk5
            return zi

        if config.verbose>0:
            print(f'sources.get_wtzip_index: Extracting info from {wtzipfile}')
        name=[]; glat=[]; glon=[]
        for filename in wtzip.namelist():
            if filename=='index.pkl': continue
            with wtzip.open(filename) as file:
                wtd = pickle.load(file, encoding='latin1')
                l,b = wtd['source_lb']
                name.append(Path(filename).name.split('_weights.pkl')[0].replace('_',' ').replace('p','+') )
                glon.append(l)
                glat.append(b)
        zip_index = dict(name=np.array(name),
                         glon=np.array(glon), glat=np.array(glat),
                         #coord= SkyCoord(glon, glat, unit='deg', frame='galactic'),
                        )

        ### write to temp file, insert back into the zip
        ### SHould be a way to just stream
        pickle.dump(zip_index, open('/tmp/wtfile_index.pkl', 'wb'))
        with zipfile.ZipFile(wtzipfile, mode='a') as wtzip:
            wtzip.write('/tmp/wtfile_index.pkl', 'index.pkl')

        zip_index['coord'] = SkyCoord(zip_index['glon'], zip_index['glat'], unit='deg', frame='galactic').fk5
    return zip_index

# %% ../nbs/03_sources.ipynb 9
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
            with open(datapath/filename, 'rb') as inp:
                wtd =  pickle.load(inp, encoding='latin1')

        elif (datapath/'weight_files.zip').exists():
            with  zipfile.ZipFile(datapath/'weight_files.zip') as wtzip:
                wtd =  pickle.load(wtzip.open(filename), encoding='latin1')

        else:
            raise Exception(f'No weight info found for {nickname}')

        self.update(wtd)
        self.__dict__.update(wtd)
        self.filename=filename
        self.config = config

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
        if not getattr(self.config, 'keep_pixels', False):
            photons.drop(['nest_index'], axis=1, inplace=True)
            if self.config.verbose>2:
                print('Keeping pixels')
        noweight = np.isnan(photons.weight.values)
        if self.config.verbose>1:
            print(f'\tremove {sum(noweight):,} events without weight')

        ret = photons[np.logical_not(noweight)]
        assert ret is not None
        return ret

# %% ../nbs/03_sources.ipynb 11
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

# %% ../nbs/03_sources.ipynb 12
def findsource(*pars, gal=False):
    """
    Return a SkyCoord, looking up a name, or interpreting args

    Optional inputs: len(pars) is 1 for a source name or Jxxx, or 2 for coordinate pair
    - name -- look up the name, return None if not found
    - Jxxxx.x+yyyy -- intrepret  to get ra, dec, then convert
    - ra,dec -- assume frame=fk5
    - l,b, gal=True -- assume degrees, frame=galactic
    """

    import astropy.units as u
    if len(pars)==1:
        name = tname= pars[0]
        if name.startswith('4FGL '):
            # SkyCoord.from_name not updated for DR3. Use as Jname, ignoring that 'c'
            name = name[4:].replace('c','')
        if name.startswith('J') and (len(name)==10 or len(name)==12) and ('+' in name or '-' in name):
            # parse the name for (ra,dec)
            if name[5]!='.':
                tname = name[:5]+'.0'+name[5:]
            ra=(tname[1:3]+'h'+tname[3:7]+'m')
            dec = (tname[7:10]+'d'+tname[10:12]+'m')
            try:
                (ra,dec) = map(lambda a: float(Angle(a, unit=u.deg).to_string(decimal=True)),(ra,dec))
                skycoord = SkyCoord(ra, dec, unit='deg', frame='fk5')
            except ValueError:
                print(f'Attempt to parse {name} failed: expect "J1234.5+6789" or "J1234.5678"', file=sys.stderr)
                return None
#         elif name.startswith('J'):
#             print(f'The name "{name}" starts with a "J", but cannot be parsed for ra,dec', file=sys.stderr)
#             return None

        
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
    return skycoord.galactic if gal else skycoord

# %% ../nbs/03_sources.ipynb 15
class WTSkyCoord(SkyCoord):
    def __repr__(self):
        ra,dec = self.fk5.ra.deg, self.fk5.dec.deg
        return f'({ra:.3f},{dec:.3f})'

# %% ../nbs/03_sources.ipynb 16
import warnings
class FermiCatalog():

    def __init__(self,config=None, max_sep=0.1):
        from astropy.io import fits
        self.max_sep=max_sep
        if config is None:
            config = Config()
        fail=False
        # look for 4FGL catalog file, gll_psc_v28.fit currently
        if config.catalog_file is None:
                fail = True
        elif Path(config.catalog_file).is_file():
            self.catalog_file = Path(config.catalog_file).expanduser()
        else: fail=True

        if fail:
            warnings.warn('There is no link to 4FGL catalog file: set "catalog_file" in your config.yaml'
                  ' or specify if in the Config() call')
        else:
            # make this optional
            config.catalog_file = self.catalog_file
            with fits.open(self.catalog_file) as hdus:
                self.data = data = hdus[1].data

            cname = lambda n : [s.strip() for s in data[n]]
            cvar = lambda a: data[a].astype(float)
            ivar = lambda a: data[a].astype(int)
            name = list(map(lambda x: x.strip() , data['Source_Name']))

            self.skycoord = WTSkyCoord(data['RAJ2000'], data['DEJ2000'], unit='deg', frame='fk5')

            self.df = pd.DataFrame(dict(

                skycoord = self.skycoord,
                significance = cvar('Signif_Avg'),
                variability = cvar('Variability_Index'),
                assoc_prob  = cvar('ASSOC_PROB_BAY'), # for Bayesian, or _LR for likelihood ratio
                assoc1_name = cname('ASSOC1'),
                class1      = cname('CLASS1'),
                flags       = ivar('FLAGS'),
                # ....
            ))
            self.df.index = name
            self.df.index.name = 'name'

    def __repr__(self):
        return f'4FGL file {self.catalog_file.name} with {len(self)} entries'

    def field(self, colname):
        ## Todo: check type
        return self.data[colname]

    def __call__(self, skycoord):
        """select an entry by skydir return entry"""
        try:
            idx, sep2d, _= skycoord.match_to_catalog_sky(self.skycoord)
        except Exception as msg:
            print(f'Fail skycoord.match, {skycoord}: {msg}', file=sys.stderr)
            return
        csep = sep2d.deg[0]
        return self.df.iloc[idx] if csep < self.max_sep else None

    def __getitem__(self, name):   return self.df.loc[name]
    def __len__(self): return len(self.df)

# %% ../nbs/03_sources.ipynb 18
class SourceLookup():
    """ Use lists of the pointlike and catalog sources to check for correspondence of a name or position
    """

    max_sep = 0.1

    def __init__(self, config):
        from astropy.io import fits
        import pandas as pd
        self.config=config

        zip_index = get_wtzip_index(config)
        if zip_index is None:
            raise Exception('Expected zip file weight_files.zip')

        self.pt_dirs=zip_index['coord']
        self.pt_names = zip_index['name']

        catalog_file = Path(config.catalog_file).expanduser()

        # make this optional
        with fits.open(catalog_file) as hdus:
            data = hdus[1].data
        self.cat_names = list(map(lambda x: x.strip() , data['Source_Name']))
        self.cat_dirs = SkyCoord(data['RAJ2000'], data['DEJ2000'], unit='deg', frame='fk5')

    def check_folder(self, *pars):
        if len(pars)>1: return None
        name = pars[0]
        filename = self.config.datapath/'weight_files'/(name.replace(' ','_').replace('+','p')+'_weights.pkl')
        if not filename.is_file():
            return None
        with open(filename, 'rb') as inp:
            wd = pickle.load(inp, encoding='latin1')

        #print(wd.keys(), wd['source_name'], wd['source_lb'])
        self.skycoord = SkyCoord( *wd['source_lb'], unit='deg', frame='galactic')
        self.check_4fgl()
        return name

    def __call__(self, *pars, **kwargs):
        """
        Search the catalog lists. Options are:

        * name of a pointlike source
        * name of a source found by astropy, or a coordinate, which is close to a source in the pointlike list
        * a coordinate pair (ra,dec), or (l,b, gal=True)

        Returns the pointlike name.
        """
        self.psep=self.csep=99 # flag not found
        self.cat_name = None

        # first, is the name in the weidht_files folder?
        name = self.check_folder(*pars)
        if name is not None:  return name

        # then check pointlike list
        if self.pt_names is not None and len(pars)==1 and pars[0] in self.pt_names:
            idx = list(self.pt_names).index(pars[0])
            skycoord = self.pt_dirs[idx]
        else:
            # get coord either by known catalog name, or explict coordinate pair
            try:
                skycoord = findsource(*pars, **kwargs)
            except TypeError as err:
                print(err)
                return None
            except ValueError as err:
                print(err, file=sys.stderr)
                return None
            if skycoord is None:
                error = f'*** Name "{pars}" not found by astropy, and is not in the pointlike list'
                print(error, file=sys.stderr)
                return None

        self.psep=0
        # print(f'looking up: {skycoord} ', end='')
        # print(f'in {self.pt_dirs}')
        try:
            idx, sep2d, _= skycoord.match_to_catalog_sky(self.pt_dirs)
        except Exception as msg:
            print(f'Fail SkyCoord catalog lookup: {self.pt_dirs}, {msg}', file=sys.stderr)
            return None

        self.psep = sep = sep2d.deg[0]
        pt_name =  self.pt_names[idx]
        if sep > self.max_sep:
            error = f'*** Name "{pars}" is {sep:.2f} deg > {self.max_sep} from pointlike source {pt_name}'

            print(error, file=sys.stderr)

            return None
        self.skycoord=skycoord
        self.check_4fgl()
        return pt_name

    def check_4fgl(self):
        # check for 4FGL correspondence, set self.cat_name, self.csep
        if len(self.cat_dirs)==0:
            self.cat_name=None
            return

        idx, sep2d, _= self.skycoord.match_to_catalog_sky(self.cat_dirs)
        self.csep = sep2d.deg[0]
        self.cat_name = self.cat_names[idx] if self.csep < self.max_sep else None

# %% ../nbs/03_sources.ipynb 21
class PointSource():
    """Manage the position and name of a point source
    """
    def __init__(self, *pars, **kwargs):
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
            name=f'{"fk5" if not gal else "gal"} ({pars[0]},{pars[1]}) '
        self.name = name
        gal = self.skycoord.galactic
        self.l, self.b = gal.l.deg, gal.b.deg
        self.cat_name = lookup.cat_name
        # override 4FGL if identifiec pulasr
        if self.nickname.startswith('PSR ') and self.cat_name is not None and self.cat_name.startswith('4FGL '):
            self.name = self.nickname
        # use cat name if used J-name to find it
        elif self.name.startswith('J') and self.cat_name is not None and self.cat_name.startswith('4FGL'):
            self.name = self.cat_name

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
                  xlim=(x.min()/1e3,x.max()/1e3),
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
            #print('WARNING: check e0!')
            n0,gamma,cutoff,b=self.pars
            return n0*(self.e0/e)**gamma*np.exp(-(e/cutoff)**b)
