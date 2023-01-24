# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/02_weights.ipynb.

# %% auto 0
__all__ = ['get_wtzip_index', 'WeightMan']

# %% ../nbs/02_weights.ipynb 11
import os, sys, pickle, healpy, zipfile
from pathlib import Path
import numpy as np
import pandas as pd

from scipy.integrate import quad
from astropy.coordinates import SkyCoord, Angle
from .config import *

# %% ../nbs/02_weights.ipynb 12
def get_wtzip_index(config, update=False):
    """
    Return a dict with weight table info
    
    """

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

        ## update the index
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

# %% ../nbs/02_weights.ipynb 17
class WeightMan(dict):
    """ Weight Management

    * Load weight tables
    * Assign weights to photons
    """

    def __init__(self, config, source):
        """
        source -- a PointSource object
        
        Lookup is via the "nickname" property of the source object.
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

# %% ../nbs/02_weights.ipynb 20
# ?? 
class WTSkyCoord(SkyCoord):
    def __repr__(self):
        ra,dec = self.fk5.ra.deg, self.fk5.dec.deg
        return f'({ra:.3f},{dec:.3f})'
