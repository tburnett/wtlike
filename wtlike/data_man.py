# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/01_data_man.ipynb (unless otherwise specified).

__all__ = ['get_ft1_data', 'get_ft2_info', 'WeeklyData', 'get_data_files', 'GTI']

# Cell
import os
from astropy.io import fits

import healpy
from pathlib import Path
import numpy as np
import pandas as pd
pd.set_option('display.precision', 2)
import pickle
from .config import *

# Cell
def get_ft1_data( config, ft1_file):

        """
        Read in a photon data (FT1) file, bin in energy and position to convert to a compact DataFrame

        - `ft1_file` -- A monthly file generated by J. Ballet, or a weekly file from GSFC

        Depends on config items
        - `theta_cut, z_cut` -- selection criteria
        - `ebins, etypes` -- define band index
        - `nside, nest` -- define HEALPix binning

        Returns a tuple with

        - `tstart`, the start MET time

        - DataFrame  with columns
           - `band` (uint8):    energy band index*2 + 0,1 for Front/Back
           - `nest_index`  if nest else `ring_index` (uint32): HEALPIx index for the nside
           - `run_diff`   Run number difference from previous entry. (Saved as spase arary)
           - `time` (float32):    the elapsed time in s from header value TSTART in the FT1 file
           - `rtime` (float32): MET time relative to run id.

        - gti times as an interleaved start, stop array.

        For the selected events above 100 MeV, this represents 9 bytes per photon, vs. 27.
        """

        ebins = config.energy_edges
        etypes = config.etypes
        nside = config.nside
        nest = config.nest
        z_cut =config.z_max
        theta_cut = np.degrees(np.arccos(config.cos_theta_max))
        verbose = config.verbose

        with  fits.open(ft1_file) as ft1:
            tstart = ft1[1].header['TSTART']

            ## GTI - setup raveled array function to make cut
            gti_data= ft1['GTI'].data
            # extract arrays for values of interest
            data =ft1['EVENTS'].data

        a,b = sorted(gti_data.START), sorted(gti_data.STOP)

        gti_times = np.ravel(np.column_stack((a,b)))
        if np.any(np.diff(gti_times)<0):
            raise Exception(f'Non-monatonic GTI found')

        def apply_gti(time):
            x = np.digitize(time, gti_times)
            return np.bitwise_and(x,1).astype(bool)

        # apply  selections

        sel =  ((data['ENERGY'] > ebins[0]) &
                (data['ZENITH_ANGLE'] < z_cut) &
                (data['THETA'] < theta_cut))

        dsel = data[sel]

        # get the columns for output
        glon, glat, energy, et, z, theta, time, ec =\
             [dsel[x] for x in 'L B ENERGY EVENT_TYPE ZENITH_ANGLE THETA TIME EVENT_CLASS'.split()]

        # generate event_type masks
        et_mask={}
        for ie in etypes:
            et_mask[ie]= et[:,-1-ie]


        if verbose>1:
            total = sum(b)-sum(a)
            fraction = total/(b[-1]-a[0])

            print(  f'FT1: {ft1_file.name}, GTI range {a[0]:.1f}-{b[-1]:.1f}:  {len(data):,} photons'\
                    f'\n\tSelection E > {ebins[0]:.0f} MeV. theta<{theta_cut:.1f} and z<{z_cut} remove:'\
                    f' {100.- 100*len(dsel)/float(len(data)):.2f}%'
#                     f', GTI cut removes {sum(~gti_cut)}'
                 )


        # event class -- turn into single int for later mask
#         bits = np.array([1<<n for n in range(20)])
#         def to_bin(x):
#             return np.sum(bits[x[:20]])
#         ec = [to_bin(row[20]) for row in ec

        # pixelate direction
        hpindex = healpy.ang2pix(nside, glon, glat, nest=nest, lonlat=True).astype(np.uint32)
        hpname = 'nest_index' if nest else 'ring_index'

        # digitize energy and create band index incluing (front/back)
        band_index = (2*(np.digitize(energy, ebins, )-1) + et_mask[1]).astype(np.uint8)

        # sparcify the run, and make sparsified array of diffs
        run_id = dsel['RUN_ID']
        rundiff = np.diff(run_id, prepend=0).astype(np.int32)
        if np.any(rundiff<0):
            print(f'\t*** {sum(rundiff<0)} runs out of order: will need to sort.')
        rdspar = pd.arrays.SparseArray(rundiff)

        d ={    'band'  : band_index,
                hpname  : hpindex,
                'time'  : (time-tstart).astype(np.float32),
                'run_diff':rdspar,
                'rtime' : (time-run_id).astype(np.float32),
                ### following for testing
                #'run_id': run_id.astype(np.int32),
                #'met'   : time.astype(float),
           }
        if verbose>1:
            print(f'\tReturning tstart={tstart:.0f}, {len(dsel):,} photons.')

        return  tstart, pd.DataFrame(d), gti_times

# Cell
def get_ft2_info(config, filename,
                 gti=lambda t: True):
    """Process a FT2 file, with S/C history data, and return a summary DataFrame

    Parameters:

    * config -- verbose, cos_theta_max, z_max
    * filename -- spacecraft (FT2) file
    * gti -- GTI object that checkes for allowed intervals, in MJD units

    Returns: A DataFrame with fields consistent with GTI if specified

    * start, stop -- interval in MJD units
    * livetime -- sec
    * ra_scz, dec_scz --spaceraft direction
    * ra_zenith, dec_zenith -- local zenith
    """
    # combine the files into a DataFrame with following fields besides START and STOP (lower case for column)
    fields    = ['LIVETIME','RA_SCZ','DEC_SCZ', 'RA_ZENITH','DEC_ZENITH']
    with fits.open(filename) as hdu:
        scdata = hdu['SC_DATA'].data
        tstart, tstop = [float(hdu[0].header[field]) for field in  ('TSTART','TSTOP') ]

    if config.verbose>1:
        print(f'FT2: {filename.name}, MET range {tstart:.1f}-{tstop:.1f},', end='')# {"not" if gti is None else ""} applying GTI')

    # get times to check against MJD limits and GTI
    start, stop = [MJD(np.array(scdata.START, float)),
                   MJD(np.array(scdata.STOP, float))]

    # apply GTI to bin center (avoid edge effects?)
    in_gti = gti(0.5*(start+stop))
    if config.verbose>1:
        s = sum(in_gti)
        print(f' {len(start)} entries, {s} ({100*s/len(start):.1f}%) in GTI')

    t = [('start', start[in_gti]), ('stop',stop[in_gti])]+\
        [(field.lower(), np.array(scdata[field][in_gti],np.float32)) for field in fields ]

    sc_data = pd.DataFrame(dict(t) )

    return sc_data

# Cell
class WeeklyData(object):
    """Download and process weekly Fermi-LAT files from GSFC,

    at FTP 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly'

    * week: a mission week number
    * saveto: path to save the files

    Creates a pickled file as  dict with

    * tstart--start time (MJD)
    * photons -- DataFrame with photon data--see `get_ft1_data`
    * sc_data -- DataFrame with Spacecraft data--see `get_ft2_info`
    * gti_times -- array of interleaved start/stop times from the photon data
    """

    ftp = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly'
    tmp = Path('/tmp/from_gsfc')

    def __init__(self, config, week,  overwrite=False):
        """
        * week: a mission week number, starting at 9
        * saveto: path to save the files

        """
        import wget
        self.config= config
        self.saveto=Path(config.wtlike_data/'data_files')
        os.makedirs(self.saveto, exist_ok=True)
        assert week>8
        self.wk = week

        self.ft2_file, self.ft1_file = fits_files = [self.tmp/f'week{week:03d}_{x}.fits'
                                                     for x in 'ft2 ft1'.split()]
        os.makedirs(self.tmp, exist_ok=True)

        urls = []
        for ftype in  ['spacecraft', 'photon']:
             urls.append(f'{self.ftp}/{ftype}/lat_{ftype}_weekly_w{week:03d}_p{"305" if ftype=="photon" else "310" }_v001.fits')

        for url, fname in zip(urls, fits_files):
            if not fname.exists() or overwrite:
                if config.verbose>1:
                    print(f'{url.split("/")[-1]} -> {fname}')
                try:
                    wget.download(str(url), str(fname))
                except Exception as msg:
                    print(f'Failed to download {url}: {msg}')
                    raise
            else:
                if config.verbose>1: print(f'{fname} exists')

    def process_ft1(self):
        self.tstart, self.photon_data, self.gti_times = get_ft1_data(self.config, self.ft1_file)

    def process_ft2(self):

        def apply_gti(time): # note MJD
            x = np.digitize(time, MJD(self.gti_times))
            return np.bitwise_and(x,1).astype(bool)

        self.sc_data = get_ft2_info(self.config, self.ft2_file, apply_gti)

    def save(self):
        """process, then save aa dict"""

        self.process_ft1()
        self.process_ft2()

        d = dict(tstart = self.tstart,
                photons = self.photon_data,
                sc_data = self.sc_data,
                gti_times = self.gti_times)
        filename = self.saveto/f'week_{self.wk:03d}.pkl'
        pickle.dump(d, open(filename, 'wb'))
        if self.config.verbose>0:
            print(f'Saved to {filename}')


# Cell
def get_data_files(config):
    """
    Return a list of the pickled data files
    """
    if config.valid:
        weekly_folder = config.wtlike_data/'data_files'
        ff = sorted(list(weekly_folder.glob('*.pkl')))
        if len(ff)==0:
            print(f'No .pkl files found in {weekly_folder}', file=sys.stderr)
            return
        wk = list(map(lambda f: int(os.path.splitext(f)[0][-3:]), ff))
        lastweek = pickle.load(open(ff[-1],'rb'))
        if config.verbose>0:
            print(f'Weekly folder "{weekly_folder}" contains {len(wk)} weeks, from {wk[0]} to {wk[-1]}')
            print(f'Most recent data to UTC {UTC(MJD(lastweek["tstart"])+7)}')
        return ff
    else:
        print('Config not valid, no files found')
        return []

# Cell
class GTI(object):

    def __init__(self, gti_times):
        self.mjd_gti = MJD(gti_times)

    def __call__(self, time:'float'):
        x = np.digitize(time, self.mjd_gti)
        return np.bitwise_and(x,1).astype(bool)