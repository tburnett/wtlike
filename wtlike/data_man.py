# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/01_data_man.ipynb (unless otherwise specified).

__all__ = ['get_ft1_data', 'get_ft2_info', 'WeeklyData']

# Cell
from astropy.io import fits
import wget
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
           - `time` (float32):    the elapsed time in s from header value TSTART in the FT1 file

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
        if verbose>1:
            total = sum(b)-sum(a)
            fraction = total/(b[-1]-a[0])
            print(f'GTI times: {a[0]:.0f} - {b[-1]:.0f}'\
                  f' Total: {total:.0f}s, fraction {100*fraction:.1f}%'
                 )

        glon, glat, energy, et, z, theta, time, ec =\
             [data[x] for x in 'L B ENERGY EVENT_TYPE ZENITH_ANGLE THETA TIME EVENT_CLASS'.split()]

        # generate event_type masks
        et_mask={}
        for ie in etypes:
            et_mask[ie]= et[:,-1-ie]



        data_cut = np.logical_and(theta<theta_cut, z<z_cut)
        e_cut = energy>ebins[0]
        gti_cut = apply_gti(time)

        if verbose>1:
            print(  f'Photon (FT1) file {ft1_file}: {len(data):,} photons, {sum(e_cut):,} with E>{ebins[0]:.0f} MeV.'\
                    f'\n\ttheta<{theta_cut:.1f} and z<{z_cut} selections remove:'\
                    f' {100.- 100*sum(data_cut)/float(len(data)):.2f} %'\
                    f'\n\tGTI cut removes {sum(~gti_cut)}'
                 )
        # apply selection
        sel = e_cut & data_cut & gti_cut
        glon_sel = glon[sel]
        glat_sel = glat[sel]

        # event class -- turn into single int for later mask
#         bits = np.array([1<<n for n in range(20)])
#         def to_bin(x):
#             return np.sum(bits[x[:20]])
#         ec = [to_bin(row[20]) for row in ec[sel]]

        # pixelate direction
        hpindex = healpy.ang2pix(nside, glon_sel, glat_sel, nest=nest, lonlat=True).astype(np.int32)
        hpname = 'nest_index' if nest else 'ring_index'

        # digitize energy and create band index incluing (front/back)
        ee = energy[sel]
        band_index = (2*(np.digitize(ee, ebins, )-1) + et_mask[1][sel]).astype(np.uint8)

        # combine into a recarray to feed to pandas
        recarray = np.rec.fromarrays(
                    [ band_index,
                     hpindex,
                     (time-tstart)[sel].astype(np.float32),
                    # ec, ### temp
                    ],
                    names=['band', hpname, 'time',
                           #'ec', #### temp
                          ])
        if verbose>1:
            print(f'\tReturning tstart={tstart:.0f}, {len(recarray):,} photons.')

        return  tstart, pd.DataFrame.from_records(recarray), gti_times

# Cell
def get_ft2_info(config, filename, gti=None):
    """Process a FT2 file, with S/C history data, and return a summary DataFrame

    Parameters:

        - config -- verbose, cos_theta_max, z_max
        - filename -- spacecraft (FT2) files
        - gti -- GTI object with allowed intervals or None

     Returns:
        A DataFrame with fields
            - start, stop -- interval in MJD units
            - livetime -- sec
            - ra_scz, dec_scz --spaceraft direction
            - ra_zenith, dec_zenith -- local zenith
     """
    # combine the files into a DataFrame with following fields besides START and STOP (lower case for column)
    fields    = ['LIVETIME','RA_SCZ','DEC_SCZ', 'RA_ZENITH','DEC_ZENITH']
    if config.verbose>1:
        print(f'S/C history (FT2) file {filename}', end='')# {"not" if gti is None else ""} applying GTI')

    with fits.open(filename) as hdu:
        scdata = hdu['SC_DATA'].data

    # get times to check against MJD limits and GTI
    start, stop = [MJD(np.array(scdata.START, float)),
                   MJD(np.array(scdata.STOP, float))]

    # apply GTI to bin center (avoid edge effects?)
    in_gti = gti(0.5*(start+stop)) if gti else np.ones(len(start), bool)
    if config.verbose>1:
        gti_check = f', {sum(in_gti)} in GTI' if gti  else ''
        print(f' {len(start)} entries {gti_check}')

    t = [('start', start[in_gti]), ('stop',stop[in_gti])]+\
        [(field.lower(), np.array(scdata[field][in_gti],np.float32)) for field in fields ]

    sc_data = pd.DataFrame(dict(t) )

    return sc_data

# Cell
class WeeklyData(object):
    """Download and process weekly Fermi-LAT files from GSFC
    """

    ftp = 'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly'
    tmp = Path('/tmp/from_gsfc')

    def __init__(self, config, week, saveto):
        """
        * week: a mission week number
        * saveto: path to save the fiels

        """
        self.config= config
        self.saveto=Path(saveto)
        os.makedirs(self.saveto, exist_ok=True)
        assert week>8
        self.wk = week

        self.fits_files = [self.tmp/f'week{week:03d}_{x}.fits' for x in 'ft1 ft2'.split()]
        os.makedirs(self.tmp, exist_ok=True)

        urls = []
        for ftype in  ['photon', 'spacecraft']:
             urls.append(f'{self.ftp}/{ftype}/lat_{ftype}_weekly_w{week:03d}_p{"305" if ftype=="photon" else "202" }_v001.fits')

        for url, fname in zip(urls, self.fits_files):
            if not fname.exists():
                if config.verbose>1:
                    print(f'{url.split("/")[-1]} -> {fname}')
                wget.download(str(url), str(fname))
            else:
                if config.verbose>1: print(f'{fname} exists')

    def process_ft1(self):
        self.tstart, self.photon_data, self.gti_times = get_ft1_data(self.config, self.fits_files[0])

    def process_ft2(self):

        def apply_gti(time): # note MJD
            x = np.digitize(time, MJD(self.gti_times))
            return np.bitwise_and(x,1).astype(bool)

        self.sc_data = get_ft2_info(self.config, str(self.fits_files[1]), apply_gti)

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
