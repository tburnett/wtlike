# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/01-data_man.ipynb (unless otherwise specified).

__all__ = ['get_ft1_data', 'ArrowDataset', 'get_photons_near_source', 'get_ft2_info']

# Cell
from astropy.io import fits
import healpy
import numpy as np
import pandas as pd
pd.set_option('display.precision', 2)
import pickle
import pyarrow as pa
import pyarrow.parquet as pq

from .config import *

def get_ft1_data( config, ft1_file):

        """
        Read in a photon data (FT1) file, bin in energy and position to convert to a compact DataFrame

        - `ft1_file` -- A monthly file generated by J. Ballet

        Depends on config items
        - `theta_cut, z_cut` -- selection criteria
        - `ebins, etypes` -- define band index
        - `nside, nest` -- define HEALPix binning

        Returns a tuple with

        - `tstart`, the start MET time,  and

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

        ft1 = fits.open(ft1_file)
        tstart = ft1[1].header['TSTART']

        ## GTI - setup raveled array function to make cut
        gti_data= ft1['GTI'].data
        a,b = gti_data.START, gti_data.STOP
        gti_times = np.ravel(np.column_stack((a,b)))
        def apply_gti(time):
            x = np.digitize(time, gti_times)
            return np.bitwise_and(x,1).astype(bool)


        # extract arrays for values of interest
        data =ft1['EVENTS'].data
        glon, glat, energy, et, z, theta, time =\
             [data[x] for x in 'L B ENERGY EVENT_TYPE ZENITH_ANGLE THETA TIME'.split()]

        # generate event_type masks
        et_mask={}
        for ie in etypes:
            et_mask[ie]= et[:,-1-ie]

        data_cut = np.logical_and(theta<theta_cut, z<z_cut)
        e_cut = energy>ebins[0]
        gti_cut = apply_gti(time)

        if verbose>0:
            print(  f'FT1 file {ft1_file}:'\
                    f'\n\tFound {len(data):,} events, {sum(e_cut):,} with E>{ebins[0]:.0f} MeV.'\
                    f'\n\ttheta<{theta_cut:.1f} and z<{z_cut} selections remove:'\
                    f' {100.- 100*sum(data_cut)/float(len(data)):.2f} %'\
                    f'\n\tGTI cut removes {sum(~gti_cut)}'
                 )
        # apply selection
        sel = e_cut & data_cut & gti_cut
        glon_sel = glon[sel]
        glat_sel = glat[sel]

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
                     (time-tstart)[sel].astype(np.float32) ],
                    names=['band', hpname,'time'])
        if verbose>1:
            print(f'\tReturning tstart={tstart:.0f}, {len(recarray):,} photons.')



        return  tstart, pd.DataFrame.from_records(recarray), gti_times



# Cell
class ArrowDataset():
    """
    Encapsulate the  Parquet Arrow monthly dataset

    Assumes that the dataset folder is at `config.files.data`
    """

    def __init__(self, config):
        """
        Represent the  Parquet Arrow monthly dataset

        """
        # set up, assuming folder has a file `tstart.pkl` and a folder `dataset`
        self.dataroot = config.files.data
        assert self.dataroot.exists()
        self.tstarts = pickle.load(open( self.dataroot/'tstart.pkl', 'rb'))
        self.dataset = self.dataroot/'dataset'
        assert self.dataset.exists()
        if len(self.tstarts)==0:
            # empty. Is there a better way?
            self.pqset=None
            self.pdict={}
            self.month_numbers=[]
        else:
            self.pqset = t = pq.ParquetDataset(self.dataset)
            assert 'month' in t.partitions.partition_names
            u =t.partitions[0]
            self.pdict = u.dictionary
            self.month_numbers = sorted(list(map(lambda x: int(str(x)), self.pdict)))

    def __repr__(self):
        n = self.month_numbers
        s = f'Parquet DataSet at "{self.dataroot}" '
        if len(n)>0:
            s += f'has {len(self.pdict)} months, {n[0]}-{n[-1]}.'
        else:
             s+= 'is empty'
        return s

    def get_month(self, month ):
        """Retrieve the month's tstart, and its dataframe,
             with columns: band, time, nest_index

        """
        try:
            table = pq.read_pandas(self.dataset, filters=[f'month == {month}'.split()])
        except Exception as e:
            print(f'Failed to read month {month} from dataset {self.dataset}: {e}', file='sys.stderr')
            raise
        df = table.to_pandas()

        return self.tstarts[month],df

    def add_month(self, month, tstart, df, gti_times, overwrite=False):
        """
        Add the data for a month to the dataset

        - `month` -- month index, 1-based
        - `tstart` -- MET TSTART for the month
        - `df` -- DataFrame derived from a FT1 file. Expect columns band, hpindex (or nest_index), time
        - `gti_times` -- interleaved start/stop array
        """
        if month in self.tstarts and not overwrite:
            raise Exception(f'Month {month} already in dataset: specify "overwrite" ')
        df.loc[:, 'month'] = np.uint8(month)

        if 'nest_index' not in df.columns:
            hpname = 'ring_index'
            if hpname not in df.columns: hpname = 'hpindex'
            df.loc[:,'nest_index'] = healpy.ring2nest(config.nside, df[hpname]).astype(np.int32)


        table = pa.Table.from_pandas(df, preserve_index=False)
        pq.write_to_dataset(table, root_path=str(self.dataset), partition_cols=['month'])
        # update tstart
        self.tstarts[month] = tstart
        pickle.dump(self.tstarts, open( self.dataroot/'tstart.pkl', 'wb'))

        #update gti
        pickle.dump(gti_times, open(self.dataroot/f'gti/month_{month:}.pkl', 'wb'))


    def __getitem__(self, index):
        """
        """
        m = self.month_numbers[index]
        return self.get_month(m)

    def __len__(self):
        return len(self.month_numbers)

# Cell
#
def get_photons_near_source(config, source, dataset_part):
    """
    Select the photons near a source

    - source : a PointSource object
    - dataset_part : a partition (month) of the Arrow photon dataset

    Returns a DF with
    - `band` index,
    - `time` in MJD (added tstart and converted from MET)
    - `pixel` index, nest indexing
    - `radius` distance in deg from source direction
    """

    def _cone(config, source, nest=True):
        # cone geometry stuff: get corresponding pixels and center vector
        l,b,radius = source.l, source.b, config.radius
        cart = lambda l,b: healpy.dir2vec(l,b, lonlat=True)
        conepix = healpy.query_disc(config.nside, cart(l,b), np.radians(radius), nest=nest)
        center = healpy.dir2vec(l,b, lonlat=True)
        return center, conepix

    center, conepix = _cone(config,source)
    tzero, df = dataset_part
    allpix = df.nest_index.values

    # select by comparing high-order pixels (faster)
    shift=11
    a = np.right_shift(allpix, shift)
    c = np.unique(np.right_shift(conepix, shift))
    incone = np.isin(a,c)

    if config.verbose>2:
        a, b = sum(incone), len(allpix)
        print(f'Processing photons from month {df.month[0]}:\n\tPixel cone cut: select {a} from {b} ({100*a/b:.1f}%)')

    # cut df to entries in the cone
    dfc = df[incone]

    # distance from center for all accepted photons
    ll,bb = healpy.pix2ang(config.nside, dfc.nest_index,  nest=True, lonlat=True)
    cart = lambda l,b: healpy.dir2vec(l,b, lonlat=True)
    t2 = np.degrees(np.array(np.sqrt((1.-np.dot(center, cart(ll,bb)))*2), np.float32))
    in_cone = t2<config.radius

    if config.verbose>2:
        print(f'\tGeometric cone cut: select {sum(in_cone)}')
    # assume all in the GTI (should check)

    # times: convert to double, add to start, convert to MJD
    time = MJD(np.array(dfc.time, float)+tstart)

    # assemble the DataFrame, remove those outside the radius
    out_df = pd.DataFrame(np.rec.fromarrays(
        [np.array(dfc.band), time, dfc.nest_index, t2],
        names='band time pixel radius'.split()))[in_cone]
    return out_df

# Cell
def get_ft2_info(config, ft2_file_path, gti=None):
    """Process a set of FT2 files, with S/C history data

    Parameters:

        - config -- verbose, cos_theta_max, z_max
        - ft2_files -- list of spacecraft files
        - gti -- GTI object with allowed intervals or None

     """
    # combine the files into a DataFrame with following fields besides START and STOP (lower case for column)
    fields    = ['LIVETIME','RA_SCZ','DEC_SCZ', 'RA_ZENITH','DEC_ZENITH']
    ft2_files = list(ft2_file_path.glob('*.fits'))
    if config.verbose>1:
        print(f'Processing {len(ft2_files)} S/C history (FT2) files in {ft2_file_path}')
#         print(f'  applying cuts cos(theta) < {config.cos_theta_max},  z < {config.z_max}')
    sc_data=[]
    for filename in ft2_files:
        with fits.open(filename) as hdu:
            scdata = hdu['SC_DATA'].data
            # get times to check against MJD limits and GTI
            start, stop = [MJD(np.array(scdata.START, float)),
                           MJD(np.array(scdata.STOP, float))]
            if config.mjd_range is not None:
                a,b=  config.mjd_range
                if stop[-1]<a:
                    print(f'\tfile {filename.name}: skip, before selected range' )
                    continue
                elif start[0]>b:
                    print(f'\tfile {filename.name}: quit, beyond range')
                    break
            # apply GTI to bin center (avoid edge effects?)
            in_gti = gti(0.5*(start+stop)) if gti else np.ones(len(start), bool)
            if config.verbose>2:
                gti_check = f', {sum(in_gti)} in GTI' if gti  else ''
                print(f'\tfile {filename.name}: {len(start)} entries {gti_check}')
            t = [('start', start[in_gti]), ('stop',stop[in_gti])]+\
                [(field.lower(), np.array(scdata[field][in_gti],np.float32)) for field in fields ]

            sc_data.append( pd.DataFrame(dict(t) ) )

    return pd.concat(sc_data, ignore_index=True)