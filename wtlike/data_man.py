# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/01_data_man.ipynb (unless otherwise specified).

__all__ = ['get_ft1_data', 'get_ft2_info', 'FermiData', 'check_data', 'update_data', 'plot_week']

# Cell
import os, sys
import dateutil, datetime
from astropy.io import fits
from ftplib import FTP_TLS as FTP

import healpy
from pathlib import Path
import numpy as np
import pickle

from .config import Config, Timer, UTC, MJD

# Cell
def get_ft1_data( config, ft1_file):

        """
        Read in a photon data (FT1) file, bin in energy and position to convert to a compact table

        - `ft1_file` -- A weekly file from GSFC

        Depends on config items
        - `theta_cut, z_cut` -- selection criteria
        - `ebins, etypes` -- define band index
        - `nside, nest` -- define HEALPix binning

        Returns a dict with keys

        - `tstart`, the start MET time from the FT1 header

        - `photons`: a dict with keys and entries for each selected photon
           - `band` (uint8):    energy band index*2 + 0,1 for Front/Back
           - `nest_index`  if nest else `ring_index` (uint32): HEALPIx index for the nside
           - `run_ref` (uint8) reference to the run number, in the array `runs`
           - `trun` (unit32): time since the run id in 2 $\\mu s$ units

        - `gti_times` -- GTI times as an interleaved start, stop array.
        - `runs` -- a list of the run numbers, each a MET time. Expect 109 per week

        For the selected events above 100 MeV, this represents 10 bytes per photon, vs. 27 in the FT1 data
        """

        delta_t = config.offset_size
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
            raise Exception(f'Non-monatonic GTI found in file {ft1_file}')

        # def apply_gti(time):
        #     x = np.digitize(time, gti_times)
        #     return np.bitwise_and(x,1).astype(bool)

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


        if verbose>2:
            total = sum(b)-sum(a)
            fraction = total/(b[-1]-a[0])

            print(  f'FT1: {ft1_file.name}, GTI range {a[0]:.1f}-{b[-1]:.1f}:  {len(data):,} photons'\
                    f'\n\tSelection E > {ebins[0]:.0f} MeV. theta<{theta_cut:.1f} and z<{z_cut} remove:'\
                    f' {100.- 100*len(dsel)/float(len(data)):.2f}%'
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

        #
        run_id = dsel['RUN_ID'].astype(np.uint32)

        # save this for for run reference
        runlist = np.unique(run_id)

        # a dict with all photon info, which requires runlist
        runs = np.unique(run_id)
        photons = {'band' : band_index,
                  hpname : hpindex,
                  'run_ref': np.searchsorted(runs, run_id).astype(np.uint8),
                  'trun'  : ((time-run_id)/delta_t).astype(np.uint32),
                 }

        if verbose>1:
            print(f'FT1 data from file {ft1_file.name}: tstart={tstart:.0f} (UTC {UTC(MJD(tstart))[:-6]})'
                  f' selected {len(dsel):,}/{len(data):,} photons, in {len(runlist)} runs.')

        return  dict(tstart=tstart, #df,
                     photons=photons,
                     gti_times=gti_times,
                     runlist=runlist)

# Cell
def get_ft2_info(config, filename,
                 gti=lambda t: True):
    """Process a FT2 file, with S/C history data, and return a summary dict

    Parameters:

    * config -- verbose, cos_theta_max, z_max
    * filename -- spacecraft (FT2) file
    * gti -- GTI object that checkes for allowed intervals, in MJD units

    Returns: A dict with fields consistent with GTI if specified

    * start, stop -- interval in MJD units
    * livetime -- sec
    * ra_scz, dec_scz --spaceraft direction
    * ra_zenith, dec_zenith -- local zenith
    """
    # combine the files into a dict  with following fields besides START and STOP (lower case for column)
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


    sc_data = dict(t)

    return sc_data

# Internal Cell
def filepaths(week):
    """Returns: A tuple with two elements for the week number, each with two triplets with:
        ftp folder, ftp filename, local simple filename
    """
    urls = []
    for ftype, alias in  [('spacecraft','ft2'), ('photon','ft1')]:
         urls.append((
             f'{ftype}',
             f'lat_{ftype}_weekly_w{week:03d}_p{"305" if ftype=="photon" else "310" }_v001.fits',
             f'week{week:03d}_{alias}.fits',
            ))
    return urls

# Internal Cell
class GSFCweekly(dict):

    ftp_site = 'heasarc.gsfc.nasa.gov'
    ftp_path = 'fermi/data/lat/weekly'
    local_path  = '/tmp/from_gsfc'

    def __init__(self, config=None):
        """ Obtain lists of the weekly FT1 and FT2 files at GSFC, Set up as a dict, with
        keys= week numbers, values=mofification date strings
        """
        self.config = config or Config()
        # self.wtlike_data_file_path = Path(self.config.datapath/'data_files')
        # assert self.wtlike_data_file_path.is_dir(), 'Data path invalid'
        # os.makedirs(self.local_path, exist_ok=True)
        try:
            with FTP(self.ftp_site) as ftp:
                ftp.login()
                ftp.prot_p()
                ftp.cwd(self.ftp_path+'/photon') # or spacecraft
                # aet modification time and type for all files in folder

                parse_week = lambda fn:  int(fn.split('_')[3][1:])
                flist = ftp.mlsd(facts=['modify', 'type'])
                self.fileinfo = sorted([(parse_week(name), fact['modify']) for name,fact in flist
                                 if fact['type']=='file' and name.startswith('lat') ])
        except Exception as msg:
            print(f'FTP login to or download from {self.ftp_site} failed:\n\t--> {msg}',file=sys.stderr)
            self.valid=False
            return
        self.update(self.fileinfo)
        self.valid=True

    def download(self, week):
        """ Download the given week's files to the local folder

        week -- the mission week number, starting at 9. If negative, get recent one

        return the ft1, ft2 local filenames
        """

        if week<0:
            week = list(self.keys())[week]
        assert week in self, f'week {week} not found at FTP site'
        files = []
        with FTP(self.ftp_site) as ftp:
            ftp.login()
            ftp.prot_p()
            for ftp_folder, ftp_filename, local_filename in filepaths(week):
                ftp.cwd('/'+self.ftp_path+'/'+ftp_folder)
                if self.config.verbose>0:
                    print(f'GSFCweekly: {ftp_folder}/{ftp_filename} --> {local_filename}', flush=True)
                with open(f'{self.local_path}/{local_filename}', 'wb') as  localfile:
                    ftp.retrbinary('RETR ' + ftp_filename, localfile.write)
                files.append(local_filename)
        return files

    def week_number(self, met):
        return (met-233711940)//(7*24*3600)

    def load_week(self, week):
        """Load a pickled week summary """
        filename = self.wtlike_data_file_path/f'week_{week:03d}.pkl'
        assert filename.exists(), f'File {filename} does not exist'
        with open(filename, 'rb') as imp:
            ret = pickle.load(imp)
        return ret

    def check_week(self, week):
        """Returns True if the local week needs updating"""
        data = self.load_week(week)
        if 'file_date' not in data:
            return True
        # check that file date agrees
        return data['file_date'] != self[week]



# Cell
class FermiData(GSFCweekly):
    """ Manage the full data set in weekly chunks
    * Checking the current set of files at GSFC
    * downloading a week at a time to a local tmp
    * Converting to condensed format and saving to pickled dicts in wtlike_data
    """

    def __init__(self, config=None):
        super().__init__(config)
        self.wtlike_data_file_path = Path(self.config.datapath/'data_files')
        assert self.wtlike_data_file_path.is_dir(), 'Data path is invalid'
        os.makedirs(self.local_path, exist_ok=True)

    @property
    def local_filedate(self):
        """ the datetime object representing the last file date in local storage"""
        from dateutil.parser import parse
        weekly_folder = self.config.datapath/'data_files'
        ff = sorted(list(weekly_folder.glob('*.pkl')))
        if len(ff)==0:
            print(f'FermiData: No weekly summary files found in {weekly_folder}', file=sys.stderr)
            return None

        wk = list(map(lambda f: int(os.path.splitext(f)[0][-3:]), ff))
        lastweek = pickle.load(open(ff[-1],'rb'))
        return dateutil.parser.parse(lastweek['file_date'])

    @property
    def gsfc_filedate(self):
        return dateutil.parser.parse(list(self.values())[-1])

    def __str__(self):
        return f'FermiData: {len(self.fileinfo)} week files at GSFC, from {self.fileinfo[0]} to {self.fileinfo[-1]}'

    def in_temp(self):
        """return list of GSFC copied files in the local_path folder"""
        names = [f.name for f in Path(self.local_path).glob('*')]
        return names

    def __call__(self, week, test=False, tries_left=2):
        """ Download and convert the given week:
        * download FT1 and FT2 from GSFC to scratch space
        * convert each
        * save pickled dict summary
        * remove files

        """
        assert week in self, f'week {week} not found at FTP site'
        ff = filepaths(week)
        ft1_file = Path(self.local_path)/ff[1][2]
        ft2_file = Path(self.local_path)/ff[0][2]

        if self.config.verbose>1:
            print(f'FermiData: converting week {week}')

        while tries_left>0:
            try:
                if not (ft1_file.exists() and ft2_file.exists()):
                    self.download(week)
                # convert photon data to compact form, as a dice
                week_summary = get_ft1_data(self.config, ft1_file)
                break
            except Exception as e:
                print(f'*** ERROR *** Failed to convert {ft1_file}: {e} download it again)',
                      file=sys.stderr)
                tries_left -=1
                if tries_left==0:
                    print(f'Failed to convert week file {ft1_file}: quitting', file=sys.stderr)
                    return None
                else:
                    os.unlink(ft1_file)

        def apply_gti(time): # note MJD
            x = np.digitize(time, MJD(week_summary['gti_times']))
            return np.bitwise_and(x,1).astype(bool)
        sc_data = get_ft2_info(self.config, ft2_file, apply_gti)

        # finished with copies of FT1 and FT2 files: delete them
        for file in (ft1_file,ft2_file):
            os.unlink(file)

        # add file date and space craft summary
        week_summary['file_date'] = self[week]
        week_summary['sc_data'] = sc_data
        filename = self.wtlike_data_file_path/f'week_{week:03d}.pkl'
        if filename.exists() and self.config.verbose>0:
            print(f'FermiData: replacing existing {filename}',flush=True)
        if not test:
            with open(filename, 'wb') as out:
                pickle.dump(week_summary, out)

        if self.config.verbose>0:
            print(f'FermiData: Saved to {filename}', flush=True)

    def download_and_convert(self, week_range,  processes=None):
        """Download FT1 and FT2 files from GSFC and create summary files for the weeks

        * week_range: a (first, last+1) tumple, or a iterable

        """
        from multiprocessing import Pool

        processes = processes or self.config.pool_size
        txt = f', using {processes} processes ' if processes>1 else ''

        if self.config.verbose>0:
            print(f'\tDownloading {len(week_range)} week files{txt}\n', end='', flush=True)

        if processes>1:
            with Pool(processes=processes) as pool:
                pool.map(self, week_range)
        else:
            list(map(self,  week_range))

    def needs_update(self):
        """ Compare files on disk with the GSFC list and compile list that need to be downloaded

        Check the file date of the last one on disk and update it as well if it has a different filedate
        """
        gg =self.wtlike_data_file_path.glob('*.pkl')
        file_weeks= map(lambda n: int(n.name[5:8]), gg)
        ondisk = np.array(list(file_weeks))

        missing =  set(self.keys()).difference(set(ondisk))
        if len(ondisk)==0:
            return list(missing)

        last = sorted(ondisk)[-1]
        if self.check_week(last):
             missing.add(last)
        return sorted(list(missing))

    def check_data(self):
        """
        Return: sorted list of summary files, last week number, number of days in last week"""
        config = self.config
        if config.valid:
            weekly_folder = config.datapath/'data_files'
            ff = sorted(list(weekly_folder.glob('*.pkl')))
            if len(ff)==0:
                print(f'No .pkl files found in {weekly_folder}', file=sys.stderr)
                return
            getwk = lambda f: int(os.path.splitext(f)[0][-3:])
            wk = [getwk(f) for f in ff] #list(map(lambda f: int(os.path.splitext(f)[0][-3:]), ff))
            lastweek = pickle.load(open(ff[-1],'rb'))

            file_date = lastweek['file_date']
            gti = lastweek['gti_times'];
            days = (gti[-1]-gti[0])/(24*3600)
            if config.verbose>0:
                print(f'Weekly folder "{weekly_folder}" contains {len(wk)} weeks.'\
                      f'\n\t Last week in local dataset, #{wk[-1]}, has {days:.3f} days, ends at UTC {UTC(MJD(gti[-1]))}, filedate {file_date}' )
            #return ff, wk[-1], days
        else:
            print(f'Config not valid, {config.errors}', file=sys.stderr)
            return None

    def update_data(self):
        """Bring all of the local week data summaries up to date, downloading the missing ones from GSFC.
        If the last one is the current week, check to see if needs updating, by comparing file date, in days,
        from the last update with the current one at GSFC.
        """
        self.check_data()
        needs = self.needs_update()
        if len(needs)==0:
            print('--> Up to date!')
            return
        return self.download_and_convert(needs)

    def get_run_times(self, week):
        r = self.load_week(week)
        pdict = r['photons']
        if 'run_id' in pdict:
            runs = np.unique(pdict['run_id'])
        elif 'run_ref' in pdict:
            runs = r['runlist'][pdict['run_ref']]
        else:
            assert False
        return MJD(runs)

# Cell
def check_data(config=None, update=False):
    """
    Print current status, and update if requested
    """
    ff = FermiData(config)
    ff.check_data()
    print(f'Weeks needing download: {ff.needs_update()}')
    if update:
        ff.update_data()

def update_data( config=None):
    """Bring all of the local week data summaries up to date, downloading the missing ones from GSFC.
    If the last one is the current week, check to see if needs updating, by comparing file date, in days,
    from the last update with the current one at GSFC.
    """
    ff = FermiData(config)
    return ff.update_data()

# Internal Cell
def get_week_files(config, week_range=None):
    """Return list of week files

    - week_range [None] -- tuple with inclusive range. If None, get all
    """
    import pandas as pd
    data_folder = config.datapath/'data_files'
    data_files = sorted(list(data_folder.glob('*.pkl')))
    weeks = week_range or  config.week_range
    if week_range is not None:

        slc = slice(*week_range) if type(week_range)==tuple else slice(week_range,week_range)
        wk_table = pd.Series(data=[df for df in data_files],
                     index= [ int(df.name[-7:-4]) for df in  data_files],
                    )
        data_files = wk_table.loc[slc].values

        if config.verbose>1:
            q = lambda x: x if x is not None else ""
            print(f'LoadData: Loading weeks[{q(slc.start)}:{q(slc.stop)}:{q(slc.step)}]', end='' if config.verbose<2 else '\n')
    else:
        if config.verbose>1: print(f'LoadData: loading all {len(data_files)} weekly files')

    if len(data_files)==0:
        msg =  f'Specified week_range {week_range} produced no output. Note that week numbers are 9-'
        raise Exception(msg)

    return data_files

# Cell
def plot_week(week=None, mjd=None, nside=32, **kwargs):
    """
    Make an AIT plot of the given week's photon data

    Combine all energies for now

    - week -- the week number from 9
    - mjd -- [None] If set, derive the week from it
    - nside [32] -- HEALPix nside to project data before plotting.
    - kwargs -- args for healpix.ait_plot

    """
    import matplotlib.pyplot as plt
    from utilities import healpix as hpx
    from .config import mission_week

    assert nside & (nside-1) == 0, 'nside must be power of 2'
    config = Config()
    if not config.valid:
        print('No data to plot since config not valid', file=sys.stderr)
        return
    kw = dict(log=True, tick_labels=False, vmin=5, vmax=1e3,
             cblabel=f'counts per nside={nside} pixel')
    kw.update(kwargs)
    if mjd is not None: week=mission_week(mjd)

    file = get_week_files(config,(week,week))[0]
    with open(file, 'rb') as inp:
        u = pickle.load(inp)
    ni = u['photons']['nest_index']
    utc = UTC(MJD((u['gti_times'][0])))[:-5]
    n = 12*nside**2
    # to convert data nside (1024 normally)
    shift = int(np.log2(config.nside/nside))
    pmap,_ = np.histogram( np.right_shift(ni,2*shift), np.linspace(0,n,n+1))
    t = hpx.HPmap(pmap, f'week_{week:d}\n{utc}', nest=True)
    t.ait_plot( **kw);