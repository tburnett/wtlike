# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/00_config.ipynb.

# %% auto 0
__all__ = ['day', 'first_data', 'mission_start', 'Cache', 'Config', 'MJD', 'UTC', 'UTCnow', 'mission_week', 'mjd_range',
           'FermiInterval', 'FermiMonth', 'FermiWeek', 'bin_size_name', 'decorate_with', 'Timer']

# %% ../nbs/00_config.ipynb 3
import os, sys, warnings, pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple
import numpy as np

# %% ../nbs/00_config.ipynb 4
class Cache(dict):
    """
    Manage a file cache

    - `path` -- string or `filepath` object <br> This is the folder where the index and data files are saved.
    - `clear` -- set True to clear the cache on initialization

    This uses pickle to save objects, associated with a hashable key, which is used to index the
    filename in a file `index.pkl` in the same folder.

    The `__call__` function is a convenient way to use it, so one call may either store a new entry or retrieve an existing one.

    """

    def __init__(self, path, clear:bool=False):

        self.path = Path(path) if path else None
        if self.path is None: return
        if not self.path.exists() :
            print(f'Warning: cache Path {self.path} does not exist, cache disabled ',file=sys.stderr)
            self.path=None
            return

        self.index_file = self.path/'index.pkl'

        if self.path.exists():
            if clear:
                print('Clearing cache!')
                self.clear()
            else:
                self._load_index()

    def _dump_index(self):
        with open(self.index_file, 'wb') as file:
            pickle.dump(self, file)

    def _load_index(self):
        if not self.index_file.exists():
            self._dump_index()
            return
        with open(self.index_file, 'rb') as file:
            self.update(pickle.load(file))

    def add(self, key, object,  exist_ok=False):
        if not self.path: return
        assert type(key)==str, f'Expect key to be a string, got {key}'
        if key  in self:
            if not exist_ok:
                print(f'Warning: cached object for key "{key}" exists', file=sys.stderr)
            filename = self[key]
        else:
            filename = f'cache_file_{hex(key.__hash__())[3:]}.pkl'
            self[key] = filename
            self._dump_index()

        with open(self.path/filename, 'wb') as file:
            pickle.dump(object, file )


    def get(self, key):
        if key not in self:
            return None
        filename = self[key]
        # temp fix
        if str(filename)[0]=='/': filename=self.path/filename.name

        if not (self.path/filename).exists():
            # perhaps deleted by another instance?
            print(f'File for Cache key {key} not found, removing entry', file='sys.stderr')
            selt.pop(key)
            return None
        with open(self.path/filename, 'rb') as file:
            ret = pickle.load(file)
        return ret

    def clear(self):
        if not self.path: return
        for f in self.path.iterdir():
            if f.is_file:
                f.unlink()
        super().clear()

        self._dump_index()

    def remove(self, key):
        """remove entry and associated file"""
        if not self.path: return
        if key not in self:
            print(f'Cache: key {key} not found', file=sys.stderr)
            return
        filename = self.path/self[key]
        try:
            filename.unlink()
        except:
            print(f'Failed to unlink file {filename}', file=sys.stderr)
            raise
        super().pop(key)
        self._dump_index()


    def __call__(self, key, func, *pars, description='', overwrite=False, **kwargs,
                ):
        """
        One-line usage interface for cache use

        - `key` -- key to use, usually a string. Must be hashable <br>
            If None, ignore cache and return the function evaluation
        - `func` -- user function that will return an object that can be pickled
        - `pars`, `kwargs` -- pass to `func`
        - `description` -- optional string that will be printed
        - `overwrite` -- if set, overwrite previous entry if exists

        Example:
        <pre>
        mycache = Cache('/tmp/thecache', clear=True)

        def myfun(x):
            return x

        result = mycache('mykey', myfun, x=99,  description='My data')

        </pre>

        """

        if key is None or self.path is None:
            return func(*pars, **kwargs)


        if description:
            print(f'{description}: {"Saving to" if key not in self or overwrite else "Restoring from"} cache', end='')
            print('' if key == description else f' with key "{key}"')
        ret = self.get(key)
        if ret is None or overwrite:
            ret = func(*pars, **kwargs)
            self.add(key, ret, exist_ok=overwrite)
        return ret

    def show(self, starts_with=''):
        import datetime
        if not self.path: return 'Cache not enabled'
        if len(self.items())==0: return f'Cache at {self.path} is empty\n'
        title = 'Cache contents' if not starts_with else f'Cache entries starting with {starts_with}'
        s = f'{title}\n {"key":30}   {"size":>10}  {"time":20} {"name"}, folder {self.path}\n'
        for name, value in self.items():
            # temporary: override file's path if set
            if str(value)[0]=='/':
                value = self.path/value.name
            if name is None or not name.startswith(starts_with) : continue
            try:
                stat = (self.path/value).stat()
                size = stat.st_size
                mtime= str(datetime.datetime.fromtimestamp(stat.st_mtime))[:16]
                s += f'  {name:30s}  {size:10}  {mtime:20} {value}\n'
            except Exception as msg:
                s += f'{name} -- file not found\n'
                raise
        return s

    def __str__(self):
        return self.show()

# %% ../nbs/00_config.ipynb 9
class Config():
    defaults=\
    """
        verbose         : 1 # set to zero for no output
        warnings        : ignore #

        datapath        : ~/wtlike-data # where to find data--must be set
        cachepath       : ~/.cache/wtlike #

        # Expect 4FGL FITS file, e.g.,  gll_psc_v28.fit
        catalog_file    :

        # multiprocessing
        pool_size       : 1 # number of pool processes to use

        # data cuts, processing
        radius          : 4
        cos_theta_max   : 0.4
        z_max           : 100
        offset_size     : 2.e-06  # scale factor used for event time

        # binning -- actually determined by weight run
        energy_edge_pars : [2,6,17] # pars for np.logspace
        etypes          : [0, 1] # front, back
        nside           : 1024
        nest            : True

        # data selection for cell creation
        week_range      : []  # default all weeks found
        time_bins       : [0, 0, 7] # full MJD range, 7-day cells
        exp_min         : 5    # threshold for exposure per day, in cm^2 Ms units.

        # cell fitting
        use_kerr        : True  # Use the Kerr power-law exposure weighting
        likelihood_rep  : poisson
        poisson_tolerance : 0.2

    """

    def __init__(self, **kwargs):
        import yaml, warnings
        from yaml import SafeLoader

        # parameters: first defaults, then from ~/.config/wtlike/config.yaml, then kwars
        pars = yaml.load(self.defaults, Loader=SafeLoader)
        dp = Path('~/.config/wtlike/config.yaml').expanduser()
        if dp.is_file():
            with open(dp, 'r') as inp:
                userpars = yaml.load(inp, Loader=SafeLoader)
            pars.update(userpars)
            #print(f'update from user file {dp}: {userpars}')
        pars.update(kwargs)

        self.__dict__.update(pars)

        # set warnings filter if requested
        if self.warnings is not None and self.warnings != 'None':
            #print(f'*** set warnings filter to "{self.warnings}" ***')
            warnings.simplefilter(self.warnings)

        self.energy_edges = ee=np.logspace(*self.energy_edge_pars)
        self.energy_bins = np.sqrt(ee[1:] * ee[:-1])
        if not self.week_range:
            self.week_range = (None, None)

       # set up, check files paths
        self.error_msg=''
        if self.datapath is None:
            self.error_msg+='\ndatapath must be set to a folder with wtlike data'

        else:
            self.datapath = df = Path(self.datapath).expanduser()
            if not (self.datapath.is_dir() or  self.datapath.is_symlink()):
                self.error_msg+=f'\ndata_folder "{df}" is not a directory or symlink'
            subs = 'aeff_files weight_files data_files'.split()
            if self.error_msg=='':
                for sub in subs:
                    if not ( (df/sub).is_dir() or  (df/sub).is_symlink()) :
                        self.error_msg+=f'\n{df/sub} is not a directory or symlink'

        self.cachepath =  Path(self.cachepath).expanduser()
        os.makedirs(self.cachepath, exist_ok=True)
        if not self.cachepath.is_dir():
            self.error_msg +=f'cachepath {self.cachepath} is not a folder.'

        # look for 4FGL catalog file, gll_psc_v28.fit currently
        fail = False
        if self.catalog_file is None or self.datapath is None:
            t = Path(self.datapath).expanduser()
            u = sorted(list(t.glob('gll_psc_v*.fit')))
            if len(u)>0:
                self.catalog_file = u[-1]
            else:
                fail = True
        elif Path(self.catalog_file).expanduser().is_file():
            self.catalog_file = Path(config.catalog_file).expanduser()
        else: fail=True

        if fail:
            warnings.warn('There is no link to 4FGL catalog file: set "catalog_file" in your config.yaml'
                  ' or specify if in the Config() call')


    @property
    def cache(self):
        if not hasattr(self, '_cache'):
            self._cache = Cache(self.cachepath, clear=False)
        return self._cache

    @property
    def valid(self):
        if len(self.error_msg)==0: return True
        print(f'wtlike configuration is invalid:\n{self.error_msg}',file=sys.stderr)
        return False

    def __str__(self):
        s = 'Configuration parameters \n'
        for name, value in self.__dict__.items():
            if name=='files' or name.startswith('_'): continue
            s += f'  {name:15s} : {value}\n'
        return s

    def __repr__(self): return str(self)
    def get(self, *pars): return self.__dict__.get(*pars)

# %% ../nbs/00_config.ipynb 13
day = 24*3600.
first_data=54683
#mission_start = Time('2001-01-01T00:00:00', scale='utc').mjd
# From a FT2 file header
# MJDREFI =               51910. / Integer part of MJD corresponding to SC clock S
# MJDREFF =  0.00074287037037037 / Fractional part of MJD corresponding to SC cloc
mission_start = 51910.00074287037
from datetime import datetime
from astropy.time import Time

def MJD(arg):
    """ convert MET or UTC to MJD
    """

    if type(arg)==str:
        if arg=='now':
            return Time(datetime.utcnow()).mjd
        while len(arg.split('-'))<3:
            arg+='-1'
        return Time(arg, format='iso').mjd
    return (mission_start + arg/day  )

def UTC(mjd):
    " convert MJD value to ISO date string"
    t=Time(mjd, format='mjd')
    t.format='iso'; t.out_subfmt='date_hm'
    return t.value

def UTCnow():
    """ current UTC """
    t=datetime.utcnow()
    return f'UTC {t.year}-{t.month:02d}-{t.day} {t.hour:02d}:{t.minute:02d}'

def mission_week(mjd):
    """ return the mission week number for a MJD value
    (Note that week #0 starts on UTC Thursday 2008-05-29 00:00, the first data is in week 9, and that week 525 is missing)
    """
    return (mjd-54615)//(7)

def mjd_range(start,stop, make_round=True):
    """Return a tuple with a valid MJD range
    If either is greater than first_data, interpret as MJD
    otherwise offset from first_data, or now

    So 0,0 is full range, -7, 0 is last 7 days, 0,7 is first 7 days

    """
    a,b = start, stop
    # extract MJD range
    now = MJD('now')
    if a<0: a = now-a
    elif a < first_data : a+=first_data
    if b <= 0 :  b = now-b
    elif b < first_data: b += first_data

    return (round(a),round(b)) if make_round else (a,b)

# %% ../nbs/00_config.ipynb 17
class FermiInterval():
    """For iteration thru (start,stop) tuples in the Fermi data set
    """
    def __init__(self, interval=30, offset=0):
        from wtlike.config import  MJD, first_data
        self.interval=interval
        a,b = first_data+offset, MJD('now')
        self.mm = np.arange(a,round(b),interval)

    def __len__(self):
        return len(self.mm)-1

    def __getitem__(self, k):
        """Return a (start, end) MJD tuple for the kth interval, k=0...N-1. 
        if k<0, return the range for the full interval.
        """
        return (self.mm[k], self.mm[k+1],) if k>=0 else (self.mm[k-1],self.mm[k])

    def __call__(self, k):
        """ Return a timeinterval tuple
        """
        a,b = self[k]
        return (a,b,self.interval)

class FermiMonth(FermiInterval):
    def __init__(self):
        super().__init__(interval=30)

class FermiWeek(FermiInterval):
    def __init__(self):
        super().__init__(interval=7)

# %% ../nbs/00_config.ipynb 19
def bin_size_name(bins):
    """Provide a nice name, e.g., 'day' for a time interval
    """
    if np.isscalar(bins) :
        binsize = bins
    else:
        binsize = np.mean(bins)

    def check_unit(x):
        unit_table = dict(year=1/365.25, week=1/7, day=1, hour=24, min=24*60) #, min=24*60, s=24*3600)
        for name, unit in unit_table.items():
            t = x*unit
            r = np.mod(t+1e-9,1)

            if r<1e-6 : #or t>1:
                return t, name
        return x, 'day'
    n, unit =  check_unit(binsize)
    nt = f'{n:.0f}' if np.mod(n,1)<1e-3 else f'{n:.1f}'
    return f'{nt}-{unit}'# if n>1 else f'{unit}'

# %% ../nbs/00_config.ipynb 21
def decorate_with(other_func):
    def decorator(func):
        func.__doc__ += other_func.__doc__
        return func
    return decorator

# %% ../nbs/00_config.ipynb 22
import time

class Timer():
    """Usage:
    ```
    with Timer() as t:
        time.sleep(5)
    print(t)
    ```
    """

    def __init__(self):
        self.t=time.time()
        self.exit_time=1e6

    def __enter__(self):
        return self
    def __exit__(self, *pars):
        self.exit_time = time.time()-self.t
    def __repr__(self):
         return  f'elapsed time: {self.elapsed:.1f}s ({self.elapsed/60:.1f} min)'
    @property
    def elapsed(self):
        return min(time.time()-self.t, self.exit_time)
