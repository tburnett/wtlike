# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/00_config.ipynb (unless otherwise specified).

__all__ = ['Cache', 'Files', 'Config', 'MJD', 'UTC', 'day', 'first_data', 'PointSource']

# Cell
from astropy.time import Time
from astropy.coordinates import SkyCoord
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple
import os, sys
import numpy as np
import pickle

# Cell
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
        if not self.path: return
        self.path.mkdir(exist_ok=True)
        assert self.path.exists()
        self.index_file = self.path/'index.pkl'

        if self.path.exists():
            if clear:
                print('Clearing cache!')
                self.clear()
            else:
                self._load_index()
        else:
            self.path.mkdir(exist_ok=True)

    def _dump_index(self):
        with open(self.index_file, 'wb') as file:
            pickle.dump(self, file)

    def _load_index(self):
        if not self.index_file.exists():
            self._dump_index()
            return
        with open(self.index_file, 'rb') as file:
            self.update(pickle.load(file))

    def add(self, key, object, exist_ok=False):
        if not self.path: return
        if key  in self:
            if not exist_ok:
                print(f'Warning: cached object for key "{key}" exists', file=sys.stderr)
            filename = self[key]
        else:
            filename = self.path/f'cache_file_{hex(key.__hash__())[3:]}.pkl'
            self[key] = filename
            self._dump_index()

        with open(filename, 'wb') as file:
            pickle.dump(object, file )


    def get(self, key):
        if key not in self:
            return None
        filename = self[key]
        if not filename.exists():
            # perhaps deleted by another instance?
            print(f'File for Cache key {key} not found, removing entry', file='sys.stderr')
            selt.pop(key)
            return None
        with open(filename, 'rb') as file:
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
        filename = self[key]
        filename.unlink()
        super().pop(key)
        self._dump_index()


    def __call__(self, key, func, *pars, description='', overwrite=False, **kwargs,
                ):
        """
        One-line usage interface for cache use

        - `key` -- key to use, usually a string. Must be hashable <br>
            If none, ignore cache and return the function evaluation
        - `func` -- user function that will return an object that can be pickled
        - `pars`, `kwargs` -- pass to `func`
        - `description` -- optional string that will be printed

        Example:
        <pre>
        mycache = Cache('/tmp/thecache', clear=True)

        def myfun(x):
            return x

        result = mycache('mykey', myfun, x=99,  description='My data')

        </pre>

        """

        if key is None:
            return func(*pars, **kwargs)

        ret = self.get(key)
        if description:
            print(f'{description}: {"Saving to" if key not in self else "Restoring from"} cache with key "{key}"')

        if ret is None or overwrite:
            ret = func(*pars, **kwargs)
            self.add(key, ret)
        return ret

    def __str__(self):
        import datetime
        if not self.path: return 'Cache not enabled'
        s = f'Cache contents\n {"key":20}   {"size":>10}  {"time":20} {"name"}, in folder {self.path}\n'
        for name, value in self.items():
            if name is None: continue
            stat = value.stat()
            size = stat.st_size
            mtime= str(datetime.datetime.fromtimestamp(stat.st_mtime))[:16]
            s += f'  {name:20s}  {size:10}  {mtime:20} {value.name}\n'
        return s

# Cell
@dataclass
class Files:
    """ paths to the various files that we need"""

    data:str = '$HOME/data'
    ft2: str = '$HOME/work/lat-data/ft2'
    gti: str = '$HOME/work/lat-data/binned/'
    aeff:str = '$HOME/work/lat-data/aeff'
    weights: str = '$HOME/onedrive/fermi/weight_files'

    # expand -- not implemented in Path
    def __post_init__(self):
        d = self.__dict__
        for name, value in d.items():
            d[name] = Path(os.path.expandvars(value))

    @property
    def valid(self):
        """assume all files ok if aeff"""
        return self.aeff.exists()

    def __repr__(self):
        s = 'File paths for light curves\n'
        for name, value in self.__dict__.items():
            s += f'  {name:10s} : {value}\n'
        return s

# Cell
@dataclass
class Config:
    """Default light curve configuration parameters"""
    verbose : int = 3
    files :'' =  None

    # cache
    cachepath: str = '/tmp/cache'

    # photon selection
    mjd_range : Tuple=None
    radius: float = 5
    cos_theta_max:float=0.4
    z_max : float=100

    # binning
    energy_edges = np.logspace(2,6,17)
    time_interval: int = 1
    use_uint8: bool=False  # for weights

    # healpix data representation used by data
    nside : int=1024
    nest: bool=True

    # exposure calculation
    bins_per_decade: int=4
    base_spectrum: str='lambda E: (E/1000)**-2.1'
    energy_range: Tuple = (100.,1e6)

    # analysis
    likelihood_rep: str='poisson'

    def __post_init__(self):
        if self.files is None: self.files=Files()

    @property
    def cache(self):
        if not hasattr(self, '_cache'):
            self._cache = Cache(self.cachepath, clear=False)
        return self._cache

    @property
    def valid(self):
        return self.files.valid

    def __str__(self):
        s = 'Configuration parameters \n'
        for name, value in self.__dict__.items():
            if name=='files' or name.startswith('_'): continue
            s += f'  {name:15s} : {value}\n'
        return s

    def __repr__(self): return str(self)

# Cell

day = 24*3600.
first_data=54683

def MJD(met):
    "convert MET to MJD"
    #mission_start = Time('2001-01-01T00:00:00', scale='utc').mjd
    # From a FT2 file header
    # MJDREFI =               51910. / Integer part of MJD corresponding to SC clock S
    # MJDREFF =  0.00074287037037037 / Fractional part of MJD corresponding to SC cloc
    mission_start = 51910.00074287037
    return (mission_start + met/day  )

def UTC(mjd):
    " convert MJD value to ISO date string"
    t=Time(mjd, format='mjd')
    t.format='iso'; t.out_subfmt='date_hm'
    return t.value

# Cell
class PointSource():
    """Manage the position and name of a point source
    """
    def __init__(self, name, position=None):
        """position: (l,b) tuple or None. if None, expect to be found by lookup
        """
        self.name=name
        if position is None:
            skycoord = SkyCoord.from_name(name)
            gal = skycoord.galactic
            self.l,self.b = (gal.l.value, gal.b.value)
        else:
            self.l,self.b = position
            skycoord = SkyCoord(self.l,self.b, unit='deg', frame='galactic')
        self.skycoord = skycoord
    def __str__(self):
        return f'Source "{self.name}" at: (l,b)=({self.l:.3f},{self.b:.3f})'
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
        return self.name.replace(' ', '_').replace('+','p')

    @classmethod
    def fk5(cls, name, position):
        """position: (ra,dec) tuple """
        ra,dec = position
        sk = SkyCoord(ra, dec, unit='deg',  frame='fk5').transform_to('galactic')
        return cls(name, (sk.l.value, sk.b.value))