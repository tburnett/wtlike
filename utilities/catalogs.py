""" Access catalogs as dataframes 
"""
import os, sys, glob
import numpy as np
import pylab as plt
import pandas as pd
from pathlib import Path
from astropy.io import fits 
from astropy.table import Table
from astropy.coordinates import SkyCoord


class CatDF():
    """Methods common to the following classes
    """
    @property
    def skycoord(self):
        return SkyCoord(self.ra.values, self.dec.values, unit='deg', frame='fk5')
    
    def match(self, other):
        """other -- another CatDF object
        """
        idx, _delta, _ = self.skycoord.match_to_catalog_sky( other.skycoord )
        return idx, _delta.deg
    
    def hist_delta(self, other, name='', ax=None,  xmin=1e-3,xmax=4):
        fig, ax = plt.subplots() if ax is None else (ax.figure, ax)
        idx, delta = self.match(other)
        ax.hist(delta.clip(xmin, xmax), np.logspace(np.log10(xmin), np.log10(xmax), 40), 
        histtype='stepfilled', log=True);
        ax.set(xscale='log', ylim=(0.5,None), title=name);
        ax.grid(alpha=0.5)
        
    def add_match(self, other, prefix):
        idx, delta = self.match(other)
        self.loc[:,f'{prefix}_name'] = other.index[idx]
        self.loc[:,f'{prefix}_delta'] = delta
    
    def select_cone(self, other, cone_size=5, query=''):
        """ Return a DataFrame of this catalog's entries within the cone_size of the given source
        
        - other -- name | SkyCoord
        - cone_size -- size in degrees abut the other
        - query -- an optional query string to apply to the result
        
        """
        if type(other)==str:
            try:
                skycoord = SkyCoord.from_name(other)
            except Exception as err:
                print({err}, file=sys.stderr)
                return
        elif isinstance(other, SkyCoord):
            skycoord = other
        else:
            raise Exception(f'Expected {other} to be either the name of a source, or a SkyCoord')
            
        sep = self.skycoord.separation(skycoord).deg; 
        cone = sep < cone_size
        incone = self.loc[cone, :].copy()
        incone['sep'] = sep[cone]
        return incone if query=='' else incone.query(query)

    def find_nearest(self, other, cone_size=0.1):
        try:
            return self.select_cone(other, cone_size=cone_size).iloc[0].name
        except Exception as msg:
            print(f'Failed to find a source near "{other}" : {msg}', file=sys.stderr)
            

   
class LATpsr(CatDF, pd.DataFrame):
    
    """Return the LAT pulsar list as a DataFrame
    index is the Source_Name field
    """

    def __init__(self, *args, path='$FERMI/catalog/obj-pulsar-lat_v1*', **kwargs):
        from astropy.table import Table
        lcat = glob.glob(os.path.expandvars(path))[-1]
        filename= os.path.split(lcat)[-1]
        version = filename.split('_')[-1][:-5]

        print (f'Loaded LAT pulsars {version}', end='')
        super().__init__(Table.read(lcat, hdu=1).to_pandas(), *args, **kwargs)
        self.version = filename.split('.')[0][-4:]                 
        self['msec']=msec = np.array([code.decode().find('m')>-1 for code in self.CHAR_Code ], bool) 
        print (': {} entries, {} millisecond pulsars'.format(len(self), sum(msec)) )
        self.index = map(lambda name: name.strip().decode(), self.Source_Name.values)
        self.rename(columns=dict(RAJ2000='ra', DEJ2000='dec'), inplace=True)
        self.index.name = 'name'
        del self['Source_Name']
        
        self.__dict__.update(version=version, )

        
class BigFile(CatDF, pd.DataFrame):
    
    """Return BigFile as a DataFrame.
    Index is the PSRJ name
    """

    def __init__(self, path='$FERMI/catalog/Pulsars_BigFile_*.fits',  **kwargs):
        ff = sorted(glob.glob(os.path.expandvars(path)))
        filename = ff[-1]
        version =  filename.split('_')[-1][:-5] 
        print(f'Loaded BigFile {version}', end='')
        with fits.open(filename) as t:
            super().__init__( t[1].data, **kwargs)
        print(f': {len(self)} entries')
        self.version =version        
        names=[t.strip() for t in self.NAME.values]
        jnames=[t.strip() for t in self.PSRJ.values]
        psrnames = map(lambda s:'PSR '+s, jnames)
        self.rename(columns=dict(RAJD='ra', DECJD='dec'), inplace=True)
        self.index = psrnames
        self.index.name='name'

    
class UWcat(CatDF, pd.DataFrame):
    
    def __init__(self, model='uw1216'):
        
        filename = Path(os.path.expandvars('$FERMI'))/f'skymodels/sources_{model}.csv'
        assert filename.exists(), f'File {filename} not found'
        super().__init__(pd.read_csv(filename, index_col=0))
        print(f'Loaded UW model {model}: {len(self)} entries')


class Fermi4FGL(CatDF, pd.DataFrame):

    def __init__(self, path='$FERMI/catalog/'):
 
        t=Path(os.path.expandvars(path)).expanduser(); 
        assert t.is_dir(), f'path {path} is not a directory'
        filename = sorted(list(t.glob('gll_psc_v*.fit')))[-1]
        print(f'Loaded Fermi 4FGL-DR3 {filename.name}', end='')
        with fits.open(filename) as hdus:
            data = hdus[1].data

        cname= lambda n : [s.strip() for s in data[n]]
        cvar = lambda a: data[a].astype(float)
        ivar = lambda a: data[a].astype(int)
        name = list(map(lambda x: x.strip() , data['Source_Name']))

        super().__init__( dict(
            ra          = cvar('RAJ2000'),
            dec         = cvar('DEJ2000'), 
            glat        = cvar('GLAT'),
            glon        = cvar('GLON'),
            r95         = cvar('Conf_95_SemiMajor'),
            spectype    = cname('SpectrumType'),
            eflux       = cvar('Energy_Flux100'), # erg cm-2 s-1
            significance= cvar('Signif_Avg'),
            variability = cvar('Variability_Index'),
            assoc_prob  = cvar('ASSOC_PROB_BAY'), # for Bayesian, or _LR for likelihood ratio
            assoc1_name = cname('ASSOC1'),
            assoc2_name = cname('ASSOC2'),
            class1      = cname('CLASS1'),
            class2      = cname('CLASS2'),
            flags       = ivar('FLAGS'),
            # ....
        ))
        print( f': {len(self)} entries' )
        self.__dict__.update(data=data, filename=filename.name,)
        self.index = name
        self.index.name = 'name'
        self.fitscols = data.columns
        
    def field(self, name):
        return self.data.field(name)
