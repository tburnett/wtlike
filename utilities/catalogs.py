""" Access catalogs as dataframes 
"""
import os, sys, glob
import numpy as np
import pylab as plt
import pandas as pd
from pathlib import Path
from astropy.io import fits 
from astropy.table import Table
from astropy.coordinates import SkyCoord, Angle

from . spectral_functions import *

def make_jname(skycoord):
    """Return a name in the format Jhhmm.m+ddmm
    note that last digits are truncated, not rounded
    http://cds.u-strasbg.fr/vizier/Dic/iau-spec.htx#S3.2.1
    """
    import numpy as np
    sc = skycoord.fk5
    ra, dec = sc.ra.deg, sc.dec.deg
    mm = np.mod(ra*4,1440) # RA in minutes 
    ss = np.mod(mm*60,60) #seconds
    HH,MM = int(mm/60), int(mm%60)
    m = int(ss/6) # prescription for .1 min digit
    sign= '+' if dec>=0 else '-'
    dem = int(abs(dec)*60) #abs( DEC) in minutes, truncated
    return 'J' +   '{:02d}{:02d}.{:1d}'.format(HH,MM,m)\
            + sign+'{:02d}{:02.0f}'.format(int(dem/60),dem%60)

def parse_jname(name):
    """ convert a "J-name" to a SkyCoord
    """
    tname = name[:5]+'.0'+name[5:] if name[5]!='.' else name
    ra = (tname[1:3]+'h'+tname[3:7]+'m')
    dec = (tname[7:10]+'d'+tname[10:12]+'m')
    try:
        (ra,dec) = map(lambda a: float(Angle(a, unit='deg').to_string(decimal=True)),(ra,dec))
        return SkyCoord(ra, dec, unit='deg', frame='fk5')
    except ValueError as err:
        print(f'Attempt to parse "{name}" failed ({err}): expect "J1234.5+6789" or "J1234.5678"', file=sys.stderr)
        return None

class MySkyCoord(SkyCoord):
    """Subclass that overrides __repr__ to return "(ra, dec)" 
        Also make pickle-able
    """
    def __repr__(self):
        repr = lambda lon,lat: f'({lon:7.3f}, {lat:+7.3f})'
        ra,dec = self.fk5.ra.deg, self.fk5.dec.deg
        if not hasattr(ra, '__iter__'):
            return repr(ra,dec)
        return np.array([repr(lon,lat) for lon,lat in zip(ra,dec)]).__repr__()
    # For pickle: represent with RA, Dec dict
    def __getstate__(self):
        ra,dec = self.fk5.ra.deg, self.fk5.dec.deg
        return dict(ra=ra,dec=dec)
    def __setstate__(self, state):
        # tricky: make a new object, copy its dict to self
        obj = SkyCoord(state['ra'], state['dec'], unit='deg', frame='fk5')
        self.__dict__.update(obj.__dict__)


class CatDF():
    """Methods common to the following classes
    """

    @property
    def skycoord(self):
        return MySkyCoord(self.ra.values, self.dec.values, unit='deg', frame='fk5')
    
    def match(self, other):
        """other -- another CatDF object
        """
        idx, _delta, _ = self.skycoord.match_to_catalog_sky( other.skycoord )
        return idx, _delta.deg
    
    def hist_delta(self, other, name='', ax=None,  xmin=1e-3,xmax=4):
        _, ax = plt.subplots() if ax is None else (ax.figure, ax)
        _, delta = self.match(other)
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
        - cone_size -- size in degrees about the other
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
            
    def catalog_entry(self, skycoord,  cone_size=0.5):
        """ return the entry  that is closest to the skycoord
        If none within cone_size, returns None
        Add fk5 and galactic LonLat objects to the Series for display convenience
        """
        near = self.select_cone(skycoord, cone_size=cone_size)
        if len(near)==0: return None #f'No source within {cone_size} deg'
        ndf = near.sort_values('sep')
        info =  ndf.iloc[0].copy() 
        info.loc['fk5'] = LonLat(info.ra, info.dec)
        if 'glon' in info:
            info.loc['galactic'] = LonLat(info.glon, info.glat)
        else: #klugy
            info.loc['galactic'] = None
        
        return info
   
class LATpsr(CatDF, pd.DataFrame):
    
    """Return the LAT pulsar list as a DataFrame
    index is the Source_Name field
    """

    def __init__(self, *args, path='$FERMI/catalog/srcid/cat/obj-pulsar-lat_v1*', **kwargs):
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
    
    def __init__(self, model='uw1410', filename=None):

        if model is not None and filename is None:
            filename = Path(os.path.expandvars('$FERMI'))/f'skymodels/sources_{model}.csv'
        elif filename is not None:
            filename = Path(filename)
        else:
            pass
        
        assert filename.exists(), f'File {filename} not found'
        uwdf = pd.read_csv(filename, index_col=0)

        sf = []
        for n, v in uwdf.iterrows():    
            sf.append( specfun_dict[v.modelname](np.array(v.pars[1:-1].split(), float), e0=v.e0))
        uwdf.loc[:,'specfunc'] = sf

        super().__init__(uwdf)
        print(f'Loaded UW model {model}: {len(self)} entries')
        self.__dict__.update(name=model)

class LonLat():
    def __init__(self, lon,lat):
        self.lon, self.lat=lon,lat
    @property
    def as_tuple(self):
        return (self.lon, self.lat)
    def __repr__(self):
        return f'({self.lon:7.3f},{self.lat:+7.3f})'


class Fermi4FGL(CatDF, pd.DataFrame):

    PRIMARY_HEADER_KEYS = ('TELESCOP', 'INSTRUME', 'ORIGIN', 'DATE', 'VERSION')
    CATALOG_IDENTIFYING_KEYS = ('VERSION', 'CDS-NAME', 'EXTNAME', 'HDUCLAS1', 'HDUCLAS2', 'HDUVERS')

    class FlagBits():
       
        """From  Table 4 in the 4FGL DR3 paper https://arxiv.org/abs/2201.11184
            1 TS < 25 with other model or analysis
            2 Moved beyond 95% error ellipse
            3 Flux changed with other model or analysis
            4 Source/background ratio < 10%
            5 Confused
            6 Interstellar gas clump (c sources)
            9 Localization flag from pointlike
            10 Bad spectral fit quality
            12 Highly curved spectrum
            13 TS < 25 at 12 yr
            14 Soft Galactic Unassociated (§ 6.2)
            """
        def __init__(self, f):
            self.f = f
        def __getitem__(self, n):
            return self.f & 2**(n-1)>0
        def __repr__(self):
            r = ''
            for n in range(1, 16):
                if (self.f & 2**(n-1)) >0:
                    r+= f'{n},'
            return '{}' if r=='' else '{'+r[:-1]+'}'


    def __init__(self, version=None, *, path='$FERMI/catalog/', reset_index=False):
        """
        - version: "vnn" for release, or DR3 or DR4 for 4FGL data release
                    "uwnnnn" for an internal version

        - path: either a folder containing the catalog FITS files for use by the version lookup, 
                or the absolute path to an appropriate FITS file             
        """
 
        t=Path(os.path.expandvars(path)).expanduser(); 
        if version is not None:
            if not version.startswith('uw'):
                version = {'DR4':'v32', 'DR3':'v28', 'FL16Y':'v36'}.get(version, version)
        if t.is_file():
            filename = t
        elif t.is_dir():
            if version.startswith('uw'):
                u = list(t.glob(f'gll_pscP305{version}*.fits')); 
                assert len(u)>0, f'No files for version {version} found.'
                filename=u[-1]

            else:
                filename = sorted(list(t.glob('gll_psc_v*.fit')))[-1] if version is None else \
                    Path(f'{t}/gll_psc_{version}.fit')
                
        else:
            raise Exception( f'path {path} is not a directory or a FITS file')
        if not Path(filename).is_file():
            raise Exception(f'File {filename} does not exist.')
        print(f'Loaded Fermi 4FGL {filename.name}', end='')
        resolved_filename = Path(filename).expanduser().resolve()
        with fits.open(filename) as hdus:
            data = hdus[1].data
            header0 = hdus[0].header
            header1 = hdus[1].header

        primary_header_values = {k: header0.get(k) for k in self.PRIMARY_HEADER_KEYS}
        catalog_identifying_values = {k: header1.get(k) for k in self.CATALOG_IDENTIFYING_KEYS}

        cname= lambda n : [s.strip() for s in data[n]]
        cvar = lambda a: data[a].astype(np.float32)
        ivar = lambda a: data[a].astype(int)
        name = list(map(lambda x: x.strip() , data['Source_Name']))
        energy_factor = 1 if data.columns['Energy_Flux100'].unit.startswith('erg') else 1.602178e-06

        # calculate these first
        funcs = self.specfuncs(data)


        cat_subset =  dict(
            ra          = cvar('RAJ2000'),
            dec         = cvar('DEJ2000'), 
            glat        = cvar('GLAT'),
            glon        = cvar('GLON'),
            r95         = cvar('Conf_95_SemiMajor'),
            specfunc    = funcs,
            pivot       = cvar('Pivot_Energy'),
            eflux       = cvar('Energy_Flux100') * energy_factor, # to erg
            significance= cvar('Signif_Avg'),
            curvature   = 2*cvar('LP_Beta'),
            curv_unc    = 2*cvar('Unc_LP_Beta'),

            
            # ....
        )
        if 'Variability_Index' in data.columns.names:
            cat_subset['variability'] = cvar('Variability_Index')

        release_format =  'ASSOC1' in data.columns.names
        if release_format:
            # release format
            cat_subset.update(dict(
                assoc_prob  = cvar('ASSOC_PROB_BAY'), # for Bayesian, or _LR for likelihood ratio
                assoc1_name = cname('ASSOC1'),
                class1      = cname('CLASS1'),
                flags       = list(map(self.FlagBits, ivar('FLAGS'))),
            ))
        elif 'Passoc' in data.columns.names:
            # internal association format
            cat_subset.update(dict(
                assoc_prob  = cvar('Passoc'), 
                assoc1_name = cname('assoc_new'),
                class1      = cname('class_new'),
                nickname    = cname('NickName'),
                ts          = cvar('Test_Statistic'),
            ))
        else:
            # internal no assocation
            cat_subset.update(dict(
                nickname    = cname('NickName'),
                ts          = cvar('Test_Statistic'),
            ))

        super().__init__( cat_subset)

        print( f': {len(self)} entries' )
        # extract actual version
        v = header1['VERSION'] 
        v = {'v28':'DR3', 'v32':'DR4'}.get(v,v)
        vname = version if version=='FL16Y' else f'4FGL-{v}'
        index_name = str(header1.get('CDS-NAME', vname)).strip()
        if len(index_name) == 0:
            index_name = vname
        self.__dict__.update(
            data=data,
            filename=filename.name,
            catalog_file=str(resolved_filename),
            name=vname,
            primary_header_values=primary_header_values,
            catalog_identifying_values=catalog_identifying_values,
        )

        #index by source names except now with FL16Y
        if not reset_index:
            # index with Source_Name or NickName
            if release_format:
                self.index = name
            else:
                self.index = self.nickname
                self['jname'] = name 
            self.index.name = index_name
        else:
            # ordinal indexing
            self['name'] = name
        self.fitscols = data.columns

    def info(self):
        """Print catalog file path and selected FITS header identification keys."""
        print(f'Catalog file: {self.catalog_file}')

        print('Primary header keys:')
        for key in self.PRIMARY_HEADER_KEYS:
            print(f'  {key}: {self.primary_header_values.get(key)}')

        print('\nCatalog-identifying keys:')
        for key in self.CATALOG_IDENTIFYING_KEYS:
            print(f'  {key}: {self.catalog_identifying_values.get(key)}')
       
    def specfuncs(self, data):
        """ Return a list of spectral functions
        """
        # special version that
        specfun_dict = dict(
            PowerLaw=PowerLaw,
            LogParabola=LogParabola,
            PLSuperExpCutoff=PLSuperExpCutoff4,
            PLSuperExpCutoff4=PLSuperExpCutoff4,

            )
        cvar = lambda a: data[a].astype(float)
        cname= lambda n : [s.strip() for s in data[n]]
        def par_array(colnames):
            # return a transposed table of the columns 
            return np.array(list(map(cvar, colnames))).T

        pardict = dict(
            LogParabola=par_array('LP_Flux_Density LP_Index LP_beta Pivot_Energy'.split()),
            PowerLaw    =par_array('PL_Flux_Density PL_Index'.split()),
            PLSuperExpCutoff=par_array('PLEC_Flux_Density PLEC_IndexS PLEC_ExpfactorS PLEC_Exp_Index'.split()),
            PLSuperExpCutoff4=par_array('PLEC_Flux_Density PLEC_IndexS PLEC_ExpfactorS PLEC_Exp_Index'.split())
                    )

        pivot = cvar('Pivot_Energy')
        spec = []                    
        for i,name in enumerate(cname('SpectrumType')):
            pars = pardict[name][i]
            spec.append( specfun_dict[name]( pars, e0=pivot[i]))
        return spec

    def field(self, name):
        return self.data.field(name)

    def get_specfunc(self, src_name, func_name='default', ):
        """ Return a spectral function 

        - src_name - The 4FGL name
        - func_name - one of PLEC4, LP, PL, default
        """

        from utilities.spectral_functions import PLSuperExpCutoff4, LogParabola, PowerLaw 
        specs = dict(
            PLEC4=(PLSuperExpCutoff4,'PLEC_Flux_Density PLEC_IndexS PLEC_ExpfactorS PLEC_Exp_Index'.split()),
            LP = (LogParabola, 'LP_Flux_Density LP_Index LP_beta Pivot_Energy'.split(), ),
            PL = (PowerLaw, 'PL_Flux_Density PL_Index'.split(),), 
        )
        fgl_names =list(self.index) 
        row = self.data[fgl_names.index(src_name)]
        default = dict(PLSuperExpCutoff4='PLEC4', PLSuperExpCutoff='PLEC4', #allow for mistake?
                       LogParabola='LP', PowerLaw='PL')[row['SpectrumType']]
        func, parnames = specs.get(func_name if func_name !='default'  else default)    
        pars = [row[name] for name in parnames] 
        return func(pars, e0=row['Pivot_Energy']) if pars is not None else None  

    def get_series(self, name):
        """If name is a source:
                Return the full record for the source as a Pandas Series object
            else
                Return the column 
        """
        if name not in self.index:
            # it is a field name
            return pd.Series(self.field(name), index=self.index)
        index = list(self.index).index(name)
        assert index>=0, f'Failed to find {name}'
        s = self.data[index]
        cols = self.data.columns
        return pd.Series(dict( (c, s.field(c)) for c in cols.names),  name=name) 
    
    def flux_band(self, name):

        s = self.get_series(name)
        flux = s['Flux_Band']
        unc = s['Unc_Flux_Band']
        return pd.Series(dict(
                flux=flux[1:], 
                low=unc[1:,0], 
                high=unc[1:,1],
                islimit=pd.isna(unc[1:,0]),
                nufnu=s['nuFnu_Band'][1:],
                ),  name=name)

    def band_plot(self, name, ax=None, emax=None, specfun=True, ms=0, **kwargs):
        
        fb = self.flux_band(name)
        e_bins = np.array([ 0.1, 0.3, 1, 3, 10, 30, 100, 1000])
        energy = np.sqrt(e_bins[:-1]* e_bins[1:])
        e_err = np.array([energy-e_bins[:-1],e_bins[1:]-energy])
        
        lim = fb.islimit
        bar =~lim
        err = np.array([-fb.low, fb.high, ])
        fig,ax=plt.subplots(figsize=(8,5)) if ax is None else (ax.figure, ax)
        scale = 1e9 * energy # / 624.150648
        ax.errorbar(x=energy[bar], xerr=e_err[:,bar], y=(scale*fb.flux)[bar], 
                    yerr=(scale*err)[:,bar], fmt='o', ms=ms);
        if sum(lim)>0:
            _,caps,_ = ax.errorbar(x=energy[lim], xerr=e_err[:,lim], y=(y:=(scale* fb.high)[lim]), 
                        uplims=True*sum(lim), yerr=y/2,  fmt='.',color='grey')
            caps[0].set_markersize(ms)
        if specfun:
            sf = self.get_specfunc(name, 'LP')
            sf.sed_plot(ax=ax,plot_kw=dict(color='orange'),
                    )#label='4FGL-DR3')        
        # ax.scatter(energy, 1e12*fb.nufnu, marker='o', facecolor='none', edgecolor='w', s=400 )
        kw = dict(xscale='log', xlim=(0.1, 100), xlabel='Energy (GeV)', 
            yscale='log', ylim=(0.02), yticks=[0.1,1,10, 100], 
                yticklabels='0.1 1 10 100'.split()
                )
        kw.update(kwargs)
        ax.set(**kw)
        # ax.set_title(fb.name, fontsize=10) 
        # ax.text(0.5,0.92, fb.name, fontsize=12, ha='center', transform=ax.transAxes) # ylabel=r'$\nu F_{\nu}\ \ [\mathrm{10^{-12}\ erg\ cm^{-2}\ s^{-1}}]$',
        return fig   
