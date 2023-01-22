"""
process.py: Analysis of the time behaviour of a source

Do run an individual source: examine_source(name)
where name can be:
* any source identifier recognized by SIMBAD
* the name of a uw source in uw1216, specifically in the wtlike table of weights
"""
import datetime, inspect
from wtlike import *
from utilities.ipynb_docgen import *
from utilities.catalogs import  *

from wtlike.sources import SourceFinder

uwname = 'uw1400'
pars = sys.argv
if len(pars)>1:
    if pars[1].startswith('uw'):
        uwname = pars[1]
    else:
        print( f'Unrecognized parameter(s): {pars[1:]}')

def show(*pars):
    return display_markdown(*pars)

plt.rc('font', size=12)
pd.set_option('display.float_format', lambda f: f'{f:.3f}' if (abs(f)>1e-3 and abs(f)<1e5 ) else f'{f:.2e}')

with capture_hide(f'Catalogs: {uwname} and 4FGL-DR3') as catalog_setup:
    if os.environ.get('FERMI', None) is None:
        os.environ['FERMI'] ='.'
        print(f'Setting env var FERMI to ".". Expect to find folders catatag and skymodels')
    cat4 = Fermi4FGL()
    uwcat = UWcat(uwname)

sk = uwcat.skycoord.galactic
uwcat.loc[:,'glat'] = sk.b.deg
uwcat.loc[:,'glon'] = sk.l.deg
uwcat.loc[:, 'r95'] = np.sqrt( (uwcat.a * 2.64)**2 + (0.00791)**2) # 

src_finder = SourceFinder()

defaults = dict(neighbor=None, interval=30, nyquist=24, max_sep=0.5, tsmin=25, info_name='Other info')

class SourceAnalyzer():

    def __init__(self, name, 
            neighbor=None,  
            interval=30, 
            nyquist=24, 
            tsmin=50,
            **kwargs):
        """
        """
        self.log=''
        self.wtl=None
        self.max_sep = kwargs.pop('max_sep', defaults['max_sep'])
        SourceFinder.max_sep = self.max_sep #klugy
        if not self.setup(name, neighbor, interval, nyquist, tsmin): 
            self.fig = None # flag that no output
            return

        if self.wtl is not None:
        # OK: make figure with 4 subplots
            self.fig = plt.figure(figsize=(12,5))
            gs = plt.GridSpec(2, 2,  width_ratios=[4,1], wspace=0.2,hspace=0.5, top=0.95)
            (ax1,ax2,ax3,ax4) =  [self.fig.add_subplot(g) for g in gs]

            self.plot_bb(ax1)        
            self.plot_sed(ax2, self.max_sep)
            self.plot_periodogram(ax3, ax4)
        else:
            # no analysis since low TS. Make the plots anyway
            self.fig, ax = plt.subplots(figsize=(3,3))
            self.plot_sed(ax, self.max_sep)
        
    def __repr__(self):
        return f'{self.__class__.__name__}("{self.name}")'

    def setup(self, name, neighbor, interval, nyquist, tsmin):
        pt_name = self.pt_name =  src_finder.find(name, tol=self.max_sep)
        if pt_name is None: 
            self.log =  f'<font color="red">Failed lookup:</font> <br>{monospace(src_finder.log)} ' 
            #raise Exception(self.log)
            self.printout=''
            self.name='(not found)'
            self.skycoord = None
            return False
        self.skycoord = src_finder.skycoord
        self.name = name

        with capture_hide(f'{name} ({pt_name}) analysis printout') as self.printout:
            print(src_finder.log)
            if pt_name not in uwcat.index:
                print(f'UW source {pt_name} not found in {uwname}??')
                self.wtl = self.bb=self.px = None
                self.pointlike_info = FileNotFoundError
                return False
            ts = uwcat.loc[pt_name,'ts'] 

            self.pointlike_info = pi = uwcat.loc[pt_name]
            if ts < tsmin:
                self.log = f'UW source has TS={ts:.1f} < {tsmin}'
                self.wtl = None
                return True
            # pi.loc[:, 'r95'] = pt[loc]
            self.wtl = WtLike(PointSource(pt_name)) 
            if str(neighbor).strip() =='': neighbor
            if neighbor is not None and str(neighbor).strip():
                print(f'Reweighting with {neighbor}')
                self.neighbor_bb = WtLike(neighbor).bb_view()
                self.wtl = self.wtl.reweighted_view(self.neighbor_bb)
            self.bb = self.wtl.view(interval).bb_view()
            self.px = self.wtl.periodogram( 1/(4*nyquist) )
        return True

    def get_catalog_info(self, cat, select=None, show=False):
        """
        select - list of column names
        """
        cone_size = SourceFinder.max_sep
        info = cat.catalog_entry( self.skycoord, cone_size=cone_size)
        if info is None:
            return f'No {cat.name} source within {cone_size} deg'
        name = info.name
        if select is not None:
            info = info[select]

        sep, r95 = info.sep, info.r95
        color = 'red' if sep>r95 else ''
        with capture(
                f'{cat.name} ({name}) <font color={color}>(sep, r95 ={sep*60:.1f}, {r95*60:.1f} arcmin )</font>',
                show=show) as out:
            print(info)
        return out

    def display_4fgl_info(self):
        try:
            return self.get_catalog_info(cat4,
            'fk5 galactic specfunc pivot eflux significance flags variability assoc_prob class1 assoc1_name  r95 sep'.split() )
        except:
            return 'No 4FGL-DR3 info'

    def display_pointlike_info(self):
        if self.pt_name is None:
            return ''
        if self.pt_name not in uwcat.index:
            return f'<font color="red"> Source {self.pt_name} not in {uwname}</font>'
        
        return self.get_catalog_info(uwcat,
            select='fk5 galactic specfunc e0 eflux100 ts fitqual locqual aprob r95  sep'.split(),
            )

    def _repr_html_(self):
        return self.printout

    def plot_sed(self, ax, cone_size=1):

        c4entry = cat4.catalog_entry(self.skycoord, cone_size=cone_size)
        uwentry = uwcat.catalog_entry(self.skycoord, cone_size=cone_size)
        if c4entry is not None:
            c4entry.specfunc.sed_plot(ax=ax, e0=c4entry.pivot, label=cat4.name)
        if uwentry is not None:
            uwentry.specfunc.sed_plot(ax=ax, e0=uwentry.e0, label=uwcat.name,
                    yticks = [0.1,1,10], xticks=[0.1,1, 10],
                    yticklabels='0.1 1 10'.split(), xticklabels='0.1 1 10 '.split()
        )

    def plot_bb(self, ax):

        def g2fit(cell):
            from wtlike.poisson import Poisson
            from wtlike.loglike import LogLike, Gaussian2dRep
            
            ts = Poisson.from_function(LogLike(cell)).ts
            r = dict(t=cell.t, tw=cell.tw, ts=round(ts,1))            
            if ts<4:
                r.update(flux=0, counts=cell.n, beta=np.nan, sig_beta=np.nan)
            else:
                r.update(Gaussian2dRep(LogLike(cell)).fit)
            return r
        try:
            self.bb.plot(ax=ax, source_name=self.wtl.source.name, UTC=True)
            self.bb_table = self.bb.fluxes['t tw flux errors'.split()]
            self.beta_table = pd.DataFrame.from_dict( 
                dict((i, g2fit(cell)) for i,cell in self.bb.cells.iterrows()) ,orient='index')\
                    ['t flux beta sig_beta'.split()]
        except Exception as e:
            print(f'oops: bb failed: {e}')
            self.wtl.plot(ax=ax)
            self.bb_table = None
            self.beta_table = None

    def display_beta_table(self):
        """ return (hidden) markdown string"""
        bt = self.beta_table
        if bt is None or len(bt)<2: return ''
        test = np.max(bt.beta/bt.sig_beta)
        color = 'red' if test>2 else ''
        with capture_hide(f'beta values for BB intervals (<font color={color}>max sig = {test:.1f}</font> )') as cp:
            print(bt)
        return cp

    def plot_periodogram(self, ax3, ax4):
        
        def hist_peak_power(px, pname='p1', ax=None, xlim=(0,40), thresh=25,  **kwargs):
            """
            """
            p1 = px.find_peaks(pname).query('f>0.1')[pname].clip(*xlim)
            fig, ax = plt.subplots(figsize=(3,2)) if ax is None else (ax.figure, ax)
            kw = dict(xlabel=f'Power {pname}', xlim=xlim, title='' )
            kw.update(kwargs); ax.set(**kw)
            ax.hist(p1, np.linspace(*xlim,), log=True, histtype='step', lw=2);
            ax.hist(p1[p1>thresh], np.linspace(*xlim,), log=True, histtype='stepfilled', color='red');
            ax.grid(alpha=0.5)
            return fig
        # self.px = self.wtl.periodogram( 1/(4*nyquist) )
        self.px.power_plot(ax=ax3, pmax=50)
        hist_peak_power(self.px, ax=ax4, title='Peak distribution', xlabel='')

    def display_fft_peaks(self, query='p1>30 & f>0.05', show=False):
        df = self.px.find_peaks('p1').query(query)
        if len(df)==0:
            return monospace(f'No peaks satisfying {query}')
        with capture(f'{len(df)} FFT peaks satisfying {query}') as out:
            print( df)
        return out

    def get_nearby(self, radius_cut=3, var_cut=30):
        near_df = cat4.select_cone(self.skycoord, cone_size=5, query=f'variability>{var_cut}& sep<{radius_cut}')
        near_df['name'] = near_df.index
        return near_df['sep glon glat significance variability'.split()].sort_values('sep')

    def display_nearby(self, show=False, radius_cut=3, **kwargs):
        near_df = self.get_nearby(radius_cut=radius_cut, **kwargs)
        if len(near_df)==0: return f' No nearby (within {radius_cut} deg) variable sources.'
        with capture(f'{len(near_df) } nearby variable sources', show=show) as out:
            print(near_df)
        return out

# access to the SourceAnalyzer object from the notebook
proc=None
def get_proc(): return proc

@ipynb_doc
def examine_source(name, info=None, text='',  **kwargs): 

    """## {name}
    <br>
    {text}
    {log}
    {printout}
    {neighbor_plot}
    {pinfo}
    {ginfo}
    {fig}
    {nearby}
    {beta}
    {fft_peaks}
    {other_info}
    """
    global proc
    ginfo=beta=fig=fft_peaks=nearby=neighbor_plot=pinfo=other_info=''
    kw = defaults.copy() 
    kw.update(kwargs)
    max_sep = kw.pop('max_sep', None)
    info_name = kw.pop('info_name', 'Other info')
    if max_sep is not None: SourceFinder.max_sep=max_sep

    # the additional info, a dict-like object, can update default analysis parameters
    if info is not None:
        for k,v in info.items():
            if k in defaults and str(v)!='nan':
                # print(f'setting {k} from {kw[k]} to {v}')
                kw[k]=v

        with capture_hide(f'{info_name}') as other_info:
            print(info)
    else: other_info=''


    neighbor = kw.get('neighbor',None)
    if not str(neighbor).strip(): neighbor=None
    neighbor_info = '' if neighbor is None else f'Reweighted with {neighbor}'
    
    try:
        proc = self = SourceAnalyzer(name, **kw)
    except Exception as ex:
        printout = f'<font color="red">Failed: {ex}</font>\n'
        raise
        return locals()
    log = '' if not self.log else f'<font color="red">{self.log}</font>'
    # proc = self = SourceAnalyzer(name, **kw)
    printout=self.printout
    if self.skycoord is None:
        return locals()
    pinfo = self.display_pointlike_info()
    ginfo = self.display_4fgl_info()
    fig = self.fig
    
    nearby = self.display_nearby()

    if self.wtl is not None:
        beta = self.display_beta_table()
        fft_peaks = self.display_fft_peaks()

    if neighbor is None:
        neighbor_plot = ''
    elif self.wtl is not None: 
        neighbor_plot, ax = plt.subplots(figsize=(10,3))
        self.neighbor_bb.plot(ax=ax)
        neighbor_plot.summary= f'Reweighted with {neighbor} (click for its light curve)'

    return locals()

def process_df(df): #, max_sep=None, tsmin=50):
    with capture('Default parameter values') as parout:
        print(pd.Series(defaults))

    show(f"""\
        \nProcessing {len(df)} sources,  Start at {str(datetime.datetime.now())[:16]}
        \n{parout}
        \n---
        """)
    for name, info in df.iterrows():
        if type(name) !=str:
            print('Bad entry', file=sys.stderr)
            continue
        examine_source(name, info=info, )  
    show(f"""\
        \n---
        \n# Finish at {str(datetime.datetime.now())[:16]}""")


def process_excel(filename, 
        source_name_column,
        title='Excel analysis',
        query='',
        **kwargs):

    defaults.update(kwargs)

    spreadsheet = Path(filename) #'AMXP scorecard (2).xlsx')
    assert spreadsheet.is_file(), f'File "{filename}" not found'
    
    show(f"""\
        # {title}

        {catalog_setup}
        Read spreadsheet "{spreadsheet.absolute()}", 
        dated {str(datetime.datetime.fromtimestamp(spreadsheet.stat().st_mtime))[:16]}
        """)
    # make the index the source name, remove columns without a "1" (or anything?) in the first row
    df = pd.read_excel(spreadsheet)
    assert source_name_column in df.columns, f'did not find column "{source_name_column}" to use as index'
    df.index = df.loc[:,source_name_column]
    to_drop = pd.isna(df.iloc[0])
    df = df.drop(columns=df.columns[to_drop])[1:]
    if query:        
        df = df.query(query)

    with capture_hide(f'Spread sheet: {query}') as ss:
        print(df)
    show(f'{ss}')
    
    process_df(df)