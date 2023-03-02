"""
process.py: Analysis of the time behaviour of a source

To run an individual source: examine_source(name)
where name can be:
* any source identifier recognized by SIMBAD
* the name of a uw source in uw1410, specifically in the wtlike table of weights
* a "j-name": format J1234+5678 or J1234.5-6789
"""
import datetime
from wtlike import *
from utilities.ipynb_docgen import ( show, DataFrameWrapper, ipynb_doc, capture_hide, capture, monospace)
from utilities.catalogs import*

from wtlike.sources import SourceFinder
__all__ = ['show', 'uwcat', 'cat4', 'get_proc', 'Summarize', 'SourceAnalyzer', 'load_source_spreadsheet',
    'setup_excel', 'catalog_setup','get_fermi_info', 'examine_source','src_finder', 
    'process_excel', 'process_df']

uwname = 'uw1410'
pars = sys.argv
if len(pars)>1:
    if pars[1].startswith('uw'):
        uwname = pars[1]
    else:
        print( f'Unrecognized parameter(s): {pars[1:]}')

# def show(*pars):
#     return display_markdown(*pars)

plt.rc('font', size=12)

def float_format(f):
    if f==0: return '0'
    return f'{f:.3f}' if (abs(f)>1e-3 and abs(f)<1e5 ) else f'{f:.2e}'

pd.set_option('display.float_format', float_format)

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

defaults = dict(neighbor=None, interval=30, nyquist=24, max_sep=0.5, tsmin=16, info_name='Other info',
            fft_query='p1>25 & f>0.05')

class Summarize():
    """ base class for call-back class"""

    def add(self, sa): 
        pass

class SourceAnalyzer():

    def __init__(self, name, 
            neighbor=None,  
            interval=30, 
            nyquist=24, 
            tsmin=25,
            make_figs=True,
            summarizer=Summarize(),
            **kwargs):
        """
        """
        self.log=''
        self.wtl=None
        self.max_sep = kwargs.pop('max_sep', defaults['max_sep'])
        self.porb = float(kwargs.pop('porb', np.nan))
        self.fluxlim = kwargs.pop('fluxlim', None)
        self.kwargs = kwargs 
        self.fig = None
        assert isinstance(summarizer, Summarize), f'Expect {summarizer} to be subclass of Summarize'
        
        SourceFinder.max_sep = self.max_sep #klugy
        if not self.setup(name, neighbor, interval, nyquist, tsmin): 
            return
        if make_figs:
            self.make_plots()
        summarizer.add(self)

    def plot_lc(self, ax1=None,ax2=None):
        if self.wtl is  None: return

        if ax1 is None:
            # OK: make figure with 2 subplots
            fig = plt.figure(figsize=(12,2.5))
            gs = plt.GridSpec(1, 2,  width_ratios=[4,1], wspace=0.2,hspace=0.5, top=0.95)
            (ax1,ax2) =  [fig.add_subplot(g) for g in gs]
        else: fig=None
        self.plot_bb(ax1)        
        self.plot_seds(ax2, self.max_sep)
        return fig

    def make_plots(self):

        if self.wtl is not None:
            self.fig=plt.figure(figsize=(12,5))
            gs = plt.GridSpec(2, 2,  width_ratios=[4,1], wspace=0.2,hspace=0.5, top=0.95)
            (ax1,ax2,ax3,ax4) =  [self.fig.add_subplot(g) for g in gs]
            self.plot_lc(ax1=ax1, ax2=ax2)
            self.plot_FFT(ax1=ax3, ax2=ax4)
        else:
            # no analysis since low TS. Only make the sed plots
            self.fig, ax = plt.subplots(figsize=(3,3))
            self.plot_seds(ax, self.max_sep)
        
    def __repr__(self):
        return f'{self.__class__.__name__}("{self.name}")'

    def setup(self, name, neighbor, interval, nyquist, tsmin, **kwargs):
        tstart = self.kwargs.get('tstart',0)
        tstop = self.kwargs.get('tstop',0)
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

        with capture_hide(f'Analysis printout for {name} ({pt_name})') as self.printout:
            print(src_finder.log)
            if pt_name not in uwcat.index:
                print(f'UW source {pt_name} not found in {uwname}??')
                self.wtl = self.bb=self.px = None
                self.pointlike_info = FileNotFoundError
                return False
            self.pointlike_info = pi = uwcat.loc[pt_name]
            ts = pi.ts
            if ts < tsmin:
                self.log = f'UW source has TS={ts:.1f} < {tsmin}'
                self.wtl = None
                return True

            self.wtl = WtLike(PointSource(pt_name), time_bins=(tstart, tstop, interval)) 

            if neighbor is not None and str(neighbor).strip():
                print(f'Reweighting with {neighbor}')
                self.neighbor_bb = WtLike(neighbor).bb_view()
                self.wtl = self.wtl.reweighted_view(self.neighbor_bb)
            self.bb = self.wtl.view(interval).bb_view()
            
            self.px = self.wtl.periodogram( 1/(4*nyquist), tstart=tstart, tstop=tstop)

            if not(np.isnan(self.porb)) and self.porb is not None:
                # generate orbial phase view to plot later (no periastron yet)
                bins = self.kwargs.get('phase_bins', 8)
                # create a phase view
                self.orbital_phase = self.wtl.phase_view(period=self.porb, nbins=bins)
                phase_df= self.orbital_phase.fits.query('t<1')
                self.orbital_ts = -2 * np.sum([f(1) for f in phase_df.fit.values])
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
            select='jname fk5 galactic specfunc e0 eflux100 ts fitqual locqual aprob r95  sep'.split(),
            )

    def _repr_html_(self):
        if self.fig is None:
            self.fig = make_plots()
            return FigureWrapper(self.fig)._repr_html_()
        return self.printout

    def plot_seds(self, ax=None, cone_size=0.5, **kwargs):
        _, ax = plt.subplots(figsize=(4,4)) if ax is None else (ax.figure, ax)
        c4entry = cat4.catalog_entry(self.skycoord, cone_size=cone_size)
        uwentry = uwcat.catalog_entry(self.skycoord, cone_size=cone_size)

        if uwentry is not None:
            uwentry.specfunc.sed_plot(ax=ax, e0=uwentry.e0, label=uwcat.name)
        if c4entry is not None:
            c4entry.specfunc.sed_plot(ax=ax, e0=c4entry.pivot, label=cat4.name)
            
        kw = dict( yticks=[0.1,1,10], xticks=[0.1,1,10],
            yticklabels='0.1 1 10'.split(), xticklabels='0.1 1 10'.split(),
                ylabel=r'$E^2\ dN/dE\ \mathrm{(eV\ cm^{-2}\ s^{-1})}$')
        kw.update(kwargs)
        ax.set(**kw)
        
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
            if self.fluxlim is not None:
                ax.set(ylim=self.fluxlim) ###
            self.bb.plot(ax=ax, source_name=self.wtl.source.name, UTC=True)

            df_bb = self.bb.fluxes['t tw ts flux errors'.split()]
            df_beta = pd.DataFrame.from_dict( 
                dict((i, g2fit(cell)) for i,cell in self.bb.cells.iterrows()) ,orient='index')\
                    ['flux beta sig_beta'.split()]
            self.bb_table = pd.concat([df_bb, df_beta], axis=1)

            
        except Exception as e:
            print(f'oops: bb failed: {e}')
            self.wtl.plot(ax=ax)
            self.bb_table = None
            # self.beta_table = None

    def display_bb_table(self):
        """ return (hidden) markdown string"""
        bt = self.bb_table
        if bt is None or len(bt)==0: return ''
        test = np.max(bt.beta/bt.sig_beta)
        
        warning_text = f'<font color=red>Check beta:  max beta/sig_beta is {test:.1f}</font>' if test>2 else ''
        # with capture_hide(f'BB fits {warning_text}') as cp:
        #     print(bt)
        # return cp
        dw = DataFrameWrapper(bt, summary=f'{len(bt)} BB fits {warning_text}', index=False)
        return dw._repr_html_()

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

    def display_fft_peaks(self, show=False):
        query = defaults['fft_query']
        df = self.px.find_peaks('p1').query(query)
        if len(df)==0:
            return monospace(f'No FFT peaks satisfying {query}')
        df['period'] = 1/df.f
        # with capture(f'{len(df)} FFT peaks satisfying {query}: max(p1)={max(df.p1):.1f}') as out:
        #     with pd.option_context('display.precision', 6, 'display.float_format',None):
        #         print( df)
        with pd.option_context('display.precision', 6, 'display.float_format',None):
            dw = DataFrameWrapper(df,
                    summary=f'{len(df)} FFT peaks satisfying {query}: max(p1)={max(df.p1):.1f}',
                    index=False, show=show)
            return dw._repr_html_()

    def get_nearby(self, radius_cut=3, var_cut=30):
        near_df = cat4.select_cone(self.skycoord, cone_size=5, query=f'variability>{var_cut}& sep<{radius_cut}')
        near_df['name'] = near_df.index
        return near_df['sep glon glat significance variability'.split()].sort_values('sep')

    def display_nearby(self, show=False, radius_cut=3, **kwargs):
        near_df = self.get_nearby(radius_cut=radius_cut, **kwargs)
        if len(near_df)==0: return f' No nearby (within {radius_cut} deg) variable sources.'
        dw = DataFrameWrapper(near_df,
            summary=f'{len(near_df) } nearby variable sources'
            )
        return dw._repr_html_()
        # with capture(f'{len(near_df) } nearby variable sources', show=show) as out:
        #     print(near_df)
        # return out

    def plot_FFT(self, ax1=None, ax2=None, peak_query='f>0.01'):
        """ Generate a Kerr Periodogram and a histogram of the peak values
        
        px -- TimeSeries object
        
        Return the figure with two Axes plots if had to create
        """
        def hist_peak_power(px, pname='p1', ax=None, xlim=(0,40), thresh=25,  **kwargs):
            """
            """
            p1 = px.find_peaks(pname).query(peak_query)[pname].clip(*xlim)
            _, ax = plt.subplots(figsize=(3,2)) if ax is None else (ax.figure, ax)
            kw = dict(xlabel=f'Power {pname}', xlim=xlim, title='' )
            kw.update(kwargs); ax.set(**kw)
            ax.hist(p1, np.linspace(*xlim,), log=True, histtype='step', lw=2);
            ax.hist(p1[p1>thresh], np.linspace(*xlim,), log=True, histtype='stepfilled', color='red');
            ax.grid(alpha=0.5)
        
        if ax1 is None:
            fig = plt.figure(figsize=(12,2.5))
            gs = plt.GridSpec(1, 2,  width_ratios=[4,1], wspace=0.2,hspace=0.5, top=0.95)
            (ax1,ax2) =  [fig.add_subplot(g) for g in gs]
        else: fig=None

        px = self.px
        px.power_plot(ax=ax1, pmax=50)
        hist_peak_power(px, ax=ax2,  xlabel='')
        ax2.set( xlabel=r'$P_1$ peak values')
        return fig

    def plot_orbital_phase(self, ax=None, **kwargs):
        if self.porb is None:
            print('No orbital period available', file=sys.stderr)
            return
       
        fig,ax = plt.subplots(figsize=(5,3)) if ax is None else (ax.figure, ax)
        kw = dict(xlim=(0,1), xticks=np.arange(0,1.01,0.25)); 
        kw.update(kwargs)
        self.orbital_phase.plot(ax=ax, **kw)
        ax.axvline( 1.0, ls='--', color='lightgrey')
        return fig

    def plot_phase(self, ax=None, **kwargs):
        """A phase plot with perhaps different porb, bins, or ref
        """
        period = kwargs.pop('period', self.porb)
        nbins  = kwargs.pop('nbins', 10)
        ref    = kwargs.pop('ref', '2008')

        fig,ax = plt.subplots(figsize=(6,3)) if ax is None else (ax.figure, ax)
        kw = dict(xlim=(0,2),ylim=None); kw.update(kwargs)
        
        self.wtl.phase_view(period=period, nbins=nbins, reference=ref).plot(ax=ax, **kw)
        
        ax.axvline( 1.0, ls='--', color='lightgrey')
        return fig
        
    def plot_phase(self, ax=None, **kwargs):

        period = kwargs.pop('period', self.porb)
        nbins  = kwargs.pop('nbins', 10)
        ref    = kwargs.pop('ref', '2008')        

        fig,ax = plt.subplots(figsize=(6,3)) if ax is None else (ax.figure, ax)
        kw = dict(xlim=(0,1),ylim=(0,None),); kw.update(kwargs)
        
        self.wtl.phase_view(period=period, nbins=nbins, reference=ref).plot(ax=ax, **kw)
        
        ax.axvline( 1.0, ls='--', color='lightgrey')
        return fig

# access to the SourceAnalyzer object from the notebook
proc=None
def get_proc(): return proc

@ipynb_doc
def examine_source(name, info=None, text='',  **kwargs): 

    """### {name}
    {other_info}
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
    {phase_plot}
    
    """
    global proc
    ginfo=beta=fig=fft_peaks=nearby=neighbor_plot=pinfo=other_info=phase_plot=''
    kw = defaults.copy() 
    kw.update(kwargs)
    if info is not None: kw.update(info)
    max_sep = kw.pop('max_sep', None)
    info_name = kw.pop('info_name', 'Source info')
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
        beta = self.display_bb_table()
        fft_peaks = self.display_fft_peaks()

    if neighbor is None:
        neighbor_plot = ''
    elif self.wtl is not None: 
        neighbor_plot, ax = plt.subplots(figsize=(10,3))
        self.neighbor_bb.plot(ax=ax)
        neighbor_plot.summary= f'Reweighted with {neighbor} (click for its light curve)'

    if getattr(self, 'orbital_phase', None) is not None:
        phase_plot = self.plot_orbital_phase()
        phase_plot.summary=f'Phase plot (period {self.porb:.6f} d, TS={self.orbital_ts:.0f} )'

    return locals()

def process_df(df, **kwargs): #, max_sep=None, tsmin=50):

    defaults.update(kwargs)
    with capture('Default parameter values') as parout:
        print(pd.Series(defaults))

    show(f"""\
        ## Processing {len(df)} sources, from {str(datetime.datetime.now())[:16]}
        {parout}
        ---
        """)
    for name, info in df.iterrows():
        if type(name) !=str:
            print(f'Bad entry: {name}', file=sys.stderr)
            continue
        examine_source(name, info=info, **kwargs ) 

    show(f"""\
        \n---
        \n## Finish at {str(datetime.datetime.now())[:16]}""")

class WTSkyCoord(SkyCoord):
    def __repr__(self):
        ra,dec = self.fk5.ra.deg, self.fk5.dec.deg
        return f'({ra:7.3f},{dec:+7.3f})'


def load_source_spreadsheet(
            filename='AMXP scorecard (2).xlsx',
            source_name='AMXP name',
            ra_name=None, dec_name=None,
            flag_row=True):
    """Return a dataframe derived from the table in a spreadsheet
    row 1 is column names
    if flag_row is True row 2 has non-blank to include the column in the datafrane
    if ra_name (and dec_name) are not None, rename columns to "ra" and "dec"
    The index is set to the contents of source_nane
    """
    spreadsheet = Path(filename) #'AMXP scorecard (2).xlsx')
    assert spreadsheet.is_file(), f'File "{filename}" not found'
    df = pd.read_excel(spreadsheet)
    assert source_name in df.columns, f'did not find column "{source_name}" to use as index'
    # take care of trailing non-break space
    df.index = list(map(lambda s:str(s).strip().replace(u"\u2014","-"),df.loc[:,source_name]))
    if flag_row:
        to_drop = pd.isna(df.iloc[0])
        to_drop[source_name]=True
        df = df.drop(columns=df.columns[to_drop])[1:]
    else:
        df = df.drop(columns=[source_name])
    if ra_name is not None:
        df.rename(columns={ra_name: 'ra', dec_name: 'dec'}, inplace=True)
    return df

def get_fermi_info(source_names, max_sep=0.5):
    """
    Return a dataframe with information on UW and 4FGL sources near the locations of a list of source names
    """
    
    src_finder = SourceFinder()
    SourceFinder.max_sep = max_sep
    
    class SourceInfo(dict):

        class SkyCoord(SkyCoord):
            # subclass that displays (ra,dec)
            def __repr__(self):
                ra,dec = self.fk5.ra.deg, self.fk5.dec.deg
                return f'({ra:6.3f},{dec:+6.3f})'

        def __init__(self, name):
            #self.name=name.replace('_',' ')
            src_finder.log=''

            # get the nearest pointlike name from the weight tables
            self.pt_name = pt_name = src_finder(name) 

            try:
                self['skycoord'] = self.SkyCoord(src_finder.skycoord)
            except ValueError as err:
                print(f'SourceInfo: Fail to recognize name "{name}"', file=sys.stderr)
                return

            if pt_name is None: return

            # grab pointlike info
            uw = self.get_uwinfo(pt_name)
            if len(uw)==0: 
                self.update(uw_name = pt_name) 
                return
            self.update(uw_name=pt_name,
                        uw_pos=self.SkyCoord(uw.ra,uw.dec,unit='deg',frame='fk5'), 
                        uw_sep = uw.sep,
                        uw_r95 = uw.r95,
                        uw_ts = uw.ts,
                       )
            dr = self.get_4fgl()
            if len(dr)==0: return
            self.update(dr_name=dr.name,
                        dr_pos=self.SkyCoord(dr.ra,dr.dec,unit='deg', frame='fk5'),
                        dr_sep = dr.sep,
                        dr_r95 = dr.r95,
                        dr_class1=dr.class1)
        def get_uwinfo(self, pt_name):
            #return a dict 
            if pt_name not in uwcat.index:
                # print(f'{self.name}: UW source {pt_name} not found in {uwname}??', file=sys.stderr)
                return pd.Series(name=pt_name, dtype=object)
            info = uwcat.catalog_entry( self['skycoord'], cone_size=SourceFinder.max_sep)
            return pd.Series(info)

        def get_4fgl(self):
            info = cat4.catalog_entry( self['skycoord'], cone_size=SourceFinder.max_sep)
            if info is None: return pd.Series(dtype=object)
            return pd.Series(info)

        def __repr__(self):
            return f'{self.__class__.__name__}({self.name}):\n{pd.Series(self, dtype=object)}'
    
    source_names = np.atleast_1d(source_names)
    return pd.DataFrame([SourceInfo(name) for name in source_names], index=source_names )

def setup_excel(filename, 
        source_name_column,
        ra_name=None, dec_name=None,
        title='Excel analysis',
        query='',
        select = (None,),
        **kwargs):

    defaults.update(kwargs)

    spreadsheet = Path(filename) #'AMXP scorecard (2).xlsx')
    assert spreadsheet.is_file(), f'File "{filename}" not found'
    if len(title)>1 and title[0] != '#' : title = '# '+title
    show(f"""\
        {title}

        {catalog_setup}
        Read spreadsheet "{spreadsheet.absolute()}", 
        dated {str(datetime.datetime.fromtimestamp(spreadsheet.stat().st_mtime))[:16]}
        """)
    # make the index the source name, remove columns without a "1" (or anything?) in the second row
    df = pd.read_excel(spreadsheet)
    assert source_name_column in df.columns, f'did not find column "{source_name_column}" to use as index'
    df.index = df.loc[:,source_name_column]
    to_drop = pd.isna(df.iloc[0])
    if type(select) != tuple: 
        # assume selecting a single row
        select = (select,select+1)
    df = df.drop(columns=df.columns[to_drop])[1:][slice(*select)]
    if ra_name is not None and dec_name is not None:
        sc = SkyCoord(df[ra_name], df[dec_name], unit='deg', frame='fk5')
        df.loc[:,'skycoord'] = CatDF.SkyCoord(sc)
    if query:        
        df = df.query(query)
    return df

def process_excel(filename, source_name_column,  **kwargs):
    df = setup_excel(filename, source_name_column,  **kwargs)
    return process_df(df, **kwargs)