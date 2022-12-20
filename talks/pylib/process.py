import datetime
from wtlike import *
from utilities.ipynb_docgen import *
from utilities.catalogs import  *
from wtlike.poisson import *
from wtlike.loglike import *
from wtlike.sources import SourceLookup

def show(*pars):
    return display_markdown(*pars)

plt.rc('font', size=12)
pd.set_option('display.float_format', lambda f: f'{f:.3f}' if abs(f)>1e-3 else f'{f:.2e}')
SourceLookup.max_sep=0.25

with capture_hide('Catalogs setup') as catalog_setup:
    cat4 = Fermi4FGL()
    uwname = 'uw1216'
    uwcat = UWcat(uwname)

sk = uwcat.skycoord.galactic
uwcat.loc[:,'glat'] = sk.b.deg
uwcat.loc[:,'glon'] = sk.l.deg
uwcat.loc[:, 'r95'] = np.sqrt( (uwcat.a * 2.64)**2 + (0.00791)**2) # 

class SourceFinder(SourceLookup):
    def __init__(self, config=None):
        super().__init__(config or Config())
        self.log = ''
           
    def check(self, name, tol=0.25):    
        self.log = f'Source name, "{name}" '
        pt_name = self(name)

        if pt_name is None: 
            return
        elif self.psep<0.01:
            self.log += f'\n\tcorresponds to {uwname} source "{pt_name}".'
        elif self.psep<tol:
            self.log += f'\n\tis closest to {uwname} source "{pt_name:}", {self.psep:.2f} deg away'
        else:
            self.log += f'\n\tis  {self.psep:.2f} > {tol} deg  away from nearest {uwname} source.'
            return None
        return pt_name

src_finder = SourceFinder()

class Processor():

    def __init__(self, name,  interval=30, nyquist=24):

        if not self.setup(name, interval, nyquist): 
            self.fig = None # flag that no output
            return
 
        # OK: make figure with 4 subplots
        self.fig = plt.figure(figsize=(12,5))
        gs = plt.GridSpec(2, 2,  width_ratios=[4,1], wspace=0.2,hspace=0.5, top=0.95)
        (ax1,ax2,ax3,ax4) =  [self.fig.add_subplot(g) for g in gs]

        self.plot_bb(ax1)        
        self.plot_sed(ax2)
        self.plot_periodogram(ax3, ax4)
        

    def setup(self, name, interval, nyquist):
        pt_name = self.pt_name =  src_finder.check(name)
        if pt_name is None: 
            self.printout = src_finder.log 
            return False
        self.skycoord = src_finder.skycoord

        with capture_hide(f'{name} ({pt_name}) analysis printout') as self.printout:
            print(src_finder.log)
            self.pointlike_info = pi = uwcat.loc[pt_name]
            # pi.loc[:, 'r95'] = pt[loc]
            self.wtl = WtLike(PointSource(pt_name)) 
            self.bb = self.wtl.view(interval).bb_view()
            self.px = self.wtl.periodogram( 1/(4*nyquist) )
        return True

    def catalog_entry(self, cat,  cone_size=0.5):
        """ return the entry from the given catalog that is closest
        If none within cone_size, returns None
        """
        near = cat.select_cone(self.skycoord, cone_size=cone_size)
        if len(near)==0: return None #f'No source within {cone_size} deg'
        ndf = near.sort_values('sep')
        return ndf.iloc[0] 

    def get_catalog_info(self, cat, select=None, show=False, cone_size=0.5):
        """
        select - list of column names
        """

        info = self.catalog_entry(cat, cone_size=cone_size)
        if info is None:
            return f'No {cat.name} source within {cone_size} deg'
        name = info.name
        if select is not None:
            info = info[select]

        sep, r95 = info.sep, info.r95
        color = 'red' if sep>r95 else ''
        with capture(
                f'{cat.name} {name} info <font color={color}>(sep, r95 ={sep:.3f}, {r95:.3f} deg )</font>',
                show=show) as out:
            print(info)
        return out

    def display_4fgl_info(self):
        return self.get_catalog_info(cat4)

    def display_pointlike_info(self):

        return self.get_catalog_info(uwcat,
            select='ra dec glon glat ts eflux100 aprob r95 specfunc sep'.split(),
            )

    def _repr_html_(self):
        return self.printout

    def plot_sed(self, ax):
        # self.wtl.source.sed_plot(ax=ax, ylim=(0.05, 20), title='SED',
        #             yticks = [0.1,1,10], xticks=[0.1,1, 10],
        #         yticklabels='0.1 1 10'.split(), xticklabels='0.1 1 10 '.split());
        c4entry = self.catalog_entry( cat4)
        uwentry = self.catalog_entry( uwcat)
        if c4entry is not None:
            c4entry.specfunc.sed_plot(ax=ax, label=cat4.name)
        uwentry.specfunc.sed_plot(ax=ax, label=uwcat.name,
                    yticks = [0.1,1,10], xticks=[0.1,1, 10],
                    yticklabels='0.1 1 10'.split(), xticklabels='0.1 1 10 '.split()
        )

    def plot_bb(self, ax):

        def g2fit(cell):
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

    def beta_table_mkd(self):
        """ return (hidden) markdown string"""

        with capture_hide(f'beta values for BB intervals (max: {self.beta_table.beta.max():.1f})') as cp:
            print(self.beta_table)
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

    def get_fft_peaks(self, query='p0>25 & f>0.05', show=False):
        df = self.px.find_peaks().query(query)
        if len(df)==0:
            return f'No peaks satisfying {query}'
        with capture(f'{len(df)} FFT peaks satisfying {query}') as out:
            print( df)
        return out

@ipynb_doc
def process_source(name, info=None, nyquist=24):
    """## {name}

    {printout}
    {pinfo}
    {ginfo}
    {fig}
    {beta}
    {fft_peaks}
    {other_info}
    """
    self = Processor(name, nyquist=nyquist)
    printout=self.printout
    if self.fig is not None: 
        fig = self.fig
        pinfo = self.display_pointlike_info()
        ginfo = self.display_4fgl_info()
        beta = self.beta_table_mkd()
        fft_peaks = self.get_fft_peaks()

    else:
        pinfo=ginfo=beta=fig=fft_peaks=''

    if info is not None:
        with capture_hide('other info') as other_info:
            print(info)
    else: other_info=''

    return locals()

