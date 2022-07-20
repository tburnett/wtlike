from wtlike import *
from utilities.catalogs import Fermi4FGL
from pathlib import Path
from utilities.ipynb_docgen import *
import sys

def interpret(lcr):
    def extract_value(x):
        # expect "< num" or "num" -- set negative if limit
        y = x.split()
        return -float(y[-1]) if y[0]=='<'  else float(y[0])
    def extract_error(x):
        # nan if limit
        return np.nan if x=='-'  else float(x)

    flux_vals = lcr.iloc[:,4].apply(extract_value).values
    flux_errs = lcr.iloc[:,5].apply(extract_error).values

    gamma = lcr.iloc[:,6].apply(extract_error).values
    t = MJD(lcr.iloc[:,2])
    return pd.DataFrame(dict(t=t,
                             flux=flux_vals,
                             flux_unc=flux_errs,
                             ts=lcr.TS, 
                             Gamma=gamma))

# def old_interpret(lcr):
    
#     def extract_value(x):
#         y = x.split()
#         if x[0]=='<': return -float(y[1])
#         elif x[0]=='-': return np.nan
#         return float(y[0])

#     def extract_error(x):
#         y = x.split()
#         return float(y[2]) if x[0]!='<' else np.nan
#     flux_vals = lcr.iloc[:,4].apply( extract_value)
#     flux_errs = lcr.iloc[:,4].apply( extract_error)
#     gamma = lcr.iloc[:,5].apply(extract_value).values
#     t = MJD(lcr.iloc[:,2])
#     return pd.DataFrame(dict(t=t,flux=flux_vals,flux_unc=flux_errs,
#                                        ts=lcr.TS, Gamma=gamma))


class LCR():
    """
    Load an LCR light curve and compare with wtlike
    """

    base_url='https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/source.html?source_name=' 
    

    
    def __init__(self, srcname, csv_filename=None, lcr_dir='.', config=None, key='', pct=80):
        self.name=srcname

        # look up the 4FGL name
        fermi_name = srcname if srcname.startswith('4FGL ') else Fermi4FGL().find_nearest(srcname)

        # make up the filename
        if csv_filename is None:
            csv_filename = list(Path(lcr_dir).glob(f'{fermi_name.replace(" ","_")}*.csv'))[-1]
        else:
            csv_filename = Path(csv_filename)

        assert csv_filename.is_file(), f'{csv_filename} is not a file.' 

        lcr = pd.read_csv(csv_filename)
        print(f'Read LCR light curve file {csv_filename} with {len(lcr)} entries to {lcr.iloc[-1,0]}.')

        # make a dataframe, and set index to rounded time, accounting for 0.000754 d = 65 s offset
        self.lcr_df = interpret(lcr) #if len(lcr.columns)==11  else old_interpret(lcr)
        self.lcr_df.index, self.lcr_df.index.name = self.lcr_df.t.round(1), 'tindex'

        # get a mean corresponding to 
        # self.lcr_mean = lcr_mean = self.lcr_df.query('t<59765').flux.clip(0,1e-4).mean()
        # print(f'\tMean LCR flux: {lcr_mean:.2e}')

        # Wtlike weekly light curve
        self.wtl = WtLike(srcname, config=config or Config(verbose=1), key=key)
        self.wtl_df = self.wtl.fluxes.copy()
        self.wtl_df.loc[:,'flux_unc'] = self.wtl_df.errors.apply(lambda x: 0.5*(x[1]-x[0]))
        
        self.df=df = self.wtl.fluxes.copy()
        df.index, df.index.name = df.t, 'tindex' # will use time for index

        # combine with some LCR values -- depends on index for each being consistent
        # note adjust flux to be mean, for comparison
        df.loc[:,'lcr_flux'] = self.lcr_df.flux
        df.loc[:,'lcr_unc'] = self.lcr_df.flux_unc
        df.loc[:,'unc'] = df.errors.apply(lambda x: 0.5*(x[1]-x[0]))
        df.loc[:,'lcr_ts'] = self.lcr_df.ts
        # df.loc[:,'lcr_Gamma'] = self.lcr_df.Gamma

        # now adjust the LCR fits to the same mean as the top 20% wtlike
        fx_cut = np.percentile(self.wtl_df.flux, pct)
        subset = df[df.flux>fx_cut]
        lcrfx = subset.lcr_flux
        wtlfx = subset.flux
        r = lcrfx/wtlfx
        norm = self.lcr_mean = r.mean()
        print(f'Set normalization for LCR fits: {norm:.2e}')
        df.lcr_flux /= norm
        df.lcr_unc /= norm
        
        self.norm = 1. #norm = np.mean(ratio) 
        self.cdf=cdf = df[(df.lcr_flux>0) & (df.ts>9) & ~ np.isnan(df.lcr_unc) ]
        
        # stats
        print(f'Compare {len(self.lcr_df)} week fits for TS>4 measurements')
        df_both=df.query('ts<4 & lcr_flux<0 ')
        print(f'\t{"Both missing":15} {len(df_both):5}')

        df_lcr_only = df.query('ts<4 & lcr_flux>0 '); 
        print(f'\t{"LCR, no wtlike":15} {len(df_lcr_only):5}')

        df_wtl_only = df.query('ts>4 & lcr_flux<0'); 
        print(f'\t{"wtlike, no LCR":15} {len(df_wtl_only):5}')


    def plot(self, fignum=1, figsize=(12,6), **kwargs):
        fig, (ax1,ax2) = plt.subplots(2,1, figsize=figsize, sharex=True, num=fignum)
        plt.subplots_adjust( hspace=0, top=0.92 )
        
        df = self.lcr_df[ (self.lcr_df.flux>0) & ~ np.isnan(self.lcr_df.flux_unc) ]
        ax1.errorbar(x=df.t, y=df.flux/self.lcr_mean, 
                yerr=df.flux_unc/self.lcr_mean,fmt='.' , label=f'rel to {self.lcr_mean:.2e}')
        #ax1.plot(df.t, df.flux, '.');
        dfx = self.lcr_df[self.lcr_df.flux<0]
        ax1.plot(dfx.t, -dfx.flux/self.lcr_mean, 'v', ms=5,color='orange')
        
        ax1.grid(alpha=0.5)
        ax1.set(xlabel='MJD', ylabel='flux', **kwargs)
        ax1.text(0.5,0.85, 'LCR',ha='center', color='blue',
                 size='large', transform=ax1.transAxes)
        ax1.legend()

        self.wtl.plot(ax=ax2,ts_bar_min=9, source_name='', **kwargs);
        fig.suptitle(f'{self.name}: Weekly light curves by wtlike and LCR', fontsize=16);
        ax2.text(0.5,0.85, 'wtlike', ha='center', color='BLUE',
                 size='large', transform=ax2.transAxes)
        return fig
        
    def flux_comparison(self, ax=None, fignum=2, fmax=15, lim=None, **kwargs):

        cdf =self.cdf
        norm = self.norm
        fig, ax = plt.subplots( figsize=(10,10), num=fignum) if ax is None else (ax.figure, ax)
        
        lim = (0,fmax) if lim is None else lim
        y,x = cdf.flux.clip(*lim), (cdf.lcr_flux/norm).clip(*lim)
        yerr, xerr,= cdf.unc, cdf.lcr_unc/norm
        ax.errorbar(x,y, xerr=xerr, yerr=yerr,  fmt='.', color='grey');
        #ax.plot(cdf.flux,  cdf.lcr_flux/norm, '.');
        ax.plot(lim,lim, '--', color='orange')
        ax.grid(alpha=0.5)
        kw = dict(ylabel='wtlike relative flux', xlim=lim,ylim=lim, xlabel='LCR normalized flux', 
              title='Flux comparison')
        kw.update(kwargs)
        ax.set(**kw);
        return fig

    def unc_comparison(self, ax=None, fignum=3, umax=0.5, **kwargs):
        cdf =self.cdf
        lim2 = (0,umax)
        y = (cdf.unc/cdf.flux).clip(*lim2)
        x = (cdf.lcr_unc/cdf.lcr_flux).clip(*lim2)

        fig, ax = plt.subplots( figsize=(5,5), num=fignum) if ax is None else (ax.figure, ax)
        ax.plot(x, y, '.');
        ax.plot(lim2, lim2, '--', color='orange')
        ax.grid(alpha=0.6)
        kw = dict(ylabel='wtlike relative flux unc', xlim=lim2, ylim=lim2,
                xlabel='LCR relative flux unc', aspect=1,
          title = 'Uncertainty comparison')
        kw.update(kwargs)
        ax.set(**kw)
        fig.set_facecolor('white')
        return fig
        
    def ratio_comparison(self, fignum=4):
        fig, (ax, ax2) = plt.subplots(2,1,figsize=(12,6), sharex=True, num=fignum)
        plt.subplots_adjust(hspace=0, top=0.92)
        cdf = self.cdf
        norm = self.norm
        ax.plot(cdf.flux,  cdf.lcr_flux/norm/cdf.flux, '.')
        ax.grid(alpha=0.5)
        ax.axhline(1, color='grey')
        ax.set(ylim=(0,3), xscale='log', xlabel='wtlike relative flux',   ylabel='Flux',)    
        ax2.plot(cdf.flux, (cdf.lcr_unc/norm/cdf.unc).clip(0.5,2), '.');
        ax2.axhline(1, color='grey'); ax2.grid(alpha=0.5)
        ax2.set(xscale='log', ylim=(0.5,2) ,xlabel='wtlike relative flux',ylabel='Uncertainty');
        fig.suptitle('LCR/wtlike flux, flux uncertainty ratios')
        fig.set_facecolor('white')
        return fig
        
    def lcr_index(self, fignum=4):
        fig, ax =plt.subplots(figsize=(8,4), num=fignum)
        cdf = self.cdf
        ax.semilogx( cdf.lcr_flux.clip(1e-8,1e-5), cdf.lcr_Gamma.clip(1,4), '.');
        ax.set(ylabel='Photon index', xlabel='Flux', title='LCR weekly fits to 3C 273')
        ax.axhline(2.6, color='orange', label='12-year index')
        ax.grid(alpha=0.5); ax.legend();
        fig.set_facecolor('white')
        return fig

    def ts_comparison(self, ax=None, lim=(4,4000), fignum=5):

        df = self.df
        y,x = df.ts.clip(*lim), df.lcr_ts.clip(*lim)

        lmin=lim[0]
        lcr_only=y==lim[0]
        wtlike_only= x==lim[0]
        both = ~lcr_only & ~wtlike_only

        fig, ax = plt.subplots(figsize=(6,6), num=fignum) if ax is None else (ax.figure, ax)
        ax.loglog(x,y, '.', label=f'both fit [{sum(both)}]');
        ax.plot(x[lcr_only], y[lcr_only], 'o', color='violet', label=f'LCR only [{sum(lcr_only)}]')
        ax.plot(x[wtlike_only], y[wtlike_only], 'o', color='green', label=f'wtlike only [{sum(wtlike_only)}]')

        ax.plot(lim,lim)
        ax.set(ylabel='wtlike', xlabel='LCR', xlim=lim, ylim=lim,title='TS comparison');
        ax.grid(); ax.legend(fontsize=12);
        return fig

    def compare_all(self, fmax=10, umax=0.5, title=None):
    
        fig, axx = plt.subplots(ncols=3, figsize=(18,5))
        plt.subplots_adjust(wspace=0.25, top=0.86)
        self.ts_comparison(ax=axx[0])
        self.flux_comparison(ax=axx[1], fmax=fmax)
        self.unc_comparison(ax=axx[2], umax=umax);
        if title: fig.suptitle(title)
        return fig

    def flux_profile(self, ax=None, 
                    lim=(0,1.4), 
                    bins=(0.675,1.375, 14), **kwargs):

        from utilities.uplots import make_profile
        from numpy.polynomial import Polynomial
        
        df = self.df # the compaison
        # screen out bad points
        tdf = self.df[~np.isnan(df.lcr_ts)]

        pf = make_profile(tdf.flux,tdf.lcr_flux, np.linspace(*bins))

        fitfun = Polynomial.fit(pf.y, pf.x, deg=1)

        fig, ax = plt.subplots(figsize=(6,6)) if ax is None else (ax.figure, ax)
        ax.errorbar(pf.y, pf.x, pf.xerr, pf.yerr, ls='none', lw=2, marker='o',  )

        kw = dict( ylabel='wtlike', xlabel='LCR',  xlim=lim, ylim=lim, aspect=1,
                    title='')
        kw.update(kwargs)
        ax.set(**kw)    
        ax.grid(alpha=0.5)
        plotdom = np.linspace(*lim, num=25)
        ax.plot( plotdom, plotdom, '--', color='grey', label='expected');
        ax.plot( plotdom, fitfun(plotdom),'-', color='orange', label='observed');
        ax.legend();
        return fig
    
@ipynb_doc
def show_profile(name, **kwargs):
    """### Flux comparison profile plot for {name}
    {fig}
    """
    with capture_hide('setup') as out:
        lcr =LCR(name,lcr_dir='/home/burnett/work/lcr_data')
    kw = dict(lim=(0.5,1.5))
    kw.update(kwargs)
    fig = lcr.flux_profile( **kw);
    fig.width=400
    return locals()
    
@ipynb_doc
def lcr_check(name, fmax=10, umax=0.4):
    """
    #### {name}: wtlike vs. LCR
    {fig1}
    {fig2}
    """
    with capture_hide('setup output') as out1:
        lcr = LCR(name, key='', lcr_dir='/home/burnett/work/lcr_data')
    fig1 = lcr.plot()
    fig2 = lcr.compare_all(fmax=fmax, umax=umax)   
    
    return locals()