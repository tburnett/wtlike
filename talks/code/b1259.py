from wtlike import *
# from astropy.coordinates import SkyCoord
# from wtlike.skymaps import *
# from cache_decorator import Cache
from utilities.ipynb_docgen import *
plt.rc('font', size=14)

weekly=None # the WtLike object set by header()  
pg = None # TimeSeries 
clear = False

def get_global(which):
    return globals().get(which, None)

class B1259Periastron(WtLike):
    """
    """    
    tp, period = 55544.694, 1236.7243
    
    def __init__(self, **kwargs):
        super().__init__('PSR B1259-63',**kwargs)
        pd.set_option('display.precision', 3)#, 'display.colheader_justify','left')
        pd.set_option('display.width',100)
        plt.rc('font', size=14)

    @property
    def mjd_dates(self):
        return [self.tp+n*self.period for n in range(4)]
         
    def date_info(self):
        return pd.DataFrame([dict(MJD=round(d,3), UTC=UTC(d)) for d in self.mjd_dates ])
    
    def stacked_plots(self, fignum=3, ylim=(2,200), ts_min=4, p0=0.5, xlim=(-50,200)):
        self.sixhr= dailies = [self.view(round(tz)+xlim[0], round(tz)+xlim[1], 0.25) for tz in self.mjd_dates]
        self.bb_views = [daily.bb_view(p0=p0,) for daily in dailies]
        
        fignum=1
        fig, axx = plt.subplots(4,1, sharex=True, sharey=True, figsize=(15,10), num=fignum)
        plt.subplots_adjust(hspace=0, top=0.92)
        for i, (bbv, ax, tzero) in enumerate(zip(self.bb_views, axx.flatten(), self.mjd_dates)):
            tzero = round(tzero)
            bbv.plot(ax=ax, show_flux=True, tzero=tzero, xlim=xlim,
                         log=True, source_name=f'{UTC(tzero)[:10]}', ylabel='',
                        legend_loc='upper left' if i==0 else 'none',
                       error_pixsize=0);
            ax.axvline(0, ls='--', color='grey')
        
        ax.set(xlabel='Days past periastron')
        fig.suptitle('PSR B1259-63 periastron behavior')
        fig.text(0.06, 0.5,'Count flux ($\mathrm{10^{-6}\ cm^{-2}\ s^{-1}}$)', rotation='vertical', va='center' );

        fig.width=600
        return fig


@ipynb_doc
def header():
    """
    # B1259-63 Analysis
    <h6 align="right"> {date}</h6>
    > Create a B1259-63 (aka PSR J1302-6350) light curve using Bayesian Blocks  

    * [HESS high-energy results](https://arxiv.org/pdf/astro-ph/0506280.pdf)

    * [Previous *Fermi* paper](https://arxiv.org/pdf/1912.05868.pdf)

    {outp}
    """
    global weekly
    with capture_hide('Output from B1259-63 setup: create cells, fit each, run BB partition, fit partitions.') as outp:
        if weekly is None or clear:
            weekly = B1259Periastron(time_bins=(0,0,7), clear=clear)#.bb_view()
        else:
            print('(already done)')
    return locals()
header()

def data_check():
    from wtlike.data_man import update_data
    clear = update_data()
    print(('Will' if clear else 'Will not'),'update')
    with Timer() as t:
        weekly = self=B1259Periastron( clear=clear)
        print(t)
    return clear




@ipynb_doc
def B1259( clear=False):
    r"""
    ## Fit to all data

    <p style="text-align:right;">{date}</p> 
    
    Create a `WtLike` object with all the data
    
    {outp}
 
    ## The full weekly-interval light curve, showing the BB partitions
    {out2}
    {fig2}
    Table of fits (note that the "flux" column is relative to the 12-year count flux measurement
    6.7e-9 cm**-2 s**-1.)
    {bbf}
    
    ## Expand about each periastron
 
    #### Periastron dates

    Assuming {period:.2f}-day ({period_yr:.2f} year) orbital period, the MJD and UTC values are:
    
    {utc}
 
    Expand the above, with 6-hour bins, following those dates:


    {fig3}
    """
    global weekly #  make availlable for follow-up cells

    with capture_hide('Output from analysis: create cells, fit each, run BB partition, fit partitions.') as outp:
        if weekly is None or clear:
            weekly = B1259Periastron(time_bins=(0,0,7), clear=clear)#.bb_view()
        else:
            print('(already done)')
        bb = weekly.bb_view(0.5)
        period = weekly.period
        period_yr = period/365.25
    
    bbf = monospace(str(bb.fluxes), 'BB fit table')#, open=False)

    # fig 2: full light curve
    fig2 = figure( bb.plot(fignum=2, figsize=(18,4), UTC=True, title='Full weekly light curve'),
                  width=800)
    
    # fig 3 -- stack light curves
    utc = monospace(str(weekly.date_info())) 
    with capture_hide('Output from periastron analyses') as out2:
        fig3 = figure( weekly.stacked_plots(  fignum=3), width=800)
        
    return locals()

@ipynb_doc
def recent():
    """
       
    ## Recent details: 1-day bins for 2 weeks, 6-hr bins for 1 week 
    {out4}{fig4}
    
    {out5}{fig5}
    """
    global weekly
    # fig 4 -- daily, last two weeks
    with capture_hide(f'Output, with light curve table, from selecting the last two weeks') as out4:
        recent_wk = weekly.view(-14,0,1)
        print(recent_wk.fluxes)
        fig4=figure(
            recent_wk.plot(show_flux=True, ylim=(-.1,2), title='Daily bins for last 2 weeks',
                       xlabel=f'MJD, until {UTC(weekly.stop)[:10]}',fignum=4),
            caption=None, width=600)
  
    # fig 5 -- last week, 6-hr bins
    with capture_hide(f'Last week, 6-hr bins') as out5:
        lastwk = weekly.sixhr[-1].view(-7,0,0.25, exp_min=0.3)
        fig5= figure(
            lastwk.plot_with_exposure(show_flux=True,ylim=(0,2),  
                                      title=f'Last week, 6-hr bins, ending {UTC(weekly.stop)[:10]}'),
           caption=None, width=600 )
        print(lastwk.fluxes)
        
    return locals()

def lowfreqplot(pgm, ax=None, over=None, pticks=None,
                query='f<0.0208', power='p1', **kwargs):
    
    fig, ax = plt.subplots(figsize=(5,3)) if ax is None else (ax.figure, ax)
    yr = 365.25
    df = pgm.power_df.query(query).copy()
    # fig2, ax = plt.subplots(figsize=(4,3))
    x = df.f*yr
    ax.plot(x, df[power], '.-', label='periodogram', color='cornflowerblue')
    kw = dict(xlabel='Frequency (1/yr)', ylabel=f'Power {power}')
    kw.update(kwargs)
    ax.set(**kw)
    ax.grid()
    
    if over is not None:
        a,b = over.lim()
        xx = x[(x>a) & (x<b)]
        ax.plot(xx, over(xx), '--r', label=f'sinc at {1/(over.freq):.2f} yr')

    # ax.legend(fontsize=10)
    if pticks is None: return
    a,b = np.array(ax.get_xlim())
    x2 = lambda p: (1/p-a)/(b-a)
    ax.twiny().set(xlabel='Period (yr)',xlim=(0,1),
            xticks=x2(np.array(pticks)), 
            xticklabels=[ f'{t}' for t in pticks])  
    ax.legend(fontsize=12)   

@ipynb_doc
def b1259_fft():
    """### B1259-63 Periodogram
    
    {fig1}
    
    ### Low frequency details
    {fig2}
    
    The orange lines denote the orbital frequency and harmonics.
    """
    global pg
    with capture_hide('FFT analysis output') as out1:
        pg = weekly.periodogram()
    
    fig1, ax1 = plt.subplots(figsize=(12,4))
    pg.power_plot(ax=ax1, pmax=50, xlim=(0,6));
    
    fig2, ax = plt.subplots(figsize=(12,4), ) 
    forb =365.25/weekly.period
    for i in range(1,10):
        ax.axvline(i*forb, ls='--', color='orange') #, label='B1259 orbital frequency');
    lowfreqplot(pg,ax=ax, power='p0', pticks=[0.3,0.4,0.5,1,1.5,2,4], xlim=(0,3))
    return locals()    

@ipynb_doc
def peak_study(q='p1>19 & f>1.1'):
    """
    ## Peaks
    
    Run the peak-finding algorithm on the periodogram, for frequencies above 1.1. 
    {peaks}
    
    ### Distribution of peak values for frequencies >1.1
    {fig1}
    
    ### Look at the two with power > 24
    {fig2}
    """
    global pg
    if pg is None:
        with capture_hide('FFT analysis output') as out1:
            pg = weekly.periodogram()

    pd.set_option('display.precision', 4)
    peakdf = pg.find_peaks('p1')
    with capture_show(f'Peaks which satisfy {q}') as peaks:
        print(peakdf.query(q))
        
    p1 = peakdf.query('f>1.1').p1.values
    fig1, ax = plt.subplots(figsize=(5,3))
    hkw = dict(bins= np.linspace(0,35,71), histtype='stepfilled',log=True)
    ax.hist(p1, color='cornflowerblue', **hkw);
    ax.hist(p1[p1>24], color='red', label='P1>24', **hkw)
    ax.set(xlabel='Power p1',xlim=(0,35));
    ax.grid(); ax.legend();
        
    fig2, (ax1,ax2) = plt.subplots(ncols=2, figsize=(8,3), sharey=True)
    plt.subplots_adjust(wspace=0)
    f1 = 4.3419; delf=1.5e-3
    pg.power_plot(ax=ax1, pmax=50, xlim=(f1-delf,f1+delf), ylim=(-15,40))
    f2 = 5.6586# bpk.f[-1]
    pg.power_plot(ax=ax2, pmax=50, xlim=(f2-delf,f2+delf),ylim=(-15,40))
    ax2.set(ylabel='');
    return locals()

@ipynb_doc
def spectogram():
    """
    ## Spectogram

    Finally here is a spectogram

    {fig}
    """
    pg = get_global('pg');
    fig, ax = plt.subplots(figsize=(8,8))
    sg = pg.spectogram(interval=30)
    sg.sgplot(ax=ax)
    return locals()

@ipynb_doc
def examine_peaks():
    """ ### Examine the peaks
    
    Here we subdivide the data into two 2500-day intervals, perform FFT analysis on each,
    and compare the two results.
    
    {fig1}
    {fig2}
    
    The lower frequency one has support in both intervals. Since it is clearly not 
    associated with B1259, we'll look for a neighboring source.
    """
    
    ts= get_global('pg')
    if ts is None:
        ts = get_global('weekly').periodogram()
    sg = ts.spectogram( interval=2500)

    # sg.sgplot(imshow_kw=dict(vmin=0,vmax=10))
    fpk = [4.3419, 5.6586]
    d = 0.003
    figs=[]
    for fignum,f in enumerate(fpk):
        fig, (ax1,ax2) = plt.subplots(ncols=2, figsize=(10,3), num=fignum)
        fig.suptitle(f'Examine peak at {f}')
        pk = sg.loc[:, f-d:f+d]
        pk.sgplot(ax=ax1, imshow_kw=dict(vmin=0,vmax=20), grid='orange')
        pk.iloc[0,:].plot(ax=ax2)
        pk.iloc[1,:].plot(ax=ax2)
        plt.grid(); ax2.axvline(f, ls=':');
        figs.append(fig)

    fig1, fig2= figs
    return locals()