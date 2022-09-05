from wtlike import *
from wtlike.data_man import DataView
from wtlike.skymaps import HPmap
from utilities.ipynb_docgen import *
from pathlib import Path

from cache_decorator import Cache
names = ['PKS 0208-512', 'PKS 0454-234', 'S5 0716+714', 'OJ 014','PG 1553+113', 'PKS 2155-304',]
periods= [(2.6,0.2), (3.6,0.4), (2.7,0.4), (4.1,0.5), (2.2,0.2), (1.7,0.1)]

@Cache(cache_dir='cache')
def setup_pgms(names):
    pgms = []
    for name in names:
        pgm = WtLike(name).periodogram(tsamp=1)
        pgm.power_spectrum()
        pgms.append(pgm)
    return pgms

@Cache(cache_dir='cache')
def make_count_map(nside=128, bmin=8):
    return DataView().count_map(nside=nside, bmin=bmin)

@ipynb_doc
def intro():
    """
    # Blazar periodicity with Kerr periodograms
    <h6 align="right">{date}</h6>

    This study was motivated by the analysis summarized in this [WAM walkthrouth](https://confluence.slac.stanford.edu/display/SCIGRPS/Weekly+Analysis+Meeting+-+2022?preview=/332664096/349291687/walkthrough_6_24_2022_ppenil.pdf), 
    for a paper "Evidence of Periodic Variability in gamma-ray Emitting Blazars with Fermi-LAT", P. Peñil et al. It covers data to 
    UTC 2020-12-10, or MJD 59193, a total of 12.3 yr. 

    They claim significant periodicity for the six sources PKS 0208512, PKS 0454234, S5 0716+714, OJ 014, PG 1553+113 and PKS 2155304.
    
    Their locations:
    {fig}

    I use data to 2022-07-02, MJD 59762, 13.9 yr total. 
    The analysis uses [wtlike](https://github.con/tburnett/wtlike), displaying a weekly light curve and the 
    low-frequency portion of the Kerr periodogram extracted using hourly time samples. 
    Shown on the periodograms are the measurements from Peñil+ of the periods, derived at least in part from Lomb-Scargle periodograms.
    

    """
    with capture_hide('processing output') as txt:
        pgms = setup_pgms(names)

    nside=128
    cntmap = make_count_map(nside) 
    skymap = HPmap(cntmap,  unit=f'counts per nside={nside} pixel') 

    def six(afig):
        from astropy.coordinates import SkyCoord

        for name in names: #[ '4FGL J1555.7+1111', 'Cyg X-3', ]:
            sk = SkyCoord.from_name(name)
            afig.plot( sk, 'D', color='red', ms=10)
            afig.text( sk, ' '+ name, color='red')


    fig=skymap.ait_plot( log=True, tick_labels=False, pixelsize=0.1, figsize=(12,5), colorbar=False,
         pctlim=(50,98), annotator=six, cmap='Greys', alpha=0.75, title='The six Peñil et al. blazars');  
    
    return locals()


@ipynb_doc
def blazar_variability(name, period, time_bins=(0,0,7), tsamp=1/24, query='p>100 & 0.01 >f>2e-4'):
    """
    #### {name}
    {fig1}
    
    """
    with capture_hide('Setup output') as txt:
        wtl = WtLike(name, time_bins=time_bins)
        pgm = wtl.periodogram(tsamp=tsamp)
        pgm.power_spectrum()
        
    def left(ax):
        wtl.plot(ax=ax, UTC=True)    
    
    def right(ax):
        df = pgm.power_df.query('5e-4< f<0.003').copy()
        df.loc[:,'period'] = 1/df.f/365.25
        # fig2, ax = plt.subplots(figsize=(4,3))
        ax.plot(df.period, df.p1, '.-', label='periodogram')
        ax.set(xlabel='Period (yr)', ylabel='Power $p_1$',
            xscale='log', xticks=[1,2,3,4], xticklabels='1 2 3 4'.split())
        a,b = period
        ax.axvspan(a-b/2,a+b/2, color='lightgrey', label='Peñil et al.')
        ax.grid()
        ax.legend(fontsize=10)
        
    fig1 = plt.figure(figsize=(12,3))
    fig1.subplots_adjust(wspace=0)
    gs = plt.GridSpec(1,10)
    left(  fig1.add_subplot(gs[:6]) )
    right( fig1.add_subplot(gs[7:]) )
         
    # pks= pgm.find_peaks('p1').query(query)
    # period = 1/pks.f/365.25
    # df = pd.DataFrame.from_dict(dict(period=period.round(1), power=pks.p.round()))
    
    return locals()

@ipynb_doc
def power_plots():
    """
    ## Power plots with  $P_b$
    
    The Kerr periodogram analysis checks the variation of the background to the source, and
    produces a "profile" estimate of the source power $P_1$ and the background power $P_b$.
    These plots, patterned after the Kerr paper, show both for low frequency.

    {fig}
    The vertical dashed lines are the frequencies identified by Peñil+.
    """    
    fig, axx = plt.subplots(nrows=2, ncols=3, figsize=(12,9))
    plt.subplots_adjust(wspace=0.45, hspace=0.3)
    with capture_hide('processing output') as txt:
        for name, period, ax in zip(names, periods, axx.flatten()):
            pg = WtLike(name).periodogram(tsamp=1)
            pg.power_plot(ax=ax, xlim=(0,.002))
            ax.text(0.95, 0.90, name, fontsize=12, transform=ax.transAxes, ha='right')
            ax.set(xticks=[0, 0.001, 0.002])
            ax.axvline(1/period[0]/365.25, color='grey', ls='--')
    return locals()

@ipynb_doc
def summary():
    """
    ### Discussion
    
    Note the dramatic discrepancy with the Peñil+ result for PKS 0208-512. But that might be explained by the large
    background component that implies a nearby variable source. The last periodogram shows a pretty clear 1/f noise component.

    A small, but perhaps significant difference for PG 1553+113 may be hard to account for.
    
    The paper shows, in Figure 1, an example Lomb-Scargle periodogram in period units, for PKS 0454-234, which looks quite similar to mine. 
    {penil_fig1}
    
    """
    penil_fig1 = image('penil_fig1.png')
    return locals()

@ipynb_doc
def compare_with_LCR(name='4FGL J0210.7-5101'):
    """
    ## LCR-wtlike comparison for {name}
    
    This uses a light curve file downloaded from the 
    [_Fermi_ Light Curve Repository](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository),
    to compare with the wtlike analysis of the same source.
    {fig}
    """
    gname = name.replace(' ','_')+'*.csv'
    files = list(Path('.').glob(gname))
    assert len(files)>0, f'Did not find a LCR cvs file for {gname}'
    assert len(files)==1, f'Multiple files found for {name}: {files}'
    lcr_file = files[0]
    df = pd.read_csv(lcr_file)
    with capture_hide('WtLike setup output') as txt:
        wtl = WtLike('4FGL_J0210.7-5101')
    
    df.loc[:,'time'] = MJD(df.MET)
    df.loc[:,'flux'] = list(map( lambda x: np.nan if x[0]=='<' else float(x), df.iloc[:,4]) )
    
    fig, (ax1,ax2) = plt.subplots(nrows=2, figsize=(12,6), sharex=True)
    plt.subplots_adjust(hspace=0)
    ax1.plot(df.time, df.flux, '.');
    ax1.grid();
    wtl.plot(ax=ax2, UTC=True);
    
    fig.caption = f'Upper plot:  {df.columns[4]} vs time from the LCR file {lcr_file}; '\
        f'Lower plot: wtlike light curve.'  
    
    return locals()

@ipynb_doc
def phase_plots(reference='2008'):
    """
    ## Phase plots
    Plots of relative flux vs. phase with the Peñil+ measured periods. 
    The phase is measured relative to UTC {reference}.
    {fig}
    """
    def make_phase(index=0, ax=None, ):
        wtl = WtLike(names[index])
        pv = wtl.phase_view(period=periods[index][0]*365.25, reference=reference)
        pv.plot(ax=ax)

    fig, axx = plt.subplots(nrows=2, ncols=3, figsize=(12,8), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0, hspace=0)
    with capture_hide('Setup output') as txt:
        for i, ax in enumerate(axx.flatten()):
            make_phase(i, ax=ax)
            ax.set(xlabel='phase', xlim=(0,1), ylabel='Relative flux' if i%3==0 else '')
    return locals()

@ipynb_doc
def check_nearby(name, radius_cut=4, var_cut=50):
    r""" 
    #### Check {name} for nearby variable sources
    
    Display all within ${radius_cut}^\circ$ with variability >{var_cut}
    {df}
    """
    from utilities.catalogs import Fermi4FGL
    #global near_df
    with capture_hide('load 4FGL output') as out:
        cat =Fermi4FGL();
    near_df = cat.select_cone(name, cone_size=radius_cut, query=f'variability>{var_cut}')
    near_df['name'] = near_df.index
    df = near_df['name sep significance variability'.split()].sort_values('sep')

    return locals()


def skyplot(skymap):
    """
    need to generate a skymap
    ```
    @Cache()
    dvall = DataView()
    def make_count_map(nside=128, bmin=8):
        return dvall.count_map(nside=nside, bmin=bmin)

    cntmap = make_count_map() 
    skymap = HPmap(cntmap,  unit=f'counts per nside={nside} pixel') 
    ```
    """
    def six(afig):

        for name in names: #[ '4FGL J1555.7+1111', 'Cyg X-3', ]:
            sk = SkyCoord.from_name(name)
            afig.plot( sk, 'D', color='red', ms=10)
            afig.text( sk, ' '+ name, color='red')


    return skymap.ait_plot( log=True, tick_labels=False, pixelsize=0.1, figsize=(12,5), colorbar=False,
             vlim=vlim, annotator=six, cmap='Greys', alpha=0.75, title='The six periodic blazars');  



class Sinc():
    def __init__(self, A, freq, delf):
        """
        * A amplitude
        * freq frequency
        * delf -- frequency delta = 1/T 
        """
        self.A, self.freq, self.delf = A,freq, delf
        self.sincsq =  lambda x: A*np.sinc((x-freq)/delf)**2 
    def __call__(self, x):
        return self.sincsq(x)
    
    def lim(self, f=2):
        return  (self.freq-f*self.delf, self.freq+f*self.delf) 
    
    def plot(self, width=2, ax=None, pticks=None, **kwargs):
        xlim = self.lim(width) 
        x = np.linspace(*xlim)
        fig, ax = plt.subplots(figsize=(4,3)) if ax is None else (ax.figure, ax)
        ax.plot(x, self.sincsq(x),'-')
        ax.set(xlim=xlim, **kwargs)
        ax.grid(0.5)
        
        if pticks is None: return
        a,b = np.array(xlim)
        x2 = lambda p: (1/p-a)/(b-a)
        ax.twiny().set(xlabel='Period',
                xticks=x2(np.array(pticks)), 
                xticklabels=[ f'{t}' for t in pticks])  
                
def lowfreqplot(pgm, ax=None, over=None, pticks=None,
                query='f<0.0208', penil=None, **kwargs):
    
    fig, ax = plt.subplots(figsize=(5,3)) if ax is None else (ax.figure, ax)
    yr = 365.25
    df = pgm.power_df.query(query).copy()
    # fig2, ax = plt.subplots(figsize=(4,3))
    x = df.f*yr
    ax.plot(x, df.p1, '.-', label='periodogram', color='blue')
    kw = dict(xlabel='Frequency (1/yr)', ylabel='Power $p_1$')
    kw.update(kwargs)
    ax.set(**kw)
    ax.grid()
    
    if penil is not None:
        a,b = np.array(penil)
        span = 1/np.array([a+b/2, a-b/2])
        ax.axvspan(*span, color='lightgrey', label='Peñil et al.')
    
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
def sinc_overlay(power='p1'):
    """
    ### Periodograms with sinc overlays

    These plots, in frequency units, show the result of overlaying a sinc squared distribution at the
    position of the nearest observed peak to the frequency determined by the Peñil et al. analysls. Its width
    is set by the full time interval, {T:.0f} days

    {fig}
    """
    with capture_hide('processing output') as txt:
        pgms = setup_pgms(names)
        
    T = pgms[0].tspan
    yr=365.25
    fig, axx = plt.subplots(nrows=2, ncols=3, sharex=True, figsize=(15,6))
    plt.subplots_adjust(wspace=0.25, hspace=0.)

    for k,(pgm, name, period, ax) in enumerate(zip(pgms, names, periods, axx.flatten())):
        # pgm.power_plot(ax=ax, xlim=(0,.002))
        ax.text(0.95, 0.90, name, fontsize=12, transform=ax.transAxes, ha='right')
         
        # find a peak close to the Penil estimate
        penil_f = 1/(yr*period[0])
        df = pgm.find_peaks(power=power).query('f<3e-3')

        close = np.abs(penil_f-df.f)< 1/T
        df.loc[:,'close'] = close
        # print(name, penil_f, sum(close), df)
        if sum(close)>0: # 'Not found'
            peak = df[close].iloc[0]
            sinc = Sinc(peak[power], peak.f*yr, yr/T) 
        else: sinc=None
        lowfreqplot(pgm,ax=ax, over=sinc, penil=period,  pticks=[1,1.5,2,3,5,10] if k<3 else None,
                    ylim=(0,None), ylabel='',xlim=(0,1));
    fig.text(0.07, 0.5, f'Power {power}', rotation='vertical', va='center')
    fig.caption="""Low frequency periodograms for the six claimed Peñil et al. periodic blazars. 
        The overlaid red plots are the
            expected sinc functions plotted at the closest peak to the claimed frequency, 
            the uncertainty range of which is denoted by a grey bar."""
                
    return locals()   

@ipynb_doc
def study_0208():
    """ <h2 align="center"> A look at PKS 0208-512</h2>
        
    <h6 align="right">{date}</h6>
    
    Since the Kerr periodogram does not support the existence of oscillations with a period $2.2\pm 0.2$, 
    and the behavior of the low frequency spectrum suggested 1/f noise, and further there is a suggestion
    from the phase plot that the period is actually half that, I looked in more detail to try to 
    check both hypotheses.

    
    ### The analysis
    I redid the last plot in the previous analysis, but expanded the frequency scale, made it semilogy,
    and overplotted an empirical look at the low energy behavior
    {fig}
    A peak is indeed close to the hypothesized half frequency, but the low frequency behavior appears to be exponential, with 
    a 0.25 yr constant, rather than 1/f.
    
    In fact, there are two other peaks equidistant in frequency. To be investigated, 
    especially to see if there is a similar behavior with other blazars.
    """
    with capture_hide('setup printout') as txt:
        names=['PKS 0208-512']
        pgm = setup_pgms(names)[0]
    
    period=[2.6,0.2]
    power='p1'; T= pgm.tspan
    yr=365.25
    penil_f = 1/(yr*period[0])
    df = pgm.find_peaks(power=power).query('f<3e-3')

    close = np.abs(2*penil_f-df.f)< 1/T
    df.loc[:,'close'] = close
    # print(name, penil_f, sum(close), df)
    if sum(close)>0: # 'found'
        peak = df[close].iloc[0]
        sinc = Sinc(peak.p1, peak.f*yr, yr/T) 
    else: sinc=None

    fig, ax = plt.subplots(figsize=(12,5))
    x=np.linspace(0.0,2)
    ax.plot(x, 9e3*(np.exp(-x*4)), '--g',label=r'$\propto e^{-4f}$');
    lowfreqplot(pgm, ax=ax, yscale='log', ylim=(10,None), xlim=(0,3),
               pticks=[0.3,0.4,0.5,1,1.5,2,4], penil=period,
               over=sinc)
    return locals()


def main():

    intro()
    for name, period in zip(names, periods): 
        blazar_variability(name, period)

    power_plots()
    phase_plots()
    sinc_overlay(fignum=5)
    summary()
    study_0208()