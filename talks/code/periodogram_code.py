
from utilities.catalogs import *
#----------------------------------------------------------------------   
from wtlike import *
from utilities.ipynb_docgen import *

from utilities.catalogs import *
from wtlike.time_series import power_spectrum_plot

# a dict of periodogram tables -- global for these functions
pgs={}

def set_pgs( x):
    global pgs
    pgs = x
    
def get_pgs():
    return pgs
    
@ipynb_doc
def introduction():
    """
    # Pulsar Variability Study
    
    ## Introduction
    
    This study is inspired by Matthew Kerr, applying his periodogram analysis as documented 
    [here](https://arxiv.org/pdf/1910.00140.pdf), in which see Section 5, and implemented in [wtlke](https://tburnett.github.io/wtlike/). 
    
    The idea is to use the periodogram analysis to find variability in pulsars: are there any besides J2021+4026?
    
    Matthew's GI proposal, "SEARCHING FOR FAST(ER) VARIABILITY IN γ-RAY PULSARS" notes that previous
    searches are only sensitive at > 1 y timescales, and suggests looking at > 1 min time scales. 
    
    Also see the [talk by Jeff Scargle](https://confluence.slac.stanford.edu/display/SCIGRPS/Virtual+Collaboration+Meeting+Mar+21-25%2C+2022?preview=/332672681/341257592/lat_collaboration.pdf).
    """
    return locals()

    
def load_periodograms(names, cache_file='cache/periodograms.pkl' ):
    """Manage a cache of periodograms, as a dict of DataFrames 
    
    Saved format is actually np.record, avoiding problem with pandas versioning
    """
    import pickle
    from wtlike import WtLike, Config
    
    global pgs
    
    if Path(cache_file).is_file():
        with open(cache_file, 'rb') as infile:
            t = pickle.load(infile)
        pgs =  dict((key, pd.DataFrame.from_records(value)) for (key,value) in t.items())
        print(f'Loaded from file {cache_file}')
        set_pgs(pgs)
        return
    
    # create dict of dataframes to return
    print(f'Regenerating list of {len(names)} periodograms')
    spect = lambda name: WtLike(name, config=Config(verbose=0)).periodogram().power_spectrum()   
    pgs = dict( (name,  spect(name)) for name in names)
    
    # save it as records
    with open(cache_file, 'wb') as outfile:  
        pickle.dump( dict((k, v.to_records(index=False)) for (k,v) in pgs.items() ),  outfile)
        print(f'Dumped to file {cache_file}')
    set_pgs(pgs)  
    
@ipynb_doc
def get_bright_pulsars(N=20):
    """ 
    ## Setup 
    
    ### Extract a list of 4FGL sources identified as pulsars
    Select the top 20 by flux, excluding the Crab.
    {out1}
    The list: {names}
    
    ### Run the FFT analysis
    
    {out2}
    """

    with capture_show('output') as out1:
        fcat = Fermi4FGL()
        fcatpsr = fcat.query('class1=="PSR" | class1=="MSP"').copy()
        fcatpsr.loc[:,'catname'] = fcatpsr.index
        fcatpsr.index = fcatpsr.assoc1_name
        fcatpsr.index.name='psr_name'
        print(f'\t{len(fcatpsr)} identified as pulsars')

    bright = fcatpsr.sort_values('eflux', ascending=False)[:N+1]
    bright_names =  list(filter( lambda n: n!='PSR J0534+2200', bright.assoc1_name.values))
    names = sorted(bright_names)

    with capture_show('printout') as out2:
        load_periodograms(bright_names)

    return locals()

@ipynb_doc
def full_spectra(pmax=50, utc_end='UTC 2022-05-26 09:53'):
    """
    ## Plot full periodograms
    <h6 align="right"> {date} </h6>
    
    The analysis is of {N} brightest pulsars, excluding the Crab.
    The dataset has all data until {utc_end}, nearly 14 years. The subsequent FFT analysis uses 1-hour
    sampling, for a Nyquist frequency of 6 per day.
    {fig}
    """

    N = len(pgs)    
    fig, axx = plt.subplots(ncols=2, nrows=N//2, figsize=(20,20), sharex=True, sharey=True)
    plt.subplots_adjust(hspace=0, top=0.95, wspace=0.05)
    fig.suptitle('Pulsar periodograms', fontsize=20)
    fig.text(0.05, 0.5, r'$\leftarrow P_b \ \ \ \ P_1 \rightarrow $', rotation='vertical', fontsize=16)

    for (name, p), ax in zip(pgs.items(),  axx.flatten()):
        ax.text(0.05, 0.9, name,  transform=ax.transAxes, fontsize=10)
        power_spectrum_plot(p, ax=ax, pmax=pmax)
        ax.set(ylabel='', xlim=(0,6))

    return locals()

@ipynb_doc
def low_spectra(ncols=5, xlim=(1e-4, 0.03), ylim=(-62,62), profile=True):
    """
    ## Low-frequency power spectra
    
    <h6 align="right"> {date} </h6>
    
    The following plots show {plottype}.
    {fig}
    The vertial dashed lines are for the following periods that can affect data acquisition:
    {taglist}
    """

    plottype = 'the profile P_1, with P_b in orange' if profile else 'P0 only'
    tags = np.array([1/365.25, 4/365.25, 1/53.05])
    tagnames = ['year', 'quarter year', 'Fermi precession']
    taglist = pd.DataFrame(dict(description=tagnames, freq=tags, period=(1/tags)))
    N = len(pgs)
    nrows=(N+ncols-1)//ncols
    fig, axx = plt.subplots(ncols=ncols, nrows=nrows, figsize=(1+ncols*3,1+nrows*2.5 if profile else 1+nrows*2), 
                            sharex=True, sharey=True)
    plt.subplots_adjust(hspace=0, top=0.94, left=0.1, wspace=0.15, bottom=0.1)
    fig.suptitle('Low-frequency pulsar power spectra', fontsize=20)
    fig.text(0.5, 0.05, r'$\mathrm{ Frequency\ (d^{-1})}$', ha='center', fontsize=16)
    fig.text(0.05, 0.5, r'$\leftarrow P_b \ \ \ \ P_1 \rightarrow $' if profile else \
                      r'$ P_0 $', rotation='vertical')
    for (name, power_df), ax in zip(pgs.items(), axx.flatten()):
        power_spectrum_plot(power_df, ax=ax, pmax=50, profile=profile,
                            xlim=xlim, xscale='log', ylim=ylim)
        ax.text(0.08, 0.90, name,  transform=ax.transAxes, fontsize=10)
        ax.set(ylabel='',xlabel='')
        ax.set_xticks([1e-3,1e-2], labels=['$10^{-3}$', '$10^{-2}$'])
        for tag in tags:
            ax.axvline(tag, color='lightgreen', ls='--')
  
    return locals()

@ipynb_doc
def low_freq_selected(choose=(3,5), profile=False, **kwargs):
    """
    ### Look at the J2021 guys
    
    profile= {prof}
    {fig}
    This compares the known variable J2021+4026 with the nearby J2021+3651.
    
    The latter shows no variability in the light curve, yet also has apparently significant low 
    frequencies, beside the 1-year peak.
    """
    prof = profile
    pgs = get_pgs()
    keys = list(pgs.keys())
    fig, axx = plt.subplots(nrows=2,ncols=1, figsize=(12,6), sharey=True)
    fig.suptitle('Low-frequency periodograms')
    plt.subplots_adjust(hspace=0.01)
    tags = np.array([1/365.25, 4/365.25, 1/53.05])
    
    kw = dict(ylabel='', #xticks=[1e-3, 1e-2], xticklabels='0.001 0.01'.split(),
             xlim=(0.9e-3, 0.03), xscale='log', ylim=(0,75),)
    kw.update(kwargs)
    
    for k, ax in zip(choose,axx.flatten()):
        name = keys[k]
        ax.text(0.5, 0.9, name, transform=ax.transAxes, ha='center')
        power_spectrum_plot(pgs[name], ax=ax,  profile=profile, **kw)
        for tag in tags:
            ax.axvline(tag, color='lightgreen', ls='--')
        ax.xaxis.set_major_formatter(
          lambda val,pos: { 0.001:'0.001', 0.01:'0.01', 0.1:'0.1', 1.0:'1', 10.0:'10',}.get(val,''))
    return locals()

def peak_finder(x, y, cut=0):
    """
    Determine positions and values of the peaks in an array

    - x -- equally-spaced independent variable
    - y -- dependent variable
    - cut -- filter so peaks above this

    return a DataFrame with columns x and y for the positions and values of peaks in the array
    """
 
    deltax = x[1]-x[0]

    # get gradient 
    g = np.gradient(y)

    # index of first point before peak, where gradient changes sign
    qp = (g[:-1]>0) & (g[1:]<0)
    qi = np.arange(len(g)-1)[qp]

    # find peaks, interpolating the gradient to zero
    g1,g2 = g[qi],g[qi+1]
    xp = x[qi] + g1/(g1-g2)*deltax

    # estimate peak value with largest of two points
    yp = np.max(np.vstack([y[qi],y[qi+1]]),axis=0)   
    return pd.DataFrame.from_dict(dict(x=xp, y=yp)).query(f'y>{cut}')

def find_peaks(power_df, power='p0'):
    """
    Determine positions and values of the peaks in the given power spectrum

    - power: Select p0, p1, or pb

    return a DataFrame with columns f and p for the frequencies and power values
    """
    df = power_df #self.power_spectrum()
    expect = 'p0 p1 pb'.split()
    assert power in expect, f'TimeSeries.find_peaks: {power} not one of expected {expect}'
    y = df[power].values
    x = df.f.values
    df = peak_finder(x,y)
    df.loc[:,'period'] = 1/df.x
    df.rename(columns=dict(x='freq', y='power'), inplace=True)
    return df

@ipynb_doc
def check_peak(name, pmax=50):
    """## Check a peak
    
    For source {name}, check the peak, with power {ppeak:.1f} at frequency {fpeak:.4f}.
    {fig1}
    """
    pgs = get_pgs()
    z = pgs[name]; 
    deltaf = z.f[0]
    peaks = peak_finder(z.f, z.p1).query('y>25& x>0.1')
    fpeak = peaks.iloc[0].x
    ppeak = peaks.iloc[0].y

    # ax.plot(z.f, z.p1, '.-', );
    def left(ax):
        power_spectrum_plot(z,ax=ax, pmax=pmax)
        ap = dict(arrowstyle='->',color='k', lw=3)
        ax.annotate('', xy=(fpeak, 0.85*pmax), xytext=(fpeak, pmax),# transform=ax.transData,
                        arrowprops=ap);
        
    def right(ax):
        power_spectrum_plot(z,ax=ax, pmax=pmax, 
                            xlim=(fpeak-50*deltaf, fpeak+50*deltaf ))
        
    fig1  = figure(plt.figure(figsize=(12,3)))
    fig1.suptitle(name)
    plt.subplots_adjust(wspace=0)
    gs = plt.GridSpec(1,10)
    left(  fig1.add_subplot(gs[:7]) )
    right( fig1.add_subplot(gs[8:]) )


    return locals()

@ipynb_doc
def peak_distribution(fmin=0.01, fNyquist=6):
    r""" ## Peak distribution
    
    This combines the distributions of the power of the peaks of the periodograms
    of {N} pulsars, those with the highest flux excluding Geminga, Vela, and the Crab.
    Frequencies were selected from {fmin:.2f} to {fNyquist} cycles/day.
    
    {fig}
    
    The overlaid function is empirical, and perhaps can be derived for the distribution of 
    peak power values. 

    The largest peak's power found is {maxp:.1f} for {pname}. The prediction from the 
    empirical distribution above 30 is {over30:.1f}, consistent with the detected number, {N30}. 
    However, it predicts {over25:.1f} above 25, vs. the observed {N25}.

    <h6 align="right"> Generated {date} </h6>
    """
    pgs = get_pgs()
    pname = ''
    p = np.empty(0)
    N = 0; maxp = 0
    for i, (k,v) in enumerate(pgs.items()):
        if i<2: continue
        dfpk = peak_finder(v.f, v.p1).query('x>0.1')
        t = dfpk.query(f'y>{maxp}')
        if len(t)>0: #print(k,t )
            pname = k
            maxp = max(t.y)            
        p = np.append(p, dfpk.y.values)
        N+=1
    # distribution
    xi = 9/4
    N10 = sum(p>10)
    N25 = sum(p>25)
    N30 = sum(p>30)
    fn =lambda x: N10* np.exp(-(x-10)/xi)/xi
    Fn = lambda x: xi* fn(x)
    over30 = Fn(30); over25=Fn(25)
    
    fig, ax = plt.subplots(figsize=(6,4))
    fig.caption=f'Distribution of the profile power P_1 for the peaks in the {N} periodograms.'
    kw = dict(bins=np.linspace(0,35,36), histtype='stepfilled',
              facecolor='lightgrey',edgecolor='blue', log=True, lw=2)
    ax.hist(p, label=f'{N} pulsars', **kw);
    ax.grid();
    x = np.linspace(0,35)
    ax.set(ylim=(0.6,None),  xlabel=r'$\mathrm{Power}\ P_1$')

    ax.plot(x, fn(x), label=r'$\propto \exp(-\frac{4}{9} P_1)$');
    ax.legend();
    
    return locals()