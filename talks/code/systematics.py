import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from utilities import profile_plot
from utilities.ipynb_docgen import *
from wtlike.poisson import Poisson
from wtlike.config import  MJD

def measure_systematic(wtl, use_sigma=True, pfactor=1, wfactor=1, ts_min=9):
    
    def adjust(self, pfactor, wfactor):
        p  = self.poiss.p[:] # copy to avoid changing original
        p[0] *= pfactor  # move peak  for the factor
        p[1] /= wfactor**2  # reduce equivalent counts by square of K
        return Poisson(p).sigma_dev(pfactor)
    
    fits = wtl.fits.copy()

    pull = fits.fit.apply(lambda x: adjust(x, pfactor, wfactor))
    ts   = fits.fit.apply(lambda x: x.ts)
    cut  = (ts>ts_min) & ( np.abs(pull)<4 )
    sdev = pull[cut]

    return sdev

def check_systematics(wtl, title=None, **kwargs):
    
#     if not use_sigma:
#         flux = wtl.fluxes.flux.values
#         sig = wtl.fluxes.errors.apply(lambda x: (x[1]-x[0])/2).values
#         pull = (flux-1)/sig
#         cut = ( np.abs(pull)<4 )
#     else:
#         # use the equivalent Gaussian sigma, pehaps modified
#         def pull_kluge(x):
#             p  = x.poiss.p[:] # copy to avoid changing original
#             p[1] /= kfactor**2  # reduce equivalent counts by square of K
#             return Poisson(p).sigma_dev(1)
#         pull = wtl.fits.fit.apply(pull_kluge)
#         ts   = wtl.fits.fit.apply(lambda x: x.ts)
#         cut  = (ts>ts_min) & ( np.abs(pull)<4 )
        
#     # evaluate the pull distribution
#     sdev = pull[cut]
    sdev = measure_systematics(wtl, **kwargs)
    
    label=f'offset {100*sdev.mean():.1f}%\nsystematic {100*(sdev.std()-1):.1f}%'
    
    fig, (ax,ax2,ax3,ax4) = plt.subplots(1,4,figsize=(15,3))
    plt.subplots_adjust(wspace=0.3)
    sk = wtl.source.skycoord.fk5
    ra, dec = sk.ra.deg, sk.dec.deg
    title= title or f' {wtl.time_bins[2]}-day, flux: {wtl.src_flux:.1e}, fk5: ({ra:.1f}, {dec:.1f})'
    fig.suptitle(f'{wtl.source_name} {title}', fontsize=14)
    
    ax.grid(alpha=0.4); xlim=(-4,4)
    ax.hist(pull, bins=np.linspace(*xlim,num=51), histtype='step', lw=2,  label=label)
    ax.set(xlim=xlim, xlabel='pull')
    ax.axvline(0, color='grey')
    ax.legend(loc='lower left', fontsize=12)
    
    profile_plot( wtl.fluxes.e, pull, bins=10, ax=ax2,
        ylim=(-5,5), ylabel='pull', xlabel='exposure/cell', )
    ax2.axhline(0,color='grey')

    try:
        profile_plot(wtl.fluxes.ctm.values, pull, bins=np.linspace(0.4,1.0,7), ax=ax3,
                xlabel='<cos theta>', ylabel='pull', ylim=(-5,5), xlim=(0.4,1.0))
        ax3.axhline(0, color='grey')
    except:
        print('fail pull vs. ctm plot')
    
    profile_plot(wtl.fluxes.e, wtl.fluxes.ctm,  bins=10, ax=ax4,
                 ylabel='<cos theta>', xlabel='exposure/cell', ylim=(0.4,1.0) )
    fig.set_facecolor('white')
    return fig

def systematic(wtl, k=1, ts_min=9):
    """ wtl: a WtLike object
    Return the mean and std of the max L fits 
    k is a kluge factor to adjust the individual likelihood functions in order to "blow up" the 
    estimated errors by the factor.
    """
    def kluge(x):
        p  = x.poiss.p[:] # copy to avoid changing original
        p[1] /= k**2  # reduce equivalent counts by square of K
        return Poisson(p).sigma_dev(1)
    pull = wtl.fits.fit.apply(kluge)
    ts   = wtl.fits.fit.apply(lambda x: x.ts)
    cut  = (ts>ts_min) & ( np.abs(pull)<4 )
    sdev = pull[cut]
    return sdev

def pull_plot(wtl, title='', ax=None, **kwargs):
    """use the sigma_dev poisson function 
    """
    sdev = comupute_systematics(wtl, **kwargs)

    fig, ax = plt.subplots(1,1, figsize=(5,3)) if ax is None else (ax.figure,ax)

    ax.hist(sdev, np.linspace(-5,5,25), histtype='step', lw=2, density=True,
            label=f'mean  {sdev.mean():.2f}\nwidth {sdev.std():.3f}')
    dom = np.linspace(-5,5,101)
    ax.plot(dom, stats.norm.pdf(dom), '--', lw=2, label='Normal(0,1)')
    ax.grid(alpha=0.5)
    ax.set(xlabel='normalized deviation', title=title )
    ax.legend(loc='upper left', fontsize=10);
    txt = f'{wtl.source_name}'
    if kfactor!=1:
        txt += f'\n(k={kfactor})'
    ax.text(0.98,0.97, txt, transform=ax.transAxes, ha='right', va='top', fontsize=12)

def pull_fun(x, k=1):
    p  = x.poiss.p[:] # copy to avoid changing original
    p[1] /= k**2  # reduce equivalent counts by square of K
    return Poisson(p).sigma_dev(1)

def make_phase_df(wtl, tw, period, nbins=25):

    view = wtl if wtl.time_bins[2]==tw else wtl.view(0,0,tw)
    df = view.fits['t tw n e fit'.split()].copy()

    pull = df.fit.apply(pull_fun).values
    df.loc[:,'pull'] = pull
    df.loc[:,'sigma'] =df.fit.apply(lambda f: (f.errors[1]-f.errors[0])/2)
    return phases(df, period, nbins)

def phases(df, period, nbins=50):
    """
    Input: a DF with columns t, pull, sigma 
    Returns: a phased DF with x, y, yerr to make phase plot
    """
    dfc = df.copy()['t pull sigma'.split()]    

    phase = np.mod(dfc.t-MJD(2008), period).values/period
    bins =  np.linspace(0,1,nbins+1)
    phase_bin = np.digitize(phase, bins)
    counts, _ = np.histogram(phase, bins)

    dfc.loc[:,'phase_bin'] = phase_bin
    
    # create a pivot table on the phase bin, with averages
    dfp = dfc.pivot_table(values=['pull','sigma'], index='phase_bin')
    return pd.DataFrame(dict(x = (bins[1:]+bins[:-1])/2,
                             y = (dfp.pull) * dfp.sigma,
                             yerr= dfp.sigma/np.sqrt(counts),)
                       )
@ipynb_doc
def correlation_vs_offset(wtl, interval=1):
    """
    ## Correlation vs. offset
    The adjacent day correlation, clearly indicates problems with the exposure. Here I look
    at the correlation coefficient vs. the offset. 
    {fig1}
    
    And here is an FFT of the correlation vs. offset
    {fig2}
    The vertical dashed lines are at periods of {tpr} days and half that.
    """
    pull = wtl.view(0,0,interval).fits.fit.apply(pull_fun).values
    N = len(pull)
    M = min(N,800)
    r = range(1,M)

    tpr = 53.05
    corr = np.array([np.sum(pull[:N-i]*pull[i:])/(N-i) for i in r])

    def correlation_plot(xmax=400):
        fig,ax = plt.subplots(figsize=(12,3), num = 1)
        ax.plot( r, 100*corr, '-');
        ax.grid(alpha=0.5);
        ax.set(xlabel='offset [days]', xlim=(0,xmax), ylabel='correlation (%)', ylim=(-5,15)) 
        ax.axhline(0, ls='-', color='grey');
        ax.axvline(365.25, ls='--', lw=2,color='orange', label='1 year')
        ax.legend()
        fig.set_facecolor('white')
        return fig

    def fft_plot( ):
        yr = 365.25/interval
        df = 1/M * yr
        output_array = np.fft.rfft(corr)
        power =np.square(np.absolute(output_array))# / norm
        fig,ax=plt.subplots(figsize=(12,3), num=2)
        ax.plot(df*np.arange(len(power)), power, 'D', ms=7);
        ax.set(xlim=(0,0.05*yr), xlabel='Frequency $(\\mathrm{year^{-1}})$',
              yscale='log', ylim=(1e-2,None), ylabel='Power',
              xticks=np.linspace(0,20,11));
        ax.grid(alpha=0.75);
        
        for x in (yr/tpr, 2*yr/tpr):
            ax.axvline(x, ls='--', lw=2, color='orange')
        fig.set_facecolor('white')
        return fig

    fig1=figure(correlation_plot())
    fig2=figure(fft_plot())
    return locals()


def pull_plot(wtl, ax=None,  title='',nbins=25, wfactor=1, **kwargs):
    """use the sigma_dev poisson function 
    """
#     def adjust(self, pfactor, wfactor):
#         p  = self.poiss.p[:] # copy to avoid changing original
#         p[0] *= pfactor  # move peak  for the factor
#         p[1] /= wfactor**2  # reduce equivalent counts by square of K
#         return Poisson(p).sigma_dev(pfactor)
    
#     fits = wtl.fits.copy()

#     pull = fits.fit.apply(lambda x: adjust(x, pfactor, wfactor))
#     ts   = fits.fit.apply(lambda x: x.ts)
#     cut  = (ts>ts_min) & ( np.abs(pull)<4 )
#     sdev = pull[cut]
    sdev = measure_systematic(wtl, wfactor=wfactor, **kwargs)

    fig, ax = plt.subplots(1,1, figsize=(5,3)) if ax is None else (ax.figure, ax)

    ax.hist(sdev, np.linspace(-5,5,nbins+1), histtype='step', lw=2, density=True,
            label=f'mean  {sdev.mean():.3f}\nwidth {sdev.std():.3f}')
    dom = np.linspace(-5,5,101)
    ax.plot(dom, stats.norm.pdf(dom), '--', lw=2,)# label='Normal(0,1)')
    ax.grid(alpha=0.5)
    ax.set(xlabel='normalized deviation', title=title )

    ax.text(0.05, 0.97,f'mean  {sdev.mean():.3f}\nwidth {sdev.std():.3f}',
        transform=ax.transAxes, ha='left', va='top', fontsize=10 )
    #ax.legend(loc='lower center', fontsize=10);
    txt = f'{wtl.source_name}'
    if wfactor!=1:
        txt += f'\n(wfactor={wfactor})'
    ax.text(0.98,0.97, txt, transform=ax.transAxes, ha='right', va='top', fontsize=12)
    
@ipynb_doc
def pull_plots(views, sizes=[0.5,3,7]):
    """
    Plots of the normalized deviation, or pull, for various cell sizes:
    {fig1}
    """
    fig, axx = plt.subplots(1,len(sizes), figsize=(12,3), sharex=True)
    fig.set_facecolor('white')
    fig1 = figure(fig)

    for ax, w, cell_size  in zip(axx.flatten(), 
                            (views),
                            sizes):
            pull_plot(w, ax=ax, title=f'{cell_size} day cells')
            ax.set(ylim=(0,0.5))
    return locals()

@ipynb_doc
def phase_plot(wtl,  period, nbins=50):
    """### Phase analysis

    #### Plot the phases
    {fig1}
    This periodicity in the correlations suggests that
    it could be corrected.

    """

    # Now make a DF, with daily bins, with fit objects, and add pull and sigma
    pdf = make_phase_df(wtl, 0.5, period, nbins)

    fig,ax = plt.subplots(figsize=(6,3))
    ax.errorbar(x=pdf.x, y=100*pdf.y, yerr=100*pdf.yerr,  fmt=' ', marker='o');
    ax.set(xlabel=f'phase (period={period} days)', xlim=(0,1), ylim=(-3,3),
           ylabel='flux systematic (%)',)
    ax.axhline(0, color='grey')
    ax.grid(alpha=0.5)

    fig1=figure(fig, caption=f'Flux systematic for {wtl.source_name}')

    return locals()

# %run systematics.py
# sdevs = [measure_systematic(wtl_dict[name]) for name in bright_names]; 

# rms_vals = [x.std() for x in sdevs]; 
# xdf= pd.DataFrame(dict(G100=eflux, rms_vals=np.array(rms_vals).round(2)))
# xdf.sort_values('G100', ascending=True); 