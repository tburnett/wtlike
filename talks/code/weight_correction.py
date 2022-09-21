#----------------------------------------------------------------------   
# from tkinter import FALSE
from wtlike import *
from utilities.ipynb_docgen import *

from utilities.catalogs import *
from cache_decorator import Cache

wtl1, wtl2, wtlx, bb2 = None,None,None, None
near_df = None

def get_global(name):
    return globals()[name]


# S1, S2='PSR J0007+7303', '4FGL J0019.6+7327'

S1, S2, interval = 'PSR J2032+4127', 'Cyg X-3', 30
# S1, S2, interval = 'PSR J1957+5033', '4FGL J1955.4+5132', 30
# S1, S2, interval = 'PSR J2028+3332', '4FGL J2025.3+3341', 30
# S1, S2, interval = 'PSR J2021+3651', 'P88Y5611', 30


# @Cache(cache_dir='cache')
def setup_wtlike(interval=interval, clear=False):
    
    global wtl1, wtl2, bb2
    if wtl1 is None:
        time_bins=(0,0,interval)
        wtl1 = WtLike(S1, time_bins=time_bins, clear=clear)
        wtl2 = WtLike(S2, time_bins=time_bins, clear=clear)
        bb2  = wtl2.bb_view()
    assert wtl1 is not None
    return wtl1, wtl2, bb2

def get_wtl():
    return get_global('wtl1')
    # global wtl1,wtl2, bb2
    # if wtl1 is None:
    #     setup_wtlike()
    # return wtl1, wtl2, bb2

@ipynb_doc
def comment(text):
    """{text}"""
    return locals()

@ipynb_doc
def wc_intro():
    r"""# Weight correction factor
    <div align="right">{date}</div>
    For a presumably constant source $S_1$ with a variable neighbor $S_2$, I've
    [derived](https://tburnett.github.io/wtlike/sources.html#Accounting-for-variations-from-neighboring-sources) the following:
    
    $$w'_1 = \frac{w_1}{1+\alpha_2\ w_2}\ \ ,   $$
    where  $w_1$ and $w_2$ are the corresponding weights for a common pixel/band, and $\alpha_2$ is the source flux factor for $S_2$. 
    
    We will test this, with $S_1$ as the pulsar {name1} and $S_2$ the variable source {name2}. 
    The separation is ? deg.

    {fig}
    """
    name1, name2=S1,S2
    global wtl1, wtlx, bb2
    with capture_hide('Wtlike and BB processing output') as out1:
        wtl1,wtl2,bb2 = setup_wtlike()
        wtlx = wtl1.reweighted_view(bb2)
        
    fig, (ax1,ax2) = plt.subplots(nrows=2, figsize=(12,5), sharex=True)
    fig.caption=f'Upper plot: the pulsar {S1}; lower plot: Bayesian block output of the variable source.'
    plt.subplots_adjust(hspace=0)
    wtl1.plot(ax=ax1)
    bb2.plot(ax=ax2)
    return locals()


@ipynb_doc
def check_nearby(name, radius_cut=5, var_cut=50):
    r""" 
    ### Check {name} for nearby variable sources
    {out}
    
    Display all within {radius_cut} deg of {name} with variability >{var_cut}
    {df}
    """
    global near_df
    with capture_hide('load 4FGL output') as out:
        cat =Fermi4FGL();
    near_df = cat.select_cone(name, query=f'variability>{var_cut}')
    near_df['name'] = near_df.index
    df = near_df['name sep significance variability'.split()].sort_values('sep')
    return locals()


@ipynb_doc
def variability_plots(name, time_bins=(0,0,interval), ylim=None, tsamp=1/24, pmax=None ):
    """
    #### Long-term variability for {name}
    {fig1}
    
    """
    with capture_hide('Setup output') as txt:
        wtl = WtLike(name, time_bins=time_bins) if type(name)==str else name
        name = wtl.source.name
        pgm = wtl.periodogram(tsamp=tsamp)
        pgm.power_spectrum()
        
    def left(ax):
        wtl.bb_view().plot(ax=ax, UTC=True, ylim=ylim)   
    
    def right(ax):
        df = pgm.power_df.query( '5e-4< f<0.003').copy()
        df.loc[:,'period'] = 1/df.f/365.25

        ax.plot(df.period, df.p1, '.-', label='$P_1$')
        ax.plot(df.period, df.p0, '.-', label='$P_0$')
        ax.set(xlabel='Period (yr)', xscale='log', 
            xlim=(1,5),ylabel='Power', ylim=(0,pmax),
            xticks=[1,2,3,4], xticklabels='1 2 3 4'.split())
        ax.grid()
        ax.legend(fontsize=10)
        
    fig1 = plt.figure(figsize=(12,3))
    fig1.subplots_adjust(wspace=0)
    gs = plt.GridSpec(1,10)
    left(  fig1.add_subplot(gs[:6]) )
    right( fig1.add_subplot(gs[7:]) )
    return locals()

def wc_demo(ylim=(0.5,3)):
    wc_intro()
    comment("""### Demonstrate the effect of applying the correction.

#### Before:

""")
    variability_plots(wtl1, ylim=ylim, pmax=50)
    comment( f"""

#### After:

Applied correction to the pulsar {S1} from the nearby variable source {S2}
""")
    wtlx = wtl1.reweighted_view(bb2)
    variability_plots(wtlx, ylim=ylim, pmax=50)
    comment(
    """So both $P_0$ and $P_1$ are now smaller, and the relative size seems to be reversed.
    """)

if __name__ == '__main__':
    wc_demo()