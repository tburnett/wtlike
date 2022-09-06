from wtlike import *
from utilities.ipynb_docgen import *

def lowfreqplot(pgm, ax=None, over=None, 
                pticks=None,
                query='f<0.0208',**kwargs):
    
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
    
    if over is not None:
        a,b = over.lim()
        xx = x[(x>a) & (x<b)]
        ax.plot(xx, over(xx), '--r', label=f'sinc at {1/(over.freq):.2f} yr')

    if pticks is None: return
    a,b = np.array(ax.get_xlim())
    x2 = lambda p: (1/p-a)/(b-a)
    ax.twiny().set(xlabel='Period (yr)',xlim=(0,1),
            xticks=x2(np.array(pticks)), 
            xticklabels=[ f'{t}' for t in pticks])  
    ax.legend(fontsize=12) 
    
def simple_power_plot(pgm, ax=None, query=None, power='p1', pticks=None, **kwargs):
    fig, ax = plt.subplots(figsize=(5,3)) if ax is None else (ax.figure, ax)
    yr = 365.25
    df = pgm.power_df #.query(query).copy()
    # fig2, ax = plt.subplots(figsize=(4,3))
    x = df.f #*yr
    ax.plot(x, df[power], '.-', label='periodogram', color='blue')
    kw = dict(xlabel=f'Frequency (1/d)', ylabel=f'Power ${power}$')
    kw.update(kwargs)
    ax.set(**kw)
    ax.grid()
    
    if not pticks: return
    a,b = np.array(ax.get_xlim())
    x2 = lambda p: (1/p-a)/(b-a)
    ax.twiny().set(xlabel='Period (d)',xlim=(0,1),
            xticks=x2(np.array(pticks)), 
            xticklabels=[ f'{t}' for t in pticks])  
    

@ipynb_doc
def title():
    """

    # Study the binaries XSS J12270-4859 and 1FGL J1018.6-5856
    
    <h5 align="right">{date}</h5>
    
    The first is [PSR J1227-4853](https://arxiv.org/pdf/1412.4735.pdf), a 1.69 ms MSP, 
    the original black widow, orbital period 0.28788 d and the second is a 
    [Fermi-discovered](https://ntrs.nasa.gov/citations/20120016459) gamma-ray binary 
    with a 16.6-day period.
    """

    return locals()

@ipynb_doc
def j1227(name='PSR J1227-4853'):
    """
    ## {name}

    
    ### Bayesian-Block Light curve
    {fig1}
    The BB analysis identiview three blocks.
    
    ### Periodogram

    {fig2}
    This is featureless except for the low frequency reflection of the transitions.    
    
    {fig4}
    
    ### Look at the expected orbital frequency:
    {fig3}
    
    There is a certainly a peak consistent with that frequency, but at a level that is not
    significant.
    {printout}
    """
    
    with capture_hide('printout') as printout:
        wtlx = WtLike(name) #XSS 12270-4859')

        fig1, ax = plt.subplots(figsize=(12,4))
        wtlx.view(30).bb_view().plot(ax=ax);

        px = wtlx.periodogram()
        fig2, ax = plt.subplots(figsize=(12,3))
        px.power_plot(ax=ax, pmax=50)

        orbital_freq = 1/0.287887

        fig3,ax = plt.subplots(figsize=(6,3))
        simple_power_plot(px, ax=ax, xlim=(3.4724, 3.475), power='p0', pticks=[0.290,0.292, 0.294, 0.296], ylim=(0,20))

        # px.power_plot(ax=ax, xlim=(orbital_freq-0.005,orbital_freq+0.005),ylim=(0,30));
        ax.axvline(orbital_freq, ls='--', color='orange', label='known orbital period')
        # ax.set(title='1FGL J1018.6-5856');
        ax.legend();

        # px.power_plot(xlim=(0,0.005))

        fig4, ax = plt.subplots(figsize=(6,3))
        lowfreqplot(px, ax=ax, yscale='log', ylim=(1,None), xlim=(0,2),pticks=[0.3,0.4,0.5,1,1.5,2,4],);
        
    return locals()

@ipynb_doc
def j1018(name='1FGL J1018.6-5856'):
    r"""## {name}
    
    ###  Bayesian-Block Light curve and periodogram:
    {fig1}
    {fig2}
    
    Except for a fluctuation that should be looked at, the lightcurve is featureless. The periodogram shows
    the expected orbital frequency and three harmoiics
    
    ### Blow-up around orbital frequency:
    {fig3}
    The peak corresponds to a period of ${period:.2f} \pm 0.02\ \mathrm{d}$.
    
    {out}
    """
    with capture_hide('printout') as out:
        wtly = WtLike(name)

        fig1, ax = plt.subplots(figsize=(12,4)) 
        wtly.view(14).bb_view().plot(ax=ax);

        py = wtly.periodogram()
        fig2, ax = plt.subplots(figsize=(12,3))
        py.power_plot(ax=ax,xlim=(0,0.5));
        dfp = py.find_peaks().query('p0>200'); 
        period = 1/dfp.f.values[0]
        print(f'Source period: {period:.2f} d')
        fig3, ax = plt.subplots(figsize=(6,3))


        simple_power_plot(py, ax=ax, xlim=(0.0595,0.0615), pticks=[16.4,16.5, 16.6])
    
    return locals()

def main():
    title()
    j1227()
    j1018()
    