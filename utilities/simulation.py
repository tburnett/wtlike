import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from utilities.ipynb_docgen import *
from wtlike.simulation import Simulation 
from pylib.show_stuff import *

def poiss_fit(w, Nsrc=None):
    """
    Return Poisson fit of weights w, 
    Nsrc -  source events expected or None, in which case use sum(w)
    """
    from wtlike.loglike import LogLike, PoissonRep
    S = np.sum(w) if Nsrc is None else Nsrc
    return PoissonRep(LogLike( dict( n=len(w), w=w, S=S, B=len(w)-Nsrc )))   

def cumsimpson(y, dx):
    """Return the cumulative Simpson integral for an odd number>2 of equally-spaced evaluations
     * y -- array of integrand values, 2n+1 values for 2n intervals
     * dx -- interval size

     Returns array of size n, the cumulative integral estimate for pairs of intervals.

     >>> x, dx = np.linspace(0,1, 5), 1/4
     >>> print(cumsimpson(x, dx))
     >>> [0.125 0.5  ]
    """
    v = np.empty((3*(len(y)-1)//2))
    v[0::3] = y[0:-1:2]
    v[1::3] = 4*y[1::2]
    v[2::3] = y[2::2]
    return (np.cumsum(v)*dx/3)[2::3]

class FunctionGenerator:
    r"""Encapsulate a function or a table to:
    
    * define a continous function via interpolation
    * Implement functions for the intergral and inverse integral 
    * implement a Monte Carlo generation of the functions distribution
    
    """
    def __init__(self, func, limits=(0,1), n=1000,rng=None):
        from wtlike.exposure import cumsimpson
        if isinstance(rng, np.random.Generator):
            self.rng = rng
        else:
            self.rng = np.random.default_rng(rng)
        a,b = self.limits= limits 
        self.n =n
        dx = (b-a)/n
        # evalueate on bin edges
        self.xp = np.linspace(a,b,n+1)
        
        if callable(func):
            self.yp = func(self.xp)
        elif hasattr(func, '__iter__'):
            self.yp = func
        else:
            raise ValueError('Interep: Parameter is not a function or an array')
        if hasattr(func, 'rng'):
            self.rng = func.rng
            
        if np.any(self.yp<0) or np.sum(self.yp)==0:
            raise ValueError('Interp: Function or array is not positive definite')

        # the integral function
        self.cy = np.insert(cumsimpson(self.yp, dx), 0,0)
        self.cx = self.xp[::2]
        norm = self.cy[-1]
        # mean = <w> and quality=<w**2> / <w>
        self.mean = cumsimpson(self.xp*self.yp, dx)[-1]/norm
        self.quality = cumsimpson(self.xp**2 * self.yp, dx)[-1]/norm/self.mean
        
    def __call__(self, x):
        return np.interp(x, self.xp, self.yp)
    def inverse(self, y): # expects monatonic increasing??
        return np.interp(y, self.yp, self.xp)
    def cumulative(self, x):
        return np.interp(x, self.cx, self.cy)
    def inverse_integral(self, x):
        return np.interp(x, self.cy/self.cy[-1], self.cx)
    
    def generate(self, N,):
        from scipy import stats
        return self.inverse_integral(stats.uniform.rvs(size=int(N), random_state=self.rng))
    
    def plot(self, ax=None, figsize=(5,3), **kwargs):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize =figsize) if ax is None else (ax.figure, ax)
        ax.plot(self.xp, self.yp, '-')
        ax.set(**kwargs)
        return fig
    
    @classmethod
    def test(cls, n=100, limits=(0,1), rng=None):
  
        # A test function to digitize
        testfun = lambda x: x**2 
        
        def run_tests(fun):
            f = cls(fun, limits=limits, n=n, rng=rng)
            x,y= f.xp, f.yp
            # and its integral
            F = lambda x: f.cumulative(x)
            test2 = np.max(np.abs(x**3/3 - F(x)))
            assert test2<1e-4, f'Fail test 2: {test2:.2e}>1e-4'
            
            # now inverse
            FI = lambda y: f.inverse_integral(y)
            h = f.generate(1e7)
            test3= abs(h.mean()-3/4) + abs((h**2).mean()-3/5)
            assert test3<1e-3, f'Fail test 3: {test3:.2e}>1e-3'
        
        run_tests(testfun)
        run_tests(testfun(np.linspace(*limits, n+1)))

FunctionGenerator.test() 



class WeightModel:
    r"""### class WeightModel

    Input: signal and background functions  $f_S$ and $f_B$, $\alpha$. 

    Implements 
    * $f(x) = f_B(x) + (1+\alpha)\ f_S(s)$.
    * $w(x) = (1+\alpha)\ f_S(x) / f(x) $

    Here $x$ is an internal variable that is uniformly distributed from 0 to 1.
    """
    def __init__(self, fs, fb,  rng=None):
        """ Manage a weighting model

        fs, fb : functions on [0,1] representing signal and background


        """
        self.alpha = 0
        self.fb=fb
        self.fs=fs
        if isinstance(rng, np.random.Generator):
            self.rng = rng
        else:
            self.rng = np.random.default_rng(rng)        

        # constants corresonding to <w> and 1/sqrt(<w**2>)
        def integral(fun):
            from scipy.integrate import quad
            return quad(fun, 0,1)[0]
        B = integral(self.fb)
        S = integral(self.fs)
        self.norm = T=B+S

        self.source_fraction = W = S/T
        WW = integral(lambda x: self.weight(x)*fs(x))/T
        self.sensitivity = 1/np.sqrt(WW)
        self.variance = WW/W

    def __repr__(self) -> str:
        return f'{self.__class__.__name__}: signal/noise: {100*self.source_fraction:.1f}%, variance/event: {self.variance:.3f}'
 
    def __call__(self, x):
        return  (self.fb(x) + (1+self.alpha)*self.fs(x))/self.norm
    
    def weight(self, x):
        s = self.fs(x) 
        b = self.fb(x)
        return s/(s+b)
    
    def generate(self, Nsrc=10000, alpha=0, time_generator=None):
        """Generate, and return a list of weights
        * Nsrc the number of expected source events to generate if alpha is zero
        * alpha
        * time_generator : None | function that will generate a list poisson  times

        Determine how many, based on alpha and Nsrc, to generate.
        Return that many if time_generator is None, else use it to make a list of random
        times, 0 to 1, and return that and a list of weights with the same length

        Sets self.ngen with the number it tried to generate
        """
        self.alpha=alpha
        self.ngen = Ngen= int(Nsrc * (1/self.source_fraction + alpha))

        if time_generator is None:
            return self.weight( FunctionGenerator(self, rng=self.rng).generate(int(Ngen)))
        
        tt = time_generator(0,1, Ngen, rng=self.rng)
        w = self.weight( FunctionGenerator(self, rng=self.rng).generate(len(tt), ))
        return tt,w
  
    def _make_events(self, start, stop, Nsrc, wmin=1e-5):
        """ For times between start and stop, 
            generate a DF with times and weights """
        # This is pretty complicated
        #
        # evaluate the function to predict the expected distribution of weights
        t = (start+stop)/2
        self.alpha = alpha = self.source(t)-1
        
        # Generate a list of times, number determined by detection efficiency 
        Ngen = int(Nsrc * (1/wtmod.source_fraction + alpha))
        tt = time_generator(start, stop, count=Ngen, rng=self.rng)

        # Now make a list of of weights to correspond with each time
        w = wtmod.weight( FunctionGenerator(wtmod, rng=self.rng).generate(len(tt), ))
        
        df =   pd.DataFrame(dict(time=tt, weight=w.astype(np.float32)))
        return df[w>wmin] if wmin>0 else df
    
    def plot_functions(self, ax=None,  **kwargs):
        fig, ax = plt.subplots(figsize=(4,2)) if ax is None else (ax.figure, ax)
        x = np.linspace(0,1)
        for fn, label,color, lw, in zip((self.fs,self.fb, self),
                             f'$f_S$ $f_B$ total'.split(),
                             'red grey violet'.split(),
                             [1,1,2]):
            ax.plot(x, fn(x), '-', label=label, color=color, lw=lw)
        kw=dict(xlim=(0,1), ylim=(1e-3,None), xlabel='x', title='Functions', yscale='log', )
        kw.update(kwargs)
        ax.set(**kw)
        ax.grid(alpha=0.5)
        ax.legend(fontsize=10)
        return fig
    def plot_weight(self, ax=None, **kwargs):
        fig, ax = plt.subplots(figsize=(3,2)) if ax is None else (ax.figure, ax)
        x = np.linspace(0,1)
        ax.plot(x, self.weight(x), '-')
        ax.set(**kwargs)
        ax.grid(alpha=0.5)
        return fig
    def plot_both(self):
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8,2), sharex=True)
        plt.subplots_adjust(wspace=0.3)
        self.plot_functions(ax=ax1)
        self.plot_weight(ax=ax2, title='weight')
        return fig
    

    @classmethod
    def example(cls, S=0.1, rng=None):
        """
        ### WeightModel.example:

        Returns `WeightModel(fs,fb)` with:
        ```
        fs = lambda x: S * x
        fb = lambda x: (1-x)**2
        ```
        `S` (default 0.1) adjusts the relative size of the nominal signal
        """
        return cls(fs= lambda x: S * x,
                   fb= lambda x: (1-x)**2,
                   rng=rng,
                   )
    @classmethod
    def signal_spec(cls, sn, tau=None, rng=None):
        """
        Specify the signal to noise ratio
        """
        fs = lambda x: x
        fb = lambda x: (1-x)**2 if tau is None else lambda x:  np.exp(-x/tau)
        obj = cls(fs, fb)
        r0 = obj.source_fraction
        a = (1-r0)*sn/(1-sn)/r0
        fsp = lambda x: a*x
        return cls(fsp, fb, rng=rng )
        
    @classmethod
    def signal_tau(cls, sn, tau=0.1, rng=None):
        """
        Specify the signal to noise ratio
        """
        fs = lambda x: 1 if np.isscalar(x) else np.ones(len(x))
        fb =  lambda x:  np.exp(-x/tau)
        obj = cls(fs, fb)
        r0 = obj.source_fraction
        a = (1-r0)*sn/(1-sn)/r0
        fsp = lambda x: a* fs(x)
        return cls(fsp, fb, rng=rng )
    @classmethod
    def test(cls):
        wm = cls.example()
        wgen = wm.generate(1e7)
        W = wgen.mean(); WW = np.mean(wgen**2)
        test1 = np.abs(W/wm.source_fraction -1)
        test2 = 0 #np.abs(WW/W /wm.sensitivity -1)
        # show(pd.Series(dict(W_meas=W,  
        #             W_int = wm.source_fraction, #integral(weight),
        #             WW_meas=WW, 
        #             WW_int= wm.sensitivity*wm.source_fraction,
        #             test1=test1,
        #             test2=test2,
        #            ),
        #        name='Compare integrals with MC sums'))
        assert test1<2e-3 and test2<1e-3, f'WeightModel failed a test :{test1:.2e}, {test2:.2e}'

WeightModel.test()        
  

class PeriodicFunc:

    def __init__(self, flux):
        """ Superclass saves absolute flux"""
        self.flux = flux

    def __call__(self, t):
        """ relative flux is 1"""
        if np.isscalar(t): return 1
        return np.ones(len(t))
    
    def __repr__(self):
        return self.__class__.__name__
    
    @property
    def name(self): return self.__class__.__name__
    
    def info(self):
        """ return dict of public parameters"""
        d = self.__dict__
        return dict( (key,d[key]) for key in d.keys() if not key.startswith('_') )

    def plot(self, ax=None,  **kwargs):
        """
        return a plot of the function
        """
        kw = dict( xlabel='t', xlim=(0, self.period), ylabel=self.latex, ylim=(0,None))
        kw.update(kwargs)
        t = np.arange(*(kw['xlim']), step=1e-2*self.period)    

        fig,ax = plt.subplots(figsize=(3,1.5)) if ax is None else (ax.figure,ax)
        ax.plot(t,self(t),color='orange', lw=2)
        ax.grid(alpha=0.5)
        ax.set(**kw)
        return fig

    @classmethod
    def demo(cls, tlim=(0,1), **kwargs):

        cf = cls(**kwargs)

        show(f"""## {cf.__class__.__name__} demo/test""")
        show(cls.__doc__)
        vars = cf.__dict__
        hidden = np.array( list(
            map(lambda n: n.startswith('_'), vars.keys()) 
        ))
        pars = pd.Series(vars, name='value')
        with pd.option_context('display.precision', 3):
            show(pd.Series(pars[~hidden], summary='Parameters'))

        show(cf.plot())


class CosFunc(PeriodicFunc):
    """Manage a cosinusoidal signal modulation
    """
    def __init__(self, flux=1, A=1, period=1, tzero=0):
        super().__init__(flux)
        assert A>=0 and  A<=1, 'limit amplitude modulation'
        self.A = A
        self.period = period
        self._w = 2*np.pi/period
        self.tzero = tzero
  
    def __call__(self, t):
        # t = np.atleast_1d(t)
        return 1+ self.A*np.cos(self._w*(t-self.tzero ) )
    @property
    def latex(self):
        return rf'$1 + {self.A:.1f}\ \cos(2\pi f_S t)$'

        
class HatFunc(PeriodicFunc):
    """
    A periodic step function
    """
    def __init__(self, flux=1, high=1.5, low=0, period=1,  tzero=0):
        super().__init__(flux)
        assert high>1 and low>=0 and low<1, 'Required'
        self._dc= period*(1-low)/(high-low)
  
        self.period = period
        self.tzero = tzero
        self.high = high
        self.low = low

    @property
    def latex(self):
        return 'step function'
    def __call__(self, t):
        # t = np.atleast_1d(t)
        tau = np.mod(t-self.tzero, self.period)
        return np.where(tau<self._dc, self.high, self.low)

sec_per_day = 86400
def time_generator(start, stop, count, rng=None):
    """ Generate a list of times, distributed randomly

    - start, stop: times
    - count : expected number to generate with rate=count/(stop-start)
    - rng : random state generator

    returns : list of times between start and stop. Note that the actual number is Poisson-distributed
    """
    scale = (stop-start)/count
    if rng is None: rng = np.random.default_rng()
    elif not isinstance(rng, np.random.Generator):
        rng = np.random.default_rng(rng)
    ## equivalent to below, but **much** slower!
    # tt =[]
    # t = start
    # while True:
    #     t += rng.exponential(scale =scale)
    #     if t>stop: break
    #     tt.append(t)
    ##
    # Generate 4 sigma more than needed, insert start at front, do cumsum, truncate at stop
    gensize = int(count+4*np.sqrt(count))
    tx = np.insert(rng.exponential(scale=scale, size=gensize),0,start)
    tt = np.cumsum( tx ).astype(np.float32)
    
    return tt[tt<stop]

def make_exposure(fexp, start, stop, interval=300/sec_per_day):
    """
    - fexp -- exposure in cm^2, a value or a function of time in day units
    - start, stop -- range of time in day units
    - interval [300] --(day units) default 5-min interval (fermi data is 30 s)

    Returns: a DataFrame with start, stop, exp, latter in cm^2 s
    """
    if np.isscalar(fexp):
        fval = fexp
        fexp = lambda t: fval
    
    nbins = int( (stop-start) / interval )
    edges = np.linspace(start, start+nbins*interval, nbins+1)
    starts, stops = edges[:-1], edges[1:]
    exp = fexp(starts) * interval *sec_per_day
    
    return pd.DataFrame.from_dict(
                dict(start=starts, stop=stops,  exp=exp,)
                                )


class WtSim(Simulation):
    """Simulate weighted photon data

    - name : required name
    - exposure : DataFrame-like object with columns (start, stop, exp)
    - source : function of time (d) returning flux in cm^-2 s^-1. 
            has a property `flux` with nominal value for evaluation of relative flux
    - weight_model : instance of WeightModel which manages weight generation
    - rng : random generator instance, or integer seed

    """
    def __init__(self, 
                 name, 
                 exposure,
                 source,
                 weight_model, 
                 rng=None,
                 time_generator=None):

        self.name = name
        self.source = source
        self.time_generator = time_generator

        self.exposure = exposure #make_exposure(efun, tstart, tstop, interval=86400)
        self.filename = None #flag that not a regular source
    
        if isinstance(rng, np.random.Generator):
            self.rng = rng
        else:
            self.rng = np.random.default_rng(rng)
        self.weight_model = weight_model
        self.weight_model.rng = self.rng

    def __repr__(self):
        exp = self.exposure
        return f"""Simulation {self.name}, {len(exp)} intervals from {exp.start[0]} to {exp.stop.values[-1]}
            """

    def _make_events(self, start, stop, Nsrc, wmin=1e-5):
        """ For times between start and stop, 
            generate a DF with times and weights """

        t = (start+stop)/2
        self.alpha = alpha = self.source(t)-1
    
        tt,w = self.weight_model.generate(Nsrc, alpha, time_generator)
        tt = start + tt * (stop-start)
        
        df =   pd.DataFrame(dict(time=tt, weight=w.astype(np.float32)))
        return df[w>wmin] if wmin>0 else df

    def run(self):
        """ Create dataframe `self.photons`, using the exposure and source, and the weight model for the latter
        """
        simdata = []
        deltat = 0 # the time from the last event to end of the (start,stop) interval
        for start, stop, exp in self.exposure.itertuples(index=False, name=None):
            Nsrc = exp * self.source.flux
            events = self._make_events( start-deltat, stop, Nsrc, ) 
            simdata.append(events )
            deltat = stop
            if len(events)>0: deltat -= events.time.values[-1]
        
        self.photons =  pd.concat(simdata)

 
def show_weight_stat(weights, text=''):

    N = len(weights)
    W = np.sum(weights)
    U = np.sum(weights**2)
    show(f"""{text}
    |  N  | mean | quality|
    |-----|--------|-------|
    | {N:,d} | {W/N:.2e} | {U/W:.3f}|
    """)

def show_hist(x, bins=100, hkw={}, **kwargs):
    fig, ax =plt.subplots(figsize=kwargs.pop('figsize', (4,2)))
    summary = kwargs.pop('summary',None)
    ax.set(**kwargs)
    hkwt = dict(histtype='step', lw=2, density=True)
    hkwt.update(hkw)
    h = ax.hist(x, bins, **hkwt)
    show(fig, summary=summary)

def weight_analysis(weights, name, nbins=100):
    bins = np.linspace(1/nbins,1,nbins) #skip first bin
    fig, ax =plt.subplots(figsize=(4,2))
    whist,_,_= ax.hist(weights.clip(0,10), bins=bins, log=False,histtype='step',
                       density=False, lw=2)

    ax.set(xscale='linear', xlabel='weight', title=name)
    ax.grid(alpha=0.5)
    show(fig)
    return whist

 
class RunSim:
    def __init__(self, src_info, exp_info, tsamp=1/36, phase_bins=8):
        self.exp_info = exp_info
        self.periodic = src_info.pop('func')(**src_info)
        self.exposure = make_exposure(**exp_info) 
        show(f"""---
            ## Run {self.periodic} simulation 
            """)

        sim = WtSim(f'{self.periodic} Simulation', 
                  exposure=self.exposure,
                  source = self.periodic,
                  weight_model=WeightModel.example(S=0.005),
                  rng=42) 
        with capture('Analysis output') as self.output:
            with Timer() as elapsed:
                self.wtl = WtLike(sim, time_bins=(1,0,1))
                self.px = self.wtl.periodogram(tsamp=tsamp)
                self.phase_view = self.wtl.phase_view(self.periodic.period, phase_bins)
            print(elapsed)
        show(self.output)
        
    def show_all(self):
        show_dict(self.exp_info,  summary='Exposure parameters',show=True) 

        show(self.exposure.head(), summary=f"""Exposure ( {len(self.exposure)} intervals from 
             {self.exp_info['start']} to {self.exp_info['stop']}.)""")

        show_dict(self.periodic.info(), summary=self.periodic.name + ' parameters')
        show_source_flux(self.exposure, self.periodic)
        show("""### WtLike analysis""")

        show(self.wtl.plot(), figsize=(15,3), summary='light curve')
        show(self.wtl.fluxes.head(), summary='Flux df head()' )
        
    def show_phase_plot(self, xlim=(0,1)):
        show('### Phase light curve')
        fig, ax = plt.subplots(figsize=(3,2))
        self.phase_view.plot(ax, xlim=xlim,) 
        self.periodic.plot(ax=ax)
        show(fig)

    def show_phase_data(self):
        df = self.phase_view.fluxes
        show(df.iloc[:len(df)//2], summary=f'Phase view data')
    
    def show_power_plot(self, xlim=(0,9), ylim=(-10,None)):

        self.px.power_plot(pmax=50, xlim=xlim, ylim=ylim)
        show(plt.gcf(), figsize=(8,2))
        show_fft_peaks(self.px);     
