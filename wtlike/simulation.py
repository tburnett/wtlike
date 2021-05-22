# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/04_simulation.ipynb (unless otherwise specified).

__all__ = ['generate_times', 'sec_per_day', 'WeightFunction', 'make_exposure', 'Simulation']

# Cell
import os
import numpy as np
import pandas as pd
from scipy import stats

from .config import Config, PointSource


# Cell
import numbers

class _Sampler():
    """ Sample an arbitrary function or histogram

    - func -- a function, a histogram, or a fixed value<br>
        If a function, must be positive definite.<br>
        Assume histogram bins are 0 to 1.
    - a,b  -- limits (default 0,1)
    - n    -- table size (ignored if a histogram or value)

    """

    def __init__(self, func, limits=(0,1), n=100):

        a,b = limits
        self.x = np.linspace(a,b,n+1) # bin edges
        dx = (b-a)/(n)/2
        self.deltafun=None

        if callable(func):
            # A function
            # evaluate at bin centers
            y = np.array([func(t-dx) for t in self.x])
            if np.any(y<0) or np.sum(y)==0:
                raise ValueError('Function is not positive definite')
        elif isinstance(func, numbers.Number):
            # a single value, or delta function
            self.deltafun = func
            if  func<0 or func>1:
                raise ValueError('Value not in range [0,1]')
            self.mean=func
            return
        else:
            n = len(func)
            self.x = np.linspace(a,b,n)
            y = func
        cy = np.cumsum(y)
        d = cy[-1]-cy[0]
        self.sy = (cy-cy[0])/d

        self.mean = np.sum( (self.x-dx) * y) / d

    def _evaluate(self, r):
        """evaluate inverse integral. expect 0<r<1 """
        return np.interp(r, self.sy, self.x)

    def __call__(self, size):
        """Generate `size` values
        """
        if self.deltafun: return np.full(size, self.deltafun)

        return self._evaluate(stats.uniform.rvs(size=size))

# Cell
sec_per_day = 24*3600

def generate_times(start, stop, count):
    """ Generate a list of times, distributed randomly

    - start, stop: times
    - count : expected number to generate with rate=count/(stop-start)

    returns : list of times between start and stop. Note that the actual number is Poisson-distributed
    """
    # note: can speed this up by making groups of random calls

    tt =[]
    t = start
    scale = (stop-start)/count
    while True:
        t += np.random.exponential(scale =scale)
        if t>stop: break
        tt.append(t)
    return tt

# Cell
class WeightFunction(object):

    def __init__(self, s=1,b=1, wt_signif=0.1):
        self.s = s
        self.b = b
        self.lam = wt_signif

    def __call__(self, r):
        return (self.s * np.exp(-r/self.lam)/(self.lam*(1-np.exp(-1/self.lam))) + self.b)

    def sample(self, s,b, n):
        self.s = s
        self.b = b
        return _Sampler(self, n=1000)(n);

    def weights(self, s, b, n):
        h = self.sample(s,b,n)
        return 1-b/self(h)

# Cell
def make_exposure(fexp, start, stop, interval=300):
    """
    - fexp -- exposure in cm^2, a value or a function of time in day units
    - start, stop -- range of time in day units
    - interval [300] -- 5-min interval (fermi data is 30 s)

    Returns: a DataFrame with start, stop, exp
    """
    def check_scalar( f):
        if np.isscalar(f):
            fval = f
            return lambda t: fval
        return f
    fexp = check_scalar(fexp)

    nbins = int((stop-start)*sec_per_day / interval)
    edges = np.linspace(start, start+nbins*interval/sec_per_day, nbins+1)
    starts, stops = edges[:-1], edges[1:]
    exp = fexp(starts) * interval
    return pd.DataFrame.from_dict(dict(start=starts, stop=stops, exp=exp))

# exp  = make_exposure(500, 0, 1 )
# days  = np.sum(exp.stop-exp.start); secs = days*24*3600
# exptot=np.sum(exp.exp)
# exp_text = f' average {exptot/secs:.0f} cm^2 for {secs/1e6:.1f} Ms'
# print(exp_text)
# exp.head()

# Cell
class Simulation(object):

    def __init__(self, name, src_flux, tstart, tstop, bkg_rate=1e-6,  efun=3000, wt_signif=0.1):
        """
        - src_flux : source flux, scalar or function of days, typically around 1e-7
        - tstart, tstop :(days)
        - bkg_rate : background flux, scalar or function of day, typicaly 1e-6 for 4-deg cone
        - efun : scalar, function (of time in days) of the exposure/s. Typically 3000 cm^2 for fermi

        - wt_signif : now the width of the PSF in (r/rmax)**2 coordinates

        """
        def check_scalar( f):
            if np.isscalar(f):
                fval = f
                return lambda t: fval
            return f
        self.name = name
        self.src_fun = check_scalar(src_flux)
        self.bkg_fun = check_scalar(bkg_rate)
        self.flux_fun = lambda t: src_fun(t)+bkg_fun(t)
        self.wt_signif=wt_signif

        self.exposure = make_exposure(efun, tstart, tstop)


    def run(self):
        times = []
        weights = []
        for start, stop, exp in self.exposure.itertuples(index=False,name=None):

            src = self.src_fun((start+stop)/2)
            bkg = self.bkg_fun((start+stop)/2)
            delta_t = (stop-start)*sec_per_day # tolal tim
            counts = (src+bkg) * exp #
            #print(f'From {start} to {stop}, exposure/s {exp/delta_t:.0f}, counts {counts:.0f}')
            new_times = generate_times(start, stop, counts)
            wfun = WeightFunction(wt_signif=self.wt_signif)
            new_wts = wfun.weights(s=src, b=bkg, n=len(new_times));

            assert len(new_times)==len(new_wts)
            times = np.append(times, new_times)
            weights = np.append(weights, new_wts)

        print(f'generated {len(times)} photons')
        self.photons=pd.DataFrame(dict(time=times, weight=weights.astype(np.float32)))