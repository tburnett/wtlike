# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/06_poisson.ipynb (unless otherwise specified).

__all__ = ['Poisson', 'PoissonFitter']

# Cell
import sys
import numpy as np
from numpy import polyfit
from scipy import optimize, special, stats

class Poisson(object):
    r"""This is a functor class which returns the log of a three-parameter Poisson-like function used to represent
    the flux likelihood. The parameters are, in order:

    $s_p$ : flux at peak, if positive; if negative, there is only a limit

    $e$ : normalization factor to convert flux to equivalent counts. must be >0

    $b$  : background flux: must be >=0

    This parameterization is equivalent to that described in William Tompkins' thesis
    (arxiv: astro-ph/0202141) and Nolan, et al., 2003,  ApJ 597:615:627.
    The functional form is that of a Poisson distribution.

    The log likelihood expression is for flux $s>=0$
       $$w(s | s_p,e,b) = e\ (s_p+b) \ \log[ e\ (s+b) ] - e\ s + \mathrm{const}$$
    the const is such that $w=0$ at the peak.

    A slightly more elegant expression, used in the cdf function, is to define
       $\beta = e b$,  $\mu = e s_p + \beta$, and  $x = e s$.

    Then the log likelihood is

       $$ f(x | \beta) =  \mu \ \log( x + \beta) - x + \mathrm{const}$$

    where the peak is at $x=x_p=\max(0, \mu-\beta)$, and the constant is defined so that $f(x_p)=0$.

    """

    def __init__(self, p:'array of parameters'):
        """p : array of parameters
        """
        self.p = p

    def __call__(self, dom):
        """Return the value(s) of the fit function for the given domain."""
        sp,e,b = self.p
        b = abs(b) # this is a bit of a kluge to keep log argument positive
        if b==0: b = 1e-20 #another kludge
        e = abs(e) #and another
        r = e*(dom+b)
        r_peak = e*(sp+b)
        if sp > 0:
            const = r_peak*np.log(r_peak) - r_peak
        else:
            #sp=0
            #const = r_peak*np.log(r_peak) - r_peak
            t = e*b
            const = r_peak*np.log(t) - t
        f = r_peak*np.log(r) - r
        return f - const

    def __str__(self):

        e, beta, mu = self.altpars()
        return f'Poisson: mu,beta= {mu:.1f}, {beta:.1f}'

    def __repr__(self):
        if self.flux==0:
            return 'flux is zero for source'
        t = np.array(self.errors)/self.flux-1
        relerr = np.abs(np.array(self.errors)/self.flux-1)
        return f'{self.__module__}.{self.__class__.__name__}: {self.flux:.3f}[1+{relerr[0]:.3f}-{relerr[1]:.3f}]'

    @property
    def flux(self):
        return max(self.p[0], 0)

    @property
    def errors(self):
        return self.find_delta()

    @property
    def ts(self):
        return 0 if self.flux<=0 else (self(self.flux)-self(0))*2.0


    def altpars(self):
        """ return alternate parameters: e, beta, mu """
        e = abs(self.p[1])
        beta = e * abs(self.p[2])
        mu =   e * self.p[0] + beta
        return e,beta,mu

    def find_delta(self,delta_logl=.5):
        """Find points where the function decreases by delta from the max"""
        smax = max(0,self.p[0])
        ll_max = self(smax)
        ll_zero = self(0)
        func = lambda s: ll_max-self(s)-delta_logl
        if ll_max-ll_zero<delta_logl:
            s_low = 0
        else:
            #s_low = optimize.bisect(func,0,smax,xtol=.01*smax)
            s_low = optimize.brentq(func,0,smax,xtol=1e-17)
        if smax>0:
            s_high = smax*10
        else:
            s_high = 1e-15
        while func(s_high)<0: s_high*=2
        #s_high = optimize.bisect(func,smax,s_high,xtol=.01*smax)
        s_high = optimize.brentq(func,smax,s_high,xtol=1e-17)
        if not np.all(np.isreal([s_low,s_high])):
            print('Could not find two roots!')
            return None
        return (s_low,s_high)

    def cdf(self, flux ):
        """ cumulative pdf, from flux=0, according to Bayes
        uses incomplete gamma function for integrals. (Note that the scipy function
        is a regularized incomplete gamma function)
        """
        e, beta, mu = self.altpars()
        offset = special.gammainc(mu+1, beta) # Bayes offset if beta>0
        return (special.gammainc(mu+1, beta+flux*e)-offset)/(1-offset)

    def cdfc(self,flux):
        """ complementary cumulative cdf: 1-cdf(flux)"""
        e, beta, mu = self.altpars()
        return special.gammaincc(mu+1, beta+flux*e)/special.gammaincc(mu+1, beta)

    def cdfinv(self, pval):
        """ return the inverse of the cdf
         pval : float
        """
        e,beta,mu = self.altpars()
        gbar = lambda x : special.gammainc(mu+1, beta+x)
        chatinv = lambda pv : special.gammaincinv(mu+1, pv+gbar(0)*(1-pv))-beta
        return  chatinv(pval)/e

    def cdfcinv(self, pvalc):
        """ return the inverse of cdfc = 1-cdf
        useful for limits with very low probablity
        pvalc : float
            1-pval: zero corresponds to infinite flux
        """
        e,beta,mu = self.altpars()
        if e==0: return np.nan
        gcbar = lambda x : special.gammaincc(mu+1, beta+x)
        #cchat = lambda x : gcbar(x)/gcbar(0)
        cchatinv = lambda pv : special.gammainccinv( mu+1, pv*gcbar(0) )-beta
        return cchatinv(pvalc)/e

    def percentile(self, limit=0.95):
        """Left for compatibility: use cdfinv or cdfcinv
        """
        if limit>=0.95: return self.cdfcinv(1-limit)
        return self.cdfinv(limit)

    def pts(self):
        return 0 if self.flux<=0 else (self(self.flux)-self(0))*2.0

    def zero_fraction(self):
        """ Return an estimate of the fraction of the probability that corresponds
        to negative flux. Assume only have to calculate if TS>0 or TS<16
        """
        if self.ts==0: return 1.0
        if self.ts>16: return 0.0
        #this assumes that self.cdfc(0) is 1.0, and that cdfc(-b) is the full integral
        return 1-1./self.cdfc(-self.p[2])

    def sigma_dev(self, val):
        """Express the deviation of `val` from compatibility with this PDF using a sigma measure"""
        # norm(0,1).isf is inverse survival probability for Gaussian
        return stats.norm(0,1).isf(self.cdf(val))

    @classmethod
    def from_fit(cls, counts, flux, sig_flux, tol=5):
        """
        Create a Poisson instance using the fit parameters from an analysis
        of a likelihood function.

        - counts -- number of weights
        - flux -- peak value, in flux units
        - sig_flux -- its standard deviation, square root of the variance

        This is relevant if the likelihood is not truncated by the Bayes requirement.

        If $s$ and $v$ are the signal rate and its variance, determined by an optimization for $n$ measurements,
        so $s$ is large compared with $\sqrt{v}$, then the poisson-like parameters can be determined as follows.

        We use the properties of the Poisson distribution
        $f(n;\lambda) = \exp(n \log\lambda - \lambda + \mathrm{const})$,
        where the parameter $\lambda$ is equal to the expected value of number of occurrences $n$ and
        to its variance, and that the function we want is shifted by the background $b$ and scaled by a factor
        $k$ so we use $f(k(s-b)| \lambda)$. This implies that for the expected value of $s$, $\lambda = n$,
        and $ k(s-b)= k^2 v = n$.

        """
        mu = counts
        beta =  counts-np.sqrt(counts)/(sig_flux/flux)
        assert beta>=0, f'std invalid, too small for number of counts'
        k = (mu-beta) /flux
        b = beta/k
        P = cls([flux, k , b])
        # check values at 0 and the peaks
        if P(0) > -tol**2/2:
            raise(f'Bayes violation: P(0)= {P(0):.2f} too large for'\
                f' {tol} sigma limit. ')

        assert np.abs(P(flux))<1e-3, f'Fail validity: peak {P(mean)}'
        return P

# Cell
class PoissonFitter(object):
    """ Helper class to fit a log likelihood function to the Poisson

    parameters
    ----------
    - func : function of one parameter
    - fmax : position of maximum value, or None
           if None, estimate using fmin
    - scale: float | None
        estimate for scale to use; if None, estimate from derivatime
    - tol : float
        absolute tolerance in probability amplitude for fit, within default domain out to delta L of 4
    - delta : float
        value to calculate numerical derivative at zero flux

    """

    def __init__(self, func, fmax=None, scale=None,  tol=0.20, delta=1e-4, dd=-0.1, test_mode=False):
        """
        parameters
        ----------
        func : function of one parameter
        fmax : position of maximum value, or None
               if None, estimate using fmin
        scale: float | None
            estimate for scale to use; if None, estimate from derivatime
        tol : float
            absolute tolerance in probability amplitude for fit, within default domain out to delta L of 4
        delta : float
            value to calculate numerical derivative at zero flux
        """
        self.func = func
        #first check derivative at zero flux - delta is sensitive
        self.f0 = func(0)
        s = self.wprime = (func(delta)-self.f0)/delta
        if scale is None: scale= 5/s if s>0 else 1.0
        self.smax = fmax or self.find_max(scale) if s>=0 else 0.

        self.ts = 2.*(func(self.smax) - self.f0)

        if test_mode:
            return
        # determine values of the function corresponding to delta L of 0.5, 1, 2, 4
        # depending on how peaked the function is, this will be from 5 to 8
        # The Poisson will be fit to this set of values
        dlist = np.array([0.5, 1.0, 2.0, 4.0])
        if s < dd:
            # large negative derivative: this will be just an exponential
            if s < -100: s=-100. #cut off for now
            self.dom = - dlist/s
            self._poiss= Poisson([-1, -s, 1])
            self.maxdev=0
            return #no test in this case
        else:
            dom = set()
            for delta in dlist:
                a,b = self.find_delta(delta, scale, xtol=tol*1e-2)
                dom.add(a); dom.add(b)
            self.dom = np.array(sorted(list(dom)))
            self.fit()
        self.maxdev=self.check(tol)[0]

    def __repr__(self):
        return '%s.%s : wprime=%.3e maxdev=%.2f, %s' % (self.__module__,self.__class__.__name__,
            self.wprime, self.maxdev,  str(self._poiss))

    @property
    def poiss(self):
        return self._poiss

    def __call__(self, x):
        if hasattr(x, '__iter__'):
            ret =list(map(self.func, x))
            return ret[0] if len(ret)==1 else ret
        return self.func(x)

    def find_max(self, scale):
        """Return the flux value that maximizes the likelihood.
        """
        # if self.func(0) > self.func(scale/10.) and self.wprime<0:
        #     return 0
        r= optimize.fmin(lambda s: -self.func(s), scale, ftol=0.01, xtol=0.01,
                disp=False, full_output=True, retall=True)
        t = r[0][0]
        #if t==scale:
        #    raise Exception('Failure to find max value: %s' % list(r))
        return t if t>0 else 0

    def find_delta(self, delta_logl=.5, scale=1.0, xtol=1e-5):
        """Find positive points where the function decreases by delta from the max
        """
        ll_max = self(self.smax)
        ll_zero = self(0)
        func = lambda s: ll_max-self(s)-delta_logl
        if ll_max-ll_zero<delta_logl:
            s_low = 0
        else:
            s_low = optimize.brentq(func,0, self.smax, xtol=xtol)
        if self.smax>0:
            s_high = self.smax*10
        else:
            s_high = scale
        while func(s_high)<0 and s_high<1e6:
            s_high*=2
        s_high = optimize.brentq(func,self.smax,s_high, xtol=xtol)
        if not np.all(np.isreal([s_low,s_high])):
            msg = '%s.find_delta Failed to find two roots!' % self.__class__.__name__
            print(msg)
            raise Exception( msg)
        if s_high==s_low:
            msg= '%s.find_delta Failed to find high root with delta=%.1f: %s' % (self.__class__.__name__,delta_logl,s_high)
            print(msg)
            print('wprime: %.3e' % self.wprime)
            raise Exception(msg)
        return (s_low,s_high)

    def fit(self, mu=30, beta=5):
        """Do the fit, return parameters for a Poisson constructor
        mu, beta: initial parameters for the fit if the peak is positive
        """
        smax = self.smax
        if smax>0:
            # function to fit has positive peak. Fit the drived parameters mu, beta
            cod = self(self.dom)-self.func(smax)
            #print 'smax=%.2f, w(smax)=%s, w(%s)=%s' % (smax,self.func(smax), self.dom, cod)
            def fitfunc(p):
                mu,beta=p
                e=(mu-beta)/smax; b = beta/e
                self._poiss = Poisson([smax, e,b])
                r = self._poiss(self.dom)-cod
                #print'f(%.3f,%.3f): %s' % (mu,beta,r)
                return r
            mu,beta =  optimize.leastsq(fitfunc,[mu,beta], ftol=1e-6,xtol=1e-6, maxfev=10000)[0]

            e = (mu-beta)/smax; b = beta/e
            return [smax, e, b]
        else:
            # maximum is at zero, so only limit.
            x=self.dom; y=self(x)
            # exposure factor estimated from asymptotic behavior
            big= x[-1]*1e3; e = -self(big)/big;
            # preliminary fit to the quadratic coeficients and estimate parameters from the linear and 2nd order
            pf = polyfit(x,y, 2)
            b,a = pf[:2]
            beta = -e*(e+a)/(2.*b)
            mu = beta*(1+a/e)
            smax= (mu-beta)/e; b= beta/e
            # now fit the Poisson with e fixed to the asym. estimate.
            cod = self(self.dom)-self(0)
            pinit = [smax,  b]
            def fitfunc(p):
                self._poiss = Poisson([p[0], e, p[1]])
                return self._poiss(x) -cod
            t =  optimize.leastsq(fitfunc,  pinit,  xtol=1e-6,ftol=1e-6, maxfev=10000)[0]
            return self._poiss.p

    def check(self, tol=0.05):
        offset = self(self.smax)
        dom =self.dom[1:] # ignore first on
        deltas = np.array([np.exp(self.func(x)-offset)-np.exp(self._poiss(x)) for x in dom])
        t = np.abs(deltas).max()
        if t>tol:
            print(f'PoissonFitter warning: max dev= {t:.3f} > tol= {tol}. (wprime={self.wprime:.2f})', file=sys.stderr )
            #raise Exception(f'PoissonFitter: max dev= {t:.3f} > tol= {tol}. (wprime={self.wprime:.2f})' )
        return t, deltas

    def plot(self, ax=None, xticks=True, legend=True ):
        """Return a figure showing the fit"""
        import matplotlib.pyplot as plt
        xp = self.dom
        x = np.linspace(0, xp[-1]*1.05, 25)
        if ax is None:
            fig, ax = plt.subplots(figsize=(3,3))
        else: fig = ax.figure
        pfmax = self(self.smax)
        ax.plot(x, np.exp(self(x)-pfmax), '-', label='Input function')
        ax.plot(xp, np.exp(self(xp)-pfmax), 'o', color='orange', label='points used for approx.')
        ax.plot(x, np.exp(self._poiss(x)), '--', lw=2,color='orange', label='Poisson approx.')
        if legend: ax.legend(loc='upper right', prop=dict(size=8) )
        ax.set(xlim=(0,None), ylim=(0,1.05))
        if xticks:
            ax.set_xticks([0, xp[-1]])
        ax.grid(alpha=0.4)
        fig.set_facecolor('white')
        return fig

    def normalization_summary(self, nominal=None):
        """return a dict with useful stuff for normalization check
            nominal: None or float
                if specified, calculate delta_ts and pull
        """
        poiss = self.poiss
        lower, upper = poiss.errors
        maxl = poiss.flux
        err = self.maxdev
        if nominal is not None:
            mf =self(nominal)
            delta_ts = 2.*(self(maxl) - mf )
        if lower>0:
            pull = np.sign(maxl-mf) * np.sqrt(max(0, delta_ts))\
            if nominal is not None else None
            summary  = dict(
                maxl=maxl,
                lower=lower, upper=upper,
                ts=self.ts, # poiss.ts,
                err=err,
                )
        else:
            # just an upper limit
            pull = -np.sqrt(max(0, delta_ts)) if nominal is not None else None
            summary= dict(maxl=0,lower=0, upper=poiss.cdfinv(0.05), ts=0,
                 err=err,
                )
        if nominal is not None:
            summary.update(delta_ts=delta_ts, pull=pull)
        return summary