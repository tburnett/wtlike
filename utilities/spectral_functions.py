"""
These are a subset of the spectral functions defined in uw/like/Model.py, those actually used in pointlike or xFGL 
catalogs. 
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

class SedFun:
    """Create a log-log SED version of dN/dE:
    Function of log10(e/1GeV) returns log(e**2 * dN/dE(e) * 1e6) 
    
    """
    def __init__(self, cntfun):
        self.f = cntfun
    def __call__(self, loge):        
        e = 1e3 * 10**loge         
        return np.log(e**2 * self.f(e) * 1e6)
    
    @property
    def peak(self):
        """ log10 of peak energy/1 GeV
        """
        from scipy.optimize import minimize
        return  minimize(lambda x: -self(x), [1.0], bounds=[(-1,3)] ).x[0]
    
    def plot(self, show_peak=True, ax=None, **kwargs):
        fig, ax = plt.subplots(figsize=(4,2)) if ax is None else (ax.figure, ax)
        x = np.linspace(-1, 2,)
        kw = dict(); kw.update(kwargs)
        ax.plot(x, self(x))
        ax.set(**kw)
        if show_peak:
            pk = self.peak
            ax.plot(pk, self(pk), 'o')
        return fig
 

class FluxModel():
    
    emin, emax = 1e2, 1e5

    def __init__(self, pars, e0=1000, errors=[]):
        self.pars=pars
        self.e0=e0
        self.errors=errors

    def __repr__(self):

        p = self.pars
        t = f'{self.__class__.__name__}({p[0]:.2e}, '
        return t+', '.join([(f'{x:.2f}' if x<10 else f'{x:.0f}') for x in p[1:]])+')'

    @property
    def photon_flux(self):
        """photon flux (cm-2 s-1)"""
        #return quad(self, self.emin, self.emax)[0]
        # or integrate e*f(e) over log(e) for much better precision
        return quad( lambda loge: np.exp(loge)* self(np.exp(loge)),
                     np.log(self.emin), np.log(self.emax) )[0] 

    @property
    def energy_flux(self):
        """ Energy flux (erg cm-2 s-1)
        """
        func = lambda e: self(e) * e
        return 1.60218e-6 * quad(func, self.emin, self.emax)[0]
    
    @property
    def sedfun(self):
        return SedFun(self)
    
    def sed_plot(self, ax=None, e0=None,
             figsize=(5,4), label='', plot_kw={}, **kwargs):
        """Make an SED for the source
        - plot_kw -- for the plot command (lw,ls,color, etc.) 
        - kwargs -- for the Axes object (xlim, ylim, etc.)
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=figsize) if ax is None else (ax.figure, ax)
        x =np.logspace(2,5,61)

        trans = lambda x: (x/1e3, self(x)*x**2 * 1e6)
        pkw = dict(ls='-', lw=2)
        pkw.update(plot_kw)
        lines=ax.loglog(*trans(x), label=label, **pkw)
        if e0 is not None:
            ax.plot(*trans(e0),'o', color=lines[0].get_color())
        ax.grid(alpha=0.5)
        kw = dict(xlabel='Energy (GeV)',
                  ylabel=r'$\mathrm{Energy\ Flux\ (eV\ cm^{-2}\ s^{-1})}$',
                  xlim=(x.min()/1e3,x.max()/1e3),
                  ylim=(0.1, 100),
                 )
        kw.update(kwargs)
        ax.set(**kw)
        if label!='': ax.legend(fontsize=10)

    def i_flux(self, emin=100,emax=1e5,
                e_weight=0,
                cgs=False, #error=False,two_sided=False,
                quiet=False):
        """ Return the integral flux, \int_{emin}^{emax} dE E^{e_weight} dN/dE.
            e_weight = 0 gives the photon flux (ph cm^-2 s^-1)
            e_weight = 1 gives the energy flux (MeV cm^-2 s^-1) (see kwargs)

            Optional keyword arguments:

                =========    =======================================================
                Keyword      Description
                =========    =======================================================
                emin         [100] lower bound in MeV
                emax         [1e5] upper bound in MeV
                e_weight     [0] energy power by which to scale dN/dE
                cgs          [False] if True, energy in ergs, else in MeV
                error        [False] if True, return value is a tuple with flux and estimated error
                two_sided    [False] if True, return value is a triple with flux, estimated high and low errors
                quiet        [False] set to suppress error message
                =========    =======================================================
        """

        # check for a divergent flux
        # ADW: This is very bad for hard sources...
        #if 100*self(100) <= 1e5*self(1e5): emax = min(5e5,emax)

        try:
            func    = self if e_weight == 0 else lambda e: self(e)*e**e_weight
            units  = 1.60218e-6**(e_weight) if cgs else 1. #extra factor from integral!
            epsabs = min(func(emin),func(emax))*1e-10 # needed since epsrel does not seem to work
            flux    =  units*quad(func,emin,emax,epsabs=epsabs,full_output=True)[0]
            # if error:
            #     # will silently ignore 'free' parameters without errors
            #     mask = (self.free) * (self.internal_cov_matrix.diagonal()>0)
            #     args = (emin,emax,e_weight,cgs,False)
            #     d    = self.__flux_derivs__(*args)[mask]
            #     dt  = d.reshape( (d.shape[0],1) ) #transpose
            #     try:
            #         err = (d * self.internal_cov_matrix[mask].transpose()[mask] * dt).sum()**0.5
            #     except:
            #         err=np.nan
            #     if not two_sided:
            #         return (flux,err)
            #     else: #use log transform to estimate two-sided errors
            #         log_err  = err/flux
            #         log_flux = np.log(flux)
            #         return (flux,np.exp(log_flux+log_err)-flux,flux-np.exp(log_flux-log_err))

            return flux
        except Exception as msg:
            if not quiet:
                print( f'Encountered a numerical error, "{msg}", when attempting to calculate integral flux.',
                out=sys.stderr)
            return np.nan if not error else ([flux, np.nan,np.nan] if two_sided else [flux, np.nan])
        
    def curvature(self, e=None):
        """return estimate of the curvature, using numerical derivatives
        This is exactly the parameter beta for a LogParabola 
        e: float
            energy in Mev, default e0
        """
        from scipy.misc import derivative
        def dfun(x):
            def fun(x):
                """flux as function of x=ln(e0/e)"""
                e = self.e0 * np.exp(-x)
                return self(e)
            return derivative(fun,x, dx=0.01)/fun(x)
        x=0 if e is None else np.log(self.e0/e)

        return -0.5*derivative(lambda x: dfun(x), x, dx=0.01)


class PowerLaw(FluxModel):

    def __call__(self, e):
        n0,gamma = self.pars
        return n0*(self.e0/e)**gamma

class LogParabola(FluxModel):

    def __call__(self, e):
        n0,alpha,beta,e_break=self.pars
        x = np.log(e_break/e)
        y = (alpha - beta*x)*x
        return n0*np.exp(y)

class PLSuperExpCutoff(FluxModel):

    def __call__(self,e):
        #print('WARNING: check e0!')
        n0,gamma,cutoff,b=self.pars
        return n0*(self.e0/e)**gamma*np.exp(-(e/cutoff)**b)

class PLSuperExpCutoff4(FluxModel):
    """ adaptation of Matthew's code in uw/like/Models.py
        to just implement the function

        >>> m = PLSuperExpCutoff4([1e-11,2,0.8,1], e0=1e3)
        >>> expected = np.asarray([5.54576203e-09, 3.25605721e-10, 
        >>>     1.00000000e-11, 4.71063799e-16])
        >>> np.allclose(m(10**np.arange(1,5)),expected)
        >>> True
    """
    def __call__(self,e):
        n0,gamma,d,b = self.pars
        x = e*(1./self.e0)
        return n0*x**(-gamma+d/b)*np.exp((d/b**2)*(1-x**b))

# external access: by name, ...
specfun_dict = dict(
    PowerLaw=PowerLaw,
    LogParabola=LogParabola,
    PLSuperExpCutoff=PLSuperExpCutoff,
    PLSuperExpCutoff4=PLSuperExpCutoff4,
    )

# or just run the appropriate construcctor
def spectral_function(name, *pars, **kwargs):
    return specfun_dict[name](*pars, **kwargs)


class MultiSED:
    """ Generate a table of SED plots
    - tooltips
    - caption
    """    
    def __init__(self, n, start=0, ncols=10, caption='',
                 tooltips=None, **kwargs):
        """
        * n -- Number of plots
        * start -- for annotating with sequential index
        * ncols -- number of columns
        * caption -- use if invoked with `show`
        * tooltips -- relevant for show
        
        """

        nrows = (n+ncols-1)//ncols
        self.caption = caption
        self.tooltips = tooltips
        self.fig, axx = plt.subplots(ncols=ncols, nrows=nrows, 
                figsize=(12,1.2*nrows), sharex=True, sharey=True,
                gridspec_kw=dict(
                    left= 0.01,  right=0.99, hspace=0.05,
                    top = 0.95, bottom=0.01, wspace=0.05))
        self.axx = axx.flatten()

        # keywords for all the Axes objects
        self.ax_kw = dict(xticks=[], yticks=[], xlabel='', ylabel='', ylim=(0.05,5),)
        self.ax_kw.update(kwargs)
        
        # this could be flexible, perhaps use names
        self.annotate(start)
        self.axx[0].set(**self.ax_kw)
        
    def annotate(self, start=0):
        
        # add id's to top right of all subplots -- here just an index
        for j,ax in enumerate(self.axx):
            ax.text(0.95,0.95, str(start+j), transform=ax.transAxes, 
                fontsize=10, ha='right',va='top')
            
    def plots(self, specfuns, **kwargs ):
        """
        * specfuns -- a list of specral functions to apply to set of Axes
        * kwargs -- plotting options, like color, lw, ls        
        """
        for ax, f in zip(self.axx, specfuns):
            f.sed_plot(ax, plot_kw=dict(kwargs))  
            ax.set(**self.ax_kw)
        return self
    
    def _repr_html_(self):
        # for ipython display
        from utilities.ipynb_docgen import FigureWrapper
        return FigureWrapper(
            self, 
            tooltips=self.tooltips,
            caption=self.caption,
        )._repr_html_()