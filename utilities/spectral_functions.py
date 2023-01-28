"""
These are a subset of the spectral functions defined in uw/like/Model.py, those actually used in pointlike or xFGL 
catalogs. 
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


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
    
    def sed_plot(self, ax=None, e0=None,
             figsize=(5,4), label='', plot_kw={}, **kwargs):
        """Make an SED for the source

        - kwargs -- for the Axes object (xlim, ylim, etc.)
        """
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=figsize) if ax is None else (ax.figure, ax)
        x =np.logspace(2,5,61)

        trans = lambda x: (x/1e3, self(x)*x**2 * 1e6)
        lines=ax.loglog(*trans(x), '-', lw=2, label=label, **plot_kw)
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
        if label!='': ax.legend()

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