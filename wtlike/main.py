# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/90_main.ipynb (unless otherwise specified).

__all__ = ['WtLike']

# Cell
import numpy as np
from .bayesian import get_bb_partition
from .lightcurve import fit_cells, LightCurve, flux_plot
from .cell_data import partition_cells
from .config import MJD


class WtLike(LightCurve):
    """
    Summary
    ---------
    There are three layers of initialization, implemented in superclasses,
    each with parameters. The classnames, associated parameters and data members set:

    SourceData -- load photons and exposure
        parameters:
          - source : name, a PointSource object, or a Simulation object
          - config [Config()] : basic configuration
          - week_range [None] : range of weeks to load
          - key [''] : the cache key: '' means construct one with the source name, None to disable
          - clear [False] : if using cache, clear the contents first
        sets:
          - photons
          - exposure

    CellData -- create cells
        parameters:
          - time_bins [Config().time_bins] : binning: start, stop, binsize
        sets:
          - cells

    LightCurve -- likelihood analysis of the cells
        parameters:
          - e_min [10] -- threshold for exposure (cm^2 units)
          - n_min [2]  -- likelihood has trouble with this few
          - lc_key [None] -- possible cache for light curve
        sets:
          - fits, fluxes

    WtLike (this class) -- no parameters (may add BB-specific ones)
        Implements:  bb_view, plot_BB
        sets:
          - bb_flux  (only if bb_view invoked)

    """
    def bb_view(self, p0=0.05, key=None, clear=False):
        """Return a view with the BB analysis applied

        - p0 -- false positive probability parameter

        Its `plot` function will by default show an overplot on the parent's data points.
        """
        #  a new instance
        r = self.view()

        # bb analysis on this to make new  set of cells and poisson fits
        bb_edges  = get_bb_partition(self.config, self.fits,  p0=p0, key=key, clear=clear)
        r.cells = partition_cells(self.config, self.cells, bb_edges)

        r.fits = fit_cells(self.config, r.cells, )
        r.isBB = True
        r.bayes_p0 = p0
        return r

    def plot(self, *pars, **kwargs):
        # which view type is this?
        if getattr(self, 'isBB', False):
            return self.plot_bb(*pars, **kwargs)
        elif getattr(self, 'is_phase', False):
            return self.plot_phase(*pars, **kwargs)
        else:
            return super().plot(*pars, **kwargs)

    def plot_bb(self, ax=None, **kwargs):
        """Plot the light curve with BB overplot
        """
        import matplotlib.pyplot as plt
        self.check_plot_kwargs(kwargs)
        figsize = kwargs.pop('figsize', (12,4))
        fignum = kwargs.pop('fignum', 1)
        ts_min = kwargs.pop('ts_min',-1)
        source_name =kwargs.pop('source_name', self.source_name)
        fig, ax = plt.subplots(figsize=figsize, num=fignum) if ax is None else (ax.figure, ax)


        colors = kwargs.pop('colors', ('lightblue', 'wheat', 'blue') )
        flux_plot(self.parent.fits, ax=ax, colors=colors, source_name=source_name,
                  label=self.step_name+' bins', **kwargs)
        flux_plot(self.fits, ax=ax, step=True,
                  label=f'BB (p0={100*self.bayes_p0:.0f}%)', zorder=10,**kwargs)
        ax.grid(alpha=0.5)
        fig.set_facecolor('white')
        return fig

    def phase_view(self, period, nbins=25, reference='2008'):
        """ Return a "phase" view, in which the cell time binning is according to phase.

        * reference -- a UTC data for aligning the bins.
        """
        ref = 0 if not reference else MJD(reference)

        # helper function that returns the bin number as a float in [0,period)
        binner = lambda t: np.mod(t-ref,period)/period * nbins

        # adjust start to correspond to edge of bin

        # create a view with nbins per period and get the cells
        st = self.start # the start of data taking
        self.reference_bin =strefbin = binner(st)
        stnew = st +np.mod(-strefbin,1)*period/nbins
        view = self.view(stnew, 0, period/nbins)
        cells = view.cells
        bw = 1/nbins

        def concat(pcells, t):
            newcell = dict(t=t, tw=bw)
            for col in 'n e S B'.split():
                newcell[col] = pcells[col].sum()
            newcell['w'] = np.concatenate(list(pcells.w.values))
            return newcell

        # concatenate all in the same phase bin--note cyclic rotation
        k = int(strefbin)
        fcells = [concat(cells.iloc[((ibin-k)%nbins):-1:nbins], (ibin+0.5)*bw)  for ibin in range(nbins) ]

        view.cells = pd.DataFrame(fcells)
        view.update()
        view.is_phase = True
        view.period = period
        return  view

    def plot_phase(self, ax=None, **kwargs):
        """Plot a phase lightcurve

        """
        kw = dict(ylim=(0.975, 1.025), xlim=(0,1) )
        kw.update(kwargs)
        fig, ax = plt.subplots(figsize=(10,5)) if ax is None else (ax.figure, ax)
        fig = super().plot(ax=ax, xlabel=f'phase for {self.period}-day period');
        ax.set(**kwargs );
        ax.axhline(1.0, color='grey');
        return fig