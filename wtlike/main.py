# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/90-main.ipynb (unless otherwise specified).

__all__ = ['WtLike']

# Cell
from .bayesian import get_bb_partition
from .lightcurve import fit_cells, LightCurve, flux_plot
from .cell_data import partition_cells

class WtLike(LightCurve):
    """
    Inherit from `LightCurve` and add Bayesian Block capability
    """
    def bb_view(self, key=None, clear=False, **kwargs):
        """Return a view with the BB analysis applied

        Its `plot` function will by default show an overplot on the parent's data points.
        """
        #  a new instance
        r = self.view()

        # bb analysis on this to make new  set of cells and poisson fits
        bb_edges  = get_bb_partition(self.config, self.fits,  key=key, clear=clear)
        r.cells = partition_cells(self.config, self.cells, bb_edges)

        r.fits = fit_cells(self.config, r.cells, )
        r.isBB = True
        return r

    def plot(self, *pars, **kwargs):

        if getattr(self, 'isBB',  None) is None:
            return super().plot(*pars, **kwargs)
        else:
            return self.plot_BB(*pars, **kwargs)

    def plot_BB(self, ax=None, **kwargs):
        """Plot the light curve with BB overplot
        """
        import matplotlib.pyplot as plt

        figsize = kwargs.pop('figsize', (12,4))
        fignum = kwargs.pop('fignum', 1)
        ts_min = kwargs.pop('ts_min',-1)
        source_name =kwargs.pop('source_name', self.source_name)
        fig, ax = ig, ax = plt.subplots(figsize=figsize, num=fignum) if ax is None else (ax.figure, ax)


        colors = kwargs.pop('colors', ('lightblue', 'wheat', 'blue') )
        flux_plot(self.config, self.parent.fits, ax=ax, colors=colors, source_name=source_name,
                  label=self.step_name+' bins', **kwargs)
        flux_plot(self.config, self.fits, ax=ax, step=True,
                  label='BB overlay', zorder=10,**kwargs)

        fig.set_facecolor('white')
        return fig
