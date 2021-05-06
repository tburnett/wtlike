# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/07_cell_data.ipynb (unless otherwise specified).

__all__ = ['CellData', 'concatenate_cells', 'partition_cells']

# Cell
import os
import numpy as np
import pandas as pd
from .config import *
from .source_data import *
from .loglike import LogLike, PoissonRep

# Cell
class CellData(SourceData):
    """Manage a set of cells generated from a data set

        Invoke superclass to load photon data and exposure for the source.

        * time_bins, default config.time_bins

        """

    def __init__(self, *pars, **kwargs):
        """

        """
        bins = kwargs.pop('bins', kwargs.pop('time_bins', Config().time_bins))
        #  load phots
        super().__init__(*pars, **kwargs )
        self.use_uint8  = self.config.use_uint8
        self.rebin(bins)

    def rebin(self, newbins):
        """bin, or rebin
        """
        photon_data = self.p_df
        self.cell_edges = edges = time_bin_edges(self.config, self.e_df, newbins)
        if self.config.verbose>0:
            step = newbins[2]
            self.step_name = 'orbit-based' if step<=0 else bin_size_name(step)
            print(f'Bin photon data into {int(len(edges)/2)} {self.step_name}'\
                  f' bins from {edges[0]:.1f} to {edges[-1]:.1f}')

        # note need to take care of interleave
        expose = self.binned_exposure( edges ) [0::2]
        self.fexposure=(expose/self.exptot).astype(np.float32)
        self.get_cells()

    def get_cells(self):
        """
        Generate the cell DataFrame
        """


        # restrict photons to range of bin times
        photons = self.photons.query(f'{self.cell_edges[0]}<time<{self.cell_edges[-1]}')

        # use photon times to get indices into photon list
        edges = np.searchsorted(photons.time, self.cell_edges)

        wts = photons.weight.values
        start,stop = self.cell_edges[0::2], self.cell_edges[1::2]
        center = (start+stop)/2
        width = (stop-start)
        cells = []
        ek = np.append(edges[0::2], edges[-1])

        for k, (t, tw, e) in enumerate( zip(center, width, self.fexposure) ):
            w = wts[ek[k]:ek[k+1]]
            n = len(w)
            cells.append(dict(t=t, tw=tw,
                              e=e,
                              n=n,
                              w=w,
                              S=e*self.S,
                              B=e*self.B,
                             )
                        )
        self.df = self.cells =  pd.DataFrame(cells)
        return self.df

    def update(self): pass # virtual

    def view(self, newbins):
        """Return a "view": a new instance of this class with a different set of cells

        """
        import copy
        if self.config.verbose>1:
            print(f'Making a view of the class {self.__class__}')
        r = copy.copy(self)
        # check to see if new binning is contained
        r.rebin(newbins)

        r.update()
        return r


    def __repr__(self):
        return f'''{self.__class__}:
        {len(self.fexposure)} intervals from {self.cell_edges[0]:.1f} to {self.cell_edges[-1]:.1f} for source {self.source_name}
        S {self.S:.2f}  B {self.B:.2f} '''


    def concatenate( self ):
        """
        Combine this set of cells to one
        Return a dict with summed n, S, B, and concatenated w
        """

        cells = self.df

        newcell = dict()

        if 't' in cells:
            ca, cb =cells.iloc[0], cells.iloc[-1]
            newcell.update(dict(t= 0.5*(ca.t-ca.tw/2 + cb.t+cb.tw/2), tw=cb.t-ca.t ))

        for col in ' n S B'.split():
            newcell[col] = cells[col].sum()
        newcell['w'] = np.concatenate(list(cells.w.values))
        return newcell


    def full_likelihood(self ):
        """Concatentate all the cells, return a LogLike object
        """
        return LogLike(self.concatenate())

    def plot_concatenated(self, fignum=1, **kwargs):
        """Likelihood function, with fit for concatenated data
        """
        import matplotlib.pyplot as plt
        lka = self.full_likelihood()
        fig,ax = plt.subplots(figsize=(4,2), num=fignum)
        lka.plot(ax=ax, **kwargs)
        return fig

# Cell
def concatenate_cells( cells):
    """
    Combine a group of cells to one
    - cells: dataframe with cells containing  n, w, S, B<br>
            Optionally, if $t$ is present, generate t and tw
    Return a dict with summed n, S, B, and concatenated w
    """
    newcell = dict()
    if 't' in cells:
        ca, cb =cells.iloc[0], cells.iloc[-1]
        newcell.update(dict(t= 0.5*(ca.t-ca.tw/2 + cb.t+cb.tw/2), tw=cb.t-ca.t ))

    for col in ' n S B'.split():
        newcell[col] = cells[col].sum()
    newcell['w'] = np.concatenate(list(cells.w.values))
    return newcell

# Cell
def partition_cells(config, cells, edges):
    """ Partition a set of cells
     - cells -- A DataFrame of cells
     - edges  -- a list of edge times delimiting boundaries between cells

    Returns a DataFrame of combined cells, with times and widths adjusted to account for missing cells

    """
    # get indices of  cell idexes just beyond each edge time
    ii = np.searchsorted(cells.t, edges)

    # Get the appropriate boundary times to apply to combined cells
    # this is complicated by missing cells, need to put boundary in gaps if ncessary
    ileft = ii[:-1]
    cleft = cells.iloc[ileft ]
    tleft =  (cleft.t - cleft.tw/2).values
    iright = ii[1:]-1
    cright = cells.iloc[iright ]
    tright = (cright.t+cright.tw/2).values
    betweens = 0.5*(tleft[1:] + tright[:-1])
    tboundary = np.append(np.insert(betweens, 0, tleft[0]), tright[-1])

    # now combine the cells,
    newcells = []
    for k in range(len(ii)-1):
        a,b = ii[k:k+2]
        subset = cells.iloc[a:b];

#         ca, cb = subset.iloc[0], subset.iloc[-1]
#         newcell = dict(t= 0.5*(ca.t-ca.tw/2 + cb.t+cb.tw/2)  )
        tl, tr = tboundary[k:k+2]
        newcell = dict(t=0.5*(tl+tr), tw=tr-tl)

        for col in 'e n S B'.split():
            newcell[col] = subset[col].sum()
        newcell['e'] /= len(subset)
        newcell['w'] = np.concatenate(list(subset.w.values)) #np.array(w, np.uint8)
        newcells.append(newcell)
    return pd.DataFrame(newcells)