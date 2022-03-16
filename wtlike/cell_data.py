# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/08_cell_data.ipynb (unless otherwise specified).

__all__ = ['time_bin_edges', 'contiguous_bins', 'CellData', 'concatenate_cells', 'partition_cells']

# Cell
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rc('font', size=14)
plt.rc('figure', facecolor='white')
from .config import *
from .sources import *
from .source_data import *
from .loglike import LogLike, PoissonRep, poisson_tolerance

# Cell
def time_bin_edges(config, exposure, tbin=None):
    """Return an interleaved array of time bin, or cell start/stop values

    - exposure -- the weighted exposure table derived from the spacecraft info and the source. Used only
        to establish nominal start/stop
    - tbin: an array (a,b,d), default config.time_bins to specify binning

        interpretation of a, b, d:

        a:  if > 50000, interpret as MJD for start
            if < 0, back from stop
            otherwise, offset from exposure start

        b:  if > 50000, interpret MJD value for stop
            if > 0, increment from start
            otherwise, offset from exposure stop

        d : if positive, the day bin size
            if 0; return contiguous bins


    """
    # nominal total range, MJD edges
    start = np.round(exposure.start.values[0])
    stop =  np.round(exposure.stop.values[-1])

    a, b, step = tbin if tbin is not None else config.time_bins

    if a>50000: start=a
    elif a<0: start = stop+a
    else : start += a

    if b>50000: stop=b
    elif b>0: stop = start+b
    else: stop += b

    if step<=0:
        return contiguous_bins(exposure.query(f'{start}<start<{stop}'),)

    # adjust stop
    nbins = int((stop-start)/step)
    assert nbins>0, 'Bad binning: no bins'
    stop = start+(nbins)*step
    u =  np.linspace(start,stop, nbins+1 )

    # make an interleaved start/stop array
    v = np.empty(2*nbins, float)
    v[0::2] = u[:-1]
    v[1::2] = u[1:]
    return v

# Cell
def contiguous_bins(exposure, min_gap=20, min_duration=600):

    """ return a start/stop interleaved array for contiguous intervals

    """

    stop = exposure.stop.values
    start = exposure.start.values

    # interleave  the starts and stops
    ssint = np.empty(2*len(start))
    ssint[0::2] = start
    ssint[1::2] = stop

    # Tag the (stpp,start) pairs < 10 sec as  not adjacent
    not_adjacent = np.diff(ssint)[1::2] > min_gap/(24*3600) ;
    #print(f'{sum(not_adjacent)} (start,stop) pairs are not closer than {min_gap} s')

    # make a mask, keep ends
    mask = np.empty(2*len(start), bool)
    mask[0] = mask[-1] = True
    #

    # insert into mask -- keep only the (stop,start) pairs  which are not adjacent
    mask[1:-2:2] = not_adjacent
    mask[2:-1:2] = not_adjacent

    # apply mask, split into start and stop
    keep = ssint[mask]
    return keep

# Cell
class CellData(SourceData):
    """Manage a set of cells generated from a data set

        Invoke superclass to load photon data and exposure for the source.

        * time_bins, default config.time_bins

        The `e` cell entry is the weighted exposure for the cell in units $cm^2\ Ms$.
        """

    def __init__(self, *pars, **kwargs):
        """
        """
        bins = kwargs.pop('bins', kwargs.pop('time_bins', Config().time_bins))

        #  load source data
        super().__init__(*pars, **kwargs )

        self.rebin(bins)
        self.parent = None

    def rebin(self, newbins):
        """bin, or rebin
        """
        photon_data = self.photons
        self.time_bins = newbins
        self.cell_edges = edges = time_bin_edges(self.config, self.exposure, newbins)
        if self.config.verbose>0:
            step = newbins[2]
            self.step_name = 'orbit-based' if step<=0 else bin_size_name(step)
            print(f'CellData.rebin: Bin photon data into {int(len(edges)/2)} {self.step_name}'\
                  f' bins from {edges[0]:.1f} to {edges[-1]:.1f}')
        # note need to take care of interleave
        #self.binexp = self.binned_exposure( edges ) [0::2]

        self.make_cells()

    def get_exposure_per_cell(self, exposure_factor=1e-6):
        """
        Return a dict of arrays per cell:
        - exp -- exposure, in cm^2 Ms units, if exposure_factor==1e-6
        - costh -- mean cos theta per cell
        - exp_energy if exp_fract in the exposure DF, set exposure energy

        Note: total exposure vs. energy is:
            t =df.apply(lambda r: np.array(r.exp * np.array(r.exp_fract), float), axis=1).values
            u = np.vstack(t)
        """
        exp = self.exposure.exp.values
        costh = self.exposure.cos_theta.values
        # the cell index list
        eci = np.searchsorted(self.exposure.stop, self.cell_edges).reshape(len(self.cell_edges)//2,2)
        cell_exp = np.array([exp[slice(*ecx)].sum()*exposure_factor for ecx in eci], np.float32) #np.float32)
        cell_costh =np.array([costh[slice(*ecx)].mean() for ecx in eci], np.float32) #np.float32)

        ef = self.exposure.get('exp_fract', None)
        if ef is not None:
            efa = np.array([np.array(x,float) for x in ef])
            efe = efa.T * exp
            cee = np.array([efe.T[slice(*ecx)].sum(axis=0)*exposure_factor for ecx in eci], np.float32)
        else:
            cee = None
        return dict(exp=cell_exp, costh=cell_costh, exp_energy= cee)

    def get_weights_per_cell(self):
        """
        Return a list of arrays of the weights per cell
        """
        wts = self.photons.weight.values.astype(np.float32)
        # use photon times to get cell index range into photon list
        photon_cell = np.searchsorted(self.photons.time, self.cell_edges).reshape(len(self.cell_edges)//2,2)
        return  [wts[slice(*cell)] for cell in photon_cell]

    def make_cells(self, exposure_factor=1e-6):
        """
        Generate a "cells" DataFrame, binning according to self.cell_edges

        - exposure_factor --  recast exposure as cm^2 * Ms if `exposure_factor`==1e-6`

        Thus the `e` cell entry is the weighted exposure for the cell in units $\mathrm{cm^2 Ms}$.
        """
        # ncells = len(self.cell_edges)//2
        start,stop = self.cell_edges[0::2], self.cell_edges[1::2]
        center = (start+stop)/2
        width = (stop-start)
        etot = self.exptot*exposure_factor
        Sk, Bk = self.S/etot, self.B/etot

        # invoke the helper functions to make lists to incorporate
        weights = self.get_weights_per_cell()
        ec = self.get_exposure_per_cell(exposure_factor)
        cell_exp, cell_costh = ec['exp'], ec['costh']

        self.cells = pd.DataFrame.from_dict(dict(
            t=center,
            tw=width,
            e = cell_exp,
            ctm = cell_costh,
            n = [len(w) for w in weights],
            w = weights,
            S = cell_exp*Sk,
            B = cell_exp*Bk,
        ))

    def update(self): pass # virtual

    def view(self, *pars, exp_min=None, no_update=False):
        """
        Return a "view", a copy of this instance with a perhaps a different set of cells

        - pars -- start, stop, step  to define new binning. Or start, step, or just step
           start and stop are either MJD values, or offsets from the start or stop.
           step -- the cell size in days, or if zero, orbit-based binning

        - exp_min [None] -- If specified, a different minimum exposure, in cm^2 Ms units to use for fitting
            from.

        - no_update -- avoid fitting the cells if invoked by LightCurve or WtLike
        """
        import copy
        if self.config.verbose>1:
            print(f'Making a view of the class {self.__class__}')
        r = copy.copy(self)
        if exp_min is not None: r.exp_min=exp_min

        if len(pars)==3:
            newbins = pars
        elif len(pars)==2: # new limits, same interval
            newbins = (pars[0], pars[1], self.time_bins[0])
        elif len(pars)==1:
            if type(pars[0])==tuple:
                newbins = pars[0]
            else:
                newbins = (self.time_bins[0], self.time_bins[1], pars[0])
        else:
            newbins=None

        if newbins is not None:
            r.rebin(newbins)
            if not no_update:
                r.update()

        r.parent = self
        return r

    def __repr__(self):
        return f'''{self.__class__}:
        {len(self.exposure):,} intervals from {self.cell_edges[0]:.1f} to {self.cell_edges[-1]:.1f} for source {self.source_name}
        S {self.S:.2f}  B {self.B:.2f} '''


    def concatenate( self ):
        """
        Combine this set of cells to one
        Return a dict with summed n, S, B, and concatenated w
        """

        cells = self.cells

        newcell = dict()

        if 't' in cells:
            ca, cb =cells.iloc[0], cells.iloc[-1]
            newcell.update(dict(t= 0.5*(ca.t-ca.tw/2 + cb.t+cb.tw/2), tw=cb.t-ca.t ))

        for col in ' n e S B'.split():
            newcell[col] = cells[col].sum()
        newcell['w'] = np.concatenate(list(cells.w.values))
        return newcell


    def full_likelihood(self ):
        """Concatentate all the cells, return a LogLike object
        """
        return LogLike(self.concatenate())

    def plot_concatenated(self, fignum=1, ax=None, **kwargs):
        """Likelihood function, with fit for concatenated data
        """
        lka = self.full_likelihood()
        fig,ax = plt.subplots(figsize=(4,2), num=fignum) if ax is None else (ax.figure, ax)
        lka.plot(ax=ax, **kwargs)
        return fig

    def x_view(self, corr_func):
        """ Return a new view with a modified exposure and same binning
        - corr_func -- exposure correction factor function of (MJD) time
        """
        import copy
        if self.config.verbose>0:
            print(f'CellData.x_view: Making a view of the class {self.__class__.__name__} with adjusted exposure')
        r = copy.copy(self)
        # need a copy of this object to change the exp field?
        r.exposure = copy.copy(self.exposure)
        t = r.exposure.stop
        mod_exp = r.exposure.exp * corr_func(t)
        r.exposure.drop('exp', axis=1)
        r.exposure.loc[:,'exp'] = mod_exp

        r.rebin(self.time_bins)
        r.update()
        r.parent = self

        return r

    def phase_view(self, period, nbins=25, reference='2008'):
        """ Return a "phase" view, in which the cell time binning is according to phase.

        * period -- 'year' | 'precession' | float
        * reference -- a UTC date for aligning the bins.
        """
        ref = 0 if not reference else MJD(reference)

        period = dict(precession=53.05, year=365.25).get(period, period)
        assert np.isreal(period), f'The specified period, "{name}", is not a real number'

        if self.config.verbose>0:
            print(f'CellData.phase_view: Create phase view, {nbins} bins with period {period} days.')

        # helper function that returns the bin number as a float in [0,period)
        binner = lambda t: np.mod(t-ref,period)/period * nbins

        # adjust start to correspond to edge of bin

        # create a view with nbins per period and get the cells
        st = self.start # the start of data taking
        self.reference_bin =strefbin = binner(st)
        stnew = st +np.mod(-strefbin,1)*period/nbins
        view = self.view(stnew, 0, period/nbins, no_update=True)
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
        view.update() # does fit
        view.is_phase = True # tag to choose proper plot
        view.period = period
        return  view


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
    # get indices of  cell indexes just beyond each edge time
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
        check = cells.iloc[a:b]
        subset = check[~pd.isna(check.n)]

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