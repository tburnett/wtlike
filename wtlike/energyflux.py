# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/11_energyflux.ipynb.

# %% auto 0
__all__ = ['pull_plot', 'EnergyCells', 'FluxFixer']

# %% ../nbs/11_energyflux.ipynb 11
import sys, os
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from .source_data import SourceData
from .config import *
from .loglike import LogLike
from .poisson import Poisson
from .effective_area import EffectiveArea
from .lightcurve import flux_plot
plt.rc('font', size=12)

# %% ../nbs/11_energyflux.ipynb 16
def contiguous_idx(exposure, min_gap=20, min_duration=600):

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
    # keep = ssint[mask]
    r = (np.arange(len(mask))//2)[mask]
    
    # reshape as (start,stop) pairs X entries in exposure
    return r.reshape(len(r)//2,2)

# %% ../nbs/11_energyflux.ipynb 17
def pull_plot(pulls, ax=None, nbins=25, **kwargs):
    """ Histogram of a set of "pulls", presumed to be normal (0,1) distributed, 
    the function for which is overplotted.
    """
    from scipy import stats
    x = pulls[~np.isnan(pulls)]
    
    fig, ax = plt.subplots(1,1, figsize=(5,3)) if ax is None else (ax.figure, ax)
    if len(x)<2:
        ax.text(0.5, 0.5, '(no data)', color='red',
                transform=ax.transAxes, ha='center', va='center', )
        ax.set(xticks=[], yticks=[])
        return fig
    
    kw = dict(ylim=(0.,0.45))
    kw.update(kwargs)
    ax.set(**kwargs)
    ax.hist(pulls, np.linspace(-5,5,nbins+1), histtype='step', lw=2, density=True, )
    dom = np.linspace(-5,5,101)
    ax.plot(dom, stats.norm.pdf(dom), '--', lw=2,)# label='Normal(0,1)')
    ax.grid(alpha=0.5)
    
    ax.text(0.95, 0.97,f'{len(x)}\n{np.nanmean(x):.3f}\n'+r'$\pm$'+f'{np.nanstd(x):.3f}',
            transform=ax.transAxes, ha='right', va='top',  
            fontdict=dict(family='monospace'), fontsize=8 )
    return fig

# %% ../nbs/11_energyflux.ipynb 18
class LivetimeHistory:
    """
    Manage the livetime history
    The spacecraft history is a sequence of 30-s intervals each with a livetime<30s and a value for z=$\cos\\theta$.
    
    This class makes a table for each contiguous interval, usually an orbit of the summed livetime for bins in z.
    
    It creates a table of 
    """
    ctbins=np.linspace(0.4,1.0,13)

    def __init__(self, wtl ):
        """
        - wtl : a WtLike object with SC history,
              the exposure DataFrame has columns start, stop, livetime, cos_theta
        """
        self.source_name = wtl.source_name
        edf = wtl.exposure
        def lt_table(edf, istart, iend):
            df = edf.iloc[istart:iend+1,]
            h,_ = np.histogram(df.cos_theta, bins=self.ctbins, weights=df.livetime)
            return  h.astype(np.float32)
        
        # makes  a table of the start/stop indices for contiguous chunks
        idx = contiguous_idx(edf)
        
        self.start = edf.start.values[idx[:,0]]
        self.stop  = edf.stop.values[idx[:,1]]
        lt    = np.array([lt_table(edf, a,b) for (a,b) in idx], np.float32 )
        self.cum_lt  = np.insert(np.cumsum(lt, axis=0), 0, np.zeros(12), axis=0);
        
    def __repr__(self):
        
        return f'Spacecraft history for source {self.source_name}: a table of {self.cum_lt.shape} entries'
        
    def __call__(self, index_range=None):
        """Return the live time array for the range of orbit indices

        """    
        if index_range is None:
            return self.cum_lt[-1]

        return self.cum_lt[index_range[1]] - self.cum_lt[index_range[0]]
    
    #---------------------------------------------------
    #        plot routines
    #--------------------------------------------------
 
    def lt_vs_ct(self, ax=None):
        """ Plot of livetime vs $\cos\\theta$.
        """

        fig, ax = plt.subplots(figsize=(5,3)) if ax is None else (ax.figure, ax)
        cumx  = self.cum_lt[-1] 
        ax.stairs(100* cumx/sum(cumx), self.ctbins, lw=2)
        ax.set(xlim=(0.4,1.0), xlabel=r'$\cos\theta$', ylabel='Fraction (%)')
        ax.grid(alpha=0.5)
        
    def lt_vs_time(self, ax=None, tbinsize=120):
        """ Plot of livetime fraction vs time.
        """

        # sum over cos theta to get cumulative total livetime 
        cuml = np.sum(self.cum_lt,axis=1)

        # the time interval array of bin edges
        tbins = np.arange(self.start[0], self.stop[-1], tbinsize)
        # index of the time dimension 
        cum_idx  = np.searchsorted(self.start, tbins)

        fig, ax = plt.subplots(figsize=(12,2)) if ax is None else (ax.figure, ax)
        ax.stairs(np.diff(cuml[cum_idx])/(tbinsize*24*3600), tbins, lw=2);
        ax.set(xlabel='Time (MJD)', ylabel='Livetime fraction');
        ax.grid(alpha=0.5);

# %% ../nbs/11_energyflux.ipynb 19
class EnergyCells(object):
    r"""
    Create a set of cells using energy-dependent livetime
    
    Input is an instance of the `SourceData` class for a given source, containing photon and spacecraft info 
    
    Processing creates:
     
    * Orbit-based exposure vs $\cos\theta$, implemented by `LivetimeHistory` and using the effective area
            
    * The basic time-defined cells, implemented by `FermiInterval`
    
    Each cell has 8 sub-cells for the 4/decade energies from 100 MeV to 10 GeV. Then the likelihood function
    for each sub-cell is determined using the tools `LogLike` and `Poisson`.  
    
    """    
        
    def __init__(self, wtl:SourceData, interval,
                adjust_slope=-0.55):
                    
  
        if not isinstance(wtl, SourceData):
            raise TypeError(f'Expect {wtl} to be an instance of SourceData')
            
        self.interval = interval 
        lth = LivetimeHistory(wtl)
        self.config = wtl.config
        self.photons = wtl.photons
        self.Aeff =  EffectiveArea(file_path=self.config.datapath/'aeff_files').tabulate()
        self.adjust_slope = adjust_slope
        # orbit info from the LivetimeHistory object
        self.source = wtl.source
        self.source_name = wtl.source_name
        self.lth = lth
        self.orbit_times = 0.5*(lth.start+lth.stop)
        
        # setup all the cells/subcells with front/back mask
        self.cells = self.make_cells()
        # self.ss, self.bb = self._get_S_and_B()
        self.ftable = {}


    def make_cells(self):
        """
        create a cell DF corresponding to the interval
        """
        # get interval info
        self.fi = fi =  FermiInterval(self.interval)
        self.ncells = len(fi)
        cell_edges = fi.mm
        # get the orbit indeces for each cell 
        ci = cell_indices = np.append(np.searchsorted(cell_edges, self.orbit_times), 99999)

        # differences show transition (with guard at end)
        cidiff = np.append(np.insert(np.diff(ci),0,0),99);
        # array of where changes are
        cichange = (cidiff>0)[:-1]

        # then the corresponding starting orbit indices and cell numbers
        t = np.arange(len(ci))[cichange]

        # orbit indices for start and stop in the cell
        # effective times from those
        oi_start = np.arange(len(ci))[t[:-1]]
        oi_stop  = np.arange(len(ci))[t[1:]]
        tstart   = self.lth.start[oi_start]
        tstop    = self.lth.stop[oi_stop-1]
        # dt       = tstop-tstart

        # reconstruct cell numbers, allowing for missing ones in cell_indices array
        cell_number = np.cumsum(cidiff[:-1][cichange])-1 
        xdf = pd.DataFrame.from_dict(dict(
                    orbit_index_range = [(a,b) for a,b in zip(oi_start, oi_stop)],
  
                             )
                        )
        xdf.index=cell_number[:-1]
        xdf.loc[:,'time'] = cell_edges[xdf.index] + self.interval/2
        
        event_times = self.photons.time.values
        # 

        # make an interleaved (start,stop) array of the actual cell edges
        cell_edges = np.empty(2*len(xdf))
        cell_edges[::2] = tstart
        cell_edges[1::2] = tstop

        # get ranges of eventindices per cell
        event_cell_indices = np.searchsorted(event_times, cell_edges ).reshape(len(xdf),2)

        xdf.loc[:,'event_index_range'] = [tuple(t) for t in event_cell_indices]
        return xdf['time orbit_index_range event_index_range'.split()]

    def _get_Fbar(self):
        """Return the normalized cell exposure matrices as a (Ncells x 16 x 2) numpy array 
        """
        L = np.array([ cell['livetime'] for cell in self]).T
        F_cells = np.dot(self.Aeff, L)
        F_full = F_cells.sum(axis=2).T
        F_bar = F_cells.T / F_full

        return F_bar #.shape, F_bar[0]

    # def _get_wsum(self):
    #     """Return the summed weights as a ( 16 x 2) numpy array
    #     """ 
    #     Wsum = np.zeros(32)
    #     for k, v in  self.photons.groupby('band'):
    #         Wsum[k] = v.weight.sum()
    #     return Wsum.reshape((16,2) )

    def _get_S_and_B(self, fbmask=3):
        """Return the expected sum(w) and sum(1-w) as (Ncell X 16) arrays
        
        - fbmask-- 1,2,3 for F, B, F+B
        """ 

        vmask = [[0,0], [1,0], [0,1], [1,1]][fbmask]
        
        # sum the weights according to energy/event type index
        W = np.zeros(32)
        N = np.zeros(32)
        for k, v in  self.photons.groupby('band'):
            W[k] = v.weight.sum()
            N[k] = len(v.weight)
        # reshape to enable explicit F,B index
        W = W.reshape((16,2))
        N = N.reshape((16,2))

        # now use the normalized exposure matrix, then sum over front, back
        F_bar = self._get_Fbar() 
        # S = (F_bar * W).sum(axis=2)
        # B = (F_bar * (N-W)).sum(axis=2) 
        Sk = np.dot( F_bar * W, vmask) 
        Bk = np.dot( F_bar * (N-W), vmask) 
        
        # ad-hoc adjustment per cell depending on energy distribution

        a = self.energy_adjustment()
        return Sk * a, Bk * a

    def energy_adjustment(self):
        # either 1 or a column vector of factor for each cell
        # preferred? -0.55
        slope = getattr(self, 'adjust_slope', None)
        if slope is None: return 1
    
        ew = self.energy_moment()
        return (1 + slope*ew)[:,None]
        

    def __getitem__(self, k):
        """ #### Index: return a dict for cell # k
        It is compatible with the LogLike class
        
        * time: central time (MJD)
        * livetime array 
        
        For the events in the cell, arrays of:
        * bands: 'bands'
        * weights: 'w' 
        """
        t = dict( self.cells.iloc[k,:])
        r = dict(time=t['time'])
        oir = t['orbit_index_range']
        eir = t['event_index_range']
        r['livetime'] = self.lth(oir)
        r['bands'] = self.photons.band.values[slice(*eir)]
        r['w'] = w = self.photons.weight.values[slice(*eir)].astype(np.float32)
        r['n'] = len(w)
        if hasattr(self, 'ss'):
            # after setup, need these
            r['Sk'] = self.ss[k].astype(np.float32)
            r['Bk'] = self.bb[k].astype(np.float32)
            # sums for normalization flux
            r['S'] = np.sum(r['Sk'])
            r['B'] = np.sum(r['Bk'])
        return r

    def __len__(self):
        return self.ncells
    
    def __repr__(self):
        return f'EnergyCells: manage energy-dependent cells for source {self.source_name} with {self.ncells} {self.interval}-day cells'
    
   
    def flux_table(self, fbmask=3,  monitor=False):
        """Return a Ncells x 8 DataFrame of Poisson objects
        index is time, column name is energy index
        
        """
        ft = self.ftable.get(fbmask, None)
        if ft is not None:
            return ft
            
        def cell_w(cell, fbmask=3, nbands=8):
            """
            For a given cell, return a dict of weights for nbands, with Front/Back selection
            """
            photons = pd.DataFrame.from_dict(
                dict( eband=cell['bands']//2, 
                     evtbit=cell['bands']%2+1, 
                     weight=cell['w'])
                )
            if fbmask<3: #possible values: 1 (Front), 2(back), or 3(both)
                photons = photons[np.bitwise_and(photons.evtbit, fbmask)>0]

            ecells = {}
            for k, v in photons.groupby('eband'):    
                ecells[k] = v.weight.values
            return ecells

        def subcell_poisson(cell, Sk,Bk, fbmask=3, nbands=8):

            cw = cell_w(cell, fbmask, nbands)

            pr = [] 
            for k in range(nbands):
                w = cw.get(k,[])
                if len(w)<2:
                    # non-valid for little of no info
                    p = Poisson((0,np.nan,0))
                else:
                    subcell = dict( w=w, n=len(w), S=Sk[k], B=Bk[k],)
                    p = Poisson.from_function(LogLike(subcell)) 
                pr.append(p)
            return pr

        subs = {}
        S, B = self._get_S_and_B(fbmask)

        for i, cell in enumerate(self):
            if monitor: print(f'\rCell index: {i:4d} /{n-1:4d} ', end='')
            Sk,Bk = S[i], B[i]
            subs[cell['time']] = subcell_poisson(cell, Sk,Bk, fbmask)
        ft =  pd.DataFrame.from_dict(subs, orient='index')
        self.ftable[fbmask]=ft
        return ft
    
    
    
    def count_flux(self, fbmask=None):
        """Return a frame with the relative count flux Poisson-like log-likelihoods,  combining the energy-dependent
        likelihoods.
        """
        
        
        if fbmask is None:
            # add all Sk
            Sk,Bk, = self._get_S_and_B()
            fits = []
            for cell, S,B in zip(self, Sk, Bk):
                cell['S'] = S.sum()
                cell['B'] = B.sum()
                fits.append( Poisson.from_function(LogLike(cell) ))
            t=np.array([cell['time'] for cell in self]); width = np.diff(t)[0]; 
            
        else:
            # fbmask specified: make table add log-likelihoods 
            ft = self.flux_table(fbmask)
            valid = lambda poiss: poiss.p[1]>1
            fits = []

            for i, plist in ft.T.iteritems(): 
                filtered = list(filter(valid,plist))
                if len(filtered)==0: # completely empty
                    t = plist[0]
                elif len(filtered)==1: # only one
                    t = filtered[0]
                else:
                    t = Poisson.from_list(filtered) 
                fits.append(t) 
            t=np.array(ft.index); width = np.diff(t)[0]; 
        
        df = pd.DataFrame(dict(t=t, 
                               tw=np.full(len(t),width), 
                                fit=fits)
                         )
        return df
    
    def energy_moment(self):
        """ Return a measure of the exposure energy distribution.  
        """
        F = self._get_Fbar()[:,:8].sum(axis=2)    
        return (np.arange(8) * F).sum(axis=1) / F.sum(axis=1) - 3.5
    


    def spectral_index(self, fbmask=3, tol=0.3):
        """
        
        """
        
        def add_loglikes(plist, scale_sigs=None):
            """ 
            - plist -- list of (presumably 8) of Poisson objects to combine
            - scale_sigs -- None | list of factors to scale up sigma 
            """

            isvalid = lambda poiss: poiss.p[1]>1
            adjusted = lambda p, c : Poisson([p[0], p[1]/c**2, p[2]]) 

            if scale_sigs is not None:
                # if rescaling, replace the list 
                plist = [adjusted(poiss.p,c) for poiss,c in zip(plist, scale_sigs)]
            # remove any invalid ones
            filtered = list(filter(isvalid,plist))
            if len(filtered)==0: return plist[0]# completely empty
            if len(filtered)==1: return filtered[0] # only one

            return Poisson.from_list(filtered, tol=tol) 
        
        ft = self.flux_table(fbmask)
        ck = np.log(10)/4 * (3.5-np.arange(8))

        fits = [ add_loglikes(plist, scale_sigs=ck)  for t, plist in ft.iterrows()]
        t=np.array(ft.index)
        return pd.DataFrame(pd.Series(fits, index=t))

    #=========================================================================
    #           plot functions
    #--------------------------------------------------------------------------
    
    def plot_count_flux(self, ax=None, **kwargs):
        """Plot of the relative count flux vs. time. 

        """
        return flux_plot( self.count_flux(),
                         ax=ax, title='Count relative flux using energy', **kwargs);
   
    def plot_fractional_livetime(self, ax=None):
        """Fractional live time per cell
        """

        clt  = self.lth.cum_lt
        cell_live_time = self.cells.orbit_index_range.apply(lambda oir: sum(clt[oir[1]]-clt[oir[0]]))  

        fig, ax = plt.subplots(figsize=(12,3)) if ax is None else (ax.figure, ax)
        ax.plot( self.cells.time.values,  cell_live_time  /(self.interval*24*3600), '+g')
        ax.set(xlabel='MJD', ylabel='Livetime fraction', ylim=(0,None), )
        return fig
    
    def plot_lt_vs_w(self, ax=None, **kwargs):
        """ The total livetime per cell vs the sum of weights
        """
        fig, ax = plt.subplots(figsize=(8,3)) if ax is None else (ax.figure, ax)
        q = np.array([(sum(c['w']), sum(c['livetime'])) for c in self])
        ax.plot(q[:,0], q[:,1]/1e3, '.')
        kw = dict(xlabel='sum of weights', ylabel='live time (ks)')
        kw.update(kwargs); ax.set(**kw)
        ax.grid()
        return fig
   
    def plot_wsum_vs_ssum(self, ax=None, **kwargs):
        """Comparison of the ratio of the sum of weights per cell to the expected value, assuming a constant source.
        """
        wsum = np.array([sum(c['w']) for c in self])
        
        S, B = self._get_S_and_B()
        
        ssum = S.sum(axis=1)
        
        fig, ax = plt.subplots(figsize=(8,3)) if ax is None else (ax.figure, ax)
        ax.plot(ssum, (wsum/ssum).clip(0.8,1.2), '.');
        kw = dict(ylabel=r'\sum w\ /\ S$', xlabel=r'$S$', ylim=(0.8,1.2))
        kw.update(kwargs)
        ax.set(**kw)
        # ax.plot([0,wsum.max()], [0,wsum.max()], ':', color='orange')
        ax.axhline(1, color='grey')
        ax.grid();
        return fig

    def plot_ecell_likelihoods(self, icell, nbands=8, ax=None):
        """ Plots of the likelihoods for energy bands of a given cell
        """
        SS,BB = self._get_S_and_B()
        cell = self[icell]

        # make energy subcells
        ecells = {}
        photons = pd.DataFrame.from_dict(dict(eband=cell['bands']//2, weight=cell['w']))
        for k, v in photons.groupby('eband'):    
            ecells[k] = v.weight.values

        # sums = np.array([sum(w) for w in ecells.values()],)

        ll = [] 
        for ic in range(nbands):
            subcell = dict( w= ecells[ic], n=len(ecells[ic]), S=SS[icell][ic], B=BB[icell][ic])
            ll.append(LogLike(subcell))

        fig, axx = plt.subplots(nrows=2, ncols=4, figsize=(12,4),sharex=True, sharey=True) if ax is None else (ax.figure, ax)
        plt.subplots_adjust(hspace=0, wspace=0)
        fig.suptitle(f'cell #{icell}')
        for ell,ax in zip(ll, axx.flatten()):
            ell.plot(ax=ax, ylabel='')
            ax.axvline(1, color='grey', ls='--')
            
    def plot_livetime_per_band(self, tmax=2e4,   **kwargs):
        """Scatter plot of the livetime, in ks, for each energy band. Plots are capped at {tmax} s. 
        Points with zero are red.
        """
        
        lt = np.array([cell['livetime'] for cell in self])
        df = pd.DataFrame(lt)
        fig, axx = plt.subplots(nrows=8, figsize=(12,8), sharex=True, sharey=True)
        plt.subplots_adjust(top=0.94, hspace=0)
        fig.suptitle(f'Live time per band ({self.interval}-day interval)')
        t = self.cells.time.values
        kw = dict(xlabel='MJD')
        kw.update(kwargs)
        for i, ax in enumerate(axx):
            x = df.iloc[:, i].clip(0, tmax)*1e-3
            ax.set(**kw)
            ax.plot(t, x, '.');
            ax.plot(t[x==0], x[x==0], '.r')
        return fig
            
    def plot_fluxes(self,  **kwargs):
        """Flux light curves for each energy band.
        """
        from wtlike.lightcurve import flux_plot
        flux_table = self.flux_table()

        fig, axx = plt.subplots(nrows =8, figsize=(15,12), sharex=True, sharey=True)
        fig.suptitle(f'Energy-band lightcurves for {self.source_name}')
        plt.subplots_adjust(hspace=0, top=0.95)
        width = np.diff(flux_table.index)[0]
        for iband, ax in enumerate(axx.flatten()):
            lc = pd.DataFrame.from_dict(dict(t=flux_table.index, tw=width, fit=flux_table.iloc[:,iband]))
            ax.axhline(1, color='grey')
            flux_plot(lc, ax=ax, source_name=f'Energy {iband}', axhline=dict(color='grey'), **kwargs)
        return fig
    
    def plot_poisson_parameters(self, **kwargs):
        """Scatter plot of the width and background poisson parameters.
        """

        ft = self.flux_table()
        f1 =  np.vectorize(lambda x: x.p[1])
        f2 =  np.vectorize(lambda x: x.p[2])

        z =   ft.T #to_numpy().T

        fig, ax = plt.subplots(figsize=(10,6))
        kw = dict(xlim=(0,None), xlabel='sqrt(count equivalence)',
               ylim=(-0.1,1.25), ylabel='relative background',
              title='Poisson parameters')
        kw.update(kwargs)
        ax.set(**kw)
        for i, (x,y) in enumerate(zip(f1(z),f2(z))):
            t = np.sqrt(np.where(x>0, x, 0) )
            ax.plot(t, y, '.', label=f'band {i}');
        ax.legend();ax.grid();
        return fig


# %% ../nbs/11_energyflux.ipynb 22
def pull_plot(pulls, ax=None, nbins=25, **kwargs):
    """ Histogram of a set of "pulls", presumed to be normal (0,1) distributed, 
    the function for which is overplotted.
    """
    from scipy import stats
    x = pulls[~np.isnan(pulls)]
    
    fig, ax = plt.subplots(1,1, figsize=(5,3)) if ax is None else (ax.figure, ax)
    if len(x)<2:
        ax.text(0.5, 0.5, '(no data)', color='red',
                transform=ax.transAxes, ha='center', va='center', )
        ax.set(xticks=[], yticks=[])
        return fig
    
    kw = dict(ylim=(0.,0.45))
    kw.update(kwargs)
    ax.set(**kwargs)
    ax.hist(pulls, np.linspace(-5,5,nbins+1), histtype='step', lw=2, density=True, )
    dom = np.linspace(-5,5,101)
    ax.plot(dom, stats.norm.pdf(dom), '--', lw=2,)# label='Normal(0,1)')
    ax.grid(alpha=0.5)
    
    ax.text(0.95, 0.97,f'{len(x)}\n{np.nanmean(x):.3f}\n'+r'$\pm$'+f'{np.nanstd(x):.3f}',
            transform=ax.transAxes, ha='right', va='top',  
            fontdict=dict(family='monospace'), fontsize=8 )
    return fig

class FluxFixer():

    def __init__(self, fluxfits, emom):
        self.fluxfits = fluxfits
        fluxes = fluxfits.apply(lambda x: x.flux)
        x, y = emom, fluxes-1
        self.x=x; self.y = y
        
        class LSQ():
            def __init__(self, x, y):
                self.a, self.b = np.linalg.lstsq(np.vstack([x, np.ones(len(x))]).T, y, rcond=None)[0]
            def __call__(self, x):
                return self.a*x + self.b
            def __repr__(self):
                return rf'$y={self.a:.2f}x + {self.b:.3f}$'
        self.fitfun = LSQ(x,y)
        
    def scatter_plot(self, ax=None, xlim=np.array([-0.12,0.2]), ylim=(-1,1)):
        
        fig, ax = plt.subplots(figsize=(10,3)) if ax is None else (ax.figure, ax)
        
        ax.plot(xlim, self.fitfun(xlim), '--', color='orange', label=f'least-squares fit:\n{self.fitfun}')
        ax.set(ylim=ylim, xlim=xlim);
        ax.grid(alpha=0.5)
        ax.scatter(self.x, self.y, marker='.',); 
        ax.legend(loc='upper right')
        ax.axhline(0, color='grey')
        
    def pull_std(self):
        return self.fluxfits.apply(lambda x: x.sigma_dev(1)).std()
        
    def pull_hist(self, ax=None, fixit=False, **kwargs):
        
        if fixit:
            fits = self.fluxfits.values
            adjust = self.fitfun(self.x)
            pulls = np.array([fit.sigma_dev(yy+1) for fit, yy in zip(fits, adjust)])
        else:
            pulls = self.fluxfits.apply(lambda x: x.sigma_dev(1))
            
        fig, ax = plt.subplots(figsize=(2,1.5)) if ax is None else (ax.figure, ax)
        kw=dict(yticks=[], ylim=(0,0.5))
        kw.update(kwargs)
        pull_plot(pulls, ax=ax, **kw)
        
    def plot_fix(self):

        fig, (ax1,ax2, ax3) = plt.subplots(ncols=3, figsize=(8,1.5), 
                                           gridspec_kw=dict(width_ratios=[3,1,1]))
        plt.subplots_adjust(wspace=0.1)
        self.scatter_plot(ax=ax1)
        self.pull_hist(ax=ax2)
        self.pull_hist(ax=ax3, fixit=True)
 
