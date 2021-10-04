# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/03_exposure.ipynb (unless otherwise specified).

__all__ = ['BaseAeff', 'KerrAeff', 'SourceAeff', 'time_bin_edges', 'sc_data_selection', 'binned_exposure']

# Cell
import pandas as pd
import numpy as np

from .config import (Config, UTC, MJD)
from .effective_area import EffectiveArea

# Internal Cell
def _sc_process(config, source, sc_data):

    """
    - source -- contains ra, dec in degrees
    - sc_data -- DF constructed from spacecraft data (FT2).

    Return: a DF with the S/C data for the source direction, wtih cos theta and zenith cuts

    columns:
    - start, stop, livetime -- from the FT2 info
    - cos_theta -- angle between bore and direction
    """

    # calculate cosines with respect to sky direction
    ra_r,dec_r = np.radians(source.ra), np.radians(source.dec)
    sdec, cdec = np.sin(dec_r), np.cos(dec_r)

    def cosines( ra2, dec2):
        ra2_r =  np.radians(ra2.values)
        dec2_r = np.radians(dec2.values)
        return np.cos(dec2_r)*cdec*np.cos(ra_r-ra2_r) + np.sin(dec2_r)*sdec

    cos_thetas = cosines(sc_data.ra_scz,    sc_data.dec_scz)
    zcosines = cosines(sc_data.ra_zenith, sc_data.dec_zenith)
    # mask out entries too close to zenith, or too far away from ROI center
    mask =   (cos_thetas >= config.cos_theta_max) & (zcosines>=np.cos(np.radians(config.z_max)))
    if config.verbose>1:
        print(f'\tFound {len(mask):,} S/C entries:  {sum(mask):,} remain after zenith and theta cuts')
    dfm = sc_data.loc[mask,:]
    livetime = dfm.livetime.values

    #return livetime, cos_thetas[mask]

    return  pd.DataFrame(
        dict(
            start=sc_data.start[mask],
            stop=sc_data.stop[mask],
            livetime=livetime,
            cos_theta=cos_thetas[mask],
        )
    )

# Cell
# don't use now?
from abc import abstractmethod

class BaseAeff(object):
    """
    Base class for calculating the weighted effective area

    It is a functor of `costh`, returning  the effective area weighted by a spectrum.
    The exposure is then a sum over intervals of this weighted effective area times the livetime per interval.

    A subclass must define a `setup` method to define the spectral function and minimum energy
    to include Back events.
    """

    # default energies to use for Simpson's integral: 8/decade from 100 MeV to 10 GeV
    edom = np.logspace(2,4,41)

    def __init__(self, config, source):

        self.source = source
        self.Aeff = EffectiveArea(file_path=config.datapath/'aeff_files')
        self.setup()

    @abstractmethod
    def setup(self):
        pass

    def __call__(self, cos_theta) : #sc_data):
        """
        Return the weighted effective area as a function of $\\cos( \\theta) $

        """
        # as set by self.setup -- also self.back_min
        edom = self.edom
        wts = self.wts
        cos_theta = np.atleast_1d(cos_theta)

        # a table of the weights for each pair in spectral flux and cos_theta arrays
        rvals = np.empty([len(wts),len(cos_theta)])

        for i, (en,wt) in enumerate(zip(edom,wts)):
            faeff, baeff = np.array(self.Aeff( [en], cos_theta ))
#             if en >= self.back_min:
#                 # both front and back
#                 rvals[i] = (faeff+baeff) * wt
#             else: # front only
#                 rvals[i] = faeff * wt
            rvals[i] =  faeff * wt if en<self.back_min else  (faeff+baeff) * wt

        from scipy.integrate import simpson
        aeff = simpson(rvals, edom,axis=0) / simpson(wts,edom)

        return aeff

class KerrAeff(BaseAeff):
    """
    Kerr implementation, from godot, of the weighted exposure, using a $E^{-2.1}$ powerlaw spectrum
    """

    bins_per_decade: int=10
    base_spectrum: str='lambda E: (E/1000)**-2.1'
    energy_range: tuple =(100,1e4)# (100.,1e6)

    def setup(self ):

        """set up energy domain, evaluate fluxes
           This is the Kerr version, using a wired-in powerlaw spectrum
        """

        emin,emax = self.energy_range
        loge1=np.log10(emin); loge2=np.log10(emax)
        self.edom= np.logspace(loge1, loge2, int((loge2-loge1)*self.bins_per_decade+1))

        spectrum = eval(self.base_spectrum) #lambda E: (E/1000)**-2.1

        self.wts = spectrum(self.edom)

        # the threshold for including Back events
        # should use default, this is what I was usinq
        self.back_min=0

class SourceAeff(BaseAeff):
    """
    BaseAeff subclass that uses the source spectrum applied only to used bands
    """

    def setup(self):
        f = self.source.spectral_model
        self.wts  = f(self.edom)
        self.back_min=10**2.5

# Cell
def time_bin_edges(config, exposure, tbin=None):
    """Return an interleaved array of start/stop values

    tbin: an array (a,b,d), default config.time_bins

    interpretation of a, b:

    if > 50000, interpret as MJD

    if <0, back from stop

    otherwise, offset from start

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
def sc_data_selection(config, source, sc_data):

    """
    Return a DataFrame with the S/C data for the source direction, wtih cos theta and zenith cuts

    columns:
    - start, stop, livetime -- from the FT2 info
    - cos_theta -- angle between bore and direction
    - exp -- the exposure: effective area at angle wighted by a default spectral function, times livetime

    """

    sc_df = _sc_process(config, source, sc_data)
    cos_theta = sc_df.cos_theta.values
    livetime = sc_df.livetime.values

    # now get appropriate weighted effective area, multipy by livetime

    if config.use_kerr:
        sc_df.loc[:,'exp'] = KerrAeff(config, source)(cos_theta) * livetime
    else:
        sc_df.loc[:,'exp'] = SourceAeff(config, source)(cos_theta) * livetime

    return sc_df

# Cell
def binned_exposure(config, exposure, time_edges):
    """Bin the exposure into cells

    - exposure -- A DataFrame derived from FT2
    - time_bins: list of edges, an interleaved start/stop array


    Returns:

    An array of exposure integrated over each time bin. Assumes that the time bins
    are contained within the exposure.

    it is interleaved, client must apply [0::2] selection.

    """

    # get exposure calculation
    exp   =exposure.exp.values
    estart= exposure.start.values
    estop = exposure.stop.values

    # determine bins,

    #use cumulative exposure to integrate over larger periods
    cumexp = np.concatenate(([0],np.cumsum(exp)) )

    # get index into tstop array of the bin edges
    edge_index = np.searchsorted(estop, time_edges)

    # return the exposure integrated over the intervals
    cum = cumexp[edge_index]

    # difference is exposure per interval
    bexp = np.diff(cum)
#     if config.verbose>1:
#         print(f'exposure per bin:\n{pd.Series(bexp).describe(percentiles=[])}')
    return bexp