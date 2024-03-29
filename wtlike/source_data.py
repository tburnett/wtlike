# AUTOGENERATED! DO NOT EDIT! File to edit: ../nbs/05_source_data.ipynb.

# %% auto 0
__all__ = ['SourceData']

# %% ../nbs/05_source_data.ipynb 3
import os, sys
import numpy as np
import pandas as pd
from .config import *
from .sources import PointSource
from .load_data import load_source_data, binned_exposure
from .simulation import *

# %% ../nbs/05_source_data.ipynb 4
class SourceData(object):
    """ Load the photon data near the source and the associated exposure.
    --or--
    Use a Simulation object to generate data

    Either from:
      1. `config.wtlike_data/'data_files'`, the Path to folder with list of pickle files
      2. the cache, with key `{source.name}_data`

    * source : name, PointSource, or Simulation
    * `config` : basic configuration
    * `source` : PointSource object if specified
    * `clear` : if set, overwrite the cached results

    Calculate the values for

    * S, B : sums of w and 1-w
    * exptot : total associated exposure

    """

    def __init__(self, source, config=None,  clear=False,
                 week_range=None, key=''):
        """

        """

        self.config = config if config else Config()
        self.verbose = self.config.verbose
        self.simulated=False
        self.used_key = None

        # if self.verbose>0:
        #     print(f'SourceData: week_range: {week_range}')
        ## source is either a name, a PointSource object, or a Simulation
        if type(source)==str:

            try:
                self.source = PointSource(source, config=self.config)
            except Exception as e:
                print(f'{e}', file=sys.stderr)
                raise
            # if a string, use it as the name
            self.source_name = source

        elif isinstance(source, PointSource):
            self.source = source # do I need this?
            self.source_name = source.name

        elif isinstance(source, Simulation):
            self.simulated=True
            self.source=None
            self.source_name = getattr(source, 'source_name', '')
            # can put this into cache
            source.run()
            self.photons = source.photons
            self.exposure = source.exposure

        else: # assume really PointSource
            self.source = source # do I need this?
            self.source_name = source.name

        # not sure why
        # if self.source is not None:
        #     key = f'{self.source.filename}_data' if key=='' else key
        #     self.source.data_key = key
        # else: # no cache for sim, yet
        #     key=None


        if not self.simulated:
            # either load from data, or from a chache--also key used to retrieve data
            if self.config.verbose>1:
                print(f'Loading source data, week_range={week_range}, key={key}')

            ret =load_source_data( self.config, self.source, week_range, key, clear)
            self.photons, self.exposure = ret[:2]
            if len(ret)>2: self.used_key = ret[2]


        else: #TODO
            pass

        # make range of MJD or days available
        self.start = self.exposure.start[0]
        self.stop =  self.exposure.stop.values[-1]
        self.exptot = self.exposure.exp.sum()

        # estimates for signal and background counts in total exposure
        w = self.photons.weight
        self.S = np.sum(w)
        self.B = np.sum(1-w)

        if self.verbose>0:
            print(SourceData.__repr__(self))

    def rates(self):
        print(f'Average fluxes for {self.source_name}: signal {self.S/self.exptot:.2e}/s, background {self.B/self.exptot:.2e}/s')

    def __repr__(self):
        time = self.photons.time.values

        exp = self.exposure
        days  = np.sum(exp.stop-exp.start); secs = days*24*3600
        exp_text = f' average effective area {self.exptot/secs:.0f} cm^2 for {secs/1e6:.1f} Ms'

        if not self.simulated:
            photon_text = f'photons from {UTC(time[0])[:10]} to {UTC(time[-1])[:10]}'
        else:
            photon_text = f'simulated photons over {days:.1f} days.'

        r = f'SourceData: Source {self.source_name} with:'\
            f'\n\t data:     {len(self.photons):9,} {photon_text}'\
            f'\n\t exposure: {len(self.exposure):9,} intervals, {exp_text}'

        self.src_flux, self.bkg_flux = self.S/self.exptot,  self.B/self.exptot
        r+= f'\n\t rates:  source {self.src_flux:.2e}/s, background {self.bkg_flux:.2e}/s,'
        if not self.simulated and 'ts' in self.source.fit_info :
            r+= f' TS {self.source.fit_info["ts"]:.1f}'
#             f' S/N ratio {self.src_flux/self.bkg_flux:.2e}'

        return r

    def binned_exposure(self, time_edges):
        """Bin the exposure

        - time_bins: list of edges.
        """
        return binned_exposure(self.config, self.exposure,  time_edges)

    def binned_cos_theta(self, time_bins=None):
        """ Calculate average cosine of angle with respect to bore axis, per time bin
        """
        if time_bins is None:
            time_bins = get_default_bins(self.config, self.exposure)
        df = self.exposure.copy()
        estop =df.stop.values
        df.loc[:,'tbin'] =np.digitize(estop, time_bins)
        ct = df.groupby('tbin').mean()['cos_theta']
        return ct, time_bins

    def weight_histogram(self, nbins=1000, key=''):
        """ return a weight distribution
        """
        def doit(nbins):
            return np.histogram(self.p_df.weight.values, np.linspace(0,1,nbins+1))[0]

        key = f'{self.source_name}_weight_hist' if key=='' else key
        description = f'Weight histogram for {self.source_name}' if self.config.verbose>0 else ''
        return self.config.cache(key, doit, nbins, description=description)

    def plot_data(self):
        import matplotlib.pyplot as plt
        if self.simulated:
            print(f'Simulated!')
            fig, (ax1, ax4) = plt.subplots(1,2, figsize=(8,4))
            ax1.hist(self.photons.time.values, 500, histtype='step');
            ax1.set(xlabel='Time (MJD)')

            ax4.hist(self.photons.weight, 100, histtype='step')
            ax4.set(xlabel='weight')

        else:
            fig, (ax1,ax2, ax3,ax4) = plt.subplots(1,4, figsize=(15,4))
            ax1.hist(self.photons.time.values, 100, histtype='step');
            ax1.set(xlabel='Time (MJD)')
            ax2.hist(self.photons.radius.values**2, 100, histtype='step', log=True);
            ax2.set(xlabel='Radius**2 (deg**2)', ylim=(100, None));

            ax3.hist(self.photons.band, 32, histtype='step', log=True);
            ax3.set(xlabel='Band index')
            ax4.hist(self.photons.weight, 100, histtype='step')
            ax4.set(xlabel='weight')
