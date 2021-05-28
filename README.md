# Making FermiLAT gamma-ray Light curves with wtlike 
> Quickly create a light curve for any 4FGL source, on any time scale, with optional Bayesian Block analysis 


#### github links:  
[this document](https://tburnett.github.io/wtlike/),   [the repository](https://github.com/tburnett/wtlike)

## Introduction
`wtlike`(Pronounced "DUB-Tee-like"), is a library optimized for interactive exploration in a [Jupyter notebook](https://jupyter.org/) with access to all Fermi-LAT data, and to analyze the time dependence of any source in the
4FGL catalog on any time scale, with the option of performing a [Bayesian Block](https://arxiv.org/pdf/1207.5578.pdf)  partition to select optimal time intervals. The source can be identified by the 4FGL name, or any equivalent common name.

Here is a minimal demo:

```python
from wtlike import *
wtl = WtLike('3C 273', clear=True)
```

    SourceData: photons and exposure for 3C 273: Saving to cache with key "3C 273_data"
    	Assembling photon data and exposure for source 3C 273 from folder "/home/burnett/wtlike_data/data_files",
    	 with 669 files, last file:  week_678.pkl: loading all files
    ....................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................

    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py:204: RuntimeWarning: Mean of empty slice.
      b = np.array(healpy.pix2ang(nside, shifted_pix, nest=True, lonlat=True)).mean(axis=1).round(1)
    /home/burnett/miniconda3/lib/python3.7/site-packages/numpy/core/_methods.py:163: RuntimeWarning: invalid value encountered in true_divide
      ret, rcount, out=ret, casting='unsafe', subok=False)



    ---------------------------------------------------------------------------

    Exception                                 Traceback (most recent call last)

    <ipython-input-2-482849d3de1f> in <module>
          1 from wtlike import *
    ----> 2 wtl = WtLike('3C 273', clear=True)
    

    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/lightcurve.py in __init__(self, *pars, **kwargs)
        335         self.n_min = kwargs.pop('n_min', 2)
        336         self.lc_key = kwargs.pop('lc_key', None)
    --> 337         super().__init__(*pars, **kwargs)
        338         self.update()
        339 


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/cell_data.py in __init__(self, *pars, **kwargs)
         30         bins = kwargs.pop('bins', kwargs.pop('time_bins', Config().time_bins))
         31         #  load source data
    ---> 32         super().__init__(*pars, **kwargs )
         33 
         34 


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/source_data.py in __init__(self, source, config, clear, week_range, key)
        422                             _load_from_weekly_data, self.config, self.source, week_range,
        423                             overwrite=clear,
    --> 424                             description=f'SourceData: photons and exposure for {self.source_name}')
        425             photons, self.exposure = r[:2]
        426             self.runs = r[2] if len(r)==3 else None


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/config.py in __call__(self, key, func, description, overwrite, *pars, **kwargs)
        143 
        144         if ret is None or overwrite:
    --> 145             ret = func(*pars, **kwargs)
        146             self.add(key, ret, exist_ok=overwrite)
        147         return ret


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/source_data.py in _load_from_weekly_data(config, source, week_range)
        331         exposure = _calculate_exposure_for_source(config, source, week )
        332         if photons is not None:
    --> 333             add_weights(config, photons, source)
        334             pp.append(photons)
        335             runs.append( add_exposure_to_events(config, exposure, photons,)


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py in add_weights(config, photon_data, source)
        117     ## NEW
        118     wtman = WeightMan(config, filename=weight_file)
    --> 119     photon_data = wtman.add_weights(photon_data)
        120 
        121     ## OLD


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py in add_weights(self, photons)
        269         """
        270         if self.format==0:
    --> 271             self._old_format(photons)
        272         else:
        273             photons = self._new_format(photons)


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py in _old_format(self, photons)
        204             b = np.array(healpy.pix2ang(nside, shifted_pix, nest=True, lonlat=True)).mean(axis=1).round(1)
        205 
    --> 206             raise Exception(f'There was no overlap of the photon data at {b} and the weights at {a}')
        207         shifted_pix[bad] = 12*nside**2 # set index to be beyond pixel indices
        208 


    Exception: There was no overlap of the photon data at [nan nan] and the weights at [290.   64.3]


```python
wtl.plot();
```

    SourceData: photons and exposure for 3C 273: Saving to cache with key "3C 273_data"
    	Assembling photon data and exposure for source 3C 273 from folder "/home/burnett/wtlike_data/data_files",
    	 with 669 files, last file:  week_678.pkl: loading all files
    ....................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................................

    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py:204: RuntimeWarning: Mean of empty slice.
      b = np.array(healpy.pix2ang(nside, shifted_pix, nest=True, lonlat=True)).mean(axis=1).round(1)
    /home/burnett/miniconda3/lib/python3.7/site-packages/numpy/core/_methods.py:163: RuntimeWarning: invalid value encountered in true_divide
      ret, rcount, out=ret, casting='unsafe', subok=False)



    ---------------------------------------------------------------------------

    Exception                                 Traceback (most recent call last)

    <ipython-input-5-73bc6ad3422a> in <module>
          1 from wtlike import *
    ----> 2 wtl = WtLike('3C 273', clear=True)
          3 wtl.plot();


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/lightcurve.py in __init__(self, *pars, **kwargs)
        335         self.n_min = kwargs.pop('n_min', 2)
        336         self.lc_key = kwargs.pop('lc_key', None)
    --> 337         super().__init__(*pars, **kwargs)
        338         self.update()
        339 


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/cell_data.py in __init__(self, *pars, **kwargs)
         30         bins = kwargs.pop('bins', kwargs.pop('time_bins', Config().time_bins))
         31         #  load source data
    ---> 32         super().__init__(*pars, **kwargs )
         33 
         34 


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/source_data.py in __init__(self, source, config, clear, week_range, key)
        422                             _load_from_weekly_data, self.config, self.source, week_range,
        423                             overwrite=clear,
    --> 424                             description=f'SourceData: photons and exposure for {self.source_name}')
        425             photons, self.exposure = r[:2]
        426             self.runs = r[2] if len(r)==3 else None


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/config.py in __call__(self, key, func, description, overwrite, *pars, **kwargs)
        143 
        144         if ret is None or overwrite:
    --> 145             ret = func(*pars, **kwargs)
        146             self.add(key, ret, exist_ok=overwrite)
        147         return ret


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/source_data.py in _load_from_weekly_data(config, source, week_range)
        331         exposure = _calculate_exposure_for_source(config, source, week )
        332         if photons is not None:
    --> 333             add_weights(config, photons, source)
        334             pp.append(photons)
        335             runs.append( add_exposure_to_events(config, exposure, photons,)


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py in add_weights(config, photon_data, source)
        117     ## NEW
        118     wtman = WeightMan(config, filename=weight_file)
    --> 119     photon_data = wtman.add_weights(photon_data)
        120 
        121     ## OLD


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py in add_weights(self, photons)
        269         """
        270         if self.format==0:
    --> 271             self._old_format(photons)
        272         else:
        273             photons = self._new_format(photons)


    /mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py in _old_format(self, photons)
        204             b = np.array(healpy.pix2ang(nside, shifted_pix, nest=True, lonlat=True)).mean(axis=1).round(1)
        205 
    --> 206             raise Exception(f'There was no overlap of the photon data at {b} and the weights at {a}')
        207         shifted_pix[bad] = 12*nside**2 # set index to be beyond pixel indices
        208 


    Exception: There was no overlap of the photon data at [nan nan] and the weights at [290.   64.3]


```python
debug
```

    > [0;32m/mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py[0m(206)[0;36m_old_format[0;34m()[0m
    [0;32m    204 [0;31m            [0mb[0m [0;34m=[0m [0mnp[0m[0;34m.[0m[0marray[0m[0;34m([0m[0mhealpy[0m[0;34m.[0m[0mpix2ang[0m[0;34m([0m[0mnside[0m[0;34m,[0m [0mshifted_pix[0m[0;34m,[0m [0mnest[0m[0;34m=[0m[0;32mTrue[0m[0;34m,[0m [0mlonlat[0m[0;34m=[0m[0;32mTrue[0m[0;34m)[0m[0;34m)[0m[0;34m.[0m[0mmean[0m[0;34m([0m[0maxis[0m[0;34m=[0m[0;36m1[0m[0;34m)[0m[0;34m.[0m[0mround[0m[0;34m([0m[0;36m1[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    205 [0;31m[0;34m[0m[0m
    [0m[0;32m--> 206 [0;31m            [0;32mraise[0m [0mException[0m[0;34m([0m[0;34mf'There was no overlap of the photon data at {b} and the weights at {a}'[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    207 [0;31m        [0mshifted_pix[0m[0;34m[[0m[0mbad[0m[0;34m][0m [0;34m=[0m [0;36m12[0m[0;34m*[0m[0mnside[0m[0;34m**[0m[0;36m2[0m [0;31m# set index to be beyond pixel indices[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    208 [0;31m[0;34m[0m[0m
    [0m


    ipdb>  ll


    [1;32m    190 [0m    [0;32mdef[0m [0m_old_format[0m[0;34m([0m[0mself[0m[0;34m,[0m [0mphotons[0m[0;34m)[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    191 [0m        [0;32mif[0m [0;32mnot[0m [0mself[0m[0;34m.[0m[0mconfig[0m[0;34m.[0m[0mnest[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    192 [0m            [0;31m# data are RING[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    193 [0m            [0mphoton_pix[0m [0;34m=[0m [0mhealpy[0m[0;34m.[0m[0mring2nest[0m[0;34m([0m[0mconfig[0m[0;34m.[0m[0mnside[0m[0;34m,[0m [0mphotons[0m[0;34m.[0m[0mpixel[0m[0;34m.[0m[0mvalues[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    194 [0m        [0;32melse[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    195 [0m            [0mphoton_pix[0m [0;34m=[0m [0mphotons[0m[0;34m.[0m[0mpixel[0m[0;34m.[0m[0mvalues[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    196 [0m        [0mnside[0m [0;34m=[0m [0mself[0m[0;34m.[0m[0mnside_wt[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    197 [0m        [0mto_shift[0m [0;34m=[0m [0;36m2[0m[0;34m*[0m[0mint[0m[0;34m([0m[0mnp[0m[0;34m.[0m[0mlog2[0m[0;34m([0m[0mself[0m[0;34m.[0m[0mconfig[0m[0;34m.[0m[0mnside[0m[0;34m//[0m[0mself[0m[0;34m.[0m[0mnside_wt[0m[0;34m)[0m[0;34m)[0m[0;34m;[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    198 [0m        [0mshifted_pix[0m [0;34m=[0m   [0mnp[0m[0;34m.[0m[0mright_shift[0m[0;34m([0m[0mphoton_pix[0m[0;34m,[0m [0mto_shift[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    199 [0m        [0mbad[0m [0;34m=[0m [0mnp[0m[0;34m.[0m[0mlogical_not[0m[0;34m([0m[0mnp[0m[0;34m.[0m[0misin[0m[0;34m([0m[0mshifted_pix[0m[0;34m,[0m [0mself[0m[0;34m.[0m[0mwt_pix[0m[0;34m)[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    200 [0m        [0;32mif[0m [0mself[0m[0;34m.[0m[0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m>[0m[0;36m0[0m [0;34m&[0m [0msum[0m[0;34m([0m[0mbad[0m[0;34m)[0m[0;34m>[0m[0;36m0[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    201 [0m            [0mprint[0m[0;34m([0m[0;34mf'\tApplying weights: {sum(bad)} / {len(bad)} photon pixels are outside weight region'[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    202 [0m        [0;32mif[0m [0msum[0m[0;34m([0m[0mbad[0m[0;34m)[0m[0;34m==[0m[0mlen[0m[0;34m([0m[0mbad[0m[0;34m)[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    203 [0m            [0ma[0m [0;34m=[0m [0mnp[0m[0;34m.[0m[0marray[0m[0;34m([0m[0mhealpy[0m[0;34m.[0m[0mpix2ang[0m[0;34m([0m[0mnside[0m[0;34m,[0m [0mself[0m[0;34m.[0m[0mwt_pix[0m[0;34m,[0m [0mnest[0m[0;34m=[0m[0;32mTrue[0m[0;34m,[0m [0mlonlat[0m[0;34m=[0m[0;32mTrue[0m[0;34m)[0m[0;34m)[0m[0;34m.[0m[0mmean[0m[0;34m([0m[0maxis[0m[0;34m=[0m[0;36m1[0m[0;34m)[0m[0;34m.[0m[0mround[0m[0;34m([0m[0;36m1[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    204 [0m            [0mb[0m [0;34m=[0m [0mnp[0m[0;34m.[0m[0marray[0m[0;34m([0m[0mhealpy[0m[0;34m.[0m[0mpix2ang[0m[0;34m([0m[0mnside[0m[0;34m,[0m [0mshifted_pix[0m[0;34m,[0m [0mnest[0m[0;34m=[0m[0;32mTrue[0m[0;34m,[0m [0mlonlat[0m[0;34m=[0m[0;32mTrue[0m[0;34m)[0m[0;34m)[0m[0;34m.[0m[0mmean[0m[0;34m([0m[0maxis[0m[0;34m=[0m[0;36m1[0m[0;34m)[0m[0;34m.[0m[0mround[0m[0;34m([0m[0;36m1[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    205 [0m[0;34m[0m[0m
    [0;32m--> 206 [0;31m            [0;32mraise[0m [0mException[0m[0;34m([0m[0;34mf'There was no overlap of the photon data at {b} and the weights at {a}'[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[1;32m    207 [0m        [0mshifted_pix[0m[0;34m[[0m[0mbad[0m[0;34m][0m [0;34m=[0m [0;36m12[0m[0;34m*[0m[0mnside[0m[0;34m**[0m[0;36m2[0m [0;31m# set index to be beyond pixel indices[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    208 [0m[0;34m[0m[0m
    [1;32m    209 [0m        [0;31m# find indices with search and add a "weights" column[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    210 [0m        [0;31m# (expect that wt_pix are NEST ordering and sorted)[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    211 [0m        [0mweight_index[0m [0;34m=[0m [0mnp[0m[0;34m.[0m[0msearchsorted[0m[0;34m([0m[0mself[0m[0;34m.[0m[0mwt_pix[0m[0;34m,[0m[0mshifted_pix[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    212 [0m        [0mband_index[0m [0;34m=[0m [0mnp[0m[0;34m.[0m[0mfmin[0m[0;34m([0m[0;36m31[0m[0;34m,[0m [0mphotons[0m[0;34m.[0m[0mband[0m[0;34m.[0m[0mvalues[0m[0;34m)[0m [0;31m#all above 1 TeV into last bin[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    213 [0m[0;34m[0m[0m
    [1;32m    214 [0m        [0;31m# final grand lookup -- isn't numpy wonderful![0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    215 [0m        [0mphotons[0m[0;34m.[0m[0mloc[0m[0;34m[[0m[0;34m:[0m[0;34m,[0m[0;34m'weight'[0m[0;34m][0m [0;34m=[0m [0mself[0m[0;34m.[0m[0mwts[0m[0;34m[[0m[0mtuple[0m[0;34m([0m[0;34m[[0m[0mband_index[0m[0;34m,[0m [0mweight_index[0m[0;34m][0m[0;34m)[0m[0;34m][0m[0;34m[0m[0;34m[0m[0m
    [1;32m    216 [0m[0;34m[0m[0m
    


    ipdb>  p to_shift


    8


    ipdb>  p photon_pix


    array([], dtype=int32)


    ipdb>  p config.nside


    *** NameError: name 'config' is not defined


    ipdb>  p self.config.nside


    1024


    ipdb>  p self.config.nest


    True


    ipdb>  p len(photons)


    0


    ipdb>  up


    > [0;32m/mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py[0m(271)[0;36madd_weights[0;34m()[0m
    [0;32m    269 [0;31m        """
    [0m[0;32m    270 [0;31m        [0;32mif[0m [0mself[0m[0;34m.[0m[0mformat[0m[0;34m==[0m[0;36m0[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m--> 271 [0;31m            [0mself[0m[0;34m.[0m[0m_old_format[0m[0;34m([0m[0mphotons[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    272 [0;31m        [0;32melse[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    273 [0;31m            [0mphotons[0m [0;34m=[0m [0mself[0m[0;34m.[0m[0m_new_format[0m[0;34m([0m[0mphotons[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m


    ipdb>  p photons


    Empty DataFrame
    Columns: [band, time, pixel, radius]
    Index: []


    ipdb>  ll


    [1;32m    264 [0m    [0;32mdef[0m [0madd_weights[0m[0;34m([0m[0mself[0m[0;34m,[0m [0mphotons[0m[0;34m)[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    265 [0m        """
    [1;32m    266 [0m        [0mget[0m [0mthe[0m [0mphoton[0m [0mpixel[0m [0mids[0m[0;34m,[0m [0mconvert[0m [0mto[0m [0mNEST[0m [0;34m([0m[0;32mif[0m [0;32mnot[0m [0malready[0m[0;34m)[0m [0;32mand[0m [0mright[0m [0mshift[0m [0mthem[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    267 [0m        [0madd[0m [0;34m'weight'[0m[0;34m,[0m [0mremove[0m [0;34m'band'[0m[0;34m,[0m [0;34m'pixel'[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    268 [0m[0;34m[0m[0m
    [1;32m    269 [0m        """
    [1;32m    270 [0m        [0;32mif[0m [0mself[0m[0;34m.[0m[0mformat[0m[0;34m==[0m[0;36m0[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [0;32m--> 271 [0;31m            [0mself[0m[0;34m.[0m[0m_old_format[0m[0;34m([0m[0mphotons[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[1;32m    272 [0m        [0;32melse[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    273 [0m            [0mphotons[0m [0;34m=[0m [0mself[0m[0;34m.[0m[0m_new_format[0m[0;34m([0m[0mphotons[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    274 [0m[0;34m[0m[0m
    [1;32m    275 [0m        [0;31m# don't need these columns now (add flag to config to control??)[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    276 [0m        [0mphotons[0m[0;34m.[0m[0mdrop[0m[0;34m([0m[0;34m[[0m[0;34m'pixel'[0m[0;34m][0m[0;34m,[0m [0maxis[0m[0;34m=[0m[0;36m1[0m[0;34m,[0m [0minplace[0m[0;34m=[0m[0;32mTrue[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    277 [0m[0;34m[0m[0m
    [1;32m    278 [0m        [0;32mif[0m [0mself[0m[0;34m.[0m[0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m>[0m[0;36m1[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    279 [0m            [0mprint[0m[0;34m([0m[0;34mf'\t{sum(np.isnan(photons.weight.values)):,} events without weight'[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    280 [0m        [0;32mreturn[0m [0mphotons[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    281 [0m[0;34m[0m[0m
    


    ipdb>  up


    > [0;32m/mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/weights.py[0m(119)[0;36madd_weights[0;34m()[0m
    [0;32m    117 [0;31m    [0;31m## NEW[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    118 [0;31m    [0mwtman[0m [0;34m=[0m [0mWeightMan[0m[0;34m([0m[0mconfig[0m[0;34m,[0m [0mfilename[0m[0;34m=[0m[0mweight_file[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m--> 119 [0;31m    [0mphoton_data[0m [0;34m=[0m [0mwtman[0m[0;34m.[0m[0madd_weights[0m[0;34m([0m[0mphoton_data[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    120 [0;31m[0;34m[0m[0m
    [0m[0;32m    121 [0;31m    [0;31m## OLD[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [0m


    ipdb>  p photon_data


    Empty DataFrame
    Columns: [band, time, pixel, radius]
    Index: []


    ipdb>  ll


    [1;32m    104 [0m[0;32mdef[0m [0madd_weights[0m[0;34m([0m[0mconfig[0m[0;34m,[0m  [0mphoton_data[0m[0;34m,[0m [0msource[0m[0;34m)[0m[0;34m:[0m [0;31m# nbins=50):[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    105 [0m    """ add weights for the source to the photon data
    [1;32m    106 [0m[0;34m[0m[0m
    [1;32m    107 [0m    [0;34m-[0m [0mphoton_data[0m [0;34m-[0m[0;34m-[0m [0mDataFrame[0m [0;32mwith[0m [0mphoton[0m [0mdata[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    108 [0m[0;34m[0m[0m
    [1;32m    109 [0m    [0;34m-[0m [0msource[0m [0;34m-[0m[0;34m-[0m[0;31m [0m[0;31m`[0m[0mPointSource[0m[0;31m`[0m [0mobject[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    110 [0m[0;34m[0m[0m
    [1;32m    111 [0m    [0mReturn[0m [0mthe[0m [0mweight[0m [0mvalue[0m [0mhistogram[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    112 [0m    """
    [1;32m    113 [0m    [0mweight_file[0m [0;34m=[0m  [0mcheck_weights[0m[0;34m([0m[0mconfig[0m[0;34m,[0m  [0msource[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    114 [0m    [0;32mif[0m [0mweight_file[0m [0;32mis[0m [0;32mNone[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    115 [0m        [0;32mraise[0m [0mException[0m[0;34m([0m[0;34mf'Weight file not found for {source}'[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    116 [0m[0;34m[0m[0m
    [1;32m    117 [0m    [0;31m## NEW[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    118 [0m    [0mwtman[0m [0;34m=[0m [0mWeightMan[0m[0;34m([0m[0mconfig[0m[0;34m,[0m [0mfilename[0m[0;34m=[0m[0mweight_file[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0;32m--> 119 [0;31m    [0mphoton_data[0m [0;34m=[0m [0mwtman[0m[0;34m.[0m[0madd_weights[0m[0;34m([0m[0mphoton_data[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[1;32m    120 [0m[0;34m[0m[0m
    


    ipdb>  up


    > [0;32m/mnt/c/users/thbur/OneDrive/work/wtlike/wtlike/source_data.py[0m(333)[0;36m_load_from_weekly_data[0;34m()[0m
    [0;32m    331 [0;31m        [0mexposure[0m [0;34m=[0m [0m_calculate_exposure_for_source[0m[0;34m([0m[0mconfig[0m[0;34m,[0m [0msource[0m[0;34m,[0m [0mweek[0m [0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    332 [0;31m        [0;32mif[0m [0mphotons[0m [0;32mis[0m [0;32mnot[0m [0;32mNone[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m--> 333 [0;31m            [0madd_weights[0m[0;34m([0m[0mconfig[0m[0;34m,[0m [0mphotons[0m[0;34m,[0m [0msource[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    334 [0;31m            [0mpp[0m[0;34m.[0m[0mappend[0m[0;34m([0m[0mphotons[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[0;32m    335 [0;31m            runs.append( add_exposure_to_events(config, exposure, photons,)
    [0m


    ipdb>  ll


    [1;32m    289 [0m[0;32mdef[0m [0m_load_from_weekly_data[0m[0;34m([0m[0mconfig[0m[0;34m,[0m [0msource[0m[0;34m,[0m [0mweek_range[0m[0;34m=[0m[0;32mNone[0m[0;34m)[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    290 [0m    """
    [1;32m    291 [0m    [0mGenerate[0m [0mcombinded[0m [0mDataFrames[0m [0;32mfrom[0m [0ma[0m [0mlist[0m [0mof[0m [0mpickled[0m [0mfiles[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    292 [0m    [0mEither[0m [0mweekly[0m [0;32mor[0m [0mmonthly[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    293 [0m[0;34m[0m[0m
    [1;32m    294 [0m    [0mkwargs[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    295 [0m    [0;34m-[0m [0mweek_range[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    296 [0m    """
    [1;32m    297 [0m[0;34m[0m[0m
    [1;32m    298 [0m    [0;31m# check weights[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    299 [0m    [0mweight_file[0m [0;34m=[0m  [0mcheck_weights[0m[0;34m([0m[0mconfig[0m[0;34m,[0m  [0msource[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    300 [0m    [0;32massert[0m [0mweight_file[0m [0;32mis[0m [0;32mnot[0m [0;32mNone[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    301 [0m[0;34m[0m[0m
    [1;32m    302 [0m    [0mdata_folder[0m [0;34m=[0m [0mconfig[0m[0;34m.[0m[0mwtlike_data[0m[0;34m/[0m[0;34m'data_files'[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    303 [0m    [0mdata_files[0m [0;34m=[0m [0msorted[0m[0;34m([0m[0mlist[0m[0;34m([0m[0mdata_folder[0m[0;34m.[0m[0mglob[0m[0;34m([0m[0;34m'*.pkl'[0m[0;34m)[0m[0;34m)[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    304 [0m    [0miname[0m [0;34m=[0m [0mdata_folder[0m[0;34m.[0m[0mname[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    305 [0m[0;34m[0m[0m
    [1;32m    306 [0m    [0;32mif[0m [0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m>[0m[0;36m0[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    307 [0m        print(f"\tAssembling photon data and exposure for source {source.name} from"\
    [1;32m    308 [0m              [0;34mf' folder "{data_folder}",\n\t with {len(data_files)} files,'[0m[0;31m\[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    309 [0m              f' last file:  {data_files[-1].name}: ', end='')
    [1;32m    310 [0m[0;34m[0m[0m
    [1;32m    311 [0m    [0mw1[0m[0;34m,[0m[0mw2[0m [0;34m=[0m [0mweek_range[0m [0;32mor[0m  [0mconfig[0m[0;34m.[0m[0mweek_range[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    312 [0m    [0;32mif[0m [0mw1[0m [0;32mis[0m [0;32mnot[0m [0;32mNone[0m [0;32mor[0m [0mw2[0m [0;32mis[0m [0;32mnot[0m [0;32mNone[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    313 [0m        [0;32mif[0m [0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m>[0m[0;36m0[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    314 [0m            [0mprint[0m[0;34m([0m[0;34mf'\tLoading weeks {w1}:{w2}'[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    315 [0m        [0mdata_files[0m[0;34m=[0m [0mdata_files[0m[0;34m[[0m[0mw1[0m[0;34m:[0m[0mw2[0m[0;34m][0m[0;34m[0m[0;34m[0m[0m
    [1;32m    316 [0m    [0;32melse[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    317 [0m        [0;32mif[0m [0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m>[0m[0;36m0[0m[0;34m:[0m [0mprint[0m[0;34m([0m[0;34m'loading all files'[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    318 [0m[0;34m[0m[0m
    [1;32m    319 [0m    [0mverbose[0m[0;34m,[0m [0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m=[0m[0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m,[0m [0;36m0[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    320 [0m[0;34m[0m[0m
    [1;32m    321 [0m    [0;31m# list of data framees[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    322 [0m    [0mpp[0m [0;34m=[0m [0;34m[[0m[0;34m][0m[0;34m[0m[0;34m[0m[0m
    [1;32m    323 [0m    [0mee[0m [0;34m=[0m [0;34m[[0m[0;34m][0m[0;34m[0m[0;34m[0m[0m
    [1;32m    324 [0m    [0mruns[0m [0;34m=[0m [0;34m[[0m[0;34m][0m[0;34m[0m[0;34m[0m[0m
    [1;32m    325 [0m    [0;32mfor[0m [0mf[0m [0;32min[0m [0mdata_files[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    326 [0m        [0mprint[0m[0;34m([0m[0;34m'.'[0m[0;34m,[0m [0mend[0m[0;34m=[0m[0;34m''[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    327 [0m        [0;32mwith[0m [0mopen[0m[0;34m([0m[0mf[0m[0;34m,[0m [0;34m'rb'[0m[0;34m)[0m [0;32mas[0m [0minp[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    328 [0m            [0mweek[0m [0;34m=[0m [0mpickle[0m[0;34m.[0m[0mload[0m[0;34m([0m[0minp[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    329 [0m[0;34m[0m[0m
    [1;32m    330 [0m        [0mphotons[0m [0;34m=[0m [0m_get_photons_near_source[0m[0;34m([0m[0mconfig[0m[0;34m,[0m [0msource[0m[0;34m,[0m [0mweek[0m [0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    331 [0m        [0mexposure[0m [0;34m=[0m [0m_calculate_exposure_for_source[0m[0;34m([0m[0mconfig[0m[0;34m,[0m [0msource[0m[0;34m,[0m [0mweek[0m [0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    332 [0m        [0;32mif[0m [0mphotons[0m [0;32mis[0m [0;32mnot[0m [0;32mNone[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [0;32m--> 333 [0;31m            [0madd_weights[0m[0;34m([0m[0mconfig[0m[0;34m,[0m [0mphotons[0m[0;34m,[0m [0msource[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [0m[1;32m    334 [0m            [0mpp[0m[0;34m.[0m[0mappend[0m[0;34m([0m[0mphotons[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    335 [0m            runs.append( add_exposure_to_events(config, exposure, photons,)
    [1;32m    336 [0m                       )
    [1;32m    337 [0m        [0mee[0m[0;34m.[0m[0mappend[0m[0;34m([0m[0mexposure[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    338 [0m[0;34m[0m[0m
    [1;32m    339 [0m    [0mprint[0m[0;34m([0m[0;34m''[0m[0;34m)[0m[0;34m;[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    340 [0m    [0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m=[0m[0mverbose[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    341 [0m    [0;31m# concatenate the two lists of DataFrames[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    342 [0m    [0mp_df[0m [0;34m=[0m [0mpd[0m[0;34m.[0m[0mconcat[0m[0;34m([0m[0mpp[0m[0;34m,[0m [0mignore_index[0m[0;34m=[0m[0;32mTrue[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    343 [0m    [0me_df[0m [0;34m=[0m [0mpd[0m[0;34m.[0m[0mconcat[0m[0;34m([0m[0mee[0m[0;34m,[0m [0mignore_index[0m[0;34m=[0m[0;32mTrue[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    344 [0m[0;34m[0m[0m
    [1;32m    345 [0m    [0;31m# process exposure to find runs, and add tau column to photons[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    346 [0m    [0;31m### runs = add_exposure_to_events(config, e_df, p_df)[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    347 [0m[0;34m[0m[0m
    [1;32m    348 [0m    [0;32mif[0m [0mconfig[0m[0;34m.[0m[0mverbose[0m[0;34m>[0m[0;36m1[0m[0;34m:[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    349 [0m        [0mtimes[0m [0;34m=[0m [0mp_df[0m[0;34m.[0m[0mtime[0m[0;34m.[0m[0mvalues[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    350 [0m        [0mprint[0m[0;34m([0m[0;34mf'Loaded {len(p_df):,} photons from {UTC(times[0])} to  {UTC(times[-1])} '[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    351 [0m        [0mprint[0m[0;34m([0m[0;34mf'Calculated {len(e_df):,} exposure entries'[0m[0;34m)[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    352 [0m[0;34m[0m[0m
    [1;32m    353 [0m    [0;31m# add weights to photon data[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    354 [0m    [0;31m### add_weights(config, p_df, source)[0m[0;34m[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    355 [0m[0;34m[0m[0m
    [1;32m    356 [0m    [0;32mreturn[0m [0mp_df[0m[0;34m,[0m [0me_df[0m[0;34m,[0m [0mruns[0m[0;34m[0m[0;34m[0m[0m
    [1;32m    357 [0m[0;34m[0m[0m
    


    ipdb>  p photons


    Empty DataFrame
    Columns: [band, time, pixel, radius]
    Index: []


    ipdb>  p week


    {'tstart': 560912295.0, 'photons':         band  nest_index       time
    0          4    12295332       7.32
    1          3    11914943       8.45
    2          9    10811538       8.62
    3          0    10525286      10.21
    4          1    10574926      10.95
    ...      ...         ...        ...
    322690     0     4917176  604775.19
    322691     0     4813633  604775.38
    322692     1     4806884  604775.81
    322693     0    12451316  604776.62
    322694     1     5017066  604776.75
    
    [322695 rows x 3 columns], 'sc_data':           start      stop  livetime  ra_scz  dec_scz  ra_zenith  dec_zenith
    0      58402.04  58402.04     27.32  357.45   -43.93     329.87       -0.37
    1      58402.04  58402.04     27.39  358.79   -43.21     331.57        0.45
    2      58402.04  58402.04     27.44    0.14   -42.49     333.28        1.26
    3      58402.04  58402.04     27.43    1.47   -41.79     334.99        2.08
    4      58402.04  58402.04     27.48    2.80   -41.09     336.70        2.89
    ...         ...       ...       ...     ...      ...        ...         ...
    17104  58409.04  58409.04     27.56  303.58   -47.22     274.56       -3.93
    17105  58409.04  58409.04     27.51  304.96   -46.47     276.27       -3.12
    17106  58409.04  58409.04     27.53  306.34   -45.72     277.98       -2.31
    17107  58409.04  58409.04     27.52  307.71   -44.98     279.68       -1.50
    17108  58409.04  58409.04     18.80  309.06   -44.24     281.39       -0.69
    
    [17109 rows x 7 columns], 'gti_times': array([5.60912296e+08, 5.60918009e+08, 5.60918023e+08, 5.60923714e+08,
           5.60923728e+08, 5.60929420e+08, 5.60929434e+08, 5.60935125e+08,
           5.60935139e+08, 5.60940831e+08, 5.60940845e+08, 5.60943458e+08,
           5.60944087e+08, 5.60949234e+08, 5.60950613e+08, 5.60955087e+08,
           5.60956701e+08, 5.60960951e+08, 5.60962724e+08, 5.60966819e+08,
           5.60968698e+08, 5.60972993e+08, 5.60974657e+08, 5.60979147e+08,
           5.60980610e+08, 5.60985222e+08, 5.60986559e+08, 5.60991268e+08,
           5.60992406e+08, 5.60997308e+08, 5.60997787e+08, 5.61003590e+08,
           5.61003604e+08, 5.61009295e+08, 5.61009309e+08, 5.61015001e+08,
           5.61015015e+08, 5.61020706e+08, 5.61020720e+08, 5.61026411e+08,
           5.61026425e+08, 5.61029039e+08, 5.61029950e+08, 5.61034840e+08,
           5.61036270e+08, 5.61040696e+08, 5.61042342e+08, 5.61046560e+08,
           5.61048354e+08, 5.61052442e+08, 5.61054325e+08, 5.61058662e+08,
           5.61060282e+08, 5.61064797e+08, 5.61066234e+08, 5.61070864e+08,
           5.61072183e+08, 5.61076908e+08, 5.61077979e+08, 5.61082949e+08,
           5.61083275e+08, 5.61089170e+08, 5.61089184e+08, 5.61094876e+08,
           5.61094890e+08, 5.61100581e+08, 5.61100595e+08, 5.61106287e+08,
           5.61106301e+08, 5.61111992e+08, 5.61112006e+08, 5.61114620e+08,
           5.61115703e+08, 5.61120447e+08, 5.61121922e+08, 5.61126305e+08,
           5.61127980e+08, 5.61132170e+08, 5.61133983e+08, 5.61138089e+08,
           5.61139951e+08, 5.61144327e+08, 5.61145907e+08, 5.61150445e+08,
           5.61151858e+08, 5.61156505e+08, 5.61157798e+08, 5.61162548e+08,
           5.61163526e+08, 5.61168616e+08, 5.61168755e+08, 5.61174751e+08,
           5.61174765e+08, 5.61180456e+08, 5.61180470e+08, 5.61186162e+08,
           5.61186176e+08, 5.61191867e+08, 5.61191881e+08, 5.61197573e+08,
           5.61197587e+08, 5.61200214e+08, 5.61201414e+08, 5.61206053e+08,
           5.61207571e+08, 5.61211915e+08, 5.61213617e+08, 5.61217779e+08,
           5.61219611e+08, 5.61223769e+08, 5.61225576e+08, 5.61229988e+08,
           5.61231531e+08, 5.61236091e+08, 5.61237482e+08, 5.61242146e+08,
           5.61243414e+08, 5.61248188e+08, 5.61249040e+08, 5.61254626e+08,
           5.61254640e+08, 5.61260331e+08, 5.61260345e+08, 5.61266037e+08,
           5.61266051e+08, 5.61271742e+08, 5.61271756e+08, 5.61277448e+08,
           5.61277462e+08, 5.61283153e+08, 5.61283167e+08, 5.61285812e+08,
           5.61287093e+08, 5.61291661e+08, 5.61293217e+08, 5.61297524e+08,
           5.61299254e+08, 5.61303389e+08, 5.61305239e+08, 5.61309454e+08,
           5.61311202e+08, 5.61315645e+08, 5.61317156e+08, 5.61321736e+08,
           5.61323106e+08, 5.61327787e+08, 5.61329024e+08, 5.61333828e+08,
           5.61334552e+08, 5.61340206e+08, 5.61340220e+08, 5.61345912e+08,
           5.61345926e+08, 5.61351617e+08, 5.61351631e+08, 5.61357322e+08,
           5.61357336e+08, 5.61363028e+08, 5.61363042e+08, 5.61365732e+08,
           5.61365908e+08, 5.61368733e+08, 5.61368747e+08, 5.61371419e+08,
           5.61372759e+08, 5.61377270e+08, 5.61378860e+08, 5.61383133e+08,
           5.61384889e+08, 5.61389000e+08, 5.61390866e+08, 5.61395132e+08,
           5.61396827e+08, 5.61401299e+08, 5.61402780e+08, 5.61407380e+08,
           5.61408730e+08, 5.61413427e+08, 5.61414605e+08, 5.61419468e+08,
           5.61420043e+08, 5.61425786e+08, 5.61425800e+08, 5.61431492e+08,
           5.61431506e+08, 5.61437197e+08, 5.61437211e+08, 5.61442902e+08,
           5.61442916e+08, 5.61448608e+08, 5.61448622e+08, 5.61451240e+08,
           5.61451992e+08, 5.61457025e+08, 5.61458420e+08, 5.61462879e+08,
           5.61464501e+08, 5.61468743e+08, 5.61470520e+08, 5.61474611e+08,
           5.61476493e+08, 5.61480805e+08, 5.61482451e+08, 5.61486951e+08,
           5.61488404e+08, 5.61493022e+08, 5.61494353e+08, 5.61499067e+08,
           5.61500185e+08, 5.61505108e+08, 5.61505532e+08, 5.61511366e+08,
           5.61511380e+08, 5.61517072e+08])}


    ipdb>  q


This assumes that the name for the source, in this case the historically famous first [quasar](https://en.wikipedia.org/wiki/Quasar#Background) to be discovered, can be associated with a 4FGL catalog source. The plot shows, as a function of the MJD time, weekly measurements of deviations of the flux relative to the average of the 12-year interval used to define the 4FGL-DR3 catalog.

The first stage, extracting data for the source, takes ~10 min, but, using an included [cache system](https://tburnett.github.io/wtlike/config.html#Cache), only has to be done once.

## Overview

This package has code that was developed with the [nbdev](https://nbdev.fast.ai/) code/tests/documentation environment from the [github package lat-timing](https://github.com/tburnett/lat-timing) to generate light curves of Fermi-LAT sources.  
It is based on a [paper](https://arxiv.org/pdf/1910.00140.pdf) by Matthew Kerr, which derives the [weighted likelihood formalism](https://tburnett.github.io/wtlike/loglike#The-Kerr-likelihood-formula) used here, specifically with
the [Bayesian Block](https://arxiv.org/pdf/1207.5578.pdf) to detect and characterize variability of a gamma-ray source.

There are several innovative design features that significantly improve the speed and portability.

* Condensed photon and spacecraft data. 
* Weight tables
* A cache to improve interactivity
* Likelihood functions fit to a modified Poisson function
* Unbinned likelihood 
* A simple user interface


### How it works

Historically, gamma-ray source measurements have used two methods:
1. For a fixed time interval, where the energy and 
position are used, with a model including all potential sources and a model of the detector response to 
define the likelihood. This has been the only way to study weak sources. A light curve must apply this to
each time interval.
2. Or, for very bright flares, for example GRBs, one can simply count the number of photons within a
circular region, that is, aperture photometry.

Matthew Kerr [introduced](https://arxiv.org/pdf/1910.00140.pdf) a third method, basically counting photons but using information from a static
likelihood analysis to assign a "weight" to each photon, the probability for being from the source in question, then optimizing this likelihood. This assumes that the only thing changing is the flux of
the source.  He calls it "retrospective", since the analysis for the full time is then applied back to the individual photons.

### Individual photon Likelihood ("unbinned")
We use  a version of the Kerr likelihood where the fundamental "cell" is a single photon. The likelihood for any group of photons 
is easily determined by adding the photons.

Assumptions:
* Source photons are completely contained in the dataset boundaries (the "ROI").
* The instrument response is constant with time.
* The background pattern is constant. (Clearly violated if a neighboring source varies!)

For a photon $i$ with weight $w_i$ and exposure $\tau_i$,

{% raw %}
$$ \displaystyle\log\mathcal{L}_i(\alpha) = \log (1 + \alpha \ w_i ) - \alpha \ w_i\ R\ \tau_i $$
{% endraw %}

where:
* $\alpha$ -- fractional deviation of the source flux from the average as measured by the catalog analysis 
* $w_i$ -- probability, for the nominal source rate, that photon $i$ is from the source.
* $R$ -- expected source rate in $\mathrm{cm^{-2}\ s^{-1}}$ units. 
* $\tau_i$ -- integration of the exposure rate for the live time preceding this detection in $\mathrm{cm^2} s$ units. 
Live time, the time since the start of a run, or the previous photon. It is equivalent to time.  
This behaves like a time, but accounts for variation of the exposure rate, often rapid with respect to event rates.

(A note about "exposure": It has units $\mathrm{cm^2\ s}$ and is the integral of the "exposure rate" over a time interval.
For *Fermi*, the rate is typically $\mathrm{3000 cm^2}$ but varies as much as a factor of three over a single orbit.)

This is evaluated in the module  [loglike](https://tburnett.github.io/wtlike/loglike).

### Photon Data

*Fermi* data are retrieved from the Fermi-LAT weekly data files extracted from the [GSFC FTP server](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly), 
with subfolders for the photon data, `photon` and spacecraft data, `spacecraft`. It is [described here](http://fermi.gsfc.nasa.gov/ssc/data/access/http://fermi.gsfc.nasa.gov/ssc/data/access/). The files are organized into individual weeks, starting UTC midnight each Thursday. Files for the most recent week are updated daily.

We convert each pair of photon and spacecraft files to two DataFrame tables with the minimum information needed to evaluate the likelihood, as compressed as possible. Particularly, the energy and position are binned. Details can be seen in the module [data_man](https://tburnett.github.io/wtlike/data_man).  

The entire data set in this format occupies <2 GB.

### Select Data for a source

All further processing uses a subset of the photons within a cone, currently $4^\circ$, about the selected source, and 
evaluates the exposure during 30-s intervals for this direction.  In the class
[`SourceData`](https://tburnett.github.io/wtlike/source_data#SourceData), implemented in the module 
[`source_data`](https://tburnett.github.io/wtlike/source_data) we
1. Extract photons
2. Evaluate exposure, using the effective area tables and a spectral assumption.
3. Determine the weight for each photon, using the table for the source. See the module [weights](https://tburnett.github.io/wtlike/weights) for details.
4. For each photon's livetime, determine the exposure $\tau$. 

The result is a photon DataFrame, containing for each photon, the time $t$ in MJD units, $w$,  and $\tau$.

This class is a superclass of the user interface class `WtLike` introduced above.

### Partition into cells

A "cell", the terminology used by Kerr, the set of photons in a specific time interval. The class 
[CellData](https://tburnett.github.io/wtlike/cell_data#CellData)
in the module [cell_data](https://tburnett.github.io/wtlike/cell_data), a subclass of SourceData manages this.

This class is instantiated with a tuple to define the binning in time. Denoted by `(a, b, c)`, the elements are:
* `a`--start time 
* `b`--stop time
* `c`-- bin size in days, but 0 means orbit-based, intervals are contiguous eposure.

For the start and stop, values > 50000 are interpreted as MJD. Otherwise they are relative to start if positive
or stop if negative for the full dataset, both rounded to a full day. Zero means actual start for `a` and stop for `b`.
The default binning is simply `(0, 0, 7)` for weekly bins with the full dataset. Hourly for the most recent week would be `(-7, 0, 1/24)`.

A DataFrame table of the cells is created as a data member `cells`, with content

* `t` -- MJD time at cell center
* `tw` -- cell width in days
* `e` -- cell exposure, for reference 
* `n` -- the number of photons
* `w` -- a list of `n` weights
* `S` -- expected number source photons, the nominal source rate times the sum of $\tau$ values.
* `B` -- expected number of background photons (unused)

### Views

`CellData` implements a method `view`, which accepts a binning tuple as an argument, returns a *copy* of the current object, which can be a subclass, assigning to it the binning. Thus the view has all the attributes of its parent, but
with a different set of cells. 

So the following creates a new WtLike object that we generated above, rebins a copy with 25-day bins in the first 100 days, generates a list of the cells, then removes it since it wasn't assigned a reference variable.


```python
wtl.view(0,100,25).cells
```

    CellData: Bin photon data into 4 4-week bins from 54683.0 to 54783.0
    LightCurve: select 4 cells for fitting with e>1 & n>2





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>t</th>
      <th>tw</th>
      <th>e</th>
      <th>n</th>
      <th>w</th>
      <th>S</th>
      <th>B</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>54695.5</td>
      <td>25.0</td>
      <td>1573.86</td>
      <td>653</td>
      <td>[0.48898458, 0.6999376, 0.11305099, 0.19166616...</td>
      <td>184.38</td>
      <td>541.88</td>
    </tr>
    <tr>
      <th>1</th>
      <td>54720.5</td>
      <td>25.0</td>
      <td>1974.69</td>
      <td>1575</td>
      <td>[0.434471, 0.052078784, 0.05434912, 0.60623187...</td>
      <td>231.33</td>
      <td>679.88</td>
    </tr>
    <tr>
      <th>2</th>
      <td>54745.5</td>
      <td>25.0</td>
      <td>1533.28</td>
      <td>1283</td>
      <td>[0.26372972, 0.31009516, 0.83747405, 0.0952251...</td>
      <td>179.62</td>
      <td>527.90</td>
    </tr>
    <tr>
      <th>3</th>
      <td>54770.5</td>
      <td>25.0</td>
      <td>2039.27</td>
      <td>1328</td>
      <td>[0.0763341, 0.40165585, 0.07539207, 0.59247845...</td>
      <td>238.90</td>
      <td>702.12</td>
    </tr>
  </tbody>
</table>
</div>



### Evaluate Likelihoods and make light curve plots

The class [`LightCurve`](https://tburnett.github.io/wtlike/lightcurve), implemented in the module [`lightcurve`](https://tburnett.github.io/wtlike#LightCurve) is a subclass of `SourceData`.
An instance invokes its superclass to generate the set of cells, then evaluates the likelihoods.

#### Poisson-like Likelihood
We fit the likelihood for each cell, using only a few evaluations, to a [3-parameter Poisson-like formula](https://tburnett.github.io/poisson). Two advantages of this are:
* efficiency -- for large numbers of photons, this is much faster
* convenience -- the [`Poisson`](https://tburnett.github.io/wtlike/poisson.html#Poisson) object implements functions that return the TS, 95% CL, maximum, and uncertainties, using the [incomplete gamma function](https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.gammainc.html).

This results in a DataFrame `fits` containing
* `t, tw, n`, from the cell, and
* `fit`, the `Poisson` object

#### Light curve plots

The function `Likelihood.plot` actually invokes [`flux_plot`](https://tburnett.github.io/lightcurve#flux_plot)

### Apply Bayesian Blocks
The class [`WtLike`](https://tburnett.github.io/wtlike/90-main.html#WtLike), 
implememted in the module [`main`](https://tburnett.github.io/wtlike/90-main.html), adds an implementation of the Bayesian Block (BB) algorithm, from the module [`bayesian`](https://tburnett.github.io/wtlike)  using likelihood instead of counts.
There we have two subclasses of `astropy.stats.bayesian_blocks`, [`CountFitness`](https://tburnett.github.io/wtlike/bayesian.html#CountFitness) and the default [`LikelihoodFitness`](https://tburnett.github.io/wtlike/bayesian.html#LikelihoodFitness).
    
This code creates partitions between boundaries of a set of cells. Usage is via a special view, 
[bb_view`](https://tburnett.github.io/wtlike/main#WtLike.bb_view)
                     

```python
bb = wtl.bb_view()
bb.plot();
```

    LightCurve: select 656 cells for fitting with e>1 & n>2
    Bayesian Blocks: using penalty 0.05
    Partitioned 656 cells into 91 blocks, using LikelihoodFitness 
    LightCurve: Loaded 91 / 91 cells for fitting



![png](docs/images/output_14_1.png)


As you see, this made 91 blocks from the 656 weeks, fit each, and overplotted in on the weekly light curve.

### Simulation

Finally, a simulation option is available. See the [tutorial](https://tburnett.github.io/wtlike/tutorial/)  for an example

## Installation

Note that this is in beta mode. 
It requires:  matplotlib pandas scipy astropy healpy

To install from pyPI:

```
pip install wtlike
```
Data requirements: There are three sets of files:

- **photon data**<br> 
These are a set of weekly pickled python `dict` objects with compressed condensed photon and spacecraft data extracted from the GSFC FTP site. They contain every photon above 100 MeV, and less than $100^\circ$ 
from the zenith.

- **weight tables**<br>
Each source to be analyzed needs a table defining the photon weight as a function of position, energy, and event type. These are currently generated by pointlike. (A `fermipy`-generated version would be preferable.)

- **effective area**<br>
A standard *Fermi* instrument response file (IRF) defining the effective area as a function of detector angle and energy. 

A set of these is available as a 1.6 GB zip file.

## Input data

There are three data sources which `wtlike` needs to function:


-	The photon/spacecraft data
-	A table of weights for each source
-	An effective area IRF table 

These must be found under a folder, which by default is `~/wtlike_data`. In that folder there must be (perhaps links to) three folders named `data_files`, `weight_files`, `aeff_files`.  A copy of what I'm using is at `/afs/slac/g/glast/users/burnett/wtlike_data`
