# The `wtlike` package
> Code for generating fermi-LAT light curves. <br>


## Background

This package has code that is being adapted to the [nbdev](https://nbdev.fast.ai/) environment from [github package lat-timing](https://github.com/tburnett/lat-timing) to manage light curves of Fermi-LAT sources.  

As pointed out there, it is based on a [paper](https://arxiv.org/pdf/1910.00140.pdf) by Matthew Kerr, 

And at the same time, I've ported some code from  my [jupydoc](https://github.com/tburnett/jupydoc) documention package to allow enhanced documentation combining Markdown and code.

## Installation
Currently in pre-alpha, and must be cloned. This brings up the `nbdev` stuff as well.




## Module summary


### `config` -- configuration parameters, time and source  

### `lightcurve` -- The light curve

The light curve can be represented by plots of flux vs. time. The time range,
limited by 'config.mjd_range`, if set. The actual livetime, during which *Fermi* is
collecting data, is further limited by the GTI, for good-time interval. This is a list
of start,stop pairs.

### `gti`

The module [gti](/light_curves/lgti.html) defines `get_gti`.

During the valid times, a the flux, or rate, is estimated by counting the number 
of photons and dividing by the exposure.

The source is defined by instantiating a [PointSource(/light_curves/config#PointSource) object, defined in 


### `exposure` -- the exposure

The exposure for the specified source is calculated from the  [exposure](/light_curves/exposure.html) module,
which implements `get_exposure`. It depends on:

- Space craft info (FT2)
The FT2 file(s) contain spacecraft position and orientation as a function of time.

### `effective_area` -- Efffective Area
The module [effective_area](light_curves/effective_area.html) defines the functor class
`EffectiveArea`, needed to calculate the exposure


### `photon_data`

### `weights`



### `cells`

A "cell" represents a time interval to measure the flux.

### `bayesian`

### `simulation`




Dependencies:
```
- simulation: loglike, exposure, lightcurve

- bayesian: lightcurve, cells

- lightcurve: loglike, cells

- loglike:  poisson

- cells: photon_data, weigths, exposure

- weights: photon_data

- photon_data: gti

- exposure: gti, effective_area
```

