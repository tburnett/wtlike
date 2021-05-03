# The `wtlike` package
> Code for generating fermi-LAT light curves.


### GitHub Links

- [this document](https://tburnett.github.io/wtlike/)
-  [repository](https://github.com/tburnett/wtlike)

## Context

This package has code that was adapted to the [nbdev](https://nbdev.fast.ai/) code/tests/documentation environment from the [github package lat-timing](https://github.com/tburnett/lat-timing) to manage light curves of Fermi-LAT sources.  
It is based on a [paper](https://arxiv.org/pdf/1910.00140.pdf) by Matthew Kerr, which derives the [weighted likelihood formalism](https://tburnett.github.io/wtlike/loglike#The-Kerr-likelihood-formula) used here, specifically with
the [Bayesian Block](https://arxiv.org/pdf/1207.5578.pdf) to detect and characterize variability of a gamma-ray source.

Also, I've ported some code from  my [jupydoc](https://github.com/tburnett/jupydoc) documentation package supporting enhanced documentation combining Markdown and code, such that the 
Markdown reflects execution of the code.

## Installation
Currently in pre-alpha, and must be cloned. This brings up the `nbdev` stuff as well. 

## Demo

The following code cell loads the data for the BL Lac blazar, and plots by default, a daily light curve for the full *fermi* mission

```python
from wtlike import *
wtl = WtLike('BL Lac')
wtl.plot_flux();
```

    photons and exposure for BL Lac: Restoring from cache with key "BL Lac__data"
    BBanalysis: Source BL Lac with:
    	 data:       404,789 photons from   2008-08-04 to 2021-04-28
    	 exposure: 3,171,944 intervals from 2008-08-04 to 2021-04-28
    Rebin with 4650 1-day from 54683.0 to 59333.0
    renormalizing exposure to fraction of total: total, sum: 72355436 1
    Loaded 4522 / 4522 cells with exposure > 1e-06 for light curve analysis



![png](docs/images/output_2_1.png)


The variable `wtl` has lots of capabilities.
To examine a subset of the data at the end of the current data, we create a new `WtLike` object and plot it.

```python
wtl2 = wtl.rebinned_copy((-5,0, 1/24)) # for the last 5 days, 1-hour bid
wtl2.plot_flux(fmt='o');
```

    Rebin with 120 1-hour from 59328.0 to 59333.0
    renormalizing exposure to fraction of total: total, sum: 3974348 1
    The new range is within old--bin exposure factor 18.2
    Loaded 107 / 107 cells with exposure > 1e-06 for light curve analysis



![png](docs/images/output_4_1.png)


## Input data

There are three data sources that `wtlike` needs to function:

-	An effective area IRF table 
-	The photon/spacecraft data
-	A table of weights 

These must be found under a folder, which by default is `~/wtlike_data`. In that folder there must be (perhaps links to) three folders named `data_files`, `weight_files`, `aeff_files`. 

## Module summary

### Configuration [config](https://tburnett.github.io/wtlike/config)
Implements basic configuration information, [Config](https://tburnett.github.io/wtlike/config#Config), a cache system [Cache](https://tburnett.github.io/wtlike/config#Cache), point source info [PointSource](https://tburnett.github.io/wtlike/config#PointSource), and [time conversion](https://tburnett.github.io/wtlike/config#Time-conversion)

### Photon and Spacecraft Data  [data_man](https://tburnett.github.io/wtlike/data_man)
This module manages conversion of the weekly FT1 (photons) and FT2 (spacecraft) files, downloaded from  [GSFC](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly), to a folder containing  pickled files, each with tables of photons, space craft data, and a list of GTI times derived from the FT1 file. The total size of this is 2.8 GB. A class [WeeklyData](https://tburnett.github.io/wtlike/data_man#WeeklyData) exports the results.

### Source data  [source_data](https://tburnett.github.io/wtlike/source_data)
The module depends on a specific source. It extracts the photons within a disk, and calculates the exposure for this direction. It assumes that a weigtht analysis has been done for this source, which it uses to apply a weight to each photon. This is handled by the class [SourceData](https://tburnett.github.io/wtlike/source_data#SourceData). 

### Cell data [cell_data](https://tburnett.github.io/wtlike/cell_data)
The next step is to define a set of time bins. This module, implementing the class [CellData(SourceData)](https://tburnett.github.io/wtlike/cell_data#CellData), creates a set of cells.

### The light-curve  [light_curve](https://tburnett.github.io/wtlike/lightcurve)
The the class [LightCurve(CellData)](https://tburnett.github.io/wtlike/lightcurve#LightCurve) uses the set of cell defined by its super class, and evaluates the likelihood for each. This function is represented by a Poisson-like function for further analysis. It creates a table with this information for plotting a light curve.

### Bayesian Blocks [bayesian](https://tburnett.github.io/wtlike/bayesian) 
THis module defines the class [BBanalysis(LightCurve)](https://tburnett.github.io/wtlike/bayesian#BBanalysis). INheriting from `LightCurve`, it adds Bayesian block capability.

### Simulation [simulation](https://tburnett.github.io/wtlike/simulation)
A light curve can be also generated with a simulation.
