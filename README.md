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
Currently in pre-alpha, and must be cloned. This includes the `nbdev` support as well. 

## Demo

The following code cell loads the data for the BL Lac blazar, and plots by default, a daily light curve for the full *fermi* mission

```python
from wtlike import *
wtl = WtLike('BL Lac', bins=(0,0,7)) # how to define 7-day bins for the full dataset. Default is (0,0,1).
wtl.plot_flux();
```

    photons and exposure for BL Lac: Restoring from cache with key "BL Lac__data"
    BBanalysis: Source BL Lac with:
    	 data:       404,789 photons from   2008-08-04 to 2021-04-28
    	 exposure: 3,171,944 intervals from 2008-08-04 to 2021-04-28
    Bin photon data into 664 1-week bins from 54683.0 to 59331.0
    Loaded 656 / 656 cells with exposure > 1e-06 for light curve analysis



![png](docs/images/output_2_1.png)


The variable `wtl` has lots of capabilities.
To examine a subset of the data at the end of the current data, we create a new `WtLike` object and plot it.

```python
wtl2 = wtl.rebinned_copy((-5,0, 1/24)) # for the last 5 days, 1-hour bins
wtl2.plot_flux(fmt='o'); # Accepts plt.plot args, e.g. xlim, ylim, etc.
```

    Bin photon data into 120 1-hour bins from 59328.0 to 59333.0
    The new range is within old--bin exposure factor 127.4
    Loaded 107 / 107 cells with exposure > 1e-06 for light curve analysis



![png](docs/images/output_4_1.png)


Or, to do a Bayesian Block partiton with these 1-hour bins, perform fits, and overplot the result, just run the following.

```python
wtl2.plot_BB(fmt='o');
```

    BL Lac__bb_edges: Restoring from cache
    Partitioned 107 cells into 8 blocks, using LikelihoodFitness 
    Loaded 8 / 8 cells for fitting



![png](docs/images/output_6_1.png)


Finally, let's look at the values plotted above:

```python
wtl2.bb_table()
```




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
      <th>n</th>
      <th>flux</th>
      <th>ts</th>
      <th>errors</th>
      <th>limit</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>59328.73</td>
      <td>1.46</td>
      <td>312</td>
      <td>0.33</td>
      <td>89.9</td>
      <td>(-0.05, 0.051)</td>
      <td>0.42</td>
    </tr>
    <tr>
      <th>1</th>
      <td>59330.02</td>
      <td>1.12</td>
      <td>362</td>
      <td>0.78</td>
      <td>236.7</td>
      <td>(-0.079, 0.082)</td>
      <td>0.92</td>
    </tr>
    <tr>
      <th>2</th>
      <td>59330.83</td>
      <td>0.50</td>
      <td>290</td>
      <td>1.83</td>
      <td>439.6</td>
      <td>(-0.155, 0.162)</td>
      <td>2.10</td>
    </tr>
    <tr>
      <th>3</th>
      <td>59331.17</td>
      <td>0.17</td>
      <td>229</td>
      <td>7.32</td>
      <td>727.3</td>
      <td>(-0.563, 0.588)</td>
      <td>8.33</td>
    </tr>
    <tr>
      <th>4</th>
      <td>59331.31</td>
      <td>0.12</td>
      <td>103</td>
      <td>3.67</td>
      <td>242.3</td>
      <td>(-0.445, 0.475)</td>
      <td>4.50</td>
    </tr>
    <tr>
      <th>5</th>
      <td>59331.67</td>
      <td>0.58</td>
      <td>236</td>
      <td>1.01</td>
      <td>200.0</td>
      <td>(-0.116, 0.121)</td>
      <td>1.22</td>
    </tr>
    <tr>
      <th>6</th>
      <td>59332.35</td>
      <td>0.79</td>
      <td>213</td>
      <td>0.48</td>
      <td>88.4</td>
      <td>(-0.076, 0.079)</td>
      <td>0.62</td>
    </tr>
    <tr>
      <th>7</th>
      <td>59332.88</td>
      <td>0.25</td>
      <td>101</td>
      <td>1.56</td>
      <td>130.9</td>
      <td>(-0.222, 0.238)</td>
      <td>1.97</td>
    </tr>
  </tbody>
</table>
</div>



## Input data

There are three data sources which `wtlike` needs to function:


-	The photon/spacecraft data
-	A table of weights for each source
-	An effective area IRF table 

These must be found under a folder, which by default is `~/wtlike_data`. In that folder there must be (perhaps links to) three folders named `data_files`, `weight_files`, `aeff_files`.  A copy of what I'm using is at `/afs/slac/g/glast/users/burnett/wtlike_data`

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
THis module defines the class [BBanalysis(LightCurve)](https://tburnett.github.io/wtlike/bayesian#BBanalysis). INheriting from `LightCurve`, it adds Bayesian block capability. It is the class returned by `from wtlike import WtLike'

### Simulation [simulation](https://tburnett.github.io/wtlike/simulation)
A light curve can be also generated with a simulation.
