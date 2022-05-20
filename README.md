# Making FermiLAT gamma-ray light curves with wtlike 
> Quickly create a light curve for any 4FGL source, on any time scale, with optional Bayesian Block analysis 


#### github links:  
[this document](https://tburnett.github.io/wtlike/),   [the repository](https://github.com/tburnett/wtlike)

## Introduction
`wtlike`(Perhaps pronounced "DUB-Tee-like"), is a library optimized for interactive exploration in a [Jupyter notebook](https://jupyter.org/) with access to all Fermi-LAT data, and to analyze the time dependence of any source in the
4FGL catalog on any time scale, with the option of performing a [Bayesian Block](https://arxiv.org/pdf/1207.5578.pdf)  partition to select optimal time intervals. The source can be identified by the 4FGL name, or any equivalent common name.

Here is a minimal demo:

```python
from wtlike import *
config = Config()
if config.valid:
    with Timer() as t:
        wtl = WtLike('3C 273',week_range=(9,3*52+9))
        print(t)
```

    SourceData: week_range: (9, 165)
    SourceData:  3C 273: Restoring from cache with key "P88Y3157_weeks_9-165"
    SourceData: Source 3C 273 with:
    	 data:        35,895 photons from 2008-08-04 to 2011-08-03
    	 exposure:   713,330 intervals,  average effective area 2783 cm^2 for 21.3 Ms
    	 rates:  source 2.05e-07/s, background 4.00e-07/s, TS 30341.0
    CellData.rebin: Bin photon data into 156 1-week bins from 54683.0 to 55775.0
    LightCurve: select 156 cells for fitting with e>35 & n>2
    elapsed time: 0.8s (0.0 min)


This created a `WtLike` object, loading the first 3 years of data, by specifying weeks from #9, the first data-taking week.
The reason to specify only the first three years here is to avoid the 10 min or so that extracting the furl 13+ years would take for this demo. If that had already been done, then using an included [cache system](https://tburnett.github.io/wtlike/config.html#Cache), it only has to be done once.

Now ask it to make a plot:

```python
if config.valid: 
    wtl.plot(UTC=1);
```


    
![png](docs/images/output_3_0.png)
    


This assumes that the name for the source, in this case the historically famous first [quasar](https://en.wikipedia.org/wiki/Quasar#Background) to be discovered, can be associated with a 4FGL catalog source. The plot shows, as a function of UTC (or MJD if desired) time, weekly measurements of deviations of the flux relative to the average of the 12-year interval used to define the 4FGL-DR3 catalog.



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
2. Or, for very bright flares, for example GRBs or AGN flares, one can simply count the number of photons within a
circular region, that is, aperture photometry.

Matthew Kerr [introduced](https://arxiv.org/pdf/1910.00140.pdf) a third method, basically counting photons but using information from a static
likelihood analysis to assign a "weight" to each photon, the probability for being from the source in question, then optimizing this likelihood. This assumes that the only thing changing is the flux of
the source.  He calls it "retrospective", since the analysis for the full time is then applied back to the individual photons. Another way of looking at it is to make the assumption that the time dependence of a source's photon flux factorizes according to the energy and spatitial dependences.  

### Likelihood evaluation

In a significant modification from Kerr's implemetation as described in that paper, we evaluate weights for each photon by a table lookup.

Assumptions:
* Source photons are completely contained in the dataset boundaries (the "ROI").
* The instrument response is constant with time.
* The background pattern is constant. (Clearly violated if a neighboring source varies!)

The likelihood evaluation is implemented in the module  [loglike](https://tburnett.github.io/wtlike/loglike).

### Photon Data

*Fermi* data are retrieved from the Fermi-LAT weekly data files extracted from the [GSFC FTP server](https://heasarc.gsfc.nasa.gov/FTP/fermi/data/lat/weekly), 
with subfolders for the photon data, `photon` and spacecraft data, `spacecraft`. It is [described here](http://fermi.gsfc.nasa.gov/ssc/data/access/http://fermi.gsfc.nasa.gov/ssc/data/access/). The files are organized into individual weeks, starting UTC midnight each Thursday. Files for the most recent week are updated daily.

We convert each pair of photon and spacecraft files to two DataFrame tables with the minimum information needed to evaluate the likelihood, as compressed as possible. Particularly, the energy and position are binned. Details can be seen in the module [data_man](https://tburnett.github.io/wtlike/data_man).  

The entire data set (SOURCE class, energy>100 MeV) in this format occupies ~2 GB.

### Select Data for a source

All further processing uses a subset of the photons within a cone, currently $4^\circ$, about the selected source, and 
evaluates the exposure during 30-s intervals for this direction.  In the class
[`SourceData`](https://tburnett.github.io/wtlike/source_data#SourceData), implemented in the module 
[`source_data`](https://tburnett.github.io/wtlike/source_data) we
1. Extract photons
2. Evaluate exposure, using the effective area tables and a spectral assumption.
3. Determine the weight for each photon, using the table for the source. See the module [weights](https://tburnett.github.io/wtlike/weights) for details.

The result is a photon DataFrame, containing for each photon, the time $t$ in MJD units, $w$.

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
* `S` -- expected number source photons, the nominal source rate times the total exposure.
* `B` -- expected number of background photons (unused)

### Views

`CellData` implements a method `view`, which accepts a binning tuple as an argument, returns a *copy* of the current object, which can be a subclass, assigning to it the binning. Thus the view has all the attributes of its parent, but
with a different set of cells. 

So the following creates a new WtLike object that we generated above, rebins a copy with 25-day bins in the first 100 days, generates a list of the cells, then removes it since it wasn't assigned a reference variable.


```python
if config.valid:
    print(wtl.view(0,100,25).cells)
```

    CellData.rebin: Bin photon data into 4 25-day bins from 54683.0 to 54783.0
    LightCurve: select 4 cells for fitting with e>125 & n>2
             t    tw            e       ctm     n  \
    0  54695.5  25.0  1488.251709  0.670376   616   
    1  54720.5  25.0  1868.208862  0.681914  1523   
    2  54745.5  25.0  1450.841553  0.679041  1281   
    3  54770.5  25.0  1929.609741  0.682303  1258   
    
                                                       w           S           B  
    0  [4.8901367e-01, 0.7885742, 0.11303711, 0.19165...  305.076996  595.216980  
    1  [0.4345703, 0.6064453, 0.0690918, 0.062561035,...  382.964508  747.178467  
    2  [0.33911133, 0.3100586, 0.6225586, 0.06994629,...  297.408295  580.255005  
    3  [0.09112549, 0.58251953, 0.07537842, 0.3457031...  395.551086  771.735352  


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
if config.valid:
    bb = wtl.bb_view()
    bb.plot(UTC=1);
```

    Bayesian Blocks: partitioning 156 cells using LikelihoodFitness with penalty 5%
    	found 43 / 156 blocks.
    LightCurve: Loaded 43 / 43 cells for fitting



    
![png](docs/images/output_14_1.png)
    


As you see, this made 94 blocks from the 656 weeks, fit each, and overplotted in on the weekly light curve.
<details class="description">
    <summary>Code details ...</summary>
    
```python
# collapse-hide
from utilities.ipynb_docgen import *
@ipynb_doc
def pgram(wtl):
    """
    ### Kerr Periodogram

    One can analyze the frequency spectrum by generating a 
    [Kerr periodogram](https://arxiv.org/pdf/1910.00140.pdf). 
    A `WtLike` object has a `periodogram` function that returns a `TimeSeries` object.
    {out}
    {fig}
    """
    with capture_hide('Setup output') as out:
        p =wtl.periodogram()
        p.power_plot(pmax=100 );
    fig = plt.gcf()
    return locals()

if config.valid: pgram(wtl)
```

</details>


### Kerr Periodogram

One can analyze the frequency spectrum by generating a 
[Kerr periodogram](https://arxiv.org/pdf/1910.00140.pdf). 
A `WtLike` object has a `periodogram` function that returns a `TimeSeries` object.
<details  class="nbdoc-description" >  <summary> Setup output </summary>  <div style="margin-left: 5%;"><pre>CellData.rebin: Bin photon data into 26256 1-hour bins from 54683.0 to 55777.0<br>TimeSeries: creating power spectra, 26,256 samples, size 1.00 h: Nyquist is 6.0 /d<br></pre></div> </details>
<figure style="margin-left: 5%" title="Figure 7">   <img src="data:image/png;base64,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
" alt="Figure 7"  <br> </figure>



### Simulation

Finally, a simulation option is available. See the [tutorial](https://tburnett.github.io/wtlike/tutorial/)  for an example

## Installation

Note that this is in beta mode. 
It depends on the packages:  matplotlib, pandas, scipy, astropy, and healpy, which pip should install if needed.

To install from pyPI:

```
pip install wtlike
```

or
```
pip install -U wtlike
```
to upgrade.

## Input data

There are three data sources which `wtlike` needs to function:

-	The photon/spacecraft data, one file per week
-	A table of weights for each source
-	An effective area table, a grid in energy/angle for "front" and "back" events 

These must be found in a folder, which by default is `~/wtlike_data`. In that folder there must be (perhaps links to) three folders named `data_files`, `weight_files`, `aeff_files`, and a zip file, `weight_file.zip`.  

## Google colab setup

An easy way to test this system is to use a Jupyter notebook hosted by Google. 
