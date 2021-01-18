# The `wtlike` package
> Code for generating fermi-LAT light curves.


### GitHub Links

- [this document](https://tburnett.github.io/wtlike/)
-  [repository](https://github.com/tburnett/wtlike)

## Context

This package has code that was adapted to the [nbdev](https://nbdev.fast.ai/) code/tests/documentation environment from the [github package lat-timing](https://github.com/tburnett/lat-timing) to manage light curves of Fermi-LAT sources.  
It is based on a [paper](https://arxiv.org/pdf/1910.00140.pdf) by Matthew Kerr, which derives the [weighted likelihood formalism](/loglike.html#The-Kerr-likelihood-formula) used here, specifically with
the [Bayesian Block](https://arxiv.org/pdf/1207.5578.pdf) to detect and characterize variability of a gamma-ray source.

Also, I've ported some code from  my [jupydoc](https://github.com/tburnett/jupydoc) documentation package supporting enhanced documentation combining Markdown and code, such that the 
Markdown reflects execution of the code.

## Installation
Currently in pre-alpha, and must be cloned. This brings up the `nbdev` stuff as well.



## Module summary

### Primary ingredients:

- **Photon data**: time, energy, and position. This is managed as a [Parquet](https://parquet.apache.org/) database, derived from the public _Fermi_-LAT data. Modules are [`photon_data`](/wtlike/photon_data) and `data_man`.

- **Exposure**: Need to normalize expected rates as a function of time, position, and energy. Modules: [`exposure`](/wtlike/exposure.html), [`effective_area`](/wtlike/effective_area.html), and [`gti`](/wtlike/gti.html).

- **Weights**: for each photon, this is the probability that it was from a specific source. This requires an analysis of all the souces in the region about the source, involving `gtlike` or `pointlike`. Module `weights` finds the data for a the source, and adds it to the photon list.

### The light-curve 
The module `cells` creates a list of "cells", time intervals with a list of the photon weights, and an estimation from the exposure about the expected distribution. 
Each cell is fit using the modules `loglike` and `poisson`.
Finally, `lightcurve` makes a list of the fits.

### Simulation
A light curve can be also generated with a simulation, see `simulation`.

### Bayesian Blocks
The `bayesian` module creates a partition, and refits the blocks.




