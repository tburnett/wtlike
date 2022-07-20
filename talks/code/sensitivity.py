import os, glob
import numpy as np
import pylab as plt
import pandas as pd

import healpy
from utilities import healpix as hp
from utilities.ipynb_docgen import *

def histit(zmap, ax=None, label=''):
    fig, ax = plt.subplots(figsize=(4,3)) if ax is None else (ax.figure, ax)
    hkw = dict( bins=np.logspace(-1,1,40),histtype='step', lw=2, label=label)
    ax.hist(zmap,log=False,**hkw);

    ax.set(xscale='log');
    return ax


def eflux_vs_b(lt95, psrs, nside=256):

    # get b for the table
    import healpy
    plt.rc('font', size=14)
    glon, glat = healpy.pix2ang(nside, range(nside*nside*12), lonlat=True)

    tbl=lt95.map #*1e12 #tconv
    good = ~np.isnan(tbl); 
    if sum(~good)>0: print( f'bad points: {sum(~good)}')
    glat_fun = lambda x: np.sign(x)*abs(x)**0.65

    sglat = np.array( [glat_fun(x) for x in glat] ) #sign(glat)*abs(glat)**0.75
    ssb = np.arange(glat_fun(-90), glat_fun(90), 1.0)

    # generate range
    x=[]; y = []; dy = []; low=[]; high=[]
    for i in range(len(ssb)-1):
        cut = (sglat>ssb[i])*(sglat<ssb[i+1])*good
        x.append(sglat[cut].mean())
        subset = tbl[cut] * 1.7 ###########3.25
        y.append( subset.mean())
        dy.append( subset.std())
        low.append(np.percentile(subset,5))
        high.append(np.percentile(subset,95))
    # plt.plot(x,low, '.-');
    # plt.plot(x,high,'+-');
    # plt.gcf().set_facecolor('white');

    G100 = psrs.eflux * 1e12
    # get the glat
    gb = psrs.glat

    fig, ax = plt.subplots(figsize=(8,6))
    if False:
        ax.plot(x, y, '-b')
        ax.plot(x, low, '--b')
        ax.plot(x, high, '--b')
    ax.plot(list(map(glat_fun, gb)), G100, 'o')
    ax.set( yscale='log', xlim=(glat_fun(-80),glat_fun(80)), 
           ylim=(0.2,2000),
        xlabel = r'$\mathsf{ Galactic\ latitude\ (degrees)}$',
        ylabel = r'$\mathsf{ G100\ (10^{-12}\ erg\ cm^{-2}\ s^{-1})}$',)
    xticks = np.array([-60, -30, -10, -5, -2, 0,2, 5, 10, 30, 60])
    ax.set_xticks(glat_fun(xticks))
    ax.set_xticklabels(xticks)
    ax.set_yticks((1,10,100,1000))
    ax.set_yticklabels(['1','10','100', '1000'])
    # 
    xr = np.array(x)[::-1]
    highr = np.array(high)[::-1]
    t=ax.fill(np.hstack([x, xr]), np.hstack([low, highr]), 'grey', label='limit 5 to 95% range')
    t[0].set_alpha(0.5)

    ax.axvline(0, color='grey')
    leg=ax.legend()
    fig.set_facecolor('white');


class FluxModel():
    emin, emax = 1e2, 1e5
    def __init__(self, pars, e0=1000):
        self.pars=pars
        self.e0=e0

    def photon_flux(self):
        from scipy.integrate import quad
        return quad(self, self.emin, self.emax)[0]

    def energy_flux(self):
        func = lambda e: self(e) * e**2
        return quad(func, self.emin, self.emax)[0]

class LogParabola(FluxModel):

    def __call__(self, e):
        n0,alpha,beta,e_break=self.pars
        x = np.log(e_break/e)
        y = (alpha - beta*x)*x
        return n0*np.exp(y)

class PLSuperExpCutoff(FluxModel):

    def __call__(self,e):
        #print('WARNING: check e0!')
        n0,gamma,cutoff,b=self.pars
        return n0*(self.e0/e)**gamma*np.exp(-(e/cutoff)**b)