"""useful plot utilities
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats 


   

def profile_plot(x, y, bins, ax=None, grid=0.5, **kwargs):
    """ Make a "profile" plot 
      Show <y> with rms errors for each bin in x
      
      - **kwargs -- apply ms, ls, lw, color to errorbar; rest apply to the Axes object
    """
    # useful kwargs for errorbar
    ms=   kwargs.pop('ms', 'o')    # marker symbol
    ls =  kwargs.pop('ls', 'none') # line style
    lw =  kwargs.pop('lw', 2)      # line wiedth
    color=kwargs.pop('color', None)# color
    
    # get <y> and <y**2> for each bin in  x; calculate rms
    (my, myy), bins,_ = stats.binned_statistic(x, [y, y*y], bins=bins)
    std = np.sqrt(myy-my**2)
    # my, std, bins = profile_stats(x,y,bins)
    good = ~pd.isna(my) * ~pd.isna(std)
    x, xerr = ((bins[:-1]+bins[1:])/2)[good], (bins[1]-bins[0])/2
    y, yerr  = my[good], std[good]
    
    # make the errorbar plot
    fig, ax = plt.subplots() if ax is None else (ax.figure, ax)
    ax.errorbar(x, y, yerr, xerr,  ls=ls, marker=ms, color=color, lw=lw )
    ax.set(**kwargs)
    if grid: ax.grid(alpha=grid)
    
    