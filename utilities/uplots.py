"""useful plot utilities
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats 

def make_profile(x,y,bins):
    from scipy import stats 
    # get <y> and <y**2> for each bin in  x; calculate rms
    (my, myy), bins,_ = stats.binned_statistic(x, [y, y*y], bins=bins)
    std = np.sqrt(myy-my**2)
    # my, std, bins = profile_stats(x,y,bins)
    good = ~pd.isna(my) * ~pd.isna(std)
    x, xerr = ((bins[:-1]+bins[1:])/2)[good], (bins[1]-bins[0])/2
    y, yerr  = my[good], std[good]
    return pd.DataFrame(dict(x=x, y=y, xerr=xerr, yerr=yerr))
   

def profile_plot(x, y, bins, ax=None, grid=0.5, **kwargs):
    """ Make a "profile" plot 
      Show <y> with rms errors for each bin in x
      
      - **kwargs -- apply ms, ls, lw, color to errorbar; rest apply to the Axes object
    """
    # process the data
    pf = make_profile(x,y,bins)
    
    # useful kwargs for errorbar
    ms=   kwargs.pop('ms', 'o')    # marker symbol
    ls =  kwargs.pop('ls', 'none') # line style
    lw =  kwargs.pop('lw', 2)      # line wiedth
    color=kwargs.pop('color', None)# color
    
    # make the errorbar plot
    fig, ax = plt.subplots() if ax is None else (ax.figure, ax)
    ax.errorbar(pf.x, pf.y, pf.yerr, pf.xerr,  ls=ls, marker=ms, color=color, lw=lw )
    ax.set(**kwargs)
    if grid: ax.grid(alpha=grid)
  
  
class MultiHist:
    
    def __init__(self, 
                 ax_info,
                 histkw = {},
                 ncols=4,
                figsize=(12,3),
                ):
    
        hkw_default =  dict(histtype='step', hatch='////', lw=2)
        hkw = hkw_default; hkw.update(histkw)
        n = len(ax_info)
        ncols = min(ncols,n)
        nrows = int((n + ncols-1)//ncols)

        self.fig, axx = plt.subplots(nrows=nrows,
                                          ncols=ncols, figsize=figsize)
        plt.subplots_adjust(wspace=0.3, hspace=0.35)
        self.axx = axx.flatten()
        self.names=[]
        self.hkw = []
        for ax, (name,bins, kw) in zip(self.axx, ax_info):
            self.names.append(name)
            t = hkw.copy(); t.update(bins=bins)
            self.hkw.append(t)
            ax.set(**kw)
            ax.grid(alpha=0.5)
        for ax in self.axx[n:]: ax.set_visible(False)
        
    def fill(self, df):
        for ax, name,  kw in zip(self.axx, self.names, self.hkw):
            ax.hist(df[name], **kw)
