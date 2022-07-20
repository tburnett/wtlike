from wtlike import *
from astropy.coordinates import SkyCoord
from wtlike.skymaps import *
from cache_decorator import Cache
from utilities.ipynb_docgen import *
plt.rc('font', size=16)

@Cache()
def get_count_map(nside=256, bmin=8):
    from wtlike.data_man import DataView
    return DataView().count_map(nside=nside, bmin=bmin)
# set cache if not already
get_count_map()

@ipynb_doc
def hpmap_demo():
    """
    ## HEALPix projections
    
    These are both invoked by member functions of the class `HPmap`, which manages a HEALPix array.
    ### All-sky: AIT
    {fig1}
    
    ### A region: ZEA
    {fig2}
    """
    cntmap = get_count_map() 
    hpmap = HPmap(cntmap,  unit=f'counts per pixel') 

    aitkw = dict(cmap='jet', pctlim=(50,99.75), #vlim = (10,1000), 
                 log=True, tick_labels=False, pixelsize=0.1, figsize=(20,10), 
            colorbar=True, cb_kw=dict(shrink=0.7),unit='counts/pixel')
    fig1 = hpmap.ait_plot(**aitkw);
    
    fig2 = hpmap.zea_plot(SkyCoord(180,0,unit='deg',frame='galactic'), log=True, cb_kw=dict(shrink=0.7), size=45,
                          title='Anti-center', figsize=8, vmin=15, vmax=2e3)

    return locals()