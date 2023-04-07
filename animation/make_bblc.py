"""
Manage a Bayesian-block Light Curve Repository
"""
import numpy as np
from zipfile import ZipFile
from pathlib import Path
import matplotlib.pyplot as plt


def load_vars(vmin=100):
    """Return 4FGL DF with variable sources, sorted by energy flux, and corespondit SkyCoord with positions
    """
    from utilities.catalogs import Fermi4FGL
    cat4 = Fermi4FGL()

    var_cut = (cat4.variability>vmin) & (cat4.class1!='PSR') & (cat4.class1!='HMB') 
    cv = cat4[var_cut]

    # size to use in the scatter plot
    s = np.sqrt(cv.eflux/4e-15)
    scut = s>30 
    df = cv[scut].sort_values('eflux', ascending=False)
    print(f'Selected {len(df)} variable sources for BBLC generation')
    return df


class BBLCR(ZipFile):
    """ Access the Bayesian-block 4FGL light curves in a zip file
    """
    def __init__(self, lczipfile='files/lcr.zip'):
        from pathlib import Path
        assert Path(lczipfile).is_file()
        super().__init__(lczipfile)
        
    def __len__(self):
        return len(self.filelist)
    
    @property
    def index(self):
        """  list of source names derived from filenames"""
        return np.array([ Path(f.filename).name.replace('_',' ')[:-4] for f in self.filelist])
    
    def __getitem__(self, idx):
        """Return a tuple (name, bb light curve)

        where name is the 4FGL name, and the light curve is a Series
        including t: MJD time; tw: width in days; and flux: relative flux
    
        """
        import pickle
        fname = self.filelist[idx]
        srcname = Path(fname.filename).name[:-4].replace('4FGL_','4FGL ')
        lc = pickle.load(self.open(fname))

        # now expand the BB light curve data frame entry to weekly
        cumwk = np.cumsum(lc.tw.values)/7
        ix  = np.arange(int(cumwk[-1]))
        return srcname, lc.flux[np.searchsorted(cumwk, ix)]
        


    @property
    def asdict(self):
        """ Return all weekly light curves as a dict"""
        return dict(x for x in self)
    

class LightCurves:
    
    def __init__(self, bblcr, cat):
        self.bb = bblcr
        self.cat = cat
        print(f'Loaded zip LC repository: {len(bblcr)} entries.')
        srclist = self.bb.index
        catlist = self.cat.index
        a, b = set(srclist), set(catlist)
        diff = a.difference(b)
        assert len(diff)==0 , f'Sources {diff} in BBLC not found in 4FGL-DR3'
        missing = b.difference(a)
        if len(missing)>0:
            print(f'*** Missing {len(missing)} souurces')

    def __len__(self): return len(self.bb)

    def __getitem__(self, idx):
        name, lc = self.bb[idx]
        eflux = (self.cat.loc[name].eflux/1e-12).astype(np.float32)
        return name, lc.values*eflux


    def plot_static(self):
        """Make an Aitoff plot of the 
        
        """

        from wtlike.skymaps import AitoffFigure
        from astropy.coordinates import SkyCoord
        
        s = np.sqrt(self.cat.eflux/4e-15)
        v = self.cat.variability
        logc = np.log10(v)
        skd  = SkyCoord(self.cat.ra, self.cat.dec, unit='deg', frame='fk5')
        
        
        fig = plt.figure(figsize=(15,6))
        fig.suptitle(f'{len(self.cat)} variable sources in 4FGL-DR3', fontsize=16, ha='right');
        plt.subplots_adjust(top=0.95)

        afig = AitoffFigure(fig) 

        afig.ax.set_facecolor('lavender')

        sc = afig.scatter(skd, c=logc, s=s, 
                        alpha=0.9, marker='o', cmap='Reds', edgecolor='grey')

        cb = plt.colorbar(sc, ax=afig.ax)    
        cb.set_label(r'$\mathrm{log_{10}(variability)}$') 

        return fig
    
    def plot_lc_grid(self):
        
        fig, axx = plt.subplots(ncols=5, nrows=10, figsize=(20,20), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0, hspace=0)
        for (name, lc), ax  in zip(self, axx.flatten()):
            ax.plot(np.arange(len(lc)), lc, color='maroon');
            ax.text(0.025,0.92, self.cat.loc[name].assoc1_name,  
                    transform=ax.transAxes, va='top')
        return fig



    

def run(N=0, repository='files/lcr'):

    """Create a set of light curves in the repository folder
    """
    import os   
    from wtlike import WtLike, Config

    df = load_vars()
    os.makedirs(repository, exist_ok=True)

    for name in df.index:
        fn = Path(repository)/(name.replace(' ', '_')+'.pkl')
        if fn.is_file(): continue
        print(f'Saving BBLC to {fn}')
        bb = WtLike(name, key=None, config=Config(verbose=0)).bb_view().fluxes
        with open(fn, 'wb') as out:
            bb.to_pickle( out)
        N-=1
        if N==0: break

    