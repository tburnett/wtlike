"""
"""
import sys, pickle
import numpy as np
import pandas as pd
from pathlib import Path
from wtlike import FermiInterval

def draw_map(mapper, time_range, name='', show_sun=False, sigma=0,
                fig=None, figsize=(9,4), label='', fignum=1, **kwargs):

    """Draw an aitoff map.

    - mapper : the HEALPix map, or a function of a (start,stop) time range that returns one
    - time_range : (start,stop) tuple
    - name ['']: name to give to the map which will be drawn in UL corner
    - label ['']: text to draw in LL corner for time interval, e.g., month number
    - sigma [0]: smoothing parameter (degrees)
    - show_sun [False] : flag to show the sun path
    - utc_flage [True] : show progress bar in UTC, else MJD
    - ** kwargs : pass to the ait_plot function, e.g. vmin, vmax log, colorbar, title
    """
    from wtlike.config import first_data, MJD
    from wtlike.skymaps import HPmap
    import matplotlib.pyplot as plt

    def sun_dir(mjd):
        """The sun direction in galactic coordinates
        """
        from astropy.coordinates import get_sun, SkyCoord
        from astropy.time import Time
        s =  get_sun(Time(mjd, format='mjd'))
        return SkyCoord(s.ra, s.dec, frame='fk5').galactic

    # set up the figure and an aitoff Axes
    fig = plt.figure(figsize=figsize, num=fignum) if fig is None else fig
    time_bar = kwargs.pop('show_time_bar', True)
    ax1 = fig.add_axes([0.10,0.2,0.9,0.9], projection='aitoff')
    kw = dict(log=False, tick_labels=False, vmin=None, vmax=None, title='', colorbar=True)
    kw.update(**kwargs)
    hpm=HPmap( mapper(time_range) if callable(mapper) else mapper, name, sigma=sigma)
    hpm.ait_plot(  ax=ax1,  **kw);

    if show_sun:
        sd = sun_dir(np.linspace(*time_range, num=100))
        l = sd.l.radian
        l[l>np.pi] -= 2*np.pi
        sdf = pd.DataFrame(dict(l=-l, b=sd.b.radian)).sort_values('l')
        ax1.plot(sdf.l, sdf.b, '-', lw=4, color='yellow')

    # set up a second Axes object to show the time range
    if time_range!=(0,0) and time_bar:
        ax2 = fig.add_axes([0.35, 0.1, 0.55, 0.11])
        yrs = list(range(2008,2024,2))

        ax2.axvspan(*time_range, color='orange')
        ax2.axvspan(MJD('2008'), first_data, color='lightgrey')
        yrs = [str(yr) for yr in range(2008,2025, 2)]
        ax2.set( xlim=(first_data, MJD('now')), ylim=(0,120), yticks=[], aspect=1,
                xticks=[MJD(yr) for yr in yrs],)#        xticklabels=yrs,)
        ax2.set_xticklabels(yrs, fontsize=10)
        if label: fig.text(0.33, 0.14, label,fontsize=10, ha='right')
    return fig



class Tranche:

    def __init__(self, size, total):
        self.start = np.append(np.arange(0, total, size), total)
    
    def __getitem__(self, k):
        return  list(range(self.start[k], self.start[k+1]))
    
    def __len__(self): 
        return len(self.start)-1

    def __repr___(self):
        return f'Tranche(0, {self.start[-1]}, {self.start[1]})'

def week_number(interval):
    
    from wtlike.config import first_data
    return int((interval[0]-first_data)/7)+1


class ExposureRatio:
    def __init__(self, filename):
        assert Path(filename).is_file()
        with open(filename, 'rb') as inp:
            self.edict = pickle.load(inp)
        self.average = self.edict[0]/(len(self.edict)-1)

    def __call__(self,interval):
        if interval==(0,0):
            return self.edict[0]
        week = week_number(interval)
        ratio = self.edict[week]/self.average
        ratio[np.isnan(ratio)]=0
        return ratio
    def __len__(self): return len(self.edict)-1

class FluxDifference:

    def __init__(self,filename):
        assert Path(filename).is_file()
        with open(filename, 'rb') as inp:
            self.fdict = pickle.load(inp)
        self.average = self.fdict[0]

    def __call__(self, interval):
        if interval==(0,0):
            return self.fdict[0]

        week = week_number(interval)
        diff = self.fdict[week]-self.average
        diff[np.isnan(diff)]=0
        return diff
    
    def __getitem__(self, week):
        diff = self.fdict[week]-self.average
        diff[np.isnan(diff)]=0
        return diff

    def __len__(self): return len(self.fdict)-1

class FluxPlot:
    def __init__(self, fdict, **kwargs):
        assert Path(fdict).is_file()
        self.mapper = FluxDifference(fdict)
        self.kw = dict(colorbar=False,  grid_color='0.25', vmin=1, vmax=10,)
        self.kw.update(kwargs)

    def __len__(self): return len(self.mapper)
        
    def __call__(self, week, fig=None, show_time_bar=True):        
        interval = FermiInterval(7)[week-1]  if week>0 else (0,0)
        label = f'week {week:3d}' if week>0 else '14+ yeears'         
        fig = draw_map(self.mapper, interval, label=label, fig=fig,
                       show_time_bar=show_time_bar,   **self.kw)           
        return fig

class ExposurePlot:
    def __init__(self,edict, **kwargs):
        assert Path(edict).is_file()
        self.mapper = ExposureRatio(edict)
        self.kw = dict(colorbar=False,  grid_color='0.75', vmin=0, vmax=2,
             cmap='coolwarm', show_sun=True)
        self.kw.update(kwargs)


    def __len__(self): return len(self.mapper)
        
    def __call__(self, week, fig=None, show_time_bar=True):        
        interval = FermiInterval(7)[week-1] if week>0 else (0,0) 
        label = f'week {week:3d}' if week>0 else '14+ yeears'        
        fig = draw_map(self.mapper, interval, label=label, fig=fig, 
                       show_time_bar=show_time_bar, **self.kw)
        return fig


class PlottingJob:

    def __init__(self, jobname, tranche_size=200, job_time='1:00:00'):
        self.jobname = jobname
        self.job_time = job_time
        kind=None
        for test in 'flux exposure'.split():
            if jobname.find(test)>=0: kind=test
        if kind is None: raise Exception(f'Bad jobname: {jobname}')
        p = Path('files') #weekly_expos
        map_dict_file = sorted(list(Path('files').glob(f'weekly_{kind}*.pkl')))[-1]
            
        self.outpath = p/jobname
        self.outpath.mkdir(exist_ok=True)
        (self.outpath/'logs').mkdir(exist_ok=True)
        (self.outpath/'plots').mkdir(exist_ok=True)
        self.plotter = ExposurePlot(map_dict_file) if kind=='exposure' else FluxPlot(map_dict_file)
        self.tranche = Tranche(tranche_size, len(self.plotter))
  
    def __len__(self):
        return len(self.tranche)

    def __call__(self, tid):

        for idx in self.tranche[tid-1]:
            self.plot_week(idx+1)

    def plot_week(self, week):
        
        outfile = self.outpath/f'plots/week_{week:03d}.png'
        if outfile.is_file(): return
        fig = self.plotter(week)
        fig.savefig(outfile, bbox_inches='tight')
        fig.clear()

    def submit(self):
        from simple_slurm import Slurm
        myfile = Path(sys.argv[0]).stem
        module = 'pylib.'+myfile
        output=self.outpath/'logs/%A_%a.log'
        # print(f'log output to {output} ')
        print(f'Running module {module}')
        slurm = Slurm(
            job_name=self.jobname,
            output=output,
            time=self.job_time,
            array = f'1-{len(self)}',
            ntasks=len(self),
            )
        jobid =  slurm.sbatch(f'python -m {module}')
        return jobid

        