# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/04_skymaps.ipynb (unless otherwise specified).

__all__ = ['valid', 'HPmap', 'AitoffFigure', 'ait_plot', 'plot_week', 'SquareWCS', 'sun_dir', 'make_ltmap']

# Cell
import os, sys
import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

import healpy
from healpy.rotator import Rotator

from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
#from utilities import healpix as hpx

from .config import *
from .data_man import DataView, get_week_map
from .sources import findsource

valid = Config().valid;

# Cell
class HPmap(object):
    """
    Manage HEALPix array
    """
    def __init__(self,
            hpmap:'HEALPix array',
            name='',
            cblabel=None,
            unit='',
            sigma=None,
            nest=False):
        """create from a HEALPix array
        """
        self.name = name
        self.cblabel = cblabel if cblabel is not None else unit
        self.unit = unit
        self.nside = healpy.get_nside(hpmap)
        # reorder as defaut RING if nest is set
        self.map = hpmap if not nest else healpy.reorder(hpmap, n2r=True)

        if sigma is not None:
            self.smooth(sigma)

    def __str__(self):
        return f'<{self.__class__.__name__}>, name "{self.name}" nside {self.nside} colorbar label "{self.cblabel}"'
    def __repr__(self): return str(self)

    def __call__(self, sc:'SkyCoord') -> 'value[s]':
        """
        Return value of corresponding pixel
        """
        skyindex = healpy.ang2pix(self.nside, sc.l.deg, sc.b.deg, lonlat=True)
        return self.map[skyindex]

    def get_all_neighbors(self, pix):
        return healpy.get_all_neighbours(self.nside, pix)

    def smooth(self, sigma):
            self.map = healpy.smoothing(self.map, np.radians(sigma))

    @classmethod
    def from_FITS(cls, filename, *pars, **kwargs):
        with fits.open(filename) as hdus:
            header, data = hdus[1].header, hdus[1].data
        kw = dict(unit=header.get('TUNIT1', ''), name=header['TTYPE1'])
        kw.update(**kwargs)
        return cls(data.field(0), *pars, **kw)

    def ait_plot(self,  **kwargs):
        """
        Return a AitoffFigure object, which can be used to add source positions and text
        """
        kw= dict(label=self.name, cblabel=self.cblabel,)
        kw.update(**kwargs)
        return AitoffFigure(ait_plot(self, **kw))

    def to_FITS(self,  filename=''):
        """return a HDUlist object with one skymap column

        - filename [''] write to the file if set
        """

        column = fits.Column(name=self.name, format='E', array=self.map, unit=self.unit)

        nside = self.nside
        cards = [fits.Card(*pars) for pars in [
            ('FIRSTPIX', 0,             'First pixel (0 based)'),
            ('LASTPIX',  12*nside**2, 'Last pixel (0 based)'),
            ('MAPTYPE',  'Fullsky' , ''  ),
            ('PIXTYPE',  'HEALPIX',      'Pixel algorithm',),
            ('ORDERING', 'RING',         'Ordering scheme'),
            ('NSIDE' ,   nside,    'Resolution Parameter'),
            ('ORDER',    int(np.log2(nside)),   'redundant'),
            ('INDXSCHM', 'IMPLICIT' ,''),
            ('OBJECT' ,  'FULLSKY', ''),
            ('COORDSYS', 'GAL', ''),
        ]]
        table = fits.BinTableHDU.from_columns([column],header=fits.Header(cards))
        table.name = 'SKYMAP'

        hdus = fits.HDUList([
                    fits.PrimaryHDU(header=None),
                    table,
                ])
        if filename:
            hdus.writeto(filename, overwrite=True)
            print(f'Wrote FITS file {filename}')
        return hdus

# Cell
class AitoffFigure():

    """ Implement plot and text converting from (l,b) in degrees, or a SkyCoord.

    """
    def __init__(self, fig):
        self.fig = fig
        self.ax = fig.axes[0]
        assert self.ax.__class__.__name__ == 'AitoffAxesSubplot', 'expect figure to have aitoff axes instance'

    def _trans(self, *args):
        if len(args)==2 and isinstance(args[0], SkyCoord):
            sc = args[0].galactic
            l, y, a = sc.l.rad, sc.b.rad, args[1]
            x = -(l if l<=np.pi else l- 2*np.pi)
        elif len(args)==3:
            l, b, a = args
            x = -np.radians( (l if l<=180 else l-360 ))
            y = np.radians(b)
        else:
            raise Exception('Expect positional parameters l,b,a, or skycoord,a')
        return x,y, a

    def plot(self, *args, **kwargs):
        self.ax.plot(*self._trans(*args), **kwargs)

    def text(self, *args, **kwargs):
        self.ax.text(*self._trans(*args), **kwargs)

def ait_plot(mappable,
        pars=[],
        label='',
        title='',
        fig=None, ax=None, fignum=1, figsize=(15,8),
        pixelsize:'pixel size in deg'=1,
        projection='aitoff',
        cmap='jet',
        vmin=None, vmax=None,
        log=False,
        colorbar=True,
        cblabel='',
        cb_kw={},
        axes_pos=111,
        axes_kw={},
        tick_labels=True,
        alpha:'apply to pcolormesh'=None,
        ):
    """
    """
    #
    # code inspired by https://stackoverflow.com/questions/46063033/matplotlib-extent-with-mollweide-projection

    # make a mesh grid
    nx, ny = int(360/pixelsize), int(180/pixelsize)
    lon = np.linspace(180,-180, nx) # note reversed
    lat = np.linspace(-90., 90, ny)
    Lon,Lat = np.meshgrid(lon,lat)

    #  an arrary of values corresponding to the grid
    # dirs = SkyDir.from_galactic(Lon, Lat)
    dirs = SkyCoord(Lon, Lat, unit='deg', frame='galactic')
    arr = mappable(dirs, *np.atleast_1d(pars))

    if ax:
        fig = ax.figure
        assert ax.__class__.__name__.startswith('AitoffAxes'), 'Require that be a AitoffAxes object'
    else:
        fig = plt.figure(figsize=figsize, num=fignum) if fig is None else fig
        # this needs to be more flexible
        ax = fig.add_subplot(axes_pos, projection=projection, **axes_kw)

    # reverse longitude sign here for display
    if log:
        norm = colors.LogNorm(vmin=vmin,vmax=vmax)
        vmin=vmax=None
    else:
        norm = None
    im = ax.pcolormesh(-np.radians(Lon), np.radians(Lat), arr,  shading='nearest',
        norm=norm, cmap=cmap,  vmin=vmin, vmax=vmax, alpha=alpha)

    if tick_labels:
        ff = lambda d: d if d>=0 else d+360
        ax.set_xticklabels([f'${ff(d):d}^\degree$' for d in np.linspace(150,-150, 11).astype(int)])
    else:
        ax.set(xticklabels=[], yticklabels=[], visible=True)

    if colorbar:
        ticklabels = cb_kw.pop('ticklabels', None)
        cb_kw.update(label=cblabel,)
        cb = plt.colorbar(im, ax=ax, **cb_kw)
        if ticklabels is not None:
            cb.ax.set_yticklabels(ticklabels)
    ax.grid(color='grey')
    if label:
        ax.text( 0., 0.97, label, transform=ax.transAxes)
    if title:
        plt.suptitle(title, fontsize=12)

    return fig

# Cell
def plot_week(week=None, mjd=None, nside=32, **kwargs):
    """
    Make an AIT plot of the given week's photon data

    Combine all energies for now

    - week -- the week number from 9
    - mjd -- [None] If set, derive the week from it
    - nside [32] -- HEALPix nside to project data before plotting.
    - kwargs -- args for healpix.ait_plot

    """
    from utilities import healpix as hpx
    from .config import mission_week
    if mjd is not None: week=mission_week(mjd)

    pmap = get_week_map(week, nside=nside)

    kw = dict(log=True, tick_labels=False, vmin=5, vmax=1e3,
             cblabel=f'counts per nside={nside} pixel')

    t = hpx.HPmap(pmap, f'week_{week:d}', #\n{utc}',
                  **kwargs)#f'week_{week:d}\n{utc}', nest=True)
    t.ait_plot(**kw)
    return plt.gcf()

# Cell

class SquareWCS(WCS):
    """
    Create and use a WCS object

    - center : a SkyCoord that will be the center
    - size   : width and height of the display
    - pixsize [0.1] : pixel size
    - frame [None] : The frame is taken from the center SkyCoord, unless specified here --  only accept "galactic" or "fk5"
    - proj ["ZEA"] : projection to use
    """

    def __init__(self, center, size, pixsize=0.1, frame=None, proj='ZEA'):
        """

        """
        assert isinstance(center, SkyCoord), 'Expect SkyCoord'

        frame = frame or center.frame.name
        if frame=='galactic':
            lon, lat = center.galactic.l.deg, center.galactic.b.deg
            lon_name,lat_name = 'GLON','GLAT'
            self.axis_labels='$l$', '$b$'
        elif frame=='fk5':
            lon,lat = center.fk5.ra.deg, center.fk5.dec.deg
            lon_name, lat_name = 'RA--', 'DEC-'
            self.axis_labels = 'RA', 'Dec'
        else:
            raise Exception(f'Expect frame to be "galactic" or "fk5", not {frame}')

        nx=ny=naxis = int(size/pixsize) | 1 # make odd so central pixel has source in middle
        self.center = center
        self.frame=frame
        self.galactic = frame=='galactic'
        super().__init__(
                         dict(
            NAXIS1=nx, CTYPE1=f'{lon_name}-{proj}', CUNIT1='deg', CRPIX1=nx//2+1, CRVAL1=lon, CDELT1=-pixsize,
            NAXIS2=ny, CTYPE2=f'{lat_name}-{proj}', CUNIT2='deg', CRPIX2=ny//2+1, CRVAL2=lat, CDELT2=pixsize, )
              )

    def _make_grid(self):
        # get coordinates of every pixel`
        nx, ny = self.array_shape
        pixlists = list(range(1,nx+1)),list(range(1,ny+1))
        cgrid = self.pixel_to_world(*np.meshgrid(*pixlists) )
        if not self.galactic:
            cgrid = cgrid.galactic
        lon, lat = (cgrid.l.deg, cgrid.b.deg)
        return lon, lat

    def plot(self, hmap, log=False, cmap='jet', colorbar=True,
             cblabel='', vmin=None, vmax=None, cb_kw={},
             annotator=None, title=None):
        """
        - hmap -- a HEALPix map

        """

        import healpy as hp
        from matplotlib import colors

        wcs = self
        grid = self._make_grid();
        nside = hp.get_nside(hmap)

        # lon, lat = grid.l.deg, grid.b.deg
        ipix = hp.ang2pix(nside, *grid, lonlat=True)

        fig = plt.figure(figsize=(6,6))
        fig.add_subplot(111, projection=self)
        ax = fig.axes[0]

        if log:
            norm = colors.LogNorm(vmin=vmin,vmax=vmax)
            vmin=vmax=None
        else:
            norm = None
        ipix = hp.ang2pix(nside, *grid, lonlat=True)
        im = ax.imshow(hmap[ipix], cmap=cmap, origin='lower', norm=norm, vmin=vmin);

        nx, ny = wcs.array_shape
        ax.set(xlabel=self.axis_labels[0], xlim=(-0.5, nx-0.5),
               ylabel=self.axis_labels[1], ylim=(-0.5, ny-0.5),
              title= title)
        ax.grid();
        if colorbar:
            ticklabels = cb_kw.pop('ticklabels', None)
            cb_kw.update(label=cblabel,)
            cb = plt.colorbar(im, ax=ax, **cb_kw)
            if ticklabels is not None:
                cb.ax.set_yticklabels(ticklabels)
        if annotator is not None:
            annotator(ax, self.frame)
        return fig

# Cell
def sun_dir(mjd):
    """The sun direction in galactic coordinates
    """
    from astropy.coordinates import get_sun, SkyCoord
    from astropy.time import Time
    s =  get_sun(Time(mjd, format='mjd'))
    return SkyCoord(s.ra, s.dec, frame='fk5').galactic

def make_ltmap(time_range, sigma=1, show_sun=True,
               figsize=(9,4), utc_flag=True,fignum=1,
               **kwargs):
    """Make a livetime map.
    """

    from utilities.healpix import HPmap

    from .data_man import DataView
    from .config import first_data, MJD
    # make the healpux livetime map for the time range
    dv  =DataView(time_range)
    ltmap = dv.livetime_map(nside=128, sigma=sigma)

    # set up figure, make the plot and maybe colorbar
    fig = plt.figure(figsize=figsize, num=fignum) #(12,5))
    ax1 = fig.add_axes([0.15,0.2,0.95,0.95], projection='aitoff')
    kw = dict(log=False, tick_labels=False, vmin=None, vmax=None, title='', colorbar=True)
    kw.update(**kwargs)
    HPmap(ltmap, '').ait_plot( ax=ax1, **kw)

    if show_sun:
        sd = sun_dir(np.linspace(*time_range, num=30))
        l = sd.l.radian
        l[l>np.pi] -= 2*np.pi
        ax1.plot(-l, sd.b.radian, '.', color='orange')

    # set up ases object to show time
    ax2 = fig.add_axes([0.35, 0.1, 0.55, 0.11])
    kw2 = dict( xlim=(first_data, MJD('now')), ylim=(0,100),yticks=[], aspect=1 )
    ax2.axvspan(*dv.time_range, color='orange')
    if utc_flag:
        # even years if True else interpret as int
        cnt = 2 if type(utc_flag)==bool else utc_flag
        yrs = [str(yr) for yr in range(2008,2024, cnt)] #get this from kwarg maybe
        yrkw = dict( xticks=[MJD(yr) for yr in yrs], xticklabels=yrs,)#  xlabel='UTC',)
        kw2.update(**yrkw)
    ax2.set(**kw2)
    return fig
