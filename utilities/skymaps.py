"""

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.io import fits
# https://astropy-healpix.readthedocs.io/en/latest/healpy_compat.html

# import healpy
from astropy_healpix import healpy

def get_nside(v):
    return healpy.npix_to_nside(len(v))

#| export
def _transform_plot_args(*args):
    """ 
    Helper for map displays position specification
    Expect args to be one of:
    * a SkyCoord object (which is perhaps a list of positions)
    * lists of l, b in degrees
    * pandas DataFrame or Series object with `glon` and `glat` columns
    * The name of a source known to SkyCoord  

    Returns a SkyCoord object using the first or first two args, and the remaining args
    Subsequently expect to convert to radians for AIT or pixels for ZEA. 
    """
        
    if len(args)==0: raise ValueError('No args')
    arg, rest = args[0], args[1:]
    if isinstance(arg, pd.DataFrame) or isinstance(arg, pd.Series):
        df = arg
        if 'glon' in df: 
            sc  = SkyCoord(df.glon, df.glat, unit='deg', frame='galactic')
        elif 'GLON' in df:
            sc  = SkyCoord(df.GLON, df.GLAT, unit='deg', frame='galactic')
        elif 'ra' in df:
            sc = SkyCoord(df.ra, df.dec, unit='deg', frame='fk5').galactic
        else:
            raise ValueError('DataFrame must have GLON,GLAT or glon,glat or ra,dec')
    elif isinstance(arg, SkyCoord):
        sc = arg.galactic
    elif type(arg)==str:
        sc = SkyCoord.from_name(arg).galactic
    elif len(args)>1:
        l,b = args[:2]
        rest = args[2:]
        sc = SkyCoord(l,b, unit='deg', frame='galactic')
    else:
        raise ValueError('Expect DataFrame with glat,glon or SkyCoord or name or l,b')
    return sc, rest

def _to_radians(*args):
 
    sc, rest = _transform_plot_args(*args)
    # convert to radians 
    l, b = sc.l.deg, sc.b.deg
    x  = -np.radians(np.atleast_1d(l))
    x[x<-np.pi] += 2*np.pi # equivalent to mod(l+pi,2pi)-pi I think
    y = np.radians(np.atleast_1d(b))
    return [x,y] + list(rest) 

def _to_pixels(*args, wcs):
    sc, rest = _transform_plot_args(*args)
    return wcs.world_to_pixel(sc) + tuple(rest) 

class SkyPlotMixin():
    """
    A mix-in for sky plotting classes, implementing common features.
    """

    def plot(self, *args, **kwargs):
        self.ax.plot(*self._transform_plot_args(*args), **kwargs)
        return self

    def text(self, *args, **kwargs):
        pars = self._transform_plot_args(*args)
        print(pars)
        if len(np.atleast_1d(pars[0]))==1:
            self.ax.text(*pars, **kwargs)
        elif len(pars)==3:
            for x,y,text in zip(*pars):
                self.ax.text(x,y,text, **kwargs)
        else:
            raise Exception(f'text parameters not understood')
        return self
    
    def axes_text(self, *args, **kwargs):
        """
        Add text to the plot with coordinates relative to the axes.

        :param args: Positional arguments for the text (x, y, text).
        :param kwargs: Additional keyword arguments for styling the text.
        """
        self.ax.text(*args, transform=self.ax.transAxes, **kwargs)
        return self
    
    def scatter(self, *args, **kwargs):
        """
        """
        x,y = self._transform_plot_args(*args)
        self.mappable = self.ax.scatter(x ,y, **kwargs)
        return self
    
    def colorbar(self, **kwargs):
        if hasattr(self, 'mappable'):
            self.cbar=self.ax.figure.colorbar(self.mappable, **kwargs)
            delattr(self, 'mappable') # prevent a second call to mappable
        else:
            raise Exception('No mappable to make colorbar from')
        return self
    
    def apply(self, func, *pars, **kwargs):
        """Pass self, which behaves like an Axis, to a user-supplied function and return self."""
        func(self, *pars, **kwargs)
        return self   
    
    def legend(self, *pars, **kwargs):
        self.ax.legend(*pars, **kwargs)
        return self
    
    def title(self, *pars, **kwargs):
        self.ax.set_title(*pars, **kwargs)
        return self
    def show(self):   
        plt.show()
        return None
    
    
class AITfigure(SkyPlotMixin):

    def __init__(self, fig=None, *, figsize=(10,5), grid_color='grey',
                 show_ticklabels=False, 
                 **kwargs):
        """
        fig -- [None] 
        figsize
        grid_color -- Suppress grid if None

        """
        if isinstance(figsize, ( int, float)): figsize=(figsize, figsize/2)
        self.figure = fig or plt.figure(figsize=figsize)
        if len(self.figure.axes)==0:
            ax=self.figure.add_subplot(111, projection='aitoff')
            if not show_ticklabels: 
                ax.set(xticklabels=[], yticklabels=[], visible=True)
            else:
                # reverse order and add 360 to negative ones
                # expect each label ends with degree symbol
                xtl = ax.get_xticklabels()[::-1]
                for txt in xtl:
                    t = txt.get_text()    
                    x = int(t[:-1])
                    if x>=0: continue
                    txt.set_text(str(x+360)+t[-1])
                ax.set_xticklabels(xtl)

            ax.grid(color=grid_color)
        self.ax = self.figure.axes[0]
        self._transform_plot_args = _to_radians

        assert self.ax.__class__.__name__.startswith('Aitoff'), 'expect figure to have aitoff Axes instance'
        if grid_color is not None:
            self.ax.grid(color=grid_color)
        self.ax.set(**kwargs)
        
    def __getattr__(self, name):
        # pass everything else to self.ax
        return self.ax.__getattribute__( name)
        
    def imshow(self, X, **kwargs):
        """ X: either a HEALPix array or a 2-D grid appropriate for Axes.imshow
        """
        try:
            return self.healpix_fill(X, **kwargs)
        except TypeError:
            # perhaps there is a valid way to display an array
            return self.ax.imshow(X, **kwargs)
           
    def healpix_fill(self, hparray, 
            pixelsize:float=1,
            colorbar:bool=False,
            cb_kw:dict={},
            unit:str='',
            **kwargs):
        """Fill with the values from the HEALPix array hparray
        """

        ax = self.ax

        nside = get_nside(hparray)
        
        # code inspired by https://stackoverflow.com/questions/46063033/matplotlib-extent-with-mollweide-projection
        # make a mesh grid for lon,lat in degrees
        Lon,Lat = np.meshgrid(np.linspace(180,-180, int(360/pixelsize)), # note reversed
                            np.linspace(-90., 90, int(180/pixelsize)))

        # get the pixels and look up values
        pixels = healpy.ang2pix(nside, Lon, Lat, lonlat=True)
        values = hparray[pixels]

        # plot them (reverse Lon here too)
        im =ax.pcolormesh(-np.radians(Lon), np.radians(Lat), values,  shading='nearest', **kwargs)

        if colorbar:
            ticklabels = cb_kw.pop('ticklabels', None)
            cb_kw.update(label=unit,)
            cb = plt.colorbar(im, ax=ax, **cb_kw)
            if ticklabels is not None:
                cb.ax.set_yticklabels(ticklabels)
        
        # A kluge that will show the grid on top of image
        ax.set_axisbelow(False) 
        self.mappable= im # for processing by colorbar
        return self

def ait_plot(mappable,
        pars=[],
        label='',
        title='',
        fig=None, ax=None, fignum=1, figsize=(20,8),
        pixelsize:'pixel size in deg'=1.0,
        projection='aitoff',
        cmap='inferno',
        vmin=None, vmax=None,  vlim=None, pctlim=None,
        log=False,
        colorbar:bool=True,

        unit='',
        grid_color='grey',
        cb_kw={},
        axes_pos=111,
        axes_kw={},
        tick_labels=False,
        alpha:'apply to pcolormesh'=None,
        ):

    """
    """
    from matplotlib import colors
    #
    # code inspired by https://stackoverflow.com/questions/46063033/matplotlib-extent-with-mollweide-projection
    # make a mesh grid
    nx, ny = int(360/pixelsize), int(180/pixelsize)
    lon = np.linspace(180,-180, nx) # note reversed
    lat = np.linspace(-90., 90, ny)
    Lon,Lat = np.meshgrid(lon,lat)

    #  an arrary of values corresponding to the grid
    # dirs = SkyDir.from_galactic(Lon, Lat)
    dirs = SkyCoord(Lon, Lat, unit='deg')
    arr = mappable(dirs, *np.atleast_1d(pars))

    if ax is not None:
        fig = ax.figure
        assert ax.__class__.__name__.startswith('AitoffAxes'), 'Require that be a AitoffAxes object'
    else:
        if isinstance(figsize, ( int, float)): figsize=(figsize, figsize/2)
        fig = plt.figure(figsize=figsize, num=fignum) if fig is None else fig
        # this needs to be more flexible
        ax = fig.add_subplot(axes_pos, projection=projection, **axes_kw)

    # reverse longitude sign here for display
    if pctlim is not None:
        vlim=np.percentile(mappable.map, np.array(pctlim)) #[40, 99.98])).round(),

    if vlim is not None:
        vmin, vmax = vlim

    if log:
        norm = colors.LogNorm(vmin=vmin,vmax=vmax)
        vmin=vmax=None
    else:
        norm = None

    
    im = ax.pcolormesh(-np.radians(Lon), np.radians(Lat), arr,  shading='nearest',
        norm=norm, cmap=cmap,  vmin=vmin, vmax=vmax, alpha=alpha)

    if tick_labels:
        ff = lambda d: d if d>=0 else d+360
        ax.set_xticklabels([r'${{ff(d):d}}^\degree$' for d in np.linspace(150,-150, 11).astype(int)])
    else:
        ax.set(xticklabels=[], yticklabels=[], visible=True)

    if colorbar:
        ticklabels = cb_kw.pop('ticklabels', None)
        cb_kw.update(label=unit,)
        cb = plt.colorbar(im, ax=ax, **cb_kw)
        if ticklabels is not None:
            cb.ax.set_yticklabels(ticklabels)
    if grid_color: 
        ax.grid(color=grid_color)
        ax.set_axisbelow(False) # klugy thing to show grid on top of image
    if label:
        ax.text( 0., 0.97, label, transform=ax.transAxes)
    if title:
        plt.suptitle(title, fontsize=12)
    return AitoffFigure(fig)


class ZEAfigure(WCS, SkyPlotMixin):
    """
    Create a WCS image rectangle

    - center : a SkyCoord that will be the center, or (l,b) pair, or a source name
    - size   : (width,height) of the display (deg) or square if singlet
    - pixelsize [0.1] : pixel size (deg)
    - frame [None] : The frame is taken from the center SkyCoord, 
            unless specified here --  only accept "galactic", "fk5" or "icrs"
    - proj ["ZEA"] : projection to use
    - axes_visible [True] : whether to show the axes ticks and labels. If False,
            the axes ticks and labels are hidden and the grid is turned off.
    
    To get the WCS properties from the generated Axes object (actually WCSAxesSubplot):
         ax.wcs.wcs.crval for (l,b)
    """
    def __init__(self, center, size, *, figsize=(5,5), fig=None, axpos=111,
                 pixelsize:float=0.1, axes_visible=True,
                 frame=None, proj='ZEA', unit='', major_formatter='d.d', **kwargs):
        if isinstance(figsize, (int,float)): figsize=(figsize, figsize)
        frame = 'galactic' if frame is None else frame.lower()
        if isinstance(center, SkyCoord):
            pass
        elif type(center)==tuple:
            center = SkyCoord(*center,unit='deg', frame=frame)
        elif type(center)==str:
            center = SkyCoord.from_name(center).transform_to(frame)
        else: raise Exception( 'Expect center to be: SkyCoord, name, or (l,b) tuple')
        # size is single float or tuple(floats) for width, height in pixelsize units
        size = np.atleast_1d(size).astype(float)
        if len(size)==1: size=np.full(2, size[0])
        assert len(size)==2, 'Expect size to be single float, or (float, float) tuple'
            
        frame = frame or center.frame.name
        if frame=='galactic':
            lon, lat = center.galactic.l.deg, center.galactic.b.deg
            lon_name,lat_name = 'GLON','GLAT'
            xlabel, ylabel='$l$', '$b$'
        elif frame=='fk5' or frame=='icrs':
            lon,lat = center.fk5.ra.deg, center.fk5.dec.deg
            lon_name, lat_name = 'RA--', 'DEC-'
            xlabel, ylabel = 'RA', 'Dec'
        else:
            raise ValueError(f'Expect frame to be "galactic" "fk5" "icrs",  not {frame}')
        self.galactic = frame=='galactic'
        self.frame = frame
        self.center = center
        self.unit=unit
        ny, nx = (size/pixelsize).astype(int)
        super().__init__( dict(
            NAXIS1=nx, CTYPE1=f'{lon_name}-{proj}', CUNIT1='deg', CRPIX1=ny//2+1, CRVAL1=lon, CDELT1=-pixelsize,
            NAXIS2=ny, CTYPE2=f'{lat_name}-{proj}', CUNIT2='deg', CRPIX2=nx//2+1, CRVAL2=lat, CDELT2=pixelsize, )
             )
        if isinstance(figsize, (int,float)): figsize=(figsize, figsize)
        if fig is None: fig =  plt.figure(figsize=figsize)
        self.figure = fig
        self.ax = ax =  fig.add_subplot(axpos, projection=self)
        nx, ny = self.array_shape
        ax.set(xlim=(-0.5, nx-0.5), xlabel=xlabel, 
                  ylim=(-0.5, ny-0.5), ylabel=ylabel)# title=title)
        lon_axis = ax.coords[0]
        lat_axis = ax.coords[1]
        lon_axis.set_format_unit('deg')
        lat_axis.set_format_unit('deg')
        lon_axis.set_major_formatter(major_formatter)
        lat_axis.set_major_formatter(major_formatter)
        kw = dict(xlim=(-0.5, nx-0.5), xlabel=xlabel, 
                  ylim=(-0.5, ny-0.5), ylabel=ylabel)
        kw.update(kwargs)
        ax.set(**kwargs)
        if not axes_visible:
            for axis in ax.coords:
                axis.set_ticklabel_visible(False)
                axis.set_ticks_visible(False)
            ax.grid(False)

        self._transform_plot_args = _transform_plot_args

    def imshow(self, dmap, *, cmap='jet', 
               log=False, vmax=None, vmin=None, norm=None,
               colorbar=False, tick_labels=None, cb_kw={},
               grid=True, grid_kw=None, **kwargs
               ):
        
        """ par --a HEALPix array, ring indexing, or a function that takes a SkyCoord and returns values, 
            or a 2-D array to display with imshow
        """
        from matplotlib import colors
        ax = self.ax

        if log:
            norm = colors.LogNorm(vmin=vmin,vmax=vmax)
            vmin=vmax=None
        else:
            norm = None
       
        if callable(dmap):
            # get the coordinates of the pixels and evaluate the function on them
            ny, nx = self.array_shape
            pixel_y, pixel_x = np.indices((ny, nx))
            pixcoord = self.pixel_to_world(pixel_x, pixel_y)
            self.mappable = self.ax.imshow(dmap(pixcoord),
                                           cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, **kwargs)

        elif isinstance(dmap, np.ndarray) and dmap.ndim==2:
            self.mappable = self.ax.imshow(dmap, cmap=cmap, norm=norm, vmin=vmin, vmax=vmax, **kwargs)
         # otherwise expect a HEALPix array, and get the pixel coordinates to look up the
        
        elif isinstance(dmap, np.ndarray) and dmap.ndim==1:
            def _process_healpix_array(): 

                # make a meshgrid of the coordinates of the pixels to display
                nside = get_nside(dmap)
                nx, ny = self.array_shape
                pixlists = list(range(1,nx+1)),list(range(1,ny+1))
                cgrid = self.pixel_to_world(*np.meshgrid(*pixlists) )
                
                # get the corresponding pixel indices
                if not self.galactic:   cgrid = cgrid.galactic
                lon, lat = (cgrid.l.deg, cgrid.b.deg)
                ipix = healpy.ang2pix(nside, lon,lat, lonlat=True)
                # return the values
                return dmap[ipix]    
            
            self.mappable = self.ax.imshow(_process_healpix_array(),
                        cmap=cmap,  norm=norm,  vmin=vmin, vmax=vmax, **kwargs
                        )
        else:
            raise ValueError('Expect dmap to be a 2-D array, a HEALPix array, or a function that takes a SkyCoord and returns values')  
        
        if tick_labels:
            ff = lambda d: d if d>=0 else d+360
            ax.set_xticklabels([r'${{ff(d):d}}^\degree$' for d in np.linspace(150,-150, 11).astype(int)])
        else:
            ax.set(xticklabels=[], yticklabels=[], visible=True)

        if grid:
            grid_kw = dict(color='white', alpha=0.5) if grid_kw is None else grid_kw
            ax.grid(**grid_kw)

        if colorbar:
            ticklabels = cb_kw.pop('ticklabels', None)
            cb_kw.update(label=self.unit,)
            cb = plt.colorbar(self.mappable, ax=self.ax, **cb_kw)
            if ticklabels is not None:
                cb.ax.set_yticklabels(ticklabels)
        return self

    def scatter(self, *args, **kwargs):
        self.mappable = self.ax.scatter_coord(_transform_plot_args(*args)[0], **kwargs)
        return self
    
    def plot(self, *args, **kwargs):
        sc, rest = _transform_plot_args(*args)
        if rest:  self.ax.plot_coord(sc, rest, **kwargs)
        else:     self.ax.plot_coord(sc, **kwargs)
        return self


class SquareWCS(WCS):
    """
    Create and use a WCS object

    - center : a SkyCoord that will be the center
    - size   : width and height of the display (deg)
    - pixelsize [0.1] : pixel size
    - frame [None] : The frame is taken from the center SkyCoord, 
            unless specified here --  only accept "galactic" or "fk5"
    - proj ["ZEA"] : projection to use
    
    To get the WCS properties from the generated Axes object (actually WCSAxesSubplot):
         ax.wcs.wcs.crval for (l,b)
    """

    def __init__(self, center, size, pixsize=0.1, frame=None, proj='ZEA', unit=''):
        """

        """
        if type(center)==str:
            center = SkyCoord.from_name(center)
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

        nx=ny=naxis = int(size/pixsize) | 1 # make odd so central pixel is in the middle
        self.center = center
        self.frame=frame
        self.unit=unit
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

    def plot_map(self, hmap, fig=None, axpos=111, figsize=(8,8),
             log=False, cmap='jet', colorbar=True, 
             unit='', vmin=None, vmax=None, cb_kw={},
             title=None, **kwargs):
        """
        Plot a map
        
        - hmap -- a HEALPix map
        - fig  [None] -- a Figure
        - axpos [111] 


        Return a WCSAxes object
        """


        from matplotlib import colors

        wcs = self
        grid = self._make_grid();
        nside = healpy.get_nside(hmap)

        ipix = healpy.ang2pix(nside, *grid, lonlat=True)

        fig = fig or plt.figure(figsize=figsize)
        if type(axpos)==tuple:
            ax = fig.add_subplot(*axpos, projection=self)
        else:
            ax = fig.add_subplot(axpos, projection=self)
        ax = fig.axes[-1] ####  KLUGE
        if log:
            norm = colors.LogNorm(vmin=vmin,vmax=vmax)
            vmin=vmax=None
        else:
            norm = None
        ipix = healpy.ang2pix(nside, *grid, lonlat=True)
        im = ax.imshow(hmap[ipix], cmap=cmap, origin='lower', 
                       norm=norm,  vmin=vmin, vmax=vmax);

        nx, ny = wcs.array_shape
        ax.set(xlabel=self.axis_labels[0], xlim=(-0.5, nx-0.5),
               ylabel=self.axis_labels[1], ylim=(-0.5, ny-0.5),
               title= title)
        ax.grid()
        if colorbar:
            ticklabels = cb_kw.pop('ticklabels', None)
            cb_kw.update(label=self.unit,)
            cb = plt.colorbar(im, ax=ax,  **cb_kw)
            if ticklabels is not None:
                cb.ax.set_yticklabels(ticklabels)
        # set the default transform for world
        ax.transAxes = ax.get_transform(self.frame)
            # annotator(ax, self.frame)
        return ax
    

class ZEAaxis:
    """
    Wraps WCSaxis to intercept plot, scatter and text

    """
    def __init__(self, ax):
        # ax is a WCSaxis object
        self.ax= ax

    def _to_pixel(self, *args):
        return _to_pixels(*args, wcs=self.ax.wcs)
        # return self.ax.wcs.world_to_pixel(sc) + tuple(rest) 
    
    def plot(self, *args, **kwargs):
        return self.ax.plot(*self._to_pixel(*args), **kwargs)
    
    # def scatter(self, *args, **kwargs):
    #     return self.ax.scatter(*self._to_pixel(*args), **kwargs)
    # def text(self, *args, **kwargs):
    #     return self.ax.text(*self._to_pixel(*args),  **kwargs)

    # def __getattr__(self, name):
    #     # pass everything else to self.ax
    #     return self.ax.__getattribute__( name)
    def scatter(self, *args, **kwargs):
        self.ax.scatter(*self._to_pixel(*args), **kwargs)
        return self
    
    def text(self, *args, **kwargs):

        pars = self._to_pixel(*args)
        if len(np.atleast_1d(pars[0]))==1:
            self.ax.text(*pars, **kwargs)
        elif len(pars)==3:
            crpix = self.wcs.wcs.crpix
 
            for x,y,text in zip(*pars): #, pars[1], pars[2]):
                # print(x,y, text)
                if x<0 or x>= crpix[0] or  y<0 or y>=crpix[1]: continue
                self.ax.text(x,y,text, **kwargs)
        else:
            raise Exception(f'text parameters not understood')
        # self.ax.text(*self._to_pixel(*args),  **kwargs)
        return self

    def __getattr__(self, name):
        # pass everything else to self.ax
        return self.ax.__getattribute__( name)
    
    def apply(self, func, *pars, **kwargs):
        """Pass the axis to a user-supplied function and return self."""
        func(self, *pars, **kwargs)
        return self    

class HPmap(object):
    """
    Manage HEALPix array
    """
    def __init__(self,
            hpmap:'HEALPix array',
            name='',
            cblabel=None,
            unit='',
            sigma:'smooth parameter'=None,
            nest=False):
        """create from a HEALPix array
        """
        self.name = name
        self.cblabel = cblabel if cblabel is not None else unit
        self.unit = unit
        self.nside = get_nside(hpmap)
        # reorder as defaut RING if nest is set
        self.map = hpmap if not nest else healpy.reorder(hpmap, n2r=True)

        if sigma is not None:
            self.smooth(sigma)

    def __str__(self):
        return f'<{self.__class__.__name__}>, name "{self.name}" nside {self.nside} unit "{self.unit}"'
    
    def __repr__(self): return str(self)

    def __call__(self, sc:'SkyCoord') -> 'value[s]':
        """
        Return value of corresponding pixel(s)
        """
        sp = sc.spherical # avoid dependence on specific frame
        skyindex = healpy.ang2pix(self.nside, sp.lon.deg, sp.lat.deg, lonlat=True)
        return self.map[skyindex]

    def get_all_neighbors(self, pix):
        return healpy.get_all_neighbours(self.nside, pix)

    def smooth(self, sigma):
            self.map = healpy.smoothing(self.map, np.radians(sigma))

    def ait_plot(self, **kwargs):
        """
        Invoke the function ait_plot to draw a representation

        """
        kw= dict(label=self.name, cblabel=self.unit,)
        kw.update(**kwargs)
        return ait_plot(self, **kw)
        
    def zea_plot(self, *args, size=10, pixelsize=0.1, 
                fig=None, figsize=(8,8), **kwargs):
        """
        Create a SquareWCS object that defines its WCS as a ZEA square centered on the position given by `args`, 
        which can be a source name.
        Return the ZEAaxis object for plotting
        """
        sc, _ = _transform_plot_args(*args)
        swcs = SquareWCS(sc, size, pixelsize)
        return ZEAaxis(
            swcs.plot_map(self.map, unit=self.unit, fig=fig, figsize=figsize, **kwargs)
        )


    def convolve(self, beam_window=None, sigma=0):
        """Convolve the map with a "beam", or PSF

        - beam_window: Legendre coefficients or None
        - sigma  Gaussian sigma in degrees

        return a HEALPix array
        """

        return healpy.alm2map(
            healpy.smoothalm( healpy.map2alm(self.map),
                              sigma=np.radians(sigma),
                              beam_window=beam_window
                            ),
            nside=self.nside
            )

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
            ('COORDSYS', 'GALACTIC', ''),
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
    
    @classmethod
    def from_FITS(cls, filename, field=0, *pars, **kwargs):
        with fits.open(filename) as hdus:
            header, data = hdus[1].header, hdus[1].data
        kw = dict(unit=header.get(f'BUNIT', ''), name=header.get(f'TTYPE{field+1}', ''))
        kw.update(**kwargs)
        return cls(data.field(field), *pars, **kw)