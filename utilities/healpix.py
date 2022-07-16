import os, sys
from matplotlib import units
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib import colors
from astropy.io import fits
from astropy.coordinates import SkyCoord
import healpy

from . skydir import SkyDir


class HPmap(object):
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

    def __call__(self, skydir:'SkyDir') -> 'value[s]':
        """
        """
        skyindex = skydir.to_healpix(nside=self.nside)
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

    @classmethod
    def from_HPmap(cls, 
            other:'an existing HPmap object to resample', 
            nside,
            name:'override name of other'='',
            cblabel=''):
        skydir = SkyDir.from_healpix(range(12*nside**2))
        return cls(other(skydir), name=name or other.name, cblabel=cblable or other.cblabel)


    def ait_plot(self,  **kwargs):
        kw= dict(label=self.name, cblabel=self.cblabel,)
        kw.update(**kwargs)
        return ait_plot(self, **kw)

    def __truediv__(self, other ):
        return self.map / other.map

    def __mul__(self, other):
        return self.map * other.map

    def to_FITS(self,  filename=''):
        """return a HDUlist object with one skymap column
        
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


class HPcube(HPmap):
    """Implement a spectral cube with logarithmic interpolation
    """

    def __init__(self, spectral_cube:'2-d array, or list of 1-d arrays', 
                    energies:'corresponding energies',
                    cblabel='', 
                    energy:'default value'=1000,
                    unit:'units'='',
                    sigma:'somothing parameter (deg)'=0, 
                     ):
        sc = spectral_cube
        if type(sc)==list:
            assert len(sc)==len(energies), 'Input was a list of wrong length'
            spectral_cube = np.vstack(sc)
        sh = spectral_cube.shape
        assert len(sh)==2, 'Expect 2-d array'
        # transpose if needed               
        self.spectra = spectral_cube.T if sh[0]>sh[1] else spectral_cube
        
        try:
            self.nside = healpy.get_nside(self.spectra[0])
        except TypeError:
            print(f'shape {sh} not compatible with HEALPix', outfile=sys.stderr)
            raise
        if sigma>0:
            # use the healpy smoothing function
            self.spectra = np.array(list(map(lambda m: healpy.smoothing(m,np.radians(sigma),
                    verbose=False),  self.spectra)))
        self.maps = list(map(HPmap, self.spectra))
        self.energies = energies
        self.unit = unit
        if not cblabel and unit:
            unit = unit.replace(' ', '\\ ')
            self.cblabel = fr'$\mathrm{{ {unit} }}$'
        else:
            self.cblabel=cblabel
        self.loge = np.log(energies) # for logarithmin interpolation
        self.set_energy(energy)
        if np.std(np.diff(self.energies))<1:
            print('Warning: this is not a spectral cube, which should use logarithmic interpolation')
     
    def __str__(self):
        ee = self.energies
        return f'<{self.__class__.__name__}>: nside {self.nside}, '\
               f' {len(ee)} energies {ee[0]:.2e}-{ee[-1]:.2e} MeV '

    def __getitem__(self, index):
        return self.maps[index]
    def __len__(self):
        return len(self.maps)

    def __iter__(self): # iterate over all keys
        for index in range(len(self)):
            yield self[index]

    @property
    def spectral_cube(self):
        return np.array([layer.map for layer in self])

    def ait_plot(self, index, **kwargs):
        if type(index)==int:
            if index<len(self):
                # small index: just get plane by index
                self[index].ait_plot(  **kwargs)
                return
        # interpret as an energy
        self.set_energy(float(index))
        kw = dict(log=True, cblabel=getattr(self, 'cblabel', ''))
        kw.update(**kwargs)
        ait_plot(self, self.energy, **kw)        

    def set_energy(self, energy): 
        # set up logarithmic interpolation

        self.energy=energy
        # get the pair of energies
        if energy< self.energies[0]: i=0
        elif energy>self.energies[-1]: i= len(self.energies)-2
        else:
            i = np.where(self.energies>=energy)[0][0]-1
         
        a,b = self.loge[i], self.loge[i+1]
        self.energy_index = i #= max(0, min(int(r), len(self.energies)-2))
        self.energy_interpolation = (np.log(energy)-a)/(b-a)

        self.eplane1 = self[i].map
        self.eplane2 = self[i+1].map

    def __call__(self, skydir:'SkyDir', energy=None) -> 'value[s]':
        """
        """
        skyindex = skydir.to_healpix(nside=self.nside)
        if energy is None:
            ee = [self.energy]
        else:
            ee = np.atleast_1d(energy)
        ret = []
        for energy in ee:    
            self.set_energy(energy)    
            a = self.energy_interpolation
            u, v = self.eplane1[skyindex], self.eplane2[skyindex]
            # avoid interpolation if close to a plane
            if np.abs(a) < 1e-2:  # or v<=0 or np.isnan(v):
                w = u
            elif np.abs(1-a)< 1e-2: # or u<=0 or np.isnan(u):
                w = v
            else:
                w = np.exp( np.log(u) * (1-a) + np.log(v) * a  )
            ret.append(w)
        return ret[0] if len(ee)==1 else np.array(ret)

    def hpmap(self, energy, label=''):
        """ return an HPmap for the given energy
        """
        self.set_energy(energy)
        a = self.energy_interpolation
        if a<0.002:
            map = self.eplane1
        elif a>0.998:
            map = self.eplane2
        else:
            map = np.exp( np.log(self.eplane1) * (1-a) 
                + np.log(self.eplane2) * a   )
        return HPmap(map, label=label, cblabel=self.cblabel)

    @classmethod
    def from_FITS(cls, filename):
        try:
            hdus = fits.open(filename)
        except Exception  as msg:
            raise Exception('FITS: Failed to open {}: {}'.format(filename, msg))
        try:
            if len(hdus)==2:
                energies = []
                print(f'No energy table in file {filename}')
            else:
                if hdus[2].columns[0].name=='CHANNEL':
                    # binned format: assume next 2 columns are min, max and use geometric mean
                    emin,emax = [hdus[2].data.field(i) for i in (1,2)]
                    energies = np.sqrt(emin*emax)
                else:
                    energies = hdus[2].data.field(0)
            hdu1 = hdus[1]
            data = hdu1.data
            vector_mode = len(hdu1.columns)==1
            
            if vector_mode:
                # one vector column, expect 2d array with shape (12*nside**2, len(energies))
                spectral_cube = hdus[1].data.field(0).T
            else:
                # one column per energy: expect len(energies) columns
                spectral_cube = np.vstack([col.array for col in data.columns])
 
            nside = int(np.sqrt(spectral_cube.shape[0]/12.))
            assert spectral_cube.shape[0]==len(energies), 'shape inconsistent with number of energies'
     
            unit = hdu1.header.get('BUNIT', '')

            assert hdu1.header.get('ORDERING','RING')=='RING', 'Wrong ordering'
            assert hdu1.header.get('COORDSYS', 'GAL')=='GAL', 'Wrong coordsys'
                
        except Exception as msg:
            print(f'bad file or unexpected FITS format, file {filename}: {msg}')
            raise
        hdus.close()
        return cls(spectral_cube, energies, unit=unit) 

    @classmethod
    def from_cube(cls, 
            other:'anothter HPcube',
            nside:'new nside or get from other'=None, 
            sigma:'smoothing sigma in deg'=0,
            energies:'energies or get from other'=None,
            ):
        """create from an existing HPcube to change nside, energies, and/or smooth
        """
        assert isinstance(other, HPcube), 'Expected an HPcube object'
        
        nside = nside or other.nside
        energies = energies or other.energies
        sd_list = sd = SkyDir.from_healpix(range(12*nside**2), nside=nside)
        resampled = other(sd_list, energies)
        return cls(resampled, energies, sigma=sigma) 





    def to_FITS(self, outfile, 
            vector_format:'if True, one-column format'=False, 
            overwrite=True):

        def spectral_table(self):
            
            el = self.energies

            if vector_format:
                columns = [fits.Column(name='spectra', format=f'{len(el)}E',
                            unit=getattr(self, 'unit',''), 
                            array=np.vstack([x.map for x in self]).T)
                            ]
            else:
                columns = [fits.Column(name=f'CHANNEL{i:02d}', format='E', unit='', array=plane.map)
                    for i, plane in enumerate(self)]
                
            table = fits.BinTableHDU.from_columns(columns)
            table.name = 'SKYMAP' 
            # add HEALPix and energy info to the header 
            nside = self.nside

            emin, deltae= el[0], np.log(el[1]/el[0])
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
                    ('NRBINS',   len(el),      'Number of energy bins'),
                    ('COORDSYS', 'GAL', ''),
                    ('EMIN',     emin,           'Minimum energy'  ),
                    ('DELTAE',   deltae,         'Step in energy (log)'),
                    ('BUNIT',    getattr(self,'unit', ''), ''),
                ]]
            for card in cards: table.header.append(card)
            return table

        def energy_table(self):
            column= fits.Column( name='energy', format='E', unit='MeV', array=self.energies)
            table = fits.BinTableHDU.from_columns([column])
            table.name='ENERGIES'
            return table

        def primary(self):
            cards = [fits.Card('COMMENT', 'Written by utilities/healpix.py '),
                    ]
            return fits.PrimaryHDU(header=fits.header.Header(cards))

        # def hdus(self):
        #     return [ primary(self),  #primary
        #             spectral_table(self),           # this table
        #             energy_table(self),
        #         ]
      

        hdu_list = fits.HDUList( 
                [   primary(self),  #primary
                    spectral_table(self),           # this table
                    energy_table(self),
                ])
        hdu_list.writeto(outfile, overwrite=overwrite)
        hdu_list.close()
        print( f'\nwrote {"vector format" if vector_format else ""} '\
            f'FITS Skymap file, nside={self.nside}, {len(self.energies)} energies, to {outfile}')


    # def ratio(self, other, skydir, energy, label=''):
    #     """return the ratio of this with another cube at the energy
    #     """
    #     return self(skydir,energy) / other(skydir,energy)

class HPcubeOp(HPcube):

    """ Base class for operations 
    """
    def __init__(self, hpmap1, hpmap2):
        self.map1=hpmap1
        self.map2=hpmap2
        self.set_energy(1000.)

    def set_energy(self, energy):
        self.energy=float(energy)


class HPratio(HPcubeOp):
    """Ratio of two HpMap obects
    """    
    def __call__(self, skydir:'SkyDir', 
                energy) -> 'value[s]':
        return self.map1(skydir,float(energy)) / self.map2(skydir,float(energy))


class HPproduct(HPcubeOp):
    """ Product of two HPmap objecs
    """

    def __call__(self, skydir:'SkyDir', 
                energy) -> 'value[s]':
        self.set_energy(energy)
        return self.map1(skydir,self.energy) * self.map2(skydir,self.energy)

class AitoffFigure():

    """ Implement plot and text converting from (l,b) in degrees, or a SkyCoord.

    """    
    def __init__(self, fig):
        self.ax=fig.axes[0]
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
        fig=None, ax=None, fignum=1, figsize=(12,5),
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
    dirs = SkyDir.from_galactic(Lon, Lat)
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


def ait_multiplot(
        mappable:'either an HPCube, or list of HPmap objects', 
        energies:'energies to evaluate the HPcube'=None, 
        labels:'list of labels'=[], 
        fignum=None,
        sizex=14,
        sizey=None,
        cmap='jet', 
        vmin=None, vmax=None, 
        log=False,
        colorbar=True,
        cb_shrink=0.6,
        nx:'number of columns'=2,
        title:'used for suptitle'='',         
        **kwargs ):

    if energies is not None:
        n = len(energies)
    else:
        mappable = np.atleast_1d(mappable)
        n = len(mappable)
    if len(labels)==0:
        labels=['']*n    
 
    ny = (n+nx-1)//nx
    sizey = sizey or  (0.5 + 3.5*ny)  #[4.0, 7.5, 11, 14.5, 18][ny-1]
    kw = dict(cmap=cmap, vmin=vmin, vmax=vmax, log=log, colorbar=colorbar,cb_kw=dict(shrink=cb_shrink), **kwargs)
    
    fig, axx = plt.subplots(ny,nx, figsize=(sizex, sizey), num=fignum,
            gridspec_kw={'wspace':0.05, 'hspace':0.0,'left':0.2, 'top':0.9},
            subplot_kw=dict(projection='aitoff',visible=False ))
    
    #fig.tight_layout()
    
    if energies is not None:
        for ax, energy, label in zip(axx.flatten(), energies, labels):
            ait_plot(mappable,  energy,  ax=ax,  label=label, **kw)
            # cbtext = kwargs.get('cbtext', kwargs.get('cblabel', '')) 

    else:
        for m, ax, label in zip(mappable, axx.flatten(), labels):
            ait_plot(m,  ax=ax,  label=label, **kw)

    if title: fig.text(0.4, 0.92, title, fontsize=16, ha='center');
    fig.set_facecolor('white')
    return fig
        

def ait_multiplot_bands(mappable, fignum=1, **kwargs):
    """Special for the 8 UW energy bands
    """
    band_energies = np.logspace(2.125, 3.875, 8)
    band_labels = list(map(lambda x: f'{x/1e3:.2f} GeV', band_energies))
    fs =mpl.rcParams['font.size']
    plt.rc('font', size=10)
    ret= ait_multiplot( mappable,  band_energies,  labels=band_labels,
            nx=4, cb_shrink=0.5, sizex=20,sizey=5.,  
             **kwargs)
    mpl.rcParams['font.size']= fs
    return ret


class Polyfit(object):
    """ Manage a log polynormal fit to each pixel
    """
    def __init__(self, 
            cubefile:'either a filename or a spectral cube',
            sigsfile=None, start=0, stop=8, deg=2, limits=(0.5,25)):
        """
        """
        
        if type(cubefile)==str:
            m = HPcube.from_FITS(cubefile)
        else:
            m = cubefile
        meas = m.spectral_cube
        self.nside=m.nside
        self.limits=limits
        
        if sigsfile:
            msig = HPcube.from_FITS(sigsfile)
            sig  = msig.spectral_cube # if msig else np.ones(meas.shape)#np.array([msig[i].map for i in range(8)])
            weights = 100./sig[start:,:] #from percent
            self.wtmean = weights.mean(axis=1)
        else:
            # actual numbers 
            self.wtmean = np.array([ 66.4 ,92.33, 129.35, 121.26,  93.3 ,  64.26,  41.42,  25.74])

        self.planes = np.array(range(start,stop)) # plane numbers
        self.values = meas[start:,:]

        self.fit, self.residuals, self.rank, self.svals, self.rcond =\
            np.polyfit(self.planes,self.values, deg=deg, full=True, w=self.wtmean)
            
        labels= 'intercept slope curvature'.split()   
        
    @classmethod
    def from_product(cls, cf1, cf2, 
                nside=64, 
                log_energy_range=(2.125, 3.875, 8)):
        """ Product of two PolyFit objects
        """
        # make a spectral cube of the two PolyFit guys, assumog nside and energies
        nside=64
        coords  = SkyDir.from_healpix(range(12*nside**2), nside)
        energies = np.logspace(*log_energy_range)
        cube = np.vstack([cf1(coords, energy)*cf2(coords,energy) for energy in energies ])

        return cls(HPcube(cube, energies=energies))

    def __getitem__(self, i):
        return self.fit[i]

    def energy_index(self, energy):
        """convert energy in MeV to correspond to np.logspace(2.125, 3.875, 8)
        """
        energy = np.atleast_1d(energy)
        return 4*np.log10(energy/100)-0.5
    

    def __call__(self, 
                coord: 'SkyDir', 
                energy:'energies in MeV',
                factor_cap=25,
                
                )->'interpolated list of factors':

        """
        Implements the sky function using the logparabola fit
        Note that either arg can be multivalued, but not both.
        """
        energy = np.atleast_1d(energy)
        pix = np.atleast_1d(coord.to_healpix(self.nside) )        
        x = self.energy_index(energy)   
        fit= self.fit[:,pix]; 
        xx= np.vstack([x**2, x, np.ones(len(x))] )
        ret = np.matmul(fit.T, xx).T
        ret = np.clip(ret, *self.limits)
        return ret[0] if ret.shape[0]==1 else ret

    def parabolic_fit(self, 
            x:'energy index', 
            pix:'pixel index'):
        """Evaluate the parabolic fit in energy index
        """
        # if not hasattr(x, '__iter__'):
        #     x = np.array([x])
        x = np.atleast_1d(x)
        fit= self.fit[:,pix]; 
        t =fit.reshape(3,1)
        return ( t * np.vstack([x**2, x, np.ones(len(x))] )).sum(axis=0)
    
    def ang2pix(self, glon, glat):
        return healpy.ang2pix(64, glon, glat, lonlat=True)
        
    def get_fit(self, pix):
             
        y = self.values[:,pix]
        yerr = 1/self.wtmean
        fn = lambda xx : self.parabolic_fit(xx, pix)
        return y, yerr, fn
    
    def plot_fit(self, glon, glat, ax=None, axis_labels=True):

        pix = self.ang2pix(glon, glat)
        y, yerr, fn = self.get_fit(pix)

        fig, ax =plt.subplots() if ax is None else (ax.figure, ax)
        npl = len(self.planes)
        xx = np.linspace(self.planes[0]-0.5,self.planes[-1]+0.5,2*npl+1)

        ax.errorbar(self.planes, y, yerr=yerr, fmt='o', ms=8, label='measured values' if axis_labels else '');
        ax.plot(xx, fn(xx), '-', lw=2, label='parabolic fit' if axis_labels else '');
        ax.text(0.05,0.9,'({:3.0f}, {:+2.0f})'.format(glon, glat), transform=ax.transAxes)
        if axis_labels:
            ax.set(ylabel='flux factor', xlabel='Energy (GeV)')
            ax.set(xticks=[-0.5, 3.5, 7.5], xticklabels='0.1 1.0 10'.split())
            ax.legend()
        else:
            ax.set_xticks(self.planes[::2])  
        ax.axhline(1.0, color='grey')
        ax.grid(alpha=0.3)
        
    def multiplot(self, glons, glats, grid_shape=(4,5), title=''):
 
        fig, axx = plt.subplots(grid_shape[0],grid_shape[1], figsize=(12,12), sharex=True, sharey=True,
                            gridspec_kw=dict(left=0.05, right = 0.95,top=0.95, wspace=0, hspace=0)  )
        for glon, glat, ax in zip(glons, glats, axx.flatten()):
            self.plot_fit( glon, glat, ax, axis_labels=False)
        fig.suptitle(title, fontsize=16); 
        fig.set_facecolor('white')
           
    def ait_plots(self, fignum=1, title=''):
        fig, axx = plt.subplots(2,2, figsize=(14, 7.5), num=None,
                gridspec_kw={'wspace':0.05, 'hspace':0.0,'left':0.0, 'top':0.9},
                subplot_kw=dict(projection='aitoff') )#, visible=False), )

        for m, ax, label in zip(self, axx.flatten(), 
                'curvature  slope intercept '.split()):
           # ax.set_visible(True)
            ait_plot(HPmap(m, label=label) ,ax=ax, label=label,cb_kw=dict(shrink=0.7) )
        ait_plot(HPmap( self.residuals,label='residuals'),  ax=axx[1,1],label='residuals')
        if title: 
            fig.text(0.4, 0.92, title, fontsize=16, ha='center')
        fig.set_facecolor('white')
        return fig