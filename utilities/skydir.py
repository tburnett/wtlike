"""
"""
import numpy as np
from astropy.coordinates import SkyCoord 
import healpy  


class SkyDir(object):
    """ Provide convenient interface to astropy's SkyCoord:
    * specific to fk5 and galacic use
    * added HEALPix capability
    """
    nside = None
    nest = False
    
    def __init__(self, ra, dec, frame='fk5'):
        self.coord = SkyCoord(ra,dec, unit='deg', frame=frame)
    
    def __str__(self):
        return str(self.coord)
    def __repr__(self): return str(self)

    def __len__(self): return len(self.coord)
    def __getitem__(self, index): return self.coord[index]
    
    @property
    def fk5(self):
        t=self.coord.fk5
        return t.ra.value, t.dec.value
    @property
    def galactic(self):
        t = self.coord.galactic
        return t.l.value, t.b.value

    @classmethod
    def from_healpix(cls, ipix, nside=None) -> 'SkyDir':
        if nside is None:
            nside = SkyDir.nside
        assert nside, 'Expect nside to be specified, or a class variable'
        t = healpy.pix2ang(nside, ipix, nest=SkyDir.nest, lonlat=True)
        return cls(t[0], t[1], frame='galactic')

    @classmethod
    def from_galactic(cls, l, b) -> 'SkyDir':
        return cls(l, b, frame='galactic')
    @classmethod
    def gal(cls, l, b) -> 'SkyDir':
        return cls(l, b, frame='galactic')

    def to_healpix(self, nside=None) ->'List of pixels':
        if nside is None:
            nside = SkyDir.nside
        assert nside, 'Expect nside to be specified, a class variable'
        l,b = self.galactic
        return healpy.ang2pix(nside, l,b, nest=SkyDir.nest, lonlat=True)

    def match(self, catalog:'another SkyDir with more than one position'
            )->'idx, sepdeg':

        idx, sep, dist = self.coord.match_to_catalog_sky(catalog.coord)
        return idx, sep.deg

def test():
    SkyDir.nside=64

    r = np.array(range(10))
    t = SkyDir.from_healpix(r); 
    u = SkyDir.from_galactic(*t.galactic);

    v = u.to_healpix()
    assert np.all(r==v), 'Failed healpy test'
    
