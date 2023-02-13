import numpy as np
import healpy
from wtlike.skymaps import AitoffFigure


class Clusterer:
    
    def cluster(self, indices, min_size=2):

        def grow(indeces):

            def neighbors(i,j):
                iv = healpy.get_all_neighbours(self.nside, i)
                return j in iv
            i = 0
            cluster =list(indeces[0:1])
            remainder = indeces[1:]
            while i<len(cluster):
                newremain = []
                for j in remainder:
                    if neighbors(cluster[i],j): cluster.append(j)
                    else: newremain.append(j)
                i += 1
                remainder = newremain
            return cluster, remainder

        if not self.quiet:
            print( f'Clustering {len(indices)} pixels...')
        ret = []
        rem = indices
        while len(rem)>0:
            clu, rem = grow(rem)
            if len(clu)> min_size:
                ret.append(np.array(clu))
        if not self.quiet:
            print( f'Found {len(ret)} clusters')
        return ret

    def __init__(self, cmap, min_size=2,  threshold=0, quiet=False):
        """ Make a list of clusters from a HEALPix map of pixels

        cmap -- the map, RING ordering
        min_size [2] -- minimum cluster size
        theeshold [0] 

        """
        self.nside = healpy.get_nside(cmap)
        self.cmap = cmap 
        self.quiet = quiet
        
        self.indices = np.arange(len(cmap))[cmap>threshold]
        self.clusters = self.cluster(self.indices, min_size=min_size)
        
    def plotit(self):
        
        afig = AitoffFigure(figsize=(10,5))

        lon = []; lat=[]; s=[]
        for c in self.clusters:
            w = self.cmap[c]
            vec =np.array(healpy.pix2vec(self.nside, c))
            wvec = (vec*w).mean(axis=1)/w.mean()
            l,b = healpy.vec2dir(wvec, lonlat=True)
            lon.append(l); lat.append(b)
            s.append(sum(w))
            
        afig.scatter(lon,lat, 4*np.sqrt(np.array(s)), alpha=0.6, color='r')
        