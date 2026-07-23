import numpy as np
import healpy
from wtlike.skymaps import AitoffFigure


def grouper(points, radius, min_size=1, return_indices=False):
    """Group points into connected clusters using a distance threshold.

    Parameters
    ----------
    points : array-like, shape (n_points, n_dims)
        Point coordinates.
    radius : float
        Maximum pairwise distance for two points to be considered neighbors.
    min_size : int, default=1
        Minimum number of points required to keep a cluster.
    return_indices : bool, default=False
        If True, return index arrays for each cluster; otherwise return points.

    Returns
    -------
    list[np.ndarray]
        A list of clusters, sorted by descending cluster size.
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2:
        raise ValueError("points must be a 2D array with shape (n_points, n_dims)")
    if radius <= 0:
        raise ValueError("radius must be > 0")
    if min_size < 1:
        raise ValueError("min_size must be >= 1")

    n_points = len(pts)
    if n_points == 0:
        return []

    r2 = float(radius) ** 2
    unvisited = np.ones(n_points, dtype=bool)
    clusters = []

    # Build connected components where edges join points within `radius`.
    for seed in range(n_points):
        if not unvisited[seed]:
            continue

        stack = [seed]
        unvisited[seed] = False
        cluster = []

        while stack:
            i = stack.pop()
            cluster.append(i)

            remaining = np.flatnonzero(unvisited)
            if len(remaining) == 0:
                continue

            d2 = ((pts[remaining] - pts[i]) ** 2).sum(axis=1)
            neighbors = remaining[d2 <= r2]
            if len(neighbors):
                unvisited[neighbors] = False
                stack.extend(neighbors.tolist())

        if len(cluster) >= min_size:
            clusters.append(np.array(cluster, dtype=int))

    clusters.sort(key=len, reverse=True)
    if return_indices:
        return clusters
    return [pts[idx] for idx in clusters]


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
        