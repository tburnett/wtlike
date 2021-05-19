# AUTOGENERATED! DO NOT EDIT! File to edit: nbs/03_weights.ipynb (unless otherwise specified).

__all__ = ['check_weights', 'load_weights', 'add_weights']

# Cell
import os, sys,  pickle, healpy
import numpy as np
from .config import *
#from wtlike.photon_data import *



# Cell
def check_weights(config, source):
    """
    Check that weights for the source are available: if so, return the weight file name

    - source -- A PointSource object with information on source location

    Returns the filepath to the file if successful, otherwise, print a message abount available files
    """
    weight_files = config.wtlike_data/'weight_files'
    assert weight_files.is_dir(), f'Expect {weight_files} to be a directory'
    weight_file = weight_files/ (source.filename+'_weights.pkl')
    if not weight_file.exists():
        available = np.array(list(map(lambda p: p.name[:p.name.find('_weights')],
                          weight_files.glob('*_weights.pkl'))))
        print(f'{source} not found in list of weight files at\n\t {weight_files}.\n Available:\n{available}',
             file = sys.stderr)
        return None
    return weight_file

# Cell
def load_weights(config, filename, ):
    """Load the weight informaton

    filename: pickled dict with map info

    """
    # load a pickle containing weights, generated by pointlike
    assert os.path.exists(filename),f'File {filename} not found.'
    with open(filename, 'rb') as file:
        wtd = pickle.load(file, encoding='latin1')
    assert type(wtd)==dict, 'Expect a dictionary'
    test_elements = 'energy_bins pixels weights nside model_name radius order roi_name'.split()
    assert np.all([x in wtd.keys() for x in test_elements]),f'Dict missing one of the keys {test_elements}'
    if config.verbose>0:
        print(f'Load weights from file {os.path.realpath(filename)}')
        pos = wtd['source_lb']
        print(f'\tFound: {wtd["source_name"]} at ({pos[0]:.2f}, {pos[1]:.2f})')
    # extract pixel ids and nside used
    wt_pix   = wtd['pixels']
    nside_wt = wtd['nside']

    # merge the weights into a table, with default nans
    # indexing is band id rows by weight pixel columns
    # append one empty column for photons not in a weight pixel
    # calculated weights are in a dict with band id keys
    wts = np.full((32, len(wt_pix)+1), np.nan, dtype=np.float32)
    weight_dict = wtd['weights']
    for k in weight_dict.keys():
        t = weight_dict[k]
        if len(t.shape)==2:
            t = t.T[0] #???
        wts[k,:-1] = t
    return wts , wt_pix , nside_wt

# Cell
def _add_weights(config, wts, wt_pix, nside_wt, photon_data):
    """ get the photon pixel ids, convert to NEST (if not already) and right shift them
        add 'weight', remove 'band', 'pixel'
    """
    if not config.nest:
        # data are RING
        photon_pix = healpy.ring2nest(config.nside, photon_data.pixel.values)
    else:
        photon_pix = photon_data.pixel.values
    to_shift = 2*int(np.log2(config.nside/nside_wt));
    shifted_pix =   np.right_shift(photon_pix, to_shift)
    bad = np.logical_not(np.isin(shifted_pix, wt_pix))
    if config.verbose>0 & sum(bad)>0:
        print(f'\tApplying weights: {sum(bad)} / {len(bad)} photon pixels are outside weight region')
    if sum(bad)==len(bad):
        a = np.array(healpy.pix2ang(nside_wt, wt_pix, nest=True, lonlat=True)).mean(axis=1).round(1)
        b = np.array(healpy.pix2ang(nside_wt, shifted_pix, nest=True, lonlat=True)).mean(axis=1).round(1)

        raise Exception(f'There was no overlap of the photon data at {b} and the weights at {a}')
    shifted_pix[bad] = 12*nside_wt**2 # set index to be beyond pixel indices

    # find indices with search and add a "weights" column
    # (expect that wt_pix are NEST ordering and sorted)
    weight_index = np.searchsorted(wt_pix,shifted_pix)
    band_index = np.fmin(31, photon_data.band.values) #all above 1 TeV into last bin

    # final grand lookup -- isn't numpy wonderful!
    photon_data.loc[:,'weight'] = wts[tuple([band_index, weight_index])]

    # don't need these columns now (add flag to config to control??)
    photon_data.drop(['band', 'pixel'], axis=1)

    if config.verbose>1:
        print(f'\t{sum(np.isnan(photon_data.weight.values))} events without weight')


# Cell
def add_weights(config,  photon_data, source, nbins=50):
    """ add weights for the source to the photon data

    - photon_data -- DataFrame with photon data

    - source -- `PointSource` object

    Return the weight value histogram
    """
    weight_file =  check_weights(config,  source)
    if weight_file is None:
        raise Exception(f'Weight file not found for {source}')

    wts, wt_pix, nside_wt = load_weights(config, weight_file)
    _add_weights(config, wts, wt_pix, nside_wt, photon_data)

    return np.histogram(photon_data.weight.values, np.linspace(0,1,nbins+1))[0]

# Cell
# def get_weight_hist(config,  source, nbins=50, key=''):
#     """ return a weight distribution

#     - photon_data -- DataFrame with photon data
#     - source -- `PointSource` object

#     Uses `add_weights`.
#     """
#     def doit(nbins):
#         weight_file =  check_weights(config,  source)
#         if weight_file is None:
#             raise Exception(f'Weight file not found for {source}')
#         photon_data = get_photon_data(config,  source )
#         return add_weights(config, photon_data, source, nbins=nbins)

#     key = f'weight_hist_{source.name}' if key=='' else key
#     description = f'Weight histogram for {source.name}' if config.verbose>1 else ''
#     return config.cache(key, doit, nbins, description=description)
