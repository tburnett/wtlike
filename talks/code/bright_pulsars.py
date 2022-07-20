from wtlike import *
from astropy.coordinates import SkyCoord
from wtlike.skymaps import *
from cache_decorator import Cache
from utilities.ipynb_docgen import *

%run ../tmp//exposure/systematics.py
%run ../../../pulsars/code/sensitivity.py

from wtlike.time_series import power_spectrum_plot