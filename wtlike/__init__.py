"""
# wtlike--make FermiLAT gamma-ray Light curves

see https://tburnett.github.io/wtlike/tutorial
"""

__version__ = "0.6.4"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import Config, UTC, MJD, Timer, FermiInterval
from .data_man import check_data, update_data, DataView
from .sources import PointSource, FermiCatalog
from .loglike import PoissonRep
from .simulation import Simulation
from .time_series import TimeSeries
from .main import WtLike
