"""
# wtlike--make FermiLAT gamma-ray Light curves

see https://tburnett.github.io/wtlike/tutorial
"""

__version__ = "0.9.1"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path

from .config import Config, UTC, MJD, Timer, FermiInterval
from .data_man import check_data, update_data, DataView
from .sources import PointSource, SourceFinder
from .exposure import WeightedAeff
from .loglike import PoissonRep
from .simulation import Simulation
from .time_series import TimeSeries
from .main import WtLike
