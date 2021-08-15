"""
# wtlike--make FermiLAT gamma-ray Light curves

see https://tburnett.github.io/wtlike/tutorial
"""

__version__ = "0.3.4"

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import Config, UTC, MJD
from .data_man import check_data
from .sources import PointSource
from .loglike import PoissonRep
from .simulation import Simulation
from .main import WtLike