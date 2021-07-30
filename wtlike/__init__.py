"""
# wtlike--make FermiLAT gamma-ray Light curves

see https://tburnett.github.io/wtlike/tutorial
"""

__version__ = "0.2.4"


from .config import Config, UTC, MJD
from .data_man import check_data
from .sources import PointSource
from .loglike import PoissonRep
from .simulation import Simulation
from .main import WtLike