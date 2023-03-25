import warnings

import matplotlib as mpl

from pymsaviz.config import get_msa_testdata
from pymsaviz.msaviz import MsaViz

warnings.filterwarnings("ignore")

__all__ = [
    "MsaViz",
    "get_msa_testdata",
]

__version__ = "0.4.0"

# Setting matplotlib rc(runtime configuration) parameters
# https://matplotlib.org/stable/tutorials/introductory/customizing.html
mpl_rc_params = {
    # Legend
    "legend.loc": "upper left",  # Default: best
    "legend.frameon": False,  # Default: True
    "legend.handlelength": 1,  # Default: 2.0
    "legend.handleheight": 1,  # Default: 0.7
    # Savefig
    "savefig.bbox": "tight",  # Default: None
    "savefig.pad_inches": 0.5,  # Default: 0.1
    # SVG
    "svg.fonttype": "none",
}
mpl.rcParams.update(mpl_rc_params)
