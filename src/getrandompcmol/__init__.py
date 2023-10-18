"""
GetRandomPCMol
=======

Dummy command line tool to square a number.
"""

from .__version__ import __version__
from .evaluate_conf import eval_conf_ensemble
from .miscelleanous import bcolors, chdir, checkifinpath, create_directory
from .qmcalc import xtbopt
from .randommols import console_entry_point, main
