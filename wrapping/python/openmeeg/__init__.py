from . import _distributor_init

from .openmeeg import *
from ._version import __version__
from ._make_geometry import make_geometry, make_nested_geometry, read_geometry
from ._utils import get_log_level, set_log_level, use_log_level

set_log_level("warning")
