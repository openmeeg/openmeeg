from . import _distributor_init

# Here we import as few things as possible to keep our API as limited as
# possible
from ._openmeeg_cxx import HeadMat, Sensors, Integrator, Head2EEGMat, Head2MEGMat
from ._openmeeg_cxx import DipSourceMat, DipSource2MEGMat
from ._openmeeg_cxx import (
    GainEEG,
    GainMEG,
    GainEEGadjoint,
    GainMEGadjoint,
    GainEEGMEGadjoint,
)
from ._openmeeg_cxx import SurfSourceMat, SurfSource2MEGMat
from ._openmeeg_cxx import Matrix, SymMatrix
from ._version import __version__
from ._make_geometry import make_geometry, make_nested_geometry, read_geometry
from ._utils import get_log_level, set_log_level, use_log_level

set_log_level("warning")
