from . import _distributor_init

# Here we import as few things as possible to keep our API as limited as
# possible

from ._openmeeg_cxx import (
    Matrix,
    SymMatrix,
    SparseMatrix,

    HeadMat,
    Sensors,

    CorticalMat,
    CorticalMat2,
    DipSource2InternalPotMat,
    DipSource2MEGMat,
    DipSourceMat,
    EITSourceMat,
    Forward,
    GainEEG,
    GainMEG,
    GainEEGadjoint,
    GainMEGadjoint,
    GainEEGMEGadjoint,
    GainEITInternalPot,
    GainInternalPot,
    Head2ECoGMat,
    Head2EEGMat,
    Head2MEGMat,
    SurfSourceMat,
    SurfSource2MEGMat,
    Surf2VolMat,

    Integrator,
)

from ._version import __version__
from ._make_geometry import make_geometry, make_nested_geometry, read_geometry
from ._utils import get_log_level, set_log_level, use_log_level

set_log_level("warning")
