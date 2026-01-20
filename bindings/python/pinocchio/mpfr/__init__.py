#
# Copyright (c) 2021 INRIA
#
# ruff: noqa: F401, F403, F405

# Manually register submodules
import sys

from .. import utils
from ..explog import exp, log
from ..pinocchio_pywrap_mpfr import *
from ..pinocchio_pywrap_mpfr import __raw_version__, __version__

sys.modules["pinocchio.mpfr.rpy"] = rpy
sys.modules["pinocchio.mpfr.cholesky"] = cholesky

if WITH_COLLISION:
    import coal
    from coal import (
        CachedMeshLoader,
        CollisionGeometry,
        CollisionResult,
        Contact,
        DistanceResult,
        MeshLoader,
        StdVec_CollisionResult,
        StdVec_Contact,
        StdVec_DistanceResult,
    )

    # Deprecated, should be removed in next major release
    hppfcl = coal
