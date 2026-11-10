"""Shared simulation settings."""

from __future__ import annotations

from functools import partial

from mxlpy.integrators import Native
from mxlpy.runners import KineticRunner

NATIVE = partial(Native, method="CVODE_BDF")

# Native's autocorrelation oscillation detector flags monotonically settling
# trajectories of the PETC and metabolism models as oscillatory (false positives,
# ADR 0004), so steady-state scans run without it.
NATIVE_RUNNER = KineticRunner(integrator=NATIVE, oscillation_detector=None)
