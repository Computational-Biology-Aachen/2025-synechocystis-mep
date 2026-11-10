"""Shared simulation settings."""

from __future__ import annotations

from functools import partial

from mxlpy.integrators import Native, Sundials
from mxlpy.runners import KineticRunner

NATIVE = partial(Native, method="CVODE_BDF")

# Native's autocorrelation oscillation detector flags monotonically settling
# trajectories of the PETC and metabolism models as oscillatory (false positives,
# ADR 0004), so steady-state scans run without it.
NATIVE_RUNNER = KineticRunner(integrator=NATIVE, oscillation_detector=None)

# Native's steady-state search gives up on very slowly relaxing systems (e.g. the
# metabolism layer at near-zero supply), where Sundials converges. Used as fallback.
SUNDIALS_RUNNER = KineticRunner(integrator=Sundials, oscillation_detector=None)


def steady_state_scan(model, to_scan):
    """Steady-state scan with the native runner, retrying failures with Sundials."""
    import pandas as pd
    from mxlpy import scan

    res = scan.steady_state(
        model, to_scan=to_scan, parallel=False, runner=NATIVE_RUNNER
    )
    args, fluxes = (
        res.get_args().reset_index(drop=True),
        res.fluxes.reset_index(drop=True),
    )
    failed = args.isna().any(axis=1).to_numpy()
    if failed.any():
        retry = scan.steady_state(
            model,
            to_scan=to_scan[failed].reset_index(drop=True),
            parallel=False,
            runner=SUNDIALS_RUNNER,
        )
        args.loc[failed] = retry.get_args().to_numpy()
        fluxes.loc[failed] = retry.fluxes.to_numpy()
    return pd.DataFrame(args), pd.DataFrame(fluxes)
