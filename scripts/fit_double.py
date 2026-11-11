"""Fit the `double` growth model to the batch growth data (ADR 0004).

    uv run python scripts/fit_double.py [--pigment-regulation none|coupled]

Data: data/growth_rates_by_density.csv, 48 evenly spaced points over 96 h (checked:
the integral of the reported growth rates matches ln(N_end/N_0) within 3 %). Each
simulation starts from the steady state at the pre-culture density, taken as the
density at index 24 of the data (as in the legacy notebook), then dilutes to the first
data point.
Fitted (log-space): petc_scale, pbs_scale. Writes docs/reports/double-fit[-coupled].json.
"""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import Simulator
from scipy.optimize import least_squares

from synmep.data import load_legacy_inputs
from synmep.models.double_growth import get_model
from synmep.runners import NATIVE

warnings.filterwarnings("ignore")
FIT_BY_VARIANT = {
    "none": ("petc_scale", "pbs_scale"),
    "coupled": ("petc_scale", "pbs_scale"),
}
SIGMA_MU = 0.005  # [h^-1]
SIGMA_LN_N = 0.05


def simulate(
    ms, md, pars: dict[str, float], n_pre: float, n0: float, t: np.ndarray
) -> pd.DataFrame | None:
    for m in (ms, md):
        m.update_parameters(pars)
    y0 = ms.get_initial_conditions() | {"N": n_pre}
    # steady mode relaxes every state at k_relax = 1 min^-1, so one long time course
    # reaches the fixed point; Native's steady-state search needed its full 1000
    # iterations here (~16 s per evaluation)
    ss = Simulator(ms, y0=y0, integrator=NATIVE).simulate(1e4).get_result()
    if not hasattr(ss.value, "variables"):
        return None
    y0 = ss.value.variables.iloc[-1][ms.get_variable_names()].to_dict()
    rhs = ms.get_right_hand_side(y0).to_numpy()
    scale = np.abs(np.array(list(y0.values()))) + 1e-9
    if np.max(np.abs(rhs) / scale) > 1e-6:
        return None
    y0 |= {"N": n0}
    tc = Simulator(md, y0=y0, integrator=NATIVE).simulate_time_course(t).get_result()
    return tc.value.get_args() if hasattr(tc.value, "variables") else None


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pigment-regulation",
        default="none",
        choices=list(FIT_BY_VARIANT),
    )
    variant = parser.parse_args().pigment_regulation
    fit_names = FIT_BY_VARIANT[variant]
    suffix = "" if variant == "none" else f"-{variant}"
    inp = load_legacy_inputs()
    data = inp.exp_growth_rate_by_density
    t = np.linspace(0, 96 * 60, len(data))
    ms = get_model(inp, "steady", pigment_regulation=variant)
    md = get_model(inp, "dynamic", pigment_regulation=variant)
    p0 = np.array([md.get_parameter_values()[n] for n in fit_names])
    n0, n_pre = float(data.index[0]), float(data.index[24])

    def residuals(logp: np.ndarray) -> np.ndarray:
        res = simulate(
            ms, md, dict(zip(fit_names, p0 * np.exp(logp), strict=True)), n_pre, n0, t
        )
        if res is None or len(res) != len(t):
            return np.full(2 * len(t), 1e3)
        r_mu = (res["mu"].to_numpy() * 60 - data.to_numpy()) / SIGMA_MU
        r_n = (np.log(res["N"].to_numpy()) - np.log(data.index.to_numpy())) / SIGMA_LN_N
        return np.concatenate([r_mu, r_n])

    x0 = np.zeros(len(fit_names))  # log of the factor relative to the defaults
    t0 = time.perf_counter()
    cost0 = 0.5 * np.sum(residuals(x0) ** 2)
    fit = least_squares(residuals, x0, x_scale=1.0, diff_step=1e-3, max_nfev=200)
    wall = time.perf_counter() - t0
    best = dict(zip(fit_names, p0 * np.exp(fit.x), strict=True))
    res = simulate(ms, md, best, n_pre, n0, t)
    assert res is not None
    mu_sim = res["mu"].to_numpy() * 60
    report = {
        "variant": variant,
        "fitted": best,
        "cost_initial": float(cost0),
        "cost_final": float(fit.cost),
        "nfev": int(fit.nfev),
        "wall_s": wall,
        "rmse_mu_per_h": float(np.sqrt(np.mean((mu_sim - data.to_numpy()) ** 2))),
        "rmse_ln_N": float(
            np.sqrt(np.mean((np.log(res["N"].to_numpy()) - np.log(data.index)) ** 2))
        ),
        "trajectory": {
            "t_h": (t / 60).tolist(),
            "mu_sim": mu_sim.tolist(),
            "mu_data": data.tolist(),
            "N_sim": res["N"].tolist(),
            "N_data": data.index.tolist(),
            "Chl": res["Chl"].tolist(),
            "Car": res["Car"].tolist(),
            "SQ": res["SQ"].tolist(),
        },
    }
    Path(f"docs/reports/double-fit{suffix}.json").write_text(
        json.dumps(report, indent=1)
    )
    print(json.dumps({k: v for k, v in report.items() if k != "trajectory"}, indent=1))


if __name__ == "__main__":
    main()
