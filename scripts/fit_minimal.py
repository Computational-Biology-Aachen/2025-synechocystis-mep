"""Fit `minimal-m0` (ADR 0005) to the batch growth data.

    uv run python scripts/fit_minimal.py [--capacity]

Same protocol as `fit_g_growth.py`: 48 points of data/growth_rates_by_density.csv over
96 h; each simulation starts from the steady state at the pre-culture density (index 24
of the data), then dilutes to the first data point. The enzyme cost doesn't act in the
wild type (no overexpression), so the fit uses the cost-free model.
Fitted (log-space, relative to the defaults): petc_scale, pbs_scale, pigment_coupling
(exponent of the growth coupling of the pigment-branch enzymes) (+ capacity_K,
capacity_tau), within `BOUNDS`; least squares from the 3 best points of a coarse
start grid. Writes docs/reports/minimal-m0-fit[-cap].json.
"""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import Simulator
from scipy.optimize import least_squares

from synmep.data import load_legacy_inputs
from synmep.models import minimal
from synmep.runners import NATIVE

warnings.filterwarnings("ignore")
SIGMA_MU = 0.005  # [h^-1]
SIGMA_LN_N = 0.05
BATCH = {"light": 80.0, "sample_depth": 0.0135}  # MC-1000 batch (paper Fig. 9)
# plausible ranges: pbs_scale within the PETC training domain (ADR 0004), PETC
# efficiency within a factor ~5 of pfennig2024
BOUNDS = {"petc_scale": (0.2, 5.0), "pbs_scale": (0.25, 4.0),
          "pigment_coupling": (0.1, 10.0), "capacity_K": (0.01, 100.0),
          "capacity_tau": (30.0, 6000.0)}  # fmt: skip


def simulate(ms, md, pars: dict[str, float], n_pre: float, n0: float,
             t: np.ndarray) -> pd.DataFrame | None:  # fmt: skip
    ms.update_parameters(pars | {"N": n_pre})
    md.update_parameters(pars)
    # steady mode: slow states sped up 100x, so a long time course reaches the fixed point
    ss = Simulator(ms, integrator=NATIVE).simulate(1e5).get_result()
    if not hasattr(ss.value, "variables"):
        return None
    y0 = ss.value.variables.iloc[-1][ms.get_variable_names()].to_dict()
    rhs = ms.get_right_hand_side(y0).to_numpy()
    if np.max(np.abs(rhs) / (np.abs(np.array(list(y0.values()))) + 1e-9)) > 1e-5:
        return None
    y0 |= {"N": n0}
    tc = Simulator(md, y0=y0, integrator=NATIVE).simulate_time_course(t).get_result()
    return tc.value.get_args() if hasattr(tc.value, "variables") else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capacity", action="store_true")
    capacity = parser.parse_args().capacity
    fit_names = ["petc_scale", "pbs_scale", "pigment_coupling"]
    if capacity:
        fit_names += ["capacity_K", "capacity_tau"]
    inp = load_legacy_inputs()
    data = inp.exp_growth_rate_by_density
    t = np.linspace(0, 96 * 60, len(data))
    kw = dict(cost="none", capacity=capacity, use_fit=False, **BATCH)
    ms = minimal.get_model(inp, growth_mode="steady", **kw)
    md = minimal.get_model(inp, growth_mode="dynamic", **kw)
    p0 = np.array([md.get_parameter_values()[n] for n in fit_names])
    n0, n_pre = float(data.index[0]), float(data.index[24])

    def residuals(logp: np.ndarray) -> np.ndarray:
        pars = dict(zip(fit_names, p0 * np.exp(logp), strict=True))
        res = simulate(ms, md, pars, n_pre, n0, t)
        if res is None or len(res) != len(t):
            return np.full(2 * len(t), 1e3)
        r_mu = (res["mu"].to_numpy() * 60 - data.to_numpy()) / SIGMA_MU
        r_n = (np.log(res["N"].to_numpy()) - np.log(data.index.to_numpy())) / SIGMA_LN_N
        return np.concatenate([r_mu, r_n])

    t0 = time.perf_counter()
    cost0 = 0.5 * np.sum(residuals(np.zeros(len(fit_names))) ** 2)
    # coarse grid over petc_scale x pbs_scale as start (the defaults can sit where the
    # finite-difference Jacobian vanishes)
    grid = [
        np.array([np.log(a), np.log(b), np.log(c), *np.zeros(len(fit_names) - 3)])
        for a in (0.5, 1.0, 2.0, 4.0)
        for b in (0.5, 1.0, 2.0, 4.0)
        for c in (0.3, 1.0, 3.0)
    ]
    lo = np.log([BOUNDS[n][0] for n in fit_names] / p0)
    hi = np.log([BOUNDS[n][1] for n in fit_names] / p0)
    grid = [np.clip(x, lo, hi) for x in grid]
    starts = sorted(grid, key=lambda x: float(np.sum(residuals(x) ** 2)))[:3]
    fit = min(
        (
            least_squares(residuals, x0, bounds=(lo, hi), diff_step=1e-3, max_nfev=200)
            for x0 in starts
        ),
        key=lambda f: f.cost,
    )
    wall = time.perf_counter() - t0
    best = dict(zip(fit_names, p0 * np.exp(fit.x), strict=True))
    res = simulate(ms, md, best, n_pre, n0, t)
    assert res is not None
    mu_sim = res["mu"].to_numpy() * 60
    report = {
        "capacity": capacity,
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
    out = Path(f"docs/reports/minimal-m0-fit{'-cap' if capacity else ''}.json")
    out.write_text(json.dumps(report, indent=1))
    print(json.dumps({k: v for k, v in report.items() if k != "trajectory"}, indent=1))


if __name__ == "__main__":
    main()
