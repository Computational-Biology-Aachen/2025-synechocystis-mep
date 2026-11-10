"""Benchmark v0: mxlpy default integrator (Sundials CVODE) vs native runner.

    uv run python scripts/benchmark_v0.py

Workloads mirror the legacy reference runs (scripts/legacy/make_reference.py), so the
legacy modelbase/assimulo wall times in data/reference/timings.json are a third
column. Accuracy is checked against the legacy reference outputs. Results are
written to docs/reports/benchmark-v0.json.
"""

from __future__ import annotations

import json
import os
import sys
import time
import warnings
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import Simulator, scan
from mxlpy.integrators import Native, Sundials
from mxlpy.parallel import parallelise
from mxlpy.runners import KineticRunner

from synmep.data import DAY_IN_MIN, load_legacy_inputs
from synmep.models.v0_faithful import get_model, sanitise

warnings.filterwarnings("ignore")
REF = Path("data/reference")
OUT = Path("docs/reports/benchmark-v0.json")
REPEATS = 3

NATIVE = partial(Native, method="CVODE_BDF")
INTEGRATORS = {"sundials": Sundials, "native": NATIVE}


def _ref_y0_ss() -> dict[str, float]:
    return {
        sanitise(k): v for k, v in json.loads((REF / "y0_ss.json").read_text()).items()
    }


def _timed(fn, *args, **kwargs):
    t = time.perf_counter()
    res = fn(*args, **kwargs)
    return time.perf_counter() - t, res


def steady_state(model, y0, integrator) -> pd.Series:
    res = Simulator(model, y0=y0, integrator=integrator).simulate_to_steady_state()
    return res.get_result().unwrap_or_err().variables.iloc[-1]


def fig09(model, y0, integrator) -> pd.DataFrame:
    s = Simulator(model, y0=y0, integrator=integrator).simulate(96 * 60, steps=96 * 4)
    return s.get_result().unwrap_or_err().get_combined()


def fig11_point(t_ind: float, model, y0, integrator, t_end: float) -> pd.Series:
    s = Simulator(model, y0=y0, integrator=integrator)
    if t_ind > 0:
        s.simulate(t_ind)
    s.update_parameter("CrtE3Kcat", 0.01)
    s.simulate(t_end)
    return s.get_result().unwrap_or_err().get_combined().iloc[-1]


def main() -> None:
    inp = load_legacy_inputs()
    y0_ss = _ref_y0_ss()
    dens0 = inp.exp_growth_rate_by_density.index[24] / 1e6
    results: dict[str, dict] = {"cpu_count": os.cpu_count(), "workloads": {}}
    W = results["workloads"]

    # --- 0. model construction + native compile (cold) ------------------------------
    t_build, ms = _timed(get_model, inp, "dshc", "steady")
    y0 = ms.get_initial_conditions() | {"Cell_density": dens0}
    t_cold, _ = _timed(steady_state, ms, y0, NATIVE)
    W["build"] = {"model_construction_s": t_build, "native_first_call_s": t_cold}
    print(W["build"], flush=True)

    # --- 1. single steady state -----------------------------------------------------
    W["steady_state"] = {}
    for name, integ in INTEGRATORS.items():
        times, last = [], None
        for _ in range(REPEATS):
            dt, last = _timed(steady_state, ms, y0, integ)
            times.append(dt)
        dev = max(abs(last[k] - v) / abs(v) for k, v in y0_ss.items() if abs(v) > 1e-12)
        W["steady_state"][name] = {
            "median_s": float(np.median(times)),
            "max_rel_dev": dev,
        }
    print(W["steady_state"], flush=True)

    # --- 2. Fig 9 batch growth time course ------------------------------------------
    md = get_model(inp, "dshc", "dynamic")
    fig9_pars = {"light": 80, "sample_depth_m": 0.0135, "ShcKcat": 0, "Shc_ksynth": 0}
    fig9_pars |= {"DXS_inh_str": 0}
    # ShcKcat doesn't exist in dshc (unused parameter, dropped in v0)
    md.update_parameters(
        {k: v for k, v in fig9_pars.items() if k in md.get_parameter_names()}
    )
    y0_9 = y0_ss | {"Cell_density": 54562873.0217531 / 1e6}
    ref9 = pd.read_csv(REF / "fig09_timecourse.csv", index_col=0)
    _timed(fig09, md, y0_9, NATIVE)  # compile the dynamic variant outside the timing
    W["fig09"] = {}
    for name, integ in INTEGRATORS.items():
        times, res = [], None
        for _ in range(REPEATS):
            dt, res = _timed(fig09, md, y0_9, integ)
            times.append(dt)
        dev = {
            col: float(
                np.max(np.abs(res[col].to_numpy() - ref9[col].to_numpy()))
                / np.max(np.abs(ref9[col].to_numpy()))
            )
            for col in ("Growth_rate", "Cell_density", "SQ")
        }
        W["fig09"][name] = {
            "median_s": float(np.median(times)),
            "max_dev_rel_to_max": dev,
        }
    print(W["fig09"], flush=True)

    # --- 3. Fig 10 steady-state scan over cell density ------------------------------
    m10 = get_model(inp, "dshc", "steady")
    m10.update_parameters({"light": 80, "sample_depth_m": 0.02})
    dens = inp.exp_growth_rate_by_density.index[::2]
    to_scan = pd.DataFrame({"Cell_density": np.asarray(dens) / 1e6})
    ref10 = pd.read_csv(REF / "fig10_growth_by_density.csv", index_col=0).squeeze()
    W["fig10"] = {}
    for name, integ, parallel in (
        ("sundials_parallel", Sundials, True),
        ("native_serial", NATIVE, False),
        ("native_parallel", NATIVE, True),
    ):
        dt, res = _timed(
            scan.steady_state,
            m10,
            to_scan=to_scan,
            y0=m10.get_initial_conditions(),
            parallel=parallel,
            runner=KineticRunner(integrator=integ),
        )
        mu = res.get_args()["Growth_rate"].to_numpy() * 60
        dev = float(np.nanmax(np.abs(mu - ref10.to_numpy()) / np.abs(ref10.to_numpy())))
        W["fig10"][name] = {"wall_s": dt, "max_rel_dev": dev, "n": len(to_scan)}
    print(W["fig10"], flush=True)

    # --- 4. Fig 11 induction-time protocol ------------------------------------------
    t_end = 14 * DAY_IN_MIN
    pts = np.arange(0, t_end, DAY_IN_MIN, dtype=float)
    ref11 = pd.read_csv(REF / "fig11_induction.csv", index_col=0)
    W["fig11"] = {}
    for name, integ, parallel in (
        ("sundials_parallel", Sundials, True),
        ("native_serial", NATIVE, False),
        ("native_parallel", NATIVE, True),
    ):
        worker = partial(
            fig11_point, model=md_base(inp), y0=y0_ss, integrator=integ, t_end=t_end
        )
        dt, res = _timed(
            parallelise,
            worker,
            list(zip(pts, pts, strict=True)),
            parallel=parallel,
            disable_tqdm=True,
        )
        df = pd.DataFrame(dict(res)).T.sort_index()
        dev = {
            col: float(
                np.max(
                    np.abs(df[col].to_numpy() - ref11[col].to_numpy())
                    / np.abs(ref11[col].to_numpy())
                )
            )
            for col in ("Cell_density", "SQ")
        }
        W["fig11"][name] = {"wall_s": dt, "max_rel_dev": dev, "n": len(pts)}
    print(W["fig11"], flush=True)

    W["legacy_modelbase_assimulo"] = json.loads((REF / "timings.json").read_text())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(results, indent=1))


def md_base(inp):
    return get_model(inp, "dshc", "dynamic")


if __name__ == "__main__":
    sys.exit(main())
