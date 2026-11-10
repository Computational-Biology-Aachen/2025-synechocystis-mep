"""Generate reference outputs of the legacy (modelbase) model.

Runs ONLY in the old pixi environment:

    pixi run python scripts/legacy/make_reference.py [--quick]

Writes to data/reference/:
- surrogate_poly.npz            exact coefficients of the PETC surrogate polynomial
- surrogate_check.csv           surrogate inputs + sklearn predictions (exactness test)
- parameters.json / y0.json     legacy parameters and initial conditions (dshc, wt)
- rhs_points.csv                RHS + fluxes at several states (strongest port check)
- y0_ss.json                    steady state used as initial condition for the figures
- fig09_timecourse.csv          Fig 9 batch growth trajectory
- fig10_growth_by_density.csv   Fig 10 steady-state growth-rate scan
- fig11_induction.csv           Fig 11 at a subset of induction time points
- timings.json                  legacy wall times (benchmark baseline)
"""

from __future__ import annotations

import json
import sys
import time
import warnings
from functools import partial
from pathlib import Path

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
from legacy_setup import (  # type: ignore
    BIOMASS_REGRESSION_MODEL,
    DAY_IN_MIN,
    EXP_GROWTH_RATE_BY_DENSITY,
    get_parameterised_model,
    unwrap,
)
from modelbase.ode import Simulator

from parallel import parallelise

OUT = Path("data/reference")
OUT.mkdir(parents=True, exist_ok=True)
QUICK = "--quick" in sys.argv
TIMINGS: dict[str, float] = {}


def _timed(name: str):
    def deco(fn):
        def inner(*a, **k):
            t = time.perf_counter()
            r = fn(*a, **k)
            TIMINGS[name] = time.perf_counter() - t
            print(f"{name}: {TIMINGS[name]:.1f} s", flush=True)
            return r

        return inner

    return deco


def dump_surrogate() -> None:
    pipe = BIOMASS_REGRESSION_MODEL
    poly = pipe.named_steps["polynomialfeatures"]
    reg = pipe.named_steps["multioutputregressor"]
    coef = np.vstack([e.coef_ for e in reg.estimators_])
    intercept = np.array([e.intercept_ for e in reg.estimators_])
    np.savez(
        OUT / "surrogate_poly.npz",
        powers=poly.powers_,
        coef=coef,
        intercept=intercept,
        feature_names=np.array(pipe.feature_names_in_).astype(str),
        output_names=np.array(["ATP", "NADPH", "3PGA", "Fd_red"]).astype(str),
    )
    # Exactness check samples: spread around a typical operating point
    rng = np.random.default_rng(0)
    base = {
        "complex_abs_ps1": 5e3, "complex_abs_ps2": 5e3, "complex_abs_pbs": 5e3,
        "light_ocp": 10.0, "pigment_phycocyanin": 1.5,
        "pigment_allophycocyanin": 0.7, "pigment_beta_carotene": 0.03,
    }  # fmt: skip
    names = list(pipe.feature_names_in_)
    X = pd.DataFrame(
        {n: base.get(n, 1.0) * 10 ** rng.uniform(-1, 1, 200) for n in names}
    )
    Y = pipe.predict(X)
    df = X.copy()
    for i, o in enumerate(["ATP", "NADPH", "3PGA", "Fd_red"]):
        df[f"out_{o}"] = Y[:, i]
    df.to_csv(OUT / "surrogate_check.csv", index=False)


def dump_structure() -> None:
    for strain in ("dshc", "wt"):
        for mode in ("steady", "dynamic"):
            m, y0 = get_parameterised_model(strain=strain, growth_mode=mode)  # type: ignore
            tag = f"{strain}_{mode}"
            (OUT / f"parameters_{tag}.json").write_text(
                json.dumps(
                    {k: float(v) for k, v in m.get_parameters().items()}, indent=1
                )
            )
            (OUT / f"y0_{tag}.json").write_text(
                json.dumps({k: float(v) for k, v in y0.items()}, indent=1)
            )


def _rhs_rows(m, states: dict[str, dict[str, float]]) -> pd.DataFrame:
    rows = {}
    names = m.get_compounds()
    for label, y in states.items():
        yv = np.array([y[n] for n in names])
        rhs = dict(zip(names, np.asarray(m._get_rhs(0, yv)), strict=True))
        fl = m.get_fluxes_dict(y, 0)
        full = m.get_full_concentration_dict(y, 0)
        row = {f"y:{k}": v for k, v in y.items()}
        row |= {f"dydt:{k}": float(v) for k, v in rhs.items()}
        row |= {f"flux:{k}": float(np.atleast_1d(v)[0]) for k, v in fl.items()}
        row |= {
            f"full:{k}": float(np.atleast_1d(v)[0])
            for k, v in full.items()
            if k not in y
        }
        rows[label] = row
    return pd.DataFrame(rows).T


@_timed("steady_state_y0")
def steady_state_y0() -> dict[str, float]:
    ms, y0 = get_parameterised_model(strain="dshc", growth_mode="steady")
    return unwrap(
        Simulator(ms)
        .initialise(y0 | {"Cell_density": EXP_GROWTH_RATE_BY_DENSITY.index[24] / 1e6})
        .simulate_to_steady_state_and()
        .get_new_y0()
    )


def dump_rhs_points(y0_ss: dict[str, float]) -> None:
    rng = np.random.default_rng(1)
    for strain, mode in (("dshc", "steady"), ("dshc", "dynamic"), ("wt", "dynamic")):
        m, y0 = get_parameterised_model(strain=strain, growth_mode=mode)  # type: ignore
        y0 = y0 | {"Cell_density": EXP_GROWTH_RATE_BY_DENSITY.index[24] / 1e6}
        states = {"y0": y0, "ss": y0_ss}
        for i in range(3):
            states[f"rand{i}"] = {
                k: v * 10 ** rng.uniform(-0.5, 0.5) for k, v in y0_ss.items()
            }
        _rhs_rows(m, states).to_csv(OUT / f"rhs_points_{strain}_{mode}.csv")


@_timed("fig09")
def fig09(y0_ss: dict[str, float]) -> pd.DataFrame:
    md, _ = get_parameterised_model(strain="dshc", growth_mode="dynamic")
    md.update_parameters(
        {"light": 80, "sample_depth_m": 0.0135, "ShcKcat": 0, "Shc_ksynth": 0, "DXS_inh_str": 0}
    )  # fmt: skip
    s = Simulator(md).initialise(y0_ss | {"Cell_density": 54562873.0217531 / 1e6})
    s.simulate(96 * 60, steps=96 * 4)
    return unwrap(s.get_full_results_df())


def _fig10_worker(density: float, model, y0) -> float:
    try:
        df = (
            Simulator(model)
            .initialise(y0 | {"Cell_density": density / 1e6})
            .simulate_to_steady_state_and()
            .get_full_results_df()
        )
        if df is not None:
            return float(df.iloc[-1]["Growth_rate"] * 60)
    except Exception:  # noqa: BLE001
        pass
    return float("nan")


@_timed("fig10")
def fig10() -> pd.Series:
    ms, y0 = get_parameterised_model(strain="dshc", growth_mode="steady")
    ms.update_parameters({"light": 80, "sample_depth_m": 0.02})
    dens = EXP_GROWTH_RATE_BY_DENSITY.index[::2]
    res = parallelise(
        partial(_fig10_worker, model=ms, y0=y0),
        inputs=list(zip(dens, dens, strict=True)),
        disable_tqdm=True,
    )
    return pd.Series(res, dtype=float, name="growth_rate_per_h").sort_index()


def _fig11_worker(t_ind: float, y0_ss, t_end: float) -> pd.Series:
    md, _ = get_parameterised_model(strain="dshc", growth_mode="dynamic")
    s = Simulator(md).initialise(y0_ss)
    if t_ind > 0:
        s.simulate(t_ind)
    s.update_parameter("CrtE3Kcat", 0.01)
    s.simulate(t_end)
    c = s.get_full_results_df()
    return pd.Series({}) if c is None else c.iloc[-1]


@_timed("fig11")
def fig11(y0_ss) -> pd.DataFrame:
    t_end = 14 * DAY_IN_MIN
    pts = np.arange(0, t_end, DAY_IN_MIN * (2 if QUICK else 1), dtype=float)
    res = parallelise(
        partial(_fig11_worker, y0_ss=y0_ss, t_end=t_end),
        list(zip(pts, pts, strict=True)),
        disable_tqdm=True,
    )
    return pd.DataFrame(res).T.sort_index()


def main() -> None:
    dump_surrogate()
    dump_structure()
    y0_ss = steady_state_y0()
    (OUT / "y0_ss.json").write_text(
        json.dumps({k: float(v) for k, v in y0_ss.items()}, indent=1)
    )
    dump_rhs_points(y0_ss)
    fig09(y0_ss).to_csv(OUT / "fig09_timecourse.csv")
    fig10().to_csv(OUT / "fig10_growth_by_density.csv")
    fig11(y0_ss).to_csv(OUT / "fig11_induction.csv")
    (OUT / "timings.json").write_text(json.dumps(TIMINGS, indent=1))


if __name__ == "__main__":
    main()
