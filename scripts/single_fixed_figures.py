"""Figs 9 and 10 of single-legacy and single-fixed against the batch growth data.

    uv run python scripts/v1_figures.py

Fig 9: pre-culture steady state at the density of data row 24, diluted to the first
data point, 80 µE, 13.5 mm, 96 h. Fig 10: steady-state growth rate vs. cell density at
80 µE, 20 mm (as in the notebook). Writes docs/reports/figures/single-fig9-fig10.png and
docs/reports/single-fixed-figures.json.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path

import matplotlib as mpl
import numpy as np
from mxlpy import Simulator
from mxlpy.integrators import Sundials

from synmep.data import load_legacy_inputs
from synmep.models import single_fixed, single_legacy
from synmep.runners import NATIVE

mpl.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
OUT = Path("docs/reports")


def steady_state(model, y0):
    for integ in (NATIVE, Sundials):
        res = (
            Simulator(model, y0=y0, integrator=integ)
            .simulate_to_steady_state(oscillation_detector=None)
            .get_result()
        )
        if hasattr(res.value, "variables"):
            return res.value
    return None


def fig9(module, inp) -> dict | None:
    data = inp.exp_growth_rate_by_density
    ms = module.get_model(inp, "dshc", "steady")
    y0 = ms.get_initial_conditions() | {"Cell_density": data.index[24] / 1e6}
    ss = steady_state(ms, y0)
    if ss is None:
        return None
    y0 = ss.variables.iloc[-1][ms.get_variable_names()].to_dict()
    md = module.get_model(inp, "dshc", "dynamic")
    pars = {"light": 80, "sample_depth_m": 0.0135, "Shc_ksynth": 0, "DXS_inh_str": 0}
    md.update_parameters(
        {k: v for k, v in pars.items() if k in md.get_parameter_names()}
    )
    t = np.linspace(0, 96 * 60, len(data))
    res = (
        Simulator(md, y0=y0 | {"Cell_density": data.index[0] / 1e6}, integrator=NATIVE)
        .simulate_time_course(t)
        .get_result()
    )
    if not hasattr(res.value, "variables"):
        return None
    a = res.value.get_args()
    mu, n = a["Growth_rate"].to_numpy() * 60, a["Cell_density"].to_numpy() * 1e6
    return {
        "t_h": (t / 60).tolist(),
        "mu": mu.tolist(),
        "N": n.tolist(),
        "rmse_mu_per_h": float(np.sqrt(np.mean((mu - data.to_numpy()) ** 2))),
        "rmse_ln_N": float(np.sqrt(np.mean((np.log(n) - np.log(data.index)) ** 2))),
    }


def fig10(module, inp) -> dict:
    data = inp.exp_growth_rate_by_density
    m = module.get_model(inp, "dshc", "steady")
    m.update_parameters({"light": 80, "sample_depth_m": 0.02})
    dens = data.index[::2]
    mus = []
    for d in dens:
        ss = steady_state(m, m.get_initial_conditions() | {"Cell_density": d / 1e6})
        mus.append(
            float(ss.get_args().iloc[-1]["Growth_rate"] * 60) if ss else float("nan")
        )
    mus = np.array(mus)
    ref = data.iloc[::2].to_numpy()
    ok = np.isfinite(mus)
    return {
        "density": list(dens),
        "mu": mus.tolist(),
        "n_converged": int(ok.sum()),
        "rmse_mu_per_h": float(np.sqrt(np.mean((mus[ok] - ref[ok]) ** 2))),
    }


def main() -> None:
    legacy, fixed_inp = load_legacy_inputs(), single_fixed.load_inputs()
    data = legacy.exp_growth_rate_by_density
    results = {
        name: {"fig9": fig9(mod, inp), "fig10": fig10(mod, inp)}
        for name, mod, inp in (
            ("single-legacy", single_legacy, legacy),
            ("single-fixed", single_fixed, fixed_inp),
        )
    }
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(13, 3.6), layout="constrained")
    t = np.linspace(0, 96, len(data))
    a1.plot(t, data.to_numpy(), "o", ms=3, color="k", label="data")
    a2.semilogy(t, data.index, "o", ms=3, color="k", label="data")
    a3.semilogx(data.index, data.to_numpy(), "o", ms=3, color="k", label="data")
    for name, r in results.items():
        if r["fig9"]:
            a1.plot(r["fig9"]["t_h"], r["fig9"]["mu"], label=name)
            a2.semilogy(r["fig9"]["t_h"], r["fig9"]["N"], label=name)
        a3.semilogx(r["fig10"]["density"], r["fig10"]["mu"], label=name)
    a1.set(xlabel="t [h]", ylabel=r"growth rate [h$^{-1}$]", title="Fig 9: batch")
    a2.set(xlabel="t [h]", ylabel=r"cells mL$^{-1}$", title="Fig 9: density")
    a3.set(
        xlabel=r"cells mL$^{-1}$",
        ylabel=r"steady-state growth rate [h$^{-1}$]",
        title="Fig 10",
    )
    for ax in (a1, a2, a3):
        ax.legend()
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "figures" / "single-fig9-fig10.png", dpi=110)
    summary = {
        name: {
            "fig9_rmse_mu_per_h": r["fig9"]["rmse_mu_per_h"] if r["fig9"] else None,
            "fig9_rmse_ln_N": r["fig9"]["rmse_ln_N"] if r["fig9"] else None,
            "fig10_rmse_mu_per_h": r["fig10"]["rmse_mu_per_h"],
            "fig10_n_converged": r["fig10"]["n_converged"],
        }
        for name, r in results.items()
    }
    (OUT / "single-fixed-figures.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
