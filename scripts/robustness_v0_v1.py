"""Robustness of v0 vs v1 under parameter perturbations (docs/reports/02-v1-robustness.md).

    uv run python scripts/robustness_v0_v1.py

Same perturbation groups as the ADR 0002 probe of the legacy model: each sample scales
every parameter of a group by an independent log-uniform factor in [0.1, 10] (the
"density" group scales the initial cell density). Steady state (strain dshc, steady
growth mode, density of data row 24) with the native runner (CVODE_BDF, no oscillation
detection), retried with Sundials on failure. The same random factors are used for both
models. Writes docs/reports/v1-robustness.json.
"""

from __future__ import annotations

import copy
import json
import time
import warnings
from pathlib import Path

import numpy as np
from mxlpy import Simulator
from mxlpy.integrators import Sundials

from synmep.data import load_legacy_inputs
from synmep.models import v0_faithful, v1_fixed
from synmep.runners import NATIVE

warnings.filterwarnings("ignore")
N = 24
SEED = 0
GROUPS = ("kcat", "ksynth", "out_metabolism", "light", "density")


def group_parameters(names: list[str], group: str) -> list[str]:
    suffix = {"kcat": "Kcat", "ksynth": "_ksynth", "out_metabolism": "_out_metabolism"}
    if group == "light":
        return ["light"]
    if group == "density":
        return []
    return [n for n in names if n.endswith(suffix[group])]


def steady_state(model, y0) -> tuple[str, float | None]:
    for label, integ in (("native", NATIVE), ("sundials", Sundials)):
        res = (
            Simulator(model, y0=y0, integrator=integ)
            .simulate_to_steady_state(oscillation_detector=None)
            .get_result()
        )
        if hasattr(res.value, "variables"):
            mu = float(res.value.get_args().iloc[-1]["Growth_rate"])
            if np.isfinite(mu) and mu > 0:
                return label, mu
    return "failed", None


def run(name: str, base, y0: dict[str, float]) -> dict:
    rng = np.random.default_rng(SEED)
    names = base.get_parameter_names()
    out = {}
    for group in GROUPS:
        pars = group_parameters(names, group)
        factors = 10 ** rng.uniform(-1, 1, (N, max(len(pars), 1)))
        counts = {"native": 0, "sundials": 0, "failed": 0}
        mus, t0 = [], time.perf_counter()
        for i in range(N):
            m = copy.deepcopy(base)
            yy = dict(y0)
            if group == "density":
                yy["Cell_density"] = y0["Cell_density"] * factors[i, 0]
            else:
                values = m.get_parameter_values()
                m.update_parameters(
                    {p: values[p] * f for p, f in zip(pars, factors[i], strict=True)}
                )
            label, mu = steady_state(m, yy)
            counts[label] += 1
            if mu is not None:
                mus.append(mu * 60)
        out[group] = {
            "n": N,
            "n_parameters": len(pars),
            "native_ok": counts["native"] / N,
            "any_ok": (counts["native"] + counts["sundials"]) / N,
            "wall_s_per_sample": (time.perf_counter() - t0) / N,
            "mu_per_h_range": [min(mus), max(mus)] if mus else None,
        }
        print(name, group, out[group], flush=True)
    return out


def main() -> None:
    data = load_legacy_inputs().exp_growth_rate_by_density
    results = {}
    for name, module, inp in (
        ("v0", v0_faithful, load_legacy_inputs()),
        ("v1", v1_fixed, v1_fixed.load_inputs()),
    ):
        base = module.get_model(inp, "dshc", "steady")
        y0 = base.get_initial_conditions() | {"Cell_density": data.index[24] / 1e6}
        results[name] = run(name, base, y0)
    Path("docs/reports/v1-robustness.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
