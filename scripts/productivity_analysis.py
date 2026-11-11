"""Maximal squalene productivity of a culture: the squalene/pigment trade-off.

    uv run python scripts/productivity_analysis.py [--strains] [variant ...]

With --strains, the strains of the MCA overexpression sequence
(docs/reports/strain-optimisation.json) are evaluated instead of the squalene synthase
scan → docs/reports/productivity-strains.json.

Question: what is the maximal production rate of squalene per culture volume, given
that diverting FPP into squalene takes precursors away from the pigments needed for
growth?

Setting: continuous culture (turbidostat) at steady state, Δshc strain, 100 µmol m^-2
s^-1, 1 cm light path. Squalene leaves the cells only by growth dilution, so at steady
state the volumetric productivity is

    P [mg L^-1 h^-1] = v_sq [mM min^-1] · 60 · N [cells mL^-1] · 1e3 · V_cell [L] · 410.73

Scanned: the squalene synthase level (the "pull" on FPP, 0.1-1e4x wild type) × the
cell density the turbidostat holds. Reported per variant: the trade-off between growth
and specific squalene production, pigment content, and the optimum of P.

Thread g additionally gets a stoichiometric bound: all 3PGA supplied at that density
going to squalene (12 3PGA per squalene via 6 DXS reactions).

Variants: v0, v1 (empirical growth law μ ∝ v_Chl), g-none and g-coupled (pigments
follow synthesis, so the trade-off acts through light absorption and self-shading),
g-photoacclimation and g-capacity (pigment contents fixed by construction: controls
without the pigment side of the trade-off). Writes docs/reports/productivity.json.
"""

from __future__ import annotations

import copy
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

from synmep.data import load_legacy_inputs
from synmep.runners import NATIVE_RUNNER, SUNDIALS_RUNNER

warnings.filterwarnings("ignore")
OUT = Path("docs/reports/productivity.json")
OUT_STRAINS = Path("docs/reports/productivity-strains.json")
STRAINS = Path("docs/reports/strain-optimisation.json")
M_SQ = 410.73
SQS_FACTORS = np.geomspace(0.1, 1e4, 16)
DENSITIES = np.geomspace(3e6, 3e10, 17)  # [cells mL^-1] (3e9 ~ OD 30, practical limit)
PGA_PER_SQ = 12  # 6 C5 units, each from one DXS reaction (G3P + PYR = 2 3PGA)


def legacy_adapter(name: str):
    if name == "v0":
        from synmep.models import v0_faithful as mod

        inp, runner, vcell = load_legacy_inputs(), NATIVE_RUNNER, 2e-15
    else:
        from synmep.data_v1 import load_inputs
        from synmep.models import v1_fixed as mod

        inp, runner, vcell = load_inputs(), SUNDIALS_RUNNER, 4e-15
    base = mod.get_model(inp, "dshc", "steady")
    base.update_parameters({"light": 100.0, "sample_depth_m": 0.01})
    p0 = base.get_parameter_values()

    def configure(folds: dict[str, float], n: float):
        m = copy.deepcopy(base)
        m.update_parameters({k: p0[k] * f for k, f in folds.items()})
        m.update_variable("Cell_density", n / 1e6)  # [cells nL^-1]
        return m

    keys = {"v_sq": "v_Sqs", "mu": "Growth_rate", "chl": "chlorophyll_phytol",
            "sq": "SQ"}  # fmt: skip
    return configure, runner, vcell, keys, None, "Sqs_ksynth"


def g_adapter(variant: str):
    from synmep.models import g_coupled

    base = g_coupled.get_model(load_legacy_inputs(), variant)  # type: ignore[arg-type]

    p0 = base.get_parameter_values()

    def configure(folds: dict[str, float], n: float):
        m = copy.deepcopy(base)
        m.update_parameters({k: p0[k] * f for k, f in folds.items()} | {"N": n})
        return m

    keys = {"v_sq": "v_sq", "mu": "mu", "chl": "Chl", "sq": "SQ"}
    return configure, NATIVE_RUNNER, 4e-15, keys, "S_PGA", "E_Sqs"


def v3_adapter(name: str):
    from synmep.models import v3_minimal

    base = v3_minimal.get_model(load_legacy_inputs(), **v3_minimal.variant(name))
    p0 = base.get_parameter_values()

    def configure(folds: dict[str, float], n: float):
        m = copy.deepcopy(base)
        m.update_parameters({k: p0[k] * f for k, f in folds.items()} | {"N": n})
        return m

    keys = {"v_sq": "v_sq", "mu": "mu", "chl": "Chl", "sq": "SQ"}
    return configure, NATIVE_RUNNER, 4e-15, keys, "S_PGA", "E_Sqs"


def _get(fluxes, args, key: str) -> float:
    return float(fluxes[key] if key in fluxes else args[key])


def adapter(name: str):
    if name in ("v0", "v1"):
        return legacy_adapter(name)
    if name.startswith("v3"):
        return v3_adapter(name)
    return g_adapter(name.removeprefix("g-"))


def evaluate(configure, runner, vcell, keys, supply, settings) -> list[dict]:
    """Steady states for (label, folds, density) settings, with continuation."""
    out, y0 = [], None
    for label, folds, n in settings:
        m = configure(folds, float(n))
        if y0 is not None:  # continuation from the previous steady state
            m.update_variables({k: v for k, v in y0.items() if k != "Cell_density"})
        res = runner.steady_state(m)
        row = label | {"density": float(n)}
        if hasattr(res.value, "fluxes"):
            a = res.value.get_args().iloc[-1]
            fl = res.value.fluxes.iloc[-1]
            v_sq, mu = _get(fl, a, keys["v_sq"]), _get(fl, a, keys["mu"])
            culture_volume_fraction = n * 1e3 * vcell  # L cell per L culture
            row |= {
                "mu_per_h": mu * 60,
                "v_sq": v_sq,
                "chl": _get(fl, a, keys["chl"]),
                "sq_content": _get(fl, a, keys["sq"]),
                "productivity_mg_L_h": v_sq * 60 * culture_volume_fraction * M_SQ,
            }
            if supply is not None:
                row["bound_mg_L_h"] = (
                    _get(fl, a, supply)
                    / PGA_PER_SQ
                    * 60
                    * culture_volume_fraction
                    * M_SQ
                )
            y0 = a[m.get_variable_names()].to_dict()
        out.append(row)
    return out


def run(name: str) -> dict:
    """Squalene synthase level x density."""
    configure, runner, vcell, keys, supply, sqs = adapter(name)
    settings = [
        ({"sqs_factor": float(f)}, {sqs: float(f)}, n)
        for f in SQS_FACTORS
        for n in DENSITIES
    ]
    return {"grid": evaluate(configure, runner, vcell, keys, supply, settings),
            "cell_volume": vcell}  # fmt: skip


def run_strains(name: str) -> dict:
    """Strains from the MCA overexpression sequence (strain-optimisation.json) x density.

    Round k = the first k overexpressions of that variant's sequence, each with its
    fold change; round 0 = wild type.
    """
    configure, runner, vcell, keys, supply, _ = adapter(name)
    seq = json.loads(STRAINS.read_text())[name]
    settings = []
    for r in seq["rounds"]:
        if r.get("failed"):
            continue
        folds = {p: seq["folds"][p] for p in r["overexpressed"]}
        label = {"round": r["round"], "overexpressed": list(r["overexpressed"])}
        settings += [(label, folds, n) for n in DENSITIES]
    return {"grid": evaluate(configure, runner, vcell, keys, supply, settings),
            "cell_volume": vcell}  # fmt: skip


def main() -> None:
    strains = "--strains" in sys.argv
    args = [a for a in sys.argv[1:] if a != "--strains"]
    out = OUT_STRAINS if strains else OUT
    names = args or ["v0", "v1", "g-none", "g-coupled", "g-photoacclimation",
                     "g-capacity", "v3-drain", "v3-proteome", "v3-drain-cap",
                     "v3-proteome-cap"]  # fmt: skip
    results = json.loads(out.read_text()) if out.exists() else {}
    for name in names:
        t = time.perf_counter()
        print(name, flush=True)
        res = run_strains(name) if strains else run(name)
        results[name] = res | {"wall_s": time.perf_counter() - t}
        out.write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
