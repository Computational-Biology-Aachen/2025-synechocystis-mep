"""MCA-based strain optimisation (notebook / paper Fig 12) for all model variants.

    uv run python scripts/strain_optimisation.py [variant ...]

Each round: flux response coefficients (normalised) of squalene synthesis and biomass
synthesis to every enzyme parameter at steady growth; the enzyme with the largest sum
is overexpressed; **each enzyme at most once**. The first round is the answer to "which
single overexpression is best".

Overexpression is a fold change of the enzyme (v0/v1: its synthesis rate, thread g:
its level): fold = 1 + 1.331057 kDa / (2 E0) with E0 the legacy (2 fL) proteomics
concentration [µM]. For v0 this equals the notebook's additive
`ksynth + 1.331057 kDa / (5.5 h)` exactly. Thread g's lumped MEP chain gets the
geometric mean of the chain enzymes' folds (caveat: that's several genes).

Conditions as in the notebook: steady growth (fixed cell density, 1e7 cells mL^-1),
light 100 µmol m^-2 s^-1, sample depth 1 cm. Writes docs/reports/strain-optimisation.json.
"""

from __future__ import annotations

import copy
import json
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import Model, mca

from synmep.data import load_legacy_inputs
from synmep.runners import NATIVE_RUNNER, SUNDIALS_RUNNER

warnings.filterwarnings("ignore")
OUT = Path("docs/reports/strain-optimisation.json")
ROUNDS = 10
M_SQUALENE = 410.73  # [g mol^-1]

LEGACY_ENZYMES = ("DXS", "ThiG", "DXR", "IspD", "IspE", "IspF", "IspG", "IspH",
                  "Idi", "CrtE", "Sqs", "Shc", "ChlP", "CrtB", "ChlG")  # fmt: skip
G_ENZYMES = {  # thread-g lump -> legacy enzymes it stands for
    "DXS": ("DXS",),
    "MEP": ("DXR", "IspD", "IspE", "IspF", "IspG", "IspH"),
    "CrtE": ("CrtE",),
    "Sqs": ("Sqs",),
    "CrtB": ("CrtB",),
    "Chl": ("ChlP", "ChlG"),
}


def legacy_folds() -> dict[str, float]:
    inp = load_legacy_inputs()
    poi = inp.proteins_of_interest
    folds = {}
    for e in LEGACY_ENZYMES:
        kda = poi[poi["Enzyme"] == e]["kDa"]
        if len(kda) and e in inp.enzyme_concentrations_umol:
            e0 = inp.enzyme_concentrations_umol[e]
            folds[e] = 1 + 1.331057 * float(kda.iloc[0]) / (2 * e0)
    return folds


@dataclass
class Variant:
    name: str
    model: Model
    parameters: dict[str, float]  # parameter -> overexpression fold
    flux_sq: str
    flux_bio: str
    runner: object
    sq_content: callable  # steady-state result row -> intracellular SQ [mM]


def variants(selected: list[str]) -> list[Variant]:
    folds = legacy_folds()
    out = []
    if "v0" in selected:
        from synmep.models import v0_faithful

        m = v0_faithful.get_model(load_legacy_inputs(), "dshc", "steady")
        pars = {f"{e}_ksynth": f for e, f in folds.items()
                if f"{e}_ksynth" in m.get_parameter_names()}  # fmt: skip
        out.append(Variant("v0", m, pars, "v_Sqs", "v_Biomass", NATIVE_RUNNER,
                           lambda r: r["v_Sqs"] / r["Growth_rate"]))  # fmt: skip
    if "v1" in selected:
        from synmep.data_v1 import load_inputs
        from synmep.models import v1_fixed

        m = v1_fixed.get_model(load_inputs(), "dshc", "steady")
        pars = {f"{e}_ksynth": f for e, f in folds.items()
                if f"{e}_ksynth" in m.get_parameter_names()}  # fmt: skip
        out.append(Variant("v1", m, pars, "v_Sqs", "v_Biomass", SUNDIALS_RUNNER,
                           lambda r: r["v_Sqs"] / r["Growth_rate"]))  # fmt: skip
    from synmep.models import g_coupled

    for v in ("none", "coupled", "photoacclimation", "capacity"):
        if f"g-{v}" not in selected:
            continue
        m = g_coupled.get_model(load_legacy_inputs(), v)  # type: ignore[arg-type]
        pars = g_folds(folds)
        out.append(Variant(f"g-{v}", m, pars, "v_sq", "v_bio", NATIVE_RUNNER,
                           lambda r: r["v_sq"] / r["mu"]))  # fmt: skip
    from synmep.models import v3_minimal

    for name in v3_minimal.VARIANTS:
        if name not in selected:
            continue
        m = v3_minimal.get_model(load_legacy_inputs(), **v3_minimal.variant(name))
        out.append(Variant(name, m, g_folds(folds), "v_sq", "v_bio", NATIVE_RUNNER,
                           lambda r: r["v_sq"] / r["mu"]))  # fmt: skip
    return out


def g_folds(folds: dict[str, float]) -> dict[str, float]:
    return {
        f"E_{lump}": float(np.exp(np.mean([np.log(folds[e]) for e in members])))
        for lump, members in G_ENZYMES.items()
    }


def steady_state(v: Variant, model: Model, y0: dict[str, float]) -> pd.Series | None:
    res = v.runner.steady_state(_with_y0(model, y0))
    if not hasattr(res.value, "fluxes"):
        return None
    row = pd.concat([res.value.get_args().iloc[-1], res.value.fluxes.iloc[-1]])
    return row[~row.index.duplicated()]


def _with_y0(model: Model, y0: dict[str, float]) -> Model:
    m = copy.deepcopy(model)
    m.update_variables({k: v for k, v in y0.items() if k in m.get_variable_names()})
    return m


def optimise(v: Variant) -> dict:
    model = copy.deepcopy(v.model)
    y0 = model.get_initial_conditions()
    rounds, chosen = [], []
    for k in range(ROUNDS + 1):
        ss = steady_state(v, model, y0)
        if ss is None:
            rounds.append({"round": k, "failed": True})
            break
        y0 = {n: float(ss[n]) for n in model.get_variable_names()}
        entry = {
            "round": k,
            "overexpressed": list(chosen),
            "flux_sq": float(ss[v.flux_sq]),
            "flux_bio": float(ss[v.flux_bio]),
            "growth_rate_per_h": float(ss["mu" if "mu" in ss else "Growth_rate"]) * 60,
            "sq_content_mM": float(v.sq_content(ss)),
        }
        entry["sq_content_g_per_L"] = entry["sq_content_mM"] * M_SQUALENE / 1000
        candidates = [p for p in v.parameters if p not in chosen]
        if k == ROUNDS or not candidates:
            rounds.append(entry)
            break
        rc = mca.response_coefficients(
            _with_y0(model, y0), to_scan=candidates, parallel=False, runner=v.runner
        )
        frc = rc.fluxes.T  # rows: parameters
        score = frc[v.flux_sq] + frc[v.flux_bio]
        entry["response_sq"] = frc[v.flux_sq].to_dict()
        entry["response_bio"] = frc[v.flux_bio].to_dict()
        rounds.append(entry)
        if score.isna().all():
            break
        best = str(score.idxmax())
        model.update_parameter(
            best, model.get_parameter_values()[best] * v.parameters[best]
        )
        chosen.append(best)
        print(f"  {v.name} round {k}: {best} (score {score[best]:.3f})", flush=True)
    return {"folds": v.parameters, "rounds": rounds}


def main() -> None:
    all_variants = ["v0", "v1", "g-none", "g-coupled", "g-photoacclimation",
                    "g-capacity", "v3-drain", "v3-proteome", "v3-drain-cap",
                    "v3-proteome-cap"]  # fmt: skip
    selected = sys.argv[1:] or all_variants
    results = json.loads(OUT.read_text()) if OUT.exists() else {}
    for v in variants(selected):
        t = time.perf_counter()
        print(v.name, flush=True)
        results[v.name] = optimise(v) | {"wall_s": time.perf_counter() - t}
        OUT.write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
