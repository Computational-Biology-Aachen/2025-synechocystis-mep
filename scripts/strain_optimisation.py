"""MCA-based strain optimisation (notebook / paper Fig 12) for all model variants.

    uv run python scripts/strain_optimisation.py [variant ...]

Each round: flux response coefficients (normalised) of squalene synthesis and biomass
synthesis to every enzyme parameter at steady growth; the enzyme with the largest sum
is overexpressed; **each enzyme at most once**. A step is only taken if its control is
positive: the summed response must be > 0 **and** the response of the intracellular
squalene content (state `SQ`) must be > 0 (an overexpression that lowers the squalene
content is pointless). Among those passing both, the best is taken unless the realised
overexpression (a large step, where the local control can mislead) lowers the
steady-state squalene content; then it is rejected (`rejected`) and the next one tried.
If none is left, the sequence stops (`stop`). The first round is the answer to
"which single overexpression is best".

Overexpression is a fold change of the enzyme (`single`: its synthesis rate, `double`
and `minimal`: its level): fold = 1 + 1.331057 kDa / (2 E0) with E0 the legacy (2 fL) proteomics
concentration [µM]. For single-legacy this equals the notebook's additive
`ksynth + 1.331057 kDa / (5.5 h)` exactly. The lumped MEP chain (`double`, `minimal`) gets the
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
SQ_STATE = "SQ"  # intracellular squalene [mM] (all variants)
M_SQUALENE = 410.73  # [g mol^-1]

LEGACY_ENZYMES = ("DXS", "ThiG", "DXR", "IspD", "IspE", "IspF", "IspG", "IspH",
                  "Idi", "CrtE", "Sqs", "Shc", "ChlP", "CrtB", "ChlG")  # fmt: skip
G_ENZYMES = {  # lump (`double`, `minimal`) -> legacy enzymes it stands for
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
    if "single-legacy" in selected:
        from synmep.models import single_legacy

        m = single_legacy.get_model(load_legacy_inputs(), "dshc", "steady")
        pars = {f"{e}_ksynth": f for e, f in folds.items()
                if f"{e}_ksynth" in m.get_parameter_names()}  # fmt: skip
        out.append(Variant("single-legacy", m, pars, "v_Sqs", "v_Biomass", NATIVE_RUNNER,
                           lambda r: r["v_Sqs"] / r["Growth_rate"]))  # fmt: skip
    if "single-fixed" in selected:
        from synmep.data_single_fixed import load_inputs
        from synmep.models import single_fixed

        m = single_fixed.get_model(load_inputs(), "dshc", "steady")
        pars = {f"{e}_ksynth": f for e, f in folds.items()
                if f"{e}_ksynth" in m.get_parameter_names()}  # fmt: skip
        out.append(Variant("single-fixed", m, pars, "v_Sqs", "v_Biomass", SUNDIALS_RUNNER,
                           lambda r: r["v_Sqs"] / r["Growth_rate"]))  # fmt: skip
    from synmep.models import double

    for v in ("none", "coupled"):
        if f"double-{v}" not in selected:
            continue
        m = double.get_model(load_legacy_inputs(), v)  # type: ignore[arg-type]
        pars = lump_folds(folds)
        out.append(Variant(f"double-{v}", m, pars, "v_sq", "v_bio", NATIVE_RUNNER,
                           lambda r: r["v_sq"] / r["mu"]))  # fmt: skip
    from synmep.models import minimal

    for name in minimal.VARIANTS:
        if name not in selected:
            continue
        m = minimal.get_model(load_legacy_inputs(), **minimal.variant(name))
        out.append(Variant(name, m, lump_folds(folds), "v_sq", "v_bio", NATIVE_RUNNER,
                           lambda r: r["v_sq"] / r["mu"]))  # fmt: skip
    return out


def lump_folds(folds: dict[str, float]) -> dict[str, float]:
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
        sq_content = rc.variables.T[SQ_STATE]
        score = frc[v.flux_sq] + frc[v.flux_bio]
        entry["response_sq"] = frc[v.flux_sq].to_dict()
        entry["response_bio"] = frc[v.flux_bio].to_dict()
        entry["response_sq_content"] = sq_content.to_dict()
        rounds.append(entry)
        eligible = score[(score > 0) & (sq_content > 0)].sort_values(ascending=False)
        # the overexpression is a large step, so the local control can mislead: a step
        # that lowers the realised squalene content is rejected, the next one tried
        best, rejected = None, []
        for p in map(str, eligible.index):
            trial = copy.deepcopy(model)
            trial.update_parameter(p, trial.get_parameter_values()[p] * v.parameters[p])
            ss_trial = steady_state(v, trial, y0)
            if ss_trial is not None and v.sq_content(ss_trial) > entry["sq_content_mM"]:
                best, model = p, trial
                break
            rejected.append(p)
        if rejected:
            entry["rejected"] = rejected
        if best is None:
            entry["stop"] = (
                "no enzyme with positive control on score and squalene content"
                if eligible.empty
                else "every eligible step lowered the squalene content"
            )
            print(f"  {v.name} round {k}: stop", flush=True)
            break
        chosen.append(best)
        print(f"  {v.name} round {k}: {best} (score {score[best]:.3f})", flush=True)
    return {"folds": v.parameters, "rounds": rounds}


def main() -> None:
    all_variants = ["single-legacy", "single-fixed", "double-none", "double-coupled",
                    "minimal-m0", "minimal-m0-cap"]  # fmt: skip
    selected = sys.argv[1:] or all_variants
    results = json.loads(OUT.read_text()) if OUT.exists() else {}
    for v in variants(selected):
        t = time.perf_counter()
        print(v.name, flush=True)
        results[v.name] = optimise(v) | {"wall_s": time.perf_counter() - t}
        OUT.write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
