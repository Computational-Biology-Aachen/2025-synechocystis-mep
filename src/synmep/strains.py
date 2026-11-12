"""MCA strain optimisation for any registered variant (ADR 0006 §2.2).

Paper protocol (Fig 12, report 04): steady growth at 1e7 cells mL^-1,
100 µmol m^-2 s^-1, 1 cm, warm white; each round the enzyme with the largest summed
response of squalene and biomass synthesis is overexpressed by its proteomics-based
fold, each enzyme at most once, only steps with positive control on the score and on
intracellular squalene that raise it when realised.

Unlike `scripts/strain_optimisation.py` (the original per-family script, report 04)
this works on the fitted models through `synmep.evaluate` adapters. Response
coefficients are finite differences of steady states (1 %), and the score is
R(SQ) + 2 R(mu), which at steady state equals R(v_sq) + R(v_bio) (v_sq = SQ mu,
v_bio ∝ mu) in every family, including `full`, whose fluxes are per layer.
"""

from __future__ import annotations

from functools import cache

import numpy as np

from synmep import evaluate
from synmep.data import load_legacy_inputs

LEGACY_ENZYMES = ("DXS", "ThiG", "DXR", "IspD", "IspE", "IspF", "IspG", "IspH",
                  "Idi", "CrtE", "Sqs", "Shc", "ChlP", "CrtB", "ChlG")  # fmt: skip
LUMPS = {  # lump (`double`, `minimal`, `full`) -> legacy enzymes it stands for
    "DXS": ("DXS",),
    "MEP": ("DXR", "IspD", "IspE", "IspF", "IspG", "IspH"),
    "CrtE": ("CrtE",),
    "Sqs": ("Sqs",),
    "CrtB": ("CrtB",),
    "Chl": ("ChlP", "ChlG"),
}
MCA = {"density": 1e7, "light": 100.0, "sample_depth": 0.01}
STEP = 0.01  # relative perturbation for the response coefficients
ROUNDS = 10


@cache
def legacy_folds() -> dict[str, float]:
    """fold = 1 + 1.331057 kDa / (2 E0), E0 the legacy (2 fL) proteomics level [µM]."""
    inp = load_legacy_inputs()
    poi = inp.proteins_of_interest
    folds = {}
    for e in LEGACY_ENZYMES:
        kda = poi[poi["Enzyme"] == e]["kDa"]
        if len(kda) and e in inp.enzyme_concentrations_umol:
            e0 = inp.enzyme_concentrations_umol[e]
            folds[e] = 1 + 1.331057 * float(kda.iloc[0]) / (2 * e0)
    return folds


def enzyme_folds(a: evaluate.Adapter) -> dict[str, float]:
    """Overexpressable parameter -> fold for this variant (geometric mean for lumps)."""
    folds = legacy_folds()
    if a.name.startswith("single-"):
        m = a.build("steady", light=MCA["light"], sample_depth=MCA["sample_depth"])
        names = set(m.get_parameter_names())
        return {f"{e}_ksynth": f for e, f in folds.items() if f"{e}_ksynth" in names}
    return {
        f"E_{lump}": float(np.exp(np.mean([np.log(folds[e]) for e in members])))
        for lump, members in LUMPS.items()
    }


def _steady(a, folds, pars, y0=None):
    light = {k: v for k, v in MCA.items() if k != "density"}
    m = evaluate.steady_model(a, folds, MCA["density"], pars=pars, **light)
    if y0:
        m.update_variables({k: v for k, v in y0.items() if k in m.get_variable_names()})
    row = evaluate.steady_state(a, m, growing=True)
    if row is None:
        return None, None
    return row, {k: float(row[k]) for k in m.get_variable_names()}


def optimise(a: evaluate.Adapter, pars: dict[str, float] | None = None) -> dict:
    """Overexpression sequence; `pars` are absolute parameter values (fits)."""
    pars = dict(pars or {})
    folds_all = enzyme_folds(a)
    sq, mu = a.keys["sq"], a.keys["mu"]
    chosen: dict[str, float] = {}
    rounds = []
    row, y0 = _steady(a, chosen, pars)
    for k in range(ROUNDS + 1):
        if row is None:
            rounds.append({"round": k, "failed": True})
            break
        entry = {
            "round": k,
            "overexpressed": list(chosen),
            "sq_mM": float(row[sq]),
            "mu_per_h": float(row[mu]) * 60,
        }
        candidates = [p for p in folds_all if p not in chosen]
        if k == ROUNDS or not candidates:
            rounds.append(entry)
            break
        score, sq_resp = {}, {}
        for p in candidates:
            r1, _ = _steady(a, chosen | {p: 1 + STEP}, pars, y0)
            if r1 is None:
                continue
            d_sq = np.log(float(r1[sq]) / float(row[sq])) / np.log(1 + STEP)
            d_mu = np.log(float(r1[mu]) / float(row[mu])) / np.log(1 + STEP)
            sq_resp[p], score[p] = d_sq, d_sq + 2 * d_mu
        entry["score"], entry["response_sq"] = score, sq_resp
        eligible = sorted(
            (p for p in score if score[p] > 0 and sq_resp[p] > 0),
            key=lambda p: -score[p],
        )
        best, rejected = None, []
        for p in eligible:
            trial = chosen | {p: folds_all[p]}
            r1, y1 = _steady(a, trial, pars, y0)
            if r1 is not None and float(r1[sq]) > float(row[sq]):
                best, row, y0 = p, r1, y1
                break
            rejected.append(p)
        if rejected:
            entry["rejected"] = rejected
        rounds.append(entry)
        if best is None:
            entry["stop"] = (
                "no enzyme with positive control" if not eligible
                else "every eligible step lowered the squalene content"
            )  # fmt: skip
            break
        chosen[best] = folds_all[best]
    return {"folds": folds_all, "sequence": list(chosen), "rounds": rounds}


def strains(result: dict) -> list[dict[str, float]]:
    """Strain k = the first k overexpressions (k = 0: wild type)."""
    seq, folds = result["sequence"], result["folds"]
    return [{p: folds[p] for p in seq[:k]} for k in range(len(seq) + 1)]
