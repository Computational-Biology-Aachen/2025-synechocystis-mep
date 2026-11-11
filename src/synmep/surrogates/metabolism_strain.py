"""Per-strain metabolism surrogates for the `double` family (ADR 0006 checkpoint A).

A surrogate with all six overexpressable enzymes as inputs reached only 9–20 % median /
40–100 % p95 error (`docs/reports/metabolism-surrogate_enz*.json`), the WT surrogate
3 %. Goal 1 only needs the discrete strains of the MCA sequence, so each strain (enzyme
folds) gets its own surrogate in the two supplies S_PGA, S_ATP, trained on demand from
the metabolism ODE and cached in `data/surrogates/metabolism/`.

Supplies are sampled as in `scripts/train_metabolism_surrogate.py`, from the PETC
training states (white light + LED mixtures) times a PETC efficiency 0.3–10.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from synmep.data import DATA
from synmep.models.metabolism import (
    BIOMASS,
    REFERENCE,
    SUPPLIES,
    SURROGATE_OUTPUTS,
    get_model,
    reference_parameters,
)
from synmep.runners import steady_state_scan
from synmep.surrogates.poly import fit_polynomial, predict

CACHE = DATA / "surrogates" / "metabolism"
PETC_SCALE_RANGE = (0.3, 10.0)
FLOOR = 1e-3  # [mM min^-1] supply floor, as in double_growth
DEGREES = (3, 4, 5, 6)


def key(folds: dict[str, float], coupling: float) -> str:
    text = json.dumps(
        {
            "biomass": BIOMASS,
            "coupling": coupling,
            "folds": {k: round(v, 6) for k, v in sorted(folds.items())},
        }
    )
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def sample(n: int, rng: np.random.Generator) -> pd.DataFrame:
    s = DATA / "surrogates"
    petc = pd.read_csv(s / "petc_training.csv").dropna()
    led = pd.read_csv(s / "petc_training_led.csv")
    petc = pd.concat([petc, led[petc.columns].dropna()], ignore_index=True)
    rows = petc.iloc[rng.integers(0, len(petc), n)].reset_index(drop=True)
    scale = 10 ** rng.uniform(*np.log10(PETC_SCALE_RANGE), n)
    conv = rows["chl"].to_numpy() * 0.06 * scale
    return pd.DataFrame(
        {
            "S_PGA": np.maximum(rows["J_3PGA"].to_numpy() * conv, FLOOR),
            "S_ATP": np.maximum(
                rows["J_ATP"].to_numpy()
                * conv
                * reference_parameters(REFERENCE)["atp_scale"],
                FLOOR,
            ),
        }
    )


def train(
    folds: dict[str, float], coupling: float, n: int = 3000, seed: int = 0
) -> dict:
    rng = np.random.default_rng(seed)
    x = sample(n, rng)
    m = get_model(pigment_coupling=coupling)
    p0 = m.get_parameter_values()
    m.update_parameters({k: p0[k] * f for k, f in folds.items()})
    args, fluxes = steady_state_scan(m, x)
    mu = args["mu"]
    y = pd.DataFrame(
        {
            "mu": mu,
            "Chl_star": fluxes["v_chl"] / mu,
            "Car_star": fluxes["v_car"] / mu,
            "SQ_star": fluxes["v_sq"] / mu,
        }
    )
    data = pd.concat([x, y], axis=1)
    data = data[data.notna().all(axis=1) & (y > 0).all(axis=1)].reset_index(drop=True)
    idx = rng.permutation(len(data))
    test, tr = data.iloc[idx[: len(data) // 5]], data.iloc[idx[len(data) // 5 :]]
    best = None
    for degree in DEGREES:
        poly = fit_polynomial(
            tr[list(SUPPLIES)], tr[list(SURROGATE_OUTPUTS)], degree, log_outputs=True
        )
        pred = predict(poly, test[list(SUPPLIES)])
        err = {
            o: float(np.quantile(np.abs(pred[o] / test[o] - 1), 0.95))
            for o in SURROGATE_OUTPUTS
        }
        if best is None or max(err.values()) < max(best[1].values()):
            best = (poly, err, degree)
    assert best is not None
    poly, err, degree = best
    return poly | {
        "folds": folds,
        "coupling": coupling,
        "n_converged": len(data),
        "degree": degree,
        "p95_abs_rel_err": err,
    }


def load_or_train(folds: dict[str, float], coupling: float) -> dict:
    path = CACHE / f"{key(folds, coupling)}.json"
    if path.exists():
        return json.loads(path.read_text())
    poly = train(folds, coupling)
    CACHE.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(poly, indent=1))
    tmp.replace(path)
    return poly


def path_for(folds: dict[str, float], coupling: float) -> Path:
    return CACHE / f"{key(folds, coupling)}.json"
