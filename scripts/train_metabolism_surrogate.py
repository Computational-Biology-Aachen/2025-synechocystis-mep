"""Train the metabolism QSS surrogate from the g_metabolism ODE model (ADR 0004).

    uv run python scripts/train_metabolism_surrogate.py

Supply inputs are *physically produced* pairs: (J_3PGA, J_ATP, Chl) triples from the
PETC training data (culture states pushed through the light model and pfennig2024)
times a PETC efficiency factor (0.3-10x, fitted in g_growth):
S = J × Chl × 0.06 × petc_scale. Sampling S_PGA and S_ATP independently (or with a
fixed ratio band) left real operating points outside the training data (ADR 0004).
The CrtE level is sampled 0.01-3x WT (Fig 11 knock-down); all other enzymes are at WT
level. The metabolism layer is solved to steady state with the native runner and a
polynomial in clamped log inputs is fitted to log outputs: mu and the steady-state
contents X* = v_X / mu (Chl, carotenoids, squalene).
Writes data/surrogates/metabolism_training.csv, data/surrogates/metabolism_poly.json
and docs/reports/metabolism-surrogate.json.
"""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from synmep.models.g_metabolism import (
    SUPPLIES,
    SURROGATE_OUTPUTS,
    get_model,
)
from synmep.runners import steady_state_scan
from synmep.surrogates.poly import fit_polynomial, predict

warnings.filterwarnings("ignore")
OUT = Path("data/surrogates")
N_SAMPLES = 4000
DEGREES = (3, 4, 5, 6)
PETC_SCALE_RANGE = (0.3, 10.0)
SEED = 0
S_PGA_FLOOR = 1e-3  # [mM min^-1] as in g_growth (net respiration at very low light)
S_ATP_FLOOR = 1e-3
# An 8-input surrogate (all enzyme levels 0.1-10x) reached only ~10-14 % median /
# ~50 % p95 relative error (degree 4, 6000 samples). Fig 12 (MCA over all enzymes) uses
# the coupled ODE model instead, so only CrtE (Fig 11 knock-down) is an input here.
SURROGATE_ENZYMES = {"CrtE": (0.01, 3.0)}


def sample(n: int, rng: np.random.Generator) -> pd.DataFrame:
    petc = pd.read_csv(OUT / "petc_training.csv").dropna()
    rows = petc.iloc[rng.integers(0, len(petc), n)].reset_index(drop=True)
    scale = 10 ** rng.uniform(*np.log10(PETC_SCALE_RANGE), n)
    conv = rows["chl"].to_numpy() * 0.06 * scale
    out = {
        "S_PGA": np.maximum(rows["J_3PGA"].to_numpy() * conv, S_PGA_FLOOR),
        "S_ATP": np.maximum(rows["J_ATP"].to_numpy() * conv, S_ATP_FLOOR),
    }
    out |= {
        f"E_{e}": 10 ** rng.uniform(*np.log10(r), n)
        for e, r in SURROGATE_ENZYMES.items()
    }
    return pd.DataFrame(out)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    x = sample(N_SAMPLES, rng)
    t = time.perf_counter()
    args, fluxes = steady_state_scan(get_model(), x)
    wall = time.perf_counter() - t
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
    data.to_csv(OUT / "metabolism_training.csv", index=False)
    ok = data.notna().all(axis=1) & (data[list(SURROGATE_OUTPUTS)] > 0).all(axis=1)
    data = data[ok].reset_index(drop=True)

    inputs = [*SUPPLIES, *(f"E_{e}" for e in SURROGATE_ENZYMES)]
    idx = rng.permutation(len(data))
    n_test = len(data) // 5
    test, train = data.iloc[idx[:n_test]], data.iloc[idx[n_test:]]
    report: dict = {
        "n_samples": N_SAMPLES,
        "n_converged": int(ok.sum()),
        "wall_s": wall,
    }
    best = None
    for degree in DEGREES:
        poly = fit_polynomial(
            train[inputs], train[list(SURROGATE_OUTPUTS)], degree, log_outputs=True
        )
        pred = predict(poly, test[inputs])
        entry = {"n_terms": len(poly["powers"])}
        for o in SURROGATE_OUTPUTS:
            rel = pred[o].to_numpy() / test[o].to_numpy() - 1
            entry[o] = {
                "median_abs_rel_err": float(np.median(np.abs(rel))),
                "p95_abs_rel_err": float(np.quantile(np.abs(rel), 0.95)),
            }
        report[f"degree_{degree}"] = entry
        worst = max(entry[o]["p95_abs_rel_err"] for o in SURROGATE_OUTPUTS)
        if best is None or worst < best[0]:
            best = (worst, degree, poly)
    assert best is not None
    report["chosen_degree"] = best[1]
    (OUT / "metabolism_poly.json").write_text(json.dumps(best[2], indent=1))
    Path("docs/reports/metabolism-surrogate.json").write_text(
        json.dumps(report, indent=1)
    )
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
