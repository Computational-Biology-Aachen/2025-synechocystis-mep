"""Train the PETC surrogate from the clamped, native-lowerable pfennig2024 (ADR 0004).

    uv run --extra surrogate python scripts/train_petc_surrogate.py

Samples absorbed light per complex around the operating point of the legacy model
(overall light scale 0.03-30x; complex ratios 0.5-2x, which covers the legacy runs:
PS2/PS1 is constant, PBS/PS1 and OCP/PS1 vary by < 2x), solves pfennig2024 to steady
state with the native runner and fits a polynomial in the (clamped) log inputs.
Writes data/surrogates/petc_training.csv, data/surrogates/petc_poly.json and
docs/reports/petc-surrogate.json.
"""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import scan

from synmep.runners import NATIVE_RUNNER
from synmep.petc.pfennig import LIGHT_INPUTS, OUTPUTS, get_petc_model, outputs
from synmep.surrogates.poly import fit_polynomial, predict

warnings.filterwarnings("ignore")
OUT = Path("data/surrogates")
N_SAMPLES = 4000
DEGREE = 5
SCALE = (0.03, 30.0)
RATIO = (0.5, 2.0)
SEED = 0

# median operating point of the legacy model in the Fig 9 batch (µmol mmol(Chl)^-1 s^-1)
REFERENCE = {
    "absorption_ps1": 72.47,
    "absorption_ps2": 4.25,
    "absorption_pbs": 203.6,
    "light_OCP": 2.49,
}


def sample(n: int, rng: np.random.Generator) -> pd.DataFrame:
    scale = 10 ** rng.uniform(*np.log10(SCALE), n)
    out = {"absorption_ps1": REFERENCE["absorption_ps1"] * scale}
    for name in LIGHT_INPUTS[1:]:
        ratio = 10 ** rng.uniform(*np.log10(RATIO), n)
        out[name] = REFERENCE[name] * scale * ratio
    return pd.DataFrame(out)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    x = sample(N_SAMPLES, rng)
    t = time.perf_counter()
    res = scan.steady_state(
        get_petc_model(),
        to_scan=x,
        parallel=False,  # the loky-parallel scan hung (ADR 0004); serial takes ~6 min
        runner=NATIVE_RUNNER,
    )
    wall = time.perf_counter() - t
    y = outputs(res.fluxes).reset_index(drop=True)
    data = pd.concat([x, y], axis=1)
    ok = data.notna().all(axis=1)
    data.to_csv(OUT / "petc_training.csv", index=False)
    data = data[ok].reset_index(drop=True)

    idx = rng.permutation(len(data))
    n_test = len(data) // 5
    test, train = data.iloc[idx[:n_test]], data.iloc[idx[n_test:]]
    poly = fit_polynomial(train[list(LIGHT_INPUTS)], train[list(OUTPUTS)], DEGREE)
    report = {
        "n_samples": N_SAMPLES,
        "n_converged": int(ok.sum()),
        "wall_s": wall,
        "degree": DEGREE,
        "light_scale_range": SCALE,
        "ratio_range": RATIO,
        "n_terms": len(poly["powers"]),
    }
    pred = predict(poly, test[list(LIGHT_INPUTS)])
    for o in OUTPUTS:
        err = pred[o].to_numpy() - test[o].to_numpy()
        scale = np.abs(test[o]).max()
        report[o] = {
            "rmse": float(np.sqrt(np.mean(err**2))),
            "max_abs_err": float(np.abs(err).max()),
            "max_abs_err_rel_to_range": float(np.abs(err).max() / scale),
            "r2": float(1 - np.sum(err**2) / np.sum((test[o] - test[o].mean()) ** 2)),
        }
    (OUT / "petc_poly.json").write_text(json.dumps(poly, indent=1))
    Path("docs/reports/petc-surrogate.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
