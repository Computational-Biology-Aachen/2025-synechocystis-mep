"""Closed-form PETC for the minimal model (ADR 0005), fitted to the pfennig2024 data.

    uv run python scripts/fit_petc_closed_form.py

Data: data/surrogates/petc_training.csv (pfennig2024 steady states at physically
realisable absorbed-light combinations, ADR 0004). Model:

    I      = a_PS1 + w_PS2 a_PS2 + w_PBS a_PBS          [absorbed light, PS1 units]
    J_3PGA = J_max tanh(I / I_k) - r                    [mmol mol(Chl)^-1 s^-1]
    J_ATP  = c0 + c1 J_3PGA + c2 J_3PGA^2

J_ATP is almost a function of J_3PGA in pfennig2024, so the PETC has one saturating
light response and a fixed ATP:3PGA relation. Writes data/surrogates/petc_closed_form.json
and docs/reports/petc-closed-form.json (hold-out accuracy vs. the degree-5 polynomial).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

DATA = Path("data/surrogates/petc_training.csv")
OUT = Path("data/surrogates/petc_closed_form.json")
REPORT = Path("docs/reports/petc-closed-form.json")
INPUTS = ["absorption_ps1", "absorption_ps2", "absorption_pbs"]


def j_3pga(p: dict[str, float], x: np.ndarray) -> np.ndarray:
    i = x[:, 0] + p["w_ps2"] * x[:, 1] + p["w_pbs"] * x[:, 2]
    return p["J_max"] * np.tanh(i / p["I_k"]) - p["r"]


def j_atp(p: dict[str, float], j3: np.ndarray) -> np.ndarray:
    return p["c0"] + p["c1"] * j3 + p["c2"] * j3**2


def r2(pred: np.ndarray, y: np.ndarray) -> float:
    return float(1 - np.sum((pred - y) ** 2) / np.sum((y - y.mean()) ** 2))


def main() -> None:
    d = pd.read_csv(DATA).dropna()
    x, y3, ya = d[INPUTS].to_numpy(), d["J_3PGA"].to_numpy(), d["J_ATP"].to_numpy()
    idx = np.random.default_rng(0).permutation(len(d))
    tr, te = idx[: int(0.8 * len(d))], idx[int(0.8 * len(d)) :]
    names = ["J_max", "I_k", "w_ps2", "w_pbs", "r"]

    def unpack(logp: np.ndarray) -> dict[str, float]:
        return dict(zip(names, np.exp(logp), strict=True))

    fit = least_squares(
        lambda lp: j_3pga(unpack(lp), x[tr]) - y3[tr],
        np.log([20.0, 2e4, 10.0, 15.0, 0.3]),
    )
    p = unpack(fit.x)
    c2, c1, c0 = np.polyfit(y3[tr], ya[tr], 2)
    p |= {"c0": float(c0), "c1": float(c1), "c2": float(c2)}
    pred3 = j_3pga(p, x[te])
    preda = j_atp(p, pred3)
    poly = json.loads(Path("docs/reports/petc-surrogate.json").read_text())
    report = {
        "parameters": p,
        "n_train": len(tr),
        "n_test": len(te),
        "J_3PGA": {
            "r2": r2(pred3, y3[te]),
            "rmse": float(np.sqrt(np.mean((pred3 - y3[te]) ** 2))),
        },
        "J_ATP": {
            "r2": r2(preda, ya[te]),
            "rmse": float(np.sqrt(np.mean((preda - ya[te]) ** 2))),
        },
        "polynomial_degree5": {k: poly[k]["r2"] for k in ("J_3PGA", "J_ATP")},
    }
    OUT.write_text(json.dumps({"inputs": INPUTS, **p}, indent=1))
    REPORT.write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
