"""Polynomial surrogates in clamped, scaled log inputs (ADR 0004).

Inputs are log-transformed, clamped to the training domain and scaled to [-1, 1]:
`z = (clip(log x, lo, hi) - mid) / half`. Clamping makes predictions constant outside
the domain instead of extrapolating (the legacy surrogate's failure mode, ADR 0002).
The representation is plain JSON so it can be committed and inspected.
"""

from __future__ import annotations

import itertools
from typing import Any

import numpy as np
import pandas as pd
import sympy as sp

Poly = dict[str, Any]


def _powers(n_inputs: int, degree: int) -> list[tuple[int, ...]]:
    return [
        p
        for p in itertools.product(range(degree + 1), repeat=n_inputs)
        if sum(p) <= degree
    ]


def _scaled(poly: Poly, x: pd.DataFrame) -> np.ndarray:
    logx = np.log(np.clip(x[poly["inputs"]].to_numpy(dtype=float), 1e-300, None))
    lo, hi = np.array(poly["log_lo"]), np.array(poly["log_hi"])
    return (np.clip(logx, lo, hi) - (hi + lo) / 2) / ((hi - lo) / 2)


def _design(z: np.ndarray, powers: list[tuple[int, ...]]) -> np.ndarray:
    return np.column_stack([np.prod(z ** np.array(p), axis=1) for p in powers])


def fit_polynomial(x: pd.DataFrame, y: pd.DataFrame, degree: int) -> Poly:
    logx = np.log(x.to_numpy(dtype=float))
    poly: Poly = {
        "inputs": list(x.columns),
        "outputs": list(y.columns),
        "degree": degree,
        "log_lo": logx.min(axis=0).tolist(),
        "log_hi": logx.max(axis=0).tolist(),
        "powers": [list(p) for p in _powers(x.shape[1], degree)],
    }
    a = _design(_scaled(poly, x), [tuple(p) for p in poly["powers"]])
    coef, *_ = np.linalg.lstsq(a, y.to_numpy(dtype=float), rcond=None)
    poly["coef"] = {o: coef[:, j].tolist() for j, o in enumerate(y.columns)}
    return poly


def predict(poly: Poly, x: pd.DataFrame) -> pd.DataFrame:
    a = _design(_scaled(poly, x), [tuple(p) for p in poly["powers"]])
    return pd.DataFrame(
        {o: a @ np.array(c) for o, c in poly["coef"].items()}, index=x.index
    )


def expressions(poly: Poly, inputs: dict[str, sp.Expr]) -> dict[str, sp.Expr]:
    """Sympy expressions of all outputs, given expressions for the inputs."""
    z = []
    for name, lo, hi in zip(
        poly["inputs"], poly["log_lo"], poly["log_hi"], strict=True
    ):
        logx = sp.Min(sp.Max(sp.log(inputs[name]), lo), hi)
        z.append((logx - (hi + lo) / 2) / ((hi - lo) / 2))
    monomials = [
        sp.Mul(*(zi**p for zi, p in zip(z, powers, strict=True) if p))
        for powers in poly["powers"]
    ]
    return {
        o: sp.Add(*(float(c) * m for c, m in zip(coef, monomials, strict=True)))
        for o, coef in poly["coef"].items()
    }
