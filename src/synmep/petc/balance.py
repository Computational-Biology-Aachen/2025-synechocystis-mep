"""PETC forms that respond to the PS1/PS2 excitation balance (ADR 0006 §2.4).

Inputs (pfennig2024 surrogate inputs): absorbed light per complex a_PS1, a_PS2, a_PBS
[µmol(photons) mmol(Chl)^-1 s^-1] and the OCP light. Outputs J_3PGA, J_ATP
[mmol mol(Chl)^-1 s^-1].

Common excitation model: phycobilisomes pass a share phi of their excitation to PS2 and
1 - phi to PS1, reduced by OCP quenching (blue-green light):

    b_eff = b / (1 + k_ocp L_OCP)
    e1 = a_PS1 + (1 - phi) b_eff a_PBS
    e2 = w2 (a_PS2 + phi b_eff a_PBS)

- **m0** (`minimal-m0`, reference): J_3PGA = J_max tanh(I / I_k) - r with one weighted
  sum I = a_PS1 + w_PS2 a_PS2 + w_PBS a_PBS. Light is interchangeable.
- **m1**: linear electron flow needs both photosystems,
  J_3PGA = J_max tanh(smin(e1, e2) / I_k) - r, smin(x, y) = (x^-p + y^-p)^(-1/p).
- **m2**: mini electron transport with reduced fractions q (PQ), s (P700/PC) and f (Fd):
  v2 = e2 (1 - q), v_b = k_b q (1 - s), v1 = e1 s (1 - f), v_c = V f / (K + f).
  At steady state all four equal J; eliminating q, s, f leaves one monotone equation in
  J, solved here by bisection (in a model the three fractions are fast states).

J_ATP in m1/m2: c0 + c1 J + c2 J^2 + c3 x, x = excess PS1 excitation (cyclic flow),
m1: e1 - smin(e1, e2) relative to I_k; m2: unused PS1 capacity e1 s - J.
"""

from __future__ import annotations

import numpy as np

INPUTS = ("absorption_ps1", "absorption_ps2", "absorption_pbs", "light_OCP")

M1 = ("J_max", "I_k", "w2", "b", "phi", "k_ocp", "p", "r")
M2 = ("w1", "w2", "b", "phi", "k_ocp", "k_b", "V", "K", "r")
ATP = ("c0", "c1", "c2", "c3")


def excitation(p: dict[str, float], x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    a1, a2, apbs, ocp = x.T
    b_eff = p["b"] / (1 + p["k_ocp"] * ocp)
    e1 = p.get("w1", 1.0) * (a1 + (1 - p["phi"]) * b_eff * apbs)
    e2 = p["w2"] * (a2 + p["phi"] * b_eff * apbs)
    return e1, e2


def smin(x: np.ndarray, y: np.ndarray, power: float) -> np.ndarray:
    x, y = np.maximum(x, 1e-12), np.maximum(y, 1e-12)
    return (x**-power + y**-power) ** (-1 / power)


def m1(p: dict[str, float], x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """J_3PGA and the excess-PS1 term of m1."""
    e1, e2 = excitation(p, x)
    m = smin(e1, e2, p["p"])
    return p["J_max"] * np.tanh(m / p["I_k"]) - p["r"], (e1 - m) / p["I_k"]


def m2_flux(p: dict[str, float], e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    """Linear electron flux J of m2 (bisection on the monotone steady-state residual)."""
    k_b, v, k = p["k_b"], p["V"], p["K"]
    hi = np.minimum(np.minimum(e2, v), k_b) * (1 - 1e-12)
    lo = np.zeros_like(hi)

    def residual(j: np.ndarray) -> np.ndarray:
        q = 1 - j / e2
        s = 1 - j / (k_b * q)
        f = k * j / (v - j)
        return e1 * s * (1 - f) - j  # > 0: PS1 could carry more than j

    for _ in range(60):
        mid = 0.5 * (lo + hi)
        pos = residual(mid) > 0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    return 0.5 * (lo + hi)


def m2(p: dict[str, float], x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """J_3PGA and the excess-PS1 term of m2."""
    e1, e2 = excitation(p, x)
    j = m2_flux(p, e1, e2)
    q = 1 - j / e2
    s = 1 - j / (p["k_b"] * q)
    return j - p["r"], (e1 * s - j) / p["V"]


def m0(p: dict[str, float], x: np.ndarray) -> np.ndarray:
    i = x[:, 0] + p["w_ps2"] * x[:, 1] + p["w_pbs"] * x[:, 2]
    return p["J_max"] * np.tanh(i / p["I_k"]) - p["r"]


def j_atp(c: dict[str, float], j3: np.ndarray, excess: np.ndarray) -> np.ndarray:
    return c["c0"] + c["c1"] * j3 + c["c2"] * j3**2 + c["c3"] * excess
