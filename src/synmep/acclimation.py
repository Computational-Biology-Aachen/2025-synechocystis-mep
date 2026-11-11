"""Capacity acclimation as an option for every family (ADR 0006 §2.3).

Photosynthetic capacity Cap scales the PETC supply and relaxes over hours towards
Cap* = (1 + K) L_rel / (K + L_rel), L_rel = absorbed light per PS1 Chl relative to the
Fig 9 median (`L_REF`), so Cap* = 1 at L_rel = 1 and grows with light (K, tau fitted).
`minimal` has its own implementation of the same rule (ADR 0005); `single-legacy`
stays faithful and has none.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import sympy as sp

if TYPE_CHECKING:
    from mxlpy import Model

L_REF = 72.47  # [µmol(photons) mmol(Chl)^-1 s^-1] PS1 absorption at the Fig 9 median
PARAMETERS = {"L_ref": L_REF, "capacity_K": 1.0, "capacity_tau": 600.0}  # tau [min]


def relax(x: float, target: float, rate: float) -> float:
    return rate * (target - x)


def inverse(x: float) -> float:
    return 1 / x


def target_expressions(abs_ps1: str = "complex_abs_ps1_mean") -> dict[str, sp.Expr]:
    """`L_rel` and `Cap_target` from the (depth-mean) PS1 absorption."""
    k, l_rel = sp.Symbol("capacity_K"), sp.Symbol("L_rel")
    return {
        "L_rel": sp.Symbol(abs_ps1) / sp.Symbol("L_ref"),
        "Cap_target": (1 + k) * l_rel / (k + l_rel),
    }


def add_state(m: Model, rate: str | None = None) -> None:
    """Variable Cap (start 1) relaxing to `Cap_target` at `rate`, default 1/tau.

    The steady state doesn't depend on the rate; steady modes may pass a fast one.
    """
    m.add_variable("Cap", 1.0)
    if rate is None:
        m.add_derived("capacity_rate", fn=inverse, args=["capacity_tau"])
        rate = "capacity_rate"
    m.add_reaction(
        "net_Cap", fn=relax, args=["Cap", "Cap_target", rate], stoichiometry={"Cap": 1}
    )
