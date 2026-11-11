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


# --- chromatic acclimation: PS1:PS2 ratio seeks balanced excitation (ADR 0006 §2.6)
# PBS share to PS2 among the attached PBS in pfennig2024 (PBS_PS2 / (PBS_PS1 + PBS_PS2))
PBS_TO_PS2 = 0.51 / (0.39 + 0.51)
CHROMATIC = {
    "chromatic_gain": 0.5,  # [-] exponent g of the excitation imbalance
    "chromatic_tau": 1440.0,  # [min] ~1 day, photosystem stoichiometry acclimation
    "ps_ratio_min": 2.0,  # PS1:PS2 range of the PETC training data (report 08)
    "ps_ratio_max": 12.0,
    # E2/E1 that cells keep: the imbalance at the reference state under warm white
    # (`calibrate_chromatic`), so warm white keeps the measured PS1:PS2 ratio
    "chromatic_reference": 1.0,
}


def chromatic_expressions(prefix_means: str = "complex_abs_") -> dict[str, sp.Expr]:
    """Total PS1 and PS2 excitation and the bounded target of `ps_ratio`.

    E1 = a_PS1 + (1 - phi) a_PBS, E2 = a_PS2 + phi a_PBS (linear electron flow needs
    equal throughput of both), target r* = clamp(r ((E2/E1) / I_ref)^g, r_min, r_max):
    more PS1 when PS2 is more excited than under the reference (PS2 light, e.g. orange),
    less in PS1 light (blue), as observed in cyanobacteria. More PS1 moves Chl absorption
    from PS2 to PS1 and so lowers E2/E1 (negative feedback). I_ref =
    `chromatic_reference`. (Per-complex excitation would give positive feedback: more
    PS1 complexes lower each one's excitation.)
    """
    r = sp.Symbol("ps_ratio")
    a1, a2, apbs = (sp.Symbol(f"{prefix_means}{c}_mean") for c in ("ps1", "ps2", "pbs"))
    e1 = a1 + (1 - PBS_TO_PS2) * apbs
    e2 = a2 + PBS_TO_PS2 * apbs
    imbalance = e2 / e1 / sp.Symbol("chromatic_reference")
    target = r * imbalance ** sp.Symbol("chromatic_gain")
    return {
        "chromatic_imbalance": e2 / e1 / sp.Symbol("chromatic_reference"),
        "ps_ratio_target": sp.Min(
            sp.Max(target, sp.Symbol("ps_ratio_min")), sp.Symbol("ps_ratio_max")
        ),
    }


def add_chromatic_state(m: Model, rate: str | None = None) -> None:
    """Turn the parameter `ps_ratio` into a state relaxing to `ps_ratio_target`."""
    r0 = m.get_parameter_values()["ps_ratio"]
    m.remove_parameter("ps_ratio")
    m.add_variable("ps_ratio", r0)
    if rate is None:
        m.add_derived("chromatic_rate", fn=inverse, args=["chromatic_tau"])
        rate = "chromatic_rate"
    m.add_reaction(
        "net_ps_ratio",
        fn=relax,
        args=["ps_ratio", "ps_ratio_target", rate],
        stoichiometry={"ps_ratio": 1},
    )


def imbalance(args) -> float:
    """E2/E1 (unnormalised) from derived values of a chromatic model at a state."""
    return float(args["chromatic_imbalance"]) * float(args["chromatic_reference"])
