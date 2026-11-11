"""`double` family at steady growth with the explicit metabolism ODE (ADR 0004).

For metabolic control analysis over all enzyme parameters (Fig 12; the metabolism
surrogate covers only the overexpressable lumps within their fold ranges, and MCA
needs exact local derivatives) and turbidostat scans. Cell density N is a parameter
(steady growth, e.g. turbidostat); the PETC still enters as its surrogate, the
metabolism as the full metabolism ODE, and the slow pools as

    dX/dt = v_X - mu X     (X = Chl, Car, SQ)

`coupled`: the pigment-branch enzymes scale with growth. Parameters (petc_scale,
pbs_scale) are taken from the `double` growth fit.
"""

from __future__ import annotations

import json
from pathlib import Path

import sympy as sp
from mxlpy import Model

from synmep import acclimation
from synmep import light as light_sources
from synmep.data import DATA, LegacyInputs
from synmep.exprfn import expr_fn
from synmep.models import metabolism
from synmep.models.double_growth import (
    CELL_VOLUME,
    PETC_POLY,
    PigmentRegulation,
    per_nl,
    petc_supply_expressions,
    scaled,
)


def fitted_parameters(variant: PigmentRegulation) -> dict[str, float]:
    suffix = "" if variant == "none" else f"-{variant}"
    path = Path("docs/reports") / f"double-fit{suffix}.json"
    return json.loads(path.read_text())["fitted"] if path.exists() else {}


def flux_minus_dilution(v: float, x: float, mu: float) -> float:
    return v - mu * x


def get_model(
    inp: LegacyInputs,
    variant: PigmentRegulation = "none",
    *,
    cell_density: float = 1e7,  # [cells mL^-1] legacy Fig 12 condition
    light: float = 100.0,  # legacy Fig 12 condition
    sample_depth: float = 0.01,  # legacy default
    depth_nodes: int = 5,
    spectral_step: int = 5,
    light_channels: tuple[str, ...] | None = None,
    capacity: bool = False,
    chromatic: bool = False,
) -> Model:
    m = metabolism.get_model(pigment_coupling=1.0 if variant == "coupled" else 0.0)
    m.remove_parameters(["S_PGA", "S_ATP"])
    ref = metabolism.REFERENCE
    m.add_variables({"Chl": ref.chl, "Car": ref.carotenoids, "SQ": ref.squalene})
    m.add_parameters(
        {
            "N": cell_density,
            "light": light,
            "sample_depth_m": sample_depth,
            "cell_volume": CELL_VOLUME,
            "phycocyanin_ref": 53180 * 2e-15 / CELL_VOLUME,
            "allophycocyanin_ref": 24076 * 2e-15 / CELL_VOLUME,
            "pbs_scale": 1.0,
            "ps_ratio": 5.9,
            "beta_carotene_fraction": 0.26,
            "petc_scale": 1.0,
            "S_PGA_floor": 1e-3,
        }
    )
    m.update_parameters(
        {
            k: v
            for k, v in fitted_parameters(variant).items()
            if k in m.get_parameter_names()
        }
    )
    m.add_derived("N_nL", fn=per_nl, args=["N"])
    m.add_derived("phycocyanin", fn=scaled, args=["phycocyanin_ref", "pbs_scale"])
    m.add_derived(
        "allophycocyanin", fn=scaled, args=["allophycocyanin_ref", "pbs_scale"]
    )

    if light_channels is not None:
        light_sources.add_parameters(m, light_channels)
    exprs = petc_supply_expressions(
        inp, json.loads((DATA / "surrogates" / PETC_POLY).read_text()),
        depth_nodes, spectral_step, light_channels,
    )  # fmt: skip
    chl, scale = sp.Symbol("Chl"), sp.Symbol("petc_scale")
    if chromatic:
        m.add_parameters(acclimation.CHROMATIC)
        exprs |= acclimation.chromatic_expressions()
        if "k_relax" not in m.get_parameter_names():
            m.add_parameter("k_relax", 1.0)
    if capacity:
        if "k_relax" not in m.get_parameter_names():
            m.add_parameter("k_relax", 1.0)
        m.add_parameters(acclimation.PARAMETERS)
        exprs |= acclimation.target_expressions()
        scale = scale * sp.Symbol("Cap")
    exprs["S_PGA"] = sp.Max(
        sp.Symbol("J_3PGA") * chl * 0.06 * scale, sp.Symbol("S_PGA_floor")
    )
    exprs["S_ATP"] = (
        sp.Max(sp.Symbol("J_ATP"), 0) * chl * 0.06 * scale * sp.Symbol("atp_scale")
    )
    for name, expr in exprs.items():
        fn, args = expr_fn(expr, origin=name)
        m.add_derived(name, fn=fn, args=args)

    for state, flux in (("Chl", "v_chl"), ("Car", "v_car"), ("SQ", "v_sq")):
        m.add_reaction(
            f"net_{state}",
            fn=flux_minus_dilution,
            args=[flux, state, "mu"],
            stoichiometry={state: 1},
        )
    if capacity:  # steady state only: fast relaxation, same fixed point
        acclimation.add_state(m, "k_relax")
    if chromatic:
        acclimation.add_chromatic_state(m, "k_relax")
    return m
