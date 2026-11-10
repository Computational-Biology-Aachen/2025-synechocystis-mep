"""Ground-up growth layer (hours timescale), thread "g" (ADR 0004).

Slow states: cell density `N` [cells mL^-1], chlorophyll `Chl`, carotenoids `Car` and
squalene `SQ` [mM, intracellular]. Faster layers enter as QSS surrogates:

    N, Chl, Car --self-shading--> absorbed light per complex at 5 depths
                --PETC surrogate, averaged over depth--> J_3PGA, J_ATP
                                                  [mmol mol(Chl)^-1 s^-1]
                × Chl × petc_scale --> S_PGA, S_ATP [mM min^-1]
                --metabolism surrogate (+ E_CrtE)--> mu, X* (steady-state contents)

    dN/dt = mu N,   dX/dt = mu (X* - X)   (X = Chl, Car, SQ)

`mu (X* - X)` equals the flux form `v_X - mu X` at quasi-steady state (X* = v_X / mu),
but stays bounded when growth stops (ADR 0004).

Strain: Δshc (squalene isn't converted to hopanoids). Cell volume 4 fL.
"""

from __future__ import annotations

import json
from typing import Literal

import numpy as np
import sympy as sp
from mxlpy import Model

from synmep.data import DATA, LegacyInputs
from synmep.exprfn import expr_fn
from synmep.models.g_metabolism import REFERENCE
from synmep.surrogates.petc_legacy import Symbols, light_expressions
from synmep.surrogates.poly import expressions

CELL_VOLUME = 4e-15  # [L cell^-1] (paper, ADR 0004)
GrowthMode = Literal["steady", "dynamic"]


def _load(name: str) -> dict:
    return json.loads((DATA / "surrogates" / name).read_text())


def relaxation(x: float, x_star: float, mu: float) -> float:
    return mu * (x_star - x)


def growth(n: float, mu: float) -> float:
    return mu * n


def scaled(x: float, factor: float) -> float:
    return x * factor


def per_nl(n: float) -> float:
    """[cells mL^-1] -> [cells nL^-1] (unit of the legacy light model)."""
    return n * 1e-6


def get_model(
    inp: LegacyInputs,
    growth_mode: GrowthMode = "dynamic",
    petc_poly: dict | None = None,
    metabolism_poly: dict | None = None,
    depth_nodes: int = 5,
) -> Model:
    """Growth model; `depth_nodes` = 0 evaluates the PETC at the depth-averaged light."""
    petc_poly = _load("petc_poly.json") if petc_poly is None else petc_poly
    metabolism_poly = (
        _load("metabolism_poly.json") if metabolism_poly is None else metabolism_poly
    )
    m = Model()
    m.add_variables(
        {
            "N": 5.456287e7,  # first density of the batch data
            "Chl": REFERENCE.chl,
            "Car": REFERENCE.carotenoids,
            "SQ": 0.0,
        }
    )
    m.add_parameters(
        {
            "light": 80.0,  # [µmol(photons) m^-2 s^-1] MC-1000 batch (paper Fig. 9)
            "sample_depth_m": 0.0135,  # half the culture tube diameter (paper Fig. 9)
            "cell_volume": CELL_VOLUME,
            # legacy per-cell phycobiliprotein amounts (53180 / 24076 mg L^-1 at 2 fL,
            # source unknown, docs/assumptions.md) times a fitted factor
            "phycocyanin_ref": 53180 * 2e-15 / CELL_VOLUME,  # [mg L^-1]
            "allophycocyanin_ref": 24076 * 2e-15 / CELL_VOLUME,  # [mg L^-1]
            "pbs_scale": 1.0,  # [-] fitted
            "ps_ratio": 5.9,
            "beta_carotene_fraction": 0.26,
            "petc_scale": 1.0,  # [-] fitted: efficiency of the PETC supply
            "S_PGA_floor": 1e-3,  # [mM min^-1] net respiration at very low light
            "E_CrtE": 1.0,  # [-] CrtE level relative to WT (Fig 11 knock-down)
        }
    )
    m.add_derived("N_nL", fn=per_nl, args=["N"])
    m.add_derived("phycocyanin", fn=scaled, args=["phycocyanin_ref", "pbs_scale"])
    m.add_derived(
        "allophycocyanin", fn=scaled, args=["allophycocyanin_ref", "pbs_scale"]
    )

    symbols = Symbols(
        chl="Chl",
        carotenoids="Car",
        cell_density="N_nL",
        light="light",
        phycocyanin="phycocyanin",
        allophycocyanin="allophycocyanin",
        ps_ratio="ps_ratio",
        beta_carotene_fraction="beta_carotene_fraction",
        sample_depth="sample_depth_m",
        cell_volume="cell_volume",
    )
    names = {
        "absorption_ps1": "complex_abs_ps1",
        "absorption_ps2": "complex_abs_ps2",
        "absorption_pbs": "complex_abs_pbs",
        "light_OCP": "light_ocp",
    }
    exprs: dict[str, sp.Expr] = {}
    if depth_nodes == 0:
        # PETC at the depth-averaged light (legacy)
        exprs |= light_expressions(inp, symbols)
        exprs |= expressions(
            petc_poly, {i: sp.Symbol(names[i]) for i in petc_poly["inputs"]}
        )
    else:
        # PETC at the local light of Gauss-Legendre depth nodes, outputs averaged:
        # the PETC is nonlinear, so the mean of the rates != the rate at the mean light
        nodes, weights = np.polynomial.legendre.leggauss(depth_nodes)
        j = {"J_3PGA": [], "J_ATP": []}
        for i, (x, w) in enumerate(zip(nodes, weights, strict=True)):
            pre = f"z{i}_"
            exprs[f"{pre}depth"] = float((x + 1) / 2) * sp.Symbol("sample_depth_m")
            exprs |= light_expressions(inp, symbols, depth=f"{pre}depth", prefix=pre)
            local = expressions(
                petc_poly, {k: sp.Symbol(pre + v) for k, v in names.items()}
            )
            for out, expr in local.items():
                exprs[pre + out] = expr
                j[out].append(float(w / 2) * sp.Symbol(pre + out))
        exprs |= {out: sp.Add(*terms) for out, terms in j.items()}
    chl, scale = sp.Symbol("Chl"), sp.Symbol("petc_scale")
    exprs["S_PGA"] = sp.Max(
        sp.Symbol("J_3PGA") * chl * 0.06 * scale, sp.Symbol("S_PGA_floor")
    )
    exprs["S_ATP"] = sp.Max(sp.Symbol("J_ATP"), 0) * chl * 0.06 * scale
    exprs |= expressions(
        metabolism_poly, {i: sp.Symbol(i) for i in metabolism_poly["inputs"]}
    )
    for name, expr in exprs.items():
        fn, args = expr_fn(expr, origin=name)
        m.add_derived(name, fn=fn, args=args)

    if growth_mode == "dynamic":
        m.add_reaction("growth", fn=growth, args=["N", "mu"], stoichiometry={"N": 1})
    # In steady mode (fixed N) only the fixed point X = X*(X) matters, which doesn't
    # depend on the relaxation rate; relaxing at mu (~1e-5 min^-1 in dense cultures)
    # makes steady-state searches fail, so a fast constant rate is used instead.
    if growth_mode == "steady":
        m.add_parameter("k_relax", 1.0)  # [min^-1]
    rate = "mu" if growth_mode == "dynamic" else "k_relax"
    for state in ("Chl", "Car", "SQ"):
        m.add_reaction(
            f"net_{state}",
            fn=relaxation,
            args=[state, f"{state}_star", rate],
            stoichiometry={state: 1},
        )
    return m
