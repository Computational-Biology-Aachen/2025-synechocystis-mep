"""Growth layer (hours timescale) of the `double` family (ADR 0004, ADR 0006).

Slow states: cell density `N` [cells mL^-1], chlorophyll `Chl`, carotenoids `Car` and
squalene `SQ` [mM, intracellular]. Faster layers enter as QSS surrogates:

    N, Chl, Car --self-shading--> absorbed light per complex at 5 depths
                --PETC surrogate, averaged over depth--> J_3PGA, J_ATP
                                                  [mmol mol(Chl)^-1 s^-1]
                × Chl × petc_scale --> S_PGA, S_ATP [mM min^-1]
                --metabolism surrogate (+ enzyme levels E_*)--> mu, X* (contents)

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

from synmep import acclimation
from synmep import light as light_sources
from synmep.data import DATA, LegacyInputs
from synmep.exprfn import expr_fn
from synmep.models.metabolism import REFERENCE, reference_parameters
from synmep.surrogates.petc_legacy import Symbols, light_expressions
from synmep.surrogates.poly import expressions

CELL_VOLUME = 4e-15  # [L cell^-1] (paper, ADR 0004)
GrowthMode = Literal["steady", "dynamic"]
# ADR 0004: "coupled" = pigment-branch enzymes scale with growth (separate metabolism
# surrogate). The fixed-pigment variants (photoacclimation, capacity) are retired
# (ADR 0006).
PigmentRegulation = Literal["none", "coupled"]
L_REF = acclimation.L_REF
# PETC surrogate trained on white light + LED mixtures (report 08); the metabolism
# surrogate is per strain (`synmep.surrogates.metabolism_strain`, ADR 0006)
PETC_POLY = "petc_poly_led.json"


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


def petc_supply_expressions(
    inp: LegacyInputs,
    petc_poly: dict,
    depth_nodes: int,
    spectral_step: int,
    light_channels: tuple[str, ...] | None = None,
) -> dict[str, sp.Expr]:
    """Light model + PETC surrogate -> J_3PGA, J_ATP and the mean PS1 absorption.

    Expects model names Chl, Car, N_nL, light, phycocyanin, allophycocyanin, ps_ratio,
    beta_carotene_fraction, sample_depth_m, cell_volume.
    """
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
    ch = None if light_channels is None else light_sources.channels(light_channels)
    if depth_nodes == 0:
        # PETC at the depth-averaged light (legacy)
        exprs |= light_expressions(inp, symbols, step=spectral_step, channels=ch)
        exprs |= expressions(
            petc_poly, {i: sp.Symbol(names[i]) for i in petc_poly["inputs"]}
        )
    else:
        # PETC at the local light of Gauss-Legendre depth nodes, outputs averaged:
        # the PETC is nonlinear, so the mean of the rates != the rate at the mean light
        nodes, weights = np.polynomial.legendre.leggauss(depth_nodes)
        i_ps1 = petc_poly["inputs"].index("absorption_ps1")
        a1_min = float(np.exp(petc_poly["log_lo"][i_ps1]))
        j = {"J_3PGA": [], "J_ATP": []}
        for i, (x, w) in enumerate(zip(nodes, weights, strict=True)):
            pre = f"z{i}_"
            exprs[f"{pre}depth"] = float((x + 1) / 2) * sp.Symbol("sample_depth_m")
            exprs |= light_expressions(
                inp,
                symbols,
                depth=f"{pre}depth",
                prefix=pre,
                step=spectral_step,
                channels=ch,
            )
            # Below the training domain (dark culture layers) all inputs are scaled up
            # together to its edge (ratios stay on the training manifold) and the
            # outputs are interpolated linearly to zero supply in darkness, instead of
            # clamping each input separately (off-manifold, spurious supply).
            a_ps1 = sp.Symbol(pre + names["absorption_ps1"])
            exprs[pre + "light_rel"] = a_ps1 / a1_min
            rel = sp.Symbol(pre + "light_rel")
            up = sp.Max(1, 1 / rel)
            local = expressions(
                petc_poly, {k: sp.Symbol(pre + v) * up for k, v in names.items()}
            )
            for out, expr in local.items():
                exprs[pre + out] = expr * sp.Min(1, rel)
                j[out].append(float(w / 2) * sp.Symbol(pre + out))
        exprs |= {out: sp.Add(*terms) for out, terms in j.items()}
        for c in ("ps1", "ps2", "pbs"):
            exprs[f"complex_abs_{c}_mean"] = sp.Add(
                *(
                    float(w / 2) * sp.Symbol(f"z{i}_complex_abs_{c}")
                    for i, w in enumerate(weights)
                )
            )
    if depth_nodes == 0:
        for c in ("ps1", "ps2", "pbs"):
            exprs[f"complex_abs_{c}_mean"] = sp.Symbol(f"complex_abs_{c}")
    return exprs


def get_model(
    inp: LegacyInputs,
    growth_mode: GrowthMode = "dynamic",
    petc_poly: dict | None = None,
    metabolism_poly: dict | None = None,
    depth_nodes: int = 5,
    pigment_regulation: PigmentRegulation = "none",
    spectral_step: int = 5,
    light_channels: tuple[str, ...] | None = None,
    capacity: bool = False,
    folds: dict[str, float] | None = None,
    chromatic: bool = False,
) -> Model:
    """Growth model.

    `depth_nodes` = 0 evaluates the PETC at the depth-averaged light. `light_channels`
    (`synmep.light` channel names) makes the spectrum a mixture with weights `w_<name>`.
    `capacity` adds capacity acclimation (`synmep.acclimation`). `folds` (enzyme
    parameter -> fold vs WT) selects the strain's metabolism surrogate (trained on
    first use). `spectral_step`
    [nm] coarsens the wavelength grid of the light model: 5 nm changes the absorbed
    light by a median 3e-4 (max 6.5 %, only in very dark culture layers) and makes the
    model ~5x smaller than at 1 nm (ADR 0004).
    """
    petc_poly = _load(PETC_POLY) if petc_poly is None else petc_poly
    if metabolism_poly is None:
        from synmep.surrogates.metabolism_strain import load_or_train

        coupling = 1.0 if pigment_regulation == "coupled" else 0.0
        metabolism_poly = load_or_train(dict(folds or {}), coupling)
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
            "atp_scale": reference_parameters(REFERENCE)["atp_scale"],
        }
        # enzyme levels relative to WT that the metabolism surrogate takes
        | {i: 1.0 for i in metabolism_poly["inputs"] if i.startswith("E_")}
    )
    m.add_derived("N_nL", fn=per_nl, args=["N"])
    m.add_derived("phycocyanin", fn=scaled, args=["phycocyanin_ref", "pbs_scale"])
    m.add_derived(
        "allophycocyanin", fn=scaled, args=["allophycocyanin_ref", "pbs_scale"]
    )

    if light_channels is not None:
        light_sources.add_parameters(m, light_channels)
    exprs = petc_supply_expressions(
        inp, petc_poly, depth_nodes, spectral_step, light_channels
    )
    chl, scale = sp.Symbol("Chl"), sp.Symbol("petc_scale")
    if capacity:
        m.add_parameters(acclimation.PARAMETERS)
        exprs |= acclimation.target_expressions()
        scale = scale * sp.Symbol("Cap")
    if chromatic:
        m.add_parameters(acclimation.CHROMATIC)
        exprs |= acclimation.chromatic_expressions()
    exprs["S_PGA"] = sp.Max(
        sp.Symbol("J_3PGA") * chl * 0.06 * scale, sp.Symbol("S_PGA_floor")
    )
    exprs["S_ATP"] = (
        sp.Max(sp.Symbol("J_ATP"), 0) * chl * 0.06 * scale * sp.Symbol("atp_scale")
    )
    exprs |= expressions(
        metabolism_poly, {i: sp.Symbol(i) for i in metabolism_poly["inputs"]}
    )
    targets = {s: sp.Symbol(f"{s}_star") for s in ("Chl", "Car", "SQ")}
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
    if capacity:
        acclimation.add_state(m, None if growth_mode == "dynamic" else "k_relax")
    if chromatic:
        acclimation.add_chromatic_state(
            m, None if growth_mode == "dynamic" else "k_relax"
        )
    for state in ("Chl", "Car", "SQ"):
        m.add_reaction(
            f"net_{state}",
            fn=relaxation,
            args=[state, str(targets[state]), rate],
            stoichiometry={state: 1},
        )
    return m
