"""Minimal model (v3, ADR 0005): culture squalene productivity with the pigment trade-off.

One ODE model, built from thread g (ADR 0004) with three reductions:

- PETC: closed form instead of a polynomial surrogate (`scripts/fit_petc_closed_form.py`)

      I = a_PS1 + w_PS2 a_PS2 + w_PBS a_PBS,   J_3PGA = J_max tanh(I / I_k) - r,
      J_ATP = c0 + c1 J_3PGA + c2 J_3PGA^2

  evaluated at Gauss-Legendre depths of the self-shaded light profile and averaged.
- metabolism: the thread-g metabolism ODE (`g_metabolism`) simulated directly, with
  growth-coupled pigment-branch enzymes (pigment regulation A), no surrogate.
- slow pools: dX/dt = v_X - mu X for Chl, carotenoids and squalene; dN/dt = mu N.

Additions:

- **enzyme cost** of overexpression. Delta = sum_i (E_i - 1) f_i is the extra protein
  as a fraction of the proteome (f_i: WT proteome mass fraction of enzyme i from the
  proteomics). `cost="drain"`: the extra protein is synthesised at mu Delta P_cell
  [g L^-1 min^-1] and drains 3PGA, ATP and NADPH like biomass. `cost="proteome"`: the
  extra protein displaces growth and photosynthetic capacity, which scale with
  (1 - phi_Q) / (1 - phi_Q + Delta) (phi_Q: fixed housekeeping share).
- optional **capacity acclimation** (`capacity=True`): photosynthetic capacity Cap
  relaxes over hours to (1 + K) L_rel / (K + L_rel).

Dry weight per cell volume from the measured protein (ADR 0005 §2.5): the gravimetric
dry weight of the proteomics samples (3.88 pg cell^-1, 0.97 g mL^-1 at 4 fL) contains
only 22 % protein, so it's taken to be inflated; protein per cell (0.83 pg) / 0.445 gives
1.85 pg DW cell^-1, 0.46 g mL^-1.

`growth_mode="steady"`: cell density N is a parameter (turbidostat) and the slow
equations are multiplied by `slow_speedup`, which leaves the steady state unchanged
and lets steady-state searches converge.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Literal

import numpy as np
import sympy as sp
from mxlpy import Model

from synmep.data import CELL_VOLUME as LEGACY_CELL_VOLUME
from synmep.data import DATA, PROTEIN_CONTENT, LegacyInputs, proteomics_samples
from synmep.exprfn import expr_fn
from synmep.models import g_metabolism
from synmep.models.g_growth import CELL_VOLUME, L_REF, growth, per_nl, scaled
from synmep.surrogates.petc_legacy import Symbols, light_expressions

GrowthMode = Literal["steady", "dynamic"]
Cost = Literal["none", "drain", "proteome"]
FIT = Path("docs/reports")

# WT proteome mass fractions of the lumped enzymes (see `proteome_fractions`)
LUMPS = {
    "DXS": ("DXS",),
    "MEP": ("DXR", "IspD", "IspE", "IspF", "IspG", "IspH"),
    "CrtE": ("CrtE",),
    "Sqs": ("Sqs",),
    "CrtB": ("CrtB",),
    "Chl": ("ChlP", "ChlG"),
}
VARIANTS = ("v3-drain", "v3-proteome", "v3-drain-cap", "v3-proteome-cap")
PHI_Q = 0.45  # [-] growth-independent proteome share (Scott2010, E. coli; assumption)


def variant(name: str) -> dict:
    """'v3-<cost>[-cap]' -> `get_model` keyword arguments."""
    parts = name.split("-")
    return {"cost": parts[1], "capacity": parts[-1] == "cap"}


def dry_weight_density(cell_volume: float = CELL_VOLUME) -> float:
    """[gDW L(cell)^-1] from the median measured protein per cell / protein content."""
    protein_pg = float(proteomics_samples()["protein_pg"].median())
    return protein_pg * 1e-12 / PROTEIN_CONTENT / cell_volume


def proteome_fractions(inp: LegacyInputs) -> dict[str, float]:
    """WT proteome mass fraction of each lump, f = E0 M / protein concentration.

    E0 [µM] are the legacy proteomics concentrations (gravimetric dry weight x 0.445,
    per legacy cell volume, as loaded by `load_legacy_inputs()`); dividing by the protein concentration on the same basis
    gives the (basis-free) mass fraction.
    """
    samples = proteomics_samples()
    wt = samples[samples["sample"].str.startswith("WT")]
    protein_basis = (  # [g protein L(cell)^-1] as used for E0
        float(wt["dry_weight_pg"].median())
        * 1e-12
        * PROTEIN_CONTENT
        / LEGACY_CELL_VOLUME
    )
    poi = inp.proteins_of_interest
    out = {}
    for lump, members in LUMPS.items():
        out[lump] = sum(
            inp.enzyme_concentrations_umol[e]
            * 1e-6
            * float(poi[poi["Enzyme"] == e]["kDa"].iloc[0])
            * 1e3
            / protein_basis
            for e in members
        )
    return out


def petc_closed_form() -> dict[str, float]:
    return json.loads((DATA / "surrogates" / "petc_closed_form.json").read_text())


def fitted_parameters(capacity: bool) -> dict[str, float]:
    path = FIT / f"v3-fit{'-capacity' if capacity else ''}.json"
    return json.loads(path.read_text())["fitted"] if path.exists() else {}


def flux_minus_dilution(v: float, x: float, mu: float, speed: float) -> float:
    return speed * (v - mu * x)


def relax(x: float, target: float, rate: float) -> float:
    return rate * (target - x)


def protein_drain(mu: float, delta: float, rho_dw: float, protein: float) -> float:
    """[g protein L^-1 min^-1] synthesis of the extra (overexpressed) protein."""
    return mu * delta * rho_dw * protein


def _tanh(x: sp.Expr) -> sp.Expr:
    # x >= 0 (absorbed light): written with exp(-2x), lowers natively, no overflow
    e = sp.exp(-2 * x)
    return (1 - e) / (1 + e)


def petc_expressions(
    inp: LegacyInputs, depth_nodes: int, spectral_step: int
) -> dict[str, sp.Expr]:
    """Light model + closed-form PETC -> J_3PGA, J_ATP [mmol mol(Chl)^-1 s^-1] (depth
    means) and the mean PS1 absorption `complex_abs_ps1_mean`."""
    pc = petc_closed_form()
    symbols = Symbols(
        chl="Chl", carotenoids="Car", cell_density="N_nL", light="light",
        phycocyanin="phycocyanin", allophycocyanin="allophycocyanin",
        ps_ratio="ps_ratio", beta_carotene_fraction="beta_carotene_fraction",
        sample_depth="sample_depth_m", cell_volume="cell_volume",
    )  # fmt: skip
    nodes, weights = np.polynomial.legendre.leggauss(depth_nodes)
    exprs: dict[str, sp.Expr] = {}
    means: dict[str, list[sp.Expr]] = {"J_3PGA": [], "J_ATP": [], "abs_ps1": []}
    for i, (x, w) in enumerate(zip(nodes, weights, strict=True)):
        pre = f"z{i}_"
        exprs[f"{pre}depth"] = float((x + 1) / 2) * sp.Symbol("sample_depth_m")
        exprs |= light_expressions(
            inp, symbols, depth=f"{pre}depth", prefix=pre, step=spectral_step
        )
        a1, a2, apbs = (
            sp.Symbol(f"{pre}complex_abs_{c}") for c in ("ps1", "ps2", "pbs")
        )
        exprs[f"{pre}I_abs"] = a1 + pc["w_ps2"] * a2 + pc["w_pbs"] * apbs
        j3 = sp.Symbol(f"{pre}J_3PGA")
        exprs[f"{pre}J_3PGA"] = (
            pc["J_max"] * _tanh(sp.Symbol(f"{pre}I_abs") / pc["I_k"]) - pc["r"]
        )
        exprs[f"{pre}J_ATP"] = pc["c0"] + pc["c1"] * j3 + pc["c2"] * j3**2
        for k, name in (
            ("J_3PGA", "J_3PGA"),
            ("J_ATP", "J_ATP"),
            ("abs_ps1", "complex_abs_ps1"),
        ):
            means[k].append(float(w / 2) * sp.Symbol(f"{pre}{name}"))
    exprs["J_3PGA"] = sp.Add(*means["J_3PGA"])
    exprs["J_ATP"] = sp.Add(*means["J_ATP"])
    exprs["complex_abs_ps1_mean"] = sp.Add(*means["abs_ps1"])
    return exprs


def get_model(
    inp: LegacyInputs,
    *,
    growth_mode: GrowthMode = "steady",
    cost: Cost = "proteome",
    capacity: bool = False,
    cell_density: float = 1e7,  # [cells mL^-1] (steady: parameter N; dynamic: N(0))
    light: float = 100.0,  # [µmol(photons) m^-2 s^-1]
    sample_depth: float = 0.01,  # [m]
    depth_nodes: int = 5,
    spectral_step: int = 5,
    use_fit: bool = True,
) -> Model:
    rho_dw = dry_weight_density()
    ref = g_metabolism.Reference(rho_dw=rho_dw)
    m = g_metabolism.get_model(ref, pigment_coupling=1.0)
    v_bio_wt = m.get_parameter_values()["V_bio"]
    m.remove_parameters(["S_PGA", "S_ATP"])
    m.add_variables({"Chl": ref.chl, "Car": ref.carotenoids, "SQ": ref.squalene})
    dynamic = growth_mode == "dynamic"
    if dynamic:
        m.add_variable("N", cell_density)
    m.add_parameters(
        ({} if dynamic else {"N": cell_density})
        | {
            "light": light,
            "sample_depth_m": sample_depth,
            "cell_volume": CELL_VOLUME,
            # legacy per-cell phycobiliprotein amounts (thread g), times a fitted factor
            "phycocyanin_ref": 53180 * 2e-15 / CELL_VOLUME,  # [mg L^-1]
            "allophycocyanin_ref": 24076 * 2e-15 / CELL_VOLUME,  # [mg L^-1]
            "pbs_scale": 1.0,  # [-] fitted
            "ps_ratio": 5.9,
            "beta_carotene_fraction": 0.26,
            "petc_scale": 1.0,  # [-] fitted: efficiency of the PETC supply
            "S_PGA_floor": 1e-3,  # [mM min^-1] net respiration at very low light
            "slow_speedup": 1.0 if dynamic else 100.0,  # [-] steady mode only
            "protein_content": PROTEIN_CONTENT,  # [g protein gDW^-1]
            "phi_Q": PHI_Q,
        }
        | {f"f_{k}": v for k, v in proteome_fractions(inp).items()}
    )
    if capacity:
        m.add_variable("Cap", 1.0)
        m.add_parameters({"L_ref": L_REF, "capacity_K": 1.0, "capacity_tau": 600.0})
    if use_fit:
        fit = fitted_parameters(capacity)
        m.update_parameters(
            {k: v for k, v in fit.items() if k in m.get_parameter_names()}
        )

    m.add_derived("N_nL", fn=per_nl, args=["N"])
    m.add_derived("phycocyanin", fn=scaled, args=["phycocyanin_ref", "pbs_scale"])
    m.add_derived(
        "allophycocyanin", fn=scaled, args=["allophycocyanin_ref", "pbs_scale"]
    )

    exprs = petc_expressions(inp, depth_nodes, spectral_step)
    delta = sp.Add(*((sp.Symbol(f"E_{k}") - 1) * sp.Symbol(f"f_{k}") for k in LUMPS))
    exprs["protein_extra"] = delta  # [-] extra protein / proteome
    supply_scale = sp.Symbol("petc_scale")
    if cost == "proteome":
        phi_q = sp.Symbol("phi_Q")
        exprs["capacity_share"] = (1 - phi_q) / (1 - phi_q + sp.Symbol("protein_extra"))
        supply_scale = supply_scale * sp.Symbol("capacity_share")
    if capacity:
        exprs["L_rel"] = sp.Symbol("complex_abs_ps1_mean") / sp.Symbol("L_ref")
        k_cap, l_rel = sp.Symbol("capacity_K"), sp.Symbol("L_rel")
        exprs["Cap_target"] = (1 + k_cap) * l_rel / (k_cap + l_rel)
        supply_scale = supply_scale * sp.Symbol("Cap")
    chl = sp.Symbol("Chl")
    exprs["S_PGA"] = sp.Max(
        sp.Symbol("J_3PGA") * chl * 0.06 * supply_scale, sp.Symbol("S_PGA_floor")
    )
    exprs["S_ATP"] = sp.Max(sp.Symbol("J_ATP"), 0) * chl * 0.06 * supply_scale
    if cost == "proteome":  # growth capacity is displaced too
        m.add_parameter("V_bio_WT", v_bio_wt)
        m.remove_parameter("V_bio")
        exprs["V_bio"] = sp.Symbol("V_bio_WT") * sp.Symbol("capacity_share")
    for name, expr in exprs.items():
        fn, args = expr_fn(expr, origin=name)
        m.add_derived(name, fn=fn, args=args)

    if cost == "drain":
        m.add_reaction(
            "v_protein_extra",
            fn=protein_drain,
            args=["mu", "protein_extra", "rho_DW", "protein_content"],
            stoichiometry={s: -c for s, c in g_metabolism.BIOMASS.items()},
        )
    if dynamic:
        m.add_reaction("growth", fn=growth, args=["N", "mu"], stoichiometry={"N": 1})
    for state, flux in (("Chl", "v_chl"), ("Car", "v_car"), ("SQ", "v_sq")):
        m.add_reaction(
            f"net_{state}",
            fn=flux_minus_dilution,
            args=[flux, state, "mu", "slow_speedup"],
            stoichiometry={state: 1},
        )
    if capacity:
        m.add_derived("capacity_rate", fn=_rate, args=["capacity_tau", "slow_speedup"])
        m.add_reaction(
            "net_Cap",
            fn=relax,
            args=["Cap", "Cap_target", "capacity_rate"],
            stoichiometry={"Cap": 1},
        )
    return m


def _rate(tau: float, speedup: float) -> float:
    return speedup / tau
