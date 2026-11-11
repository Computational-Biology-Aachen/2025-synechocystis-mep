"""`full` family: pfennig2024 included directly, no surrogates (ADR 0006 §2.5).

Model time is seconds (pfennig2024's unit); `mu` is reported in min^-1 like the other
families.

Structure:

- **Depth layers** i = 0..K-1 at Gauss–Legendre depths of the sample (weights w_i).
  Each layer has its local light (`light_expressions` with depth) and a full copy of
  pfennig2024 (17 states, prefix `z<i>_`), with absorbed light per complex and OCP
  light from the light model, 3PGA free (not clamped).
- **Fast pools per layer:** pfennig2024's pools, among them 3PGA, ATP/ADP and
  NADPH/NADP+, plus FPP and GGPP (relaxation < 1 min, Jacobian of the metabolism).
  **Shared, slow:** MEPi and IPP (3 and 12 min), Chl, carotenoids, squalene, N.
- **Metabolism** (`synmep.models.metabolism`) runs per layer: every reaction reads its
  layer's fast pools and the shared slow pools. It acts fully on its layer's pools and
  on shared pools with weight w_i. Concentrations in mM are converted from pfennig's
  mmol mol(Chl)^-1 with the current Chl, rates from mM min^-1 to mmol mol(Chl)^-1 s^-1.
- **Changes to pfennig2024:** the lumped ATP consumption is removed (the metabolism's
  biomass, MEP, maintenance demand replace it). The metabolism's maximal rates and
  Michaelis constants are re-derived (`metabolism.Reference`) at pfennig2024's ATP and
  NADPH levels. pfennig2024's soluble totals (adenylates, NADP, NAD) and Michaelis
  constants are per Chl. `soluble_per_cell=True` makes the totals per-cell constants
  instead; without also converting the per-Chl Michaelis constants that creates a
  positive feedback (more Chl -> less NADP per Chl -> slower CBB -> slower growth ->
  more Chl), so the default keeps pfennig2024's per-Chl totals.
- mu = sum_i w_i v_bio,i / rho_DW; pigment-branch enzymes scale with growth (as
  `minimal`, exponent `pigment_coupling`); dry weight from protein (ADR 0005).
"""

from __future__ import annotations

import dataclasses
from typing import Literal

import numpy as np
import sympy as sp
from mxlpy import Derived, Model

from synmep import light as light_sources
from synmep.data import LegacyInputs
from synmep.exprfn import expr_fn
from synmep.models import metabolism
from synmep.models.double_growth import CELL_VOLUME, per_nl, scaled
from synmep.models.minimal import dry_weight_density
from synmep.petc.pfennig import LIGHT_INPUTS, get_petc_model
from synmep.surrogates.petc_legacy import Symbols, light_expressions

GrowthMode = Literal["steady", "dynamic"]
PER_LAYER_METABOLITES = ("FPP", "GGPP")
SHARED_METABOLITES = ("MEPi", "IPP")
PFENNIG_POOLS = {"PGA": "3PGA", "ATP": "ATP", "NADPH": "NADPH"}  # metabolism -> pf
SOLUBLE_TOTALS = ("AP_tot", "NADP_tot", "NAD_tot")  # [mmol mol(Chl)^-1] in pfennig
REFERENCE_CHL = 19.9  # [mM] pfennig2024's per-Chl totals hold at this Chl content
BIOMASS_ATP = metabolism.BIOMASS_ATP  # [mmol gDW^-1] (ADR 0006 checkpoint B)
# non-growth-associated ATP maintenance [mmol ATP gDW^-1 h^-1], lower bound from protein
# turnover in an iJN678-based Synechocystis model (arXiv 2506.05916). Saturating in ATP
# (KM_NGAM), so it stops when ATP is gone (ADR 0006 checkpoint B). Without an ATP sink
# that doesn't need NADPH, dim layers lock up: ATP fully charged, ATP synthase stalls,
# lumen pH ~3, photosynthetic control stops electron transport.
NGAM = 0.13
KM_NGAM = 0.1  # [mM] (assumption, well below the ~8.5 mM adenylate pool)
# biomass synthesis' NADPH half-saturation x10 vs the metabolism reference (0.149 ->
# 1.49 mM, NADP(H) pool 0.53 mM): with the reference value biomass synthesis drains
# NADPH in dim layers until it's nearly gone, then stops, a relaxation oscillation
# against FNR (~4 s period) that made batches take 5-7 min. x3 still oscillates.
# Lowers growth at given supply; the fitted efficiency compensates (report 09).
KM_BIO_NADPH_FACTOR = 10.0


def per_chl(total_cell: float, chl: float) -> float:
    """[mmol mol(Chl)^-1] from a per-cell concentration [mM]."""
    return total_cell * 1000 / chl


def to_mM(x_per_chl: float, chl: float) -> float:
    return x_per_chl * chl / 1000


def ngam_rate(atp: float, ngam: float, rho_dw: float, km: float) -> float:
    """[mM min^-1] maintenance ATP use, ngam [mmol gDW^-1 h^-1] x rho_DW, saturating."""
    return ngam / 60 * rho_dw * atp / (km + atp)


def slow_dilution(x: float, mu: float, speed: float) -> float:
    """[mM s^-1] growth dilution of a slow pool (mu in min^-1), sped up in steady mode."""
    return speed * mu * x / 60


def growth(n: float, mu: float) -> float:
    return mu * n / 60


def _rename(names: list[str], table: dict[str, str]) -> list[str]:
    return [table.get(n, n) for n in names]


def _copy_pfennig(m: Model, pf: Model, pre: str, inputs: dict[str, str]) -> None:
    """Add a prefixed copy of pfennig2024's states, derived values, surrogates and
    reactions; parameters stay shared, the light inputs map to `inputs`."""
    local = set(pf.get_variable_names()) | set(pf.get_raw_derived())
    local |= {o for s in pf.get_raw_surrogates().values() for o in s.outputs}
    table = {n: pre + n for n in local} | inputs
    for name, var in pf.get_raw_variables().items():
        m.add_variable(pre + name, var.initial_value)
    for name, d in pf.get_raw_derived().items():
        m.add_derived(pre + name, fn=d.fn, args=_rename(d.args, table))
    for name, s in pf.get_raw_surrogates().items():
        m.add_surrogate(
            pre + name,
            dataclasses.replace(
                s,
                args=_rename(s.args, table),
                outputs=_rename(s.outputs, table),
                units={table[k]: v for k, v in s.units.items()},
                stoichiometries={
                    table.get(k, k): {table.get(a, a): b for a, b in v.items()}
                    for k, v in s.stoichiometries.items()
                },
            ),
        )
    for name, r in pf.get_raw_reactions().items():
        if name == "vATPconsumption":
            continue  # replaced by the metabolism's ATP demand
        stoich = {}
        for k, v in r.stoichiometry.items():
            if isinstance(v, Derived):
                v = Derived(fn=v.fn, args=_rename(v.args, table))
            stoich[table.get(k, k)] = v
        m.add_reaction(
            pre + name, fn=r.fn, args=_rename(r.args, table), stoichiometry=stoich
        )


def _stoich(expr: sp.Expr, origin: str) -> Derived | float:
    if not expr.free_symbols:
        return float(expr)
    fn, args = expr_fn(expr, origin)
    return Derived(fn=fn, args=args)


def get_model(
    inp: LegacyInputs,
    *,
    growth_mode: GrowthMode = "steady",
    cell_density: float = 1e7,  # [cells mL^-1]
    light: float = 100.0,  # [µmol(photons) m^-2 s^-1]
    sample_depth: float = 0.01,  # [m]
    layers: int = 3,
    spectral_step: int = 5,
    light_channels: tuple[str, ...] | None = None,
    pigment_coupling: float = 1.0,
    soluble_per_cell: bool = False,
    biomass_atp: float = BIOMASS_ATP,
    km_bio_nadph_factor: float = KM_BIO_NADPH_FACTOR,
) -> Model:
    dynamic = growth_mode == "dynamic"
    pf = get_petc_model(clamp_pga=False)
    pf_ss = _pfennig_reference()
    rho_dw = dry_weight_density()
    ref = metabolism.Reference(
        rho_dw=rho_dw,
        pools=metabolism.REFERENCE.pools
        | {
            "ATP": to_mM(pf_ss["ATP"], REFERENCE_CHL),
            "NADPH": to_mM(pf_ss["NADPH"], REFERENCE_CHL),
        },
    )
    met = metabolism.get_model(
        ref,
        pigment_coupling=pigment_coupling,
        biomass=metabolism.BIOMASS | {"ATP": biomass_atp},
    )
    nodes, weights = np.polynomial.legendre.leggauss(layers)
    w = weights / 2

    m = Model()
    pf_pars = pf.get_parameter_values()
    totals_cell = (
        {f"{t}_cell": to_mM(pf_pars[t], REFERENCE_CHL) for t in SOLUBLE_TOTALS}
        if soluble_per_cell
        else {}
    )
    skip = set(LIGHT_INPUTS) | (set(SOLUBLE_TOTALS) if soluble_per_cell else set())
    m.add_parameters({k: v for k, v in pf_pars.items() if k not in skip})
    met_pars = met.get_parameter_values()
    # ATP demand is mechanistic here (biomass, MEP chain); the surrogate families'
    # maintenance sink balanced an assumed supply and would drain pfennig2024's
    # nearly fully charged adenylate pool
    met_pars["k_ATP_maintenance"] = 0.0
    met_pars["Km_bio_NADPH"] *= km_bio_nadph_factor
    skip_met = ("S_PGA", "S_ATP", "r_NADPH", "atp_scale")
    m.add_parameters({k: v for k, v in met_pars.items() if k not in skip_met})
    m.add_parameters(
        totals_cell
        | ({} if dynamic else {"N": cell_density})
        | {
            "light": light,
            "sample_depth_m": sample_depth,
            "cell_volume": CELL_VOLUME,
            "phycocyanin_ref": 53180 * 2e-15 / CELL_VOLUME,
            "allophycocyanin_ref": 24076 * 2e-15 / CELL_VOLUME,
            "pbs_scale": 1.0,
            "ps_ratio": 5.9,
            "beta_carotene_fraction": 0.26,
            "slow_speedup": 1.0 if dynamic else 100.0,
            "NGAM": NGAM,
            "KM_NGAM": KM_NGAM,
        }
    )
    if light_channels is not None:
        light_sources.add_parameters(m, light_channels)
    ch = None if light_channels is None else light_sources.channels(light_channels)

    m.add_variables({"Chl": ref.chl, "Car": ref.carotenoids, "SQ": ref.squalene})
    m.add_variables({k: ref.pools[k] for k in SHARED_METABOLITES})
    if dynamic:
        m.add_variable("N", cell_density)
    m.add_derived("N_nL", fn=per_nl, args=["N"])
    m.add_derived("phycocyanin", fn=scaled, args=["phycocyanin_ref", "pbs_scale"])
    m.add_derived(
        "allophycocyanin", fn=scaled, args=["allophycocyanin_ref", "pbs_scale"]
    )
    if soluble_per_cell:
        for t in SOLUBLE_TOTALS:
            m.add_derived(t, fn=per_chl, args=[f"{t}_cell", "Chl"])

    symbols = Symbols(
        chl="Chl", carotenoids="Car", cell_density="N_nL", light="light",
        phycocyanin="phycocyanin", allophycocyanin="allophycocyanin",
        ps_ratio="ps_ratio", beta_carotene_fraction="beta_carotene_fraction",
        sample_depth="sample_depth_m", cell_volume="cell_volume",
    )  # fmt: skip
    exprs: dict[str, sp.Expr] = {}
    for i, x in enumerate(nodes):
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
    for name, expr in exprs.items():
        fn, args = expr_fn(expr, origin=name)
        m.add_derived(name, fn=fn, args=args)

    # per layer: pfennig2024 copy and the metabolism's fast pools
    met_reactions = met.get_raw_reactions()
    for i in range(layers):
        pre = f"z{i}_"
        inputs = {
            "absorption_ps1": f"{pre}complex_abs_ps1",
            "absorption_ps2": f"{pre}complex_abs_ps2",
            "absorption_pbs": f"{pre}complex_abs_pbs",
            "light_OCP": f"{pre}light_ocp",
        }
        _copy_pfennig(m, pf, pre, inputs)
        m.add_variables({pre + k: ref.pools[k] for k in PER_LAYER_METABOLITES})
        for mk, pk in PFENNIG_POOLS.items():
            m.add_derived(f"{pre}{mk}_mM", fn=to_mM, args=[pre + pk, "Chl"])
        bio = met.get_raw_derived()["mu"]
        table = {mk: f"{pre}{mk}_mM" for mk in PFENNIG_POOLS}
        table |= {k: pre + k for k in PER_LAYER_METABOLITES}
        m.add_derived(f"{pre}mu", fn=bio.fn, args=_rename(bio.args, table))

    fn, args = expr_fn(
        sp.Add(*(float(w[i]) * sp.Symbol(f"z{i}_mu") for i in range(layers))), "mu"
    )
    m.add_derived("mu", fn=fn, args=args)
    for name in ("E_CrtB_eff", "E_Chl_eff"):
        d = met.get_raw_derived()[name]
        m.add_derived(name, fn=d.fn, args=d.args)

    conv = sp.Integer(1000) / (
        60 * sp.Symbol("Chl")
    )  # mM min^-1 -> mmol molChl^-1 s^-1
    products = {"v_chl": "Chl", "v_car": "Car", "v_sq": "SQ"}
    for i in range(layers):
        pre = f"z{i}_"
        table = {mk: f"{pre}{mk}_mM" for mk in PFENNIG_POOLS}
        table |= {k: pre + k for k in PER_LAYER_METABOLITES}
        for name, r in met_reactions.items():
            if name.startswith(("supply_", "dilution_")):
                continue
            stoich: dict[str, Derived | float] = {}
            for k, c in r.stoichiometry.items():
                if k in PFENNIG_POOLS:
                    stoich[pre + PFENNIG_POOLS[k]] = _stoich(
                        c * conv, f"{pre}{name}_{k}"
                    )
                elif k in PER_LAYER_METABOLITES:
                    stoich[pre + k] = c / 60
                else:  # shared slow pool
                    stoich[k] = float(w[i]) * c / 60
            if name in products:  # synthesis feeds the shared pigment / squalene pools
                stoich[products[name]] = _stoich(
                    float(w[i]) * sp.Symbol("slow_speedup") / 60, f"{pre}{name}_out"
                )
            m.add_reaction(
                pre + name, fn=r.fn, args=_rename(r.args, table), stoichiometry=stoich
            )
        m.add_reaction(
            f"{pre}atp_ngam",
            fn=ngam_rate,
            args=[f"{pre}ATP_mM", "NGAM", "rho_DW", "KM_NGAM"],
            stoichiometry={pre + "ATP": _stoich(-conv, f"{pre}atp_ngam_ATP")},
        )
        for k in PER_LAYER_METABOLITES:
            m.add_reaction(
                f"{pre}dilution_{k}",
                fn=metabolism.dilution,
                args=[pre + k, "mu"],
                stoichiometry={pre + k: -1 / 60},
            )
    # per-cell synthesis rates [mM min^-1] (layer means) as derived values for the
    # evaluators: the rate laws evaluated once more per layer
    for v in ("v_chl", "v_car", "v_sq"):
        r = met_reactions[v]
        for i in range(layers):
            pre = f"z{i}_"
            table = {mk: f"{pre}{mk}_mM" for mk in PFENNIG_POOLS}
            table |= {k: pre + k for k in PER_LAYER_METABOLITES}
            m.add_derived(f"{pre}{v}_rate", fn=r.fn, args=_rename(r.args, table))
        total = sp.Add(
            *(float(w[i]) * sp.Symbol(f"z{i}_{v}_rate") for i in range(layers))
        )
        fn, args = expr_fn(total, v)
        m.add_derived(v, fn=fn, args=args)
    for k in SHARED_METABOLITES:
        m.add_reaction(
            f"dilution_{k}", fn=metabolism.dilution, args=[k, "mu"],
            stoichiometry={k: -1 / 60},
        )  # fmt: skip
    for s in ("Chl", "Car", "SQ"):
        m.add_reaction(
            f"dilution_{s}", fn=slow_dilution, args=[s, "mu", "slow_speedup"],
            stoichiometry={s: -1},
        )  # fmt: skip
    if dynamic:
        m.add_reaction("growth", fn=growth, args=["N", "mu"], stoichiometry={"N": 1})
    return m


def _pfennig_reference() -> dict[str, float]:
    """pfennig2024 steady state (clamped 3PGA, default light) for the reference levels."""
    from synmep.runners import NATIVE_RUNNER

    res = NATIVE_RUNNER.steady_state(get_petc_model())
    return res.value.variables.iloc[-1].to_dict()
