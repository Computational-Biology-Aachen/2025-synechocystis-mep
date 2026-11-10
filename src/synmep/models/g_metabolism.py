"""Ground-up metabolism layer (minutes timescale), thread "g" (ADR 0004).

Fast pools [mM]: PGA (3PGA and all sugar phosphates), ATP, NADPH, MEPi (lumped MEP
pathway intermediates DXP…HMBPP), IPP (IPP + DMAPP), FPP, GGPP.

Inputs (parameters, set by the slower layers):
- `S_PGA`, `S_ATP` [mM min^-1]: supply from the PETC layer (`J_3PGA`, `J_ATP` × Chl).
  pfennig2024 has no explicit NADPH demand, so NADPH supply is `r_NADPH · S_PGA`.
- `E_<enzyme>` [-]: enzyme level relative to the wild type reference (1 = WT);
  `E_CrtE` acts on both the FPP and the GGPP synthase step

Outputs: `mu` [min^-1] (growth rate from the biomass synthesis flux) and the product
fluxes `v_chl`, `v_car`, `v_sq` [mM min^-1] into the slow pools.

Rate laws are irreversible saturable (Michaelis-Menten products) with product
inhibition at the branch points (keeps pools bounded when growth, and with it
dilution, stops). Every pool is diluted by growth. Maximal rates and sink constants
are set by construction (`reference_parameters`): at a reference state each enzyme runs
at about half saturation while carrying the flux needed to keep the measured
pigment/product contents at the reference growth rate (the legacy "minimal kcat x 2"
idea, made explicit), and the pools sit at their reference values.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from mxlpy import Model

# biomass composition flux [mmol gDW^-1] (paper Tab. biomass_stoichiometry, Knoop2013)
BIOMASS = {"PGA": 6.934, "ATP": 0.774, "NADPH": 2.924}
POOLS = ("PGA", "ATP", "NADPH", "MEPi", "IPP", "FPP", "GGPP")
ENZYMES = ("DXS", "MEP", "CrtE", "Sqs", "CrtB", "Chl")  # CrtE makes FPP and GGPP
SUPPLIES = ("S_PGA", "S_ATP")
OUTPUTS = ("mu", "v_chl", "v_car", "v_sq")
# what the growth layer needs: mu and the steady-state contents X* = v_X / mu [mM]
# (bounded, unlike fluxes divided by a vanishing mu; ADR 0004)
SURROGATE_OUTPUTS = ("mu", "Chl_star", "Car_star", "SQ_star")

# Km / Ki [mM]: legacy blast/literature values where available (docs/assumptions.md)
KM = {
    "DXS_PGA": 0.24,  # DXS Km(G3P)
    "MEP_MEPi": 0.195,  # DXR Km(DXP), first step of the lumped chain
    "MEP_NADPH": 0.035,  # (assumption: as ChlP) cofactor saturation of the chain
    "MEP_ATP": 0.1,  # (assumption)
    "FPPS_IPP": 0.0034,  # CrtE Km(IPP)
    "GGPPS_FPP": 0.275,  # CrtE Km(FPP)
    "GGPPS_IPP": 0.0034,
    "Sqs_FPP": 0.0028,
    "Sqs_NADPH": 0.5,
    "CrtB_GGPP": 0.005,
    "Chl_GGPP": 0.06,  # ChlP Km(GGPP)
    "Chl_NADPH": 0.035,  # ChlP Km(NADPH)
}
KI = {
    "DXS_IPP": 0.065,  # IPP/DMAPP on DXS (poplar DXS, Banerjee2013)
    "DXS_MEPi": 1.0,  # (assumption) MEP intermediates on DXS
    "FPPS_FPP": 0.5,  # (assumption) product inhibition
    "GGPPS_GGPP": 0.5,  # (assumption) product inhibition
}
RHO_DW = 970.0  # [gDW L(cell)^-1] 3.88 pg DW cell^-1 (proteomics samples) / 4 fL
MU_MAX = math.log(2) / (5.5 * 60)  # [min^-1] turbidostat doubling time 5.5 h (paper)


def mm(s: float, km: float) -> float:
    return s / (km + s)


def v_biomass(pga: float, atp: float, nadph: float, vmax: float, k_pga: float,
              k_atp: float, k_nadph: float) -> float:  # fmt: skip
    return vmax * mm(pga, k_pga) * mm(atp, k_atp) * mm(nadph, k_nadph)


def growth_rate(pga: float, atp: float, nadph: float, vmax: float, k_pga: float,
                k_atp: float, k_nadph: float, rho_dw: float) -> float:  # fmt: skip
    return v_biomass(pga, atp, nadph, vmax, k_pga, k_atp, k_nadph) / rho_dw


def v_dxs(pga: float, ipp: float, mepi: float, e: float, vmax: float, km: float,
          ki_ipp: float, ki_mepi: float) -> float:  # fmt: skip
    return e * vmax * mm(pga, km) / (1 + ipp / ki_ipp + mepi / ki_mepi)


def v_mep(mepi: float, nadph: float, atp: float, e: float, vmax: float, km: float,
          km_nadph: float, km_atp: float) -> float:  # fmt: skip
    return e * vmax * mm(mepi, km) * mm(nadph, km_nadph) * mm(atp, km_atp)


def v_mm1(s: float, e: float, vmax: float, km: float) -> float:
    return e * vmax * mm(s, km)


def v_mm1_pi(s: float, p: float, e: float, vmax: float, km: float, ki: float) -> float:
    return e * vmax * mm(s, km) / (1 + p / ki)


def v_mm2(s1: float, s2: float, e: float, vmax: float, km1: float, km2: float) -> float:
    return e * vmax * mm(s1, km1) * mm(s2, km2)


def v_mm2_pi(s1: float, s2: float, p: float, e: float, vmax: float, km1: float,
             km2: float, ki: float) -> float:  # fmt: skip
    return e * vmax * mm(s1, km1) * mm(s2, km2) / (1 + p / ki)


def proportional(x: float, r: float) -> float:
    return r * x


def first_order(s: float, k: float) -> float:
    return k * s


def dilution(s: float, mu: float) -> float:
    return mu * s


def constant(x: float) -> float:
    return x


@dataclass(frozen=True)
class Reference:
    """Reference state used to set maximal rates and sinks (4 fL cells, ADR 0004)."""

    mu: float = (
        0.04 / 60
    )  # [min^-1] initial batch growth rate (growth_rates_by_density)
    chl: float = 19.9  # [mM] 7.5 mg Chl OD^-1 L^-1, 1.0548e8 cells mL^-1 OD^-1, 4 fL
    carotenoids: float = 8.15  # [mM] 2 mg OD^-1 L^-1 (M 581.6), same conversion
    squalene: float = 2.9  # [mM] 0.5 mg OD^-1 L^-1 (M 410.7), sqs strain
    j_atp_petc: float = 3.14  # [mmol mol(Chl)^-1 s^-1] PETC J_ATP at the Fig 9 median
    pools: dict[str, float] = field(
        default_factory=lambda: {
            "PGA": 1.0,  # [mM] lumped sugar phosphates (assumption)
            "ATP": 2.0,  # [mM] typical cyanobacterial ATP (Doello2018 ~1-3 mM)
            "NADPH": 0.5,  # [mM] (assumption, Cooley2001 order of magnitude)
            "MEPi": KM["MEP_MEPi"],
            "IPP": KI["DXS_IPP"],
            "FPP": KM["Sqs_FPP"],
            "GGPP": KM["CrtB_GGPP"],
        }
    )
    pga_overflow: float = 0.05  # fraction of PGA supply going to storage at reference
    nadph_overflow: float = 0.2  # fraction of NADPH supply oxidised at reference


REFERENCE = Reference()


def reference_fluxes(ref: Reference = REFERENCE) -> dict[str, float]:
    """Fluxes at the reference state, including dilution of the intermediate pools."""
    mu, p = ref.mu, ref.pools
    j = {"chl": mu * ref.chl, "car": mu * ref.carotenoids, "sq": mu * ref.squalene}
    j["ggpps"] = j["chl"] + 2 * j["car"] + mu * p["GGPP"]
    j["fpps"] = j["ggpps"] + 2 * j["sq"] + mu * p["FPP"]
    # C5 units: 3 per FPP, 1 more per GGPP, plus dilution of IPP
    j["mep"] = 3 * j["fpps"] + j["ggpps"] + mu * p["IPP"]
    j["dxs"] = j["mep"] + mu * p["MEPi"]
    j["ipp"] = j["mep"]
    j["bio"] = mu * RHO_DW  # [gDW L^-1 min^-1]
    return j


def reference_parameters(ref: Reference = REFERENCE) -> dict[str, float]:
    """Maximal rates and sink constants reproducing the reference state."""
    j, p, sat = reference_fluxes(ref), ref.pools, 0.5
    dil = {k: ref.mu * v for k, v in p.items()}

    pga_demand = BIOMASS["PGA"] * j["bio"] + 2 * j["dxs"] + dil["PGA"]
    s_pga = pga_demand / (1 - ref.pga_overflow)
    nadph_demand = (
        BIOMASS["NADPH"] * j["bio"] + 3 * j["mep"] + 3 * j["chl"] + j["sq"]
    ) + dil["NADPH"]
    s_nadph = nadph_demand / (1 - ref.nadph_overflow)
    s_atp = ref.j_atp_petc * ref.chl * 0.06
    atp_demand = BIOMASS["ATP"] * j["bio"] + 3 * j["mep"] + dil["ATP"]
    bio_sat = (ref.mu / MU_MAX) ** (1 / 3)
    return {
        "V_DXS": j["dxs"]
        / (mm(p["PGA"], KM["DXS_PGA"]) / (1 + 1 + p["MEPi"] / KI["DXS_MEPi"])),
        "V_MEP": j["mep"]
        / (sat * mm(p["NADPH"], KM["MEP_NADPH"]) * mm(p["ATP"], KM["MEP_ATP"])),
        "V_FPPS": j["fpps"]
        * (1 + p["FPP"] / KI["FPPS_FPP"])
        / mm(p["IPP"], KM["FPPS_IPP"]),
        "V_GGPPS": j["ggpps"]
        * (1 + p["GGPP"] / KI["GGPPS_GGPP"])
        / (mm(p["FPP"], KM["GGPPS_FPP"]) * mm(p["IPP"], KM["GGPPS_IPP"])),
        "V_Sqs": j["sq"] / (sat * mm(p["NADPH"], KM["Sqs_NADPH"])),
        "V_CrtB": j["car"] / sat,
        "V_Chl": j["chl"]
        / (mm(p["GGPP"], KM["Chl_GGPP"]) * mm(p["NADPH"], KM["Chl_NADPH"])),
        # capacity mu_max; each cofactor term equally saturated so that mu = mu_ref
        "V_bio": MU_MAX * RHO_DW,
        **{f"Km_bio_{s}": p[s] * (1 / bio_sat - 1) for s in BIOMASS},
        "r_NADPH": s_nadph / s_pga,
        "k_PGA_storage": ref.pga_overflow * s_pga / p["PGA"],
        "k_NADPH_oxidation": ref.nadph_overflow * s_nadph / p["NADPH"],
        "k_ATP_maintenance": max(s_atp - atp_demand, 0.0) / p["ATP"],
        "S_PGA": s_pga,
        "S_ATP": s_atp,
    }


def get_model(ref: Reference = REFERENCE) -> Model:
    m = Model()
    m.add_variables(dict(ref.pools))
    m.add_parameters(
        {f"E_{e}": 1.0 for e in ENZYMES}
        | {f"Km_{k}": v for k, v in KM.items()}
        | {f"Ki_{k}": v for k, v in KI.items()}
        | reference_parameters(ref)
        | {"rho_DW": RHO_DW}
    )
    m.add_derived("S_NADPH", fn=proportional, args=["S_PGA", "r_NADPH"])
    for s in ("PGA", "ATP", "NADPH"):
        m.add_reaction(
            f"supply_{s}", fn=constant, args=[f"S_{s}"], stoichiometry={s: 1}
        )

    bio_args = [
        "PGA", "ATP", "NADPH", "V_bio", "Km_bio_PGA", "Km_bio_ATP", "Km_bio_NADPH",
    ]  # fmt: skip
    m.add_reaction(
        "v_bio",
        fn=v_biomass,
        args=bio_args,
        stoichiometry={s: -BIOMASS[s] for s in BIOMASS},
    )
    m.add_derived("mu", fn=growth_rate, args=[*bio_args, "rho_DW"])
    m.add_reaction(
        "v_dxs",
        fn=v_dxs,
        args=["PGA", "IPP", "MEPi", "E_DXS", "V_DXS", "Km_DXS_PGA", "Ki_DXS_IPP",
              "Ki_DXS_MEPi"],
        stoichiometry={"PGA": -2, "MEPi": 1},  # G3P + PYR -> DXP
    )  # fmt: skip
    m.add_reaction(
        "v_mep",
        fn=v_mep,
        args=["MEPi", "NADPH", "ATP", "E_MEP", "V_MEP", "Km_MEP_MEPi",
              "Km_MEP_NADPH", "Km_MEP_ATP"],
        # DXR (NADPH), IspD (CTP~ATP), IspE (ATP), IspG/IspH (2 Fd each ~ 2 NADPH)
        stoichiometry={"MEPi": -1, "NADPH": -3, "ATP": -3, "IPP": 1},
    )  # fmt: skip
    m.add_reaction(
        "v_fpps",
        fn=v_mm1_pi,
        args=["IPP", "FPP", "E_CrtE", "V_FPPS", "Km_FPPS_IPP", "Ki_FPPS_FPP"],
        stoichiometry={"IPP": -3, "FPP": 1},
    )
    m.add_reaction(
        "v_ggpps",
        fn=v_mm2_pi,
        args=["FPP", "IPP", "GGPP", "E_CrtE", "V_GGPPS", "Km_GGPPS_FPP",
              "Km_GGPPS_IPP", "Ki_GGPPS_GGPP"],
        stoichiometry={"FPP": -1, "IPP": -1, "GGPP": 1},
    )  # fmt: skip
    m.add_reaction(
        "v_sq",
        fn=v_mm2,
        args=["FPP", "NADPH", "E_Sqs", "V_Sqs", "Km_Sqs_FPP", "Km_Sqs_NADPH"],
        stoichiometry={"FPP": -2, "NADPH": -1},
    )
    m.add_reaction(
        "v_car",
        fn=v_mm1,
        args=["GGPP", "E_CrtB", "V_CrtB", "Km_CrtB_GGPP"],
        stoichiometry={"GGPP": -2},
    )
    m.add_reaction(
        "v_chl",
        fn=v_mm2,
        args=["GGPP", "NADPH", "E_Chl", "V_Chl", "Km_Chl_GGPP", "Km_Chl_NADPH"],
        stoichiometry={"GGPP": -1, "NADPH": -3},  # ChlG + ChlP (3 NADPH)
    )
    for name, pool, k in (
        ("atp_maintenance", "ATP", "k_ATP_maintenance"),
        ("nadph_oxidation", "NADPH", "k_NADPH_oxidation"),
        ("pga_storage", "PGA", "k_PGA_storage"),
    ):
        m.add_reaction(name, fn=first_order, args=[pool, k], stoichiometry={pool: -1})
    for p in POOLS:
        m.add_reaction(
            f"dilution_{p}", fn=dilution, args=[p, "mu"], stoichiometry={p: -1}
        )
    return m
