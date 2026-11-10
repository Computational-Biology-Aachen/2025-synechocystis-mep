"""v1: v0 with the ADR 0002 bugs fixed, one commit per fix (fix log: ADR 0002 §5).

Reproduces `model.get_model` + `main.ipynb:get_parameterised_model` 1:1, including
the bugs listed in ADR 0002. Names are sanitised to valid identifiers (`sanitise`) and every reaction is prefixed with `v_`, because the legacy model reuses
names across compounds, parameters and reactions (e.g. `DXS`, `3PGA_out_metabolism`).
"""

from __future__ import annotations

from typing import Literal

from mxlpy import Model

from synmep import kinetics_v1 as k
from synmep.data import CAROT, CELL_VOLUME, CHLA, SQUALENE, LegacyInputs
from synmep.data_v1 import load_inputs
from synmep.exprfn import clip_zero, expr_fn
from synmep.surrogates.petc_legacy import Symbols, petc_expressions

__all__ = ["get_model", "load_inputs"]

Strain = Literal["wt", "dshc"]
GrowthMode = Literal["steady", "dynamic"]

ENZYME_NAMES = [
    "DXS", "ThiG", "DXR", "IspD", "IspE", "IspF", "IspG", "IspH1", "IspH2", "Idi",
    "CrtE1", "CrtE2", "CrtE3", "Sqs", "Shc", "ChlP", "CrtB", "ChlG",
]  # fmt: skip
COMPOUNDS = [
    # MEP pathway
    "3PGA", "ATP", "NADPH", "Fdred", "DXP", "ThiDP", "MEP", "CDP_ME", "CDP_MEP",
    "MEcPP", "HMB_PP", "DMAPP", "IPP",
    # terpenoid synthesis and products
    "GPP", "FPP", "GGPP", "SQ", "carotenoids", "Hopanoids", "phytol",
    "chlorophyll-GGPP", "chlorophyll-phytol",
    # enzymes
    "DXS", "DXS_nonactive", "ThiG", "DXR", "IspD", "IspE", "IspF", "IspG", "IspH",
    "Idi", "CrtE", "Sqs", "Shc", "ChlP", "CrtB", "ChlG",
]  # fmt: skip


def sanitise(name: str) -> str:
    """Legacy name → valid identifier."""
    return name.replace("3PGA", "PGA3").replace("-", "_").replace(" ", "_")


def reaction_name(legacy: str) -> str:
    return f"v_{sanitise(legacy)}"


def legacy_parameters(
    inp: LegacyInputs,
    strain: Strain,
    growth_mode: GrowthMode,
) -> tuple[dict[str, float], dict[str, float]]:
    """Legacy parameters and initial conditions, keyed by legacy names."""
    del strain, growth_mode  # parameters are identical for all variants
    all_compounds = [*COMPOUNDS, "G3P", "PYR"]
    pars: dict[str, float] = {
        "light": 100.0,
        "ATP_out_metabolism": 0.5,
        "3PGA_out_metabolism": 0.5,
        "NADPH_out_metabolism": 0.5,
        "Fdred_out_metabolism": 0.5,
        "Kna_DXS": 0.000004,
        "DXS_na_kdeg": 1 / 50000,
        "DXS_inh_str": 0.1,
        "ps_ratio": 5.9,
        "beta_carotene_fraction": 0.26,
        "sample_depth_m": 0.01,
        "cell_volume": CELL_VOLUME,
        "NADPHsink_kfwd": 0.1,
        "NADPHsink_krev": 0.1,
        "Fdredsink_kfwd": 0.1,
        "Fdredsink_krev": 0.1,
        "phycocyanin": 53180 * 2e-15 / CELL_VOLUME,
        "allophycocyanin": 24076 * 2e-15 / CELL_VOLUME,
    }
    ep = inp.enzyme_parameters
    for enzyme in ENZYME_NAMES:
        pars[f"{enzyme}Kcat"] = ep[f"{enzyme}Kcat"]
        pars[f"{enzyme}_Keq"] = ep[f"{enzyme}_Keq"]
        pars[f"{enzyme}_h"] = 1
        for s in all_compounds:
            for p in ("Km", "Ki", "Kd"):
                if (key := f"{enzyme}{s}{p}") in ep:
                    pars[key] = ep[key]
    for d, value in inp.enzyme_concentrations_umol.items():
        pars[f"{d}_ksynth"] = value / (5.5 * 60) * 2
        for suffix in ("KmA", "KmB", "KmC", "KmD"):
            pars[f"{d}_synth_{suffix}"] = 1
        pars[f"{d}_synth_h"] = 1
    # notebook overrides
    pars |= {"light": 100, "cell_volume": CELL_VOLUME}
    pars |= {"DXS_inh_str": 0.001, "ThiGDXPKm": 0.01}
    pars |= {
        key: v
        for key, v in inp.min_kcat_values.items()
        if key not in {"IspHKcat", "CrtEKcat"}
    }

    y0 = dict.fromkeys(COMPOUNDS, 0.01)
    y0 |= inp.enzyme_concentrations_umol
    y0 |= {
        "3PGA": 0.1,
        "ATP": 0.1,
        "NADPH": 0.1,
        "Fdred": 0.1,
        "ThiDP": 10,
        "chlorophyll-phytol": CHLA / CELL_VOLUME,
        "Cell_density": 1e1,
        "carotenoids": CAROT / CELL_VOLUME,
        "SQ": SQUALENE / CELL_VOLUME,
    }
    return {kk: float(v) for kk, v in pars.items()}, {
        kk: float(v) for kk, v in y0.items()
    }


def get_model(
    inp: LegacyInputs,
    strain: Strain = "dshc",
    growth_mode: GrowthMode = "steady",
) -> Model:
    pars, y0 = legacy_parameters(inp, strain, growth_mode)
    S = sanitise

    m = Model()
    m.add_variables({S(c): y0[c] for c in [*COMPOUNDS, "Cell_density"]})
    m.add_parameters({S(p): v for p, v in pars.items()})

    def rxn(name: str, fn, args: list[str], stoich: dict[str, float]) -> None:
        m.add_reaction(
            reaction_name(name),
            fn=fn,
            args=[S(a) for a in args],
            stoichiometry={S(c): float(v) for c, v in stoich.items()},
        )

    m.add_derived(
        "Chla_synthesis",
        fn=k.chl_synthesis_rate,
        args=[
            S(a)
            for a in [
                "chlorophyll-GGPP",
                "NADPH",
                "chlorophyll-phytol",
                "GGPP",
                "ChlPKcat",
                "ChlP",
                "ChlPchlorophyll-GGPPKm",
                "ChlPGGPPKm",
                "ChlPNADPHKm",
                "ChlPchlorophyll-phytolKm",
                "ChlPchlorophyll-GGPPKi",
                "ChlP_Keq",
                "phytol",
                "ChlGKcat",
                "ChlG",
                "ChlGphytolKm",
                "ChlGchlorophyll-phytolKm",
                "ChlG_Keq",
                "ChlGGGPPKm",
            ]
        ],
    )
    m.add_derived(
        "Growth_rate", fn=k.growth_function, args=["Chla_synthesis", "cell_volume"]
    )

    if growth_mode == "dynamic":
        rxn(
            "Cell_density_increase",
            k.mass_action_1,
            ["Cell_density", "Growth_rate"],
            {"Cell_density": 1},
        )
    for c in sorted(COMPOUNDS):
        rxn(f"{c}_out", k.mass_action_1, [c, "Growth_rate"], {c: -1})

    rxn(
        "Biomass",
        k.growth_function,
        ["Chla_synthesis", "cell_volume"],
        {
            "NADPH": -8.762903549999999 / 60,
            "ATP": -33.0939284280016 / 60,
            "Fdred": -25.97602870400234 / 60,
            "3PGA": -78.47848943300428 / 60,
        },
    )

    # PETC surrogate incl. self-shading, as explicit expressions (no model data)
    petc = petc_expressions(
        inp,
        Symbols(
            chl=S("chlorophyll-phytol"),
            carotenoids="carotenoids",
            cell_density="Cell_density",
            light="light",
            phycocyanin="phycocyanin",
            allophycocyanin="allophycocyanin",
            ps_ratio="ps_ratio",
            beta_carotene_fraction="beta_carotene_fraction",
            sample_depth="sample_depth_m",
            cell_volume="cell_volume",
        ),
    )
    for name, expr in petc.items():
        fn, args = expr_fn(expr, origin=S(name))
        m.add_derived(S(name), fn=fn, args=args)
    for c in ("ATP", "NADPH", "3PGA", "Fdred"):
        m.add_derived(S(f"v{c}_prod"), fn=clip_zero, args=[S(f"v{c}_prod_raw")])
        rxn(f"{c} production", k.unit_conversion, [f"v{c}_prod", "cell_volume"], {c: 1})

    m.add_derived("G3P", fn=k.fraction_PGA_derivates_G3P, args=[S("3PGA")])
    m.add_derived("PYR", fn=k.fraction_PGA_derivates_PYR, args=[S("3PGA")])

    for c in ("3PGA", "ATP", "NADPH", "Fdred"):
        rxn(f"{c}_out_metabolism", k.mass_action_1, [c, f"{c}_out_metabolism"], {c: -1})

    costs = inp.enzyme_costs
    sf = 1e-4
    cofactors = ("NADPH", "ATP", "Fdred", "3PGA")

    def synth_stoich(enzyme: str) -> dict[str, float]:
        st = {c: sf * -float(costs.loc[enzyme, c]) for c in cofactors}  # type: ignore[arg-type]
        return st | {enzyme: sf * 1000}

    synth_args = ["ATP", "NADPH", "Fdred", "3PGA"]
    for enzyme in costs.index.delete(0):
        rxn(
            f"{enzyme}_synth",
            k.ping_pong_4,
            [*synth_args, f"{enzyme}_ksynth"]
            + [f"{enzyme}_synth_Km{x}" for x in "ABCD"],
            synth_stoich(enzyme),
        )
    rxn(
        "DXS_synth",
        k.DXS_synth_f,
        [*synth_args, "DXS_ksynth"]
        + [f"DXS_synth_Km{x}" for x in "ABCD"]
        + ["DXS_nonactive", "DXS_inh_str"],
        synth_stoich("DXS"),
    )
    rxn(
        "DXS_monomerization",
        k.DXS_na,
        [
            "DXS",
            "Kna_DXS",
            "IPP",
            "DMAPP",
            "DXSIPPKi",
            "DXSDMAPPKi",
            "ThiDP",
            "DXSThiDPKd",
        ],
        {"DXS": -1, "DXS_nonactive": 1},
    )
    rxn(
        "DXS_degradation",
        k.mass_action_1,
        ["DXS_nonactive", "DXS_na_kdeg"],
        {"DXS_nonactive": -1},
    )
    rxn(
        "DXS",
        k.v_dxs,
        ["PYR", "G3P", "DXS", "DXSKcat", "DXSPYRKm", "DXSG3PKm", "DXSPYRKi",
         "IPP", "DMAPP", "ThiDP", "DXSIPPKi", "DXSDMAPPKi", "DXSThiDPKd"],
        {"3PGA": -2, "DXP": 1},
    )  # fmt: skip
    rxn(
        "ThiG",
        k.reversible_michaelis_menten_keq_kcat,
        ["DXP", "ThiDP", "ThiGKcat", "ThiG", "ThiGDXPKm", "ThiGThiDPKm", "ThiG_Keq"],
        {"DXP": -1, "ATP": -1, "ThiDP": 1},
    )

    def hill(name: str, a: str, b: str, p: str, enzyme: str, stoich) -> None:
        e = name
        rxn(
            name,
            k.universal_hill_2_1,
            [a, b, p, f"{e}Kcat", enzyme, f"{e}{a}Km", f"{e}{b}Km", f"{e}{p}Km",
             f"{e}_Keq", f"{e}_h"],
            stoich,
        )  # fmt: skip

    hill("DXR", "DXP", "NADPH", "MEP", "DXR", {"DXP": -1, "NADPH": -1, "MEP": 1})
    hill("IspD", "MEP", "ATP", "CDP_ME", "IspD", {"MEP": -1, "ATP": -2, "CDP_ME": 1})
    hill(
        "IspE",
        "CDP_ME",
        "ATP",
        "CDP_MEP",
        "IspE",
        {"CDP_ME": -1, "ATP": -1, "CDP_MEP": 1},
    )
    rxn(
        "IspF",
        k.reversible_michaelis_menten_keq_kcat,
        ["CDP_MEP", "MEcPP", "IspFKcat", "IspF", "IspFCDP_MEPKm", "IspFMEcPPKm",
         "IspF_Keq"],
        {"CDP_MEP": -1, "MEcPP": 1},
    )  # fmt: skip
    hill(
        "IspG",
        "MEcPP",
        "Fdred",
        "HMB_PP",
        "IspG",
        {"MEcPP": -1, "Fdred": -1, "HMB_PP": 1},
    )

    m.add_derived("IspH1_enzyme", fn=k.IspH1, args=["IspH"])
    m.add_derived("IspH2_enzyme", fn=k.IspH2, args=["IspH"])
    hill(
        "IspH1",
        "HMB_PP",
        "NADPH",
        "IPP",
        "IspH1_enzyme",
        {"HMB_PP": -1, "NADPH": -1, "IPP": 1},
    )
    hill(
        "IspH2",
        "HMB_PP",
        "NADPH",
        "DMAPP",
        "IspH2_enzyme",
        {"HMB_PP": -1, "NADPH": -1, "DMAPP": 1},
    )
    rxn(
        "Idi",
        k.reversible_michaelis_menten_keq_kcat,
        ["IPP", "DMAPP", "IdiKcat", "Idi", "IdiIPPKm", "IdiDMAPPKm", "Idi_Keq"],
        {"IPP": -1, "DMAPP": 1},
    )
    rxn(
        "CrtE1",
        k.reversible_ordered_2_comp2,
        ["DMAPP", "IPP", "GPP", "GPP", "FPP", "CrtE1Kcat", "CrtE", "CrtE1DMAPPKm",
         "CrtE1GPPKm", "CrtE1FPPKm", "CrtE1IPPKm", "CrtE1GPPKm", "CrtE1DMAPPKi",
         "CrtE1_Keq"],
        {"DMAPP": -1, "IPP": -1, "GPP": 1},
    )  # fmt: skip
    rxn(
        "CrtE2",
        k.reversible_ordered_2_comp2,
        ["GPP", "IPP", "FPP", "DMAPP", "FPP", "CrtE2Kcat", "CrtE", "CrtE2GPPKm",
         "CrtE2DMAPPKm", "CrtE2FPPKm", "CrtE2IPPKm", "CrtE2FPPKm", "CrtE2GPPKi",
         "CrtE2_Keq"],
        {"GPP": -1, "IPP": -1, "FPP": 1},
    )  # fmt: skip
    rxn(
        "CrtE3",
        k.reversible_ordered_2_comp2,
        ["FPP", "IPP", "GGPP", "DMAPP", "GPP", "CrtE3Kcat", "CrtE", "CrtE3FPPKm",
         "CrtE3DMAPPKm", "CrtE3GPPKm", "CrtE3IPPKm", "CrtE3GGPPKm", "CrtE3FPPKi",
         "CrtE3_Keq"],
        {"FPP": -1, "IPP": -1, "GGPP": 1},
    )  # fmt: skip
    rxn(
        "Sqs",
        k.v_sqs,
        ["FPP", "FPP", "NADPH", "SQ", "SqsKcat", "Sqs", "SqsFPPKm", "SqsFPPKm",
         "SqsNADPHKm", "SqsSQKm", "Sqs_Keq", "Sqs_h"],
        {"FPP": -2, "NADPH": -1, "SQ": 1},
    )  # fmt: skip
    if strain == "wt":
        rxn(
            "Shc",
            k.reversible_michaelis_menten_keq_kcat,
            ["SQ", "Hopanoids", "ShcKcat", "Shc", "ShcSQKm", "ShcHopanoidsKm",
             "Shc_Keq"],
            {"SQ": -1, "Hopanoids": 1},
        )  # fmt: skip
    rxn(
        "carotenoids",
        k.reversible_michaelis_menten_keq_kcat,
        ["GGPP", "carotenoids", "CrtBKcat", "CrtB", "CrtBGGPPKm",
         "CrtBcarotenoidsKm", "CrtB_Keq"],
        {"GGPP": -2, "carotenoids": 1},
    )  # fmt: skip
    rxn(
        "phytol_synthase1",
        k.reversible_ordered_2_comp1,
        ["GGPP", "NADPH", "phytol", "chlorophyll-GGPP", "ChlPKcat", "ChlP",
         "ChlPGGPPKm", "ChlPchlorophyll-GGPPKm", "ChlPNADPHKm", "ChlPphytolKm",
         "ChlPGGPPKi", "ChlP_Keq"],
        {"GGPP": -1, "NADPH": -3, "phytol": 1},
    )  # fmt: skip
    rxn(
        "phytol_synthase2",
        k.reversible_ordered_2_comp1,
        ["chlorophyll-GGPP", "NADPH", "chlorophyll-phytol", "GGPP", "ChlPKcat",
         "ChlP", "ChlPchlorophyll-GGPPKm", "ChlPGGPPKm", "ChlPNADPHKm",
         "ChlPchlorophyll-phytolKm", "ChlPchlorophyll-GGPPKi", "ChlP_Keq"],
        {"chlorophyll-GGPP": -1, "NADPH": -3, "chlorophyll-phytol": 1},
    )  # fmt: skip
    rxn(
        "chlorophyll_synthesis1",
        k.reversible_michaelis_menten_keq_kcat_comp1,
        ["phytol", "chlorophyll-phytol", "ChlGKcat", "ChlG", "ChlGphytolKm",
         "ChlGchlorophyll-phytolKm", "ChlG_Keq", "GGPP", "ChlGGGPPKm"],
        {"phytol": -1, "chlorophyll-phytol": 1},
    )  # fmt: skip
    rxn(
        "chlorophyll_synthesis2",
        k.reversible_michaelis_menten_keq_kcat_comp1,
        ["GGPP", "chlorophyll-GGPP", "ChlGKcat", "ChlG", "ChlGGGPPKm",
         "ChlGchlorophyll-GGPPKm", "ChlG_Keq", "phytol", "ChlGphytolKm"],
        {"GGPP": -1, "chlorophyll-GGPP": 1},
    )  # fmt: skip
    rxn(
        "NADPH_sink_Fdred",
        k.reversible_mass_action_1_1,
        ["ATP", "Fdred", "Fdredsink_kfwd", "Fdredsink_krev"],
        {"ATP": -1, "Fdred": 1},
    )
    rxn(
        "NADPH_sink_ATP",
        k.reversible_mass_action_1_1,
        ["ATP", "NADPH", "NADPHsink_kfwd", "NADPHsink_krev"],
        {"ATP": -1, "NADPH": 1},
    )
    m.remove_parameters(m.get_unused_parameters())
    return m
