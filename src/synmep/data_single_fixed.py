"""Inputs of `single-fixed`: the legacy preprocessing with the ADR 0002 data bugs fixed.

Only the parts that differ from `synmep.data` / `synmep.legacy_data` live here.
"""

# ruff: noqa: PLC0206  (function body derived from verbatim legacy code)

from __future__ import annotations

import math
from dataclasses import replace
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from synmep.data import DATA, LegacyInputs, load_legacy_inputs
from synmep.legacy_data import calculate_min_kcats

DOUBLING_TIME = 330.0  # [min] 5.5 h turbidostat doubling time (paper Results)
CELL_VOLUME = 4e-15  # [L cell^-1] (paper Methods; legacy code: 2 fL, B9)


def get_enzyme_parameters(  # B6: defaults in mM
    *,
    mode: Literal["blast", "AI"],
    blasted_parameters_csv: Path,  # enzyme_parameters_blasted.csv
    predicted_kms_csv: Path,  # predicted_Kms.csv
    equilibrium_constants_csv: Path,  # Equilibrium_constants.csv
) -> dict:
    Enzyme_names: list[str] = [
        "DXS",
        "ThiG",
        "DXR",
        "IspD",
        "IspE",
        "IspF",
        "IspG",
        "IspH1",
        "IspH2",
        "Idi",
        "CrtE1",
        "CrtE2",
        "CrtE3",
        "Sqs",
        "Shc",
        "ChlP",
        "CrtB",
        "ChlG",
    ]
    enzyme_parameters = pd.read_csv(
        blasted_parameters_csv,
        delimiter=";",
    )
    AI_kms = pd.read_csv(
        predicted_kms_csv,
        delimiter=";",
    )
    AI_kms = AI_kms[
        AI_kms["organism"].isin(
            ["Synechocystis sp. PCC6803", "Synechocystis sp. PCC 6803"]
        )
    ]

    K_eqs = pd.read_csv(
        equilibrium_constants_csv,
        delimiter=";",
        index_col="Enzyme",
    )
    K_eqs["Keq"] = K_eqs["Keq"].astype(float)
    enzyme_parameters.dropna(how="all", inplace=True)
    enzyme_parameters_dict_list = {}
    enzyme_parameters["Km [µM]"] = enzyme_parameters["Km [uM]"].astype(float)
    enzyme_parameters["Kcat [s^-1]"] = enzyme_parameters["Kcat [s^-1]"].astype(float)
    enzyme_parameters["Ki [µM]"] = enzyme_parameters["Ki [uM]"].astype(float)
    enzyme_parameters["Kd [µM]"] = enzyme_parameters["Kd [uM]"].astype(float)
    substrates = [
        "ATP",
        "NADPH",
        "Fdred",
        "G3P",
        "PYR",
        "DXP",
        "ThiDP",
        "MEP",
        "CDP_ME",
        "CDP_MEP",
        "MEcPP",
        "HMB_PP",
        "DMAPP",
        "IPP",
        "GPP",
        "FPP",
        "GGPP",
        "SQ",
        "carotenoids",
        "Hopanoids",
        "phytol",
        "chlorophyll-GGPP",
        "chlorophyll-phytol",
    ]

    # B6 (ADR 0002): the legacy defaults were applied after the µM -> mM conversion of
    # the data, so they ended up as 100 mM / 1000 mM / 1 mM. They are meant in µM.
    default_values = {
        "eq": 10**10,  # (legacy, unsourced) effectively irreversible
        "Km": 100 / 1000,  # 100 µM (legacy, unsourced)
        "at": 600,  # [min^-1] (legacy, unsourced); overwritten by min kcats
        "Ki": 1000 / 1000,  # 1000 µM (legacy, unsourced)
        "Kd": 1 / 1000,  # 1 µM (legacy, unsourced)
    }
    for name in Enzyme_names:
        for subs in substrates:
            tempdict = {}
            if mode == "AI":
                tempdict["Km" + str(subs)] = (
                    AI_kms.loc[
                        (AI_kms["Enzyme"] == name) & (AI_kms["Substrate"] == subs)
                    ]["predicted_km"]
                    .dropna(how="all")
                    .mean()
                    / 1000
                )
            elif mode == "blast":
                if np.isnan(
                    enzyme_parameters.loc[
                        (enzyme_parameters["Enzyme"] == name)
                        & (enzyme_parameters["Substrate"] == subs)
                    ]["Km [µM]"]
                    .dropna(how="all")
                    .mean()
                ):
                    tempdict["Km" + str(subs)] = (
                        AI_kms.loc[
                            (AI_kms["Enzyme"] == name) & (AI_kms["Substrate"] == subs)
                        ]["predicted_km"]
                        .dropna(how="all")
                        .mean()
                        / 1000
                    )
                else:
                    tempdict["Km" + str(subs)] = (
                        enzyme_parameters.loc[
                            (enzyme_parameters["Enzyme"] == name)
                            & (enzyme_parameters["Substrate"] == subs)
                        ]["Km [µM]"]
                        .dropna(how="all")
                        .mean()
                        / 1000
                    )

            tempdict["Ki" + str(subs)] = (
                enzyme_parameters.loc[
                    (enzyme_parameters["Enzyme"] == name)
                    & (enzyme_parameters["Substrate"] == subs)
                ]["Ki [µM]"]
                .dropna(how="all")
                .mean()
                / 1000
            )
            tempdict["Kd" + str(subs)] = (
                enzyme_parameters.loc[
                    (enzyme_parameters["Enzyme"] == name)
                    & (enzyme_parameters["Substrate"] == subs)
                ]["Kd [µM]"]
                .dropna(how="all")
                .mean()
                / 1000
            )

            enzyme_parameters_dict_list[name + subs + "Km"] = tempdict["Km" + str(subs)]
            enzyme_parameters_dict_list[name + subs + "Ki"] = tempdict["Ki" + str(subs)]
            enzyme_parameters_dict_list[name + subs + "Kd"] = tempdict["Kd" + str(subs)]
        enzyme_parameters_dict_list[name + "_Keq"] = K_eqs.at[name, "Keq"]
        enzyme_parameters_dict_list[name + "Kcat"] = (
            enzyme_parameters.loc[(enzyme_parameters["Enzyme"] == name)]["Kcat [s^-1]"]
            .dropna(how="all")
            .mean()
            * 60
        )

    # Handle enzyme-substrate combinations not included in the CSV
    for key in enzyme_parameters_dict_list:
        if np.isnan(enzyme_parameters_dict_list[key]):
            enzyme_parameters_dict_list[key] = default_values[key[-2:]]

    return pd.Series(enzyme_parameters_dict_list, dtype=float).to_dict()


def load_inputs(data: Path = DATA) -> LegacyInputs:
    inp = load_legacy_inputs(data, cell_volume=CELL_VOLUME)
    enzyme_parameters = get_enzyme_parameters(
        mode="blast",
        blasted_parameters_csv=data / "enzyme_parameters_blasted.csv",
        predicted_kms_csv=data / "predicted_Kms.csv",
        equilibrium_constants_csv=data / "Equilibrium_constants.csv",
    )
    # B10 (ADR 0002): the legacy minimal kcats used flux = content / T_d instead of
    # content * ln2 / T_d (and T_d = 300 min). Passing T_d / ln2 gives exactly
    # content * ln2 / T_d with the legacy function.
    min_kcats = calculate_min_kcats(
        enzyme_concentrations_umol=inp.enzyme_concentrations_umol,
        cell_volume=CELL_VOLUME,
        doubling_time=DOUBLING_TIME / math.log(2),
    )
    return replace(inp, enzyme_parameters=enzyme_parameters, min_kcat_values=min_kcats)
