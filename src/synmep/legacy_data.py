# ruff: noqa: B023, PERF102, PLC0206  (verbatim legacy code, see module docstring)
"""Legacy data preprocessing, moved verbatim from the legacy `model.py`.

These functions turn the raw CSVs in `data/` into enzyme costs, enzyme
concentrations [µM], kinetic parameters and minimal kcats. They are kept
unchanged for single-legacy (faithful port, ADR 0001); known issues are listed in ADR 0002
(B6: default Km/Ki/Kd applied after µM→mM; B10: min-kcat estimate).
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd


def get_enzyme_costs(
    *,
    kugler_model_csv: Path,
    aa_costs_csv: Path,
) -> pd.DataFrame:
    costs = pd.read_csv(kugler_model_csv, delimiter=",", index_col=0)
    costs.index.name = "Amino Acid"

    enzyme_seqs = pd.read_csv(aa_costs_csv, delimiter=";", index_col=2).iloc[:, 2:]
    enzyme_seqs = enzyme_seqs.dropna()

    enzyme_costs = pd.DataFrame(
        0.0,
        columns=["AA seq", "NADPH", "ATP", "Fdred", "3PGA"],
        index=enzyme_seqs.index,
    )
    enzyme_costs["AA seq"] = enzyme_seqs["Amino acid Sequence"]
    for enzyme in enzyme_costs.index:
        for aa in enzyme_costs.loc[enzyme, "AA seq"]:  # type: ignore
            enzyme_costs.loc[enzyme, "NADPH"] += costs.loc[aa, "NADPH"]  # type: ignore
            enzyme_costs.loc[enzyme, "ATP"] += costs.loc[aa, "ATP"]  # type: ignore
            enzyme_costs.loc[enzyme, "Fdred"] += costs.loc[aa, "Fdred"]  # type: ignore
            enzyme_costs.loc[enzyme, "3PGA"] += costs.loc[aa, "3PGA"]  # type: ignore

    return enzyme_costs.drop(index=["PQS"])


def get_proteomics_results(
    *,
    proteomics_diff_csv: Path,  # Proteomics_DiffAnalysis_2WayAnova_.csv
    proteomics_poi_csv: Path,  # Proteomics_Proteins_of_interest.csv
    proteomics_concentration_csv: Path,  # Proteomics_samples_proteinconcentration.xlsx
    proteomics_dry_weight_csv: Path,  # Proteomics_Dry_weights.csv
    proteomics_dilutions_csv: Path,  # Proteomics_samples_dilutions.xlsx
    proteomics_od_wt_csv: Path,  # 20240328_Spectrum_OD_WT_proteomics.csv
    proteomics_od_sqs_csv: Path,  # 20240329_Spectrum_OD_sqs_proteomics.csv
    proteomics_od_dcas9_csv: Path,  # 20240329_Spectrum_OD_dCas9_proteomics.csv
    cell_volume: float,
    protein_content: float,
) -> tuple[pd.DataFrame, float]:
    df = (
        pd.read_csv(
            proteomics_diff_csv,
            delimiter=";",
        )
        .drop([0, 1], axis=0)
        .reset_index(drop=True)
    )

    cols = df.columns.tolist()
    cols = cols[:28] + cols[40:] + cols[28:40]
    df = df[cols]
    poi = pd.read_csv(proteomics_poi_csv, delimiter=";").drop([0, 15], axis=0)
    protein_conc = pd.read_csv(
        proteomics_concentration_csv,
        delimiter=";",
    ).reset_index(drop=True)
    dry_weights = pd.read_csv(proteomics_dry_weight_csv, delimiter=";").reset_index(
        drop=True
    )
    dilutions = pd.read_csv(
        proteomics_dilutions_csv,
        delimiter=";",
        index_col=0,
        header=0,
    ).reset_index(drop=True)

    dry_weights["Dry weight of sample [mg/mL]"] = (
        dry_weights["Difference"] / dilutions["Cell dry weight volume [mL]"]
    )
    # Measured dry weight of 2 mL sample / 2ml -> mg/mL in the sample

    abundance_columns = [col for col in df.columns if "Abundances" in col]
    # Apply transformations
    for col in abundance_columns:
        # Convert the values to float and reverse the log2 transformation
        df[col] = df[col].astype(float).apply(lambda x: 2**x)
        # Calculate the total protein abundance for normalization
        total_protein = df[col].sum()
        # Get relative abundance values
        df[col] = df[col] / total_protein

    OD_df = pd.read_csv(
        proteomics_od_wt_csv,
        delimiter=";",
        index_col=0,
        header=1,
    )
    OD_df = pd.concat(
        [
            OD_df,
            pd.read_csv(
                proteomics_od_sqs_csv,
                delimiter=";",
                index_col=0,
                header=1,
            ),
        ],
        axis=1,
    )
    OD_df = pd.concat(
        [
            OD_df,
            pd.read_csv(
                proteomics_od_dcas9_csv,
                delimiter=";",
                index_col=0,
                header=1,
            ),
        ],
        axis=1,
    )

    # FIXME: why load a second time?
    dilutions = pd.read_csv(
        proteomics_dilutions_csv,
        delimiter=";",
        index_col=0,
        header=0,
    )
    OD_df2 = pd.DataFrame(
        columns=["OD measured", "OD in sample", "OD in culture"], index=dilutions.index
    )
    OD_df2.loc[:, "OD measured"] = OD_df.loc[750].dropna()

    for i in OD_df2.index:
        OD_df2.loc[i, "OD in sample"] = (
            OD_df2.loc[i, "OD measured"]
            / dilutions.loc[i]["dilution factor for OD measurement"]
        )
        OD_df2.loc[i, "OD in culture"] = (
            OD_df2.loc[i, "OD measured"]
            / dilutions.loc[i]["dilution factor for OD measurement"]
            / (35 / dilutions.loc[i]["resuspension volume [mL]"])
        )
        # Measured OD of diluted sample / dilution factor -> OD of the sample
        # OD of the sample/(Pellet volume/resuspension volume) -> OD of the original culture (Assuming no sampling loss)
    OD_df2.reset_index(drop=True, inplace=True)
    OD_df2["Cell count"] = OD_df2["OD in sample"] * 105482164.18012363
    # Converting OD in sample to cell count in sample [cells/mL]
    protein_conc["Concentration in cell sample [ug/ml]"] = (
        protein_conc["Concentration in cell sample [ug/ml]"] / 1000
    )
    protein_conc.rename(
        columns={
            "Concentration in cell sample [ug/ml]": "Concentration in cell sample [mg/ml]"
        },
        inplace=True,
    )
    # Getting the lower detection limit:
    poi_df = df.copy()
    protein_conc["relative conc [mg Protein/mg CDW]"] = (
        protein_conc["Concentration in cell sample [mg/ml]"]
        / dry_weights["Dry weight of sample [mg/mL]"]
    )
    for i, col in enumerate(abundance_columns):
        # print(dry_weights.loc[i, "Dry weight of sample [mg/mL]"])
        poi_df[col] = (
            poi_df[col].astype(float)
            * (
                dry_weights.loc[i, "Dry weight of sample [mg/mL]"]
                * protein_content  # type: ignore
                / OD_df2.loc[i, "Cell count"].astype(float)  # type: ignore
            )
            / cell_volume
        )

        # Dry weight*0.445 -> Protein content mg/mL of sample
        # Protein content [mg/mL] / Cell count [cells/mL] -> Protein content mg/cell
        # Protein content [mg/cell]/ 4 fL (4E-15 L)/cell -> Protein concentration mg/L intracellular
    # Merge poi_df with poi to get molecular weights
    poi_df = poi_df.merge(
        poi[["UniProt Accession", "kDa"]],
        left_on="T: Accession",
        right_on="UniProt Accession",
        how="left",
    )
    detection_limit = poi_df[abundance_columns].min(axis=0).min()

    # Dividing by the molecular weights to get molar concentration
    # mg/L / kDa(kg/mol==g/mmol==mg/µmol) -> µmol/L = µM
    for col in abundance_columns:
        poi_df[col] = poi_df.apply(
            lambda row: row[col] / row["kDa"] if pd.notnull(row["kDa"]) else row[col],
            axis=1,
        )
    detection_limit_IspD = detection_limit * 2 / 24.64349

    # Drop the extra column from the merge
    poi_df.drop(columns=["UniProt Accession"], inplace=True)
    # Verify the changes
    # Extract time points from the column names
    time_points = [4, 24, 48, 96]
    strains = ["Wild type", "sqs", "dCas9"]
    for strain in strains:
        for time_point in time_points:
            poi_df[f"Average_{strain}_{time_point}"] = np.nan
            poi_df[f"Std_{strain}_{time_point}"] = np.nan
    time_point_columns = {tp: [] for tp in time_points}
    # Calculate row-wise averages and standard deviations for each time point
    for strain in strains:
        for time_point, columns in time_point_columns.items():
            strain_columns = [
                col
                for col in poi_df.columns
                if strain in col and str(time_point) in col and "Average" not in col
            ]
            if np.isnan(poi_df[f"Average_{strain}_{time_point}"].iloc[0]):
                poi_df[f"Average_{strain}_{time_point}"] = poi_df[strain_columns].mean(
                    axis=1
                )
                poi_df[f"Std_{strain}_{time_point}"] = poi_df[strain_columns].std(
                    axis=1
                )
            else:
                pass

    return poi_df, detection_limit_IspD


def get_enzyme_concentrations(
    *,
    equilibrium_constants_csv: Path,  # Equilibrium_constants.csv
    abundances: pd.DataFrame,
    detection_limit: float,
) -> dict[str, float]:
    enzyme_info = pd.read_csv(equilibrium_constants_csv, delimiter=";", header=0)
    # Identify relevant columns
    Wild_type_columns = [
        col
        for col in abundances.columns
        if any(tp in col for tp in ["F4", "F16", "F28", "F5", "F17", "F29"])
    ]
    # print(Wild_type_columns)
    Wild_type_columns.append("T: Accession")
    abundances = abundances[Wild_type_columns]
    abundances = abundances.merge(
        enzyme_info[["Uniprot accession", "Enzyme"]],
        left_on="T: Accession",
        right_on="Uniprot accession",
        how="left",
    )
    # Drop the redundant 'Uniprot accession' column after merging
    abundances = abundances.drop(columns=["Uniprot accession"])
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
    enzyme_concentrations = {}
    Wild_type_columns.remove("T: Accession")
    for e in Enzyme_names:
        mean_series = (
            abundances[abundances.loc[:, "Enzyme"] == e]
            .loc[:, Wild_type_columns]
            .mean(axis=1)
        )
        if not mean_series.empty:
            enzyme_concentrations.update({e: mean_series.iloc[0]})
        else:
            enzyme_concentrations.update({e: detection_limit * 10})
    modified_enzyme_concentrations = {}
    for e in enzyme_concentrations:
        if e.endswith(("1", "2", "3")):
            modified_key = e[:-1]
            modified_enzyme_concentrations[modified_key] = enzyme_concentrations[e]
        else:
            modified_enzyme_concentrations[e] = enzyme_concentrations[e]
    enzyme_concentrations = modified_enzyme_concentrations
    enzyme_concentrations["DXS_nonactive"] = enzyme_concentrations["DXS"] * 0.2
    enzyme_concentrations["DXS"] = enzyme_concentrations["DXS"] * 0.8
    return enzyme_concentrations


def get_enzyme_parameters(
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

    default_values = {
        "eq": 10**10,  # FIXME: justify why you choose these values
        "Km": 100,  # FIXME: justify why you choose these values
        "at": 600,  # FIXME: justify why you choose these values
        "Ki": 1000,  # FIXME: justify why you choose these values
        "Kd": 1,  # FIXME: justify why you choose these values
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


def unit_conversion(value: float, cell_volume: float) -> float:
    """
    Function for converting values in [mmol mol(chl)^-1 s^-1] to [µM s^-1]
    [mmol mol(chl)^-1 s^-1]*[mol(chl)]/[L] = mM s^-1 -> µM min⁻1
    Chlorophyll [µM intracellular] -> M intracellular -> mol
    """
    # FIXME: why aren't these model parameters?
    # FIXME: where do the numbers come from??
    return (
        value * ((7.5 / 893.509 / (105482164.18012363 * 1000) / cell_volume) * 1e-3)
    ) * 60


def calculate_min_kcats(
    *,
    enzyme_concentrations_umol: dict,
    cell_volume: float,
    doubling_time: float,
) -> dict[str, float]:
    enzyme_concentrations_umol = enzyme_concentrations_umol.copy()
    # FIXME: relate to _unit_conversion
    chla: float = (
        7.5 / 893.509 / (105482164.18012363 * 1000)
    )  # mg/OD/L -> mmol/OD/L / Cells/L/OD --> mmol/cell
    carot = 2 / 581.5565 / (105482164.18012363 * 1000)
    squalene = 0.5 / 410.73 / (105482164.18012363 * 1000)

    enzyme_concentrations_umol.pop("ThiG")
    enzyme_concentrations_umol.pop("Shc")
    enzyme_concentrations_umol.pop("DXS_nonactive")

    Total_terpenoid_concentrations = (
        (chla * 1000 * 4) + (carot * 1000 * 8) + (squalene * 1000 * 6)
    ) / cell_volume  # Products converted to mols of C5 body needed
    # doubling_time = np.log(2)/0.13*60 #Approx. maximum growth rate in turbidostat | Unit: 1/min
    MEP_flux = Total_terpenoid_concentrations / doubling_time  # Unit: µM/min
    # Calculate kcat for each enzyme
    kcat_values = {
        enzyme + "Kcat": MEP_flux / Et
        for enzyme, Et in enzyme_concentrations_umol.items()
    }
    kcat_values["ThiGKcat"] = 100  # Completely unknown concentration
    kcat_values["SqsKcat"] = (
        squalene
        * 1000
        / cell_volume
        / doubling_time
        / enzyme_concentrations_umol["Sqs"]
    )  # mmol/cell -> µmol/cell --> µmol/L intracellular --> µM/min--> µM/min substrate/µM enzyme -->
    kcat_values["ChlPKcat"] = (
        chla * 1000 / cell_volume / doubling_time / enzyme_concentrations_umol["ChlP"]
    )
    kcat_values["ChlGKcat"] = (
        chla * 1000 / cell_volume / doubling_time / enzyme_concentrations_umol["ChlG"]
    )
    kcat_values["CrtBKcat"] = (
        carot * 1000 / cell_volume / doubling_time / enzyme_concentrations_umol["CrtB"]
    )
    kcat_values["IdiKcat"] = 138  # From literature
    kcat_values["IspH1Kcat"] = kcat_values["IspHKcat"]
    kcat_values["IspH2Kcat"] = kcat_values["IspHKcat"]

    kcat_values["CrtE1Kcat"] = kcat_values["CrtEKcat"]
    kcat_values["CrtE2Kcat"] = kcat_values["CrtEKcat"]
    kcat_values["CrtE3Kcat"] = kcat_values["CrtEKcat"]
    for key in kcat_values:
        kcat_values[key] = kcat_values[key] * 2  # FIXME: explain number
    return pd.Series(kcat_values, dtype=float).to_dict()
