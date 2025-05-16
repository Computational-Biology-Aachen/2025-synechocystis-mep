from __future__ import annotations

from functools import partial
from pathlib import Path
from typing import Any, Literal, Union, cast
from warnings import warn

import numpy as np
import numpy.typing as npt
import pandas as pd
from modelbase.ode import Model
from modelbase.ode import ratefunctions as rf
from scipy.integrate import simpson
from scipy.integrate import simpson as _simpson
from sklearn.pipeline import Pipeline

from kinetics import (
    DXS_na,
    DXS_synth_f,
    IspH1,
    IspH2,
    _irradiance_depth_int,
    chl_synthesis_rate,
    fraction_PGA_derivates_G3P,
    fraction_PGA_derivates_PYR,
    growth_function,
    reversible_michaelis_menten_keq_kcat,
    reversible_michaelis_menten_keq_kcat_comp1,
    reversible_ordered_2_comp1,
    reversible_ordered_2_comp2,
    universal_hill_2_1,
    v_dxs,
    v_sqs,
)

parameters = {  # Module parameters
    "PSIItot": 0.830,  # [mmol mol(Chl)^-1] total concentration of photosystem II complexes (Zavrel2023)
    "PSItot": 3.270,  # [mmol mol(Chl)^-1] total concentration of photosystem I complexes (Zavrel2023)
    "Q_tot": 13.000,  # [mmol mol(Chl)^-1] total PHOTOACTIVE PQ concentration (Khorobrykh2020)
    "PC_tot": 1.571,  # [mmol mol(Chl)^-1] total concentration of plastocyanin (PC_ox + PC_red) (Zavrel2019)
    "Fd_tot": 3.597,  # [mmol mol(Chl)^-1] total concentration of ferredoxin (Fd_ox + Fd_red) (Moal2012)
    "NADP_tot": 26.805,  # [mmol mol(Chl)^-1] total concentration of NADP species (NADP + NADPH) (Kauny2014)
    "NAD_tot": 11.169,  # [mmol mol(Chl)^-1] total concentration of NAD species (NAD + NADH) (Tanaka2021)
    "AP_tot": 430.143,  # [mmol mol(Chl)^-1] total concentration of adenosine species (ADP + ATP) (Doello2018)
    "O2ext": 55.402,  # [mmol mol(Chl)^-1] concentration of oxygen in the surrounding medium (Kihara2014)
    "CO2ext": 3.103,  # [mmol mol(Chl)^-1] saturated concentration of CO2 in 25 °C water with ~10 ‰ Cl^- ions (Li1971)
    "F": 96.485,  # [C mmol^-1] Faraday's constant (Richardson2019)
    "R": 8.300e-03,  # [J K^-1 mmol^-1] ideal gas constant (Richardson2019)
    "T": 298.150,  # [K] temperature (set)
    "bHi": 100.000,  # [unitless] buffering constant of the thylakoid lumen (estimated)
    "bHo": 1.111e03,  # [unitless] buffering constant of the cytoplasm, assumed to be 1/f_V_lumen times larger (estimated)
    "cf_lumen": 4.613e-05,  # [mol(Chl) ml^-1] conversion factor for [mmol mol(Chl)^-1] -> [mol l^-1] for the thylakoid lumen (derived)
    "cf_cytoplasm": 4.562e-06,  # [mol(Chl) ml^-1] conversion factor for [mmol mol(Chl)^-1] -> [mol l^-1] for the cytoplasm (derived)
    # 'fCin': 100.000, # [unitless] ratio of intracellular to external CO2 concentration with activity of the CCM (Hagemann2021)
    "E0_QA": -0.140,  # [V] standard electrode potential of the reduction of PS2 plastoquinone A (Lewis2022)
    "E0_PQ": 0.533,  # [V] standard electrode potential of the reduction of free plastoquinone (Lewis2022)
    "E0_PC": 0.350,  # [V] standard electrode potential of the reduction of free plastocyanin (Lewis2022)
    "E0_P700": 0.410,  # [V] standard electrode potential of the reduction of the oxidised PS1 reaction center (Lewis2022)
    "E0_FA": -0.580,  # [V] standard electrode potential of the reduction of PS1 iron-sulfur cluster A (Lewis2022)
    "E0_Fd": -0.410,  # [V] standard electrode potential of the reduction of free ferredoxin (Lewis2022)
    "E0_NADP": -0.113,  # [V] standard electrode potential of the reduction of NADP to NADPH (Falkowski2007)
    "E0_succinate/fumarate": 0.443,  # [V] standard electrode potential of the reduction of fumarate to succinate (Falkowski2007)
    "kH0": 5.000e08,  # [s^-1] rate constant of (unregulated) excitation quenching by heat (Ebenhoh2014)
    "kHst": 1.000e09,  # [s^-1] rate constant of state transition regulated excitation quenching by heat (guess)
    "kF": 6.250e08,  # [s^-1] rate constant of excitation quenching by fluorescence (Ebenhoh2014)
    "k2": 2.500e09,  # [s^-1] rate constant of excitation quenching by photochemistry (Ebenhoh2014,Bernat2009)
    "kPQred": 250.000,  # [mol(Chl) mmol^-1 s^-1] rate constant of PQ reduction via PS2 (Matuszynska2019)
    # 'kQuench': 2.500e-05, # [mol(Chl) mmol^-1 s^-1] rate constant of PS2 quenching by Q_red (manually fitted)
    # 'kUnquench': 1.250e-03, # [mol(Chl) mmol^-1 s^-1] rate constant of PS2 unquenching by Q_ox (manually fitted)
    # 'kATPconsumption': 0.300, # [s^-1] rate constant of ATP consumption by processes other than the CBB (guess)
    # 'kNADHconsumption': 10.000, # [s^-1] rate constant of NADH consumption by processes other than NDH (guess)
    "kFlvactivation": 0.500,  # [s^-1] rate constant of Flv activation by reduced Fd (manually fitted)
    "kFlvdeactivation": 0.100,  # [s^-1] rate constant of Flv deactivation by oxidised Fd (manually fitted)
    # 'kCBBactivation': 5.000e-02, # [s^-1] rate constant of CBB activation by reduced Fd (manually fitted)
    "kCBBdeactivation": 2.500e-03,  # [s^-1] rate constant of CBB deactivation by oxidised Fd (manually fitted)
    "kPCox": 2.500e03,  # [mol(Chl) mmol^-1 s^-1] rate constant of PC oxidation via PS1 (Matuszynska2019)
    "kFdred": 2.500e05,  # [mol(Chl) mmol^-1 s^-1] rate constant of Fd reduction via PS1 (Matuszynska2019)
    # 'k_F1': 1.000, # [s^-1] rate constant of excitation quenching via fluorescence (guess)
    "k_ox1": 0.211,  # [mol(Chl)^2 mmol^-2 s^-1] rate constant of oxygen reduction via bd-type (Cyd) terminal oxidases (Ermakova2016)
    # 'k_Q': 1.755e+03, # [mol(Chl) mmol^-1 s^-1] rate constant of PC reduction by the cytochrome b6f complex (Setif2020)
    "k_NDH": 3.136,  # [mol(Chl) mmol^-1 s^-1] rate constant of PQ reduction by NDH-2 (Cooley2001)
    "k_SDH": 0.430,  # [mol(Chl) mmol^-1 s^-1] rate constant of PQ reduction by SDH (Cooley2001)
    "k_O2": 35.935,  # [mol(Chl)^2 mmol^-2 s^-1] rate constant of Fd oxidation by Flv 1/3 (Ermakova2016)
    "k_FN_fwd": 3.730,  # [mol(Chl) mmol^-1 s^-1] rate constant of NADP reduction by FNR (Kauny2014)
    "k_FN_rev": 53.364,  # [mol(Chl)^3 mmol^-3 s^-1] rate constant of reverse flux through FNR in darkness (Kauny2014)
    # 'k_pass': 1.000e-02, # [mmol mol(Chl)^-1 s^-1] rate constant of protons leaking across the thylakoid membrane per delta pH (guess)
    "k_NQ": 100.000,  # [mol(Chl) mmol^-1 s^-1] rate constant of PQ reduction by NDH (Theune2021)
    # 'k_aa': 1.156, # [mol(Chl)^2 mmol^-2 s^-1] rate constant of oxygen reduction via aa3-type (COX) terminal oxidases (Ermakova2016)
    "kRespiration": 3.506e-06,  # [mol(Chl)^3 mmol^-3 s^-1] rate constant of 3PGA oxidation and fumarate reduction by glycolysis and the TCA cycle (Ermakova2016)
    "kO2out": 4.025e03,  # [s^-1] rate constant of oxygen diffusion out of the cell (Kihara2014)
    "kCCM": 4.025e03,  # [s^-1] rate constant of CO2 diffusion into the cell, assumed to be identical to kO2out (guess)
    # 'fluo_influence': {'PS2': 1.000, 'PS1': 1.000, 'PBS': 1.250}, # [unitless] factors multiplied to the calculated fluorescence (no effect at 1) (manually fitted)
    "PBS_free": 9.000e-02,  # [unitless] fraction of unbound PBS (Zavrel2023)
    "PBS_PS1": 0.390,  # [unitless] fraction of PBS bound to PSI (Zavrel2023)
    "PBS_PS2": 0.510,  # [unitless] fraction of PBS bound to PSII (Zavrel2023)
    # 'lcf': 0.500, # [excitations photons^-1] light conversion factor (manually fitted)
    # 'kATPsynth': 10.000, # [s^-1] rate constant of ATP synthesis (fit to supply CBB)
    "Pi_mol": 1.000e-02,  # [mmol mol(Chl)^-1] molar conccentration of phosphate (Matuszynska2019)
    "DeltaG0_ATP": 30.600,  # [kJ mol^-1] energy of ATP formation (Matuszynska2019)
    "HPR": 4.667,  # [unitless] number of protons (14) passing through the ATP synthase per ATP (3) synthesized (Pogoryelov2007)
    "vOxy_max": 16.059,  # [mmol mol(Chl)^-1 s^-1] approximate Rubisco oxygenation rate (Savir2010)
    "KMATP": 24.088,  # [mmol mol(Chl)^-1] order of magnitude of the michaelis constant for ATP consumption in the CBB cycle (Wadano1998,Tsukamoto2013)
    "KMNADPH": 18.066,  # [mmol mol(Chl)^-1] approxiate michaelis constant for NADPH consumption in the CBB cycle (Koksharova1998)
    "KMCO2": 72.264,  # [mmol mol(Chl)^-1] order of magnitude of the michaelis constant for CO2 consumption by cyanobacterial Rubisco (Savir2010)
    "KIO2": 240.880,  # [mmol mol(Chl)^-1] order of magnitude of the michaelis inhibition constant of O2 for CO2 consumption by cyanobacterial Rubisco, assumed equal to KMO2 (Savir2010)
    "KMO2": 240.880,  # [mmol mol(Chl)^-1] order of magnitude of the michaelis constant for O2 consumption by cyanobacterial Rubisco (Savir2010)
    "KICO2": 72.264,  # [mmol l^-1] order of magnitude of the michaelis inhibition constant of CO2 for O2 consumption by cyanobacterial Rubisco, assumed equal to KMCO2 (Savir2010)
    # 'KMPGA': 0.100, # [mmol mol(Chl)^-1] arbitrary michaelis constant limiting oxygenation reactions for low 3PGA (guess)
    "vCBB_max": 51.301,  # [mmol mol(Chl)^-1 s^-1] approximate maximal rate of the Calvin Benson Bassham cycle (Zavrel2017)
    "kPR": 1.523e-06,  # [mol(Chl)^4 mmol^-4 s^-1] rate constant of (2-phospho)glycolate recycling into 3PGA (Huege2011)
    "fCin": 1032.3825198828176,
    "k_F1": 1.0144691550897802,
    "k_Q": 1789.6877149900638,
    "k_pass": 0.0103757688087645,
    "k_aa": 1.0726342220286855,
    "fluo_influence": {
        "PS2": 1.0871604554468057,
        "PS1": 1.0580747331180056,
        "PBS": 1.2652931291892002,
    },
    "lcf": 0.4852970468572075,
    "KMPGA": 0.108119035344566,
    "kATPsynth": 10.508685365751944,
    "kATPconsumption": 0.3007216562727887,
    "kNADHconsumption": 10.774127129643528,
    "kUnquench": 0.1084895168082126,
    "KMUnquench": 0.2057753667172297,
    "kQuench": 0.0021298894016293,
    "KHillFdred": 39.138885950985646,
    "nHillFdred": 3.803733201004677,
    "kCBBactivation": 0.0368469911771887,
    "KMFdred": 0.2867372651572032,
    "kOCPactivation": 9.119349718302952e-05,
    "kOCPdeactivation": 0.0013779872861697,
    "OCPmax": 0.2894610142573665,
    "vNQ_max": 47.52413111835536,
    "KMNQ_Qox": 1.225188500783819,
    "KMNQ_Fdred": 1.3324893590389415,
}


def _ensure_single_value(value: npt.NDArray | float) -> float:
    """Extract the last item if value is a list or tuple, otherwise return the value itself."""
    return value[0] if isinstance(value, np.ndarray) else value


def get_pigment_content() -> pd.Series:
    """[mg(Pigment) mg(Chla)^-1]
    relative pigment concentrations in a synechocystis cell (Zavrel2023)
    """
    return pd.Series(
        {
            "chla": 1.000,
            "beta_carotene": 0.176,
            "allophycocyanin": 1.118,
            "phycocyanin": 6.765,
        },
        dtype=float,
    )


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
    for e in enzyme_concentrations.keys():
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
    default_values = {"eq": 10**10, "Km": 100, "at": 600, "Ki": 1000, "Kd": 1}
    for name in Enzyme_names:
        # substrates = pd.unique(enzyme_parameters[enzyme_parameters["Enzyme"] == name]["Substrate"])
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
    for key in enzyme_parameters_dict_list.keys():
        if np.isnan(enzyme_parameters_dict_list[key]):
            enzyme_parameters_dict_list[key] = default_values[key[-2:]]

    return pd.Series(enzyme_parameters_dict_list, dtype=float).to_dict()


def unit_conversion(value: float, cell_volume: float) -> float:
    """
    Function for converting values in [mmol mol(chl)^-1 s^-1] to [µM s^-1]
    [mmol mol(chl)^-1 s^-1]*[mol(chl)]/[L] = mM s^-1 -> µM min⁻1
    Chlorophyll [µM intracellular] -> M intracellular -> mol
    """

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
    # kcat_values["IspDKcat"] = 1000
    kcat_values["IspH1Kcat"] = kcat_values["IspHKcat"]
    kcat_values["IspH2Kcat"] = kcat_values["IspHKcat"]

    kcat_values["CrtE1Kcat"] = kcat_values["CrtEKcat"]
    kcat_values["CrtE2Kcat"] = kcat_values["CrtEKcat"]
    kcat_values["CrtE3Kcat"] = kcat_values["CrtEKcat"]
    for key in kcat_values.keys():
        kcat_values[key] = kcat_values[key] * 2  # FIXME: explain number
    return pd.Series(kcat_values, dtype=float).to_dict()


# Define functions
def get_pigment_association(
    *,
    ps1_ratio: float,
    ps_comp: pd.DataFrame,
    molar_masses: pd.DataFrame,
    beta_carotene_method: str = "original",
    pigments: pd.Series | None = None,
    verbose: bool = True,
) -> pd.DataFrame:
    """Create a DataFrame containing the fractions of pigments associated to the main complexes (PS1, PS2. PBS)

    Args:
        ps1_ratio (float): ratio of PS1:PS2
        beta_carotene_method (str): how to calculate the beta carotene energy donation, one of 'original' or 'stoichiometric'
        pigments (pd.DataFrame): DataFrame with the pigment amounts [mg(pigment) mg(chlorophyll a)^-1] in named columns
        verbose (bool): issue a warning if provided and calculated beta carotene amounts don't match
    """
    # Adapt the ratio tab in _ps_comp
    ps_comp.loc["ratio", "ps1"] = ps1_ratio

    # Chlorophyll a: Fraction of total pool
    chla = chla = ps_comp.loc[["ratio", "n_chla"], :].prod(axis=0)
    chla = chla / chla.sum()

    # Beta carotene: Fraction of total pool wih adjustment factor
    beta_carotene = ps_comp.loc[["ratio", "n_beta_carotene"], :].prod(axis=0)
    beta_carotene = beta_carotene / beta_carotene.sum()

    if beta_carotene_method == "original":
        # Original in Fuente 2021: Fraction of total pool but only 75% in membrane
        beta_carotene = beta_carotene.mul(ps_comp.loc["beta_carotene_in_membrane", :])

    elif beta_carotene_method == "stoichiometric":
        if pigments is not None:
            # Calculate the molar ratio of beta carotene following from the stoichiometry
            bc_stoich = ps_comp.loc[["ratio", "n_beta_carotene"], :].prod(axis=0) * (
                1 / ps_comp.loc["ratio", :].sum()
            )
            bc_stoich = bc_stoich.mul(1 / ps_comp.loc["n_chla", :]).sum()  # type: ignore

            # Compare with the measured value and calculate a proportionality factor
            bc_measured = (
                pigments.loc["beta_carotene"]
                * molar_masses["chla"]
                / molar_masses["beta_carotene"]
            )
            bc_fac = np.min(np.array([bc_stoich / bc_measured, 1]))

            if bc_stoich > bc_measured and verbose:
                warn(
                    f"The provided pigments have a molar beta carotene:chl a ratio of {np.round(bc_measured, 3)} while the photosystem stoichiometry requires {np.round(bc_stoich, 3)}"
                )

            # Also take into account that only ps1 is excited by carotenoid absorption
            beta_carotene = beta_carotene * pd.Series({"ps1": bc_fac, "ps2": 0})
        else:
            raise ValueError(
                "with beta_carotene_method 'stoichiometric', the pigments have to be provided"
            )
    else:
        raise ValueError(
            "beta_carotene_method has to be either 'original' or 'stoichiometric'"
        )

    # Construct the full association map including pbs
    ass_map = pd.DataFrame({"chla": chla, "beta_carotene": beta_carotene})

    pbs = pd.Series([1, 1], index=["allophycocyanin", "phycocyanin"], name="pbs")
    ass_map = pd.concat([ass_map, pd.DataFrame([pbs])], axis=0).fillna(0)

    return ass_map


def get_pigment_absorption(
    *,
    pigments: pd.Series,
    abs_coef: pd.DataFrame,
) -> pd.DataFrame:
    """Get the per-pigment wavelength-specific light absorption coefficients for a certain pigment composition

    Args:
        pigments (pd.DataFrame): DataFrame with the pigment amounts [mg(pigment) mg(chlorophyll a)^-1] in named columns

    Returns:
        pd.DataFrame: a DataFrame with the chlorophyll a normalised absorption coefficients [m^2_mg(chlorophyll a)^-1]
    """
    # Multiply the per-pigment abosrption coefficients with the pigment concentration
    # Also add the new unit as a index level
    abs_coef_pigm = abs_coef.droplevel(1, axis=1).mul(pigments, axis=1)
    abs_coef_pigm = pd.concat(
        {"absorption_coefficient_m2_mgChla-1": abs_coef_pigm},
        names=["description"],
        axis=1,
    ).swaplevel(0, 1, 1)

    return abs_coef_pigm


def _get_complex_absorption(
    *,
    corrected_light: pd.Series | pd.DataFrame,
    ps_ratio: float,
    pigment_content: pd.Series,
    abs_coef: pd.DataFrame,
    molar_masses: pd.DataFrame,
    ps_comp: pd.DataFrame,
) -> pd.Series:
    """Get the absorption of each photosynthetic complex for a given input light

    Args:
        light (npt.ArrayLike | pd.DataFrame | pd.Series): mapping of light wavelengths to intensities
        absorption (pd.DataFrame): DataFrame containing the cellular, wavelength-specific absorption of each pigment
        association (pd.DataFrame): DataFrame containing the association of each pigment to the photosynthetic complexes

    Returns:
        pd.Series: _description_
    """
    absorption = get_pigment_absorption(
        pigments=pigment_content,
        abs_coef=abs_coef,
    )
    association = get_pigment_association(
        ps1_ratio=ps_ratio,
        molar_masses=molar_masses,
        ps_comp=ps_comp,
        beta_carotene_method="stoichiometric",
        pigments=pigment_content,
        verbose=False,
    )

    M_chl = 893.509  # [g mol^-1] molar mass of chlorophyll a (C55H72MgN4O5)

    if isinstance(corrected_light, pd.DataFrame):
        corrected_light = corrected_light.squeeze(np.argmin(abs_coef.shape))  # type: ignore

    abs_light = absorption.mul(corrected_light, axis=0)  # type: ignore
    abs_pigment = abs_light.apply(_simpson, axis=0)  # type: ignore

    # With the association calculate the complexes light absorption
    abs_complex = association.dot(
        pd.DataFrame({"absorption_m2_mgChla-1": abs_pigment.droplevel(1)})
    )

    # [µmol(Photons) mmol(Chl)^-1 s^-1]
    return abs_complex.squeeze(1) * M_chl


def get_mean_sample_light(
    *,
    I0: Union[pd.Series, npt.ArrayLike],
    depth: float,
    absorption_coef: pd.Series,
    chlorophyll_sample: float,
    depth0: float = 0,
) -> pd.Series:
    """Calculate the light experienced by an average cell in a sample

    Args:
        I0 (Union[pd.Series, np.ArrayLike]): initial light as specific irradiance per wavelength OR a list of such
        depth (float): depth of the sample in which the irradiance is exponentially attenuated
        absorption_coef (pd.Series): specific absorption coefficients per wavelength
        chlorophyll_sample (float): chlorophyll content of the sample
        depth0 (float, optional): depth of the initial light. Defaults to 0.

    Returns:
        pd.Series: adjusted light vector OR a list of such
    """
    # Calculate the adjustment for a single irradiace or a list of such
    if isinstance(I0, pd.Series):
        I_diff = (
            _irradiance_depth_int(I0, depth, absorption_coef, chlorophyll_sample)  #
            - _irradiance_depth_int(I0, depth0, absorption_coef, chlorophyll_sample)
        )
        return I_diff / (depth - depth0)

    elif isinstance(I0, (list, np.ndarray)):
        if np.all([isinstance(i, pd.Series) for i in I0]):
            I_diff = [
                (
                    _irradiance_depth_int(
                        Icurr, depth, absorption_coef, chlorophyll_sample
                    )
                    - _irradiance_depth_int(
                        Icurr, depth0, absorption_coef, chlorophyll_sample
                    )
                )
                / (depth - depth0)
                for Icurr in I0
            ]
            return I_diff  # type: ignore

    raise ValueError("I0 must be a pd.Series or a list of such")


def light_spectra(
    *,
    lights: pd.DataFrame,
    which: str | None = None,
    intensity: float = 1.0,
) -> npt.ArrayLike:
    """Get the light spectrum of a common light source in the range of 400 nm - 700 nm with a set integrated intensity

    Args:
        lights:
        which (str): name of the spectrum to be returned. One of 'solar', 'halogen_lamp', 'cool_white_led', 'warm_white_led', 'fluorescent_lamp', 'incandescent_bulb'
        intensity (float, optional): integrated light intensity. Defaults to 1.

    Returns:
        npt.ArrayLike: an array of wavelength-specific light intensities
    """
    if which is None:
        print("possible spectra:\n", sorted(list(lights.columns)))
        return np.array([])
    try:
        return lights.loc[:, which] * intensity
    except KeyError as e:
        raise KeyError(f"'{e}' is not a provided spectrum")


def _get_model_inputs(
    *,
    cell_density: float,  # [cells ml^-1]
    chlorophyll: float,  # [µmol l^-1]
    carotenoids: float,  # [µmol l^-1]
    light_intensity: float,  # Model
    phycocyanin: float,  # [mg l^-1]
    allophycocyanin: float,  # [mg l^-1]
    sample_depth_m: float,  # [m] Assuming a cuvette with 1 cm diameter
    beta_carotene_fraction: float,  # [rel] fraction of beta-carotene of cellular carotenoids
    cell_volume: float,  # [l]
    abs_coef: pd.DataFrame,
    lights: pd.DataFrame,
):
    # Calculate the chlorophyll content

    molar_weights = {  # [g mol^-1]
        "chlorophyll": 893.509,  # source: https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:189419
        "carotenoids": 536.8726,  # molar weight beta carotene, source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7235407
    }

    chlorophyll = _ensure_single_value(chlorophyll)
    carotenoids = _ensure_single_value(carotenoids)
    phycocyanin = _ensure_single_value(phycocyanin)
    allophycocyanin = _ensure_single_value(allophycocyanin)
    cell_density = _ensure_single_value(cell_density)
    mg_chlorophyll = chlorophyll * molar_weights["chlorophyll"] / 1000  # [mg l^-1]

    # Calculate the relative contents of carotenoids, phycocyanin, and allophycocyanin
    relative_carotenoids = (
        carotenoids * molar_weights["carotenoids"] / 1000
    ) / mg_chlorophyll  # [mg mg(Chla)^-1]
    relative_phycocyanin = phycocyanin / mg_chlorophyll  # [mg mg(Chla)^-1]
    relative_allophycocyanin = allophycocyanin / mg_chlorophyll  # [mg mg(Chla)^-1]

    # Get the pigment content of the cell
    pigment_content = pd.Series(
        {
            "chla": 1,  # This is always 1
            "beta_carotene": relative_carotenoids,
            "phycocyanin": relative_phycocyanin,
            "allophycocyanin": relative_allophycocyanin,
        }
    )

    # Only beta carotene is photosynthetically active
    photopigment_content = pigment_content.copy()
    photopigment_content["beta_carotene"] = (
        photopigment_content["beta_carotene"] * beta_carotene_fraction
    )

    # Get the total cellular absorption
    absorption = get_pigment_absorption(pigments=pigment_content, abs_coef=abs_coef)
    absorption = absorption.sum(axis=1)

    # Define the input light
    light = light_spectra(
        lights=lights,
        which="warm_white_led",
        intensity=light_intensity,
    )  # The model light intensity is µmol(photons) m^-2 s^-1, correct?

    # Calculate the chlorophyll concentration in the sample [mg(Chl) m^-3]
    chlorophyll_sample = (
        chlorophyll  # [µmol l^-1] intracellular
        * cell_volume  # [l cell^-1]
        * cell_density  # [cells ml^-1]
        * 1e6  # [ml m^-3]
        * 1e-3  # [mmol µmol^-1]
        * 893.509  # [g mol^-1]
    )
    # Correct the light for absorption
    corrected_light = get_mean_sample_light(
        I0=light,
        depth=sample_depth_m,  # [m]
        absorption_coef=absorption,  # wavelength-specific absorption coefficients for whole cell, summed up from pigments
        chlorophyll_sample=chlorophyll_sample,  # [mg(Chl) m^-3]
    )
    return {
        "orig_pfd": light,
        "absorption_spectrum": absorption,
        "pfd": corrected_light,
        "pigment_content": photopigment_content,
    }


def _get_simulated_ssrates(
    *,
    corrected_light: pd.Series,
    model: Any,
    ocp_absorption: pd.DataFrame,
    ps_ratio: float,  # Value in the default model
    pigment_content: pd.Series,
    output_rates: list[str],
    abs_coef: pd.DataFrame,
    molar_masses: pd.DataFrame,
    ps_comp: pd.DataFrame,
) -> pd.Series:
    # Calculate the absorbed lights from the light input
    complex_abs = _get_complex_absorption(
        corrected_light=corrected_light,
        ps_ratio=ps_ratio,
        pigment_content=pigment_content,
        abs_coef=abs_coef,
        molar_masses=molar_masses,
        ps_comp=ps_comp,
    )
    light_params = {f"complex_abs_{k}": v for k, v in complex_abs.items()}

    # Calculate the OCP absorbed light
    light_params["light_ocp"] = (
        ocp_absorption.mul(corrected_light, axis=0).apply(simpson).to_numpy()  # type: ignore
    )

    # get the input parameters for the model
    model_input = light_params.copy()

    # Include the pigment content
    pigment_input = pigment_content.loc[
        ["phycocyanin", "allophycocyanin", "beta_carotene"]
    ].to_dict()
    pigment_input = {f"pigment_{k}": v for k, v in pigment_input.items()}
    model_input.update(pigment_input)

    # return model_input
    model_input = pd.DataFrame(model_input, index=[1])

    # Create a container for the results
    return pd.Series(
        model.predict(model_input).flatten(),
        index=pd.Index(output_rates, name="production_rate"),
    )


def get_influx_rate_estimations(
    light_intensity: float,  # Model
    cell_density: float,  # [cells ml^-1]
    chlorophyll: float,  # [µmol l^-1]
    carotenoids: float,  # [µmol l^-1]
    phycocyanin: float,  # [mg l^-1]
    allophycocyanin: float,  # [mg l^-1]
    ps_ratio: float,
    beta_carotene_fraction: float,  # [rel] fraction of beta-carotene of cellular carotenoids
    sample_depth_m: float,  # [m] Assuming a cuvette with 1 cm diameter
    cell_volume: float,  # [l]
    # partial arguments
    ocp_absorption: pd.DataFrame,
    molar_masses: pd.DataFrame,
    ps_comp: pd.DataFrame,
    biomass_regression_model: Pipeline,  # The predictor model
    abs_coef: pd.DataFrame,
    lights: pd.DataFrame,
):
    output_rates = ["ATP", "NADPH", "3PGA", "Fd_red"]

    if len(chlorophyll) > 1:  # type: ignore
        rates = []
        for i, v in enumerate(chlorophyll):  # type: ignore
            pred_input = _get_model_inputs(
                cell_density=cast(
                    float,
                    np.array([cell_density[i] * 1e6]),  # type: ignore
                ),  # Conversion from [cells nL⁻1] to [cells ml^-1]
                chlorophyll=cast(
                    float,
                    np.array([chlorophyll[i] * 1000]),  # type: ignore
                ),  # Conversion from [mmol l⁻1] to [µmol l^-1]
                carotenoids=cast(
                    float,
                    np.array([carotenoids[i] * 1000]),  # type: ignore
                ),  # Conversion from [mmol l⁻1] to [µmol l^-1]
                phycocyanin=phycocyanin,  # [mg l^-1]
                allophycocyanin=allophycocyanin,  # [mg l^-1]
                light_intensity=light_intensity,  # Model
                sample_depth_m=sample_depth_m,  # [m] Assuming a cuvette with 1 cm diameter
                beta_carotene_fraction=beta_carotene_fraction,  # [rel] fraction of beta-carotene of cellular carotenoids
                cell_volume=cell_volume,  # [l]
                abs_coef=abs_coef,
                lights=lights,
            )
            rates.append(
                _get_simulated_ssrates(
                    corrected_light=cast(pd.Series, pred_input["pfd"]),
                    pigment_content=cast(pd.Series, pred_input["pigment_content"]),
                    ps_ratio=ps_ratio,  # Value in the default model
                    model=biomass_regression_model,
                    output_rates=output_rates,
                    ocp_absorption=ocp_absorption,
                    abs_coef=abs_coef,
                    molar_masses=molar_masses,
                    ps_comp=ps_comp,
                )
            )
        rows = [[], [], [], []]
        # Iterate over each Series in 'rates'
        for series in rates:
            # Append each row's value from the current Series to the respective list in 'rows'
            for i in range(4):
                rows[i].append(series.iloc[i])
        return rows

    # Get the inputs into the predictor function
    pred_input = _get_model_inputs(
        cell_density=cell_density * 1e6,  # [cells ml^-1]
        chlorophyll=chlorophyll[0] * 1000,  # type: ignore , [µmol l^-1]
        carotenoids=carotenoids[0] * 1000,  # type: ignore , [µmol l^-1]
        phycocyanin=phycocyanin,  # [mg l^-1]
        allophycocyanin=allophycocyanin,  # [mg l^-1]
        light_intensity=light_intensity,  # Model
        sample_depth_m=sample_depth_m,  # [m] Assuming a cuvette with 1 cm diameter
        beta_carotene_fraction=beta_carotene_fraction,  # [rel] fraction of beta-carotene of cellular carotenoids
        cell_volume=cell_volume,  # [l]
        abs_coef=abs_coef,
        lights=lights,
    )
    # Predict the simulated rates
    rates = _get_simulated_ssrates(
        corrected_light=cast(pd.Series, pred_input["pfd"]),
        pigment_content=cast(pd.Series, pred_input["pigment_content"]),
        ps_ratio=ps_ratio,  # Value in the default model
        model=biomass_regression_model,
        output_rates=output_rates,
        ocp_absorption=ocp_absorption,
        abs_coef=abs_coef,
        molar_masses=molar_masses,
        ps_comp=ps_comp,
    )
    rates.loc[rates < 0] = 0
    return rates


def get_model(
    growth_mode: Literal["dynamic", "steady"],
    strain: Literal["wt", "dshc"],
    enzyme_total_cost_mol: pd.DataFrame,
    enzyme_concentrations_umol: dict[str, float],
    enzyme_parameters: dict,
    cell_volume: float,
    biomass_regression_model: Any,
    abs_coef: pd.DataFrame,
    lights: pd.DataFrame,
    ocp_absorption: pd.DataFrame,
    molar_masses: pd.DataFrame,
    ps_comp: pd.DataFrame,
    ps_ratio: float = 5.9,
) -> Model:
    enzyme_names: list[str] = [
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

    m = Model()
    m.add_compounds(
        [
            "3PGA",
            "ATP",
            "NADPH",
            "Fdred",
            "DXP",
            "ThiDP",
            "MEP",
            "CDP_ME",
            "CDP_MEP",
            "MEcPP",
            "HMB_PP",
            "DMAPP",
            "IPP",
        ]
    )  # MEP pathway
    m.add_compounds(
        [
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
    )  # Terpenoid synthesis and products
    m.add_compounds(
        [
            "DXS",
            "DXS_nonactive",
            "ThiG",
            "DXR",
            "IspD",
            "IspE",
            "IspF",
            "IspG",
            "IspH",
            "Idi",
            "CrtE",
            "Sqs",
            "Shc",
            "ChlP",
            "CrtB",
            "ChlG",
        ]
    )  # Enzymes

    enzyme_synthesis_rates = {}
    for i, (d, value) in enumerate(enzyme_concentrations_umol.items()):
        enzyme_synthesis_rates[d + "_ksynth"] = (
            value / (5.5 * 60) * 2
        )  # Vmax = doubling of enzyme concentration every 5.5h at maximum growth rate
        enzyme_synthesis_rates[d + "_synth_KmA"] = 1
        enzyme_synthesis_rates[d + "_synth_KmB"] = 1
        enzyme_synthesis_rates[d + "_synth_KmC"] = 1
        enzyme_synthesis_rates[d + "_synth_KmD"] = 1
        enzyme_synthesis_rates[d + "_synth_h"] = 1

    all_compounds = m.compounds.copy()
    all_compounds.append("G3P")
    all_compounds.append("PYR")

    possible_params = ["Km", "Ki", "Kd"]
    existing_params = {}
    for enzyme in enzyme_names:
        existing_params[enzyme + "Kcat"] = enzyme_parameters[enzyme + "Kcat"]
        existing_params[enzyme + "_Keq"] = enzyme_parameters[enzyme + "_Keq"]
        existing_params[enzyme + "_h"] = 1
        for substrates in all_compounds:
            for p in possible_params:
                try:
                    existing_params[enzyme + substrates + p] = enzyme_parameters[
                        enzyme + substrates + p
                    ]
                except:  # noqa: E722
                    pass
    m.add_and_update_parameters(
        {
            "light": 100.0,
            "ATP_out_metabolism": 0.5,
            "3PGA_out_metabolism": 0.5,
            "NADPH_out_metabolism": 0.5,
            "Fdred_out_metabolism": 0.5,
            "Kna_DXS": 0.000004,
            "DXS_na_kdeg": 1 / 50000,
            "DXS_inh_str": 0.1,
            "ps_ratio": ps_ratio,
            "beta_carotene_fraction": 0.26,
            "sample_depth_m": 0.01,
            "cell_volume": cell_volume,
            "NADPHsink_kfwd": 0.1,
            "NADPHsink_krev": 0.1,
            "Fdredsink_kfwd": 0.1,
            "Fdredsink_krev": 0.1,
            "phycocyanin": 53180
            * 2e-15
            / cell_volume,  # mg/L intracellular, adjusted from default cell volume to currently used cell volume
            "allophycocyanin": 24076
            * 2e-15
            / cell_volume,  # mg/L intracellular, adjusted from default cell volume to currently used cell volume
        }
    )
    m.add_and_update_parameters(existing_params)
    m.add_and_update_parameters(enzyme_synthesis_rates)

    m.add_derived_compound(
        "Chla_synthesis",
        chl_synthesis_rate,
        args=[
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
        ],
    )
    m.add_derived_compound(
        "Growth_rate",
        growth_function,
        args=["Chla_synthesis", "cell_volume"],
    )

    m.add_compound("Cell_density")
    if growth_mode == "dynamic":
        m.add_reaction_from_args(
            "Cell_density_increase",
            rf.mass_action_1,
            stoichiometry={"Cell_density": 1},
            args=["Cell_density", "Growth_rate"],
        )

    for compound in sorted(set(m.compounds) ^ {"Cell_density"}):
        m.add_reaction_from_args(
            str(compound) + "_out",
            rf.mass_action_1,
            stoichiometry={compound: -1},
            args=[compound, "Growth_rate"],
        )

    m.add_reaction_from_args(
        "Biomass",
        growth_function,
        stoichiometry={
            "NADPH": -8.762903549999999 / 60,
            "ATP": -33.0939284280016 / 60,
            "Fdred": -25.97602870400234 / 60,
            "3PGA": -78.47848943300428 / 60,
        },
        args=["Chla_synthesis", "cell_volume"],
    )

    function = partial(
        get_influx_rate_estimations,
        ocp_absorption=ocp_absorption,
        molar_masses=molar_masses,
        ps_comp=ps_comp,
        biomass_regression_model=biomass_regression_model,
        abs_coef=abs_coef,
        lights=lights,
    )
    # FIXME: stupid mb1 lambda thingy
    function.__name__ = "production_rate"  # type: ignore
    m.add_algebraic_module_from_args(
        "production_rates",
        function=function,  # type: ignore
        derived_compounds=[
            "vATP_prod",
            "vNADPH_prod",
            "v3PGA_prod",
            "vFdred_prod",
        ],
        args=[
            "light",
            "Cell_density",
            "chlorophyll-phytol",
            "carotenoids",
            "phycocyanin",
            "allophycocyanin",
            "ps_ratio",
            "beta_carotene_fraction",
            "sample_depth_m",
            "cell_volume",
        ],
    )

    for comp in ["ATP", "NADPH", "3PGA", "Fdred"]:
        m.add_reaction_from_args(
            f"{comp} production",
            function=unit_conversion,
            stoichiometry={comp: 1},
            args=[f"v{comp}_prod", "cell_volume"],
        )

    m.add_derived_compound(
        "G3P",
        fraction_PGA_derivates_G3P,
        args=["3PGA"],
    )
    m.add_derived_compound(
        "PYR",
        fraction_PGA_derivates_PYR,
        args=["3PGA"],
    )

    m.add_reaction_from_args(
        "3PGA_out_metabolism",
        rf.mass_action_1,
        stoichiometry={"3PGA": -1},
        args=["3PGA", "3PGA_out_metabolism"],
    )
    m.add_reaction_from_args(
        "ATP_out_metabolism",
        rf.mass_action_1,
        stoichiometry={"ATP": -1},
        args=["ATP", "ATP_out_metabolism"],
    )
    m.add_reaction_from_args(
        "NADPH_out_metabolism",
        rf.mass_action_1,
        stoichiometry={"NADPH": -1},
        args=["NADPH", "NADPH_out_metabolism"],
    )
    m.add_reaction_from_args(
        "Fdred_out_metabolism",
        rf.mass_action_1,
        stoichiometry={"Fdred": -1},
        args=["Fdred", "Fdred_out_metabolism"],
    )

    o = enzyme_total_cost_mol.index.delete(0)
    sf = 1e-4  # Scaling factor to prevent giant stoichiometries for enzymes
    for enzyme in o:
        m.add_reaction_from_args(
            enzyme + "_synth",
            rf.ping_pong_4,
            stoichiometry={
                "NADPH": sf * -enzyme_total_cost_mol.loc[enzyme, "NADPH"],  # type: ignore
                "ATP": sf * -enzyme_total_cost_mol.loc[enzyme, "ATP"],  # type: ignore
                "Fdred": sf * -enzyme_total_cost_mol.loc[enzyme, "Fdred"],  # type: ignore
                "3PGA": sf * -enzyme_total_cost_mol.loc[enzyme, "3PGA"],  # type: ignore
                enzyme: sf * 1000,
            },
            args=[
                "ATP",
                "NADPH",
                "Fdred",
                "3PGA",
                enzyme + "_ksynth",
                enzyme + "_synth_KmA",
                enzyme + "_synth_KmB",
                enzyme + "_synth_KmC",
                enzyme + "_synth_KmD",
            ],
        )
    m.add_reaction_from_args(
        "DXS_synth",
        DXS_synth_f,
        stoichiometry={
            "NADPH": sf * -enzyme_total_cost_mol.loc["DXS", "NADPH"],  # type: ignore
            "ATP": sf * -enzyme_total_cost_mol.loc["DXS", "ATP"],  # type: ignore
            "Fdred": sf * -enzyme_total_cost_mol.loc["DXS", "Fdred"],  # type: ignore
            "3PGA": sf * -enzyme_total_cost_mol.loc["DXS", "3PGA"],  # type: ignore
            "DXS": sf * 1000,
        },
        args=[
            "ATP",
            "NADPH",
            "Fdred",
            "3PGA",
            "DXS_ksynth",
            "DXS_synth_KmA",
            "DXS_synth_KmB",
            "DXS_synth_KmC",
            "DXS_synth_KmD",
            "DXS_nonactive",
            "DXS_inh_str",
        ],
    )

    m.add_reaction_from_args(
        "DXS_monomerization",
        DXS_na,
        stoichiometry={"DXS": -1, "DXS_nonactive": 1},
        args=["DXS", "Kna_DXS", "IPP", "DMAPP", "DXSIPPKi", "DXSDMAPPKi"],
    )

    m.add_reaction_from_args(
        "DXS_degradation",
        rf.mass_action_1,
        stoichiometry={"DXS_nonactive": -1},
        args=["DXS_nonactive", "DXS_na_kdeg"],
    )

    m.add_reaction_from_args(
        "DXS",
        v_dxs,
        stoichiometry={"3PGA": -2, "DXP": 1},
        args=[
            "PYR",
            "G3P",
            "DXS",
            "DXSKcat",
            "DXSPYRKm",
            "DXSG3PKm",
            "DXSPYRKi",
            "IPP",
            "DMAPP",
            "ThiDP",
            "DXSIPPKi",
            "DXSDMAPPKi",
            "DXSThiDPKd",
        ],
    )

    m.add_reaction_from_args(
        "ThiG",
        reversible_michaelis_menten_keq_kcat,
        stoichiometry={"DXP": -1, "ATP": -1, "ThiDP": 1},
        args=[
            "DXP",
            "ThiDP",
            "ThiGKcat",
            "ThiG",
            "ThiGDXPKm",
            "ThiGThiDPKm",
            "ThiG_Keq",
        ],
    )
    m.add_reaction_from_args(
        "DXR",
        universal_hill_2_1,
        stoichiometry={"DXP": -1, "NADPH": -1, "MEP": 1},
        args=[
            "DXP",
            "NADPH",
            "MEP",
            "DXRKcat",
            "DXR",
            "DXRDXPKm",
            "DXRNADPHKm",
            "DXRMEPKm",
            "DXR_Keq",
            "DXR_h",
        ],
    )
    m.add_reaction_from_args(
        "IspD",
        universal_hill_2_1,
        stoichiometry={"MEP": -1, "ATP": -2, "CDP_ME": 1},
        args=[
            "MEP",
            "ATP",
            "CDP_ME",
            "IspDKcat",
            "IspD",
            "IspDMEPKm",
            "IspDATPKm",
            "IspDCDP_MEKm",
            "IspD_Keq",
            "IspD_h",
        ],
    )
    m.add_reaction_from_args(
        "IspE",
        universal_hill_2_1,
        stoichiometry={"CDP_ME": -1, "ATP": -1, "CDP_MEP": 1},
        args=[
            "CDP_ME",
            "ATP",
            "CDP_MEP",
            "IspEKcat",
            "IspE",
            "IspECDP_MEKm",
            "IspEATPKm",
            "IspECDP_MEPKm",
            "IspE_Keq",
            "IspE_h",
        ],
    )
    m.add_reaction_from_args(
        "IspF",
        reversible_michaelis_menten_keq_kcat,
        stoichiometry={"CDP_MEP": -1, "MEcPP": 1},
        args=[
            "CDP_MEP",
            "MEcPP",
            "IspFKcat",
            "IspF",
            "IspFCDP_MEPKm",
            "IspFMEcPPKm",
            "IspF_Keq",
        ],
    )
    m.add_reaction_from_args(
        "IspG",
        universal_hill_2_1,
        stoichiometry={"MEcPP": -1, "Fdred": -1, "HMB_PP": 1},
        args=[
            "MEcPP",
            "Fdred",
            "HMB_PP",
            "IspGKcat",
            "IspG",
            "IspGMEcPPKm",
            "IspGFdredKm",
            "IspGHMB_PPKm",
            "IspG_Keq",
            "IspG_h",
        ],
    )

    m.add_derived_compound("IspH1", IspH1, args=["IspH"])
    m.add_derived_compound("IspH2", IspH2, args=["IspH"])

    m.add_reaction_from_args(
        "IspH1",
        universal_hill_2_1,
        stoichiometry={"HMB_PP": -1, "NADPH": -1, "IPP": 1},
        args=[
            "HMB_PP",
            "NADPH",
            "IPP",
            "IspH1Kcat",
            "IspH1",
            "IspH1HMB_PPKm",
            "IspH1NADPHKm",
            "IspH1IPPKm",
            "IspH1_Keq",
            "IspH1_h",
        ],
    )
    m.add_reaction_from_args(
        "IspH2",
        universal_hill_2_1,
        stoichiometry={"HMB_PP": -1, "NADPH": -1, "DMAPP": 1},
        args=[
            "HMB_PP",
            "NADPH",
            "DMAPP",
            "IspH2Kcat",
            "IspH2",
            "IspH2HMB_PPKm",
            "IspH2NADPHKm",
            "IspH2DMAPPKm",
            "IspH2_Keq",
            "IspH2_h",
        ],
    )

    m.add_reaction_from_args(
        "Idi",
        reversible_michaelis_menten_keq_kcat,
        stoichiometry={"IPP": -1, "DMAPP": 1},
        args=["IPP", "DMAPP", "IdiKcat", "Idi", "IdiIPPKm", "IdiDMAPPKm", "Idi_Keq"],
    )

    m.add_reaction_from_args(
        "CrtE1",
        reversible_ordered_2_comp2,
        stoichiometry={"DMAPP": -1, "IPP": -1, "GPP": 1},
        args=[
            "DMAPP",
            "IPP",
            "GPP",
            "GPP",
            "FPP",
            "CrtE1Kcat",
            "CrtE",
            "CrtE1DMAPPKm",
            "CrtE1GPPKm",
            "CrtE1FPPKm",
            "CrtE1IPPKm",
            "CrtE1GPPKm",
            "CrtE1DMAPPKi",
            "CrtE1_Keq",
        ],
    )
    m.add_reaction_from_args(
        "CrtE2",
        reversible_ordered_2_comp2,
        stoichiometry={"GPP": -1, "IPP": -1, "FPP": 1},
        args=[
            "GPP",
            "IPP",
            "FPP",
            "DMAPP",
            "FPP",
            "CrtE2Kcat",
            "CrtE",
            "CrtE2GPPKm",
            "CrtE2DMAPPKm",
            "CrtE2FPPKm",
            "CrtE2IPPKm",
            "CrtE2FPPKm",
            "CrtE2GPPKi",
            "CrtE2_Keq",
        ],
    )
    m.add_reaction_from_args(
        "CrtE3",
        reversible_ordered_2_comp2,
        stoichiometry={"FPP": -1, "IPP": -1, "GGPP": 1},
        args=[
            "FPP",
            "IPP",
            "GGPP",
            "DMAPP",
            "GPP",
            "CrtE3Kcat",
            "CrtE",
            "CrtE3FPPKm",
            "CrtE3DMAPPKm",
            "CrtE3GPPKm",
            "CrtE3IPPKm",
            "CrtE3GGPPKm",
            "CrtE3FPPKi",
            "CrtE3_Keq",
        ],
    )

    m.add_reaction_from_args(
        "Sqs",
        v_sqs,
        stoichiometry={"FPP": -2, "NADPH": -1, "SQ": 1},
        args=[
            "FPP",
            "FPP",
            "NADPH",
            "SQ",
            "SqsKcat",
            "Sqs",
            "SqsFPPKm",
            "SqsFPPKm",
            "SqsNADPHKm",
            "SqsSQKm",
            "Sqs_Keq",
            "Sqs_h",
        ],
    )
    if strain == "wt":
        m.add_reaction_from_args(
            "Shc",
            reversible_michaelis_menten_keq_kcat,
            stoichiometry={"SQ": -1, "Hopanoids": 1},
            args=[
                "SQ",
                "Hopanoids",
                "ShcKcat",
                "Shc",
                "ShcSQKm",
                "ShcHopanoidsKm",
                "Shc_Keq",
            ],
        )

    m.add_reaction_from_args(
        "carotenoids",
        reversible_michaelis_menten_keq_kcat,
        stoichiometry={
            "GGPP": -2,
            "carotenoids": 1,
        },  # Arbitrary reaction encompassing all carotenoids
        args=[
            "GGPP",
            "carotenoids",
            "CrtBKcat",
            "CrtB",
            "CrtBGGPPKm",
            "CrtBcarotenoidsKm",
            "CrtB_Keq",
        ],
    )
    m.add_reaction_from_args(
        "phytol_synthase1",
        reversible_ordered_2_comp1,
        stoichiometry={"GGPP": -1, "NADPH": -3, "phytol": 1},  # phytol synthesis
        args=[
            "GGPP",
            "NADPH",
            "phytol",
            "chlorophyll-GGPP",
            "ChlPKcat",
            "ChlP",
            "ChlPGGPPKm",
            "ChlPchlorophyll-GGPPKm",
            "ChlPNADPHKm",
            "ChlPphytolKm",
            "ChlPGGPPKi",
            "ChlP_Keq",
        ],
        # Km = 35 µM for GGPP, from Erwinia herbicola
    )
    m.add_reaction_from_args(
        "phytol_synthase2",
        reversible_ordered_2_comp1,
        stoichiometry={
            "chlorophyll-GGPP": -1,
            "NADPH": -3,
            "chlorophyll-phytol": 1,
        },  # phytol synthesis
        args=[
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
        ],
        # Km = 35 µM for GGPP, from Erwinia herbicola
    )
    m.add_reaction_from_args(
        "chlorophyll_synthesis1",
        reversible_michaelis_menten_keq_kcat_comp1,
        stoichiometry={"phytol": -1, "chlorophyll-phytol": 1},
        args=[
            "phytol",
            "chlorophyll-phytol",
            "ChlGKcat",
            "ChlG",
            "ChlGphytolKm",
            "ChlGchlorophyll-phytolKm",
            "ChlG_Keq",
            "GGPP",
            "ChlGGGPPKm",
        ],
    )
    m.add_reaction_from_args(
        "chlorophyll_synthesis2",
        reversible_michaelis_menten_keq_kcat_comp1,
        stoichiometry={"GGPP": -1, "chlorophyll-GGPP": 1},
        args=[
            "GGPP",
            "chlorophyll-GGPP",
            "ChlGKcat",
            "ChlG",
            "ChlGGGPPKm",
            "ChlGchlorophyll-GGPPKm",
            "ChlG_Keq",
            "phytol",
            "ChlGphytolKm",
        ],
    )
    m.add_reaction_from_args(
        "NADPH_sink_Fdred",
        rf.reversible_mass_action_1_1,
        stoichiometry={"ATP": -1, "Fdred": 1},
        args=["ATP", "Fdred", "Fdredsink_kfwd", "Fdredsink_krev"],
    )
    m.add_reaction_from_args(
        "NADPH_sink_ATP",
        rf.reversible_mass_action_1_1,
        stoichiometry={"ATP": -1, "NADPH": 1},
        args=["ATP", "NADPH", "NADPHsink_kfwd", "NADPHsink_krev"],
    )
    return m
