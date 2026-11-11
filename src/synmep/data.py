"""Inputs of the legacy model, loaded exactly as in the legacy notebook (cell 1).

Everything returned here is plain data (floats, dicts, pandas objects). Nothing is
evaluated inside the ODE right-hand side; the spectral data are turned into
explicit expressions in `synmep.surrogates.petc_legacy`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from synmep.legacy_data import (
    calculate_min_kcats,
    get_enzyme_concentrations,
    get_enzyme_costs,
    get_enzyme_parameters,
    get_proteomics_results,
)

DATA = Path(__file__).resolve().parents[2] / "data"

# Legacy constants (ADR 0002, docs/assumptions.md); values unchanged for v0.
CELL_VOLUME = 2e-15  # [L cell^-1]
PROTEIN_CONTENT = 0.445  # [g protein g(DW)^-1]
CELLS_PER_OD_ML = 105482164.18012363  # [cells mL^-1 OD^-1]
M_CHL = 893.509  # [g mol^-1]
CHLA = 7.5 / M_CHL / (CELLS_PER_OD_ML * 1000)  # [mmol cell^-1]
CAROT = 2 / 581.5565 / (CELLS_PER_OD_ML * 1000)  # [mmol cell^-1]
SQUALENE = 0.5 / 410.73 / (CELLS_PER_OD_ML * 1000)  # [mmol cell^-1]
DOUBLING_TIME = 300.0  # [min]
DAY_IN_MIN = 60 * 24
PIGMENTS = ("chla", "beta_carotene", "allophycocyanin", "phycocyanin")


@dataclass(frozen=True)
class LegacyInputs:
    enzyme_costs: pd.DataFrame
    enzyme_concentrations_umol: dict[str, float]
    enzyme_parameters: dict[str, float]
    min_kcat_values: dict[str, float]
    abs_coef: pd.DataFrame  # [m^2 mg(pigment)^-1], one column per pigment
    light_spectrum: pd.Series  # warm white LED, relative [nm^-1]
    ocp_absorption: pd.Series
    ps_comp: pd.DataFrame
    molar_masses: pd.Series
    exp_growth_rate_by_density: pd.Series  # [h^-1] indexed by [cells mL^-1]
    proteins_of_interest: pd.DataFrame
    surrogate_poly: dict[str, np.ndarray]


def _spectrum(path: Path) -> pd.Series:
    return pd.read_csv(path, index_col=0).squeeze(axis=1)


def load_legacy_inputs(
    data: Path = DATA, cell_volume: float = CELL_VOLUME
) -> LegacyInputs:
    enzyme_costs = get_enzyme_costs(
        kugler_model_csv=data / "Amino_acid_costs_Kuglermodel.csv",
        aa_costs_csv=data / "Amino_acid_costs_arnold.csv",
    )
    abundances, detection_limit = get_proteomics_results(
        proteomics_diff_csv=data / "Proteomics_DiffAnalysis_2WayAnova_.csv",
        proteomics_poi_csv=data / "Proteomics_Proteins_of_interest.csv",
        proteomics_concentration_csv=data
        / "Proteomics_samples_proteinconcentration.csv",
        proteomics_dry_weight_csv=data / "Proteomics_Dry_weights.csv",
        proteomics_dilutions_csv=data / "Proteomics_samples_dilutions.csv",
        proteomics_od_wt_csv=data / "20240328_Spectrum_OD_WT_proteomics.csv",
        proteomics_od_sqs_csv=data / "20240329_Spectrum_OD_sqs_proteomics.csv",
        proteomics_od_dcas9_csv=data / "20240329_Spectrum_OD_dCas9_proteomics.csv",
        cell_volume=cell_volume,
        protein_content=PROTEIN_CONTENT,
    )
    enzyme_concentrations = get_enzyme_concentrations(
        equilibrium_constants_csv=data / "Equilibrium_constants.csv",
        abundances=abundances,
        detection_limit=detection_limit,
    )
    enzyme_parameters = get_enzyme_parameters(
        mode="blast",
        blasted_parameters_csv=data / "enzyme_parameters_blasted.csv",
        predicted_kms_csv=data / "predicted_Kms.csv",
        equilibrium_constants_csv=data / "Equilibrium_constants.csv",
    )
    min_kcats = calculate_min_kcats(
        enzyme_concentrations_umol=enzyme_concentrations,
        cell_volume=cell_volume,
        doubling_time=DOUBLING_TIME,
    )
    abs_coef = pd.DataFrame(
        {p: _spectrum(data / "per_pigment" / f"{p}.csv") for p in PIGMENTS}
    )
    with np.load(data / "reference" / "surrogate_poly.npz") as npz:
        surrogate_poly = {k: npz[k] for k in npz.files}
    return LegacyInputs(
        enzyme_costs=enzyme_costs,
        enzyme_concentrations_umol=enzyme_concentrations,
        enzyme_parameters=enzyme_parameters,
        min_kcat_values=min_kcats,
        abs_coef=abs_coef,
        light_spectrum=_spectrum(data / "lights" / "warm_white_led.csv"),
        ocp_absorption=_spectrum(data / "OCP_absorption.csv"),
        ps_comp=pd.read_csv(data / "photosystem_composition.csv", index_col=0),
        molar_masses=_spectrum(data / "molar_masses.csv"),
        exp_growth_rate_by_density=pd.read_csv(
            data / "growth_rates_by_density.csv", index_col=0
        )["growth_rate"],
        proteins_of_interest=pd.read_csv(
            data / "Proteomics_Proteins_of_interest.csv", delimiter=";"
        ).drop([0, 15], axis=0),
        surrogate_poly=surrogate_poly,
    )


def proteomics_samples(data: Path = DATA) -> pd.DataFrame:
    """Per-sample OD, gravimetric dry weight and measured protein of the proteomics
    samples, with per-cell values (cells from OD via `CELLS_PER_OD_ML`).

    Columns: sample, od [OD750 in the sample], dry_weight, protein [g L^-1 in the
    sample], cells [cells mL^-1], dry_weight_pg, protein_pg [pg cell^-1],
    protein_fraction [g protein g(DW)^-1, measured].
    """

    def od750(name: str) -> pd.Series:
        df = pd.read_csv(data / name, delimiter=";", header=1, index_col=0)
        return df.loc[750].dropna()

    dil = pd.read_csv(data / "Proteomics_samples_dilutions.csv", delimiter=";")
    dw = pd.read_csv(data / "Proteomics_Dry_weights.csv", delimiter=";")
    prot = pd.read_csv(
        data / "Proteomics_samples_proteinconcentration.csv", delimiter=";"
    )
    od = pd.concat(
        [
            od750("20240328_Spectrum_OD_WT_proteomics.csv"),
            od750("20240329_Spectrum_OD_sqs_proteomics.csv"),
            od750("20240329_Spectrum_OD_dCas9_proteomics.csv"),
        ]
    ).to_numpy(dtype=float)
    out = pd.DataFrame(
        {
            "sample": dil.iloc[:, 0],
            "od": od / dil["dilution factor for OD measurement"].to_numpy(),
            "dry_weight": (dw["Difference"] / dil["Cell dry weight volume [mL]"]),
            "protein": prot["Concentration in cell sample [ug/ml]"] / 1000,
        }
    )
    out["cells"] = out["od"] * CELLS_PER_OD_ML
    out["dry_weight_pg"] = out["dry_weight"] * 1e-3 / out["cells"] * 1e12
    out["protein_pg"] = out["protein"] * 1e-3 / out["cells"] * 1e12
    out["protein_fraction"] = out["protein"] / out["dry_weight"]
    return out
