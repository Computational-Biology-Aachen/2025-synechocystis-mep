"""Legacy (modelbase) setup, extracted verbatim from main.ipynb cell 1.

Only runs in the old pixi environment (modelbase + assimulo + scikit-learn==1.5.0).
Kept for provenance: the reference outputs in data/reference/ were generated with it.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import os
os.chdir(ROOT)

import numpy as np
import matplotlib.pyplot as plt
from modelbase.ode import Model, Simulator, mca
import pandas as pd
import seaborn as sns
import pickle
from model import (
    get_enzyme_costs,
    get_model,
    get_enzyme_concentrations,
    get_enzyme_parameters,
    get_proteomics_results,
    calculate_min_kcats,
)
from modelbase.ode.simulators import _Simulate
from parallel import parallelise
from pathlib import Path
from sklearn.pipeline import Pipeline
from typing import Literal, cast
from functools import partial
from typing import TypeVar


T = TypeVar("T")


def unwrap(x: T | None) -> T:
    if x is None:
        raise ValueError("Unexpected None")
    return x


CELL_VOLUME: float = 2e-15  # FIXME: explain values
PROTEIN_CONTENT: float = 0.445  # FIXME: explain values
MODE = "blast"
CHLA = (
    7.5 / 893.509 / (105482164.18012363 * 1000)  # FIXME: explain values
)  # mg/OD/L -> mmol/OD/L / Cells/LLOD --> mmol/cell
CAROT = 2 / 581.5565 / (105482164.18012363 * 1000)  # FIXME: explain values
SQUALENE = 0.5 / 410.73 / (105482164.18012363 * 1000)  # FIXME: explain values
DOUBLING_TIME: float = 300.0  # FIXME: explain values
DAY_IN_MIN = 60 * 24

# Load data
DATA = Path("data")

with open(DATA / "rate_regression_model_biomass.pickle", "rb") as f:
    BIOMASS_REGRESSION_MODEL = pickle.load(f)

OCP_ABSORPTION = pd.read_csv(DATA / "OCP_absorption.csv", index_col=0)
PS_COMP = pd.read_csv(DATA / "photosystem_composition.csv", index_col=0)
MOLAR_MASSES = cast(
    pd.Series, pd.read_csv(DATA / "molar_masses.csv", index_col=0).squeeze(1)
)

PIGMENTS_FUENTE = cast(
    pd.Series,
    pd.read_csv(
        DATA / "pigment_concentrations.csv",
        index_col=0,
    ).squeeze(1),
)

ABS_COEF = pd.concat(
    {
        pig: pd.read_csv(DATA / f"per_pigment/{pig}.csv", index_col=0)
        for pig in PIGMENTS_FUENTE.index
    },
    axis=1,
)

ENZYME_COSTS = get_enzyme_costs(
    kugler_model_csv=DATA / "Amino_acid_costs_Kuglermodel.csv",
    aa_costs_csv=DATA / "Amino_acid_costs_arnold.csv",
)

ABUNDANCES, DETECTION_LIMIT = get_proteomics_results(
    proteomics_diff_csv=DATA / "Proteomics_DiffAnalysis_2WayAnova_.csv",
    proteomics_poi_csv=DATA / "Proteomics_Proteins_of_interest.csv",
    proteomics_concentration_csv=DATA / "Proteomics_samples_proteinconcentration.csv",
    proteomics_dry_weight_csv=DATA / "Proteomics_Dry_weights.csv",
    proteomics_dilutions_csv=DATA / "Proteomics_samples_dilutions.csv",
    proteomics_od_wt_csv=DATA / "20240328_Spectrum_OD_WT_proteomics.csv",
    proteomics_od_sqs_csv=DATA / "20240329_Spectrum_OD_sqs_proteomics.csv",
    proteomics_od_dcas9_csv=DATA / "20240329_Spectrum_OD_dCas9_proteomics.csv",
    cell_volume=CELL_VOLUME,
    protein_content=PROTEIN_CONTENT,
)

ENZYME_CONCENTRATIONS_UMOL: dict[str, float] = get_enzyme_concentrations(
    equilibrium_constants_csv=DATA / "Equilibrium_constants.csv",
    abundances=ABUNDANCES,
    detection_limit=DETECTION_LIMIT,
)

ENZYME_PARAMETERS: dict[str, float] = get_enzyme_parameters(
    mode=MODE,
    blasted_parameters_csv=DATA / "enzyme_parameters_blasted.csv",
    predicted_kms_csv=DATA / "predicted_Kms.csv",
    equilibrium_constants_csv=DATA / "Equilibrium_constants.csv",
)

MIN_KCAT_VALUES: dict[str, float] = calculate_min_kcats(
    enzyme_concentrations_umol=ENZYME_CONCENTRATIONS_UMOL,
    cell_volume=CELL_VOLUME,
    doubling_time=DOUBLING_TIME,
)


LIGHTS = pd.DataFrame(
    {
        path.stem: pd.read_csv(path, index_col=0).squeeze(1)
        for path in sorted((DATA / "lights").glob("*.csv"))
    }
)
EXP_GROWTH_RATE_BY_DENSITY = pd.read_csv(
    DATA / "growth_rates_by_density.csv",
    index_col=0,
)["growth_rate"]

PROTEINS_OF_INTEREST = pd.read_csv(
    DATA / "Proteomics_Proteins_of_interest.csv", delimiter=";"
).drop([0, 15], axis=0)


def get_parameterised_model(
    strain: Literal["wt", "dshc"],
    growth_mode: Literal["steady", "dynamic"],
    kcat_values: dict[str, float] = MIN_KCAT_VALUES,
    cell_volume: float = CELL_VOLUME,
    enzyme_costs: pd.DataFrame = ENZYME_COSTS,
    enzyme_concentrations_umol: dict[str, float] = ENZYME_CONCENTRATIONS_UMOL,
    enzyme_parameters: dict[str, float] = ENZYME_PARAMETERS,
    biomass_regression_model: Pipeline = BIOMASS_REGRESSION_MODEL,
    abs_coef: pd.DataFrame = ABS_COEF,
    lights: pd.DataFrame = LIGHTS,
    ocp_absorption: pd.DataFrame = OCP_ABSORPTION,
    molar_masses: pd.Series = MOLAR_MASSES,
    ps_comp: pd.DataFrame = PS_COMP,
) -> tuple[Model, dict[str, float]]:
    m = get_model(
        growth_mode=growth_mode,
        strain=strain,
        cell_volume=cell_volume,
        enzyme_total_cost_mol=enzyme_costs,
        enzyme_concentrations_umol=enzyme_concentrations_umol,
        enzyme_parameters=enzyme_parameters,
        biomass_regression_model=biomass_regression_model,
        abs_coef=abs_coef,
        lights=lights,
        ocp_absorption=ocp_absorption,
        molar_masses=molar_masses,
        ps_comp=ps_comp,
    )
    m.update_parameters(
        {
            # "ShcKcat": 0,  # default 600
            # "Shc_ksynth": 0, default 0.008034116888921783
            "light": 100,
            "ATP_out_metabolism": 0.5,
            "3PGA_out_metabolism": 0.5,
            "NADPH_out_metabolism": 0.5,
            "Fdred_out_metabolism": 0.5,
            "cell_volume": cell_volume,
            "DXS_inh_str": 0.001,
            "ThiGDXPKm": 0.01,
        }
    )
    m.update_parameters(
        {k: v for k, v in kcat_values.items() if k not in {"IspHKcat", "CrtEKcat"}}
    )

    # Update initial
    starting_conc = {}
    for i in m.compounds:
        starting_conc[i] = 0.01
    starting_conc.update(enzyme_concentrations_umol)
    starting_conc.update(
        {
            "3PGA": 0.1,
            "ATP": 0.1,
            "NADPH": 0.1,
            "Fdred": 0.1,
            "ThiDP": 10,
            "chlorophyll-phytol": CHLA / cell_volume,
            "Cell_density": 1e1,
            "carotenoids": CAROT / cell_volume,
            "SQ": SQUALENE / cell_volume,
        }
    )
    return m, starting_conc


def get_doubling_time(rel_growth_rate: float) -> float:
    return np.log(2) / (rel_growth_rate * 60)


def generate_initial_conditions() -> dict[str, float]:
    ms, y0 = get_parameterised_model(
        strain="dshc",
        growth_mode="steady",
    )
    return unwrap(
        Simulator(ms)
        # FIXME: explain why this number
        .initialise(y0 | {"Cell_density": EXP_GROWTH_RATE_BY_DENSITY.index[24] / 1e6})
        .simulate_to_steady_state_and()
        .get_new_y0()
    )

