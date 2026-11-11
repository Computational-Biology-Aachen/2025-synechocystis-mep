"""Pigment contents per cell of the proteomics batch samples (ADR 0006 §2.8 fit targets).

Samples: WT, sqs and dCas9 at 4, 24, 48 and 96 h (triplicates), each with an in vivo
absorption spectrum; the sqs samples also with a methanol-extract spectrum.

- **Methanol extract (sqs):** Chl a = 12.9447 (A665 - A720) µg mL^-1 (Ritchie 2006,
  methanol), total carotenoids = (1000 (A470 - A720) - 2.86 Chl a) / 221 µg mL^-1
  (Wellburn 1994, methanol, no Chl b). Extract → culture: / dilution factor of the
  extract × (methanol volume / extracted culture volume). The methanol volume isn't
  recorded; 1 mL is assumed (`METHANOL_VOLUME`). It gives 6.6–7.8 µg Chl mL^-1 OD^-1,
  consistent with the 7.5 mg OD^-1 L^-1 the models use.
- **In vivo spectra (all strains):** non-negative least squares on 400–700 nm with the
  Fuente2021 in vivo coefficients of Chl a, β-carotene, phycocyanin, allophycocyanin
  plus scattering (constant, λ^-1, λ^-4). Scattering dominates these dilute samples, so
  only Chl a and the PC:Chl ratio are used, and Chl a is calibrated on the sqs extracts
  (median ratio extract / in vivo). Carotenoids come from the extracts only.

Contents in the cell [mM] use cells from OD (`CELLS_PER_OD_ML`) and `CELL_VOLUME` 4 fL;
carotenoids with M = 581.6 g mol^-1 as the models' reference.
"""

from __future__ import annotations

from functools import cache

import numpy as np
import pandas as pd
from scipy.optimize import nnls

from synmep.data import CELLS_PER_OD_ML, DATA, M_CHL

CELL_VOLUME = 4e-15  # [L cell^-1]
M_CAR = 581.5565  # [g mol^-1] as the legacy/reference carotenoid content
METHANOL_VOLUME = 1.0  # [mL] assumption, see module docstring
SPECTRA = {
    "WT": "20240328_Spectrum_OD_WT_proteomics.csv",
    "sqs": "20240329_Spectrum_OD_sqs_proteomics.csv",
    "dCas9": "20240329_Spectrum_OD_dCas9_proteomics.csv",
}
FUENTE = ("chla", "beta_carotene", "phycocyanin", "allophycocyanin")


def _spectra(name: str) -> pd.DataFrame:
    df = pd.read_csv(DATA / name, delimiter=";", header=1, index_col=0)
    return df.dropna(axis=1, how="all")


def _dilutions() -> pd.DataFrame:
    return pd.read_csv(
        DATA / "Proteomics_samples_dilutions.csv", delimiter=";"
    ).set_index("a")


def _per_cell_mM(ug_per_ml: pd.Series, cells_per_ml: pd.Series, molar_mass: float):
    pg = ug_per_ml / cells_per_ml * 1e6  # [pg cell^-1]
    return pg * 1e-12 / molar_mass / CELL_VOLUME * 1e3  # [mM]


def _meta(samples: pd.Index) -> pd.DataFrame:
    parts = samples.str.split()
    return pd.DataFrame(
        {"strain": parts.str[0], "t_h": parts.str[1].astype(int)}, index=samples
    )


@cache
def methanol_extract() -> pd.DataFrame:
    """Chl a and carotenoids of the sqs samples from the methanol extracts."""
    a = _spectra("20240329_Methanolex_sqs_proteomics.csv")
    dil = _dilutions().loc[a.columns]
    od = a.columns.to_series().map(lambda c: _od(c))
    to_culture = (
        METHANOL_VOLUME
        / dil["volume for methanol extraction [mL]"]
        / dil["dilution factor of methanol extract"]
    )
    chl = 12.9447 * (a.loc[665] - a.loc[720]) * to_culture
    car = (
        1000 * (a.loc[470] - a.loc[720]) - 2.86 * 12.9447 * (a.loc[665] - a.loc[720])
    ) / 221
    car = car * to_culture
    cells = od * CELLS_PER_OD_ML
    out = _meta(a.columns)
    out["chl_ug_per_ml_od"] = chl / od
    out["chl_mM"] = _per_cell_mM(chl, cells, M_CHL)
    out["car_mM"] = _per_cell_mM(car, cells, M_CAR)
    return out


def _od(sample: str) -> float:
    strain = sample.split()[0]
    a750 = _spectra(SPECTRA[strain]).loc[750, sample]
    return float(a750) / float(
        _dilutions().loc[sample, "dilution factor for OD measurement"]
    )


@cache
def in_vivo() -> pd.DataFrame:
    """NNLS decomposition of every in vivo spectrum (relative amounts, scattering)."""
    lam = np.arange(400, 701)
    coef = [
        pd.read_csv(DATA / "per_pigment" / f"{p}.csv", comment="#", index_col=0)
        .iloc[:, 0]
        .reindex(lam)
        .to_numpy()
        for p in FUENTE
    ]
    x = np.column_stack(
        [*coef, np.ones(lam.size), (lam / 550) ** -1, (lam / 550) ** -4]
    )
    rows = []
    for name in SPECTRA.values():
        spectra = _spectra(name)
        for sample in spectra.columns:
            y = spectra[sample].reindex(lam).to_numpy(dtype=float)
            amounts, _ = nnls(x, y)
            fit = x @ amounts
            rows.append(
                {"sample": sample}
                | dict(zip(FUENTE, amounts[: len(FUENTE)], strict=True))
                | {"rmse": float(np.sqrt(np.mean((y - fit) ** 2)))}
            )
    out = pd.DataFrame(rows).set_index("sample")
    return _meta(out.index).join(out)


@cache
def pigment_contents() -> pd.DataFrame:
    """Per-sample pigment contents [mM in the cell] for all strains.

    Columns: strain, t_h, chl_mM (in vivo, calibrated on the sqs extracts), car_mM
    (sqs only, extracts), pc_per_chl (in vivo, relative units), chl_extract_mM (sqs).
    """
    iv = in_vivo()
    dil = _dilutions()
    od = pd.Series({s: _od(s) for s in iv.index})
    # in vivo amount per culture volume ~ fitted amount / OD dilution factor
    chl_rel = iv["chla"] / dil.loc[iv.index, "dilution factor for OD measurement"]
    ext = methanol_extract()
    chl_ext_ug = ext["chl_ug_per_ml_od"] * od[ext.index]
    scale = float((chl_ext_ug / chl_rel[ext.index]).median())  # [µg mL^-1 per unit]
    out = iv[["strain", "t_h"]].copy()
    out["chl_mM"] = _per_cell_mM(chl_rel * scale, od * CELLS_PER_OD_ML, M_CHL)
    out["pc_per_chl"] = iv["phycocyanin"] / iv["chla"]
    out["car_mM"] = ext["car_mM"].reindex(out.index)
    out["chl_extract_mM"] = ext["chl_mM"].reindex(out.index)
    return out


def pigment_targets() -> pd.DataFrame:
    """Mean and standard deviation per strain and time point (fit targets)."""
    g = pigment_contents().groupby(["strain", "t_h"])
    cols = ["chl_mM", "car_mM", "pc_per_chl"]
    return g[cols].mean().join(g[cols].std(), rsuffix="_sd")
