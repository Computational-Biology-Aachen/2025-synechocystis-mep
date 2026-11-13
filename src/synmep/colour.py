"""Coloured-light growth data (Zavřel et al. 2024) and the model's response to it.

WT under nine narrow-band LEDs (435-687 nm) at 25 µmol m^-2 s^-1, low-density flasks
(`data/zavrel2024_light_quality.csv`). The model is the WT steady state at OD750 0.2
(2.1e7 cells mL^-1), 1 cm, one LED at a time; growth and Chl are compared relative to
the 633 nm LED (the data's growth unit is unclear, report 11). A culture without a
growing state counts as growth 0: that is the model's prediction, not a failure.
Used as a target of the shared fit (ADR 0006 §8) and by `scripts/colour_growth.py`.
"""

from __future__ import annotations

from functools import cache

import pandas as pd

from synmep import evaluate
from synmep.data import DATA

REF_NM = 633
DENSITY = 0.2 * 1.0548e8  # [cells mL^-1] OD750 0.2 (cells per OD, assumptions.md)
LIGHT = {"light": 25.0, "sample_depth": 0.01}


@cache
def data() -> pd.DataFrame:
    """Growth and Chl per cell relative to the 633 nm LED, by peak wavelength [nm]."""
    d = pd.read_csv(DATA / "zavrel2024_light_quality.csv", comment="#", index_col=0)
    return pd.DataFrame({
        "growth": d["growth_rate"] / d.loc[REF_NM, "growth_rate"],
        "chl": d["chl_fmol_cell"] / d.loc[REF_NM, "chl_fmol_cell"],
    })  # fmt: skip


def response(a: evaluate.Adapter, pars: dict[str, float] | None = None) -> pd.DataFrame:
    """Model mu [h^-1], Chl [mM] and both relative to 633 nm, per LED."""
    peaks = [int(p) for p in data().index]
    chans = tuple(f"led_{p}" for p in peaks)
    rows, states = {}, {}
    # 633 nm first (with fallbacks), then outwards, each from its neighbour's state
    for p in sorted(peaks, key=lambda p: abs(p - REF_NM)):
        w = {f"w_{c}": float(c == f"led_{p}") for c in chans}
        m = evaluate.steady_model(
            a, {}, DENSITY, light_channels=chans, pars=dict(pars or {}) | w, **LIGHT
        )
        done = [q for q in states if (q - REF_NM) * (p - REF_NM) >= 0] or list(states)
        if done:
            y0 = states[min(done, key=lambda q: abs(q - p))]
            m.update_variables({k: v for k, v in y0.items() if k != a.n_name})
        row = evaluate.steady_state(a, m, growing=True, fallback=p == REF_NM)
        if row is None:
            rows[p] = (0.0, float("nan"))
            continue
        rows[p] = (float(row[a.keys["mu"]]) * 60, float(row[a.keys["chl"]]))
        states[p] = {k: float(row[k]) for k in m.get_variable_names()}
    out = pd.DataFrame(rows, index=["mu_per_h", "chl_mM"]).T.sort_index()
    ref = out.loc[REF_NM]
    out["growth"] = out["mu_per_h"] / ref["mu_per_h"] if ref["mu_per_h"] > 0 else 0.0
    out["chl"] = out["chl_mM"] / ref["chl_mM"]
    return out
