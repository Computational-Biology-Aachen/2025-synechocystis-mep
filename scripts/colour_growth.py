"""Relative WT growth under narrow-band LEDs vs Zavřel et al. 2024 (report 11).

    uv run --extra surrogate python scripts/colour_growth.py [variant ...]

Data: `data/zavrel2024_light_quality.csv` (nine LEDs, 435-687 nm, 25 µmol m^-2 s^-1,
low-density flasks). Model: WT steady state at 2.1e7 cells mL^-1 (OD750 0.2), 1 cm,
25 µmol m^-2 s^-1, one LED at a time (Gaussian, 20 nm FWHM), growth relative to the
633 nm LED. This tests what the shared fit can't: how the PETC and the acclimation
rules respond to light colour. Limits: continuous light instead of 14:10 h, steady
state instead of exponential flask growth. Writes docs/reports/colour-growth.json.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from synmep import evaluate, variants
from synmep.data import DATA

warnings.filterwarnings("ignore")
OUT = Path("docs/reports/colour-growth.json")
REF_NM = 633
DENSITY = 0.2 * 1.0548e8  # [cells mL^-1] OD750 0.2 (cells per OD, assumptions.md)
LIGHT = {"light": 25.0, "sample_depth": 0.01}
DEFAULT = ("minimal-m0", "minimal-m0-cap", "minimal-m0-cap-abs", "minimal-m0-cap-flux",
           "minimal-m2", "minimal-m2-cap")  # fmt: skip


def data() -> pd.Series:
    d = pd.read_csv(DATA / "zavrel2024_light_quality.csv", comment="#", index_col=0)
    return d["growth_rate"] / d.loc[REF_NM, "growth_rate"]


def model_growth(name: str, peaks: list[int]) -> dict[int, float | None]:
    a = variants.adapter(name)
    chans = tuple(f"led_{p}" for p in peaks)
    out = {}
    for p in peaks:
        w = {f"w_{c}": float(c == f"led_{p}") for c in chans}
        m = evaluate.steady_model(a, {}, DENSITY, light_channels=chans, pars=w, **LIGHT)
        row = evaluate.steady_state(a, m, growing=True)
        out[p] = None if row is None else float(row[a.keys["mu"]]) * 60
    return out


def main() -> None:
    names = sys.argv[1:] or [v + s for v in DEFAULT for s in ("", "+ca")]
    rel_data = data()
    peaks = [int(p) for p in rel_data.index]
    report = json.loads(OUT.read_text()) if OUT.exists() else {}
    report["data_relative"] = rel_data.to_dict()
    for name in names:
        if name.endswith("+ca"):
            ca = variants.adapter(name)
            ref = evaluate.calibrate_chromatic(variants.adapter(name[:-3]), ca)
            ca.fixed["chromatic_reference"] = ref
        mu = model_growth(name, peaks)
        ref = mu.get(REF_NM)
        rel = {p: (v / ref if v and ref else None) for p, v in mu.items()}
        ok = [p for p in peaks if rel[p] is not None]
        err = np.array([rel[p] - rel_data[p] for p in ok])
        report[name] = {
            "mu_per_h": mu,
            "relative": rel,
            "rmse_relative": float(np.sqrt(np.mean(err**2))) if ok else None,
            "spearman": float(spearmanr([rel[p] for p in ok],
                                        [rel_data[p] for p in ok])[0]) if ok else None,
            "n_failed": len(peaks) - len(ok),
        }  # fmt: skip
        r = report[name]
        print(f"{name}: rmse {r['rmse_relative']}, rho {r['spearman']}, "
              f"failed {r['n_failed']}", flush=True)  # fmt: skip
        print("   ", {p: (round(v, 2) if v is not None else None)
                      for p, v in rel.items()}, flush=True)  # fmt: skip
        OUT.write_text(json.dumps(report, indent=1, default=str))
    print("data", {p: round(v, 2) for p, v in rel_data.items()})


if __name__ == "__main__":
    main()
