"""Comparison of all variants: gates, ranking, goal answers (ADR 0006 §2.9, stage 7).

    uv run python scripts/summary.py

Reads docs/reports/fits/*.json and docs/reports/goals/*.json, writes
docs/reports/summary.json and figures/summary-*.png, prints the markdown tables used
in report 10.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPORTS = Path("docs/reports")
FIG = REPORTS / "figures"
# Gate 1 on the PETC level (report 08): median |ratio error| of the paired spectral
# ranking vs pfennig2024; `full` contains pfennig2024 itself
SPECTRAL_ERROR = {"single-legacy": np.nan, "single-fixed": 0.101, "double": 0.101,
                  "minimal-m0": 0.030, "minimal-m2": 0.014, "full": 0.0}  # fmt: skip
GATE_MU, GATE_LN_N = 0.0058, 0.149  # single-legacy + 10 % (ADR 0006)


def spectral_error(name: str) -> float:
    for key, v in SPECTRAL_ERROR.items():
        if name.startswith(key):
            return v
    return np.nan


def short(names: list[str]) -> str:
    return ", ".join(e.removeprefix("E_").removesuffix("_ksynth") for e in names)


def load(folder: str) -> dict[str, dict]:
    return {
        p.stem: json.loads(p.read_text())
        for p in sorted((REPORTS / folder).glob("*.json"))
    }


def main() -> None:
    fits, goals = load("fits"), load("goals")
    rows = []
    for name in sorted(set(fits) | set(goals)):
        f, g = fits.get(name, {}), goals.get(name, {})
        conv = g.get("convergence", {})
        ok = sum(c["ok"] for c in conv.values())
        n = sum(c["n"] for c in conv.values())
        strains = g.get("strains", [])
        batch = [
            s["batch"]["best_productivity_mg_L_h"] for s in strains if s.get("batch")
        ]
        titer = [s["batch"]["titer_end_mg_L"] for s in strains if s.get("batch")]
        turb = [s["turbidostat_best"]["productivity_mg_L_h"] for s in strains
                if s.get("turbidostat_best")]  # fmt: skip
        spectra = {
            k: v["best_productivity_mg_L_h"]
            for k, v in g.get("spectra", {}).items()
            if v
        }
        best_spec = max(spectra, key=spectra.get) if spectra else None
        row = {
            "variant": name,
            "k": f.get("n_parameters"),
            "rmse_mu": f.get("rmse_mu_per_h"),
            "rmse_lnN": f.get("rmse_ln_N"),
            "aic": f.get("aic"),
            "validation_lnN": f.get("validation_rmse_ln_N"),
            "fig10_mu": f.get("fig10_rmse_mu_per_h"),
            "spectral_err": spectral_error(name),
            "convergence": ok / n if n else np.nan,
            "sequence": short(g.get("mca", {}).get("sequence", [])),
            "best_strain": short(g.get("best_strain", [])),
            "batch_P": max(batch) if batch else np.nan,
            "titer96": max(titer) if titer else np.nan,
            "turbidostat_P": max(turb) if turb else np.nan,
            "best_light": best_spec,
            "best_light_P": spectra.get(best_spec, np.nan) if best_spec else np.nan,
            "warm_white_P": spectra.get("warm_white_led:1", np.nan),
        }
        row["gate1"] = bool(row["spectral_err"] <= 0.05)
        row["gate2"] = bool(
            (row["rmse_mu"] or 1) <= GATE_MU and (row["rmse_lnN"] or 9) <= GATE_LN_N
        )
        row["gate3"] = bool(row["convergence"] >= 0.95)
        row["passes"] = row["gate1"] and row["gate2"] and row["gate3"]
        rows.append(row)
    df = pd.DataFrame(rows)
    passing = df[df["passes"]]
    ranking = (
        passing.sort_values(["validation_lnN", "aic"]) if len(passing) else passing
    )
    out = {
        "table": df.to_dict("records"),
        "ranking": list(ranking["variant"]),
        "spread_passing": {
            k: [float(passing[k].min()), float(passing[k].max())]
            for k in ("batch_P", "titer96", "turbidostat_P", "best_light_P")
            if len(passing)
        },
    }
    (REPORTS / "summary.json").write_text(json.dumps(out, indent=1, default=str))
    cols1 = ["variant", "k", "rmse_mu", "rmse_lnN", "aic", "validation_lnN", "fig10_mu",
             "spectral_err", "convergence", "gate1", "gate2", "gate3"]  # fmt: skip
    cols2 = ["variant", "sequence", "best_strain", "batch_P", "titer96", "turbidostat_P",
             "best_light", "best_light_P", "warm_white_P"]  # fmt: skip
    print(df[cols1].to_markdown(index=False, floatfmt=".4g"))
    print()
    print(df[cols2].to_markdown(index=False, floatfmt=".3g"))
    print("\nranking:", out["ranking"], "\nspread:", out["spread_passing"])

    # figure: batch fits
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), layout="constrained")
    data = None
    for name, f in fits.items():
        tr = f.get("trajectory")
        if not tr:
            continue
        axes[0].plot(tr["t_h"], tr["mu_per_h"], label=name)
        axes[1].plot(tr["t_h"], tr["N"], label=name)
        axes[2].plot(tr["t_h"], tr["chl_mM"], label=name)
        data = data or tr
    from synmep.evaluate import _growth_data

    d = _growth_data()
    t = np.linspace(0, 96, len(d))
    axes[0].plot(t, d.to_numpy(), "k.", label="data")
    axes[1].plot(t, d.index, "k.")
    axes[0].set(xlabel="t [h]", ylabel="μ [h$^{-1}$]")
    axes[1].set(xlabel="t [h]", ylabel="cells mL$^{-1}$", yscale="log")
    axes[2].set(xlabel="t [h]", ylabel="Chl [mM]", yscale="log")
    axes[0].legend(fontsize=7)
    fig.savefig(FIG / "summary-fits.png", dpi=100)

    # figure: goal 2 spectra per variant (batch productivity relative to warm white)
    fig, ax = plt.subplots(figsize=(13, 4.5), layout="constrained")
    for name, g in goals.items():
        sp = {
            k: v["best_productivity_mg_L_h"]
            for k, v in g.get("spectra", {}).items()
            if v
        }
        ref = sp.get("warm_white_led:1")
        if not ref:
            continue
        keys = list(sp)
        ax.plot(range(len(keys)), [sp[k] / ref for k in keys], "o-", label=name)
    if goals:
        keys = list(next(iter(goals.values())).get("spectra", {}))
        ax.set_xticks(range(len(keys)), keys, rotation=70, fontsize=7)
    ax.axhline(1, color="grey", lw=0.5)
    ax.set(ylabel="batch productivity / warm white")
    ax.legend(fontsize=7)
    fig.savefig(FIG / "summary-spectra.png", dpi=100)


if __name__ == "__main__":
    main()
