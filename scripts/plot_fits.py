"""Batch fits (Fig 9 data) of one family, from the fit JSON files.

    uv run python scripts/plot_fits.py double|minimal

Writes docs/reports/figures/<family>-fit-batch.png and prints the summary table.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPORTS = Path("docs/reports")
FITS = {
    "double": {
        "double-none": "double-fit.json",
        "double-coupled": "double-fit-coupled.json",
        "double-photoacclimation (retired)": "double-fit-photoacclimation.json",
        "double-capacity (retired)": "double-fit-capacity.json",
    },
    "minimal": {
        "double-coupled": "double-fit-coupled.json",
        "minimal-m0": "minimal-m0-fit.json",
        "minimal-m0-cap": "minimal-m0-fit-cap.json",
    },
}


def main() -> None:
    family = sys.argv[1] if len(sys.argv) > 1 else "minimal"
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), layout="constrained")
    rows = []
    for k, (label, name) in enumerate(FITS[family].items()):
        path = REPORTS / name
        if not path.exists():
            continue
        r = json.loads(path.read_text())
        tr = r["trajectory"]
        style = "--" if not label.startswith(family) else "-"
        if k == 0:
            axes[0].plot(tr["t_h"], tr["mu_data"], "k.", label="data")
            axes[1].plot(tr["t_h"], tr["N_data"], "k.", label="data")
        axes[0].plot(tr["t_h"], tr["mu_sim"], style, label=label)
        axes[1].plot(tr["t_h"], tr["N_sim"], style, label=label)
        axes[2].plot(tr["t_h"], tr["Chl"], style, label=label)
        rows.append({
            "model": label,
            "fitted": ", ".join(f"{p} {v:.3g}" for p, v in r["fitted"].items()),
            "RMSE μ [1/h]": round(r["rmse_mu_per_h"], 4),
            "RMSE ln N": round(r["rmse_ln_N"], 3),
        })  # fmt: skip
    axes[0].set(xlabel="time [h]", ylabel="growth rate [h$^{-1}$]")
    axes[1].set(xlabel="time [h]", ylabel="cells mL$^{-1}$", yscale="log")
    axes[2].set(xlabel="time [h]", ylabel="Chl [mM]", yscale="log")
    axes[0].legend(fontsize=8)
    fig.savefig(REPORTS / "figures" / f"{family}-fit-batch.png", dpi=100)
    print(pd.DataFrame(rows).to_markdown(index=False))


if __name__ == "__main__":
    main()
