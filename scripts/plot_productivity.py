"""Figures and summary for docs/reports/productivity.json.

uv run python scripts/plot_productivity.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

RES = Path("docs/reports/productivity.json")
FIG = Path("docs/reports/figures")
ORDER = [
    "single-legacy",
    "single-fixed",
    "double-none",
    "double-coupled",
    "double-photoacclimation",
    "double-capacity",
]
# --minimal: the minimal variants next to their `double` parents, separate figures
ORDER_MINIMAL = [
    "double-coupled",
    "double-capacity",
    "minimal-m0-drain",
    "minimal-m0",
    "minimal-m0-drain-cap",
    "minimal-m0-cap",
]
MINIMAL = "--minimal" in sys.argv
PIGMENT_TRADE_OFF = {"single-legacy": True, "single-fixed": True, "double-none": True, "double-coupled": True,
                     "double-photoacclimation": False, "double-capacity": False,
                     "minimal-m0-drain": True, "minimal-m0": True, "minimal-m0-drain-cap": True,
                     "minimal-m0-cap": True}  # fmt: skip
SUFFIX = "-minimal" if MINIMAL else ""


def main() -> None:
    res = json.loads(RES.read_text())
    variants = [v for v in (ORDER_MINIMAL if MINIMAL else ORDER) if v in res]
    grids = {v: pd.DataFrame(res[v]["grid"]) for v in variants}

    fig, axes = plt.subplots(1, len(variants), figsize=(3.3 * len(variants), 3.2),
                             layout="constrained", squeeze=False)  # fmt: skip
    rows = []
    for ax, v in zip(axes[0], variants, strict=True):
        g = grids[v]
        p = g.pivot(index="sqs_factor", columns="density", values="productivity_mg_L_h")
        mesh = ax.pcolormesh(p.columns, p.index, p.to_numpy(), shading="nearest",
                             cmap="viridis")  # fmt: skip
        fig.colorbar(mesh, ax=ax, label="P [mg L$^{-1}$ h$^{-1}$]")
        best = g.loc[g["productivity_mg_L_h"].idxmax()]
        ax.plot(best["density"], best["sqs_factor"], "r*", ms=12)
        ax.set(xscale="log", yscale="log", title=v, xlabel="cells mL$^{-1}$")
        if v == variants[0]:
            ax.set_ylabel("squalene synthase [x WT]")
        wt = g[np.isclose(g["sqs_factor"], 1.0)]
        wt_best = wt.loc[wt["productivity_mg_L_h"].idxmax()]
        edge = np.isclose(best["sqs_factor"], g["sqs_factor"].max())
        rows.append(
            {
                "variant": v,
                "pigment trade-off": "yes" if PIGMENT_TRADE_OFF[v] else "no (fixed)",
                "P max [mg/L/h]": round(best["productivity_mg_L_h"], 4),
                "at Sqs [x WT]": round(best["sqs_factor"], 2)
                if not edge else f">= {best['sqs_factor']:.0f} (scan edge)",
                "at density [cells/mL]": f"{best['density']:.2e}",
                "mu at optimum [1/h]": round(best["mu_per_h"], 4),
                "P at WT Sqs [mg/L/h]": round(wt_best["productivity_mg_L_h"], 4),
                "gain vs WT Sqs": round(
                    best["productivity_mg_L_h"] / wt_best["productivity_mg_L_h"], 1
                ),
                "stoich. bound [mg/L/h]": round(best["bound_mg_L_h"], 3)
                if "bound_mg_L_h" in g and pd.notna(best.get("bound_mg_L_h")) else None,
            }
        )  # fmt: skip
    fig.savefig(FIG / f"productivity-maps{SUFFIX}.png", dpi=100)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), layout="constrained")
    for v in variants:
        g = grids[v]
        n_opt = g.loc[g["productivity_mg_L_h"].idxmax(), "density"]
        s = g[g["density"] == n_opt].sort_values("sqs_factor")
        style = "-" if PIGMENT_TRADE_OFF[v] else ":"
        axes[0].plot(s["sqs_factor"], s["mu_per_h"], style, marker="o", ms=3, label=v)
        axes[1].plot(s["sqs_factor"], s["chl"], style, marker="o", ms=3, label=v)
        axes[2].plot(s["sqs_factor"], s["productivity_mg_L_h"], style, marker="o",
                     ms=3, label=v)  # fmt: skip
        if "bound_mg_L_h" in s:
            axes[2].plot(s["sqs_factor"], s["bound_mg_L_h"], style, color="grey",
                         lw=0.8)  # fmt: skip
    axes[0].set(ylabel="growth rate [h$^{-1}$]")
    axes[1].set(ylabel="Chl [mM]", yscale="log")
    axes[2].set(ylabel="P [mg L$^{-1}$ h$^{-1}$] (grey: stoich. bound)", yscale="log")
    for ax in axes:
        ax.set(xscale="log", xlabel="squalene synthase [x WT]")
    axes[0].legend(fontsize=8)
    fig.suptitle("trade-off at each variant's optimal density (dotted: no pigment "
                 "trade-off by construction)")  # fmt: skip
    fig.savefig(FIG / f"productivity-trade-off{SUFFIX}.png", dpi=100)
    print(pd.DataFrame(rows).to_markdown(index=False))


if __name__ == "__main__":
    main()
