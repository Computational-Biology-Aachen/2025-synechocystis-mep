"""Figures and summary for docs/reports/productivity-strains.json.

uv run python scripts/plot_productivity_strains.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

RES = Path("docs/reports/productivity-strains.json")
FIG = Path("docs/reports/figures/productivity-strains.png")
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


def short(p: str) -> str:
    return p.removesuffix("_ksynth").removeprefix("E_")


def main() -> None:
    res = json.loads(RES.read_text())
    variants = [v for v in (ORDER_MINIMAL if MINIMAL else ORDER) if v in res]
    fig, axes = plt.subplots(2, len(variants), figsize=(3.3 * len(variants), 6.5),
                             layout="constrained", squeeze=False)  # fmt: skip
    rows = []
    for j, v in enumerate(variants):
        g = pd.DataFrame(res[v]["grid"])
        cmap = plt.get_cmap("viridis", g["round"].max() + 1)
        best_per_round = []
        for k, s in g.groupby("round"):
            axes[0, j].plot(s["density"], s["productivity_mg_L_h"], color=cmap(k),
                            marker="o", ms=2, label=str(k))  # fmt: skip
            ok = s.dropna(subset=["productivity_mg_L_h"])
            if len(ok):
                b = ok.loc[ok["productivity_mg_L_h"].idxmax()]
                best_per_round.append((k, b["productivity_mg_L_h"], b))
        axes[0, j].set(xscale="log", title=v, xlabel="cells mL$^{-1}$")
        ks = [b[0] for b in best_per_round]
        axes[1, j].plot(ks, [b[1] for b in best_per_round], "o-")
        axes[1, j].set(xlabel="no. overexpressed (MCA order)", ylim=(0, None))
        if j == 0:
            axes[0, j].set_ylabel("P [mg L$^{-1}$ h$^{-1}$]")
            axes[1, j].set_ylabel("best P over density")
        _, p_best, row = max(best_per_round, key=lambda b: b[1])
        wt = next(b for b in best_per_round if b[0] == 0)
        rows.append({
            "variant": v,
            "best strain (overexpressed)": ", ".join(short(p) for p in row["overexpressed"])
            or "WT",
            "P max [mg/L/h]": round(p_best, 4),
            "at density [cells/mL]": f"{row['density']:.2e}",
            "mu [1/h]": round(row["mu_per_h"], 4),
            "P WT [mg/L/h]": round(wt[1], 4),
            "gain": round(p_best / wt[1], 1),
            "P after 1 OX": round(next(b[1] for b in best_per_round if b[0] == 1), 4),
        })  # fmt: skip
    axes[0, -1].legend(title="round", fontsize=7)
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG.with_stem(FIG.stem + "-minimal") if MINIMAL else FIG, dpi=100)
    print(pd.DataFrame(rows).to_markdown(index=False))


if __name__ == "__main__":
    main()
