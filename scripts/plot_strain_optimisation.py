"""Figure and summary table for docs/reports/strain-optimisation.json.

uv run python scripts/plot_strain_optimisation.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

RES = Path("docs/reports/strain-optimisation.json")
FIG = Path("docs/reports/figures/strain-optimisation.png")
ORDER = ["v0", "v1", "g-none", "g-coupled", "g-photoacclimation", "g-capacity"]
# --v3: the minimal-model variants next to their thread-g parents, separate figures
ORDER_V3 = [
    "g-coupled",
    "g-capacity",
    "v3-drain",
    "v3-proteome",
    "v3-drain-cap",
    "v3-proteome-cap",
]
V3 = "--v3" in sys.argv


def name(p: str) -> str:
    return p.removesuffix("_ksynth").removeprefix("E_")


def main() -> None:
    res = json.loads(RES.read_text())
    variants = [v for v in (ORDER_V3 if V3 else ORDER) if v in res]
    fig, axes = plt.subplots(
        2, len(variants), figsize=(3.2 * len(variants), 7), layout="constrained",
        height_ratios=(3, 1), squeeze=False,
    )  # fmt: skip
    rows = []
    for j, v in enumerate(variants):
        rounds = [r for r in res[v]["rounds"] if not r.get("failed")]
        scores = pd.DataFrame(
            {
                r["round"]: {
                    name(p): r["response_sq"][p] + r["response_bio"][p]
                    for p in r["response_sq"]
                }
                for r in rounds
                if "response_sq" in r
            }
        )
        ax = axes[0, j]
        sns.heatmap(scores, cmap="RdBu_r", center=0, ax=ax, cbar=j == len(variants) - 1)
        for k, r in enumerate(rounds[1:]):
            gene = name(r["overexpressed"][-1])
            ax.add_patch(
                plt.Rectangle((k, list(scores.index).index(gene)), 1, 1, fill=False,
                              edgecolor="red", linewidth=2)
            )  # fmt: skip
        ax.set(title=v, xlabel="round")
        ax2 = axes[1, j]
        x = [r["round"] for r in rounds]
        ax2.plot(x, [r["sq_content_g_per_L"] for r in rounds], "o-", label="SQ [g/L]")
        ax2.set(xlabel="no. overexpressed", ylim=(0, None))
        ax3 = ax2.twinx()
        ax3.plot(x, [r["growth_rate_per_h"] for r in rounds], "s--", color="C1",
                 label="μ [1/h]")  # fmt: skip
        ax3.set(ylim=(0, None))
        if j == 0:
            ax2.set_ylabel("intracellular SQ [g L^-1]")
        if j == len(variants) - 1:
            ax3.set_ylabel("growth rate [h^-1]")
        first = name(rounds[1]["overexpressed"][-1]) if len(rounds) > 1 else "-"
        rows.append({
            "variant": v,
            "single best (round 1)": first,
            "order": ", ".join(name(p) for p in rounds[-1]["overexpressed"]),
            "SQ start [g/L]": round(rounds[0]["sq_content_g_per_L"], 3),
            "SQ after 1 [g/L]": round(rounds[1]["sq_content_g_per_L"], 3)
            if len(rounds) > 1 else None,
            "SQ end [g/L]": round(rounds[-1]["sq_content_g_per_L"], 3),
            "μ start [1/h]": round(rounds[0]["growth_rate_per_h"], 4),
            "μ end [1/h]": round(rounds[-1]["growth_rate_per_h"], 4),
        })  # fmt: skip
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG.with_stem(FIG.stem + "-v3") if V3 else FIG, dpi=100)
    print(pd.DataFrame(rows).to_markdown(index=False))


if __name__ == "__main__":
    main()
