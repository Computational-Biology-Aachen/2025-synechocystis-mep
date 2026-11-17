"""Figures of the executive report for one variant (default: the recommended one).

    uv run python scripts/executive_figures.py [variant]

Reads docs/reports/fits/<variant>.json and goals/<variant>.json, runs the Sqs scan and
the steady states for the trade-off panels, writes
docs/reports/figures/exec-{fit,story,goals}.png (goals: goal 1 only) and exec-<variant>.json (numbers).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from shared_fit import TIMES, targets

from synmep import evaluate, variants

R = Path("docs/reports")
FIG = R / "figures"
NAME = sys.argv[1] if len(sys.argv) > 1 else "minimal-m2-cap"
STEADY_DENSITY = 3e8  # [cells mL^-1] mid-batch density for the partitioning
SQS_FOLDS = (1, 2, 4, 7, 12, 20, 35, 60)
# carbon atoms per reaction event (3PGA = C3; IPP = C5 from 6 C, 1 C lost as CO2)
CARBON = {"biomass": ("v_bio", 6.934 * 3), "storage": ("pga_storage", 3),
          "squalene": ("v_sq", 30), "carotenoids": ("v_car", 40),
          "Chl (phytol)": ("v_chl", 20)}  # fmt: skip


def label(strain: list[str]) -> str:
    return " + ".join(s.removeprefix("E_") for s in strain) or "WT"


def fit_figure(fit: dict) -> None:
    tr = fit["trajectory"]
    data = evaluate._growth_data()
    t_data = np.linspace(0, 96, len(data))
    tg = targets()
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4), layout="constrained")
    ax[0].plot(t_data, data.to_numpy(), "k.", ms=3, label="data (Fig 9)")
    ax[0].plot(tr["t_h"], tr["mu_per_h"], "C0", label=NAME)
    ax[0].set(xlabel="time [h]", ylabel="growth rate [h$^{-1}$]", title="a  batch growth")
    ax[1].plot(t_data, data.index.to_numpy(), "k.", ms=3)
    ax[1].plot(tr["t_h"], tr["N"], "C0")
    ax[1].set(xlabel="time [h]", ylabel="cells mL$^{-1}$", yscale="log",
              title="b  cell density")  # fmt: skip
    ax[2].plot(TIMES, tg["chl"], "ks", label="Chl data")
    ax[2].plot(tr["t_h"], tr["chl_mM"], "C2", label="Chl model")
    ax[2].set(xlabel="time [h]", ylabel="Chl [mM in cell]", title="c  pigments")
    a2 = ax[2].twinx()
    a2.plot(TIMES, tg["pc"], "k^", mfc="none", label="PC:Chl data")
    a2.plot(tr["t_h"], tr["pc_per_chl"], "C1--", label="PC:Chl model")
    a2.set_ylabel("PC:Chl [mg mg$^{-1}$]")
    h = [*ax[2].get_legend_handles_labels()[0], *a2.get_legend_handles_labels()[0]]
    ax[2].legend(handles=h, fontsize=7, loc="upper center")
    ax[0].legend(fontsize=7)
    fig.savefig(FIG / "exec-fit.png", dpi=150)


def steady(a, folds: dict, density: float):
    m = evaluate.steady_model(a, folds, density, **evaluate.BATCH_LIGHT)
    return evaluate.steady_state(a, m)


def story_figure(a, goals: dict) -> dict:
    folds_all = goals["mca"]["folds"]
    out: dict = {"sqs_scan": {}, "partition": {}, "chl_vs_mu": {}}
    # a: Sqs level in the WT and the CrtE background, batch
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.6), layout="constrained")
    for bg, style in (({}, "-o"), ({"E_CrtE": folds_all["E_CrtE"]}, "--s")):
        rows = []
        for f in SQS_FOLDS:
            b = evaluate.batch(a, bg | {"E_Sqs": float(f)})
            if b is None:
                rows.append((f, np.nan, np.nan, np.nan))
                continue
            i48 = int(np.argmin(np.abs(np.array(b["t_h"]) - 48)))
            rows.append((f, b["best_productivity_mg_L_h"], b["chl_mM"][i48],
                         np.log(b["N"][-1] / b["N"][0]) / 96))  # fmt: skip
        r = np.array(rows, dtype=float)
        key = "CrtE background" if bg else "WT background"
        out["sqs_scan"][key] = r.tolist()
        ax[0].plot(r[:, 0], r[:, 1], style, color="C0", label=f"productivity, {key}")
        ax[0].plot(r[:, 0], r[:, 3] / r[0, 3], style, color="C3", ms=4,
                   label=f"mean μ (rel.), {key}")  # fmt: skip
        ax[0].plot(r[:, 0], r[:, 2] / r[0, 2], style, color="C2", ms=4,
                   label=f"Chl at 48 h (rel.), {key}")  # fmt: skip
    ax[0].axvline(folds_all["E_Sqs"], color="0.5", lw=0.8, ls=":")
    ax[0].set(xscale="log", xlabel="Sqs level [× WT]",
              ylabel="mg L$^{-1}$ h$^{-1}$  /  relative",
              title="a  squalene vs growth and pigments (batch)")  # fmt: skip
    ax[0].legend(fontsize=6)
    # b: carbon partitioning at a mid-batch steady state, MCA strains
    strains_ = [[]] + [goals["mca"]["sequence"][: k + 1]
                       for k in range(len(goals["mca"]["sequence"]))]  # fmt: skip
    bottom = np.zeros(len(strains_))
    shares = {k: [] for k in CARBON}
    for s in strains_:
        row = steady(a, {e: folds_all[e] for e in s}, STEADY_DENSITY)
        c_in = float(row["supply_PGA"]) * 3
        for k, (flux, c) in CARBON.items():
            shares[k].append(100 * float(row[flux]) * c / c_in)
        out["partition"][label(s)] = {k: v[-1] for k, v in shares.items()} | {
            "mu_per_h": float(row["mu"]) * 60
        }
    for k, v in shares.items():
        ax[1].bar(range(len(strains_)), v, bottom=bottom, label=k)
        bottom += np.array(v)
    ax[1].set_xticks(range(len(strains_)), [label(s) for s in strains_], fontsize=7,
                     rotation=20)  # fmt: skip
    ax[1].set(ylabel="% of fixed carbon",
              title=f"b  carbon partitioning ({STEADY_DENSITY:.0e} cells mL$^{{-1}}$)")
    ax[1].legend(fontsize=7)
    # c: Chl synthesis vs growth rate across densities and strains (paper's correlation)
    for j, s in enumerate(strains_):
        pts = []
        for n in np.geomspace(1e7, 3e9, 8):
            row = steady(a, {e: folds_all[e] for e in s}, float(n))
            if row is not None:
                # Chl synthesis per cell volume [µM h^-1]
                pts.append((float(row["mu"]) * 60, float(row["v_chl"]) * 60 * 1e3))
        p = np.array(pts)
        out["chl_vs_mu"][label(s)] = p.tolist()
        ax[2].plot(p[:, 0], p[:, 1], "o-", color=f"C{j}", ms=3, label=label(s))
    ax[2].set(xlabel="growth rate [h$^{-1}$]", ylabel="Chl synthesis [µM h$^{-1}$]",
              title="c  Chl synthesis vs growth (steady states)")  # fmt: skip
    ax[2].legend(fontsize=7)
    fig.savefig(FIG / "exec-story.png", dpi=150)
    return out


def goals_figure(goals: dict) -> None:
    """Goal 1 only; the best-light goal (goal 2) is follow-up work."""
    fig, ax = plt.subplots(figsize=(5.5, 3.8), layout="constrained")
    st = goals["strains"]
    vals = [s["batch"]["best_productivity_mg_L_h"] if s["batch"] else 0 for s in st]
    ax.bar(range(len(st)), vals, color="C0")
    ax.set_xticks(range(len(st)), [label(s["strain"]) for s in st], rotation=25,
                  fontsize=7)  # fmt: skip
    ax.set(ylabel="batch productivity [mg L$^{-1}$ h$^{-1}$]",
           title="MCA strain sequence (warm white, 96 h)")
    fig.savefig(FIG / "exec-goals.png", dpi=150)


def main() -> None:
    fit = json.loads((R / "fits" / f"{NAME}.json").read_text())
    goals = json.loads((R / "goals" / f"{NAME}.json").read_text())
    a = variants.adapter(NAME)
    fit_figure(fit)
    goals_figure(goals)
    out = story_figure(a, goals)
    (R / f"exec-{NAME}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out["partition"], indent=1))
    print(json.dumps(out["sqs_scan"], indent=1))


if __name__ == "__main__":
    main()
