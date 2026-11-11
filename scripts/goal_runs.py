"""Goal runs of one fitted variant (ADR 0006 §2.1, stage 6).

    uv run --extra surrogate python scripts/goal_runs.py <variant>

Uses the shared-fit parameters (docs/reports/fits/<variant>.json, applied by the
registry). Writes docs/reports/goals/<variant>.json:

- goal 1: MCA strain sequence (paper protocol, `synmep.strains`); per strain the batch
  (Fig 9 conditions, warm white: best-harvest mean productivity, titer after 96 h)
  and the turbidostat optimum over density (100 µmol m^-2 s^-1, 1 cm); the best strain
  under the six lamps (batch);
- goal 2: best strain in the batch under blue/orange/red LED mixtures (quarter steps)
  and the lamps, at the same photon flux; with chromatic acclimation too if the
  family has it (`+ca`, calibrated at the reference state);
- convergence counts for gate 3.
"""

from __future__ import annotations

import itertools
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

from synmep import evaluate, light, strains, variants

warnings.filterwarnings("ignore")
OUT = Path("docs/reports/goals")
DENSITIES = np.geomspace(1e7, 3e9, 9)
CHANNELS = light.LAMPS + light.LEDS


def led_mixtures(step: float = 0.25) -> list[dict[str, float]]:
    n = round(1 / step)
    out = []
    for i, j in itertools.product(range(n + 1), repeat=2):
        if i + j <= n:
            mix = {"led_blue": i * step, "led_orange": j * step,
                   "led_red": (n - i - j) * step}  # fmt: skip
            out.append({k: v for k, v in mix.items() if v > 0})
    return out


class Counter:
    def __init__(self) -> None:
        self.ok, self.n = 0, 0

    def add(self, success: bool) -> None:
        self.n += 1
        self.ok += int(success)


def batch_summary(b: dict | None) -> dict | None:
    if b is None:
        return None
    return {k: b[k] for k in ("best_harvest_h", "best_productivity_mg_L_h",
                              "titer_end_mg_L")}  # fmt: skip


def spectra_runs(a, strain, count: Counter, pars=None) -> dict:
    out = {}
    sources = [{lamp: 1.0} for lamp in light.LAMPS] + led_mixtures()
    for mix in sources:
        label = "+".join(f"{k}:{v:g}" for k, v in mix.items())
        try:
            b = evaluate.batch(
                a, strain, light_channels=CHANNELS,
                weights=light.weights(mix, CHANNELS), pars=pars,
            )  # fmt: skip
        except Exception:  # noqa: BLE001
            b = None
        count.add(b is not None)
        out[label] = batch_summary(b)
        print(f"  {label}: {out[label]}", flush=True)
    return out


def main() -> None:
    name = sys.argv[1]
    a = variants.adapter(name)
    counts = {"mca": Counter(), "batch": Counter(), "turbidostat": Counter(),
              "spectra": Counter()}  # fmt: skip
    t0 = time.perf_counter()
    mca = strains.optimise(a)
    counts["mca"].add(not any(r.get("failed") for r in mca["rounds"]))
    print(name, "sequence", mca["sequence"], flush=True)
    per_strain = []
    for s in strains.strains(mca):
        b = evaluate.batch(a, s)
        counts["batch"].add(b is not None)
        tu = evaluate.turbidostat(a, s, DENSITIES)
        ok = tu["productivity_mg_L_h"].notna() if "productivity_mg_L_h" in tu else []
        for x in ok:
            counts["turbidostat"].add(bool(x))
        best_t = None
        if len(ok) and ok.any():
            i = int(tu["productivity_mg_L_h"].idxmax())
            best_t = {"density": float(tu.loc[i, "density"]),
                      "productivity_mg_L_h": float(tu.loc[i, "productivity_mg_L_h"]),
                      "mu_per_h": float(tu.loc[i, "mu_per_h"])}  # fmt: skip
        per_strain.append({"strain": list(s), "batch": batch_summary(b),
                           "turbidostat_best": best_t})  # fmt: skip
        print(" ", list(s), per_strain[-1]["batch"], best_t, flush=True)
    scored = [p for p in per_strain if p["batch"]]
    best = (
        max(scored, key=lambda p: p["batch"]["best_productivity_mg_L_h"])
        if scored
        else None
    )
    best_strain = {k: mca["folds"][k] for k in best["strain"]} if best else {}
    report = {
        "variant": name,
        "mca": mca,
        "strains": per_strain,
        "best_strain": list(best_strain),
        "spectra": spectra_runs(a, best_strain, counts["spectra"]),
    }
    if name.startswith("minimal-"):  # chromatic acclimation available
        ca = variants.adapter(f"{name}+ca")
        ref = evaluate.calibrate_chromatic(
            a,
            ca,
            light_channels=CHANNELS,
            weights=light.weights(light.DEFAULT, CHANNELS),
        )
        report["chromatic_reference"] = ref
        report["spectra_chromatic"] = spectra_runs(
            ca, best_strain, counts["spectra"], pars={"chromatic_reference": ref}
        )
    report["convergence"] = {k: {"ok": c.ok, "n": c.n} for k, c in counts.items()}
    report["wall_s"] = time.perf_counter() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report["convergence"]), f"{report['wall_s']:.0f} s")


if __name__ == "__main__":
    main()
