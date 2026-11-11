"""Shared fit of one variant (ADR 0006 §2.8, §5).

    uv run --extra surrogate python scripts/shared_fit.py <variant> [--quick]

Targets (one protocol for every variant):
- Fig 9 batch: growth rate (48 points, sigma 0.005 h^-1) and ln N (sigma 0.05);
- WT Chl per cell at 4, 24, 48, 96 h (proteomics batch, `synmep.pigments`; sigma 10 %,
  since it is a different culture);
- WT phycocyanin : Chl mass ratio at the same times (sigma 20 %), which pins
  `pbs_scale` (the batch data only pin petc_scale x pbs_scale, report 08).

Held out: the proteomics batch's cell densities (WT, from its 4 h sample to 96 h);
Fig 10 is reported as a check. Fit factors are multiplicative on the variant's own
parameter values (`SPECS`), within bounds, in log space: least squares from the best
points of a small start grid. Writes docs/reports/fits/<variant>.json, which the
registry (`synmep.variants.fit_parameters`) then applies.
"""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from synmep import evaluate, pigments, variants
from synmep.data import proteomics_samples

warnings.filterwarnings("ignore")
OUT = Path("docs/reports/fits")
SIGMA_MU, SIGMA_LN_N, SIGMA_CHL, SIGMA_PC = 0.005, 0.05, 0.10, 0.20
TIMES = (4, 24, 48, 96)  # [h] proteomics samples

# factor name -> (model parameters it scales, lower, upper factor)
PBS = (("pbs_scale",), 0.25, 8.0)
CAP = {
    "capacity_K": (("capacity_K",), 0.01, 100.0),
    "capacity_tau": (("capacity_tau",), 0.05, 10.0),
}
SPECS = {
    "single-fixed": {
        "chl_growth_slope": (("chl_growth_slope",), 0.1, 10.0),
        "pbs_scale": (("phycocyanin", "allophycocyanin"), 0.25, 8.0),
        "atp_scale": (("atp_scale",), 0.2, 5.0),
    },
    "double": {
        "petc_scale": (("petc_scale",), 0.1, 10.0),
        # the PETC surrogate is trained for pbs_scale 0.25-4 (ADR 0004)
        "pbs_scale": (("pbs_scale",), 0.25, 4.0),
    },
    "minimal": {
        "petc_scale": (("petc_scale",), 0.1, 10.0),
        "pbs_scale": PBS,
        "pigment_coupling": (("pigment_coupling",), 0.02, 3.0),
    },
    "full": {
        "pbs_scale": PBS,
        "pigment_coupling": (("pigment_coupling",), 0.02, 10.0),
        "V_bio": (("V_bio",), 0.2, 20.0),
        "vCBB_max": (("vCBB_max",), 0.2, 5.0),
    },
}


def spec(name: str) -> dict:
    family = name.split("-")[0]
    key = "single-fixed" if name.startswith("single-fixed") else family
    s = dict(SPECS[key])
    if name.endswith("-cap"):
        s |= CAP
    return s


def targets() -> dict:
    t = pigments.pigment_targets().loc["WT"]
    s = proteomics_samples()
    s["strain"] = s["sample"].str.split().str[0]
    s["t"] = s["sample"].str.split().str[1].astype(int)
    cells = s[s["strain"] == "WT"].groupby("t")["cells"].median()
    return {
        "chl": t.loc[list(TIMES), "chl_mM"].to_numpy(),
        "pc": t.loc[list(TIMES), "pc_per_chl"].to_numpy(),
        "cells": cells.loc[list(TIMES)].to_numpy(),
    }


def at(b: dict, key: str, hours) -> np.ndarray:
    return np.interp(hours, b["t_h"], b[key])


class Problem:
    def __init__(self, name: str) -> None:
        self.name = name
        self.a = variants.adapter(name)
        self.spec = spec(name)
        m = self.a.build("dynamic", **evaluate.BATCH_LIGHT)
        p = m.get_parameter_values()
        self.base = {t: p[t] for names, _, _ in self.spec.values() for t in names}
        self.lo = np.log([v[1] for v in self.spec.values()])
        self.hi = np.log([v[2] for v in self.spec.values()])
        self.data = evaluate._growth_data()
        self.t_data = np.linspace(0, 96, len(self.data))
        self.tg = targets()
        self.cache: dict[tuple, np.ndarray] = {}

    def pars(self, z: np.ndarray) -> dict[str, float]:
        out = {}
        for (names, _, _), f in zip(self.spec.values(), np.exp(z), strict=True):
            out |= {t: self.base[t] * f for t in names}
        return out

    def residuals(self, z: np.ndarray) -> np.ndarray:
        key = tuple(np.round(z, 10))
        if key in self.cache:
            return self.cache[key]
        n_res = 2 * len(self.data) + 2 * len(TIMES)
        try:
            b = evaluate.batch(self.a, pars=self.pars(z))
        except Exception:  # noqa: BLE001 - any failure is a bad point
            b = None
        if b is None:
            r = np.full(n_res, 1e3)
        else:
            mu = at(b, "mu_per_h", self.t_data)
            n = at(b, "N", self.t_data)
            r = np.concatenate([
                (mu - self.data.to_numpy()) / SIGMA_MU,
                (np.log(n) - np.log(self.data.index.to_numpy())) / SIGMA_LN_N,
                (at(b, "chl_mM", TIMES) / self.tg["chl"] - 1) / SIGMA_CHL,
                (at(b, "pc_per_chl", TIMES) / self.tg["pc"] - 1) / SIGMA_PC,
            ])  # fmt: skip
            r = np.nan_to_num(r, nan=1e3, posinf=1e3, neginf=-1e3)
        self.cache[key] = r
        return r


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("variant")
    parser.add_argument("--quick", action="store_true", help="small grid, few steps")
    parser.add_argument("--start", help="JSON of factors to start from (no grid)")
    parser.add_argument("--max-nfev", type=int)
    args = parser.parse_args()
    prob = Problem(args.variant)
    k = len(prob.spec)
    rng = np.random.default_rng(0)
    n_grid, n_starts, max_nfev = (6, 1, 15) if args.quick else (16, 2, 40)
    t0 = time.perf_counter()
    if args.max_nfev:
        max_nfev = args.max_nfev
    grid = [np.clip(np.zeros(k), prob.lo, prob.hi)]
    if args.start:
        f0 = json.loads(args.start)
        grid = [np.log([f0[n] for n in prob.spec])]
        n_starts = 1
    else:
        grid += [
            prob.lo + rng.uniform(size=k) * (prob.hi - prob.lo) for _ in range(n_grid)
        ]
    costs = [float(np.sum(prob.residuals(z) ** 2)) for z in grid]
    print(f"{args.variant}: grid {time.perf_counter() - t0:.0f} s, best cost "
          f"{min(costs):.4g}", flush=True)  # fmt: skip
    starts = [grid[i] for i in np.argsort(costs)[:n_starts]]
    best0 = dict(zip(prob.spec, map(float, np.exp(starts[0])), strict=True))
    print(f"{args.variant}: best grid factors {json.dumps(best0)}", flush=True)
    fits = [
        least_squares(
            prob.residuals,
            np.clip(z, prob.lo + 1e-9, prob.hi - 1e-9),
            bounds=(prob.lo, prob.hi),
            diff_step=1e-2,
            max_nfev=max_nfev,
        )
        for z in starts
    ]
    fit = min(fits, key=lambda f: f.cost)
    pars = prob.pars(fit.x)
    b = evaluate.batch(prob.a, pars=pars)
    report = {
        "variant": args.variant,
        "factors": dict(zip(prob.spec, map(float, np.exp(fit.x)), strict=True)),
        "fitted": pars,
        "cost": float(fit.cost),
        "n_parameters": k,
        "wall_s": time.perf_counter() - t0,
    }
    if b is not None:
        mu = at(b, "mu_per_h", prob.t_data)
        n = at(b, "N", prob.t_data)
        n_obs = 2 * len(prob.data) + 2 * len(TIMES)
        report |= {
            "rmse_mu_per_h": float(np.sqrt(np.mean((mu - prob.data.to_numpy()) ** 2))),
            "rmse_ln_N": float(
                np.sqrt(np.mean((np.log(n) - np.log(prob.data.index.to_numpy())) ** 2))
            ),
            "chl_rel_err": (at(b, "chl_mM", TIMES) / prob.tg["chl"] - 1).tolist(),
            "pc_rel_err": (at(b, "pc_per_chl", TIMES) / prob.tg["pc"] - 1).tolist(),
            "aic": float(n_obs * np.log(2 * fit.cost / n_obs) + 2 * k),
            "trajectory": {key: b[key] for key in ("t_h", "mu_per_h", "N", "chl_mM",
                                                    "pc_per_chl", "sq_mM")},
        }  # fmt: skip
        # held out: proteomics batch (WT), from its 4 h sample
        v = evaluate.batch(prob.a, pars=pars, n0=float(prob.tg["cells"][0]), hours=92)
        if v is not None:
            pred = at(v, "N", np.array(TIMES[1:]) - TIMES[0])
            err = np.log(pred) - np.log(prob.tg["cells"][1:])
            report["validation_rmse_ln_N"] = float(np.sqrt(np.mean(err**2)))
        try:
            f10 = evaluate.fig10(prob.a, pars=pars)
            report["fig10_rmse_mu_per_h"] = f10["rmse_mu_per_h"]
            report["fig10_n_converged"] = f10["n_converged"]
        except Exception as e:  # noqa: BLE001
            report["fig10_error"] = repr(e)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{args.variant}.json").write_text(json.dumps(report, indent=1))
    print(json.dumps({k: v for k, v in report.items() if k != "trajectory"}, indent=1))


if __name__ == "__main__":
    main()
