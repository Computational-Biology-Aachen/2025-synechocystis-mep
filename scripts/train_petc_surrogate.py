"""Train the PETC surrogate from the clamped, native-lowerable pfennig2024 (ADR 0004).

    uv run --extra surrogate python scripts/train_petc_surrogate.py

Training inputs are *physically realisable* absorbed-light combinations: culture states
(cell density, Chl, carotenoids, light, sample depth, position in depth) are sampled
and pushed through the self-shading light model (local light at that depth) with the growth model's pigment parameters. Sampling the
four absorbed-light inputs independently left combinations the growth model produces
(e.g. low PBS:PS1 at 4 fL) outside the training data, where the polynomial
extrapolated to negative supply. pfennig2024 is solved to steady state with the native
runner and a polynomial in the (clamped) log inputs is fitted.
Writes data/surrogates/petc_training.csv, data/surrogates/petc_poly.json and
docs/reports/petc-surrogate.json.
"""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import scan

from synmep.data import load_legacy_inputs
from synmep.models.g_growth import get_model as get_growth_model
from synmep.petc.pfennig import LIGHT_INPUTS, OUTPUTS, get_petc_model, outputs
from synmep.runners import NATIVE_RUNNER
from synmep.surrogates.petc_legacy import evaluate_light
from synmep.surrogates.poly import fit_polynomial, predict

warnings.filterwarnings("ignore")
OUT = Path("data/surrogates")
N_SAMPLES = 4000
DEGREE = 5
SEED = 0

# culture states sampled log-uniformly (depth uniformly)
STATE_RANGES = {
    "cell_density": (1e6 * 1e-6, 3e9 * 1e-6),  # [cells nL^-1] (1e6-3e9 cells mL^-1)
    "chl": (2.0, 60.0),  # [mM]
    "carotenoids": (0.5, 30.0),  # [mM]
    "light": (3.0, 1500.0),  # [µmol(photons) m^-2 s^-1]
}
DEPTH_RANGE = (0.005, 0.03)  # [m]
PBS_SCALE_RANGE = (
    0.25,
    4.0,
)  # [-] phycobiliprotein content factor (fitted in g_growth)


def light_parameters() -> dict[str, float]:
    m = get_growth_model(load_legacy_inputs())
    m.update_parameter("pbs_scale", 1.0)
    a = m.get_args()
    keys = ("phycocyanin", "allophycocyanin", "ps_ratio", "beta_carotene_fraction",
            "cell_volume")  # fmt: skip
    return {k: float(a[k]) for k in keys}


def sample(n: int, rng: np.random.Generator) -> pd.DataFrame:
    states = pd.DataFrame(
        {k: 10 ** rng.uniform(*np.log10(r), n) for k, r in STATE_RANGES.items()}
    )
    states["sample_depth"] = rng.uniform(*DEPTH_RANGE, n)
    states["pbs_scale"] = 10 ** rng.uniform(*np.log10(PBS_SCALE_RANGE), n)
    # local light at a random depth (the growth model evaluates the PETC at depth nodes)
    states["z"] = rng.uniform(0, 1, n) * states["sample_depth"]
    pars = light_parameters()
    states["phycocyanin"] = pars.pop("phycocyanin") * states["pbs_scale"]
    states["allophycocyanin"] = pars.pop("allophycocyanin") * states["pbs_scale"]
    light = evaluate_light(load_legacy_inputs(), states, pars)
    return pd.concat([states, light.set_axis(list(LIGHT_INPUTS), axis=1)], axis=1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    samples = sample(N_SAMPLES, rng)
    x = samples[list(LIGHT_INPUTS)]
    t = time.perf_counter()
    res = scan.steady_state(
        get_petc_model(),
        to_scan=x,
        parallel=False,  # the loky-parallel scan hung (ADR 0004); serial takes ~6 min
        runner=NATIVE_RUNNER,
    )
    wall = time.perf_counter() - t
    y = outputs(res.fluxes).reset_index(drop=True)
    data = pd.concat([samples, y], axis=1)
    ok = data.notna().all(axis=1)
    data.to_csv(OUT / "petc_training.csv", index=False)
    data = data[ok].reset_index(drop=True)

    idx = rng.permutation(len(data))
    n_test = len(data) // 5
    test, train = data.iloc[idx[:n_test]], data.iloc[idx[n_test:]]
    poly = fit_polynomial(train[list(LIGHT_INPUTS)], train[list(OUTPUTS)], DEGREE)
    report = {
        "n_samples": N_SAMPLES,
        "n_converged": int(ok.sum()),
        "wall_s": wall,
        "degree": DEGREE,
        "state_ranges": STATE_RANGES,
        "depth_range": DEPTH_RANGE,
        "pbs_scale_range": PBS_SCALE_RANGE,
        "n_terms": len(poly["powers"]),
    }
    pred = predict(poly, test[list(LIGHT_INPUTS)])
    for o in OUTPUTS:
        err = pred[o].to_numpy() - test[o].to_numpy()
        scale = np.abs(test[o]).max()
        report[o] = {
            "rmse": float(np.sqrt(np.mean(err**2))),
            "max_abs_err": float(np.abs(err).max()),
            "max_abs_err_rel_to_range": float(np.abs(err).max() / scale),
            "r2": float(1 - np.sum(err**2) / np.sum((test[o] - test[o].mean()) ** 2)),
        }
    (OUT / "petc_poly.json").write_text(json.dumps(poly, indent=1))
    Path("docs/reports/petc-surrogate.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
