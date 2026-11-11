"""pfennig2024 steady states under LED mixtures and lamps (ADR 0006 §2.4, §2.6).

    uv run --extra surrogate python scripts/train_petc_led.py

Same culture-state sampling as `train_petc_surrogate.py` (absorbed light from the
self-shading light model at a random depth), but every sample gets a spectrum: a
quarter one of the six lamps, the rest a Dirichlet(0.5) mixture of the blue, orange
and red LEDs and the warm white LED (many near-monochromatic cases, i.e. strongly
imbalanced PS1/PS2 excitation). Writes data/surrogates/petc_training_led.csv.
"""

from __future__ import annotations

import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import scan

# reuse the state sampling of the white-light training
from train_petc_surrogate import (  # type: ignore[import-not-found]
    DEPTH_RANGE,
    PBS_SCALE_RANGE,
    STATE_RANGES,
    light_parameters,
)

from synmep import light
from synmep.data import load_legacy_inputs
from synmep.petc.pfennig import LIGHT_INPUTS, get_petc_model, outputs
from synmep.runners import NATIVE_RUNNER
from synmep.surrogates.petc_legacy import evaluate_light

warnings.filterwarnings("ignore")
OUT = Path("data/surrogates/petc_training_led.csv")
N_SAMPLES = 4000
SEED = 1
CHANNELS = light.LAMPS + light.LEDS
MIX = ("led_blue", "led_orange", "led_red", "warm_white_led")


def spectra(n: int, rng: np.random.Generator) -> pd.DataFrame:
    w = np.zeros((n, len(CHANNELS)))
    is_lamp = rng.uniform(size=n) < 0.25
    lamps = rng.integers(len(light.LAMPS), size=n)
    mix = rng.dirichlet(np.full(len(MIX), 0.5), size=n)
    cols = [CHANNELS.index(c) for c in MIX]
    for i in range(n):
        if is_lamp[i]:
            w[i, lamps[i]] = 1.0
        else:
            w[i, cols] = mix[i]
    return pd.DataFrame(w, columns=[f"w_{c}" for c in CHANNELS])


def sample(n: int, rng: np.random.Generator) -> pd.DataFrame:
    states = pd.DataFrame(
        {k: 10 ** rng.uniform(*np.log10(r), n) for k, r in STATE_RANGES.items()}
    )
    states["sample_depth"] = rng.uniform(*DEPTH_RANGE, n)
    states["pbs_scale"] = 10 ** rng.uniform(*np.log10(PBS_SCALE_RANGE), n)
    states["ps_ratio"] = 10 ** rng.uniform(np.log10(2.0), np.log10(12.0), n)
    states["z"] = rng.uniform(0, 1, n) * states["sample_depth"]
    states = pd.concat([states, spectra(n, rng)], axis=1)
    pars = light_parameters()
    pars.pop("ps_ratio")
    states["phycocyanin"] = pars.pop("phycocyanin") * states["pbs_scale"]
    states["allophycocyanin"] = pars.pop("allophycocyanin") * states["pbs_scale"]
    absorbed = evaluate_light(
        load_legacy_inputs(), states, pars, channels=light.channels(CHANNELS)
    )
    return pd.concat([states, absorbed.set_axis(list(LIGHT_INPUTS), axis=1)], axis=1)


def main() -> None:
    rng = np.random.default_rng(SEED)
    samples = sample(N_SAMPLES, rng)
    t = time.perf_counter()
    res = scan.steady_state(
        get_petc_model(),
        to_scan=samples[list(LIGHT_INPUTS)],
        parallel=False,  # cheap points: serial is faster (report 07)
        runner=NATIVE_RUNNER,
    )
    y = outputs(res.fluxes).reset_index(drop=True)
    data = pd.concat([samples, y], axis=1)
    data.to_csv(OUT, index=False)
    ok = int(data[["J_3PGA", "J_ATP"]].notna().all(axis=1).sum())
    print(f"{ok}/{N_SAMPLES} converged in {time.perf_counter() - t:.0f} s")


if __name__ == "__main__":
    main()
