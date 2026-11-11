"""Paired test of the PETC forms: the same culture states under different spectra.

    uv run --extra surrogate python scripts/petc_paired_spectra.py

R² over all samples is dominated by light intensity. Goal 2 needs the *ranking of
spectra*: for each culture state (sampled as in `train_petc_led.py`, same incident
photon flux) J_3PGA under each source relative to warm white. Writes
data/surrogates/petc_paired.csv.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from mxlpy import scan
from train_petc_led import CHANNELS, sample  # type: ignore[import-not-found]

from synmep import light
from synmep.data import load_legacy_inputs
from synmep.petc.pfennig import LIGHT_INPUTS, get_petc_model, outputs
from synmep.runners import NATIVE_RUNNER
from synmep.surrogates.petc_legacy import evaluate_light

warnings.filterwarnings("ignore")
OUT = Path("data/surrogates/petc_paired.csv")
N_STATES = 300
SOURCES = {
    "warm_white_led": {"warm_white_led": 1.0},
    "cool_white_led": {"cool_white_led": 1.0},
    "solar": {"solar": 1.0},
    "led_blue": {"led_blue": 1.0},
    "led_orange": {"led_orange": 1.0},
    "led_red": {"led_red": 1.0},
    "orange+red": {"led_orange": 0.5, "led_red": 0.5},
    "blue+red": {"led_blue": 0.5, "led_red": 0.5},
}


def main() -> None:
    states = sample(N_STATES, np.random.default_rng(2))
    states = states.drop(
        columns=[c for c in states if c.startswith("w_")] + list(LIGHT_INPUTS)
    )
    rows = []
    for name, mix in SOURCES.items():
        w = light.weights(mix, CHANNELS)
        s = states.assign(**w)
        pars = {}
        absorbed = evaluate_light(
            load_legacy_inputs(),
            s,
            pars | _light_pars(),
            channels=light.channels(CHANNELS),
        ).set_axis(list(LIGHT_INPUTS), axis=1)
        rows.append(s.assign(source=name, state=np.arange(N_STATES)).join(absorbed))
    data = pd.concat(rows, ignore_index=True)
    res = scan.steady_state(
        get_petc_model(), to_scan=data[list(LIGHT_INPUTS)], parallel=False,
        runner=NATIVE_RUNNER,
    )  # fmt: skip
    data = data.join(outputs(res.fluxes).reset_index(drop=True))
    data.to_csv(OUT, index=False)
    piv = data.pivot(index="state", columns="source", values="J_3PGA")
    rel = piv.div(piv["warm_white_led"], axis=0)
    print(rel.median().round(3).to_string())


def _light_pars() -> dict[str, float]:
    from train_petc_surrogate import light_parameters  # type: ignore[import-not-found]

    p = light_parameters()
    for k in ("phycocyanin", "allophycocyanin", "ps_ratio"):
        p.pop(k)
    return p


if __name__ == "__main__":
    main()
