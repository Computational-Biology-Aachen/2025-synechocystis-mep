"""Light sources as spectral channels (ADR 0006 §2.1, §2.7).

Every spectrum is a photon-flux distribution on the 1 nm grid 400–700 nm, normalised
to integral 1, so `light` [µmol(photons) m^-2 s^-1] is the total PAR flux and a mixture
sum_c w_c S_c with sum w_c = 1 keeps it. The lamp spectra are the Fuente2021 curves in
`data/lights/` (assumed to be photon distributions, as in pfennig2024).

LEDs are Gaussians (`LED_FWHM`). The model's absorption data end at 700 nm, so a far-red
LED (730 nm) would barely be absorbed in the model; it's defined but not in the default
LED set (`LEDS`), see ADR 0006.

Models take `channels` (name -> spectrum); the weights are parameters `w_<name>`.
"""

from __future__ import annotations

from functools import cache

import numpy as np
import pandas as pd
from scipy.integrate import simpson

from synmep.data import DATA

WAVELENGTHS = np.arange(400, 701)  # [nm]
LAMPS = (
    "warm_white_led",
    "cool_white_led",
    "fluorescent_lamp",
    "halogen_lamp",
    "incandescent_bulb",
    "solar",
)
LED_PEAKS = {"blue": 450.0, "orange": 620.0, "red": 680.0, "farred": 730.0}  # [nm]
LEDS = ("led_blue", "led_orange", "led_red")
LED_FWHM = 20.0  # [nm]
DEFAULT = "warm_white_led"


def _normalise(y: np.ndarray) -> np.ndarray:
    return y / simpson(y, x=WAVELENGTHS)


@cache
def lamp(name: str) -> np.ndarray:
    s = pd.read_csv(DATA / "lights" / f"{name}.csv", index_col=0).squeeze(axis=1)
    return _normalise(s.reindex(WAVELENGTHS).to_numpy(dtype=float))


def led(peak: float, fwhm: float = LED_FWHM) -> np.ndarray:
    """Gaussian LED, normalised over 400–700 nm (only the part inside counts)."""
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return _normalise(np.exp(-0.5 * ((WAVELENGTHS - peak) / sigma) ** 2))


def channels(names: tuple[str, ...] = LAMPS + LEDS) -> dict[str, np.ndarray]:
    """Spectra by channel name (lamp names, `led_<colour>`)."""
    out = {}
    for n in names:
        out[n] = (
            led(LED_PEAKS[n.removeprefix("led_")]) if n.startswith("led_") else lamp(n)
        )
    return out


def weights(source: str | dict[str, float], names: tuple[str, ...]) -> dict[str, float]:
    """Parameter values `w_<channel>` for one source name or a mixture {name: weight}."""
    mix = {source: 1.0} if isinstance(source, str) else dict(source)
    unknown = set(mix) - set(names)
    if unknown:
        msg = f"not a channel of this model: {sorted(unknown)}"
        raise ValueError(msg)
    return {f"w_{n}": float(mix.get(n, 0.0)) for n in names}


def add_parameters(model, names: tuple[str, ...]) -> None:
    """Add the channel weights `w_<name>` (default: warm white LED, else the first)."""
    source = DEFAULT if DEFAULT in names else names[0]
    model.add_parameters(weights(source, names))
