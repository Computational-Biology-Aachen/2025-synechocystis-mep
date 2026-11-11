"""Light sources as channels (ADR 0006 §2.7)."""

import numpy as np
import pytest
from scipy.integrate import simpson

from synmep import light
from synmep.data import load_legacy_inputs
from synmep.models import minimal
from synmep.runners import NATIVE_RUNNER


def test_channels_are_normalised():
    for name, s in light.channels(light.LAMPS + light.LEDS).items():
        assert simpson(s, x=light.WAVELENGTHS) == pytest.approx(1, rel=1e-6), name
        assert np.all(s >= 0)


def test_weights_reject_unknown_channels():
    with pytest.raises(ValueError, match="not a channel"):
        light.weights("led_green", light.LEDS)


def test_warm_white_channel_reproduces_the_baked_in_spectrum():
    inp = load_legacy_inputs()
    names = (light.DEFAULT, "led_red")
    baked = NATIVE_RUNNER.steady_state(minimal.get_model(inp))
    mixed = NATIVE_RUNNER.steady_state(minimal.get_model(inp, light_channels=names))
    a = baked.value.get_args().iloc[-1]
    b = mixed.value.get_args().iloc[-1]
    # equal up to the renormalisation of the channel spectra (~1e-4)
    assert b["mu"] == pytest.approx(a["mu"], rel=1e-3)
