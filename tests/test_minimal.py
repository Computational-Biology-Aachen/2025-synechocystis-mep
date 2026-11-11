"""`minimal-m0` (ADR 0005): structural checks."""

import copy

import pytest

from synmep.data import load_legacy_inputs, proteomics_samples
from synmep.models import metabolism, minimal
from synmep.runners import NATIVE_RUNNER


@pytest.fixture(scope="module")
def inp():
    return load_legacy_inputs()


def _steady(m) -> dict[str, float]:
    res = NATIVE_RUNNER.steady_state(m)
    assert hasattr(res.value, "fluxes"), "no steady state"
    return res.value.get_args(include_derived_parameters=True).iloc[-1].to_dict()


def test_dry_weight_from_protein_is_physiological():
    samples = proteomics_samples()
    # the legacy 3.88 pg DW cell^-1 is the median gravimetric dry weight
    assert samples["dry_weight_pg"].median() == pytest.approx(3.883, rel=1e-3)
    rho = minimal.dry_weight_density()  # [g L^-1]
    assert 300 < rho < 550
    assert rho < metabolism.RHO_DW / 2


def test_proteome_fractions_are_small(inp):
    f = minimal.proteome_fractions(inp)
    assert set(f) == set(minimal.LUMPS)
    assert all(0 < v < 0.01 for v in f.values())


def test_enzyme_cost_slows_growth(inp):
    base = minimal.get_model(inp, cost="proteome", use_fit=False)
    free = minimal.get_model(inp, cost="none", use_fit=False)
    wt = _steady(copy.deepcopy(base))
    for m in (base, free):
        m.update_parameter("E_Sqs", 1000.0)
    costly, no_cost = _steady(base), _steady(free)
    assert wt["protein_extra"] == 0
    assert costly["protein_extra"] == pytest.approx(
        999 * minimal.proteome_fractions(inp)["Sqs"]
    )
    assert costly["mu"] < no_cost["mu"]


def test_slow_speedup_keeps_steady_state(inp):
    m = minimal.get_model(inp, cost="proteome", use_fit=False)
    fast = _steady(copy.deepcopy(m))
    m.update_parameter("slow_speedup", 10.0)
    slow = _steady(m)
    for k in ("mu", "Chl", "SQ", "PGA"):
        assert slow[k] == pytest.approx(fast[k], rel=1e-5)


def test_chromatic_acclimation_moves_ps_ratio_the_right_way(inp):
    from synmep import evaluate, light, variants

    names = light.LAMPS + light.LEDS
    fx, ch = variants.adapter("minimal-m0"), variants.adapter("minimal-m0+ca")
    ref = evaluate.calibrate_chromatic(fx, ch, light_channels=names)
    ratios = {}
    for src in ("warm_white_led", "led_blue", "led_orange"):
        m = evaluate.steady_model(
            ch, {}, 2e8, light_channels=names, pars={"chromatic_reference": ref},
            **evaluate.BATCH_LIGHT,
        )  # fmt: skip
        m.update_parameters(light.weights(src, names))
        ratios[src] = float(evaluate.steady_state(ch, m)["ps_ratio"])
    # warm white keeps the measured ratio at the reference state; PS1 light (blue)
    # lowers it, PS2 light (orange) raises it
    assert ratios["warm_white_led"] == pytest.approx(5.9, rel=0.02)
    assert ratios["led_blue"] < ratios["warm_white_led"] < ratios["led_orange"]
