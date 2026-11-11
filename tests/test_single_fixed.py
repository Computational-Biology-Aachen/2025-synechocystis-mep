"""`single-fixed` (ADR 0002 §5): checks of the individual fixes and of the assembled model."""

import math

import pytest
from mxlpy import Simulator
from mxlpy.integrators import Sundials

from synmep import kinetics_single_fixed as k
from synmep.models import single_fixed
from synmep.runners import NATIVE


@pytest.fixture(scope="module")
def inputs():
    return single_fixed.load_inputs()


@pytest.fixture(scope="module")
def model(inputs):
    return single_fixed.get_model(inputs, "dshc", "steady")


def test_b4_monomerisation_only_of_non_active_dxs():
    # no inhibitors, saturating ThiDP: (almost) all DXS active -> ~no monomerisation
    assert k.DXS_na(1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1e9, 1.0) == pytest.approx(
        0, abs=1e-8
    )
    # IPP at its Kd: half of the enzyme is inactive
    assert k.DXS_na(1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1e9, 1.0) == pytest.approx(0.5)


def test_b5_sqs_rate_finite_without_squalene():
    v = k.v_sqs(0.01, 0.01, 1.0, 0.0, 1.0, 1.0, 0.003, 0.003, 0.5, 0.1, 1e10, 1.0)
    assert math.isfinite(v)
    assert v > 0


def test_b7_ordered_rate_laws_defined_for_knockouts():
    args = (1.0, 1.0, 0.1, 0.1, 1.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 1e3)
    assert k.reversible_ordered_2_comp1(*args) == 0.0


def test_b6_kinetic_defaults_in_mm(inputs):
    p = inputs.enzyme_parameters
    assert p["DXSIPPKi"] == pytest.approx(1.0)  # 1000 µM default
    assert p["IspDATPKm"] == pytest.approx(0.1)  # 100 µM default


def test_b9_cell_volume(model):
    assert model.get_parameter_values()["cell_volume"] == 4e-15


def test_enzymes_in_mm(inputs, model):
    y0 = model.get_initial_conditions()
    assert y0["DXS"] == pytest.approx(inputs.enzyme_concentrations_umol["DXS"] * 1e-3)


def test_steady_state(inputs, model):
    # Native's steady-state search gives up on single-fixed (all methods) where Sundials
    # converges; the native runner is still used for time courses
    y0 = model.get_initial_conditions() | {
        "Cell_density": inputs.exp_growth_rate_by_density.index[24] / 1e6
    }
    res = (
        Simulator(model, y0=y0, integrator=Sundials)
        .simulate_to_steady_state(oscillation_detector=None)
        .get_result()
        .unwrap_or_err()
    )
    args = res.get_args().iloc[-1]
    assert args["Growth_rate"] > 0
    for pool in ("PGA3", "ATP", "NADPH", "Fdred"):
        assert args[pool] > 0
    # Chl is pinned by the empirical growth law: Chl* = 1 / (slope * V)
    p = model.get_parameter_values()
    chl_star = 1 / (p["chl_growth_slope"] * p["cell_volume"])
    assert args["chlorophyll_phytol"] == pytest.approx(chl_star, rel=1e-6)


def test_native_time_course(inputs):
    m = single_fixed.get_model(inputs, "dshc", "dynamic")
    res = Simulator(m, integrator=NATIVE).simulate(600.0).get_result().unwrap_or_err()
    assert res.get_args().iloc[-1]["Growth_rate"] > 0
