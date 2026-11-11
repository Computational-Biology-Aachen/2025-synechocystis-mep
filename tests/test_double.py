"""`double` family (ADR 0004): structural checks."""

import numpy as np
import pandas as pd
import pytest
import sympy as sp
from mxlpy import Simulator

from synmep.data import load_legacy_inputs
from synmep.models import double_growth, metabolism
from synmep.runners import NATIVE
from synmep.surrogates.poly import expressions, fit_polynomial, predict


def test_metabolism_reference_is_steady_state():
    m = metabolism.get_model()
    rhs = m.get_right_hand_side()
    ref = metabolism.reference_fluxes()
    assert np.max(np.abs(rhs.to_numpy())) < 1e-12 * ref["dxs"]
    assert m.get_args()["mu"] == pytest.approx(metabolism.REFERENCE.mu, rel=1e-9)


def test_polynomial_numpy_and_sympy_agree():
    rng = np.random.default_rng(0)
    x = pd.DataFrame(
        {"a": 10 ** rng.uniform(-1, 1, 200), "b": 10 ** rng.uniform(0, 2, 200)}
    )
    y = pd.DataFrame({"y": x["a"] * np.sqrt(x["b"]) / (1 + x["a"])})
    poly = fit_polynomial(x, y, 3, log_outputs=True)
    expr = expressions(poly, {"a": sp.Symbol("a"), "b": sp.Symbol("b")})["y"]
    fn = sp.lambdify([sp.Symbol("a"), sp.Symbol("b")], expr)
    probe = pd.DataFrame(
        {"a": [0.01, 0.5, 30.0], "b": [3.0, 50.0, 1e4]}
    )  # incl. clamping
    np.testing.assert_allclose(
        [fn(a, b) for a, b in probe.to_numpy()], predict(poly, probe)["y"], rtol=1e-10
    )


def test_growth_model_runs_natively():
    m = double_growth.get_model(load_legacy_inputs(), "dynamic")
    res = Simulator(m, integrator=NATIVE).simulate(600.0).get_result().unwrap_or_err()
    args = res.get_args().iloc[-1]
    assert args["mu"] > 0
    assert args["N"] > m.get_initial_conditions()["N"]
