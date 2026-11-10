import numpy as np
import pandas as pd
import pytest
import sympy as sp

from synmep.data import load_legacy_inputs
from synmep.surrogates.petc_legacy import polynomial_expressions

REF = "data/reference/surrogate_check.csv"


@pytest.fixture(scope="module")
def inputs():
    return load_legacy_inputs()


def test_polynomial_matches_sklearn(inputs):
    poly = inputs.surrogate_poly
    names = [str(n) for n in poly["feature_names"]]
    syms = {n: sp.Symbol(n) for n in names}
    exprs = polynomial_expressions(poly, syms)
    df = pd.read_csv(REF)
    order = ["vATP_prod", "vNADPH_prod", "v3PGA_prod", "vFdred_prod"]
    fns = [sp.lambdify([syms[n] for n in names], exprs[o]) for o in order]
    for col, fn in zip(["ATP", "NADPH", "3PGA", "Fd_red"], fns, strict=True):
        got = np.array([fn(*row) for row in df[names].to_numpy()])
        want = np.maximum(df[f"out_{col}"].to_numpy(), 0)
        np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-9 * np.abs(want).max())
