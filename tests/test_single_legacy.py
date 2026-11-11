"""`single-legacy` must reproduce the legacy modelbase model (ADR 0001 §2.4)."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from synmep.data import load_legacy_inputs
from synmep.models.single_legacy import get_model, reaction_name, sanitise

REF = Path("data/reference")
CASES = [("dshc", "steady"), ("dshc", "dynamic"), ("wt", "dynamic")]


@pytest.fixture(scope="module")
def inputs():
    return load_legacy_inputs()


@pytest.fixture(scope="module")
def models(inputs):
    return {c: get_model(inputs, *c) for c in CASES}


@pytest.mark.parametrize("case", CASES)
def test_parameters_match_legacy(models, case):
    ref = json.loads((REF / f"parameters_{case[0]}_{case[1]}.json").read_text())
    mine = models[case].get_parameter_values()
    assert mine, "no parameters"
    for name, value in mine.items():
        legacy = {sanitise(k): v for k, v in ref.items()}
        assert value == pytest.approx(legacy[name], rel=1e-12), name


@pytest.mark.parametrize("case", CASES)
def test_initial_conditions_match_legacy(models, case):
    ref = json.loads((REF / f"y0_{case[0]}_{case[1]}.json").read_text())
    mine = models[case].get_initial_conditions()
    for name, value in ref.items():
        assert mine[sanitise(name)] == pytest.approx(value, rel=1e-12), name


@pytest.mark.parametrize("case", CASES)
def test_rhs_and_fluxes_match_legacy(models, case):
    m = models[case]
    ref = pd.read_csv(REF / f"rhs_points_{case[0]}_{case[1]}.csv", index_col=0)
    for label, row in ref.iterrows():
        y = {
            sanitise(c[2:]): float(v)
            for c, v in row.items()
            if c.startswith("y:") and not np.isnan(v)
        }
        rhs = m.get_right_hand_side(y)
        fluxes = m.get_fluxes(y)
        # derivatives are sums of fluxes: allow cancellation noise relative to them
        scale = max(abs(float(v)) for c, v in row.items() if c.startswith("flux:"))
        for c, v in row.items():
            if c.startswith("dydt:"):
                got, want = rhs[sanitise(c[5:])], float(v)
            elif c.startswith("flux:"):
                got, want = fluxes[reaction_name(c[5:])], float(v)
            else:
                continue
            assert got == pytest.approx(want, rel=1e-7, abs=1e-12 * scale), (label, c)
