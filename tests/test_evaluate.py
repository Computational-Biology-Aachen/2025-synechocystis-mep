"""Shared batch / turbidostat evaluation (ADR 0006 §2.1)."""

import json
from pathlib import Path

import numpy as np
import pytest

from synmep import evaluate, variants


@pytest.fixture(scope="module")
def m0():
    return variants.adapter("minimal-m0")


def test_batch_reproduces_the_fit(m0):
    b = evaluate.batch(m0)
    fit = json.loads(Path("docs/reports/minimal-m0-fit.json").read_text())
    d = evaluate._growth_data()
    t = np.linspace(0, 96, len(d))
    mu = np.interp(t, b["t_h"], b["mu_per_h"])
    rmse = np.sqrt(np.mean((mu - d.to_numpy()) ** 2))
    assert rmse == pytest.approx(fit["rmse_mu_per_h"], rel=0.02)
    assert b["titer_end_mg_L"] > 0
    assert 0 < b["best_harvest_h"] <= 96


def test_overexpression_raises_turbidostat_productivity(m0):
    dens = np.array([1e8, 1e9])
    wt = evaluate.turbidostat(m0, densities=dens)
    sqs = evaluate.turbidostat(m0, {"E_Sqs": 10.0}, densities=dens)
    assert (sqs["productivity_mg_L_h"] > wt["productivity_mg_L_h"]).all()


def test_unknown_overexpression_is_rejected():
    with pytest.raises(ValueError, match="can't overexpress"):
        evaluate.batch(variants.adapter("double-none"), {"E_Sqs": 2.0})
