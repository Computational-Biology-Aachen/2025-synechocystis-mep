"""Shared batch / turbidostat evaluation (ADR 0006 §2.1)."""

import numpy as np
import pytest

from synmep import evaluate, variants


@pytest.fixture(scope="module")
def m0():
    return variants.adapter("minimal-m0")


def test_batch_is_plausible(m0):
    # (the stored minimal-m0 fit predates the LED-refitted PETC; the shared fit,
    # ADR 0006 stage 5, replaces it)
    b = evaluate.batch(m0)
    d = evaluate._growth_data()
    t = np.linspace(0, 96, len(d))
    mu = np.interp(t, b["t_h"], b["mu_per_h"])
    assert np.sqrt(np.mean((mu - d.to_numpy()) ** 2)) < 0.03
    assert b["N"][-1] > 5 * b["N"][0]
    assert b["titer_end_mg_L"] > 0
    assert 0 < b["best_harvest_h"] <= 96


def test_overexpression_raises_turbidostat_productivity(m0):
    dens = np.array([1e8, 1e9])
    wt = evaluate.turbidostat(m0, densities=dens)
    sqs = evaluate.turbidostat(m0, {"E_Sqs": 10.0}, densities=dens)
    assert (sqs["productivity_mg_L_h"] > wt["productivity_mg_L_h"]).all()


def test_unknown_overexpression_is_rejected(m0):
    with pytest.raises(ValueError, match="can't overexpress"):
        evaluate.batch(m0, {"E_Unknown": 2.0})


def test_sq_max_matches_dry_weight() -> None:
    from synmep.models.minimal import dry_weight_density

    assert evaluate.DW_DENSITY == pytest.approx(dry_weight_density(), rel=1e-3)
