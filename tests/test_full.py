"""`full` family (ADR 0006 §2.5): structure and a dilute steady state."""

from synmep import evaluate, variants
from synmep.data import load_legacy_inputs
from synmep.models import full


def test_structure():
    m = full.get_model(load_legacy_inputs(), layers=3)
    names = m.get_variable_names()
    # 3 pfennig2024 copies with free 3PGA, per-layer FPP/GGPP, shared slow pools
    assert {"z0_3PGA", "z2_ATP", "z1_FPP", "MEPi", "IPP", "Chl", "SQ"} <= set(names)
    assert "z0_vATPconsumption" not in m.get_reaction_names()


def test_dilute_culture_grows():
    tu = evaluate.turbidostat(variants.adapter("full-k3"), densities=[1e7])
    # unfitted (shared fit pending): grows, below the data (0.04-0.06 h^-1)
    assert 0.005 < tu["mu_per_h"].iloc[0] < 0.06
