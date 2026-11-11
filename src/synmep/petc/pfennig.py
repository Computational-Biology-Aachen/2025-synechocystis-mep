"""pfennig2024 PETC model, made steady-state capable and native-lowerable (ADR 0004).

Two changes to `mxlmodels.pfennig2024_synechocystis`:

1. Its spectral light handling (model data) depends only on parameters, so the
   data-dependent derived quantities are replaced by parameters: light absorbed per
   complex [µmol(photons) mmol(Chl)^-1 s^-1] and OCP light. These are the surrogate's
   inputs. Fluorescence readouts of the measuring light are dropped. (Until mxlpy 0.68
   this was also needed because Native couldn't lower model data, MxlPy#201.)
2. 3PGA has no sink besides respiration, so the published model has no steady state
   (3PGA grows linearly). 3PGA is clamped to its default pool, and its net production
   becomes an output. Adding a mass-action biomass drain instead (as the legacy
   surrogate's training model did) made cytoplasmic and lumenal protons accumulate
   without bound.

Outputs at steady state [mmol mol(Chl)^-1 s^-1]: `J_3PGA` (net 3PGA production) and
`J_ATP` (pfennig's lumped ATP consumption, i.e. ATP available for cellular demand).
"""

from __future__ import annotations

import pandas as pd
from mxlmodels.pfennig2024_synechocystis import get_pfennig2024_synechocystis
from mxlpy import Model

LIGHT_INPUTS = ("absorption_ps1", "absorption_ps2", "absorption_pbs", "light_OCP")
OUTPUTS = ("J_3PGA", "J_ATP")
PGA_PRODUCTION = {"vCBB": 1.0, "vOxy": -2.0, "vPRsalv": 1.0, "vRespiration": -1.0}


def get_petc_model(*, clamp_pga: bool = True) -> Model:
    """pfennig2024 with absorbed light as parameters; 3PGA clamped unless `clamp_pga`
    is False (the `full` family couples 3PGA to the metabolism instead)."""
    m = get_pfennig2024_synechocystis()
    defaults = m.get_args()

    data = set(m.get_raw_data())
    tainted = set(data)
    derived = m.get_raw_derived()
    changed = True
    while changed:
        changed = False
        for name, d in derived.items():
            if name not in tainted and set(d.args) & tainted:
                tainted.add(name)
                changed = True
    for name in sorted(tainted - data):
        m.remove_derived(name)
    for name in data:
        m.remove_data(name)
    m.add_parameters({name: float(defaults[name]) for name in LIGHT_INPUTS})

    if not clamp_pga:
        return m
    # clamp 3PGA: drop it from all stoichiometries, keep it as a constant pool
    pga0 = float(defaults["3PGA"])
    for rxn, r in m.get_raw_reactions().items():
        if "3PGA" in r.stoichiometry:
            stoich = {k: v for k, v in r.stoichiometry.items() if k != "3PGA"}
            m.update_reaction(rxn, stoichiometry=stoich)
    m.remove_variable("3PGA")
    m.add_parameter("3PGA", pga0)
    return m


def outputs(fluxes: pd.DataFrame | pd.Series) -> pd.DataFrame | pd.Series:
    """PETC outputs from steady-state fluxes (rows = scan points or a single Series)."""
    j_pga = sum(c * fluxes[r] for r, c in PGA_PRODUCTION.items())
    if isinstance(fluxes, pd.Series):
        return pd.Series({"J_3PGA": j_pga, "J_ATP": fluxes["vATPconsumption"]})
    return pd.DataFrame({"J_3PGA": j_pga, "J_ATP": fluxes["vATPconsumption"]})
