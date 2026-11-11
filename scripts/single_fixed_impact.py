"""Before/after impact of the `single-fixed` fixes (ADR 0002 §5).

    uv run python scripts/single_fixed_impact.py [single-legacy|single-fixed]

Reports, for strain dshc:
- steady state at the default condition (light 100, density = data row 24): growth
  rate, Chl, carotenoids, SQ, cofactor pools
- the Fig 9 batch protocol (pre-culture steady state at data row 24, diluted to the
  first data point, 80 µE, 13.5 mm, 96 h): RMSE of growth rate and ln N vs. data
"""

from __future__ import annotations

import json
import sys
import time
import warnings

import numpy as np
from mxlpy import Simulator
from mxlpy.integrators import Sundials

from synmep.models import single_fixed, single_legacy
from synmep.runners import NATIVE

warnings.filterwarnings("ignore")

POOLS = ("PGA3", "ATP", "NADPH", "Fdred", "chlorophyll_phytol", "carotenoids", "SQ")


def _steady_state(model, y0):
    for integ in (NATIVE, Sundials):
        res = (
            Simulator(model, y0=y0, integrator=integ)
            .simulate_to_steady_state(oscillation_detector=None)
            .get_result()
        )
        if hasattr(res.value, "variables"):
            return res.value
    return None


def impact(module) -> dict:
    inp = module.load_inputs() if hasattr(module, "load_inputs") else None
    if inp is None:
        from synmep.data import load_legacy_inputs

        inp = load_legacy_inputs()
    data = inp.exp_growth_rate_by_density
    out: dict = {}

    ms = module.get_model(inp, "dshc", "steady")
    y0 = ms.get_initial_conditions() | {"Cell_density": data.index[24] / 1e6}
    t = time.perf_counter()
    ss = _steady_state(ms, y0)
    out["ss_wall_s"] = time.perf_counter() - t
    if ss is None:
        out["ss"] = None
        return out
    a = ss.get_args().iloc[-1]
    out["ss"] = {"mu_per_h": float(a["Growth_rate"] * 60)} | {
        p: float(a[p]) for p in POOLS
    }
    y0_ss = ss.variables.iloc[-1][ms.get_variable_names()].to_dict()

    md = module.get_model(inp, "dshc", "dynamic")
    pars = {"light": 80, "sample_depth_m": 0.0135, "Shc_ksynth": 0, "DXS_inh_str": 0}
    md.update_parameters(
        {k: v for k, v in pars.items() if k in md.get_parameter_names()}
    )
    tp = np.linspace(0, 96 * 60, len(data))
    y0_9 = y0_ss | {"Cell_density": data.index[0] / 1e6}
    tc = Simulator(md, y0=y0_9, integrator=NATIVE).simulate_time_course(tp).get_result()
    if not hasattr(tc.value, "variables"):
        out["fig9"] = None
        return out
    r = tc.value.get_args()
    mu = r["Growth_rate"].to_numpy() * 60
    n = r["Cell_density"].to_numpy() * 1e6
    out["fig9"] = {
        "rmse_mu_per_h": float(np.sqrt(np.mean((mu - data.to_numpy()) ** 2))),
        "rmse_ln_N": float(np.sqrt(np.mean((np.log(n) - np.log(data.index)) ** 2))),
        "N_end": float(n[-1]),
        "mu_max_per_h": float(mu.max()),
    }
    return out


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "single-fixed"
    module = {"single-legacy": single_legacy, "single-fixed": single_fixed}[which]
    print(json.dumps(impact(module), indent=1))
