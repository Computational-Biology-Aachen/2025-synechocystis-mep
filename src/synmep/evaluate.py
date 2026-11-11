"""Family-independent evaluation: batch culture and turbidostat (ADR 0006 §2.1).

Every variant is wrapped in an `Adapter` that builds its model in steady mode (cell
density fixed, turbidostat / MCA / pre-culture) or dynamic mode (dN/dt = mu N, batch)
for given enzyme folds and light, and names its outputs. On top of that:

- `batch`: Fig 9 protocol. Pre-culture steady state at the density of data row 24
  under the batch light, diluted to the first data point, simulated for 96 h.
  Titer(t) = SQ(t) N(t) V_cell M_SQ [mg L^-1]; mean productivity
  P(t) = (titer(t) - titer(0)) / t; best harvest = argmax P; and titer after 96 h.
- `turbidostat`: steady states over cell densities; P = v_sq N V_cell M_SQ.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from mxlpy import Simulator

from synmep.data import load_legacy_inputs
from synmep.runners import NATIVE, NATIVE_RUNNER, SUNDIALS_RUNNER

if TYPE_CHECKING:
    from collections.abc import Callable

    from mxlpy import Model

M_SQ = 410.73  # [g mol^-1]
BATCH_LIGHT = {"light": 80.0, "sample_depth": 0.0135}  # MC-1000 batch (Fig 9)
TURBIDOSTAT_LIGHT = {"light": 100.0, "sample_depth": 0.01}  # Fig 12 / report 05


@dataclass
class Adapter:
    """How to build and read one variant.

    `build(mode, light, sample_depth, light_channels)` returns a model; `n_name` is the
    cell density (a parameter in steady mode, a variable in dynamic mode) in units of
    `n_unit` cells mL^-1; outputs `mu` [min^-1], `sq` [mM in the cell], `v_sq`
    [mM min^-1], `chl`, `car` [mM] are names of variables, derived values or fluxes.
    `enzymes` maps the overexpressable parameters to their names in the strain lists.
    """

    name: str
    build: Callable[..., Model]
    n_name: str
    n_unit: float
    n_is_variable_in_steady: bool
    cell_volume: float  # [L cell^-1]
    keys: dict[str, str]
    runner: object = field(default_factory=lambda: NATIVE_RUNNER)
    fixed: dict[str, float] = field(default_factory=dict)  # always applied
    # mode of the steady states (turbidostat, MCA) and of the batch pre-culture, and
    # whether N is a variable in the latter (default: as in steady mode); `double`
    # runs its turbidostat on the metabolism ODE, its pre-culture on the growth layer
    steady_mode: str = "steady"
    preculture_mode: str = "steady"
    preculture_n_is_variable: bool | None = None


def _row(res) -> pd.Series:
    row = pd.concat(
        [
            res.value.get_args(include_derived_parameters=True).iloc[-1],
            res.value.fluxes.iloc[-1],
        ]
    )
    return row[~row.index.duplicated()]


def steady_model(
    a: Adapter,
    folds: dict[str, float],
    density: float,
    *,
    light_channels=None,
    mode: str | None = None,
    **light,
) -> Model:
    mode = a.steady_mode if mode is None else mode
    m = a.build(mode, light_channels=light_channels, **light)
    p0 = m.get_parameter_values()
    missing = set(folds) - set(p0)
    if missing:
        msg = f"{a.name} ({mode}): can't overexpress {sorted(missing)}"
        raise ValueError(msg)
    m.update_parameters({k: p0[k] * f for k, f in folds.items()} | a.fixed)
    n_var = a.n_is_variable_in_steady
    if mode != a.steady_mode and a.preculture_n_is_variable is not None:
        n_var = a.preculture_n_is_variable
    if n_var:
        m.update_variable(a.n_name, density / a.n_unit)
    else:
        m.update_parameter(a.n_name, density / a.n_unit)
    return m


def steady_state(a: Adapter, m: Model, rtol: float = 1e-5) -> pd.Series | None:
    """Steady state via the runners; else a long time course whose end has a relative
    RHS below `rtol` (slowly relaxing states, e.g. dense pre-cultures)."""
    for runner in (a.runner, SUNDIALS_RUNNER):
        res = runner.steady_state(m)
        if hasattr(res.value, "fluxes"):
            return _row(res)
    res = Simulator(m, integrator=NATIVE).simulate(1e5).get_result()
    if not hasattr(res.value, "variables"):
        return None
    y = res.value.variables.iloc[-1][m.get_variable_names()]
    rhs = m.get_right_hand_side(y.to_dict()).to_numpy()
    if np.max(np.abs(rhs) / (np.abs(y.to_numpy()) + 1e-9)) > rtol:
        return None
    return _row(res)


def batch(
    a: Adapter,
    folds: dict[str, float] | None = None,
    *,
    hours: float = 96.0,
    n_points: int = 193,
    light_channels: tuple[str, ...] | None = None,
    weights: dict[str, float] | None = None,
    **light,
) -> dict | None:
    """Batch culture (Fig 9 protocol) of one strain; None if it fails."""
    folds = folds or {}
    light = BATCH_LIGHT | light
    data = _growth_data()
    pars = dict(weights or {})
    ms = steady_model(
        a,
        folds,
        float(data.index[24]),
        light_channels=light_channels,
        mode=a.preculture_mode,
        **light,
    )
    ms.update_parameters(pars)
    pre = steady_state(a, ms)
    if pre is None:
        return None
    md = a.build("dynamic", light_channels=light_channels, **light)
    p0 = md.get_parameter_values()
    missing = set(folds) - set(p0)
    if missing:
        msg = f"{a.name}: no batch mode for overexpressing {sorted(missing)}"
        raise ValueError(msg)
    md.update_parameters({k: p0[k] * f for k, f in folds.items()} | a.fixed | pars)
    y0 = {k: float(pre[k]) for k in md.get_variable_names() if k in pre.index}
    y0[a.n_name] = float(data.index[0]) / a.n_unit
    t = np.linspace(0, hours * 60, n_points)
    tc = Simulator(md, y0=y0, integrator=NATIVE).simulate_time_course(t).get_result()
    if not hasattr(tc.value, "variables"):
        return None
    args = tc.value.get_args(include_derived_parameters=True)
    fl = tc.value.fluxes
    both = pd.concat([args, fl.loc[:, ~fl.columns.isin(args.columns)]], axis=1)
    n = both[a.n_name].to_numpy() * a.n_unit  # [cells mL^-1]
    sq = both[a.keys["sq"]].to_numpy()  # [mM in the cell]
    titer = sq * 1e-3 * n * 1e3 * a.cell_volume * M_SQ * 1e3  # [mg L^-1]
    th = t / 60
    with np.errstate(divide="ignore", invalid="ignore"):
        prod = np.where(th > 0, (titer - titer[0]) / th, np.nan)  # [mg L^-1 h^-1]
    best = int(np.nanargmax(prod))
    return {
        "t_h": th.tolist(),
        "N": n.tolist(),
        "mu_per_h": (both[a.keys["mu"]].to_numpy() * 60).tolist(),
        "sq_mM": sq.tolist(),
        "chl_mM": both[a.keys["chl"]].tolist(),
        "titer_mg_L": titer.tolist(),
        "productivity_mg_L_h": prod.tolist(),
        "best_harvest_h": float(th[best]),
        "best_productivity_mg_L_h": float(prod[best]),
        "titer_end_mg_L": float(titer[-1] - titer[0]),
    }


def turbidostat(
    a: Adapter,
    folds: dict[str, float] | None = None,
    densities: np.ndarray | None = None,
    *,
    light_channels: tuple[str, ...] | None = None,
    weights: dict[str, float] | None = None,
    **light,
) -> pd.DataFrame:
    """Steady states over cell densities (continuation from the previous one)."""
    folds = folds or {}
    light = TURBIDOSTAT_LIGHT | light
    densities = np.geomspace(3e6, 3e10, 17) if densities is None else densities
    rows, y0 = [], None
    for n in densities:
        m = steady_model(a, folds, float(n), light_channels=light_channels, **light)
        m.update_parameters(dict(weights or {}))
        if y0 is not None:
            m.update_variables({k: v for k, v in y0.items() if k != a.n_name})
        row = steady_state(a, m)
        out = {"density": float(n)}
        if row is not None:
            v_sq = float(row[a.keys["v_sq"]])
            out |= {
                "mu_per_h": float(row[a.keys["mu"]]) * 60,
                "sq_mM": float(row[a.keys["sq"]]),
                "chl_mM": float(row[a.keys["chl"]]),
                "productivity_mg_L_h": v_sq * 60 * n * 1e3 * a.cell_volume * M_SQ,
            }
            y0 = {k: float(row[k]) for k in m.get_variable_names()}
        rows.append(out)
    return pd.DataFrame(rows)


@cache
def _growth_data() -> pd.Series:
    return load_legacy_inputs().exp_growth_rate_by_density
