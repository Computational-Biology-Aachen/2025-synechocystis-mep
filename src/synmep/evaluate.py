"""Family-independent evaluation: batch culture and turbidostat (ADR 0006 §2.1).

Every variant is wrapped in an `Adapter` that builds its model in steady mode (cell
density fixed, turbidostat / MCA / pre-culture) or dynamic mode (dN/dt = mu N, batch)
for given enzyme folds and light, and names its outputs. On top of that:

- `batch`: Fig 9 protocol. Pre-culture steady state at the density of data row 24
  under the batch light, diluted to the first data point, simulated for 96 h.
  Titer(t) = SQ(t) N(t) V_cell M_SQ [mg L^-1]; mean productivity
  P(t) = (titer(t) - titer(0)) / t; best harvest = argmax P; and titer after 96 h.
- `turbidostat`: steady states over cell densities; P = v_sq N V_cell M_SQ.

Plausibility (report 10): a pre-culture must grow (mu >= `MU_MIN`; the steady-state
solver can land on the trivial non-growing state, from which a batch never starts), and
a batch whose squalene content exceeds `SQ_MAX` (a quarter of the dry weight, a mass
balance bound far above any reported terpenoid content) counts as failed: such runs are
surrogate artefacts, not predictions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cache
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
from mxlpy import Simulator
from mxlpy.integrators import Sundials

from synmep.data import M_CHL, load_legacy_inputs
from synmep.runners import NATIVE, NATIVE_RUNNER, SUNDIALS_RUNNER

if TYPE_CHECKING:
    from collections.abc import Callable

    from mxlpy import Model

M_SQ = 410.73  # [g mol^-1]
BATCH_LIGHT = {"light": 80.0, "sample_depth": 0.0135}  # MC-1000 batch (Fig 9)
TURBIDOSTAT_LIGHT = {"light": 100.0, "sample_depth": 0.01}  # Fig 12 / report 05
MU_MIN = 1e-3  # [h^-1] growing pre-culture (Fig 9 data: 0.008-0.060 h^-1)
DW_DENSITY = 463.7  # [g L^-1 cell] `synmep.models.minimal.dry_weight_density()`
SQ_MAX = 0.25 * DW_DENSITY / M_SQ * 1e3  # [mM] 25 % of the dry weight as squalene


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
    # modes whose `build` takes the enzyme folds itself (`double`'s per-strain
    # metabolism surrogates) instead of having them applied to parameters
    folds_in_build: tuple[str, ...] = ()
    time_per_min: float = 1.0  # model time units per minute (`full`: seconds -> 60)
    chunk: int | None = None  # integrate batches in chunks of this many output steps


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
    pars: dict[str, float] | None = None,
    **light,
) -> Model:
    mode = a.steady_mode if mode is None else mode
    m = _build(a, mode, folds, light_channels, light, pars)
    n_var = a.n_is_variable_in_steady
    if mode != a.steady_mode and a.preculture_n_is_variable is not None:
        n_var = a.preculture_n_is_variable
    if n_var:
        m.update_variable(a.n_name, density / a.n_unit)
    else:
        m.update_parameter(a.n_name, density / a.n_unit)
    return m


def _build(
    a: Adapter,
    mode: str,
    folds: dict[str, float],
    light_channels,
    light,
    pars: dict[str, float] | None = None,
) -> Model:
    """Model for `mode` with enzyme folds, `a.fixed` and parameter values `pars`
    (absolute values, applied last; names missing in this mode's model are skipped)."""
    if mode in a.folds_in_build:
        m = a.build(mode, light_channels=light_channels, folds=folds, **light)
        m.update_parameters(a.fixed)
        _set(m, pars)
        return m
    m = a.build(mode, light_channels=light_channels, **light)
    p0 = m.get_parameter_values()
    missing = set(folds) - set(p0)
    if missing:
        msg = f"{a.name} ({mode}): can't overexpress {sorted(missing)}"
        raise ValueError(msg)
    m.update_parameters({k: p0[k] * f for k, f in folds.items()} | a.fixed)
    _set(m, pars)
    return m


def _set(m: Model, pars: dict[str, float] | None) -> None:
    if pars:
        names = set(m.get_parameter_names())
        m.update_parameters({k: v for k, v in pars.items() if k in names})


def steady_state(
    a: Adapter, m: Model, rtol: float = 1e-5, *, growing: bool = False
) -> pd.Series | None:
    """Steady state via the runners; else a long time course whose end has a relative
    RHS below `rtol` (slowly relaxing states, e.g. dense pre-cultures). With `growing`,
    states with mu < `MU_MIN` are skipped (the time course starts from the model's
    initial values, i.e. growing cells)."""

    def ok(row: pd.Series) -> bool:
        return not growing or float(row[a.keys["mu"]]) * 60 >= MU_MIN

    for runner in (a.runner, SUNDIALS_RUNNER):
        res = runner.steady_state(m)
        if hasattr(res.value, "fluxes") and ok(row := _row(res)):
            return row
    res = Simulator(m, integrator=NATIVE).simulate(1e5 * a.time_per_min).get_result()
    if not hasattr(res.value, "variables"):
        return None
    y = res.value.variables.iloc[-1][m.get_variable_names()]
    rhs = m.get_right_hand_side(y.to_dict()).to_numpy()
    if np.max(np.abs(rhs) / (np.abs(y.to_numpy()) + 1e-9)) > rtol:
        return None
    row = _row(res)
    return row if ok(row) else None


def _preculture(
    a: Adapter,
    folds: dict[str, float],
    density: float,
    light_channels: tuple[str, ...] | None,
    pars: dict[str, float],
    light: dict[str, float],
    steps: int = 4,
) -> pd.Series | None:
    """Growing pre-culture steady state. Some variants are bistable (a growing and a
    non-growing branch); if the solver misses the growing one, follow the light from
    warm white to the target in `steps` (cells kept under warm white, then shifted)."""

    def solve(p: dict[str, float], y0: dict[str, float] | None) -> pd.Series | None:
        m = steady_model(a, folds, density, light_channels=light_channels,
                         mode=a.preculture_mode, pars=p, **light)  # fmt: skip
        if y0:
            m.update_variables({k: v for k, v in y0.items() if k != a.n_name})
        return steady_state(a, m, growing=True)

    pre = solve(pars, None)
    w = {k: v for k, v in pars.items() if k.startswith("w_")}
    if pre is not None or not light_channels or not w:
        return pre
    from synmep.light import DEFAULT, weights

    w0 = weights(DEFAULT, light_channels)
    names = list(steady_model(a, folds, density, light_channels=light_channels,
                              mode=a.preculture_mode, pars=pars, **light)
                 .get_variable_names())  # fmt: skip
    y0 = None
    for x in np.linspace(0, 1, steps + 1):
        p = pars | {k: (1 - x) * w0[k] + x * w[k] for k in w}
        row = solve(p, y0)
        if row is None:
            return None
        y0 = {k: float(row[k]) for k in names if k in row.index}
    return row


def batch(
    a: Adapter,
    folds: dict[str, float] | None = None,
    *,
    hours: float = 96.0,
    n_points: int = 193,
    light_channels: tuple[str, ...] | None = None,
    weights: dict[str, float] | None = None,
    pars: dict[str, float] | None = None,
    n0: float | None = None,
    **light,
) -> dict | None:
    """Batch culture (Fig 9 protocol) of one strain; None if it fails. `pars` sets
    parameter values, `n0` the inoculum [cells mL^-1] (default: Fig 9's)."""
    folds = folds or {}
    light = BATCH_LIGHT | light
    data = _growth_data()
    pars = dict(pars or {}) | dict(weights or {})
    pre = _preculture(a, folds, float(data.index[24]), light_channels, pars, light)
    if pre is None:
        return None
    md = _build(a, "dynamic", folds, light_channels, light, pars)
    y0 = {k: float(pre[k]) for k in md.get_variable_names() if k in pre.index}
    y0[a.n_name] = (float(data.index[0]) if n0 is None else n0) / a.n_unit
    t = np.linspace(0, hours * 60, n_points)  # [min]
    res = _time_course(md, y0, t * a.time_per_min, a.chunk)
    if res is None:
        return None
    args, fl = res
    both = pd.concat([args, fl.loc[:, ~fl.columns.isin(args.columns)]], axis=1)
    n = both[a.n_name].to_numpy() * a.n_unit  # [cells mL^-1]
    sq = both[a.keys["sq"]].to_numpy()  # [mM in the cell]
    if not np.all(np.isfinite(sq)) or np.max(sq) > SQ_MAX:
        return None
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
        "pc_per_chl": (pc_per_chl(a, both, md)).tolist(),
        "titer_mg_L": titer.tolist(),
        "productivity_mg_L_h": prod.tolist(),
        "best_harvest_h": float(th[best]),
        "best_productivity_mg_L_h": float(prod[best]),
        "titer_end_mg_L": float(titer[-1] - titer[0]),
    }


def _time_course(
    m: Model, y0: dict[str, float], t: np.ndarray, chunk: int | None = None
) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Args and fluxes at times t. If one integration over t fails (stiff transients,
    e.g. `full`), the time span is halved recursively from the last good state; single
    intervals get a Sundials retry."""
    names = m.get_variable_names()

    def run(y: dict[str, float], ts: np.ndarray):
        integrators = (NATIVE,) if len(ts) > 2 else (NATIVE, Sundials)
        for integ in integrators:
            tc = (
                Simulator(m, y0=y, integrator=integ)
                .simulate_time_course(ts - ts[0])
                .get_result()
            )
            if hasattr(tc.value, "variables"):
                a = tc.value.get_args(include_derived_parameters=True)
                f = tc.value.fluxes
                return [a.set_axis(ts)], [f.set_axis(ts)]
        if len(ts) <= 2:
            return None
        mid = len(ts) // 2
        first = run(y, ts[: mid + 1])
        if first is None:
            return None
        y_mid = first[0][-1].iloc[-1][names].to_dict()
        second = run(y_mid, ts[mid:])
        if second is None:
            return None
        a2 = [second[0][0].iloc[1:], *second[0][1:]]
        f2 = [second[1][0].iloc[1:], *second[1][1:]]
        return first[0] + a2, first[1] + f2

    ts = np.asarray(t, dtype=float)
    if chunk is None:
        out = run(dict(y0), ts)
        return None if out is None else (pd.concat(out[0]), pd.concat(out[1]))
    # fixed chunks of `chunk` output intervals: long single integrations of `full` hit
    # a very slow regime, short ones don't
    y, args, fluxes = dict(y0), [], []
    for start in range(0, len(ts) - 1, chunk):
        part = ts[start : start + chunk + 1]
        out = run(y, part)
        if out is None:
            return None
        a_part, f_part = pd.concat(out[0]), pd.concat(out[1])
        if args:
            a_part, f_part = a_part.iloc[1:], f_part.iloc[1:]
        args.append(a_part)
        fluxes.append(f_part)
        y = a_part.iloc[-1][names].to_dict()
    return pd.concat(args), pd.concat(fluxes)


def pc_per_chl(a: Adapter, rows: pd.DataFrame, m: Model) -> np.ndarray:
    """Phycocyanin : Chl mass ratio [mg mg^-1] (phycocyanin in mg L(cell)^-1; a
    derived value or, in `single`, a parameter)."""
    key = a.keys.get("pc", "phycocyanin")
    if key in rows:
        pc = rows[key].to_numpy()
    else:
        pc = np.full(len(rows), m.get_parameter_values()[key])
    return pc / (rows[a.keys["chl"]].to_numpy() * M_CHL)


def turbidostat(
    a: Adapter,
    folds: dict[str, float] | None = None,
    densities: np.ndarray | None = None,
    *,
    light_channels: tuple[str, ...] | None = None,
    weights: dict[str, float] | None = None,
    pars: dict[str, float] | None = None,
    **light,
) -> pd.DataFrame:
    """Steady states over cell densities (continuation from the previous one)."""
    folds = folds or {}
    light = TURBIDOSTAT_LIGHT | light
    densities = np.geomspace(3e6, 3e10, 17) if densities is None else densities
    rows, y0 = [], None
    for n in densities:
        m = steady_model(
            a, folds, float(n), light_channels=light_channels, pars=pars, **light
        )
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


FIG10_LIGHT = {"light": 80.0, "sample_depth": 0.02}  # as in the notebook


def fig10(
    a: Adapter, *, light_channels=None, pars: dict[str, float] | None = None, **light
) -> dict:
    """Fig 10: steady-state growth rate vs. cell density, against the batch data.

    Note: this reuses the Fig 9 measurements (every second point), read as quasi-steady
    growth at each density; it's a check of the steady-state assumption more than
    independent validation data.
    """
    data = _growth_data()
    light = FIG10_LIGHT | light
    dens = data.index.to_numpy()[::2]
    ref = data.to_numpy()[::2]
    tu = turbidostat(
        a, densities=dens, light_channels=light_channels, pars=pars, **light
    )
    mu = tu["mu_per_h"].to_numpy() if "mu_per_h" in tu else np.full(len(dens), np.nan)
    ok = np.isfinite(mu)
    return {
        "density": dens.tolist(),
        "mu_per_h": mu.tolist(),
        "mu_data": ref.tolist(),
        "n_converged": int(ok.sum()),
        "rmse_mu_per_h": float(np.sqrt(np.mean((mu[ok] - ref[ok]) ** 2))),
    }


CHROMATIC_REFERENCE = {"density": 2e8} | BATCH_LIGHT  # ~Fig 9 median density


def calibrate_chromatic(
    fixed: Adapter,
    chromatic: Adapter,
    pars: dict[str, float] | None = None,
    light_channels=None,
    weights: dict[str, float] | None = None,
) -> float:
    """Imbalance E2/E1 at the reference state (warm white, measured PS1:PS2) for
    `chromatic_reference` (ADR 0006 §2.6): steady state of the fixed-ratio variant,
    imbalance evaluated with the chromatic variant's expressions at that state."""
    from synmep import acclimation

    ref = dict(CHROMATIC_REFERENCE)
    n = ref.pop("density")
    # same model as the batch (and its pre-culture), so warm white keeps the ratio there
    m = steady_model(
        fixed, {}, n, light_channels=light_channels, pars=pars,
        mode=fixed.preculture_mode, **ref,
    )  # fmt: skip
    m.update_parameters(dict(weights or {}))
    row = steady_state(fixed, m)
    if row is None:
        msg = f"{fixed.name}: no reference steady state"
        raise RuntimeError(msg)
    mc = steady_model(
        chromatic, {}, n, light_channels=light_channels, pars=pars,
        mode=chromatic.preculture_mode, **ref,
    )  # fmt: skip
    mc.update_parameters(dict(weights or {}))
    y = {k: float(row[k]) for k in mc.get_variable_names() if k in row.index}
    y["ps_ratio"] = float(m.get_parameter_values()["ps_ratio"])
    return acclimation.imbalance(mc.get_args(y))
