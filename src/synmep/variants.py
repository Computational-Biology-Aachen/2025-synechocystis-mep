"""Registry of the model variants in the comparison (ADR 0006 §2.3).

`adapter(name)` returns the `synmep.evaluate.Adapter` of a variant; the suffix `+ca`
adds chromatic acclimation (`double`, `minimal`, `full`; ADR 0006 §2.6). The `double`
batch (and its pre-culture) runs on the growth layer with a per-strain metabolism
surrogate (trained on first use, ~2 min); its turbidostat and MCA on the metabolism
ODE.
"""

from __future__ import annotations

import json
from functools import cache
from pathlib import Path

from synmep.data import load_legacy_inputs
from synmep.evaluate import Adapter
from synmep.runners import NATIVE_RUNNER, SUNDIALS_RUNNER

VARIANTS = (
    "single-legacy",
    "single-fixed",
    "single-fixed-cap",
    "double-none",
    "double-none-cap",
    "double-coupled",
    "double-coupled-cap",
    "minimal-m0",
    "minimal-m0-cap",
    "minimal-m2",
    "minimal-m2-cap",
    "full-k3",
)


FITS = Path("docs/reports/fits")


def fit_parameters(name: str) -> dict[str, float]:
    """Parameters of the shared fit (ADR 0006 §2.8), if done: fits/<name>.json."""
    path = FITS / f"{name}.json"
    return json.loads(path.read_text())["fitted"] if path.exists() else {}


def _with_fit(m, name: str):
    names = m.get_parameter_names()
    m.update_parameters({k: v for k, v in fit_parameters(name).items() if k in names})
    return m


def _single(name: str) -> Adapter:
    capacity = name.endswith("-cap")
    if name == "single-legacy":
        from synmep.models import single_legacy as mod

        inp, vcell, runner = load_legacy_inputs(), 2e-15, NATIVE_RUNNER
    else:
        from synmep.data_single_fixed import load_inputs
        from synmep.models import single_fixed as mod

        inp, vcell, runner = load_inputs(), 4e-15, SUNDIALS_RUNNER

    def build(mode, *, light, sample_depth, light_channels=None):
        kw = {} if light_channels is None else {"light_channels": light_channels}
        if capacity:
            kw["capacity"] = True
        if kw and name == "single-legacy":
            msg = "single-legacy has the warm white spectrum only, no options"
            raise ValueError(msg)
        m = mod.get_model(inp, "dshc", mode, **kw)
        m.update_parameters({"light": light, "sample_depth_m": sample_depth})
        return _with_fit(m, name)

    return Adapter(
        name=name,
        build=build,
        n_name="Cell_density",
        n_unit=1e6,  # [cells nL^-1]
        n_is_variable_in_steady=True,
        cell_volume=vcell,
        keys={
            "mu": "Growth_rate",
            "sq": "SQ",
            "v_sq": "v_Sqs",
            "chl": "chlorophyll_phytol",
            "car": "carotenoids",
        },
        runner=runner,
    )


def _double(name: str) -> Adapter:
    from synmep.models import double, double_growth

    chromatic = name.endswith("+ca")
    name_fit = name.removesuffix("+ca")
    capacity = name_fit.endswith("-cap")
    variant = name_fit.removeprefix("double-").removesuffix("-cap")
    inp = load_legacy_inputs()
    fitted = fit_parameters(name_fit) or double.fitted_parameters(variant)  # type: ignore[arg-type]

    def build(mode, *, light, sample_depth, light_channels=None, folds=None):
        if mode == "turbidostat":
            m = double.get_model(
                inp, variant, light=light, sample_depth=sample_depth,
                light_channels=light_channels, capacity=capacity, chromatic=chromatic,
            )  # type: ignore[arg-type]  # fmt: skip
            names = m.get_parameter_names()
            m.update_parameters({k: v for k, v in fitted.items() if k in names})
            return m
        # batch and its pre-culture: growth layer with the metabolism surrogate
        m = double_growth.get_model(
            inp, mode, pigment_regulation=variant, light_channels=light_channels,
            capacity=capacity, folds=folds, chromatic=chromatic,
        )  # type: ignore[arg-type]  # fmt: skip
        names = m.get_parameter_names()
        m.update_parameters(
            {k: v for k, v in fitted.items() if k in names}
            | {"light": light, "sample_depth_m": sample_depth}
        )
        return m

    return Adapter(
        name=name,
        build=build,
        n_name="N",
        n_unit=1.0,
        n_is_variable_in_steady=False,
        cell_volume=4e-15,
        keys={"mu": "mu", "sq": "SQ", "v_sq": "v_sq", "chl": "Chl", "car": "Car"},
        steady_mode="turbidostat",
        preculture_n_is_variable=True,
        folds_in_build=("steady", "dynamic"),
    )


def _minimal(name: str) -> Adapter:
    from synmep.models import minimal

    inp = load_legacy_inputs()
    base, chromatic = name.split("+ca")[0], name.endswith("+ca")
    kw = minimal.variant(base) | {"chromatic": chromatic}

    def build(mode, *, light, sample_depth, light_channels=None):
        m = minimal.get_model(
            inp, growth_mode=mode, light=light, sample_depth=sample_depth,
            light_channels=light_channels, **kw,
        )  # fmt: skip
        return _with_fit(m, base)

    return Adapter(
        name=name,
        build=build,
        n_name="N",
        n_unit=1.0,
        n_is_variable_in_steady=False,
        cell_volume=4e-15,
        keys={"mu": "mu", "sq": "SQ", "v_sq": "v_sq", "chl": "Chl", "car": "Car"},
    )


def _full(name: str) -> Adapter:
    from synmep.models import full

    inp = load_legacy_inputs()
    base, chromatic = name.removesuffix("+ca"), name.endswith("+ca")
    layers = int(base.removeprefix("full-k").split("-")[0])

    def build(mode, *, light, sample_depth, light_channels=None):
        m = full.get_model(
            inp, growth_mode=mode, light=light, sample_depth=sample_depth,
            light_channels=light_channels, layers=layers, chromatic=chromatic,
        )  # fmt: skip
        return _with_fit(m, base)

    return Adapter(
        name=name,
        build=build,
        n_name="N",
        n_unit=1.0,
        n_is_variable_in_steady=False,
        cell_volume=4e-15,
        keys={"mu": "mu", "sq": "SQ", "v_sq": "v_sq", "chl": "Chl", "car": "Car"},
        time_per_min=60.0,
        chunk=8,
    )


@cache
def adapter(name: str) -> Adapter:
    if name.startswith("single-"):
        return _single(name)
    if name.startswith("double-"):
        return _double(name)
    if name.startswith("minimal-"):
        return _minimal(name)
    if name.startswith("full-"):
        return _full(name)
    msg = f"unknown variant {name!r}"
    raise ValueError(msg)
