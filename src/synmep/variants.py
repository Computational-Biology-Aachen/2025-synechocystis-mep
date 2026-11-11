"""Registry of the model variants in the comparison (ADR 0006 §2.3).

`adapter(name)` returns the `synmep.evaluate.Adapter` of a variant. The `double`
batch runs on the growth layer with the metabolism surrogate, which only takes E_CrtE;
other enzyme folds need the retrained surrogate (ADR 0006 stage 3).
"""

from __future__ import annotations

from functools import cache

from synmep.data import load_legacy_inputs
from synmep.evaluate import Adapter
from synmep.runners import NATIVE_RUNNER, SUNDIALS_RUNNER

VARIANTS = (
    "single-legacy",
    "single-fixed",
    "double-none",
    "double-coupled",
    "minimal-m0",
    "minimal-m0-cap",
)


def _single(name: str) -> Adapter:
    if name == "single-legacy":
        from synmep.models import single_legacy as mod

        inp, vcell, runner = load_legacy_inputs(), 2e-15, NATIVE_RUNNER
    else:
        from synmep.data_single_fixed import load_inputs
        from synmep.models import single_fixed as mod

        inp, vcell, runner = load_inputs(), 4e-15, SUNDIALS_RUNNER

    def build(mode, *, light, sample_depth, light_channels=None):
        kw = {} if light_channels is None else {"light_channels": light_channels}
        if kw and name == "single-legacy":
            msg = "single-legacy has the warm white spectrum only"
            raise ValueError(msg)
        m = mod.get_model(inp, "dshc", mode, **kw)
        m.update_parameters({"light": light, "sample_depth_m": sample_depth})
        return m

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


def _double(variant: str) -> Adapter:
    from synmep.models import double, double_growth

    inp = load_legacy_inputs()
    fitted = double.fitted_parameters(variant)  # type: ignore[arg-type]

    def build(mode, *, light, sample_depth, light_channels=None):
        if mode == "turbidostat":
            return double.get_model(
                inp, variant, light=light, sample_depth=sample_depth,
                light_channels=light_channels,
            )  # type: ignore[arg-type]  # fmt: skip
        # batch and its pre-culture: growth layer with the metabolism surrogate
        m = double_growth.get_model(
            inp, mode, pigment_regulation=variant, light_channels=light_channels
        )  # type: ignore[arg-type]
        names = m.get_parameter_names()
        m.update_parameters(
            {k: v for k, v in fitted.items() if k in names}
            | {"light": light, "sample_depth_m": sample_depth}
        )
        return m

    return Adapter(
        name=f"double-{variant}",
        build=build,
        n_name="N",
        n_unit=1.0,
        n_is_variable_in_steady=False,
        cell_volume=4e-15,
        keys={"mu": "mu", "sq": "SQ", "v_sq": "v_sq", "chl": "Chl", "car": "Car"},
        steady_mode="turbidostat",
        preculture_n_is_variable=True,
    )


def _minimal(name: str) -> Adapter:
    from synmep.models import minimal

    inp = load_legacy_inputs()
    kw = minimal.variant(name)

    def build(mode, *, light, sample_depth, light_channels=None):
        return minimal.get_model(
            inp, growth_mode=mode, light=light, sample_depth=sample_depth,
            light_channels=light_channels, **kw,
        )  # fmt: skip

    return Adapter(
        name=name,
        build=build,
        n_name="N",
        n_unit=1.0,
        n_is_variable_in_steady=False,
        cell_volume=4e-15,
        keys={"mu": "mu", "sq": "SQ", "v_sq": "v_sq", "chl": "Chl", "car": "Car"},
    )


@cache
def adapter(name: str) -> Adapter:
    if name.startswith("single-"):
        return _single(name)
    if name.startswith("double-"):
        return _double(name.removeprefix("double-"))
    if name.startswith("minimal-"):
        return _minimal(name)
    msg = f"unknown variant {name!r}"
    raise ValueError(msg)
