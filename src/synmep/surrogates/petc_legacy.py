"""Legacy PETC surrogate and self-shading light model as explicit sympy expressions.

Faithful re-expression of the legacy `get_influx_rate_estimations` pipeline
(pandas + scipy.simpson + pickled sklearn polynomial). Spectral data and Simpson
weights are baked in as constants, because mxlpy's native integrator couldn't lower
model data (`add_data`) before 0.69 (MxlPy#201). Now it could; the constants are kept,
since single-legacy is regression-tested against the legacy outputs in this form.

Wavelength integrals use composite Simpson weights on the 1 nm grid (400-700 nm),
exactly what `scipy.integrate.simpson(y)` (dx = 1) computes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import sympy as sp
from scipy.integrate import simpson

from synmep.data import M_CHL, PIGMENTS, LegacyInputs

M_CAROTENOIDS = 536.8726  # [g mol^-1] used for the relative carotenoid content


@dataclass(frozen=True)
class Symbols:
    """Model names the expressions are written in."""

    chl: str  # [mM] chlorophyll
    carotenoids: str  # [mM]
    cell_density: str  # [cells nL^-1]
    light: str  # [µmol(photons) m^-2 s^-1]
    phycocyanin: str  # [mg L^-1]
    allophycocyanin: str  # [mg L^-1]
    ps_ratio: str
    beta_carotene_fraction: str
    sample_depth: str  # [m]
    cell_volume: str  # [L]


def simpson_weights(n: int) -> np.ndarray:
    return simpson(np.eye(n), axis=0)


def petc_expressions(inp: LegacyInputs, s: Symbols) -> dict[str, sp.Expr]:
    """Light model (`light_expressions`) plus the raw legacy PETC surrogate outputs.

    Adds `vATP_prod_raw`, `vNADPH_prod_raw`, `v3PGA_prod_raw`, `vFdred_prod_raw`
    [mmol mol(Chl)^-1 s^-1], unclipped. The legacy scalar path clips them at zero;
    single-legacy does that with `synmep.exprfn.clip_zero`.
    """
    out = light_expressions(inp, s)
    out |= polynomial_expressions(
        inp.surrogate_poly,
        {name: sp.Symbol(name) for name in inp.surrogate_poly["feature_names"]},
    )
    return out


def light_expressions(
    inp: LegacyInputs,
    s: Symbols,
    *,
    depth: str | None = None,
    prefix: str = "",
    step: int = 1,
) -> dict[str, sp.Expr]:
    """Expressions for the self-shading light model and the absorbed light per complex.

    By default the light is the depth average over the sample (legacy). With `depth`
    (name of a depth [m]) it is the local light `I0 e^{-k z}` at that depth, for
    depth-resolved PETC evaluation (ADR 0004). `prefix` namespaces all defined names,
    so several depths can coexist in one model.

    Returned in dependency order; every entry may reference the names of earlier
    entries, so each becomes one (small) derived quantity. Ends with the absorbed
    light per complex `complex_abs_ps1/ps2/pbs` [µmol(photons) mmol(Chl)^-1 s^-1],
    `light_ocp` and the photoactive pigment contents `pigment_*` [mg mg(Chl)^-1].
    """
    S = {k: sp.Symbol(v) for k, v in s.__dict__.items()}
    # `step` [nm] subsamples the 1 nm spectra (Simpson needs an odd point count)
    sl = slice(None, None, step)
    spectrum = inp.light_spectrum.to_numpy()[sl]
    w = simpson_weights(len(spectrum)) * step
    a = {p: inp.abs_coef[p].to_numpy()[sl] for p in PIGMENTS}
    ocp = inp.ocp_absorption.to_numpy()[sl]
    out: dict[str, sp.Expr] = {}

    def define(name: str, expr: sp.Expr) -> sp.Symbol:
        out[prefix + name] = expr
        return sp.Symbol(prefix + name)

    mg_chl = S["chl"] * 1000 * M_CHL / 1000  # [mg L^-1]
    content = {
        "chla": sp.Integer(1),
        "beta_carotene": define(
            "petc_rel_carotenoids", S["carotenoids"] * M_CAROTENOIDS / mg_chl
        ),
        "phycocyanin": define("petc_rel_phycocyanin", S["phycocyanin"] / mg_chl),
        "allophycocyanin": define(
            "petc_rel_allophycocyanin", S["allophycocyanin"] / mg_chl
        ),
    }
    photo = dict(content)
    photo["beta_carotene"] = define(
        "pigment_beta_carotene",
        content["beta_carotene"] * S["beta_carotene_fraction"],
    )
    # [mg(Chl) m^-3] chlorophyll in the sample
    chl_sample = define(
        "petc_chl_sample",
        S["chl"]
        * 1000
        * S["cell_volume"]
        * S["cell_density"]
        * 1e6
        * 1e6
        * 1e-3
        * M_CHL,
    )
    d = S["sample_depth"]

    # light per wavelength: depth average, legacy form of I0 (1 - e^{-kd}) / (kd),
    # or the local light I0 e^{-kz} at depth z
    z = None if depth is None else sp.Symbol(depth)
    corrected = []
    for i in range(len(spectrum)):
        k = chl_sample * sp.Add(*(float(a[p][i]) * content[p] for p in PIGMENTS))
        i0 = float(spectrum[i]) * S["light"]
        if z is None:
            corrected.append(i0 * (sp.exp(-d * k) - 1) / (-k) / d)
        else:
            corrected.append(i0 * sp.exp(-k * z))

    def integrate(per_wl: np.ndarray) -> sp.Expr:
        return sp.Add(
            *(
                float(w[i] * per_wl[i]) * corrected[i]
                for i in range(len(corrected))
                if per_wl[i] != 0
            )
        )

    absorbed = {
        p: photo[p] * define(f"petc_int_{p}", integrate(a[p])) for p in PIGMENTS
    }

    # pigment association (legacy `get_pigment_association`, stoichiometric method)
    pc = inp.ps_comp
    r = S["ps_ratio"]
    n_chla1, n_chla2 = float(pc.loc["n_chla", "ps1"]), float(pc.loc["n_chla", "ps2"])
    n_bc1 = float(pc.loc["n_beta_carotene", "ps1"])
    n_bc2 = float(pc.loc["n_beta_carotene", "ps2"])
    ratio2 = float(pc.loc["ratio", "ps2"])
    chla_ps1 = r * n_chla1 / (r * n_chla1 + ratio2 * n_chla2)
    chla_ps2 = ratio2 * n_chla2 / (r * n_chla1 + ratio2 * n_chla2)
    bc_ps1 = r * n_bc1 / (r * n_bc1 + ratio2 * n_bc2)
    bc_stoich = (r * n_bc1 / (r + ratio2)) / n_chla1 + (
        ratio2 * n_bc2 / (r + ratio2)
    ) / n_chla2
    mm = inp.molar_masses
    bc_measured = (
        photo["beta_carotene"] * float(mm["chla"]) / float(mm["beta_carotene"])
    )
    bc_fac = define("petc_bc_fac", sp.Min(bc_stoich / bc_measured, 1))

    define(
        "complex_abs_ps1",
        (chla_ps1 * absorbed["chla"] + bc_ps1 * bc_fac * absorbed["beta_carotene"])
        * M_CHL,
    )
    define("complex_abs_ps2", chla_ps2 * absorbed["chla"] * M_CHL)
    define(
        "complex_abs_pbs",
        (absorbed["allophycocyanin"] + absorbed["phycocyanin"]) * M_CHL,
    )
    define("light_ocp", integrate(ocp))
    define("pigment_phycocyanin", photo["phycocyanin"])
    define("pigment_allophycocyanin", photo["allophycocyanin"])
    return out


def polynomial_expressions(
    poly: dict[str, np.ndarray],
    inputs: dict[str, sp.Expr],
) -> dict[str, sp.Expr]:
    """Explicit (unclipped) form of the legacy PolynomialFeatures → LinearRegression."""
    x = [inputs[str(n)] for n in poly["feature_names"]]
    monomials = [
        sp.Mul(*(xi ** int(p) for xi, p in zip(x, powers, strict=True) if p))
        for powers in poly["powers"]
    ]
    names = {"ATP": "vATP_prod", "NADPH": "vNADPH_prod", "3PGA": "v3PGA_prod"}
    names["Fd_red"] = "vFdred_prod"
    res = {}
    for j, out in enumerate(poly["output_names"]):
        expr = float(poly["intercept"][j]) + sp.Add(
            *(
                float(c) * mono
                for c, mono in zip(poly["coef"][j], monomials, strict=True)
                if c != 0
            )
        )
        res[f"{names[str(out)]}_raw"] = expr
    return res


def evaluate_light(
    inp: LegacyInputs,
    states: pd.DataFrame,
    parameters: dict[str, float],
    outputs: tuple[str, ...] = ("complex_abs_ps1", "complex_abs_ps2", "complex_abs_pbs",
                                "light_ocp"),
    step: int = 1,
) -> pd.DataFrame:  # fmt: skip
    """Vectorised evaluation of `light_expressions` for many states at once.

    `states` and `parameters` together must provide every field of `Symbols` (the
    column/key names are the field names: chl, carotenoids, cell_density [cells nL^-1],
    light, ...). With a `z` column, the local light at depth z [m] is used instead of
    the depth average.
    """
    exprs = light_expressions(
        inp,
        Symbols(**{f: f for f in Symbols.__dataclass_fields__}),
        depth="z" if "z" in states else None,
        step=step,
    )
    env: dict[str, np.ndarray | float] = {c: states[c].to_numpy(float) for c in states}
    env |= {k: v for k, v in parameters.items() if k not in env}
    for name, expr in exprs.items():
        syms = sorted(expr.free_symbols, key=str)
        fn = sp.lambdify(syms, expr, modules="numpy")
        env[name] = fn(*(env[str(s)] for s in syms)) * np.ones(len(states))
    return pd.DataFrame({o: env[o] for o in outputs}, index=states.index)
