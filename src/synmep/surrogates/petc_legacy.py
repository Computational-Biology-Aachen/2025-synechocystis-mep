"""Legacy PETC surrogate and self-shading light model as explicit sympy expressions.

Faithful re-expression of the legacy `get_influx_rate_estimations` pipeline
(pandas + scipy.simpson + pickled sklearn polynomial). Spectral data and Simpson
weights are baked in as constants, so the result can be lowered by mxlpy's native
integrator, which can't handle model data (`add_data`).

Wavelength integrals use composite Simpson weights on the 1 nm grid (400-700 nm),
exactly what `scipy.integrate.simpson(y)` (dx = 1) computes.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
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
    """Return expressions for the surrogate inputs and the four production rates.

    Keys: `complex_abs_ps1/ps2/pbs`, `light_ocp`, `pigment_*` (surrogate inputs) and
    `vATP_prod`, `vNADPH_prod`, `v3PGA_prod`, `vFdred_prod`
    [mmol mol(Chl)^-1 s^-1], clipped at zero as in the legacy scalar path.
    """
    S = {k: sp.Symbol(v) for k, v in s.__dict__.items()}
    w = simpson_weights(len(inp.light_spectrum))
    a = {p: inp.abs_coef[p].to_numpy() for p in PIGMENTS}
    spectrum = inp.light_spectrum.to_numpy()
    ocp = inp.ocp_absorption.to_numpy()

    chl_uM = S["chl"] * 1000
    mg_chl = chl_uM * M_CHL / 1000  # [mg L^-1]
    content = {
        "chla": sp.Integer(1),
        "beta_carotene": S["carotenoids"] * 1000 * M_CAROTENOIDS / 1000 / mg_chl,
        "phycocyanin": S["phycocyanin"] / mg_chl,
        "allophycocyanin": S["allophycocyanin"] / mg_chl,
    }
    photo = dict(content)
    photo["beta_carotene"] = content["beta_carotene"] * S["beta_carotene_fraction"]

    # [mg(Chl) m^-3] chlorophyll in the sample
    chl_sample = (
        chl_uM * S["cell_volume"] * S["cell_density"] * 1e6 * 1e6 * 1e-3 * M_CHL
    )
    d = S["sample_depth"]

    # depth-averaged light per wavelength, legacy form of I0 (1 - e^{-kd}) / (kd)
    corrected = []
    for i in range(len(spectrum)):
        k = chl_sample * sum(float(a[p][i]) * content[p] for p in PIGMENTS)
        i0 = float(spectrum[i]) * S["light"]
        corrected.append(i0 * (sp.exp(-d * k) - 1) / (-k) / d)

    def integrate(per_wl: np.ndarray) -> sp.Expr:
        return sp.Add(
            *(
                float(w[i] * per_wl[i]) * corrected[i]
                for i in range(len(corrected))
                if per_wl[i] != 0
            )
        )

    absorbed = {p: photo[p] * integrate(a[p]) for p in PIGMENTS}

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
    bc_fac = sp.Min(bc_stoich / bc_measured, 1)

    out: dict[str, sp.Expr] = {
        "complex_abs_ps1": (
            chla_ps1 * absorbed["chla"] + bc_ps1 * bc_fac * absorbed["beta_carotene"]
        )
        * M_CHL,
        "complex_abs_ps2": chla_ps2 * absorbed["chla"] * M_CHL,
        "complex_abs_pbs": (absorbed["allophycocyanin"] + absorbed["phycocyanin"])
        * M_CHL,
        "light_ocp": integrate(ocp),
        "pigment_phycocyanin": photo["phycocyanin"],
        "pigment_allophycocyanin": photo["allophycocyanin"],
        "pigment_beta_carotene": photo["beta_carotene"],
    }
    out |= polynomial_expressions(
        inp.surrogate_poly,
        {name: sp.Symbol(name) for name in inp.surrogate_poly["feature_names"]},
    )
    return out


def polynomial_expressions(
    poly: dict[str, np.ndarray],
    inputs: dict[str, sp.Expr],
) -> dict[str, sp.Expr]:
    """Explicit form of the legacy PolynomialFeatures → LinearRegression pipeline."""
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
        res[names[str(out)]] = sp.Max(expr, 0)
    return res
