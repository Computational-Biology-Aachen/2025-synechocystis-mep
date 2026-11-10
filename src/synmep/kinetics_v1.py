"""Rate laws of v1 (started as a copy of `kinetics_legacy`; fixes per ADR 0002 §5).

modelbase rate functions are inlined below. Units: mM, min. Enzyme concentrations are
in mM too, so the legacy hidden `Et / 1000` (µM -> mM) conversions are gone.
"""


def mass_action_1(S1, k_fwd):
    return k_fwd * S1


def reversible_mass_action_1_1(S1, P1, k_fwd, k_bwd):
    return k_fwd * S1 - k_bwd * P1


def ping_pong_4(A, B, C, D, vmax, kmA, kmB, kmC, kmD):
    return (vmax * A * B * C * D) / (
        A * B * C * D
        + (kmA * B * C * D)
        + (kmB * A * C * D)
        + (kmC * A * B * D)
        + (kmD * A * B * C)
    )


def petc_supply(value, chl):
    """PETC rate [mmol mol(Chl)^-1 s^-1] x Chl [mM] -> [mM min^-1].

    The legacy `unit_conversion` documented µM min^-1 but returned mM min^-1 (B8).
    """
    return value * chl * 1e-3 * 60


def DXS_active(DXS_f, ThiDP, KdThiDP):
    return DXS_f * (ThiDP / (KdThiDP + ThiDP))


def DXS_free(
    DXS_total, IPP, DMAPP, KdIPP, KdDMAPP
):  # Amount of free, active DXS enzyme
    return DXS_total / (1 + IPP / KdIPP + DMAPP / KdDMAPP)


def v_dxs(
    Pyr, G3P, DXS, kcat, KmPyr, KmG, Kia, IPP, DMAPP, ThiDP, KdIPP, KdDMAPP, KdThiDP
):  # Rate equation of active DXS enzyme
    numerator = (
        kcat
        * (DXS_active(DXS_free(DXS, IPP, DMAPP, KdIPP, KdDMAPP), ThiDP, KdThiDP))
        * Pyr
        * G3P
    )  # Fraction of active DXS bound to ThiDP
    denominator = Kia * KmPyr + KmG * Pyr + KmPyr * G3P + Pyr * G3P
    return numerator / denominator


def DXS_na(DXS, KnaDXS, IPP, DMAPP, KdIPP, KdDMAPP, ThiDP, KdThiDP):
    """Monomerisation of non-active DXS: k_na (E_t - E_a) (paper Methods).

    B4 (ADR 0002): the legacy code used E_t - E_f with E_f already divided by 1000
    (µM - mM), so practically all DXS was monomerised regardless of IPP/DMAPP. E_a is
    the active enzyme (free of IPP/DMAPP and ThiDP-bound).
    """
    DXS_a = DXS_active(DXS_free(DXS, IPP, DMAPP, KdIPP, KdDMAPP), ThiDP, KdThiDP)
    return KnaDXS * (DXS - DXS_a)


def universal_hill_2_1(A, B, P, kcat, Et, KmA, KmB, KmP, Keq, h):
    a = A / KmA
    b = B / KmB
    p = P / KmP
    numerator = (
        kcat
        * Et
        * a
        * b
        * (1 - ((P / (A * B)) / Keq))
        * (a + p) ** (h - 1)
        * b ** (h - 1)
    )
    denominator = 1 + (a + p) ** h * (1 + (b**h))
    return numerator / denominator


def reversible_michaelis_menten_keq_kcat(S, P, kcat, Et, KmS, KmP, Keq):
    return ((kcat * Et) * (S - P / Keq)) / (KmS * (1 + P / KmP) + S)


def reversible_michaelis_menten_keq_kcat_comp1(S, P, kcat, Et, KmS, KmP, Keq, I, Ki):
    return ((kcat * Et) * (S - P / Keq)) / (KmS * (1 + P / KmP + I / Ki) + S)


def v_sqs(
    A, B, C, P, kcat, Et, KmA, KmB, KmC, KmP, Keq, h
):  # idea: universal hill bi uni, followed by binding of NADPH and product release. Essentially rate of PSPP synthesis*NADPH binding
    a = A / KmA
    b = B / KmB
    p = P / KmP
    numerator = (
        kcat
        * Et
        * a
        * b
        * (1 - ((P / (A * B)) / Keq))  # B5: was (A * B / P) / Keq (inverted)
        * (a + p) ** (h - 1)
        * b ** (h - 1)
    )
    denominator = (1 + (a + p) ** h) * (1 + (b**h))
    return numerator / denominator * (C / (C + KmC))


def reversible_ordered_2_comp1(A, B, P, A2, KcatF, Et, KmA, KmA2, KmB, KmP, KiA, Keq):
    numerator = KcatF * Et * (A * B - (P / Keq))
    denominator = (
        A * B
        + KmA * (1 + A2 / KmA2) * B
        + KmB * A
        + (KmP + P * (1 + A / KiA)) / Keq  # B7: was kcat*Et/(kcat*Et*Keq), 0/0 at Et=0
    )
    return numerator / denominator


def reversible_ordered_2_comp2(
    A, B, P, A2, A3, KcatF, Et, KmA, KmA2, KmA3, KmB, KmP, KiA, Keq
):
    numerator = KcatF * Et * (A * B - (P / Keq))
    denominator = (
        A * B
        + KmA * (1 + A2 / KmA2 + A3 / KmA3) * B
        + KmB * A
        + (KmP + P * (1 + A / KiA)) / Keq  # B7: was kcat*Et/(kcat*Et*Keq), 0/0 at Et=0
    )
    return numerator / denominator


def chl_synthesis_rate(
    chlorophyll_GGPP,
    NADPH,
    chlorophyll_phytol,
    GGPP,
    ChlP_kcat,
    ChlP,
    ChlPkm_chlorophyll_GGPP,
    ChlPKm_GGPP,
    ChlPKm_NADPH,
    ChlPKm_chlorophyll_phytol,
    ChlP_Ki_chlorophyll_GGPP,
    ChlPKeq,
    phytol,
    ChlG_kcat,
    ChlG,
    ChlGKm_phytol,
    ChlGKm_chlorophyll,
    ChlGKeq,
    ChlGKm_GGPP,
):
    V1 = reversible_ordered_2_comp1(
        chlorophyll_GGPP,
        NADPH,
        chlorophyll_phytol,
        GGPP,
        ChlP_kcat,
        ChlP,
        ChlPkm_chlorophyll_GGPP,
        ChlPKm_GGPP,
        ChlPKm_NADPH,
        ChlPKm_chlorophyll_phytol,
        ChlP_Ki_chlorophyll_GGPP,
        ChlPKeq,
    )
    # ["chlorophyll-GGPP","NADPH", "chlorophyll-phytol","GGPP", "ChlPKcat","ChlP", "ChlPchlorophyll-GGPPKm","ChlPGGPPKm","ChlPNADPHKm","ChlPchlorophyll-phytolKm", "ChlPchlorophyll-GGPPKi","ChlP_Keq"
    V2 = reversible_michaelis_menten_keq_kcat_comp1(
        phytol,
        chlorophyll_phytol,
        ChlG_kcat,
        ChlG,
        ChlGKm_phytol,
        ChlGKm_chlorophyll,
        ChlGKeq,
        GGPP,
        ChlGKm_GGPP,
    )
    return V1 + V2


def growth_function(chl_synthesis, cell_volume, slope) -> float:
    """Empirical growth law mu = slope x Chl synthesis per cell (paper Results).

    [min^-1] = [cells mmol^-1] x [mM min^-1] x [L cell^-1]. At steady state it pins
    Chl* = 1 / (slope x cell_volume) (ADR 0002 §2.2).
    """
    return slope * chl_synthesis * cell_volume


def fraction_PGA_derivates_G3P(PGA):
    return PGA * 0.027


def fraction_PGA_derivates_PYR(PGA):
    return PGA * 0.972


def IspH1(IspH):
    return IspH * 0.9375


def IspH2(IspH):
    return IspH * 0.0625


def DXS_synth_f(A, B, C, D, kfwd, KmA, KmB, KmC, KmD, DXS_nonactive, Inh_strength):
    # return universal_hill_4_1(S1,S2,S3,S4, P,kfwd, KmA,KmB,KmC,KmD, KmP, Keq, h) * (DXS_nonactive/(DXS_nonactive+Inh_strength))
    return ping_pong_4(A, B, C, D, kfwd, KmA, KmB, KmC, KmD) * (
        DXS_nonactive / (DXS_nonactive + Inh_strength)
    )
