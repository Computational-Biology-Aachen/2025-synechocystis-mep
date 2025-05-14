import numpy as np
from modelbase.ode import ratefunctions as rf


def DXS_active(DXS_f, ThiDP, KdThiDP):
    return DXS_f * (ThiDP / (KdThiDP + ThiDP))


def DXS_free(
    DXS_total, IPP, DMAPP, KdIPP, KdDMAPP
):  # Amount of free, active DXS enzyme
    return DXS_total / 1000 / (1 + IPP / KdIPP + DMAPP / KdDMAPP)


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


def DXS_na(
    DXS, KnaDXS, IPP, DMAPP, KdIPP, KdDMAPP
):  # Amount of non-active, non-bound DXS enzyme, targeted for degradation
    DXS_f = DXS_free(DXS, IPP, DMAPP, KdIPP, KdDMAPP)
    return KnaDXS * (DXS - DXS_f)


def DXS_deg(DXS_na, kturn):
    return DXS_na * kturn


def universal_hill_2_1(A, B, P, kcat, Et, KmA, KmB, KmP, Keq, h):
    a = A / KmA
    b = B / KmB
    p = P / KmP
    numerator = (
        kcat
        * Et
        / 1000
        * a
        * b
        * (1 - ((P / (A * B)) / Keq))
        * (a + p) ** (h - 1)
        * b ** (h - 1)
    )
    denominator = 1 + (a + p) ** h * (1 + (b**h))
    return numerator / denominator


def universal_hill_2_2(A, B, P, Q, kcat, Et, KmA, KmB, KmP, KmQ, Keq, h):
    a = A / KmA
    b = B / KmB
    p = P / KmP
    q = Q / KmQ
    numerator = (
        kcat
        * Et
        / 1000
        * a
        * b
        * (1 - ((P * Q) / (A * B) / Keq))
        * (a + p) ** (h - 1)
        * (b + q) ** (h - 1)
    )
    denominator = 1 + (a + p) ** (h) * (1 + (b + q) ** (h))
    return numerator / denominator


def universal_hill_4_1(A, B, C, D, P, vmax, KmA, KmB, KmC, KmD, KmP, Keq, h):
    a = A / KmA
    b = B / KmB
    c = C / KmC
    d = D / KmD
    p = P / KmP
    numerator = (
        vmax
        / 1000
        * a
        * b
        * c
        * d
        * (a) ** (h - 1)
        * (b) ** (h - 1)
        * c ** (h - 1)
        * d ** (h - 1)
    )
    denominator = (1 + (a + p) ** (h) * (1 + (b) ** (h))) * (1 + c**h) * 1 + (d**h)
    return numerator / denominator


def reversible_michaelis_menten_keq_kcat(S, P, kcat, Et, KmS, KmP, Keq):
    return ((kcat * Et / 1000) * (S - P / Keq)) / (KmS * (1 + P / KmP) + S)


def reversible_michaelis_menten_keq_kcat_comp1(S, P, kcat, Et, KmS, KmP, Keq, I, Ki):
    return ((kcat * Et / 1000) * (S - P / Keq)) / (KmS * (1 + P / KmP + I / Ki) + S)


def random_2_1_Keq(A, B, P, Et, kcat, KmA, KmB, KmP, KiA, Keq):
    return (
        kcat
        * Et
        / 1000
        * (A * B - (P / Keq))
        / (
            A * B
            + KmB * A
            + KmA * B
            + KiA * KmB
            + (kcat * Et / 1000 * Keq) * (KmP + P * (1 + A / KiA))
        )
    )


def v_sqs(
    A, B, C, P, kcat, Et, KmA, KmB, KmC, KmP, Keq, h
):  # idea: universal hill bi uni, followed by binding of NADPH and product release. Essentially rate of PSPP synthesis*NADPH binding
    a = A / KmA
    b = B / KmB
    p = P / KmP
    numerator = (
        kcat
        * Et
        / 1000
        * a
        * b
        * (1 - (((A * B) / P) / Keq))
        * (a + p) ** (h - 1)
        * b ** (h - 1)
    )
    denominator = (1 + (a + p) ** h) * (1 + (b**h))
    return numerator / denominator * (C / (C + KmC))


def reversible_ordered_2_comp1(A, B, P, A2, KcatF, Et, KmA, KmA2, KmB, KmP, KiA, Keq):
    numerator = KcatF * Et / 1000 * (A * B - (P / Keq))
    denominator = (
        A * B
        + KmA * (1 + A2 / KmA2) * B
        + KmB * A
        + (KcatF * Et / 1000 / (KcatF * Et / 1000 * Keq) * (KmP + P * (1 + A / KiA)))
    )
    return numerator / denominator


def reversible_ordered_2_comp2(
    A, B, P, A2, A3, KcatF, Et, KmA, KmA2, KmA3, KmB, KmP, KiA, Keq
):
    numerator = KcatF * Et / 1000 * (A * B - (P / Keq))
    denominator = (
        A * B
        + KmA * (1 + A2 / KmA2 + A3 / KmA3) * B
        + KmB * A
        + (KcatF * Et / 1000 / (KcatF * Et / 1000 * Keq) * (KmP + P * (1 + A / KiA)))
    )
    return numerator / denominator


def reversible_ordered_2(A, B, P, Et, Kcatf, Kcatr, KmA, KmB, KmP, KiA, Keq):
    numerator = Kcatf * Et / 1000 * (A * B - (P / Keq))
    denominator = (
        A * B
        + KmA * B
        + KmB * A
        + (Kcatf * Et / 1000 / (Kcatr * Et / 1000 * Keq) * (KmP + P * (1 + A / KiA)))
    )
    return numerator / denominator


def reversible_ordered_2_sqs(A, P, Kcat, Et, KmA, KmB, KmP, KiA, Keq):
    numerator = Kcat * Et / 1000 * (A * A - (P / Keq))
    denominator = (
        A * A
        + KmA * A
        + KmB * A
        + (Kcat * Et / 1000 / (Kcat * Et / 1000 * Keq) * (KmP + P * (1 + A / KiA)))
    )
    return numerator / denominator


def reversible_random_2_sqs(A, P, Kcat, Et, KmA, KmB, KmP, KiA, Keq):
    numerator = Kcat * Et / 1000 * (A * A - (P / Keq))
    denominator = (
        A * A
        + KmA * A
        + KmB * A
        + (Kcat * Et / 1000 / (Kcat * Et / 1000 * Keq) * (KmP + P * (1 + A / KiA)))
    )
    return numerator / denominator


def reversible_ordered_2_2(
    A, B, P, Q, VmaxF, VmaxR, KmA, KmB, KmP, KmQ, KiA, KiB, KiP, KiQ, Keq
):
    numerator = VmaxF * VmaxR * (A * B - (P * Q) / Keq)

    denominator = (
        KiA * KmB * VmaxR
        + KmB * VmaxR * A
        + KmA * VmaxR * B
        + VmaxR * A * B
        + (KmQ * VmaxF * Q) / Keq
        + (KmP * VmaxF * Q) / Keq
        + (VmaxF * P * Q) / Keq
        + (KmQ * VmaxF * A * P) / KiA * Keq
        + (KmA * VmaxR * B * Q) / KiQ
        + (VmaxR * A * B * P) / KiP
        + (VmaxF * B * P * Q) / KiB * Keq
    )
    return numerator / denominator


def fraction(S, S_tot):
    return S / S_tot


def lightreaction_arbitrary_ATP(Chl, light, cell_volume):
    return (
        np.log((light * 1 * (Chl * 893.509 / 1000 * cell_volume) / cell_volume)) * 10
    )  # Convert Intracellular conc of chlorophyll from µmol/L to mg/L and multiply the output of the PS model (µmol/mg Chl) by that value


def lightreaction_arbitrary_NADPH(Chl, light, cell_volume):
    return np.log(
        light * 0.005 * (Chl * 893.509 / 1000 * cell_volume) / cell_volume
    )  # Convert Intracellular conc of chlorophyll from µmol/L to mg/L and multiply the output of the PS model (µmol/mg Chl) by that value


def lightreaction_arbitrary_Fdred(Chl, light, cell_volume):
    return np.log(
        light * 0.001 * (Chl * 893.509 / 1000 * cell_volume) / cell_volume
    )  # Convert Intracellular conc of chlorophyll from µmol/L to mg/L and multiply the output of the PS model (µmol/mg Chl) by that value


def lightreaction_arbitrary_3PGA(Chl, light, cell_volume):
    return np.log(
        light * 100 * (Chl * 893.509 / 1000 * cell_volume) / cell_volume
    )  # Convert Intracellular conc of chlorophyll from µmol/L to mg/L and multiply the output of the PS model (µmol/mg Chl) by that value


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


def growth_function(chl_synthesis, cell_volume) -> float:
    return (39566524308958.01) * (chl_synthesis * cell_volume)


def phytol_recycling() -> float:
    return 0.0003137085008420509


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
    return rf.ping_pong_4(A, B, C, D, kfwd, KmA, KmB, KmC, KmD) * (
        DXS_nonactive / (DXS_nonactive + Inh_strength)
    )


# Define the fraction of free PBS as a derived variable
def PBS_free(PBS_PS1, PBS_PS2):
    return 1 - PBS_PS1 - PBS_PS2


def simplified_MM(Et, kcat, S):
    return (Et * kcat * S) / S


def MM(Et, kcat, S, Km):
    return (Et * kcat * S) / (S + Km)


def _irradiance_depth_int(I0, depth, absorption_coef, chlorophyll_sample):
    """Effective irradiance"""
    return (
        I0
        * 1
        / (-absorption_coef * chlorophyll_sample)
        * np.exp(-depth * absorption_coef * chlorophyll_sample)
    )
