# Report 14: pigment synthesis in proportion to growth, PS2 signal, refit of all variants

**Status:** 2026-11-17. Decisions on report 13 §4:
- (a) pigments are made in proportion to growth;
- (b) turnover with Chl priority on NADPH is tried;
- (c) 465–520 nm are left out of the colour target;
- all three photoacclimation signals are tried.

The best combination, growth-proportional synthesis with the PS2 signal, is now the
default in `minimal` and `full`. All variants are refitted with it and their goal runs
redone. The capacity variants come back as the best fits; 555 nm now gives no growth
in every `minimal` variant.

## 1. Changes

- **Pigment synthesis modes** (`minimal.pigment_synthesis`, variant suffix
  `~enzyme` / `~growth` / `~turnover`):
  - **`growth` (new default):** Chl and carotenoid synthesis are scaled by μ / μ_ref.
    The photoacclimation factor sets the target content, and pigments can no longer
    pile up when growth stops. Both rates are near zero-order in GGPP
    (K = GGPP_ref / 10). Before this fix, the old K_m(Chl, GGPP) = 0.06 mM, far above
    GGPP_ref = 0.005 mM, let a GGPP pile-up raise Chl synthesis (Chl reached 180 mM in
    the first test). V_max is set so that the rates at the reference state are
    unchanged.
  - **`turnover`:** dX/dt = v − (μ + k_d) X for Chl and carotenoids, with
    k_d = 1 / (7 d) (Chl lifetime above one week: Vavilin & Vermaas 2007, Yao 2012).
    V_max is raised by (μ_ref + k_d) / μ_ref. K_m(Chl, NADPH) is lowered 10× to give
    Chl synthesis priority on NADPH. Implemented in `minimal` only; `full` raises
    `NotImplementedError`.
  - **`enzyme`:** the old rule (report 11), kept for comparison.
- **Photoacclimation signals** (suffix `~abs` / `~ps2` / `~flux`), each relative to its
  value at the Fig 9 median state:
  - `abs`: total absorption per Chl (old default);
  - `ps2` (new default): PS2 excitation, Chl absorption of PS2 plus the PBS share that
    reaches PS2;
  - `flux`: electron flow J_3PGA per Chl.

  `full` uses its own references for `abs` and `ps2`.
- **Colour target:** 465, 495 and 520 nm are left out (`COLOUR_EXCLUDED`). There, almost
  all absorption is by carotenoids, and how much of it reaches PS2 is outside the model
  (report 13 §1). The fits report `colour_rmse_growth` (LEDs used) and
  `colour_rmse_growth_all` (all LEDs).
- `PIGMENT_CLAMP` 0.5–1.5 and the fitted `petc_r` stay as in reports 12–13.

## 2. Screen on `minimal-m0`

All runs start from the defaults, with 465–520 nm left out:

| synthesis | signal | cost | AIC | RMSE μ [h⁻¹] | 555 nm growth (data 0.83) | g |
| --- | --- | --- | --- | --- | --- | --- |
| **growth** | **ps2** | **130** | **102** | **0.0056** | 0 | 0.26 |
| growth | abs | 135 | 106 | 0.0057 | 0 | 0.01 (bound) |
| growth | flux | 135 | 106 | 0.0057 | 0 | 0.01 (bound) |
| turnover | abs / ps2 / flux | 173–175 | 135–136 | 0.0061–0.0062 | 0.08–0.09 | 1.0–3.2 |
| enzyme | abs | 203 | 153 | 0.0060 | 0.12 | 4.0 (bound) |

- **Growth-proportional synthesis fits the white-light data clearly best.**
- **Only the PS2 signal keeps photoacclimation switched on** (g = 0.26). With total
  absorption or electron flow, the fit switches it off (g at its lower bound).
- **Turnover with Chl priority keeps 555 nm alive, but barely** (0.08–0.09), and it
  fits the white-light data worse by ~35 cost units.
- **Under `growth`, 555 nm dies** because the fit needs a higher respiration term r
  (0.16–0.20 mmol mol(Chl)⁻¹ s⁻¹, against 0.03 with turnover) to match the Fig 9
  batch. The compensation point under green light then lies above the 25 µE of the
  colour experiment.

## 3. Refit of all variants (growth, ps2)

There are 14 colour residuals now instead of 18, so AIC is comparable only within this
report. Gates and colour RMSE are comparable with report 13 (in brackets).

| variant | k | RMSE μ [h⁻¹] (rep. 13) | RMSE ln N | AIC | colour RMSE, used / all (rep. 13) | LEDs with no growth | PC:Chl error | gate 2 | gate 3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| minimal-m0 | 4 | 0.0056 (0.0062) | 0.044 | 102 | 0.40 / 0.41 (0.37) | 465–555 | −72 to −81 % | ✓ | ✓ 0.968 |
| minimal-m2 | 4 | 0.0057 (0.0064) | 0.046 | 106 | 0.40 / 0.41 (0.39) | 465–555 | −74 to −82 % | ✓ | ✓ 0.959 |
| **minimal-m0-cap** | 6 | **0.0051** (0.0094) | 0.040 | **73** | 0.42 / 0.40 (0.40) | 555 | −8 to +34 % | ✓ | ✓ 1.000 |
| **minimal-m2-cap** | 6 | **0.0050** (0.0081) | 0.039 | **69** | 0.41 / 0.39 (0.40) | 555 | −10 to +32 % | ✓ | ✓ 1.000 |
| minimal-m0-cap-abs | 6 | 0.0057 (0.0086) | 0.045 | 112 | 0.42 / 0.42 (0.39) | 465–555 | −69 to −79 % | ✓ | ✓ 0.968 |
| minimal-m0-cap-flux | 6 | 0.0052 (0.0074) | 0.041 | 108 | 0.42 / 0.42 (0.38) | 465–555 | −76 to −84 % | ✓ | ✓ 0.952 |
| full-k3 | 4 | 0.0074 (0.0073) | 0.096 | 156 | 0.37 / 0.39 (0.39) | 465–520 | +17 to +59 % | ✗ | ✗ 0.778 |

Relative growth under the LEDs (633 nm = 1):

| | 435 | 465 | 495 | 520 | 555 | 596 | 663 | 687 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| data (Zavřel 2024) | 0.62 | 0.23 | 0.39 | 0.59 | 0.83 | 0.92 | 1.01 | 0.71 |
| minimal-m0 | 0.68 | 0 | 0 | 0 | 0 | 0.53 | 0.77 | 0.50 |
| minimal-m0-cap | 0.96 | 0.22 | 0.05 | 0.06 | 0 | 0.49 | 1.01 | 0.92 |
| minimal-m2-cap | 0.94 | 0.25 | 0.05 | 0.07 | 0 | 0.49 | 1.01 | 0.90 |
| full-k3 | 0.32 | 0 | 0 | 0 | 0.12 | 0.67 | 0.82 | 0.32 |

- **Every `minimal` variant now passes gates 2 and 3.**
- **`minimal-m0-cap` and `minimal-m2-cap` fit best** by ~30 AIC units:
  - RMSE μ 0.0050–0.0051;
  - Chl within −12 to +14 % at all four sample times;
  - PC:Chl fixed for the first time, with `pbs_scale` 5.9 instead of 1.2–1.7;
  - all 93 goal runs converge;
  - 465 nm grows (0.22–0.25), although it is no longer fitted.

  But two parameters end at their lower bounds: `petc_r` at 0.02 (respiration in the
  PETC closed form close to 0) and `capacity_K` at 0.01. With these parameters,
  capacity acclimation takes over what the respiration term did. 435 nm is too high
  (0.94–0.96 against 0.62), and 596 nm is too low (0.49 against 0.92).
- **The other variants all miss PC:Chl by about −75 %.** Their fits keep `pbs_scale` at
  1.2–1.7, so the PBS level and the colour response pull in opposite directions. Only
  the PS2-flux capacity rule (`-cap`) lets PBS rise without overdriving growth.
- **555 nm gives no growth in any `minimal` variant**, not even with `petc_r` at its
  bound. The cause is not tested yet; the r argument of §2 doesn't apply to the `-cap`
  variants.
- **`full-k3` stays behind:** gate 2 ✗ (RMSE μ 0.0074), gate 3 ✗ (0.778; all orange-rich
  LED settings fail), Chl −11 to −34 %, and `V_bio` at its upper bound 20.

## 4. Goal answers

| variant | goal 1: best batch strain | productivity [mg L⁻¹ h⁻¹] | titer at 96 h [mg L⁻¹] | goal 2: best LED setting | failed settings |
| --- | --- | --- | --- | --- | --- |
| minimal-m0 | CrtE + Sqs | 0.26 | 24.6 | orange (0.46) | 3 blue-rich |
| minimal-m2 | CrtE + Sqs | 0.25 | 24.1 | orange (0.45) | 3 blue-rich |
| minimal-m0-cap | Sqs + CrtE + MEP | 0.63 | 37.5 | orange (0.69) | none |
| minimal-m2-cap | Sqs + CrtE + MEP | 0.65 | 39.4 | orange (0.69) | none |
| minimal-m0-cap-abs | Sqs | 0.39 | 37.3 | orange (0.91) | 3 blue-rich |
| minimal-m0-cap-flux | Sqs | 0.32 | 30.6 | orange (0.75) | 2 blue-rich |
| full-k3 | Sqs | 0.18 | 17.1 | warm white (0.18) | 7 orange-rich |

- **Goal 1:** for the first time, a `minimal` variant gains from MEP overexpression.
  In the two best fits, Sqs + CrtE + MEP is 1.7× better than Sqs alone, and DXS on top
  helps less (0.50–0.52). `minimal-m0` and `-m2` pick CrtE + Sqs, the variants with
  other capacity signals Sqs alone.
- **Goal 2:** pure orange in every `minimal` variant, as in reports 12–13. The orange
  LED (596 nm) is underpredicted in the colour fit (0.49 against 0.92), so the
  productivity there is probably underestimated rather than overestimated. `full-k3`
  can't answer goal 2 while its orange settings fail.

## 5. Decisions needed

1. **Selection:** `minimal-m2-cap` (AIC 69) or `minimal-m0-cap` (AIC 73, the variant
   selected at checkpoint C). Both pass the gates and give the same goal answers.
   ΔAIC = 4 favours m2 slightly; m0 is cheaper to run.
2. **Bounds:** `petc_r` (0.02) and `capacity_K` (0.01) end at their lower bounds in the
   two best fits. Widen the bounds and refit, or accept (r ≈ 0 means the capacity rule
   sets the compensation point).
3. **555 nm:** diagnose why green light gives no growth in the `-cap` variants (steady
   state at 555 nm: pigments, NADPH, PS2 input), or leave it out of the colour target
   too.
4. **`full-k3`:** turnover isn't ported, `V_bio` is at its bound, and orange settings
   fail. Keep it as a reference only, or debug the orange failures.
