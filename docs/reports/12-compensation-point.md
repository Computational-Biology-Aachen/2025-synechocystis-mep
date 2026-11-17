# Report 12: fittable compensation point, refit of `minimal` and `full-k3`

**Status:** 2026-11-16. Follow-up to report 11 and ADR 0006 §8. The respiration term of
the `minimal` PETC is now fitted, all six `minimal` variants and `full-k3` are refitted
(shared fit with the Zavřel colour target) and their goal runs redone. Fitting the
respiration term improves the white-light fit but does not restore growth under
blue-green LEDs; two further causes are identified (§3), and fixing them needs design
decisions (§5).

## 1. Changes

- **`petc_r`** (`minimal.get_model`): the respiration term r of the PETC closed form
  (m0: J_3PGA = J_max tanh(I/I_k) − r; m2: J_3PGA = v_CBB − r) is a model parameter,
  default from `petc_balance.json` (m0 0.260, m2 0.337 mmol mol(Chl)⁻¹ s⁻¹). It sets
  the light compensation point, which no fitted factor could move before.
- **Shared fit** (`scripts/shared_fit.py`, family `minimal`): new factor `petc_r`
  (0.02–2), and the pigment exponent g (`pigment_coupling`) bound widened from 1.5 to 4
  (factor 3 → 8 on g₀ = 0.5), since every fit of report 11 ended at 1.5.
- All refits start from the model defaults (earlier fits of report 11 discarded), so
  the bounds are the same for every variant.

## 2. Results

Shared fit (Fig 9 batch, Chl and PC:Chl of the proteomics batch, colour data; AIC
includes the 18 colour residuals, so it is comparable only within this report and with
the photoacclimation fits of report 11, not with report 10). Gate 2: Fig 9 RMSE μ ≤
0.0058 h⁻¹ and ln N ≤ 0.149. Gate 3: ≥ 95 % of the goal runs converge.

| variant | k | RMSE μ [h⁻¹] | RMSE ln N | AIC (rep. 11) | colour RMSE | LEDs with no growth | g | r / r₀ | gate 2 | gate 3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **minimal-m0** | 4 | **0.0055** | **0.048** | **135** (156) | 0.37 | 465, 495, 520 | 4.0 (bound) | 0.23 | ✓ | ✗ 0.946 |
| minimal-m2 | 4 | 0.0057 | 0.049 | 139 (153) | 0.38 | 465, 495, 520 | 4.0 (bound) | 0.29 | ✓ | ✗ 0.843 |
| minimal-m0-cap | 6 | 0.0075 | 0.077 | 181 (—) | 0.41 | 465–555 | 2.9 | 0.77 | ✗ | ✗ 0.860 |
| minimal-m2-cap | 6 | 0.0078 | 0.084 | 186 (—) | 0.41 | 465–555 | 3.0 | 0.88 | ✗ | ✗ 0.796 |
| minimal-m0-cap-abs | 6 | 0.0073 | 0.054 | 186 (—) | 0.40 | 465–555 | 2.8 | 0.55 | ✗ | ✗ 0.892 |
| minimal-m0-cap-flux | 6 | 0.0068 | 0.054 | 193 (262) | 0.45 | 465–555 | 4.0 (bound) | 0.02 (bound) | ✗ | ✗ 0.935 |
| full-k3 | 4 | 0.0066 | 0.086 | 162 (—) | 0.39 | 465, 495, 520 | 2.0 | — | ✗ | ✗ 0.774 |

Relative growth under the nine LEDs (633 nm = 1):

| | 435 | 465 | 495 | 520 | 555 | 596 | 663 | 687 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| data (Zavřel 2024) | 0.62 | 0.23 | 0.39 | 0.59 | 0.83 | 0.92 | 1.01 | 0.71 |
| minimal-m0 (report 11) | 0.58 | 0 | 0 | 0 | 0.20 | 0.80 | 0.85 | 0.47 |
| minimal-m0 | 0.71 | 0 | 0 | 0 | 0.07 | 0.61 | 0.93 | 0.79 |
| minimal-m2 | 0.78 | 0 | 0 | 0 | 0.04 | 0.59 | 0.97 | 0.85 |
| full-k3 | 0.33 | 0 | 0 | 0 | 0.14 | 0.64 | 0.81 | 0.31 |

- **The white-light fit improves.** `minimal-m0` and `minimal-m2` now pass gate 2
  (RMSE μ 0.0055 and 0.0057 h⁻¹, from 0.0075), with r lowered 3–4× and `pbs_scale`
  off its bound (0.73 instead of 3.5). Their Chl error at 96 h drops from +75 % to
  +57 %.
- **The capacity variants get worse than without capacity:** `capacity_K` sits at its
  lower bound in every one, and none passes gate 2. Under the colour target, capacity
  acclimation no longer pays for its two parameters (report 10's `minimal-m0-cap` won
  on white light alone).
- **The colour response does not improve.** 465–520 nm still give no growth in every
  variant, 555 nm gets worse, and the pigment exponent ends at the new bound (g = 4) in
  three variants.
- **`full-k3` fits for the first time with the energy outlets** (report 11 §4): RMSE μ
  0.0066 h⁻¹ (report 10: 0.0105), Chl within −29 to −16 %, but 3 LEDs without growth
  and 77 % convergence (turbidostats at high density, orange/red LED mixtures).
- **Goal 1 is unchanged:** every `minimal` variant picks Sqs alone as best batch strain
  (0.24–0.43 mg L⁻¹ h⁻¹ under warm white); the strains with MEP overexpression fail in
  the batch in every `minimal` variant (2 of the 5 strain batches, counted in gate 3).
  `full-k3` picks Sqs + CrtE + DXS (0.29).
- **Goal 2:** best LED setting is pure orange (`minimal-m0`, `-m2`, `-cap-abs`,
  `-cap-flux`) or orange:red 3:1 (`-cap` variants); blue-rich settings fail.

## 3. Why blue-green light still gives no growth

Tests on the report-11 `minimal-m0` fit, steady state at the colour-data conditions
(2.1e7 cells mL⁻¹, 1 cm, 25 µmol m⁻² s⁻¹):

1. **Respiration (fixed now) is not the main cause.** With r × 0.05, relative growth at
   555 nm rises from 0.20 to 0.50, but 465 and 495 nm stay at 0 and 520 nm reaches
   0.07. Under warm white the compensation point at low density is already ~5 µE.
2. **Light at 465–520 nm doesn't reach PS2.** There, almost all absorption is by
   carotenoids (Fuente2021 coefficients: β-carotene 0.054 vs Chl a 0.0017 m² mg⁻¹ at
   465 nm; Chl a falls to 0.4 % of its peak at 495 nm). In the light model
   (`petc_legacy.light_expressions`) carotenoid absorption feeds PS1 only, PS2 gets Chl
   only (its β-carotene isn't counted), and the m0 fit weights PS1 at ≈ 0 (w_PS2/w_PS1 ≈
   5e8). So the PS2 + PBS input at 465 nm is ~2 % of that at 633 nm, against 23 %
   growth in the data. The spectra themselves are fine: they reproduce the measured
   whole-cell spectra of the proteomics samples to 3 % RMS.
3. **Photoacclimation has a runaway.** In dim or carotenoid-absorbed light the rule
   raises CrtB and Chl synthesis together. Chl synthesis is NADPH-limited there, so
   GGPP goes to carotenoids: 98 mM at 520 nm, 197 mM at 465 nm, 10⁴ mM at 495 nm
   (633 nm: 3.3 mM). They screen Chl, absorption per Chl falls, the pigment factor
   rises further. When μ → 0, pigments (= synthesis / μ) have no bound, since nothing
   degrades them, and the culture stays on this non-growing branch. Tried, none rescues
   465/495 nm: the rule on Chl only (GGPP then overflows into carotenoids at 633 nm
   too), a Car:Chl feedback on CrtB, and a narrower clamp (0.5–2 or 0.25–4 instead of
   0.05–20). The data contradict the rule's direction: Chl per cell varies only
   0.9–1.8× across the LEDs and is lower at 465 nm than at 633 nm.

## 4. Status of the selection

`minimal-m0` is now the best fit (AIC 135, gate 2 ✓) and misses gate 3 by one run
(88/93; the two MEP-overexpressing batches and three blue-rich LED settings). Its
goal-1 answer (Sqs strain, 0.24 mg L⁻¹ h⁻¹, 23 mg L⁻¹ in 96 h) agrees with all
variants. Its goal-2 answer (orange) can't be trusted while 465–520 nm give no growth.

## 5. Decisions needed

1. **Light routing to PS2** (§3.2): count PS2's β-carotene (with a transfer
   efficiency), and/or let PS1 absorption contribute in m0 (refit the closed form with
   w_PS1 bounded away from 0, or use m2's PS1/PS2 balance). Changes the PETC input,
   so the PETC fits (report 08) are redone.
2. **Bound pigment synthesis** (§3.3): pigment turnover (dX/dt = v − (μ + k_d) X),
   a narrower photoacclimation range matched to the 0.9–1.8× of the data, or a
   different signal (redox or excitation balance instead of total absorption). This
   departs from ADR 0006 §8.
3. **Capacity acclimation:** keep testing the `-cap` variants, or retire them now that
   they lose to `minimal-m0` under the shared fit with the colour target.
