# Report 10: comparison of all variants and selection

**Status:** draft for checkpoint C (ADR 0006 §2.10), 2026-11-13. The selection in §5
needs the user's decisions in §6 before ADR 0006 is closed.

Numbers: `docs/reports/summary.json` (`scripts/summary.py`), per variant
`docs/reports/fits/<variant>.json` (shared fit, `scripts/shared_fit.py`) and
`docs/reports/goals/<variant>.json` (`scripts/goal_runs.py`). Figures:
`figures/summary-fits.png`, `figures/summary-spectra.png`.

## 1. Protocol (ADR 0006 §2.8–2.9, amendments §5–6)

- **Shared fit.** Every variant gets the same fit: Fig 9 growth (μ, ln N), WT Chl per
  cell and PC:Chl at 4/24/48/96 h, ≤ 5 multiplicative factors. Biomass ATP is
  50 mmol gDW⁻¹ in all families.
- **Held out.** The proteomics batch (WT cell densities) is held out; Fig 10 is a
  check.
- **Goal runs.**
  - MCA (paper protocol), then each strain in the Fig 9 batch (best-harvest mean
    productivity, titer after 96 h) and in the turbidostat (optimum over density).
  - The best strain then runs under the 6 lamps and 15 blue/orange/red LED mixtures at
    equal photon flux, without and with chromatic acclimation (`+ca`).
- **Gates.**
  1. Paired spectral ranking against pfennig2024 with median error ≤ 0.05 (report 08).
  2. Fig 9 RMSE μ ≤ 0.0058 h⁻¹ and ln N ≤ 0.149.
  3. ≥ 95 % of the goal runs converge.
- **Plausibility guards** (commit 18dd606), counted as failures in gate 3: the
  pre-culture must grow (μ ≥ 1e-3 h⁻¹), and squalene must stay ≤ 25 % of the dry
  weight.
- **Fallback solver budget.** The Sundials steady-state fallback gets a 60 s
  wall-clock budget (commit f28f2a5); before, it could stall a run for hours.

## 2. Gates

| variant | k | RMSE μ [h⁻¹] | RMSE ln N | AIC | validation ln N | Fig 10 μ | gate 1 | gate 2 | gate 3 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **minimal-m0-cap** | 5 | **0.0034** | **0.026** | **60** | 0.76 | 0.010 | ✓ 0.030 | ✓ | ✗ 0.946 |
| minimal-m0 | 3 | 0.0059 | 0.054 | 99 | 0.66 | 0.011 | ✓ 0.030 | ✗ μ +2 % | ✓ 0.952 |
| minimal-m2 | 3 | 0.0058 | 0.054 | 100 | 0.67 | 0.010 | ✓ 0.014 | ✗ μ +0.5 % | ✗ 0.918 |
| minimal-m2-cap | 5 | 0.0065 | 0.054 | 69 | **0.62** | 0.017 | ✓ 0.014 | ✗ | ✗ 0.893 |
| full-k3 | 4 | 0.0105 | 0.170 | 228 | 0.67 | 0.019 (8/24 conv.) | ✓ (is pfennig) | ✗ | ✗ 0.30 |
| double-none(-cap) | 2 (4) | 0.013 (0.008) | 0.24 (0.18) | 300 (256) | 0.80 (0.93) | 0.022 (0.041) | ✗ 0.10 | ✗ | ✗ 0.86 (0.95) |
| double-coupled(-cap) | 2 (4) | 0.015 (0.011) | 0.29 (0.23) | 328 (292) | 0.76 (0.88) | 0.022 (0.037) | ✗ 0.10 | ✗ | ✗ 0.88 (0.84) |
| single-fixed(-cap) | 3 (5) | 0.012 (0.009) | 0.53 (0.39) | 445 (385) | 0.99 (0.92) | 0.018 (0.019) | ✗ 0.10 | ✗ | not run (†) |

† `single-fixed`'s MCA hangs inside the native solver, where the time budget can't
interrupt it; `single-fixed-cap` needs > 6 h, because its MCA picks 10 enzymes. Both
fail gates 1 and 2, so their goal runs were not pursued. The partial
`single-fixed-cap` run gives batch productivities of 0.004–0.015 mg L⁻¹ h⁻¹, about 25×
below `minimal`.

**No variant passes all three gates.** Two near misses:

- `minimal-m0-cap` fails gate 3 by 5 of 93 runs. Three are turbidostats at the highest
  density (3e9 cells mL⁻¹, far past the optimum). Two are the blue-rich LED settings
  without chromatic acclimation, where the culture collapses (the non-growing branch,
  §4.2). With chromatic acclimation all its spectral runs converge.
- `minimal-m0` fails gate 2 by 2 % in μ.

The **proteomics-batch validation** barely separates variants: ln N errors are
0.62–0.99 for all of them, because that culture differs systematically from Fig 9. AIC
separates them clearly: `minimal-m0-cap` has AIC 60, the next best 69, and every
non-minimal variant ≥ 228.

## 3. Goal answers

### 3.1 Goal 1: strain and batch productivity (Fig 9 conditions, warm white)

| variant | MCA sequence | best strain (batch) | batch best P [mg L⁻¹ h⁻¹] | titer 96 h [mg L⁻¹] | WT batch P | turbidostat best P (strain) |
| --- | --- | --- | --- | --- | --- | --- |
| **minimal-m0-cap** | Sqs, CrtE, MEP, DXS | Sqs | **0.40** (harvest 96 h) | **38** | 0.076 | 3.24 (Sqs+CrtE+MEP, 1.5e9 mL⁻¹) |
| minimal-m0 | Sqs | Sqs | 0.34 | 33 | — | 1.95 |
| minimal-m2 | Sqs, CrtE | Sqs | 0.35 | 33 | — | 2.92 |
| minimal-m2-cap | Sqs, CrtE, MEP, DXS | Sqs | 0.65 | 62 | — | 1.81 |
| full-k3 | Sqs | WT (Sqs batch fails) | 0.028 (WT) | 2.7 | 0.028 | 0.19 |

- **Robust across `minimal`:**
  - Sqs overexpression is the step that matters in the batch, about 5× over WT in
    minimal-m0-cap.
  - Every variant's MCA picks Sqs first.
  - The batch is still producing at 96 h, so the best harvest is the end of the run.
- **Not robust:** the further MCA steps (CrtE, MEP, DXS), chosen at 1e7 cells mL⁻¹ and
  100 µE, raise the turbidostat optimum but *lower* the batch productivity (minimal-m0-cap:
  0.40 → 0.32 → 0.22 → 0.21). The paper's low-density MCA doesn't carry over to a
  shaded batch.
- **Spread over the minimal variants:** batch 0.34–0.65 mg L⁻¹ h⁻¹, titer 33–62 mg L⁻¹,
  turbidostat 1.8–3.2 mg L⁻¹ h⁻¹.

### 3.2 Goal 2: light composition (batch productivity of the best strain / warm white)

| light | m0 | m0-cap | m2 | m2-cap | m0 +ca | m0-cap +ca | m2 +ca | m2-cap +ca |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cool white LED | 0.54 | 0.95 | 0.51 | 0.72 | 0.56 | 0.96 | 0.53 | 0.70 |
| fluorescent | 0.96 | 1.02 | 0.96 | 1.05 | 0.96 | 1.03 | 0.96 | 1.05 |
| halogen | 0.94 | 1.23 | 0.94 | 1.14 | 1.11 | 1.23 | 1.13 | 1.14 |
| incandescent | 0.90 | 1.49 | 0.91 | 1.25 | 0.99 | 1.46 | 0.98 | 1.23 |
| solar | 0.61 | 1.23 | 0.58 | 1.01 | 0.95 | 1.21 | 0.97 | 0.97 |
| red | 0.60 | 1.98 | 0.58 | 0.94 | 1.12 | 2.61 | 1.06 | 1.51 |
| orange 0.25 + red 0.75 | 0.85 | 2.94 | 0.87 | 1.35 | 1.21 | 2.82 | 1.22 | 1.35 |
| **orange 0.5 + red 0.5** | 1.03 | **3.00** | 1.02 | 1.15 | 1.29 | 2.95 | 1.26 | 1.14 |
| orange 0.75 + red 0.25 | 1.37 | 2.71 | 1.31 | 0.98 | 1.33 | 2.70 | 1.27 | 0.98 |
| orange | 1.94 | 2.86 | 1.78 | 0.75 | 1.96 | 2.86 | 1.80 | 0.75 |
| blue 0.5 + red 0.5 | fail | 1.43 | fail | 0.67 | 0.91 | 2.26 | 0.95 | 0.79 |
| blue | fail | fail | fail | 0.49 | 0.71 | 1.25 | 0.74 | 0.56 |

(All 15 mixtures: `goals/<variant>.json`; the omitted ones lie between their
neighbours.)

**Robust in all eight minimal settings:**
- Orange–red light (620/680 nm) beats warm white; the 50/50 mixture gives 1.02–3.0×.
- Blue-rich light is worst, and it collapses the culture without chromatic
  acclimation.
- Cool white is never better than warm white.

**Not robust:** the size of the gain and the best orange:red ratio.
- Without capacity acclimation, pure orange is best (1.8–1.9×), matching pfennig2024's
  paired ranking (report 08: orange is the best single source).
- minimal-m0-cap gives ~3× for anything orange or red.
- minimal-m2-cap favours red and rates pure orange *below* warm white (0.75).

The capacity rule drives this: its target follows the PS1 absorption (`acclimation.py`),
so the light PS1 absorbs (Chl, red) raises capacity. That rule is an assumption fitted
to white light only; the spectral data can't test it.

`full-k3` can't serve as the reference for goal 2: 18 of 21 spectral runs failed (19 of 21 with `+ca`), and
its Sqs strain fails in the batch. Agreement with `full` (ranking criterion 2) is
therefore not available.

## 4. Findings

1. **Bistable pre-culture.** Pigment regulation gives a growing and a non-growing
   branch. With the Sqs strain, more than ~33 % blue in a blue/red mixture makes
   photoacclimation run away: μ → 0, Chl → 0.001 mM, SQ → 3000 mM. The solver had
   returned that branch as the pre-culture, which reported productivity 0 as a success.
   Chromatic acclimation (`+ca`) removes the collapse in all minimal variants.
2. **`double`: surrogate artefacts.** Its degree-6 polynomials (p95 error 28 % for Chl
   and SQ) produce SQ > 1000 mM, red light 13–20× warm white and μ = 0.13 h⁻¹ (twice the
   data maximum), while their inputs stay inside the training box. It fails all gates.
3. **`full-k3`: not reliable enough to be a reference.**
   - It fits Fig 9 only moderately (RMSE μ 0.0105 h⁻¹, ln N 0.17).
   - Chl comes out 8–64 % high and PC:Chl up to 43 % high.
   - Most of its batch and spectral runs fail. Of those that work, the WT reaches
     0.028 mg L⁻¹ h⁻¹, a third of `minimal`'s WT.
   - The low-light and dim-layer behaviour of the coupled pfennig2024 copies (report
     09) still needs work before it can arbitrate between the surrogate families.
4. **`single`:** the legacy PETC surrogate can't fit growth and pigments at once (Chl
   70 % low, PC:Chl 150–270 % high).
5. **Fallback solver.** The Sundials fallback evaluates the right-hand side in Python.
   Without a time budget it crawled towards t_max = 1e9 for hours (minimal-m0-cap
   under the fluorescent lamp, single-fixed's MCA). With the budget, the same full-k3
   fit took 24 min instead of failing to finish in 6 h.

## 5. Selection

**Recommended: `minimal-m0-cap`, with chromatic acclimation for goal 2.**

- It has the best fit of all variants by a wide margin (AIC 60 against ≥ 69 for
  `minimal` and ≥ 228 for the rest).
- It is the only variant passing gate 2, and it passes gate 1 (paired spectral error
  0.03).
- With `+ca` every spectral run converges. Its remaining gate-3 failures are
  turbidostat densities beyond the optimum and blue light without acclimation.
- It is cheap: about 6 min for all goal runs.

Its answers:
- **Goal 1:** the Sqs-overexpressing strain gives 0.40 mg L⁻¹ h⁻¹ and 38 mg L⁻¹ squalene
  after 96 h in the Fig 9 batch (WT 0.076). The turbidostat optimum is 3.2 mg L⁻¹ h⁻¹.
  The further MCA steps don't help in the batch.
- **Goal 2:** orange–red LED light, about 3× warm white in this variant. Across the
  minimal variants the gain spans 1.0–3.0× (structural uncertainty, dominated by the
  capacity rule).

`minimal-m2-cap` is the competitor. It has the best validation error (0.62) and a more
mechanistic PETC, but it misses gates 2 and 3 and disagrees on orange. Keep it as the
sensitivity case for goal 2.

**Export (goal 3).** `minimal-m0-cap` is the smallest variant: unrolled spectral sums
and no surrogate polynomials beyond the closed-form PETC. It needs no schema change
(ADR 0006 §2.7).

## 6. Decisions needed (checkpoint C)

1. Accept `minimal-m0-cap` although it misses gate 3 by 0.4 %? Alternatively, count the
   non-growing culture under blue light as the prediction "no growth" rather than a
   failure, which makes it pass.
2. Drop `double` and `single` from further work, keeping them as history.
3. `full`: park it, or invest in its low-light behaviour so it can arbitrate goal 2?
4. The capacity rule drives the goal-2 gain. Test alternatives (e.g. capacity driven by
   total absorbed light or by electron flow), or find data with coloured light?
5. Goal 1 follow-up: run the MCA under batch conditions (shaded, time-averaged), since
   the paper protocol picks steps that hurt the batch.

## 7. Open items carried over

- Biomass ATP 50 mmol gDW⁻¹ still needs a source (assumptions.md).
- Parallel-scan speed evaluation (report 07).
- mxlpy issue candidates: compile_fn `__module__`, `_CHECK_RTOL`. Native integrations
  can't be interrupted from Python, which is what hangs `single-fixed`'s MCA. Not
  filed; ask first.
