# Report 10 (interim): comparison of all variants — status and next steps

**Status:** interim, 2026-11-12. Fits and goal runs were stopped on request before
completion. This records what is done, what the partial results say, and how to
continue. The final report 10 (checkpoint C, ADR 0006 §2.9) replaces this file.

Numbers: `docs/reports/summary.json` (from `scripts/summary.py`), per-variant
`docs/reports/fits/<variant>.json` and `docs/reports/goals/<variant>.json`.

## 1. What is complete

| variant | shared fit | goal runs (with guards, commit 18dd606) |
| --- | --- | --- |
| single-fixed | done (after the pc_per_chl fix) | **missing** |
| single-fixed-cap | done | **missing** |
| double-none / -cap | done | done |
| double-coupled / -cap | done | done |
| minimal-m0 | done | done |
| minimal-m0-cap | done | **missing** (earlier run timed out under load; restarted run was stopped) |
| minimal-m2 / -cap | done | done |
| full-k3 | **stopped** during least squares (grid done, see §4) | **missing** |

`single-legacy` is the gate-2 reference only (no shared fit, by design).

## 2. Partial results

### Gates (ADR 0006 §2.9, §5)

| variant | k | RMSE μ [h⁻¹] | RMSE ln N | AIC | validation ln N | Fig 10 μ | gate 1 (spectral) | gate 2 (Fig 9) | gate 3 (conv.) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| minimal-m0 | 3 | 0.0059 | 0.054 | 99 | 0.66 | 0.011 | ✓ 0.030 | ✗ (μ 2 % over) | ✓ 0.95 |
| minimal-m0-cap | 5 | **0.0034** | **0.026** | **60** | 0.76 | 0.010 | ✓ 0.030 | ✓ | pending |
| minimal-m2 | 3 | 0.0058 | 0.054 | 100 | 0.67 | 0.010 | ✓ 0.014 | ✗ (at the limit) | ✗ 0.92 |
| minimal-m2-cap | 5 | 0.0065 | 0.054 | 69 | **0.62** | 0.017 | ✓ 0.014 | ✗ | ✗ 0.89 |
| double-none(-cap) | 2 (4) | 0.013 (0.008) | 0.24 (0.18) | 300 (256) | 0.80 (0.93) | 0.022 (0.041) | ✗ 0.10 | ✗ | ✗ 0.86 (0.95) |
| double-coupled(-cap) | 2 (4) | 0.015 (0.011) | 0.29 (0.23) | 328 (292) | 0.76 (0.88) | 0.022 (0.037) | ✗ 0.10 | ✗ | ✗ 0.88 (0.84) |
| single-fixed(-cap) | 3 (5) | 0.012 (0.009) | 0.53 (0.39) | 445 (385) | 0.99 (0.92) | 0.018 (0.019) | ✗ 0.10 | ✗ | pending |

No variant passes all three gates yet. `minimal-m0-cap` passes gates 1 and 2 and is
the clear favourite if its goal runs converge (≥ 95 %). `minimal-m0` misses gate 2 by
2 % in μ.

### Goal answers (best strain, Fig 9 batch conditions, warm white)

| variant | MCA sequence | best strain | batch best P [mg L⁻¹ h⁻¹] | titer 96 h [mg L⁻¹] | turbidostat P | best light (batch / warm white) |
| --- | --- | --- | --- | --- | --- | --- |
| minimal-m0 | Sqs | Sqs | 0.34 | 33 | 1.95 | orange 1.94× |
| minimal-m2 | Sqs, CrtE | Sqs | 0.35 | 33 | 2.92 | orange 1.78× |
| minimal-m2-cap | Sqs, CrtE, MEP, DXS | Sqs | 0.65 | 62 | 1.81 | orange 0.25 + red 0.75: 1.35× |
| minimal-m0-cap (earlier run, no guards) | Sqs, CrtE, MEP, DXS | Sqs | 0.40 | 38 | 3.24 | — |
| double-* | Sqs, CrtE, MEP(, DXS) | Sqs (+CrtE) | 0.17–0.47 | 11–37 | 0.6–4.1 | implausible (up to 20×, §3) |

Consistent so far across the `minimal` variants:
- **Goal 1:** Sqs overexpression is the step that counts in the batch. The further MCA
  steps (CrtE, MEP, DXS), chosen at the paper's MCA conditions (1e7 cells mL⁻¹,
  100 µE), raise the turbidostat optimum but *lower* batch productivity. In
  minimal-m0-cap the batch falls from 0.40 (Sqs) to 0.21 mg L⁻¹ h⁻¹ (all four). The
  low-density MCA protocol doesn't carry over to the shaded batch, so an MCA under
  batch conditions is a candidate follow-up.
- **Goal 2:** without capacity acclimation, orange (620 nm) is best, at about 1.8–1.9×
  warm white, as the pfennig2024 paired ranking predicts (report 08: orange 1.47 at
  the PETC level). Chromatic acclimation hardly changes the orange optimum but removes
  the blue-light failures. With capacity acclimation (m2-cap) the optimum moves to
  red-rich mixtures and the gain shrinks to 1.35×: **the goal-2 answer depends on the
  acclimation option.** That is the main structural uncertainty to resolve, ideally
  against `full`.

## 3. Findings from this session

1. **Dead pre-culture branch (fixed, commit 18dd606).** `minimal-m0` with the Sqs strain
   has a fold bifurcation along the blue fraction (at ~33 % blue in blue/red mixtures).
   Photoacclimation lowers Chl as absorption per Chl rises, which runs away to μ → 0,
   Chl ≈ 0.001 mM and SQ → 3000 mM (no dilution). The steady-state solver returned this
   state as the "pre-culture", giving productivity 0 counted as success. Now the
   pre-culture must grow (μ ≥ 1e-3 h⁻¹; it retries along a warm-white → target light
   path), and those runs count as failures. Chromatic acclimation removes this
   collapse in the minimal variants.
2. **Squalene mass-balance guard (fixed, same commit).** A batch with SQ above 25 % of
   the dry weight (282 mM) is a failure. The `double` surrogates hit SQ > 1000 mM and
   productivities up to 10¹⁵× warm white under some LED mixtures. Their inputs stay
   inside the training box, but the degree-6 polynomials (p95 error 28 % for Chl and SQ)
   misbehave in sparsely sampled corners. Even after the guard, `double-*-cap` predicts
   red light 13–20× warm white with μ = 0.13 h⁻¹ (twice the Fig 9 maximum).
   **Recommendation: drop `double` from selection.** It fails all three gates, and its
   goal-2 answers are artefacts.
3. **single-fixed** can't fit Fig 9 together with the pigment targets: Chl is 70 %
   low and PC:Chl 150–270 % high, with ln N RMSE 0.39–0.53. Expected under ADR 0006
   §4: the legacy PETC surrogate fails gate 1 anyway.
4. **Proteomics-batch validation** gives ln N errors of 0.62–0.99 for *every* variant.
   That culture differs systematically from Fig 9, so validation separates the
   variants only weakly; AIC and Fig 10 discriminate more.
5. **Load:** 11 goal runs plus fits in parallel slowed a single run from seconds per
   batch to a 3 h timeout. Run at most ~6 goal runs at once on this machine (22 cores).

## 4. `full-k3`

The grid (17 points, 27 min) found best cost 863 at factors pbs_scale 4.88,
pigment_coupling 0.025 (near its lower bound 0.02), V_bio 5.76 and vCBB_max 0.35. For
comparison, minimal-m0's best cost is about 50. The least-squares stage (max_nfev 15)
was stopped. A pigment_coupling at its bound suggests that `full`'s pigment
regulation, or the bound itself, needs a look before more fitting.

## 5. Next steps (in order)

1. **Finish the missing goal runs:** `minimal-m0-cap`, `single-fixed`,
   `single-fixed-cap` (each ~10–60 min alone):
   `OMP_NUM_THREADS=1 uv run --extra surrogate python scripts/goal_runs.py <variant>`.
   Gate 3 of `minimal-m0-cap` decides whether a variant passes all gates.
2. **full-k3 fit:** restart from the grid optimum to skip the 27 min grid:
   `uv run --extra surrogate python scripts/shared_fit.py full-k3 --max-nfev 15 --start '{"pbs_scale": 4.88, "pigment_coupling": 0.0246, "V_bio": 5.76, "vCBB_max": 0.352}'`
   (expect ~2–3 h). First check whether lowering the `pigment_coupling` bound is
   justified. Then run `goal_runs.py full-k3` (long: run it alone).
3. `uv run python scripts/summary.py`: gates, ranking (now with the fallback by number
   of gates passed), agreement with `full` (`rho_full`, Spearman over light sources),
   spread, figures `summary-fits.png` and `summary-spectra.png`.
4. Decide with the user:
   - whether `double` is dropped (recommended, §3.2);
   - whether gate 2 is relaxed by the 2 % that `minimal-m0` and `minimal-m2` miss;
   - whether capacity acclimation is kept as default, since it changes the goal-2
     answer (§2);
   - whether to add an MCA under batch conditions (§2, goal 1).
5. Write the final report 10 and update ADR 0006 (status, §7 outcome), CONTEXT.md and
   lessons-learned.md (pre-culture bistability, surrogate-corner artefacts, parallel
   load).

Deferred from earlier (unchanged): parallel-scan speed evaluation (report 07), the
source for biomass ATP 50 mmol gDW⁻¹ (assumptions.md), and the mxlpy issue
candidates (compile_fn `__module__`, `_CHECK_RTOL`), which are not to be filed without
asking.

Untracked caches that save time when continuing: `data/surrogates/metabolism/*.json`
(per-strain metabolism surrogates of `double`, retrained automatically if missing).
