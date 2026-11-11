# ADR 0006: Goals, Model Families and How the Final Model Is Selected

**Status:** Accepted (2026-11-11), implementation in progress
**Scope:** the whole repository; supersedes the goal statements of ADR 0001 §1 and
ADR 0005 §2.1, and the variant names used up to report 07

---

## 1. Context

Reports 01–07 answered a narrower question than the project needs: culture
productivity in a **turbidostat** under **warm white light**. The goals have been
restated, and the model variants had grown into a set with opaque names (v0, v1, "g",
v3) mixing versions and structures. This ADR fixes the goals, the families, the
calibration protocol and the selection criteria for the summary report.

## 2. Decision

### 2.1 Goals

1. **Maximal squalene productivity of an MCA-optimised strain in a batch culture**
   (self-shading included). Reported as
   - (a) mean volumetric productivity titer(t_h)/t_h at the best harvest time t_h,
   - (c) titer after a fixed duration (96 h, as in the data),
   - and, as a reference for continuous operation, the turbidostat optimum (report 05
     setting).

   Conditions: the Fig 9 batch (MC-1000, 1.35 cm, 80 µmol m^-2 s^-1 warm white LED,
   data inoculum). Only the strain is optimised. The five other lamps in
   `data/lights/` are run at the same intensity as a first look (the lamp comparison
   itself is follow-up work).
2. **Best light composition for squalene production**: mixtures of narrow-band LEDs
   (Gaussians at ~450, 620, 680, 730 nm) at fixed total photon flux, with the six
   lamps as baseline points. This requires the PETC part of a model to respond
   non-additively to the PS1/PS2 excitation balance.
3. **Export to mxlweb-core** (green sloth, mxl-models) is a later consideration.
   It is reported per variant (with the schema extension it would need), not used to
   select.

### 2.2 MCA protocol

The paper protocol is kept (Fig 12): steady state at 1e7 cells mL^-1,
100 µmol m^-2 s^-1, 1 cm; largest summed response of squalene and biomass synthesis,
each enzyme at most once, only steps with positive control on intracellular squalene
that raise it when realised (report 04). Every strain of the sequence is then
evaluated as a batch (2.1 a, c) and as a turbidostat; the best is reported.

### 2.3 Families and variants

| family | variants | what it is |
| --- | --- | --- |
| `single` | `single-legacy` (paper baseline, bugs included), `single-fixed` (recalibrated) | one surrogate for PETC + light; metabolism and culture intertwined (the original model) |
| `double` | `double-none`, `double-coupled` (pigment regulation) | surrogate PETC, surrogate metabolism, dynamic culture (formerly thread "g") |
| `minimal` | `minimal-m0` (formerly v3, proteome cost), `minimal-m1`, `minimal-m2` | smallest model that fits the data and answers both goals |
| `full` | `full-k3` (or `-k5`) | pfennig2024 included directly, no surrogates |

Options available to every family: **capacity acclimation** (`-cap`) and
**balance-seeking chromatic acclimation** (2.6). Retired (history only):
`double-photoacclimation` (fitted to fixed pigments), `double-capacity` (no pigment
trade-off), the drain enzyme cost (< 0.3 % different from the proteome cost).

Old reports are rewritten with these names, so the document set is consistent.
Old → new: v0 → `single-legacy`, v1 → `single-fixed`, g-* → `double-*`,
v3-proteome(-cap) → `minimal-m0(-cap)`.

### 2.4 Minimal family

- `minimal-m0`: v3 (ADR 0005), closed-form PETC on one weighted light sum. Kept as the
  reference for "no PS1/PS2 balance".
- `minimal-m1`: closed form with balance, no extra states:
  J = J_max tanh(smin(e1, e2) / I_k), e1/e2 the effective PS1/PS2 excitation, PBS
  excitation distributed between them, smin a smooth minimum.
- `minimal-m2`: mini electron-transport model with 2–3 fast states (PQ, PC/P700,
  possibly lumen pH), CBB as sink.

m1 and m2 are fitted to pfennig2024 steady states resampled to include LED mixtures;
m1 is preferred if its hold-out R² on imbalanced spectra is ≥ ~0.97.

### 2.5 Full family

- **Depth:** K copies of the pfennig2024 states, one per Gauss–Legendre depth layer
  (K = 3 vs 5 checked first), each with its local light.
- **Pools shared with metabolism:** pfennig's 3PGA, ATP and NADPH are the metabolism's
  pools; the lumped consumption reactions are replaced by the real demand. The proton
  accumulation seen with a crude 3PGA drain has to be solved, not avoided.
- **Fast per layer, slow shared:** PETC states, 3PGA, ATP/ADP, NADPH/NADP+ and their
  direct consumers exist per layer; slower pools (MEP/terpenoid intermediates,
  pigments, squalene, enzymes, μ) once per cell. Layer exports are averaged with the
  quadrature weights. The cut is set from the Jacobian (relaxation < ~1 min → per
  layer).

### 2.6 Light and pigment composition

- `double` gets its PETC surrogate retrained on the LED-mixture training data (it was
  trained on white-light states only). `single` keeps the legacy surrogate.
- Pigment composition: fixed (warm-white values) as baseline; variant with a slow
  state moving PS1:PS2 towards balanced excitation, checked against literature values
  for colour-acclimated cells (pfennig2024's sources where available).

### 2.7 Export form

Spectral, state-dependent light stays as **unrolled scalar sums** on a 5 nm grid; the
LED channel weights are parameters. No schema change needed.

**Possible follow-up:** mxl-schemas ADR 0014 only allows array subtrees that reference
data and numbers (constants). Allowing state-dependent array subtrees (array-valued
codegen in the RHS of the WAT/C backends) would let the spectra live in sidecar data
and keep the model files small. To be discussed if export size or browser compile time
becomes a problem.

### 2.8 Calibration

All variants: one fit to Fig 9 (μ, ln N) plus pigment contents per cell from the OD
spectra (WT, sqs), at most ~5 free parameters. Fig 10 is held out for validation. The
measured squalene content of the sqs strain is a check only (it is the reference scale
of `double` and `minimal`).

### 2.9 Selection (summary report)

Gates:
1. captures the PS1/PS2 balance (non-additive LED response),
2. Fig 9 RMSE ≤ `single-legacy` + 10 % (μ ≤ 0.0058 h^-1, ln N ≤ 0.149),
3. ≥ 95 % convergence in the MCA, batch and LED scans.

Ranking among those passing: (1) Fig 10 validation error and AIC, (2) agreement with
`full` on the goal answers, (3) cost of a complete goal-1 + goal-2 run. The spread of
the goal answers over all passing variants is reported as structural uncertainty.
`full` stays the reference for (2) even if it fails a gate.

### 2.10 Stages and checkpoints

1. Rename everything (code, scripts, JSON, figures, reports 01–07); delete retired
   variants.
2. Shared infrastructure: pigment data, Fig 10 validation, batch/turbidostat evaluator,
   LED inputs, acclimation options.
3. LED-mixture PETC training; `minimal-m1`/`m2`; `double` retrain. **Checkpoint A.**
4. `full` (may run in parallel to 3). **Checkpoint B.**
5. Shared fit and gates.
6. Goal runs.
7. Summary report. **Checkpoint C.**

## 3. Rationale

- Batch productivity is the standard bioprocess figure and includes the shading build-up;
  the turbidostat optimum shows what continuous operation could add.
- Goal 2 is only meaningful if light is not interchangeable in the model; hence the
  balance gate and the LED training data.
- The `full` model with the same interface as a surrogate would only measure surrogate
  error; sharing pools also captures feedback of the sink on the light reactions.
- One calibration protocol and fixed gates make the comparison reproducible instead of
  a matter of taste.

## 4. Consequences

- Every family needs a dynamic batch mode that runs reliably.
- `single` will probably fail the balance gate (legacy surrogate); that is a result.
- The `full` model is the largest risk (stiffness, pfennig changes).
- Reports 01–07 and all scripts change names; earlier JSON results are regenerated.
