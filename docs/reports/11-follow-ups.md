# Report 11: checkpoint C follow-ups — batch MCA, capacity rule, colour data, `full`

**Status:** 2026-11-13. Follow-ups decided at checkpoint C (ADR 0006 §7). Two findings
need design decisions before the goal-2 answer can be trusted (§6).

## 1. MCA under batch conditions (goal 1)

`synmep.strains.optimise_batch` and `scripts/batch_mca.py` compute the control of the
batch's own 96 h mean productivity: finite differences (1 %) of the full Fig 9 batch.
A step is taken only if the realised overexpression raises the productivity. Results
are in `docs/reports/batch-mca/`.

| variant | batch sequence | P [mg L⁻¹ h⁻¹] WT → strain | control coefficients in WT |
| --- | --- | --- | --- |
| minimal-m0-cap | **Sqs** | 0.076 → 0.40 | Sqs +1.01, CrtE −1.00, Chl −0.32, CrtB +0.14, DXS −0.09, MEP −0.01 |
| minimal-m2-cap | Sqs, DXS | 0.085 → 0.65 → 0.651 | Sqs +1.21, CrtE −1.86, Chl −0.47, CrtB +0.29 |

- **The batch optimum is Sqs alone.** DXS adds only 0.09 %. CrtB has positive local
  control, but its realised fold lowers productivity.
- **CrtE has strongly negative batch control:** it draws FPP to GGPP and pigments.
  The paper protocol (1e7 cells mL⁻¹, score R(SQ) + 2R(μ)) picks it second, because
  pigment synthesis helps growth at low density. The goal-1 answer is therefore the
  Sqs strain: 0.40 mg L⁻¹ h⁻¹ and 38 mg L⁻¹ in 96 h (minimal-m0-cap).

## 2. Capacity rule alternatives (goal 2)

`minimal.get_model(capacity_signal=...)` takes one of three signals; references are
from fitted minimal-m0-cap at the Fig 9 median state (§ `CAPACITY_REF`):
- `ps1`: PS1 absorption, the original rule;
- `abs`: total absorption;
- `flux`: electron flow J_3PGA.

| variant | AIC | RMSE μ / ln N | Fig 10 μ | gates 1/2/3 | Sqs batch P | best light / warm white | orange | red | blue |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| minimal-m0-cap (`ps1`) | **60** | 0.0034 / 0.026 | 0.010 | ✓ ✓ ✗ (0.946) | 0.40 | orange 0.5 + red 0.5: 3.0 | 2.86 | 1.98 | fail |
| minimal-m0-cap-abs | 70 | 0.0045 / 0.034 | 0.015 | ✓ ✓ ✓ (0.97) | 0.48 | orange: 4.6 | 4.58 | 1.17 | 0.50 |
| minimal-m0-cap-flux | 107 | 0.0057 / 0.038 | 0.049 | ✓ ✓ ✓ (1.00) | 1.56 (†) | orange: 1.4 | 1.40 | 0.46 | 0.15 |

† 19× the WT batch, while its turbidostat (3.1) is like the others: probably a
different capacity branch. Not checked further.

The capacity rule moves the goal-2 gain from 1.4× to 4.6× and flips red from above
to below warm white. Orange (620 nm) is the best or near-best LED under every rule.
The white-light fit can't tell the rules apart well enough (AIC 60 / 70 / 107), so
data under coloured light are needed.

## 3. Coloured-light data: Zavřel et al. 2024

**Data.** Zavřel et al. 2024 (Plant Cell Physiol, doi 10.1093/pcp/pcae062); raw data
from figshare 10.6084/m9.figshare.23330363, saved as
`data/zavrel2024_light_quality.csv`.
- Synechocystis WT under nine narrow-band LEDs (435–687 nm), 25 µmol m⁻² s⁻¹, batch
  flasks, OD750 0.2.
- Measured: growth rates, Chl, PBS, carotenoids, PSII/PSI.
- The growth column is labelled h⁻¹, but its values (0.08–0.33) only make sense as
  d⁻¹, so only relative growth is used.

**Comparison.** `scripts/colour_growth.py`, results in
`docs/reports/colour-growth.json`. Model setup: WT steady state at 2.1e7 cells mL⁻¹,
1 cm, 25 µE, one LED at a time. Growth relative to the 633 nm LED:

| | 435 | 465 | 495 | 520 | 555 | 596 | 633 | 663 | 687 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **data** | 0.62 | 0.23 | 0.39 | 0.59 | 0.83 | 0.92 | 1 | 1.01 | 0.71 |
| minimal-m0 | ✗ | ✗ | ✗ | ✗ | ✗ | 0.60 | 1 | 0.72 | ✗ |
| minimal-m0 +ca | 1.10 | ✗ | ✗ | ✗ | ✗ | 0.61 | 1 | 1.08 | 1.08 |
| minimal-m0-cap | 2.86 | ✗ | ✗ | ✗ | ✗ | ✗ | 1 | 1.91 | 1.65 |
| minimal-m0-cap-abs | 1.91 | ✗ | ✗ | ✗ | 0.57 | 0.76 | 1 | 1.58 | 0.55 |
| minimal-m0-cap-flux | 2.26 | ✗ | ✗ | ✗ | ✗ | 0.86 | 1 | 1.05 | 1.40 |
| minimal-m0-cap-flux +ca | 1.18 | ✗ | ✗ | ✗ | ✗ | 0.23 | 1 | 1.09 | 1.04 |
| minimal-m2 | ✗ | ✗ | ✗ | ✗ | ✗ | 0.59 | 1 | 0.71 | ✗ |
| minimal-m2-cap | 1.12 | ✗ | ✗ | ✗ | ✗ | ✗ | 1 | 1.12 | 0.99 |

✗ = no growing culture (μ < 1e-3 h⁻¹). The cells in the data grow under every LED.

**No variant reproduces the colour data**, under any capacity rule, with or without
chromatic acclimation. There are two failures.

1. **Dim or weakly absorbed light kills the model cells.** At 25 µE, blue-green and
   green light (465–555 nm) leave no growing state.
   - The cause is the pigment regulation. Pigment-branch enzymes scale with
     (μ/μ_ref)^c, and the shared fit drives c to **6.6–13** in every `minimal`
     variant. So 10 % less growth cuts pigment synthesis by about 60 %, which cuts
     absorption and so growth further. Below a threshold the cell runs into the trivial
     state without pigments or growth. The same mechanism made the blue-light
     pre-culture collapse in report 10.
   - Real cells do the opposite in dim light: they *raise* pigment content
     (photoacclimation).
   - The large c is what lets the fit reproduce the Fig 9 growth curve: pigment
     synthesis falls as shading slows growth.
2. **Light absorbed by Chl, which mostly feeds PS1, is overrated.** Violet (435 nm) is
   predicted at 1.1–2.9× the 633 nm rate, against 0.62 measured. Red and far-red behave
   alike. Neither m0's weighted light sum nor m2's balance captures the poor use of
   PS1-absorbed light at low intensity reported by Luimstra et al. 2018 and Zavřel
   et al. 2024.

Consequence: the goal-2 answers of reports 10 and 11 (orange best, 1.4–4.6×) are not
supported by the one dataset that tests the colour response directly. "Orange is
best" agrees with the data's 633 nm optimum, but the size of the gain and the ranking
of blue, red and violet do not.

## 4. `full`: no true steady state at batch light

Diagnosis in the WT pre-culture model (Fig 9 median density, batch light, k = 3
layers):

- **Lumen protons** (`z*_Hi`) rise linearly in every layer, by about 0.01 mM s⁻¹ in the
  surface layer: 0.2 → 60 000 mM over 6e6 s.
  - The ATP pool is saturated (≈ 428 mM of the adenylate total), so ATP synthase
    stalls while PSII and cytochrome b6f keep pumping.
  - The passive leak (`vPass`) is orders of magnitude too small.
  - Report 09 replaced pfennig's lumped ATP consumption by the real (growth and
    maintenance) demand. Since then, any excess light energy has no outlet:
    pfennig2024 has no strong photosynthetic control (ΔpH back-pressure on b6f) to
    throttle the supply.
- **CO₂** builds up to ~2 × 10⁵ mM in the surface layer.
- **Dark layer (`z2`):** NAD falls to ~0 (all NADH). Photorespiratory salvage needs NAD,
  so it stalls and PG grows without bound.

Because no true steady state exists, `full`'s steady-state solvers, Fig 10, the
turbidostat and most batch runs fail (report 10: 30 % convergence). The Sqs strain
itself is fine: its pre-culture settles to μ ≈ 0.018 h⁻¹ before the slow drift.

Infrastructure changes (commits f28f2a5, 6712d4a):
- **Fallback budget:** a 60 s wall-clock budget for the Sundials fallback.
- **Pre-culture history:** the pre-culture now follows a culture's history when the
  direct solve fails (WT under warm white, then light and folds stepped to the target).
- **Settled test:** the long simulation counts as settled on < 1e-4 relative change
  over its last tenth. This correctly rejects `full`'s drifting states.
- **Effect on `full-k3`:** its fit now takes 24 min instead of failing to finish in 6 h.

## 5. Status of the selection

`minimal-m0-cap` stays the best fit to the white-light data, and its goal-1 answer
(Sqs strain, 0.40 mg L⁻¹ h⁻¹ in the Fig 9 batch) is robust across the minimal variants
and the batch MCA. Its goal-2 answer is not trustworthy until §3's failures are
fixed.

## 6. Decisions needed

1. **Pigment regulation.** Replace growth coupling (μ/μ_ref)^c with photoacclimation:
   pigment synthesis rising at low absorbed light per Chl, or driven by a redox or
   excitation signal. Refit everything on it. This touches every family that is
   still active (`minimal`, `full`).
2. **Use the Zavřel colour data as a fit or validation target** for goal 2 (relative
   growth under nine LEDs, plus Chl and PBS per cell and PSII/PSI). Fitted, they would
   pin the capacity rule and the PS1/PS2 use. Held out, they would validate it.
3. **`full`: how to give excess energy an outlet.** Options:
   - (a) restore an ATP overflow sink (pfennig's lumped consumption, saturating);
   - (b) add photosynthetic control: lumen-pH-dependent b6f rate, or a stronger
     proton leak;
   - (c) both.

   Plus a dark-layer NADH re-oxidation (respiratory terminal oxidase path) and a check
   of the CO₂ exchange. Any of these changes pfennig2024 and is reported as such.
