# ADR 0005: Minimal Model Design (v3)

**Status:** Accepted (implemented, see §5 and `docs/reports/06-v3-minimal.md`)
**Scope:** `src/synmep/models/v3_minimal.py`; supersedes the v3 row of ADR 0001

---

## 1. Context

ADR 0001 planned v3 as "the smallest model that still explains the data", with
acceptance criteria on Figs 9–12. Since then the purpose of the work has been made
explicit: **predict the maximal squalene productivity of a culture, given the
squalene/pigment trade-off, and the strain (overexpression set) and operating point
that achieve it.** The work on v0/v1 and the ground-up thread g (ADR 0004) produced
the findings this design rests on (`docs/lessons-learned.md`, reports
`03-g-thread.md`, `04-strain-optimisation.md`, `05-productivity.md`).

## 2. Decision

### 2.1 What the minimal model must answer
1. Batch growth (Figs 9/10 data): growth-rate and density time courses.
2. Strain optimisation (Fig 12, each enzyme once): ranking of overexpressions.
3. **Culture productivity**: optimum over strain × cell density (turbidostat), with
   the squalene/pigment trade-off. This is the primary question.

### 2.2 What it must contain (and why)

| element | form | evidence |
| --- | --- | --- |
| light | self-shading with **depth resolution** (few Gauss–Legendre nodes; coarse spectrum or broadband) | depth-averaged light into a nonlinear PETC needed a 5.5x fudge factor; depth resolution removed it |
| PETC | **2 outputs** (carbon and ATP supply) from absorbed light; closed-form saturating light response, linear to 0 in darkness | PETC has effectively 2 independent outputs; pigments act only via absorption; polynomial surrogates broke off-manifold |
| growth | from the **carbon/energy balance**: μ = biomass synthesis / dry weight | gets the measured magnitude without fitting; no Chl pinning |
| **pigment → growth link** | Chl (and carotenoid) content as states, set by synthesis vs. dilution; content sets light absorption per cell and self-shading | without it there is no productivity optimum: fixed-pigment variants saturate (productivity report) |
| pigment regulation | bounded (e.g. growth-coupled pigment branch, variant A) | unregulated pigments fall into a "dark trap" in dense cultures (no steady state) |
| MEP/terpenoid branch | **DXS**, lumped MEP chain, **CrtE** (FPP and GGPP), **Sqs**, pigment sinks (CrtB, ChlP/ChlG lumped) | MCA and the strain productivity scan pick Sqs, CrtE, MEP chain, DXS in every thread-g variant; the best realistic strain is Sqs + CrtE (+ MEP, DXS) |
| feedback | DXS inhibition by IPP/DMAPP (Ki 65 µM) | without it the chain is unbounded at low growth |
| enzyme cost | protein synthesis cost of overexpressed enzymes (cofactor/carbon drain or proteome share) | thread g has free enzymes, so its free optimum (100–460x Sqs) is unbounded by cost |
| capacity acclimation | optional slow state for photosynthetic capacity | best batch fit (variant C) and partly explains the early growth rise; not needed for steady-state productivity |

### 2.3 What it can drop
- The 1410-parameter cross product (v0) and per-enzyme synthesis kinetics; enzymes are
  levels (overexpression = fold change), plus a cost term.
- Explicit cofactor pools beyond 3PGA/ATP/NADPH, the ATP⇌NADPH/Fd interconversions, Fd as
  a separate pool, ThiDP, Idi (IPP/DMAPP lumped), phytol intermediates.
- Polynomial surrogates for metabolism: the minimal model is small enough to simulate the
  metabolism directly (or with closed-form QSS) — the 8-input polynomial surrogate reached
  only ~10–14 % median error and broke at limitation switches.
- Hopanoids (Δshc strain) unless the WT is needed.

### 2.4 Acceptance criteria (revised)
- Batch fit (Figs 9/10): RMSE no worse than the best thread-g variant + 10 %
  (μ ≤ 0.0085 h^-1, ln N ≤ 0.155), refitting allowed.
- Strain optimisation: same first three overexpressions as thread g (Sqs, CrtE, MEP/DXS).
- Productivity: an interior optimum in squalene synthase level and in density, within a
  factor ~2 of thread g's optimum, and the same best realistic strain.

## 3. Rationale
- The productivity question requires the trade-off, and the trade-off requires the
  pigment → light-harvesting → growth link. The batch fit alone preferred constant
  pigments (variant B, a → 0), which would remove exactly that link; so the batch fit
  can't be the only criterion.
- Every element in §2.2 is there because removing it demonstrably broke something in the
  larger models.

## 4. Consequences
- v3 is built from thread g (not from v1), which already has the right structure; v1
  stays as the bug-fixed reference of the legacy model.
- Decided at implementation (§5): both enzyme-cost formulations and both capacity
  options are kept as variants; the dry weight is taken from the measured protein.

## 5. Implementation decisions

### 5.1 Variants
`v3-<cost>[-cap]`: enzyme cost `drain` or `proteome`, capacity acclimation off/on
(user decision: "try both" for each). With Δ = Σ_i (E_i − 1) f_i the extra protein as a
fraction of the proteome (f_i: WT proteome mass fraction of enzyme i, from the
proteomics):
- **drain**: the extra protein is synthesised at μ Δ P_cell and drains 3PGA, ATP and
  NADPH with the biomass composition (1 g protein costs like 1 g biomass).
- **proteome**: growth and photosynthetic capacity scale with (1 − φ_Q)/(1 − φ_Q + Δ),
  φ_Q = 0.45 (Scott2010; assumption).
The wild type has Δ = 0, so one batch fit per capacity option serves both costs.

### 5.2 Dry weight per cell volume from the measured protein
The gravimetric dry weight of the proteomics samples (median 3.88 pg cell^-1, i.e.
0.97 g mL^-1 at 4 fL; cells are ~0.3–0.5 g mL^-1) contains only 22 % protein (measured
on the same samples; Synechocystis: ~45–55 %). The protein measurement, the cell count
calibration and 4 fL are mutually consistent; the gravimetric dry weight isn't (most
likely residual medium salts or water in the pellets, low-OD samples have ~1 mg on a 1 g
tube). Decision: DW per cell = measured protein per cell (0.83 pg) / 0.445 = 1.85 pg,
ρ_DW = 464 g L^-1 (was 970). Per-cell volume and cells/OD stay as in the paper. Consequences:
- μ per unit supply doubles; absorbed into the fitted PETC efficiency.
- Proteomics enzyme concentrations computed with DW × 0.445 are ~2× too high; v3 only
  uses the basis-free proteome fractions f_i, and v0/v1/thread g keep their values.
- 7.5 mg Chl OD^-1 L^-1 (source unknown) now means 3.8 % Chl in DW, high for
  Synechocystis (typically 1–2 %); the Chl/OD value should be checked.

### 5.3 PETC
Closed form fitted to the pfennig2024 training states (`scripts/fit_petc_closed_form.py`):
one tanh light response in a weighted absorbed light I = a_PS1 + w_PS2 a_PS2 + w_PBS a_PBS,
and J_ATP quadratic in J_3PGA (pfennig2024's J_ATP is nearly a function of J_3PGA).
Hold-out R² 0.998 / 0.979 vs. 0.983 / 0.939 for the degree-5 polynomial, 8 parameters,
valid down to darkness (no ray extension needed).

### 5.4 Pigment regulation
Variant A (pigment-branch enzymes ∝ (μ/μ_ref)^c) with the exponent c fitted (bounded
0.1–10). Fitted c = 5–9: pigment synthesis strongly follows growth.

### 5.5 Outcome
All variants meet the batch-fit and strain-ranking criteria; with capacity acclimation
also the productivity criteria except a strictly interior density optimum (productivity
peaks at the densest steady state). The cost formulation is immaterial at realistic
overexpression. Recommended minimal model: `v3-proteome-cap`. Details: report 06.
