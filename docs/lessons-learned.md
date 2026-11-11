# Lessons learned (input for v3 minimisation)

Running list kept while working on stages 0–2. Each entry: observation → implication
for the minimal model.

## Stage 0 — audit

- **1410 parameters for 39 variables / 85 reactions.** Most are Km/Ki/Kd combinations
  for every (enzyme, compound) pair that no rate law uses. → Build parameters from the
  rate laws, not from a cross product.
- **Chl steady state is pinned by the growth law** (`Chl* = 1/(slope·V)`). The kinetics
  of chlorophyll synthesis only set μ, not Chl. → In a minimal model μ can be expressed
  directly as a function of the MEP-pathway flux into Chl. The Chl pool dynamics matter
  only for transients and self-shading.
- **The biomass drain is negligible in v0** (~0.01 % of 3PGA); cofactor pools are set by
  surrogate supply vs. arbitrary first-order `*_out_metabolism` sinks. → The four
  cofactor ODEs (ATP, NADPH, Fdred, 3PGA) mostly reproduce `supply/k_out`. They are
  candidates for algebraic elimination, but only once the biomass drain is corrected
  (that correction can make cofactors growth-limiting).
- **Fast intermediates** (HMB_PP, GPP, DMAPP, IPP: 15–150 min^-1) vs. growth
  (7.7e-4 min^-1): stiffness ratio ≈ 2.3e5. → Natural QSS candidates (supports v2a).
- **Prenyl intermediates are tiny** (µM range) and only pass flux on. → Linear-chain
  lumping candidates (DXP→…→HMB_PP; GPP→FPP).
- **Enzymes relax at the dilution rate** (all ≈ 7.7e-4 min^-1 diagonal), so they are slow
  states. For steady-state questions they are just `0.1·ksynth·sat/μ`. → Possibly replace
  with constants for every enzyme except the ones the Fig 11/12 experiments act on.
- **The surrogate costs ≈ 10 ms per RHS call** (sklearn + pandas + Simpson) and dominates
  runtime (26 s for one steady state). → Precompute everything that doesn't depend on
  state, and make the rest explicit numpy.
- **The PETC surrogate has only two independent outputs.** NADPH : 3PGA : Fdred are always
  1 : 8.96 : 2.96 (the biomass stoichiometry); ATP = that flux × 33.09/8.76 + an offset
  (~2 at zero light). → A minimal model needs two supply terms (ATP, "biomass
  precursor flux"), not four.
- **Every observed failure came from surrogate extrapolation**, not intrinsic stiffness.
  → Any surrogate (including the v2 ones) must ship with its training domain and be
  bounded/smooth outside it. Prefer low-degree or monotone forms over degree-7 polynomials.

## Stage 1 — port

- **v0 has 195 used parameters** (of 1410 legacy ones). `get_unused_parameters()` is a
  cheap first cut for any variant.
- **Native compile cost scales with expression size.** The 301-wavelength sums and
  4 × 792-term polynomials take ~20 s to compile. → A minimal model should use a small
  surrogate (few inputs, low degree) and a coarse or analytic light model; this also
  makes parallel scans cheap (each worker compiles once).
- **The legacy model is fully reproducible from the CSVs + one coefficient file**: no
  hidden state in the pickle beyond the polynomial. Every later variant can start from
  `load_legacy_inputs()`.
- **`Simulator.update_parameter` mutates the model passed to it** (expected mxlpy
  behaviour). A hand-written serial loop reusing one model object carries parameter
  changes (e.g. the Fig 11 CrtE knock-down) into later runs. → Use `scan.*` (protocols,
  time courses, steady states) for repeated runs, or a fresh model per run.

## Thread "g" (ground-up model)

- **Growth from the carbon balance works without fitting the magnitude**: PETC 3PGA supply
  per g dry weight vs. the biomass stoichiometry gives μ ≈ 0.03 h^-1, the observed batch
  range. → A minimal model can drop the empirical μ ∝ v_Chl law.
- **The PETC has effectively 2 outputs and 4 inputs, and pigment contents only act
  through absorption.** → Minimal model: one "absorbed light → carbon/ATP supply" curve.
- **Surrogate inputs must come from the physical manifold**, otherwise chained
  surrogates extrapolate at exactly the operating points that matter.
- **Predict bounded quantities** (contents X*, not fluxes) when a downstream division by
  a vanishing rate (μ) is involved.
- **Depth-averaged light into a nonlinear PETC kills growth in dense cultures.** → A
  minimal model needs at least a coarse depth resolution (or an analytic depth integral of
  a saturating light-response curve).
- **Pigment contents need regulation** to stay physiological at low growth.
- Global polynomials hit their limits at limitation switches (p95 errors 15–20 % even
  in 3D). → For the minimal model, prefer closed-form saturating functions with few
  parameters over polynomial surrogates.
- **Constant pigment contents fit the batch data best** (variant B fitted acclimation
  strength → 0). → The minimal model can drop pigment dynamics entirely (Chl, Car as
  parameters), keeping only self-shading through cell density.
- **The early growth rise is not explained by pigments.** → If the minimal model must
  capture it, it needs a slow photosynthetic-capacity state, not pigment states.
- **Depth resolution matters more than PETC detail**: it removed a 5x supply factor.
- **Strain-optimisation advice is model-dependent but robust within thread g**: v0 →
  IspG, v1 → DXS, every thread-g variant → Sqs, CrtE, MEP, DXS. → For the minimal model,
  keep Sqs, CrtE (FPP/GGPP branch point) and DXS explicit; the MEP chain can stay lumped.
- **Capacity acclimation explains more of the batch than pigment acclimation**, with a
  slow (~100 h) timescale.
- **The culture-level squalene optimum only exists in variants where pigments feed back
  on growth** (v0, v1, g-none, g-coupled: optimum at 100–460x Sqs, ~1e9 cells/mL);
  fixed-pigment variants (B, C) saturate instead. → The minimal model must keep the
  FPP/GGPP competition and a pigment → light-harvesting → growth link; pigment contents
  can't simply be parameters if the trade-off question is to be answered (contrast with
  the batch fit, which preferred constant pigments).
- **Enzyme cost is missing in thread g** and would bound the optimum more tightly.
- **MCA picks and culture productivity agree on the first 2–4 overexpressions, not on
  when to stop**: later MCA picks reduce productivity at the productive high densities.
  → Optimise strains against the culture objective directly (or use MCA only for the
  ranking and the productivity scan for the stopping point).
- **Best realistic thread-g strain: Sqs + CrtE (+ MEP chain, DXS)**, 0.4–0.65 mg/L/h,
  about half the free squalene synthase optimum.

## Minimal model v3 (ADR 0005, report 06)

- **The gravimetric dry weight was the inconsistent number**: it contains only 22 %
  protein. Dry weight from the measured protein (464 g L^-1) is consistent with 4 fL and
  the cells/OD calibration.
- **A closed-form PETC beats the polynomial** (R² 0.998 / 0.979 vs. 0.983 / 0.939, 8
  parameters): one tanh light response plus J_ATP as a function of J_3PGA.
- **Simulating the metabolism directly is cheap enough** (steady state ~2 s, strain
  optimisation 3 s) and removes the surrogate's enzyme-input limitation.
- **Strongly growth-coupled pigment synthesis (exponent 5–9) fits the batch best and
  explains the initial growth rise** — a mechanism the constant-pigment variants lacked;
  untested against pigment data.
- **Enzyme-cost formulation is immaterial at realistic overexpression** (< 0.3 % in
  productivity): pigment adjustment compensates. It would matter only at > 10^3× levels.
- **Capacity acclimation changes the productivity answer by 3×** and moves the best
  density to practical values; without it v3 says "as dense as possible".
