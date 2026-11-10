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
