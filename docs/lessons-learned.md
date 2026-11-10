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
