# `full` model: pfennig2024 without surrogates (checkpoint B, ADR 0006 stage 4)

Code: `synmep.models.full` (registry name `full-k3`), `synmep.petc.pfennig`
(`get_petc_model(clamp_pga=False)`), `synmep.models.metabolism` (`biomass=` option).
Tests: `tests/test_full.py`.

## What was built (ADR 0006 §2.5)

- K = 3 Gauss–Legendre depth layers, each with its local light and a complete copy of
  pfennig2024 (17 states, 29 derived values, 6 QSS surrogates, 23 reactions), absorbed
  light per complex and OCP light from the light model, 3PGA free.
- Per layer: pfennig2024's pools (3PGA, ATP/ADP, NADPH/NADP+, …), FPP and GGPP
  (relaxation < 1 min: 0.1–0.3 min, Jacobian of the metabolism). Shared: MEPi (3 min),
  IPP (12 min), Chl, carotenoids, squalene, N.
- The metabolism's reactions run per layer on their layer's pools (mM ↔ mmol mol(Chl)^-1
  with the current Chl, min ↔ s) and act on shared pools with the layer weight;
  mu = sum_i w_i v_bio,i / rho_DW; model time in seconds.
- 62 states, 152 parameters, 113 reactions; builds in 3.4 s, first native call
  (compile + 10 min simulation) 1.2 s, steady state at 1e7 cells mL^-1 a few seconds.

## Changes to pfennig2024 and the metabolism

1. **Lumped ATP consumption removed**; the metabolism's demand (biomass, MEP chain,
   maintenance) replaces it.
2. **Metabolism re-referenced** to pfennig2024's ATP and NADPH levels (0.55 and 0.32 mM
   instead of the assumed 2 and 0.5 mM), so its maximal rates and Michaelis constants
   are consistent with the pools they now share.
3. **Biomass ATP demand 50 instead of 0.774 mmol gDW^-1** (decision needed, below).
4. Soluble totals stay per Chl as in pfennig2024. Per-cell totals (option
   `soluble_per_cell`) without converting pfennig2024's per-Chl Michaelis constants
   create a positive feedback (more Chl → less NADP per Chl → slower CBB → slower growth
   → more Chl): Chl ran away to 166 mM.

## Findings

- **With the paper's biomass composition the coupled model collapses.** The metabolism
  then needs almost only NADPH (2.9 NADPH vs 0.77 ATP per gDW); NADPH runs out, ATP piles
  up (95–256 of 430 mmol mol(Chl)^-1), ATP synthase stalls, the proton gradient
  builds, photosynthetic control throttles electron transport, and growth falls to
  ~0 (3PGA → 0). More light, 3× CBB capacity or less respiration don't help. The
  surrogate families hid this: their NADPH supply was an assumed multiple of the 3PGA
  supply and a maintenance sink absorbed the ATP balance.
- **With a genome-scale-like ATP demand (20–50 mmol gDW^-1) it works without any fitting:**
  dilute culture (1e7 cells mL^-1, 100 µmol m^-2 s^-1) μ = 0.030–0.036 h^-1 (data
  0.04–0.06), Chl 24–30 mM (data 16–21), squalene 1.1–1.3 mM, steady state converged
  (relative RHS 1e-6).
- **Spectral ranking at 1e7 cells mL^-1** (μ, h^-1): orange 0.042 > red 0.037 > warm
  white 0.030, as pfennig2024 alone.
- **Fails at low light per cell:** dense cultures (≥ 4.5e8 cells mL^-1, so also the
  Fig 9 pre-culture and therefore the batch) and blue light. The darkest layer's CBB
  barely exceeds respiration (0.34 vs 0.31 mmol mol(Chl)^-1 s^-1), ATP is high and
  NADPH low, and the integration fails within ~3 h. pfennig2024's respiration draws
  3PGA at a near-constant rate regardless of light, and the metabolism's sinks don't
  stop at low supply.

## Open (decisions at checkpoint B)

1. Biomass ATP demand: adopt a literature value (to be sourced; GEMs ~40–60 mmol gDW^-1)
   for `full` only, or for every family (then `double`/`minimal` lose their
   maintenance-sink buffer and must be refitted)?
2. Low-light failure: candidates are making respiration and the maintenance sinks
   supply-dependent, adding a fitted PETC/CBB efficiency (as the surrogate families
   have), or stiff-solver settings. It needs work before the batch (goal 1) runs.
