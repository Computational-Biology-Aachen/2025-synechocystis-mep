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

## Update after checkpoint B

Decisions (ADR 0006 §6): biomass ATP 50 mmol gDW^-1 for all families; supply-dependent
sinks for `full`.

- **Dim layers locked up energetically**: nothing used ATP without NADPH (pfennig2024's
  lumped ATP consumption removed, maintenance sink off), so ATP charged fully (ADP
  4e-8), ATP synthase stalled, lumen pH 3.4, cytoplasm pH 15, electron transport
  stopped, NADPH 1e-6, and the solver failed. Tighter tolerances (atol 1e-14) didn't
  help.
- **Fix: non-growth-associated maintenance (NGAM) per layer**, 0.13 mmol ATP gDW^-1 h^-1
  (lower bound from protein turnover, iJN678-based model, arXiv 2506.05916), saturating
  in ATP (K_M 0.1 mM). Steady states now converge up to 4.5e8 cells mL^-1 (μ 0.009
  h^-1, Chl 70 mM); 3e9 doesn't finish. Blue light works (dilute: blue 0.032, warm
  white 0.030, red 0.036, orange 0.041 h^-1, unfitted, pbs_scale 1 while the PC:Chl data
  imply ~3).
- **The batch runs** (Fig 9 protocol), unfitted μ 0.012–0.018 h^-1 vs 0.04–0.06 in the
  data (RMSE μ 0.021, ln N 1.12). But it takes 5–7 min: one long integration fails,
  short chunks work (evaluator: chunked integration with recursive halving).
- **Cause of the slowness: sustained fast oscillations in the darkest layer** (period
  ~4 s; PC, PSII states, FNR driving force), so the solver needs tiny steps over 96 h.
  Whether pfennig2024 oscillates in dim light on its own, or only coupled, isn't
  checked yet.
