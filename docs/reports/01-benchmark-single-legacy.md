# Benchmark: single-legacy with Sundials vs native runner

**Model:** `synmep.models.single_legacy` (faithful port, 39 variables, 85 reactions,
195 parameters). **Script:** `scripts/benchmark_single_legacy.py` → `benchmark-single-legacy.json`.
**Machine:** Intel Core Ultra 7 165H, 22 logical cores, Linux. mxlpy 0.68.0,
mxlpy-native 0.4.0 (llvm backend), sundials4py 7.9.0.

| Column         | Meaning                                                              |
| -------------- | -------------------------------------------------------------------- |
| legacy         | modelbase + assimulo CVODE (old pixi env, `data/reference/timings.json`) |
| sundials       | mxlpy default integrator (`Sundials`, CVODE BDF via sundials4py)     |
| native         | mxlpy `Native(method="CVODE_BDF")`, compiled model + C solver        |

All workloads mirror the legacy reference runs; tolerances are the integrators'
defaults (1e-8).

## Results

| Workload                                      | legacy  | sundials | native            | native vs sundials | native vs legacy |
| --------------------------------------------- | ------- | -------- | ----------------- | ------------------ | ---------------- |
| one steady state                              | 20.8 s  | 4.10 s   | **0.41 s**        | 10×                | 51×              |
| Fig 9: 96 h batch time course (385 points)    | 13.7 s  | 2.63 s   | **0.027 s**       | 98×                | 512×             |
| Fig 10: 24 steady states (scan)               | 89.0 s¹ | 22.7 s¹  | **24.0 s** serial / 33.4 s¹ | 0.9× / 0.7× | 3.7×        |
| Fig 11: 14 induction protocols (14 d each)    | 64.7 s¹ | 14.3 s¹  | **9.5 s** serial / 28.4 s¹  | 1.5× / 0.5× | 6.8×        |
| one-off: model construction                   | –       | 2.0 s    | 2.0 s             |                    |                  |
| one-off: first native call (lower + compile)  | –       | –        | 9.2 s             |                    |                  |

¹ parallel over 22 cores (pebble for legacy, loky for mxlpy).

Per simulation, native is **10× (steady state) to ~100× (time course)** faster than
Sundials and 50–500× faster than the legacy setup.

## Accuracy vs. legacy reference

| Workload | sundials                       | native                         |
| -------- | ------------------------------ | ------------------------------ |
| steady state, max rel. deviation over all variables | 7e-14 | 9e-7 |
| Fig 9, max deviation / max(|value|) (μ, density, SQ) | ≤ 1.5e-7 | ≤ 1.2e-7 |
| Fig 10, max rel. deviation of μ | 8e-12 | 8e-8 |
| Fig 11, max rel. deviation (density, SQ) | ≤ 1.8e-7 | ≤ 1.6e-7 |

## Interpretation

- **Per-simulation speedup is large; scan speedup is limited by compile cost.** The
  native runner compiles once per process (9.2 s for single-legacy). loky workers in parallel
  scans each compile again, so parallel native was *slower* than serial native for
  these short scans (24 and 14 runs). Serial native already matches or beats
  22-core parallel Sundials. Native in parallel pays off once a scan has many runs per
  worker (e.g. Monte Carlo, surrogate training).
- **Steady states gain less (10×) than time courses (~100×).** `simulate_to_steady_state`
  integrates in segments and checks convergence, with Python overhead between segments.
  The time course runs in one native call.
- **Compile time was 64 s before restructuring** ([MxlPy#203](https://github.com/Computational-Biology-Aachen/MxlPy/issues/203)). It is spent in mxlpy re-parsing every
  model function into sympy, not in the backend compile. Splitting the PETC block into
  small intermediates and emitting long sums as balanced trees cut it to ~9–11 s
  (ADR 0003). The remaining cost comes from the 5 × 301-term spectral sums and the
  4 × 792-term polynomials. A smaller surrogate (later families) will shrink it further. (mxlpy 0.69: 2.2 s, report 07.)

## Pitfalls found while benchmarking

- `Simulator.update_parameter` mutates the model it was given (expected behaviour;
  `scan.*` is the tool for repeated runs). The hand-written serial Fig 11 loop reused one
  model object and produced wrong results (SQ off by 2000×) until each run got its own
  model.
- Generated model functions must be readable from worker processes (their source is
  needed for lowering). They are now written to content-hashed files (`synmep.exprfn`).

## Notebook (Figs 9–12) with single-legacy + native

`main.ipynb` was ported to mxlpy and runs end to end with single-legacy and the native runner.
Fig 12 (iterative MCA) reproduces the legacy heatmap and the first eight overexpression
choices (IspG, DXS, IspD, Sqs, ChlG, ChlP, Idi, IspF). The legacy loop then stopped
because it selected the *unused* parameter `DXS_nonactive_ksynth`, which has no
proteomics entry. single-legacy drops unused parameters, so its loop continues for two more rounds
(IspH, ThiG). This is a stopping-criterion artefact, not a difference in dynamics.
