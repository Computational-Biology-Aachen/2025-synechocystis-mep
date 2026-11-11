# Re-check of the mxlpy problems with mxlpy 0.69 (main, 3b606f09)

Dependency: `mxlpy` pinned to the git commit 3b606f09 (0.69.0, unreleased) in
`[tool.uv.sources]`. Scratch scripts for all checks: session scratchpad (not in repo).
**Speed evaluation of parallel scans is unfinished** (paused, see "Open").

## Issues filed earlier

| problem (where documented) | 0.68 | 0.69 |
| --- | --- | --- |
| #201 Native can't lower model data (ADR 0003) | `KeyError` | fixed: pfennig2024 with `add_data` runs natively (0.3 s), rel. diff to Sundials 1.2e-7 |
| #202 generated fns with `math.exp` unparseable (ADR 0003, `exprfn`) | needed `synmep.exprfn` re-homing | fixed: `add_derived`/`add_reaction` with sympy exp/log lower natively |
| #203 quadratic lowering of large fns (ADR 0003, report 01) | balanced-tree workaround | fixed, see compile times below |
| #204 `get_unused_parameters` misses derived stoichiometries | removing them broke pfennig2024 | fixed: 7 unused, `bHo` kept, model simulates after removal |

## First native call (lower + compile + one steady state)

Same machine, each in a fresh process. "old exprfn" = the repo's workaround module,
"compile_fn" = mxlpy's own `compile_fn` (what `synmep.exprfn` now wraps). Steady-state
values identical to ≤ 1e-8 relative in all cases.

| model | 0.68, old exprfn | 0.69, old exprfn | 0.69, compile_fn |
| --- | --: | --: | --: |
| single-legacy | 5.0 s | 4.0 s | 2.2 s |
| single-fixed | 4.7 s | 2.9 s | 1.9 s |
| double-coupled | 14.8 s | 2.6 s | 1.0 s |
| minimal-m0-cap | 1.1 s | 0.9 s | 0.5 s |

(single-legacy was 9.2 s in report 01 on 0.68; the difference to 5.0 s here is machine load /
caching, the relative comparison is what counts.)

## Other problems

- **Oscillation detector false positives** (ADR 0004): no false positive at default
  parameters for the PETC and minimal-m0 with the default detector. A scan-level recount
  (3648/4000 under 0.68) wasn't rerun; `NATIVE_RUNNER` keeps detection off.
- **Native steady-state search on single-fixed** (report 02, open question 5): converges at
  default parameters under 0.68 and 0.69 in this check; the robustness rerun
  (perturbed parameters) wasn't done yet.
- **Parallel native scans hung** (ADR 0004): no longer hang; results identical to serial
  (max rel. diff 0) for minimal-m0 (200 points), `double` (200) and the PETC (400). But they
  are **slower than serial** for these cheap models:

  | scan | serial | parallel (22 workers) |
  | --- | --: | --: |
  | minimal-m0, 200 pts | 4.2 s | 26.7 s |
  | double-coupled, 200 pts | 28.0 s | 35.3 s |
  | PETC, 400 pts | 6.6 s | 9.0 s |

## Findings on parallel overhead (partial)

1. **Pickling functions from `compile_fn` was slow**: they have `__module__ = None`, so
   cloudpickle's `_whichmodule` searches all of `sys.modules` for each one (~10 ms per
   function; minimal-m0 has 110 → 1.2–1.5 s per model pickle, and the pickle is repeated per
   task). Setting `fn.__module__` (done in `synmep.exprfn.expr_fn`) brings it to 0.02 s.
   Worth fixing in mxlpy's `compile_fn` (not filed yet).
2. With that fixed, a 44-point minimal-m0 scan still takes ~45 s in parallel vs ~1 s serial.
   Without the shipped precompiled worker state (`parallel._WORKER_STATE` cleared) it
   takes ~18 s: the ~70 KB state seems to make loky's sequential worker spawn block per
   worker (payload above the 64 KB pipe buffer?). Unconfirmed.
3. The remaining ~18 s pool start-up isn't explained yet (worker start + imports?).

## Native's RHS check gives false positives (found later)

Native compares its compiled right-hand side with the Python model at three perturbed
states, rtol 1e-8. In the `double` metabolism ODE the ATP balance (large, nearly
cancelling terms) differs by 2–4e-8 from rounding alone, so a correct model raised
`NativeLoweringError: Compiled model disagrees with the model at variables [1]` during a
density scan. Since 0.69 the check reruns on every compile-cache hit with the new start
state, so it fires far more often than under 0.68. Workaround: `synmep.runners` sets
`_CHECK_RTOL = 1e-6`. Candidate upstream issue: scale the tolerance with the magnitude
of the terms, or treat a mismatch on a cache hit as a miss.

## Open

- Bare pool start-up timing (trivial task, with/without heavy imports) to split (2)/(3).
- Then decide whether to report (1)–(2) upstream, and whether any script should use
  parallel scans. Until then all scans stay serial.
- Robustness rerun (report 02) and scan-level oscillation-detector recount.
