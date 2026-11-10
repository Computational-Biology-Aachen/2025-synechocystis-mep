# ADR 0001: Rewrite Plan, Stages and Model Variants

**Status:** Accepted
**Scope:** whole repository (except `paper/`, which is read-only)

---

## 1. Context

The legacy model (`model.py`, `kinetics.py`, `main.ipynb`; built with modelbase and pixi)
couples three layers (see `tmp/model-scheme.png`):

1. a **surrogate** of the pfennig2024 PETC model (pickled sklearn polynomial) that
   turns light, pigments and self-shading into ATP/NADPH/Fdred/3PGA supply,
2. a **growth model** (cell density, dilution of all pools, μ ∝ Chl synthesis, enzyme
   synthesis with cofactor costs),
3. the **inner metabolism** (MEP pathway → IPP/DMAPP → GPP/FPP/GGPP → squalene,
   hopanoids, carotenoids, chlorophyll).

It is slow (~10 ms per RHS call, ~26 s per steady state), stiff, reported to be unstable
over wide parameter ranges, has unit inconsistencies (ADR 0002, `docs/units.md`), and has
1410 parameters.

## 2. Decision

Rewrite in stages, producing one model variant per stage, each a separate module in
`src/synmep/models/`:

| Variant          | Stage | What                                                                 |
| ---------------- | ----- | -------------------------------------------------------------------- |
| `v0_faithful`    | 1     | 1:1 mxlpy port, bugs included, regression-tested against legacy outputs |
| `v1_fixed`       | 2     | bug fixes (one commit each), single unit system, sourced parameters   |
| `v2a_timescale`  | 3     | fast pools (cofactors, MEP/prenyl intermediates) → steady-state polynomial surrogate; slow pools (enzymes, pigments, products, density) stay ODEs |
| `v2b_wholecell`  | 3     | whole cell at balanced growth → polynomial surrogate; growth model reduces to dN/dt = μN with self-shading |
| `v3_minimal`     | 4     | smallest model that still passes the acceptance criteria (§2.4)       |

### 2.1 Stage 0 — audit
Unit/assumption reference (`docs/units.md`, `docs/assumptions.md`), bug list and
stiffness diagnosis (ADR 0002), `docs/lessons-learned.md` started.

### 2.2 Stage 1 — reference + port
- Before deleting the pixi env, generate **reference outputs** from the legacy model
  (steady state, Fig 9 trajectory, Fig 10 scan, Fig 11/12 results) and the explicit
  coefficients of the PETC surrogate polynomial → `data/reference/`. Scripts live in
  `scripts/legacy/` (they only run in the old env; kept for provenance).
- Switch to uv with `mxlpy[sundials,native]>=0.68` from PyPI. mxlpy, mxlpy-native and
  sundials4py are exempt from the `exclude-newer` cooldown.
- Port v0 faithfully; benchmark mxlpy default sundials vs native runner (+ legacy
  modelbase/assimulo baseline) → `docs/reports/`.

### 2.3 Stages 2–4
v1 → v2a/v2b → v3 as in the table. v2a and v2b are both built on v1 and compared.

### 2.4 Acceptance criteria
- v0: matches legacy reference (steady state / trajectories, rtol ≈ 1e-6).
- v2a/v2b: surrogate hold-out error reported; Figs 9–12 compared against v1.
- v3: refitting allowed; Fig 9/10 RMSE vs. experimental growth data ≤ v1 + 10 %;
  Figs 11/12 qualitatively the same as v1 (shape of induction-timing curve; same top ~3
  MCA overexpression targets).
- Target data are the paper's Fig 9–12 outputs only.

### 2.5 Process
- Review checkpoints after the audit, v0 + benchmark, v1 and v2a/v2b.
- Conventional commits on `main`, committed frequently; commits only ever include our own
  paths (never `paper/` or `.gitmodules`).
- `paper/` is never edited, committed in, or pushed.

## 3. Rationale

- **Faithful port first** separates "did the port change behavior" from "did the fix
  change behavior". Each fix then gets an attributable before/after.
- **Separate variant modules** instead of flags keep each variant readable and let the
  minimal model shrink without dragging options along.
- **Two v2 splits**: the timescale split is physically justified (QSS of fast pools) and
  keeps transients (induction, batch). The whole-cell split is the most aggressive
  decoupling but assumes balanced growth. Which one is good enough is an empirical
  question, so both are built and compared.

## 4. Consequences

- The legacy files (`model.py`, `kinetics.py`, `parallel.py`) stay until v0 is verified,
  then are removed (recoverable from git history).
- The pixi env disappears after reference generation, so `scripts/legacy/` can then only
  be rerun by checking out an older commit.
- The growth-law choice for v1 is deliberately deferred to after the audit (see ADR 0002).
