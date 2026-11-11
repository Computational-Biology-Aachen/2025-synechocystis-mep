# ADR 0005: Minimal Model Design (v3)

**Status:** Proposed (design; implementation pending)
**Scope:** future `src/synmep/models/v3_minimal.py`; supersedes the v3 row of ADR 0001

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
- Open before implementation: enzyme cost formulation (cofactor drain as in v0/v1 vs.
  proteome share), whether capacity acclimation is included, and the dry-weight /
  cell-volume inconsistency (affects the PETC supply scale).
