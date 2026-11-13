# Architecture context

Entry point for *why* this repository is shaped the way it is. Each ADR covers one
decision in Status/Context/Decision/Rationale/Consequences form.

Model families and variant names (`single`, `double`, `minimal`, `full`) are defined in
ADR 0006. Reference documents (not ADRs, kept current):

- [`../units.md`](../units.md) — unit system, conversion derivations, steady-state magnitudes
- [`../assumptions.md`](../assumptions.md) — inventory of every magic number and its source
- [`../lessons-learned.md`](../lessons-learned.md) — running notes feeding the minimal model
- [`../reports/`](../reports/) — reports, numbered in reading order (01 benchmark of single-legacy … 10 comparison of all variants and selection)

## ADRs

→ [ADR 0001 — Rewrite plan, stages and model variants](0001-rewrite-plan-and-stages.md)
→ [ADR 0002 — Audit findings of the legacy model](0002-legacy-audit-findings.md)
→ [ADR 0003 — uv + mxlpy toolchain and a native-lowerable surrogate](0003-toolchain-and-native-lowerable-surrogate.md)
→ [ADR 0004 — Ground-up three-timescale model (`double` family)](0004-ground-up-timescale-model.md)
→ [ADR 0005 — Minimal model design (`minimal-m0`)](0005-minimal-model-design.md)
→ [ADR 0006 — Goals, model families and selection](0006-goals-families-and-selection.md)
