# MCA-based strain optimisation (Fig 12) for all model variants

Script: `scripts/strain_optimisation.py` (→ `strain-optimisation.json`), figure:
`scripts/plot_strain_optimisation.py`.

## Procedure

As in the notebook / paper Fig 12, with **each enzyme overexpressed at most once**:

1. Steady growth at fixed cell density 1e7 cells mL^-1, light 100 µmol m^-2 s^-1, sample
   depth 1 cm (the notebook's conditions), Δshc strain.
2. Normalised flux response coefficients of squalene synthesis and of biomass synthesis
   to every enzyme parameter. Score = sum of the two.
3. Overexpress the enzyme with the highest score; repeat until every enzyme was used once
   (max 10 rounds). Round 1 answers "which *single* overexpression is best".

Overexpression is a fold change of the enzyme: v0/v1 scale the enzyme's synthesis rate,
thread g its level, by `1 + 1.331057 kDa / (2 E0)` (E0 = legacy proteomics
concentration). For v0 this is exactly the notebook's additive
`ksynth + 1.331057 kDa / (5.5 h)`. Thread g lumps DXR…IspH into one "MEP" enzyme, which
gets the geometric mean of their folds (in reality that's several genes). ChlP/ChlG are
one lump ("Chl").

Thread g uses `models/g_coupled.py`: the explicit metabolism ODE coupled to the growth
layer at fixed density, because the metabolism surrogate only has E_CrtE as an input.
Parameters come from each variant's batch fit.

## Results

![strain optimisation](figures/strain-optimisation.png)

Top: score per enzyme and round (red boxes = chosen). Bottom: intracellular squalene
(blue, g L^-1) and growth rate (orange, h^-1) after n overexpressions.

| variant | single best | order | SQ start [g/L] | after 1 | end | μ start [1/h] | μ end |
|:--|:--|:--|--:|--:|--:|--:|--:|
| v0 (legacy) | IspG | IspG, DXS, IspD, Sqs, ChlG, ChlP, Idi, IspF, IspH, ThiG | 0.196 | 0.203 | 0.773 | 0.060 | 0.078 |
| v1 (bugs fixed, not recalibrated) | DXS | DXS, Sqs, CrtB, Idi, IspD, IspF, IspE, IspH, ThiG, DXR | 0.044 | 0.108 | 0.367 | 0.011 | 0.013 |
| g, unregulated pigments | Sqs | Sqs, CrtE, MEP, DXS, Chl (last round failed) | 1.08 | 5.10 | 3.72 | 0.050 | 0.062 |
| g, A: growth-coupled pigment enzymes | Sqs | Sqs, CrtE, MEP, DXS, CrtB, Chl | 0.79 | 3.74 | 4.31 | 0.087 | 0.068 |
| g, B: photoacclimation (a ≈ 0) | Sqs | Sqs, CrtE, MEP, DXS, CrtB, Chl | 0.78 | 2.70 | 4.81 | 0.095 | 0.069 |
| g, C: capacity acclimation | Sqs | Sqs, CrtE, MEP, DXS, CrtB, Chl | 0.72 | 2.25 | 5.15 | 0.121 | 0.120 |

## Interpretation

- **v0 reproduces the notebook** (same order, start 0.2 g/L): a check of the procedure
  and of the fold-change formulation.
- **The recommendation depends on the model.** v0 picks IspG first, v1 picks DXS (fixing
  the DXS feedback bug B4 and the µM/mM Km defaults B6 moves control to DXS), and all
  thread-g variants pick **Sqs, then CrtE, then the MEP chain, then DXS**. The thread-g
  ranking is robust across the four pigment/capacity variants.
- In thread g, Sqs and CrtE carry most of the control because the reference state was
  built with Sqs saturated by FPP (Km 2.8 µM, legacy value) and CrtE feeding both FPP and
  GGPP; DXS only gains control once Sqs and CrtE are raised (round 4, score 0.2–2).
- The last rounds (CrtB, Chl) have scores ≤ 0 in most variants and reduce squalene: the
  procedure keeps going because each enzyme is used once; the useful stopping point is
  after 4 overexpressions (Sqs, CrtE, MEP, DXS) in thread g and after ~4 in v0.
- Absolute squalene levels aren't comparable between threads: thread g's reference was
  built with the sqs strain's measured squalene content (2.9 mM ≈ 1.2 g/L), the legacy
  model has no such calibration. v1 isn't recalibrated yet (it grows ~5x slower here).

## Caveats

- The lumped MEP chain in thread g stands for six genes; its overexpression isn't a
  single-gene intervention.
- Fold changes are tied to the legacy (2 fL) proteomics concentrations, as in the paper.
- Responses are steady-state, at one density/light; the batch dynamics (Fig 11) aren't
  part of this.
