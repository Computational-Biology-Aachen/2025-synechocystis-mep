# Assumptions & magic numbers (legacy model)

Inventory of every hard-coded number in the legacy model and notebook, with what it
means and whether its origin is known. **Source "?"** means it isn't documented anywhere
in the repo or the paper. These must become named, sourced parameters in single-fixed (ADR 0001).

## Physiology / conversion constants

| Value              | Meaning                                | Unit                 | Source / status                                   |
| ------------------ | -------------------------------------- | -------------------- | ------------------------------------------------- |
| 2e-15              | cell volume                            | L cell^-1            | ? — paper rate-correlation figure assumes **4 fL** |
| 0.445              | protein content of dry weight          | g g^-1               | ?                                                 |
| 1.0548216418e8     | cells per OD750                        | cells mL^-1 OD^-1    | calibration curve (paper Methods), file not in repo |
| 7.5                | Chl a per OD                           | mg L^-1 OD^-1        | ?                                                 |
| 2                  | carotenoids per OD                     | mg L^-1 OD^-1        | ?                                                 |
| 0.5                | squalene per OD                        | mg L^-1 OD^-1        | ?                                                 |
| 893.509            | M(Chl a)                               | g mol^-1             | ChEBI 189419 (duplicated 4×, also in molar_masses.csv) |
| 536.8726 / 581.5565| M(carotenoids)                         | g mol^-1             | β-carotene (NIST) vs. ? — two different values used |
| 410.73             | M(squalene)                            | g mol^-1             | standard                                          |
| 300                | doubling time for min-kcat estimate    | min                  | ? (paper: 5.5 ± 0.5 h = 330 min)                  |
| 5.5 h              | enzyme doubling time for ksynth        | h                    | paper turbidostat doubling time                   |
| 3.9566524e13       | growth-law slope                       | cells mmol(Chl)^-1   | empirical fit (paper), derivation not in repo; see units §3.2 |
| 53180, 24076       | phycocyanin, allophycocyanin per cell  | mg L^-1 at 2 fL      | ? — gives PC/Chl = 1.50, APC/Chl = 0.68 mg/mg vs. 7.0 / 2.0 in Fuente defaults |
| 0.26               | photoactive fraction of carotenoids    | –                    | ?                                                 |
| 5.9                | PSI:PSII ratio                         | –                    | ? (photosystem_composition.csv says 5)            |
| 0.01 / 0.0135 / 0.02 | sample depth                         | m                    | 13.5 mm = half the MC-1000 tube diameter (paper Fig. 9 caption); others ? |
| 0.027 / 0.972      | G3P / PYR share of 3PGA pool           | –                    | ? (sum 0.999)                                     |
| 0.9375 / 0.0625    | IspH product split IPP / DMAPP         | –                    | literature ratio 15:1 (implied), as enzyme split  |
| 0.8 / 0.2          | active / non-active DXS initial split  | –                    | ?                                                 |

## Kinetics

| Value       | Meaning                                   | Status                                       |
| ----------- | ----------------------------------------- | -------------------------------------------- |
| Keq 1e10    | default for missing Keq (IspF, IspG, ChlG, ThiG, Shc) | ? — effectively irreversible           |
| Km 100      | default Km — applied *after* µM→mM, so 100 **mM** | ? — almost certainly meant 100 µM     |
| Ki 1000, Kd 1 | default Ki / Kd (mM after the same bug) | ?                                            |
| kcat 600    | default kcat (min^-1)                     | ? (overwritten by min-kcats for most enzymes)|
| ×2          | safety factor on min kcats                | ?                                            |
| ThiG kcat 100, Idi kcat 138 | fixed kcats               | Idi "from literature", ThiG "unknown"        |
| h = 1       | Hill coefficients                         | default                                      |
| synth Km 1  | Km of enzyme synthesis for each cofactor  | ? (mM)                                       |
| 1e-4        | enzyme-cost scaling `sf`                  | numerical convenience, self-consistent       |
| Kna_DXS 4e-6, DXS_na_kdeg 2e-5, DXS_inh_str 0.1/0.001/0 | DXS regulation | ? — changed per figure without explanation |
| 0.5         | `*_out_metabolism` first-order sinks      | ? — dominant sinks (units §3.3)              |
| 0.1         | ATP↔NADPH / ATP↔Fdred interconversion     | ?                                            |
| ThiGDXPKm 0.01 | overridden in notebook                 | ?                                            |

## Biomass

| Value                                 | Meaning                                       | Status                         |
| ------------------------------------- | --------------------------------------------- | ------------------------------ |
| 8.7629 / 33.0939 / 25.9760 / 78.4785  | NADPH / ATP / Fdred / 3PGA per gDW            | Kugler GSM (paper table), mmol gDW^-1 |
| /60                                   | applied to the biomass stoichiometry          | wrong — see units §3.3         |

## Notebook / experiments

| Value                       | Where   | Status                                         |
| --------------------------- | ------- | ---------------------------------------------- |
| light 100 (default), 80     | Fig 9/10| 80 µE = experimental condition (paper Fig. 9)  |
| ShcKcat=0, Shc_ksynth=0     | Fig 9   | Δshc strain, but strain="dshc" already removes Shc reaction |
| density index 24 / 5.456e7  | y0      | pre-culture density                            |
| CrtE3Kcat → 0.01            | Fig 11  | CrtE knockdown representation                  |
| 1.331057 × kDa, /(5.5·60)   | Fig 12  | overexpression amount from rhamnose-promoter proteomics |

## Structural assumptions

- PETC (pfennig2024) is at quasi-steady state relative to the terpenoid model (s vs. h).
- Chl a is the biomass proxy; growth is driven only by Chl synthesis (empirical
  correlation, paper Results).
- All pools, including enzymes and pigments, are diluted by growth.
- The 3PGA pool stands for all sugar phosphates; G3P/PYR are fixed fractions of it.
- One lumped "carotenoids" product (2 GGPP → 1 carotenoid).

## `minimal-m0` (ADR 0005)

| Value | Meaning | Unit | Source / status |
| --- | --- | --- | --- |
| 0.83 pg cell^-1 | protein per cell (median of 36 proteomics samples) | pg | measured (`Proteomics_samples_proteinconcentration.csv`, cells from OD) |
| 464 | dry weight per cell volume = 0.83 pg / 0.445 / 4 fL | g L^-1 | derived; replaces the gravimetric 970 g L^-1 (protein only 22 % of it) |
| f_i (1e-5 to 2.4e-3) | WT proteome mass fraction of each lumped enzyme | – | proteomics (E0 × kDa / protein on the same basis) |
| 0.45 | growth-independent proteome share φ_Q (proteome-cost variant) | – | Scott2010 (E. coli); assumption for Synechocystis |
| biomass composition | cost of 1 g overexpressed protein (drain-cost variant) | mmol g^-1 | assumption: as 1 g average biomass |
| closed-form PETC (8 parameters) | J_3PGA = J_max tanh(I/I_k) − r, J_ATP quadratic in J_3PGA | – | fitted to pfennig2024 steady states (`scripts/fit_petc_closed_form.py`) |
