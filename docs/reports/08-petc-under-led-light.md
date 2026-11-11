# PETC forms under LED light (checkpoint A, ADR 0006 stage 3)

Scripts: `train_petc_led.py` (→ `data/surrogates/petc_training_led.csv`),
`fit_petc_balance.py` (→ `petc_balance.json`, `petc_poly_led.json`,
`docs/reports/petc-balance.json`), `petc_paired_spectra.py` (→ `petc_paired.csv`).
Code: `synmep.petc.balance` (m0, m1, m2), `synmep.light` (lamps, LEDs).

## Data

- **White:** the 3939 pfennig2024 steady states of the `double` surrogate training
  (warm white light, culture states through the self-shading light model).
- **LED:** 3838 more (of 4000; 74 s): same culture-state sampling plus a spectrum per
  sample (a quarter one of the six lamps, the rest Dirichlet(0.5) mixtures of blue
  450 nm, orange 620 nm, red 680 nm and warm white) and a varied PS1:PS2 ratio (2–12).
- **Paired:** 300 culture states, each under 8 sources at the same incident photon flux
  (warm/cool white, solar, blue, orange, red, orange+red, blue+red). This tests what
  goal 2 needs: the *ranking of spectra* for a given culture.

## Forms (`synmep.petc.balance`)

| form | J_3PGA | parameters |
| --- | --- | --- |
| m0 | J_max tanh((a_PS1 + w2 a_PS2 + w_pbs a_PBS)/I_k) − r | 5 + 3 (ATP) |
| m1 | J_max tanh(smin(e1, e2)/I_k) − r; e1, e2 with PBS share phi to PS2, OCP quenching | 8 + 4 |
| m2 | mini electron transport (PQ, P700/PC, Fd; saturating CBB sink), steady-state J | 9 + 4 |
| poly5 | degree-5 polynomial in the 4 log inputs (the `double` surrogate form) | 126 × 2 |

## Results

Hold-out R² (fitted to the union of both training sets, 6223 train / 1554 test):

| form | white J_3PGA / J_ATP | LED | imbalanced (one LED ≥ 80 %, n = 75) |
| --- | --- | --- | --- |
| m0 | 0.998 / 0.980 | 0.998 / 0.980 | 0.998 / 0.981 |
| m1 | 0.998 / 0.980 | 0.998 / 0.980 | 0.998 / 0.978 |
| m2 | **0.9994 / 0.991** | **0.9995 / 0.992** | **0.9994 / 0.991** |
| poly5 (LED-trained) | 0.980 / 0.938 | 0.980 / 0.939 | 0.983 / 0.946 |

R² over all samples is dominated by the light intensity (3–1500 µmol m^-2 s^-1), so the
paired test is the sharper one. J_3PGA relative to warm white at equal photon flux,
median over 162 states with clearly positive supply:

| source | pfennig2024 | m0 | m1 | m2 | poly5 LED | poly5 white (current `double`) |
| --- | --: | --: | --: | --: | --: | --: |
| blue | 0.71 | 0.71 | 0.71 | 0.74 | 0.74 | **0.17** |
| orange | 1.47 | 1.37 | 1.37 | 1.46 | 1.30 | 1.21 |
| red | 1.00 | 1.00 | 1.00 | 1.00 | 0.95 | **6.37** |
| orange + red | 1.29 | 1.21 | 1.21 | 1.26 | 1.20 | 2.96 |
| blue + red | 0.91 | 0.89 | 0.89 | 0.93 | 0.86 | 0.64 |
| cool white | 0.72 | 0.73 | 0.72 | 0.71 | 0.77 | 0.62 |
| solar | 0.89 | 0.90 | 0.89 | 0.88 | 0.88 | 0.74 |
| **median abs. ratio error (p90)** | – | 0.030 (0.10) | 0.027 (0.10) | **0.014 (0.06)** | 0.10 (0.33) | 1.35 (27) |

Findings:
1. **pfennig2024 is nearly additive in absorbed light.** Orange + red (1.29) is barely
   above the mean of orange and red (1.24). In the light model the PS1:PS2 split of
   Chl absorption is fixed by `ps_ratio`, so a spectrum only shifts Chl against PBS
   absorption, and pfennig2024 responds to that almost linearly. The additive m0
   therefore reproduces the ranking well; m1's balance term adds nothing.
2. **m2 is the most accurate form** (half the error of m0, captures the orange
   advantage: 1.46 vs 1.47).
3. **The current `double` PETC surrogate is unusable off white light** (red 6.4×,
   blue 0.17× warm white): the white-light training covers a thin manifold of the
   4 inputs, and the polynomial extrapolates wildly off it. Retrained on the LED data
   it's usable but still the least accurate form (0.10 median error).
4. **The ranking itself** (pfennig2024): orange > orange+red > red ≈ warm white >
   blue+red > solar > cool white > blue, i.e. phycobilisome light is used best, blue
   worst. The magnitudes (orange +47 %, blue −29 % at equal photons) are per culture
   state, before self-shading and acclimation.

## Calibration problem found on the way: phycobilisome content

`minimal-m0` grows faster in blue LED light than in warm white (0.119 vs 0.105 h^-1 at
the MCA condition), against the pfennig2024 ranking. Not a PETC problem (its old m0
parameters rank blue at 0.67): its batch fit set `pbs_scale` = 0.44, i.e. little PBS,
which removes exactly the light that white and orange have over blue. The in vivo
spectra (`synmep.pigments`) give **PC:Chl ≈ 3.8–5.6 mg mg^-1** (legacy table: 7),
while the models have 1.5 at `pbs_scale` = 1. So the data imply `pbs_scale` ≈ 2.5–3.7;
`minimal-m0` (0.44–0.54) and `double-coupled` (9.1) are both far off, because the batch
data only pin the product petc_scale × pbs_scale. The PC:Chl ratio should be a fit
target in the shared fit.

## Decisions needed (checkpoint A)

See the summary to the user; recorded in ADR 0006 once decided.
