# Units reference

Working reference for the unit system of the legacy model (`model.py`, `kinetics.py`,
called **v0** from here on) and the target system for v1+.
Derivations are checked numerically against a legacy steady state (dshc, steady growth,
default parameters, cell density = `growth_rates_by_density.csv` row 24).

## 1. Target unit system (v1+)

| Quantity                       | Unit                          |
| ------------------------------ | ----------------------------- |
| time                           | min                           |
| intracellular concentration    | mM (= mmol L(cell)^-1)        |
| enzymes                        | mM (converted at data load)   |
| reaction rates                 | mM min^-1                     |
| first-order constants, kcat, μ | min^-1                        |
| Km, Ki, Kd                     | mM                            |
| cell density                   | cells mL^-1                   |
| light (PPFD)                   | µmol(photons) m^-2 s^-1       |
| depth                          | m                             |

No `/1000` or `*60` inside rate laws. Every conversion factor is a named parameter
with a unit and a source (see `assumptions.md`).

## 2. Legacy (v0) unit system as implemented

v0 is **mixed**: metabolites in mM, enzymes in µM, with unit conversion hidden inside
rate laws.

| Quantity                                  | v0 unit                     | Where                              |
| ----------------------------------------- | --------------------------- | ---------------------------------- |
| metabolites (3PGA, ATP, MEP …, Chl, SQ)   | mM                          | implicit (rate laws `Et/1000`)     |
| enzymes (DXS, IspD …)                     | µM                          | proteomics → µM; `Et/1000` in laws |
| kcat                                      | min^-1                      | `Kcat [s^-1] * 60`                 |
| Km, Ki, Kd                                | mM                          | `[µM] / 1000`                      |
| Cell_density                              | cells nL^-1                 | `*1e6` → cells mL^-1 at surrogate  |
| time                                      | min                         |                                    |
| PETC surrogate outputs                    | mmol mol(Chl)^-1 s^-1       | pfennig2024 native unit            |
| `vX production` (after `unit_conversion`) | mM min^-1                   | see §3.1 (docstring says µM — wrong) |
| `Growth_rate`                             | min^-1                      | `growth_function`                  |
| phycocyanin / allophycocyanin parameters  | mg L(cell)^-1               | 53180 / 24076 mg L^-1 at 2 fL      |

## 3. Derivations

### 3.1 `unit_conversion` (surrogate rate → model rate)

```
v [mmol mol(Chl)^-1 s^-1]
 × Chl_cell [mmol cell^-1]   = 7.5 mg OD^-1 L^-1 / 893.509 g mol^-1 / (1.0548e8 cells mL^-1 OD^-1 × 1000)
                             = 7.958e-14 mmol cell^-1   (≈ 71 fg Chl cell^-1)
 / V_cell [L cell^-1]        = 2e-15            → 39.79 mM Chl
 × 1e-3 [mol mmol^-1]                           → mol(Chl) L^-1
 × 60 [s min^-1]
 = mM min^-1                 (factor 2.387 per mmol mol(Chl)^-1 s^-1)
```

**Finding:** the result is mM min^-1, not µM min^-1 as documented. And it uses a
*fixed* Chl of 39.79 mM, while the model's own Chl state settles at 12.64 mM (§3.2).
The per-Chl PETC fluxes are therefore scaled by a Chl content 3.15× higher than the
cell actually has.

### 3.2 Growth law

```
μ [min^-1] = 3.9566e13 [cell mmol^-1] × v_Chl [mM min^-1] × V_cell [L cell^-1]
           = v_Chl,cell / Chl_ref,cell,      Chl_ref,cell = 1/3.9566e13 = 2.527e-14 mmol cell^-1 (≈ 22.6 fg)
```

Combined with dilution (`dChl/dt = v_Chl − μ·Chl`), the steady state is forced to
`Chl* = 1 / (3.9566e13 · V_cell) = 12.637 mM` (at 2 fL), whatever the kinetics are. This
is confirmed numerically (legacy ss Chl = 12.636945 mM). The law acts as a proportional
controller on Chl. The paper (Results, rate correlation figure) describes the slope as
fitted assuming a **4 fL** cell volume, while the code uses 2 fL.

Note: the "obvious" alternative μ = v_Chl/[Chl] makes `dChl/dt ≡ 0` (Chl conserved, zero
eigenvalue), so it is **not** a valid replacement as-is.

### 3.3 Biomass drain

The `Biomass` reaction has rate μ and stoichiometry `-stoich/60`, with stoich in
mmol gDW^-1 (paper Tab. biomass_stoichiometry, from the Knoop/Kugler GSM). As written it
consumes μ × mmol gDW^-1 / 60, which is dimensionally not a concentration rate. Correct:

```
v_biomass,X [mM min^-1] = μ [min^-1] × stoich_X [mmol gDW^-1] × ρ_DW [gDW L(cell)^-1]
```

With ρ_DW ≈ 300–400 g L^-1 the drain is ≈ 2·10^4 × larger than in v0. In v0 the biomass
reaction consumes ~0.01 % of the 3PGA supply, and ~97 % leaves through the arbitrary
first-order `*_out_metabolism` sinks (k = 0.5 min^-1).

### 3.4 Self-shading / effective light

```
I0(λ)          = spectrum_rel(λ) [nm^-1, ∫=1] × PPFD [µmol m^-2 s^-1]
a_cell(λ)      = Σ_pigments abs_coef(λ) [m^2 mg(pig)^-1] × content [mg(pig) mg(Chl)^-1]   → m^2 mg(Chl)^-1
Chl_sample     = Chl [µM] × V_cell [L] × density [cells mL^-1] × 1e6 × 1e-3 × 893.509   → mg(Chl) m^-3
k(λ)           = a_cell(λ) × Chl_sample                                                 → m^-1
Ī(λ)           = I0(λ) (1 − e^{−k d}) / (k d)        (depth-average over [0, d])
complex_abs_X  = Σ_λ w_λ assoc_X,pig a_pig(λ) content_pig Ī(λ) × M_Chl                  → µmol(photons) mmol(Chl)^-1 s^-1
```

Dimensionally consistent. `w_λ` = Simpson weights on the 1 nm grid 400–700 nm.

### 3.5 Kinetic parameters

- `kcat [min^-1] = mean(Kcat [s^-1]) × 60`, then overwritten by `calculate_min_kcats`
  (every enzyme except ThiG/Idi): `kcat = 2 × (Σ C5 demand [µM] / T_d [min]) / E [µM]`.
  This uses `conc / T_d` instead of `conc · ln2 / T_d`, and the ×2 safety factor is
  unexplained.
- Km/Ki/Kd: `[µM]/1000 → mM`. Missing values default to Km = 100, Ki = 1000, Kd = 1 —
  but the defaults are applied *after* the /1000, so they are **mM**, i.e. 100 mM Km.
- Keq: dimension depends on reaction order (e.g. 1/mM for A+B→P) but is taken raw from
  `Equilibrium_constants.csv`. Missing (IspF, IspG, ChlG, ThiG, Shc) → default 1e10.

### 3.6 Enzyme synthesis

`ksynth = 2 E0 / (5.5 h × 60)` [µM min^-1]. The reaction (ping_pong_4, Km = 1 mM for all
four cofactors) uses stoichiometry `enzyme: 1e-4 × 1000 = 0.1` (µM) and
`cofactor: −1e-4 × cost` (mM), with cost in mol cofactor mol(enzyme)^-1. This is
self-consistent (0.1 µM enzyme = 1e-4 mM ↔ 1e-4·cost mM cofactor), but the effective
synthesis rate is 0.1·ksynth, not ksynth.

## 4. Steady-state magnitudes (legacy, default parameters)

| Pool                | value  | unit  | comment                                   |
| ------------------- | ------ | ----- | ----------------------------------------- |
| 3PGA                | 21.0   | mM    | stands for "all sugars"                   |
| ATP                 | 11.5   | mM    | ~3–10× literature                         |
| NADPH               | 3.3    | mM    |                                           |
| Fdred               | 7.3    | mM    | Fd total in vivo ≪ 1 mM                   |
| ThiDP               | 20.0   | mM    | cofactor at 20 mM is implausible          |
| MEcPP               | 2.6    | mM    |                                           |
| chlorophyll-phytol  | 12.64  | mM    | pinned by growth law (§3.2)               |
| Growth_rate         | 7.7e-4 | min^-1| T_d ≈ 15 h                                |
