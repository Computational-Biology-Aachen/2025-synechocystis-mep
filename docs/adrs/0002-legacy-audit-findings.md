# ADR 0002: Audit Findings of the Legacy Model

**Status:** Accepted (findings); fixes scheduled for v1 (ADR 0001)
**Scope:** `model.py`, `kinetics.py`, `main.ipynb` (legacy, "v0")

---

## 1. Context

Before porting, the legacy model was audited for bugs, unit consistency and the reported
instability. Details of the unit derivations are in `docs/units.md`; the inventory of
magic numbers is in `docs/assumptions.md`. All numbers below come from the legacy model
in the old pixi env (dshc, steady growth mode, default notebook parameters, density =
row 24 of `growth_rates_by_density.csv`).

## 2. Decision

v0 reproduces all findings below unchanged (faithful port). v1 fixes them, one commit per
fix, each with a before/after note appended to §5 of this ADR.

### 2.1 Bugs (fix in v1)

| #  | Finding | Location | Effect |
| -- | ------- | -------- | ------ |
| B1 | **Surrogate extrapolation + clipping.** The degree-7 polynomial is evaluated far outside its training range; outputs fall non-monotonically and are clipped `rates<0 → 0` (scalar path only, not the vectorised path) | `get_influx_rate_estimations` | **Main cause of instability** (§3) |
| B2 | PETC rates converted with a *fixed* Chl (7.5 mg OD^-1 L^-1 → 39.79 mM), but the model's Chl state settles at 12.64 mM → supply is overestimated 3.15× and doesn't respond to Chl changes | `unit_conversion` | wrong magnitude, no Chl feedback |
| B3 | Biomass drain = μ·stoich/60 with stoich in mmol gDW^-1: missing ρ_DW [gDW L^-1], spurious /60 → drain ≈ 2·10^4× too small; cofactor pools are set by arbitrary `*_out_metabolism` sinks instead | `Biomass` reaction | growth is not resource-coupled |
| B4 | `DXS_na` computes `DXS − DXS_free` with DXS in µM and DXS_free in mM (already /1000) → IPP/DMAPP feedback on DXS monomerisation is effectively disabled | `kinetics.DXS_na` | regulation missing |
| B5 | `v_sqs` mass-action ratio inverted: `(A·B/P)/Keq` instead of `P/(A·B)/Keq`; diverges as SQ → 0 | `kinetics.v_sqs` | latent singularity (masked by Keq = 7.6e27) |
| B6 | Default Km/Ki/Kd applied after the µM→mM conversion → defaults are 100 mM / 1000 mM / 1 mM instead of µM | `get_enzyme_parameters` | many Km 1000× too high |
| B7 | `reversible_ordered_2_comp*` contains `kcat·Et/(kcat·Et·Keq)` → 0/0 at Et = 0 | `kinetics` | NaN for knockouts (CrtE, ChlP) |
| B8 | `unit_conversion` docstring says µM min^-1, the result is mM min^-1 | `unit_conversion` | documentation |
| B9 | Cell volume 2 fL in code vs. 4 fL used for the empirical growth-law slope (paper) | notebook / `growth_function` | inconsistent; enzyme concentrations from proteomics also scale with 1/V |
| B10 | min-kcat estimate uses `conc/T_d` instead of `conc·ln2/T_d`; T_d = 300 min vs. paper 330 min | `calculate_min_kcats` | factor 0.63 hidden in an unexplained ×2 |

### 2.2 Modelling issues (decide per variant)

- **Growth law pins Chl.** μ = v_Chl,cell / Chl_ref forces `Chl* = 1/(slope·V)` = 12.64 mM
  at every steady state. μ = v_Chl/[Chl] is *not* a fix: it makes `dChl/dt ≡ 0`.
- **"NADPH_sink_*" reactions** are reversible ATP ⇌ NADPH and ATP ⇌ Fdred mass action —
  energy interconversion with no mechanism; misnamed.
- **Implausible pool sizes**: ATP 11.5 mM, Fdred 7.3 mM, ThiDP 20 mM.
- **1410 parameters** for 85 reactions (cross product of enzyme × compound × {Km,Ki,Kd}).
- Missing Keq for IspF, IspG, ChlG, ThiG, Shc → 1e10.
- Phycobiliprotein content (PC/Chl 1.50, APC/Chl 0.68 mg mg^-1) differs from the Fuente
  defaults used to train the surrogate (7.0 / 2.0); its source is unknown.
- The pickled surrogate is a **degree-5** polynomial (792 terms per output, no input
  scaling), while the paper's Methods describe a seventh-degree polynomial.
- In the current mxl-models pfennig2024, pigment contents only enter via the absorbed
  light per complex, so the legacy surrogate's three `pigment_*` inputs are redundant
  given `complex_abs_*` and `light_ocp`.
- The surrogate outputs NADPH : 3PGA : Fdred in the fixed ratio 1 : 8.96 : 2.96 (= the
  biomass stoichiometry), so it has only **two** independent outputs (ATP, "biomass flux").

## 3. Stiffness and instability diagnosis

**Default parameters — stiff but stable.** The numerical Jacobian at steady state has all
eigenvalues ≤ 0. The fastest modes (−157, −105, −7.6 min^-1) are IPP/HMB_PP/NADPH and
GPP/FPP; the slowest (−7.7e-4 min^-1 = μ) are enzyme/pigment dilution modes (ThiDP, ChlP,
CrtE, SQ). Stiffness ratio ≈ 2.3·10^5, plus one zero mode (Cell_density, constant in steady
mode). BDF handles this; it is not the cause of failures.

**Perturbations — failures come from the surrogate.** 12 random log-uniform (±1 decade)
perturbations per parameter group, simulated to steady state:

| group perturbed      | success | median wall time [s] | μ range [min^-1]   |
| -------------------- | ------- | -------------------- | ------------------ |
| cell density         | 100 %   | 61                   | 8.8e-5 – 9.8e-4    |
| ksynth               | 92 %    | 74                   | 1.6e-4 – 5.8e-4    |
| kcat                 | 83 %    | 68                   | 1.6e-4 – 1.2e-3    |
| `*_out_metabolism`   | 83 %    | 67                   | 4.4e-5 – 1.1e-3    |
| light                | 67 %    | 86                   | 2.8e-5 – 1.1e-3    |

Failures are CVODE corrector-convergence failures, mostly early (t < 40 min). Probing the
surrogate directly (baseline: light 100, Chl 12.64 mM):

- light 200 → 400 → 800: 3PGA supply 8.3 → 1.4 → 0 (clipped); ATP 4.8 → 8.1 → 24.4
- Chl 6 → 3 → 1 mM: ATP 5.4 → 11.8 → 1870; NADPH/3PGA/Fd clipped to 0

Mechanism: anything that lowers Chl or raises light pushes the surrogate inputs (absorbed
light per Chl, pigment ratios per Chl) out of the training range. The degree-7 polynomial
then extrapolates wildly, and clipping makes the RHS non-smooth. Because Chl feeds back
into the surrogate inputs (self-shading and the per-Chl normalisation), this forms a loop
that the solver can't follow. This confirms the gut feeling that the growth/metabolism/PETC
coupling is to blame, specifically **coupling through an extrapolating surrogate**, made
worse by B2 (wrong Chl scaling).

Runtime: one RHS call costs ≈ 10 ms (sklearn + pandas + Simpson); one steady state
≈ 26 s at default parameters, 60–90 s when perturbed.

## 4. Rationale / Consequences

- v1 must keep the surrogate within its validity domain: clamp/saturate its inputs
  smoothly, or retrain it with a range covering the scenarios (light up to ≥ 400 µE,
  Chl down to the knock-down levels), with a monotone/bounded form. The surrogate's
  training domain must be recorded with the coefficients.
- Fixing B2 (dynamic Chl) and B3 (biomass drain) changes the qualitative coupling: supply
  will scale with Chl and growth will draw on cofactors. Figs 9–12 must be re-checked
  after each of these fixes.
- **Decided at checkpoint a — growth law (v1):** keep the empirical μ ∝ v_Chl law, with
  an explicit slope parameter and one consistent cell volume; the PETC supply conversion
  uses the dynamic Chl state (fixes B2). Document that Chl* is pinned. The cell volume value
  (2 vs. 4 fL) is still open.
- **Decided at checkpoint a — surrogate (v1):** **retrain** the PETC surrogate from the
  current mxl-models `pfennig2024_synechocystis` over a wider, documented input domain with
  a bounded form (fixes B1). The pfennig steady states are computed with the mxlpy
  **native** runner. The exact legacy polynomial is kept only for v0.

## 5. Fix log (v1)

Metrics from `scripts/v1_impact.py` (dshc): **ss** = steady state at light 100 and the
density of data row 24; **Fig 9** = pre-culture steady state at row 24, diluted to the
first data point, 80 µE, 13.5 mm, 96 h, RMSE vs. `growth_rates_by_density.csv`.

| Step | ss μ [h^-1] | ss Chl [mM] | ss SQ [mM] | Fig 9 RMSE μ [h^-1] | Fig 9 RMSE ln N |
| ---- | ----------- | ----------- | ---------- | ------------------- | --------------- |
| v0 (baseline) | 0.0462 | 12.64 | 0.403 | 0.0053 | 0.135 |
| B4 DXS monomerisation | 0.0317 | 12.64 | 0.299 | 0.0059 | 0.094 |
| B5 Sqs mass-action ratio | 0.0317 | 12.64 | 0.299 | 0.0059 | 0.094 |
| B7 ordered-mechanism 0/0 | 0.0317 | 12.64 | 0.299 | 0.0059 | 0.094 |
| B6 default Km/Ki/Kd in µM | 0.0222 | 12.64 | 0.202 | 0.0286 | 1.635 |
| B10 minimal kcats (ln2, 330 min) | 0.0184 | 12.64 | 0.212 | 0.0219 | 1.314 |
| units mM/min, B8 | 0.0184 | 12.64 | 0.212 | 0.0219 | 1.314 |

- **B4** — `DXS_na` now follows the paper: k_na (E_t − E_a), E_a = free (no IPP/DMAPP)
  *and* ThiDP-bound DXS, in consistent units. A unit-only fix (E_t − 1000·E_f) made
  `DXS_nonactive` vanish (IPP Ki default 1000 mM, B6), which switched off DXS synthesis
  through `DXS_synth_f` (synthesis ∝ DXS_nonactive / (DXS_nonactive + DXS_inh_str)) and
  growth collapsed. That synthesis activation by the monomer isn't described in the
  paper; it is kept as is (open question).
- **B5** — `v_sqs` uses P/(A·B)/Keq. No measurable change (Keq = 7.6e27), but the rate
  no longer diverges as SQ → 0 (e.g. a fresh Δshc culture without squalene).
- **B7** — `reversible_ordered_2_comp*` use (KmP + P(1 + A/KiA))/Keq instead of
  kcat·Et/(kcat·Et·Keq)·(...). Identical for Et > 0; defined for knock-outs (Et = 0).
- **B6** — missing Km/Ki/Kd default to 100 µM / 1000 µM / 1 µM (`synmep.data_v1`),
  instead of 100 mM / 1000 mM / 1 mM. 1191 of 1278 kinetic constants are defaults; 18 of
  them are used by the model: DXS Ki(IPP, DMAPP, PYR), Kd(ThiDP); CrtE Ki(DMAPP, GPP,
  FPP), Km(GGPP); ChlP Ki(GGPP, Chl-GGPP), Km(Chl-GGPP, Chl-phytol); ChlG Km(Chl-GGPP);
  IspD Km(ATP); IspG Km(Fd); IspH2 Km(DMAPP, NADPH); ThiG Km(ThiDP). The largest effect
  is product inhibition of Chl synthesis (Km(Chl-phytol) 100 → 0.1 mM vs. a 12.6 mM
  pool) and stronger IPP/DMAPP feedback on DXS. Because growth follows Chl synthesis,
  batch growth now overshoots after dilution (peak 0.12 h^-1, 6.8e9 cells mL^-1 at 96 h):
  **the legacy parameter set was effectively tuned around the buggy defaults**; v1 needs
  recalibration once all fixes are in.
- **B10** — minimal kcats use the flux content·ln2/T_d with T_d = 330 min (paper: 5.5 h)
  instead of content/300 min, i.e. all min-kcats × 0.63 (the unexplained ×2 safety
  factor is kept).
- **Units, B8** — one unit system (mM, min): enzymes are converted to mM at load, the
  hidden `Et / 1000` conversions in the rate laws are gone, enzyme synthesis uses an
  explicit yield (0.1, the legacy `sf · 1000` µM artefact), `DXS_inh_str` is in mM, and
  the `unit_conversion` constants are named parameters (`chl_per_od`, `M_chl`,
  `cells_per_od` → derived `chl_ref`). Behaviour-preserving: RHS identical to the
  previous step up to 7e-16 (enzyme derivatives scaled by 1e-3).
