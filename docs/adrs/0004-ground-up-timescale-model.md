# ADR 0004: Ground-Up Three-Timescale Model (thread "g")

**Status:** Accepted, in progress
**Scope:** `src/synmep/petc/`, `src/synmep/surrogates/poly.py`, `src/synmep/models/g_*.py`,
`scripts/train_*.py`

---

## 1. Context

Besides the v0 → v1 → v2 → v3 refactoring thread (ADR 0001), a second development
thread builds the model **from the ground up**. It uses values and ideas from the legacy
model where useful, but simpler rate kinetics, and is fitted to the experimental data.
The layers are separated by timescale:

| Layer       | Timescale | Treatment                                                       |
| ----------- | --------- | --------------------------------------------------------------- |
| PETC        | seconds   | QSS surrogate of pfennig2024 (retrained, native)                |
| metabolism  | minutes   | small ODE model with simple kinetics → trained QSS surrogate    |
| growth      | hours     | ODEs: cell density, pigments, products (and slow enzyme pools)  |

## 2. Decision

### 2.1 PETC layer
- pfennig2024 from mxl-models, made native-lowerable (spectral data → parameters:
  light absorbed per complex and OCP light) and steady-state capable (3PGA clamped to its
  default pool; its net production is an output). A mass-action biomass drain, as in the
  legacy training model, made protons accumulate without bound and was rejected.
- Outputs [mmol mol(Chl)^-1 s^-1]: `J_3PGA` (net 3PGA production) and `J_ATP` (ATP
  consumption = ATP available to the cell). pfennig2024 has no explicit NADPH or Fd
  demand, so reducing-equivalent supply is tied to `J_3PGA` by the paper's biomass ratio
  (NADPH : 3PGA = 2.924 : 6.934).
- Inputs: the 4 absorbed-light quantities. Pigment contents enter pfennig2024 only
  through absorption, so they aren't separate inputs (legacy: 7 inputs).
- Surrogate: polynomial (degree 3) in clamped, scaled log inputs (`surrogates/poly.py`).
  Outside the training domain predictions are constant, not extrapolated.

### 2.2 Metabolism layer (minutes)
States (mM): 3PGA (all sugar phosphates), ATP, NADPH, IPP (IPP + DMAPP), FPP, GGPP.
Simple saturable kinetics on lumped steps (DXS, lumped MEP chain DXR…IspH, CrtE as
FPP/GGPP synthase, Sqs, CrtB, ChlP+ChlG), DXS feedback inhibition by IPP/DMAPP, and
dilution of every pool by μ.

**Growth follows from the carbon and energy balance**, not from the legacy empirical
μ ∝ v_Chl law. A saturable biomass synthesis flux draws 3PGA, ATP and NADPH with the
paper's biomass stoichiometry, and μ = v_biomass / ρ_DW. Rationale: per gram dry weight
(3.88 pg DW cell^-1 and 55 g DW per g Chl from the repo's data, independent of the cell
volume) the PETC 3PGA supply at the Fig 9 median operating point supports
μ ≈ 0.033 h^-1. The measured batch growth rates are 0.01–0.05 h^-1, so the balance gets
the magnitude right without fitting. It also removes the Chl pinning of the legacy law
(ADR 0002).

Further choices:
- pfennig2024 has no explicit NADPH demand: NADPH supply is `r_NADPH · S_PGA`, with
  `r_NADPH` set so that at the reference state 20 % of the NADPH supply is oxidised.
- ATP supply is the PETC's `J_ATP` (CBB ATP use is already internal to pfennig2024).
- Capacity of biomass synthesis: μ_max = ln2 / 5.5 h (turbidostat, paper); the three
  cofactor Km of biomass synthesis are set so that μ = μ_ref at the reference pools.
- Maximal rates of all other steps follow from a reference state (μ_ref = 0.04 h^-1, the
  initial batch growth rate; measured Chl, carotenoid and squalene contents; pools at
  their Km), so the reference is an exact steady state by construction.
- Product inhibition at branch points (DXS by IPP/DMAPP — Ki 65 µM, Banerjee2013 — and
  MEP intermediates; FPPS by FPP; GGPPS by GGPP) keeps pools bounded when μ → 0.
- CrtE is one enzyme for the FPP and GGPP synthase steps (`E_CrtE`).

### 2.3 Growth layer (hours)
ODEs for cell density N, Chl, carotenoids and squalene (`dX/dt = v_X − μX`,
`dN/dt = μN`). Self-shading (legacy light model, unrolled spectral sums) turns N, Chl and
carotenoids into the PETC inputs. Enzymes are parameters (Fig 11 knock-down, Fig 12 MCA
act on them). The metabolism enters as a QSS surrogate trained on steady states of the
metabolism ODE model.

### 2.4 Units and constants
mM, min, cells mL^-1; cell volume **4 fL** (paper); enzyme concentrations from proteomics
are recomputed for 4 fL. Every constant is a named parameter with unit and source.

### 2.5 Fitting
Target: the batch growth data in the repo (`growth_rates_by_density.csv`: growth rate
vs. density over 96 h; Figs 9/10), plus the paper's doubling times as soft targets
(5.5 h WT turbidostat, ~22 h under CrtE knock-down). Few fitted parameters (e.g. maximal
biomass synthesis capacity, PETC supply scaling); everything else fixed from data or
literature.

### 2.6 Findings while building it
- **Data inconsistency:** 3.88 pg DW per cell at 4 fL is a dry-mass density of
  0.97 g mL^-1 (cells: ~0.3–0.4). The cells-per-OD calibration, the cell volume or the
  dry weights don't fully agree. Ratios that don't involve the volume (DW : Chl) look
  plausible.
- **Paper vs. code biomass stoichiometry differ** (see ADR 0002); this thread uses the
  paper's table.
- **Native's oscillation detector gives false positives** for monotonically settling
  trajectories of pfennig2024 and the metabolism model ("Oscillatory behaviour
  detected"). With detection off (`synmep.runners.NATIVE_RUNNER`) all 4000 PETC
  training points converge (with it: 3648) in half the time.
- **Parallel scans with the native runner hung** (loky workers idle, main process busy);
  serial scans are used (4000 PETC steady states in ~1 min).
- At very low supply, steady-state pigment contents per cell rise strongly, because
  pigment synthesis doesn't slow with growth (see `docs/reports/g-thread.md`).
- **Surrogates must be trained on physically realisable inputs.** Independently sampled
  inputs left real operating points (low PBS:PS1 at 4 fL; high ATP:PGA supply ratios)
  outside the joint training distribution, where per-input clamping can't help and the
  polynomials extrapolated to negative supply or spurious growth. The PETC training
  samples culture states through the light model (incl. `pbs_scale`); the metabolism
  training samples supplies from the PETC training data × `petc_scale`.
- **The metabolism surrogate predicts steady-state contents X* = v_X/μ, not fluxes.**
  The growth layer relaxes `dX/dt = μ (X* − X)`, identical at QSS, but bounded when
  μ → 0 (the flux form ran away to 255 mM carotenoids). In steady mode the relaxation
  uses a fast constant rate (same fixed point, converges).
- Native's steady-state search fails for very slowly relaxing states where Sundials
  converges; `runners.steady_state_scan` retries failures with Sundials.

### 2.7 Fit status
Two fitted parameters (`petc_scale` 5.48, `pbs_scale` 1.26): μ RMSE 0.013 h^-1, ln N
RMSE 0.23. The acclimation rise and the peak are reproduced; growth collapses after
~50 h because the PETC is evaluated at the depth-averaged light (see the report for the
proposed depth-resolved PETC and pigment regulation).

## 3. Rationale
- QSS at each timescale boundary is physically justified (s → min → h) and keeps the
  ODE system small and non-stiff for fitting.
- Training the metabolism surrogate on a mechanistic ODE model (rather than writing
  algebraic fluxes directly) keeps the kinetics inspectable and lets us check the QSS
  assumption.

## 4. Consequences
- The retrained PETC surrogate is shared with v1 (ADR 0002 §4).
- Results are compared against v1 and the experimental data in a report per stage.
- Open: whether lumped enzyme pools need their own slow dynamics (synthesis/dilution)
  for Figs 11/12, or whether parameters suffice.
