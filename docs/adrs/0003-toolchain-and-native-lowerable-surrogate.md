# ADR 0003: uv + mxlpy Toolchain and a Native-Lowerable Surrogate

**Status:** Implemented (v0)
**Scope:** `pyproject.toml`, `src/synmep/surrogates/petc_legacy.py`, `src/synmep/exprfn.py`,
`src/synmep/models/v0_faithful.py`

---

## 1. Context

The legacy model ran on modelbase + assimulo via pixi. Its PETC surrogate was a pickled
sklearn `Pipeline(PolynomialFeatures(5) → MultiOutputRegressor(LinearRegression))`,
fed by a pandas/scipy pipeline (self-shading, Simpson integration over 301
wavelengths, pigment association) that ran on **every RHS call** (~10 ms per call).
mxlpy ≥ 0.67 provides Sundials CVODE (`mxlpy[sundials]`) and a native runner
(`mxlpy[native]`) that compiles the model, but only if the model can be lowered to
sympy.

## 2. Decision

1. **Toolchain:** uv only, `mxlpy[sundials,native]>=0.68` from PyPI. `exclude-newer = "3
   days"` stays, with `exclude-newer-package` exempting mxlpy, mxlpy-native and
   sundials4py. `mxlmodels` (needed only to retrain the surrogate) is pinned to a git
   revision as the optional extra `surrogate`, because the PyPI release (1.4.0) predates
   the pfennig2024 rewrite. Python ≥ 3.13 (required by mxlmodels). The uv pre-commit
   hook is bumped to 0.13.6; 0.7.2 couldn't parse relative `exclude-newer` durations
   and silently re-locked without any cooldown.
2. **Surrogate as explicit expressions:** the polynomial coefficients are extracted once
   from the pickle (`data/reference/surrogate_poly.npz`) and expanded into sympy
   polynomials. sklearn is no longer a dependency.
3. **Spectral integrals unrolled:** spectra, absorption coefficients and composite
   Simpson weights are baked in as numeric constants. Each wavelength integral becomes
   an explicit 301-term sum of `w_λ a(λ) I0(λ) (1 − e^{−k_λ d}) / (k_λ d)`, with `k_λ`
   depending on the Chl/carotenoid/cell-density states. Full 1 nm resolution is kept;
   coarsening is only considered if compile time becomes a problem (ADR 0001).
4. **Model data (`add_data`) is not used** in any variant that should run natively.
5. **Workaround `synmep.exprfn.expr_fn`:** sympy expressions are compiled with
   `mxlpy.meta.sympy_tools.expr_to_fn` and re-homed into a real module that imports
   `math`/`numpy`/`scipy`.

## 3. Rationale

- Explicit expressions are exact (polynomial vs sklearn: rtol 1e-9; v0 RHS/fluxes vs
  legacy at 15 checkpoints: rtol 1e-7), introspectable, and lowerable.
- Two mxlpy 0.68 limitations were found while porting (reported upstream, not patched
  here, per ADR 0001):
  - **Native can't lower model data**
    ([MxlPy#201](https://github.com/Computational-Biology-Aachen/MxlPy/issues/201)). `add_data` series used in derived functions give
    `NativeLoweringError: KeyError('<data name>')`. This affects pfennig2024 itself.
  - **Generated functions with `math.*` calls can't be re-parsed**
    ([MxlPy#202](https://github.com/Computational-Biology-Aachen/MxlPy/issues/202)). `expr_to_fn`
    output calls e.g. `math.exp`, but `source_tools._handle_call` resolves `math`
    through `inspect.getmodule(fn)`, which a generated function doesn't have
    (`KeyError: 'math'`). The function then can't be converted back to sympy, and the
    native runner can't lower any model built from such expressions.
    `add_derived(name, sympy_expr)` with `exp` is affected.

## 4. Consequences

- First native call per process lowers and compiles the model: ~9 s for v0. This was
  64 s before splitting the PETC block into small intermediates and emitting long sums
  as balanced trees. The time is spent in mxlpy re-parsing functions into sympy
  (quadratic `Add` rebuilding, `Max`/`Min` simplification over huge arguments), not in
  the backend ([MxlPy#203](https://github.com/Computational-Biology-Aachen/MxlPy/issues/203)).
  After that a steady state takes ~0.4 s. loky workers in parallel scans each compile once, so short scans run faster
  serially (see `docs/reports/benchmark-v0.md`).
- For the retrained pfennig surrogate (stage 2), pfennig2024's spectral derived
  quantities depend only on parameters, so they are replaced by parameters
  (absorbed light per complex, OCP light). This both makes pfennig2024 lowerable and
  equals the surrogate's input space. In mxl-models pfennig2024 the pigment contents
  enter only through absorption, so the new surrogate needs 4 inputs, not 7.
- If mxlpy gains data lowering or fixes the generated-function lookup, `exprfn` and
  the baked-in constants can be revisited.
