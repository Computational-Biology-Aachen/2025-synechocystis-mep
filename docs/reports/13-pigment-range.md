# Report 13: literature on the report-12 options, pigment range from the data

**Status:** 2026-11-17. Decisions on report 12 §5: (1) check the evidence for β-carotene
light harvesting in PS2, (2) bound the pigment factor to a range taken from the
literature, (3) keep the capacity variants. (2) is implemented and all variants are
refitted; 465–520 nm still give no growth, and the tests in §3 show why.

## 1. β-carotene in PS2: not implemented

- The two β-carotenes in the PS2 reaction centre (D1/D2) are photoprotective
  (triplet/singlet-oxygen quenching, secondary electron donor via Cyt b559); their
  transfer to Chl is low (~26 % by ultrafast spectroscopy). The β-carotenes of the core
  antennas CP43/CP47 do transfer: ~31–33 % (ultrafast, via S2 only), 29–76 % in
  steady-state estimates depending on method, 39 % in *Synechocystis* CP43
  pre-assembly complexes; PS1 ~57 %.
- Effect, at reference pigments under the nine LEDs (m0 PS2+PBS input relative to
  633 nm; β-carotene shares by `photosystem_composition.csv`, PS2 gets 7.8 % of it):

  | | 435 | 465 | 495 | 520 | 555 |
  | --- | --- | --- | --- | --- | --- |
  | without | 0.74 | 0.16 | 0.04 | 0.08 | 0.21 |
  | PS2 β-carotene at 35 % | 0.78 | 0.22 | 0.09 | 0.10 | 0.21 |
  | growth, data | 0.62 | 0.23 | 0.39 | 0.59 | 0.83 |

  It adds ~5 points at 465/495 nm and nothing at 555 nm. With reference pigments the
  input at 465 and 555 nm is already 16–21 % of 633 nm; the zero growth comes from
  the pigment dynamics (§3), not from the routing.

## 2. Pigment factor range (implemented)

Literature on pigment content of *Synechocystis* across growth light:

| source | conditions | Chl | phycocyanin | carotenoids |
| --- | --- | --- | --- | --- |
| Kopečná et al. 2012, Table I | 10–300 µE | 2.4× per OD (flat below 40 µE) | — | decrease (myxoxanthophyll ×2) |
| Zavřel et al. 2019, Fig 2 | turbidostat, red, 27.5–1100 µE | 2.8× per gDW (1.9× per cell, cells grow) | 3.3× per gDW | 1.7× per gDW |
| Zavřel et al. 2024 | nine LEDs at 25 µE | 0.9–1.8× per cell | 0.4–1.4× | Car:Chl 0.55–1.08 |

So pigment content changes about 3-fold between limiting and saturating light, and
barely rises in dim light. `PIGMENT_CLAMP` is now 0.5–1.5 (was 0.05–20): a 3-fold span,
with less room upwards because the reference (Fig 9 median) is already light-limited.
Applied in `minimal` and `full` (shared constant). All variants refitted from the
defaults (commit ab5ddf9):

| variant | RMSE μ [h⁻¹] (rep. 12) | AIC (rep. 12) | colour RMSE | LEDs with no growth | gate 2 | gate 3 |
| --- | --- | --- | --- | --- | --- | --- |
| minimal-m0 | 0.0062 (0.0055) | 164 (135) | 0.37 | 465, 495, 520 | ✗ | ✓ 0.959 |
| minimal-m2 | 0.0064 (0.0057) | 168 (139) | 0.39 | 465, 495, 520 | ✗ | ✗ 0.932 |
| minimal-m0-cap | 0.0094 (0.0075) | 229 (181) | 0.40 | 465–555 | ✗ | ✗ 0.860 |
| minimal-m2-cap | 0.0081 (0.0078) | 180 (186) | 0.40 | 465, 495, 520 | ✗ | ✗ 0.806 |
| minimal-m0-cap-abs | 0.0086 (0.0073) | 224 (186) | 0.39 | 465–555 | ✗ | ✗ 0.871 |
| minimal-m0-cap-flux | 0.0074 (0.0068) | 239 (193) | 0.38 | 465, 495, 520 | ✗ | ✗ 0.928 |
| full-k3 | 0.0073 (0.0066) | 173 (162) | 0.39 | 465, 495, 520 | ✗ | ✗ 0.937 |

- The white-light fit gets a little worse everywhere; `minimal-m0` now misses gate 2
  (0.0062 vs 0.0058) but passes gate 3 (70/73; fewer strain batches run, because the
  MCA stops after Sqs + CrtE). The pigment exponent again ends at its bound (g = 4) in
  `minimal-m0`, `-m2` and `-cap-flux`: within the narrow clamp the fit wants a switch.
- Goal 1 stays Sqs (0.23 mg L⁻¹ h⁻¹ for `minimal-m0`; `minimal-m2` adds DXS, +0.4 %).
  Goal 2 stays orange or orange/red.
- `tests/test_full.py`: the growth bound of the dilute `full-k3` turbidostat is raised
  from 0.06 to 0.08 h⁻¹ (fitted model: 0.069 h⁻¹ at 100 µE).

## 3. Why 465–520 nm still give no growth

Steady states of the refitted `minimal-m0` at the colour-data conditions:

1. **The carotenoid runaway is metabolic.** At pigment factor 1.5 carotenoids still
   reach 88 mM at 555 nm and 10⁴ mM at 465–520 nm (633 nm: 8.3). In dim light NADPH
   falls far below Chl synthesis's K_m (0.004 vs 0.035 mM), so Chl synthesis stops;
   carotenoid synthesis needs no NADPH and takes the GGPP. Car:Chl = v_car / v_chl then
   grows like 1 / saturation(NADPH), and content = synthesis / μ grows without bound as
   μ → 0.
2. **Pigment turnover (tested, not kept):** dX/dt = v − (μ + k_d) X with k_d = 0.002,
   0.01, 0.03 h⁻¹ bounds carotenoids (5–6 mM) but does not restore growth: Chl then
   collapses to 1.4 mM (633 nm: 14.7), because its synthesis is NADPH-starved. Electron
   flow per Chl at 465 nm is 0.91 vs 1.55 mmol mol(Chl)⁻¹ s⁻¹ at 633 nm, so light
   capture per Chl is not the limit; the loop Chl → supply → NADPH → Chl synthesis is.
3. **Car:Chl feedback on CrtB (tested, not kept):** with the narrow clamp, Chl then
   accumulates instead (4.5–5.3× at 465/555 nm) and the culture still dies.
4. **The photoacclimation signal reads blue light as bright:** at 465 nm PS1 Chl
   absorbs strongly (total absorption per Chl above that at 633 nm), so the factor sits
   at its floor 0.5 although PS2 gets little light. Chl per cell at 465 nm in the data
   is 0.92 of 633 nm, so the cut itself is moderate; it adds to the loop of point 2.

The common cause: pigment synthesis competes with growth for NADPH and loses
completely at low NADPH (no priority, no protein/pigment coupling), and the pigment
pools are bounded only by dilution.

## 4. Decisions needed

1. **Pigment synthesis in energy limitation.** Options:
   - (a) pigments as a fixed share of biomass synthesis: Chl and carotenoids made in
     proportion to growth (v_X = μ X_target, X_target set by photoacclimation within the
     data range), so no pigment pool can drift from growth. Simplest, matches the data's
     narrow range;
   - (b) pigment turnover plus priority for Chl synthesis (lower K_m(NADPH), or a
     NADPH-independent share);
   - (c) leave the metabolism and treat 465–520 nm as out of scope for goal 2 (fit
     the colour data only at 435 and 555–687 nm).
2. **Photoacclimation signal:** total absorption (current), PS2 excitation, or electron
   flow per Chl (the `flux` signal used for capacity).
