# Ryazanov (2003) / Fominov (2002) — Nb/Cu₀.₄₃Ni₀.₅₇ Validation

## Citation
V.V. Ryazanov, V.A. Oboznov, A.S. Prokofiev, V.V. Bolginov, and A.K. Feofanov,
"Coupling of Two Superconductors through a Ferromagnet: Evidence for a π Junction,"
*Physical Review Letters*, vol. 86, p. 2427, 2001.

V.V. Ryazanov, V.A. Oboznov, A.Yu. Rusanov, A.V. Veretennikov, A.A. Golubov, and J. Aarts,
"Coupling of Two Superconductors through a Ferromagnet. II. Numerical and Experimental Study,"
*Physical Review Letters* (2003).

Ya.V. Fominov, N.M. Chtchelkatchev, and A.A. Golubov,
"Nonmonotonic critical temperature in superconductor/ferromagnet bilayers,"
*Physical Review B*, vol. 66, p. 014507, 2002.

## System
Nb/Cu₀.₄₃Ni₀.₅₇ bilayer (weak ferromagnet alloy).

## Physics
Experimental Tc(d_F) data from Fig. 9 of Ryazanov et al. (2003) shows:
1. Monotonic Tc suppression from ~7.1 K to ~3 K as d_F increases from 0 to ~25 nm
2. The curve is fitted by Fominov's theory (PRB 66, 014507) with parameters:
   - h ≈ 130 K (exchange field, translating to E_ex ≈ 11.2 meV)
   - γ_B ≈ 0.3 (interface barrier parameter)

## Parameters
See `params.json` for exact values.

## Expected Data
- **tc_vs_df.csv**: Digitized from Fig. 9, showing Tc vs d_F for Nb/Cu₀.₄₃Ni₀.₅₇.
  Data represents the Fominov theoretical fit curve, not raw experimental points.

## Validation Criterion
Maximum normalized deviation |Tc_computed - Tc_expected| / Tc0 < 0.05 (5%).
This is a coarser tolerance than Buzdin because the Fominov model parameters
(h, γ_B) are approximate fits to experimental data.
