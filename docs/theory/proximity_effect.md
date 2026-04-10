# S/F Proximity Effect

The proximity effect in superconductor/ferromagnet (S/F) heterostructures describes how
Cooper pairs leak from S into F. The exchange field in the ferromagnet causes the
pair amplitude to oscillate and decay — an FFLO-like effect.

## Key Equations

**Pair amplitude in the ferromagnet:**
$$F(x) = F_0 \exp\left(-\frac{x}{\xi_F}\right) \cos\left(\frac{x}{\xi_F}\right)$$

where $\xi_F = \sqrt{\hbar D_F / E_{\text{ex}}}$ is the ferromagnet coherence length.

**Critical temperature suppression (single-mode, Buzdin):**
$$T_c(d_F) = T_{c0}\left[1 - \frac{\xi_S}{d_S} \operatorname{Re}\left(\frac{1 - e^{-2\kappa_F d_F}}{\kappa_F \xi_S}\right)\right]$$

with complex wavevector $\kappa_F = (1 + i)/\xi_F$.

## References
- Buzdin, A.I. et al., JETP Lett. **35**, 178 (1982).
- Radović, Z. et al., Phys. Rev. B **44**, 759 (1991).
- Buzdin, A.I., Rev. Mod. Phys. **77**, 935 (2005).
