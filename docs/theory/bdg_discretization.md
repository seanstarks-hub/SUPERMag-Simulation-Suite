# BdG Tight-Binding Discretization

The Bogoliubov–de Gennes (BdG) approach provides a fully microscopic description
of inhomogeneous superconductors. The Hamiltonian is discretized on a tight-binding
lattice and diagonalized to find quasiparticle spectra.

## Key Equations

**BdG Hamiltonian (1D tight-binding):**
$$\hat{H} = \sum_{i,\sigma} \left(-t c_{i\sigma}^\dagger c_{i+1,\sigma} + \text{h.c.}\right) + \sum_i \left(\Delta_i c_{i\uparrow}^\dagger c_{i\downarrow}^\dagger + \text{h.c.}\right) + \sum_i h_i (n_{i\uparrow} - n_{i\downarrow})$$

**Self-consistency:**
$$\Delta_i = -V \sum_n u_n(i) v_n^*(i) \tanh\left(\frac{E_n}{2k_B T}\right)$$

## References
- de Gennes, P.G., *Superconductivity of Metals and Alloys* (1966).
- Bagwell, P.F., Phys. Rev. B **49**, 6841 (1994).
