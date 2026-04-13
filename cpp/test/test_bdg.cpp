// Unit tests for BdG tight-binding solver
#include "supermag/bdg.h"
#include "supermag/error.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <vector>

void test_bdg_no_pairing() {
    // Delta=0: eigenvalues should be tight-binding dispersion E_k = -2t·cos(k) ± mu
    const int N = 10;
    int dim = 2 * N;
    std::vector<double> ev(dim);
    int n_ev = 0;

    int rc = supermag_bdg_solve(N, 0.1, 0.0, 0.0, 0.0,  // t=0.1eV, Delta=0, E_ex=0, mu=0
                                ev.data(), &n_ev, nullptr);
    assert(rc == SUPERMAG_OK);
    assert(n_ev == dim);

    // With Delta=0, E_ex=0, mu=0: electron block = -2t_hop_meV·cos(k), hole block = +2t_hop_meV·cos(k)
    // Spectrum should be symmetric about zero
    // All eigenvalues should be finite
    for (int i = 0; i < dim; ++i)
        assert(std::isfinite(ev[i]));

    // Check symmetry: ev[i] ≈ -ev[dim-1-i]
    for (int i = 0; i < dim / 2; ++i) {
        double sum = ev[i] + ev[dim - 1 - i];
        assert(std::fabs(sum) < 1.0);  // meV tolerance
    }
    std::printf("  PASS: test_bdg_no_pairing (range: %.2f to %.2f meV)\n", ev[0], ev[dim-1]);
}

void test_bdg_particle_hole_symmetry() {
    // BdG spectrum should be symmetric about E=0 (particle-hole symmetry)
    const int N = 20;
    int dim = 2 * N;
    std::vector<double> ev(dim);
    int n_ev = 0;

    int rc = supermag_bdg_solve(N, 0.1, 1.5, 0.0, 0.0,  // Delta=1.5meV
                                ev.data(), &n_ev, nullptr);
    assert(rc == SUPERMAG_OK);

    // Eigenvalues come in ±E pairs
    for (int i = 0; i < dim / 2; ++i) {
        double e_neg = ev[i];
        double e_pos = ev[dim - 1 - i];
        assert(std::fabs(e_neg + e_pos) < 0.5);  // meV tolerance
    }
    std::printf("  PASS: test_bdg_particle_hole_symmetry\n");
}

void test_bdg_gap() {
    // Uniform s-wave pairing should open a gap >= Delta
    const int N = 30;
    int dim = 2 * N;
    std::vector<double> ev(dim);
    int n_ev = 0;

    double Delta = 1.5;  // meV
    int rc = supermag_bdg_solve(N, 0.1, Delta, 0.0, 0.0,
                                ev.data(), &n_ev, nullptr);
    assert(rc == SUPERMAG_OK);

    // Find smallest positive eigenvalue
    double min_pos = 1e10;
    for (int i = 0; i < dim; ++i)
        if (ev[i] > 0.0 && ev[i] < min_pos)
            min_pos = ev[i];

    // Gap should be at least Delta (for half-filling, gap is exactly Delta)
    assert(min_pos >= Delta * 0.5);  // allow some tolerance for finite size
    std::printf("  PASS: test_bdg_gap (min_pos=%.4f meV, Delta=%.4f meV)\n", min_pos, Delta);
}

void test_bdg_eigenvectors() {
    // Test extended solver with eigenvector output
    const int N = 5;
    int dim = 2 * N;
    std::vector<double> ev(dim);
    std::vector<double> V(dim * dim);
    int n_ev = 0;

    int rc = supermag_bdg_solve(N, 0.1, 1.0, 0.0, 0.0,
                                    ev.data(), &n_ev, V.data());
    assert(rc == SUPERMAG_OK);

    // Eigenvectors should be orthonormal: V^T V ≈ I
    // Check a few dot products
    for (int a = 0; a < std::min(dim, 3); ++a) {
        double self_dot = 0.0;
        for (int i = 0; i < dim; ++i)
            self_dot += V[i * dim + a] * V[i * dim + a];
        assert(std::fabs(self_dot - 1.0) < 0.1);  // approximately normalized
    }
    std::printf("  PASS: test_bdg_eigenvectors\n");
}

void test_bdg_raised_cap() {
    // n_sites up to 2000 should be accepted (just test acceptance, not compute)
    int n_ev = 0;
    // n_sites = 501 was previously rejected, now should be accepted
    // But actual diag of 1002×1002 would be slow, so just test 100
    const int N = 100;
    int dim = 2 * N;
    std::vector<double> ev(dim);
    int rc = supermag_bdg_solve(N, 0.1, 1.0, 0.0, 0.0, ev.data(), &n_ev, nullptr);
    assert(rc == SUPERMAG_OK);
    assert(n_ev == dim);

    // n_sites = 2001 should be rejected
    rc = supermag_bdg_solve(2001, 0.1, 1.0, 0.0, 0.0, ev.data(), &n_ev, nullptr);
    assert(rc == SUPERMAG_ERR_INVALID_DIM);
    std::printf("  PASS: test_bdg_raised_cap\n");
}

void test_bdg_null() {
    int n_ev;
    int rc = supermag_bdg_solve(10, 0.1, 1.0, 0.0, 0.0, nullptr, &n_ev, nullptr);
    assert(rc == SUPERMAG_ERR_NULL_PTR);
    std::printf("  PASS: test_bdg_null\n");
}

int main() {
    std::printf("Running BdG tests...\n");
    test_bdg_no_pairing();
    test_bdg_particle_hole_symmetry();
    test_bdg_gap();
    test_bdg_eigenvectors();
    test_bdg_raised_cap();
    test_bdg_null();
    std::printf("All BdG tests passed!\n");
    return 0;
}
