// Transfer matrix utility for cascaded proximity layers.
// Internal header — not part of the public C API.
//
// A single layer of thickness δ and complex wave vector q has transfer matrix:
//   M = [[ cosh(qδ),  sinh(qδ)/q ],
//        [ q·sinh(qδ), cosh(qδ)  ]]
//
// Total transfer matrix for N layers: M_total = M_1 · M_2 · ... · M_N
// (ordered from S interface to vacuum boundary).
//
// Effective kernel: K_eff = M_total[1][0] / M_total[0][0]
// (ratio gives the input impedance seen from the S side).
//
// Symplectic property: det(M_i) = 1 for each layer.

#ifndef SUPERMAG_TRANSFER_MATRIX_H
#define SUPERMAG_TRANSFER_MATRIX_H

#include <complex>

namespace supermag {

struct Matrix2x2 {
    std::complex<double> m[2][2];
};

// Identity matrix
Matrix2x2 mat_identity();

// Matrix multiplication C = A * B
Matrix2x2 mat_multiply(const Matrix2x2 &A, const Matrix2x2 &B);

// Determinant
std::complex<double> mat_det(const Matrix2x2 &M);

// Build transfer matrix for a single layer with wave vector q and thickness delta.
// Uses overflow-safe sinh/cosh for large |Re(q*delta)|.
Matrix2x2 layer_transfer_matrix(std::complex<double> q, double delta);

// Extract the effective kernel from a composed transfer matrix.
// K_eff = M_total[1][0] / M_total[0][0]
std::complex<double> extract_kernel(const Matrix2x2 &M_total);

} // namespace supermag

#endif
