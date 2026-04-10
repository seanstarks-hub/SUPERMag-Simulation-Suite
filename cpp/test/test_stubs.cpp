// Compile-check that all stub headers are valid C++ and link correctly.
#include "supermag/error.h"
#include "supermag/constants.h"
#include "supermag/proximity.h"
#include "supermag/usadel.h"
#include "supermag/eilenberger.h"
#include "supermag/bdg.h"
#include "supermag/ginzburg_landau.h"
#include "supermag/josephson.h"
#include "supermag/triplet.h"
#include <cstdio>
#include <cassert>

int main() {
    std::printf("Running stub compile-check tests...\n");

    // Verify error string works
    assert(supermag_error_string(SUPERMAG_OK) != nullptr);
    assert(supermag_error_string(SUPERMAG_ERR_NO_CONVERGE) != nullptr);

    // Verify constants are non-zero
    assert(supermag_const_hbar() > 0);
    assert(supermag_const_kB() > 0);
    assert(supermag_const_mu_B() > 0);
    assert(supermag_const_e() > 0);

    // Verify stubs return SUPERMAG_ERR_NO_CONVERGE
    assert(supermag_usadel_solve(0,0,0,0,0,0,0,nullptr,nullptr) == SUPERMAG_ERR_NO_CONVERGE);
    assert(supermag_eilenberger_solve(0,0,0,0,0,0,nullptr,nullptr) == SUPERMAG_ERR_NO_CONVERGE);
    assert(supermag_bdg_solve(0,0,0,0,nullptr,nullptr) == SUPERMAG_ERR_NO_CONVERGE);
    assert(supermag_gl_minimize(0,0,0,0,0,0,nullptr,nullptr) == SUPERMAG_ERR_NO_CONVERGE);
    assert(supermag_josephson_cpr(0,0,0,0,0,nullptr,nullptr) == SUPERMAG_ERR_NO_CONVERGE);
    assert(supermag_triplet_solve(0,nullptr,nullptr,0,nullptr,nullptr) == SUPERMAG_ERR_NO_CONVERGE);

    std::printf("  PASS: all headers compile and link\n");
    std::printf("  PASS: all stubs return ERR_NO_CONVERGE\n");
    std::printf("All stub tests passed!\n");
    return 0;
}
