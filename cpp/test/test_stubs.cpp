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
#include <cstring>

int main() {
    std::printf("Running stub compile-check tests...\n");

    // Verify error string works for all codes
    assert(supermag_error_string(SUPERMAG_OK) != nullptr);
    assert(supermag_error_string(SUPERMAG_ERR_NO_CONVERGE) != nullptr);
    assert(supermag_error_string(SUPERMAG_ERR_INVALID_MODEL) != nullptr);
    assert(supermag_error_string(SUPERMAG_ERR_NOT_IMPLEMENTED) != nullptr);

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

    // Verify proximity enums and struct compile
    supermag_proximity_params_t p;
    std::memset(&p, 0, sizeof(p));
    p.model = SUPERMAG_MODEL_THIN_S;
    p.phase = SUPERMAG_PHASE_ZERO;
    assert(p.model == SUPERMAG_MODEL_THIN_S);
    assert(p.phase == SUPERMAG_PHASE_ZERO);
    (void)p;

    supermag_depairing_t dp = {0.0, 0.0, 0.0, 0.0};
    assert(supermag_depairing_total(&dp) == 0.0);

    std::printf("  PASS: all headers compile and link\n");
    std::printf("  PASS: all stubs return ERR_NO_CONVERGE\n");
    std::printf("  PASS: new error codes defined\n");
    std::printf("  PASS: proximity enums and structs work\n");
    std::printf("All stub tests passed!\n");
    return 0;
}
