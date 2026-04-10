#include "supermag/error.h"

extern "C" {

const char* supermag_error_string(int code) {
    switch (code) {
        case SUPERMAG_OK:              return "Success";
        case SUPERMAG_ERR_NULL_PTR:    return "Null pointer argument";
        case SUPERMAG_ERR_INVALID_DIM: return "Invalid dimension";
        case SUPERMAG_ERR_NO_CONVERGE: return "Did not converge (or not implemented)";
        case SUPERMAG_ERR_ALLOC:       return "Allocation failure";
        default:                       return "Unknown error";
    }
}

}
