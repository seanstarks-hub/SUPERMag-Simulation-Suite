#include "supermag/proximity.h"

extern "C" {

double supermag_depairing_total(const supermag_depairing_t *dp) {
    if (!dp) return 0.0;
    return dp->ag + dp->zeeman + dp->orbital + dp->spin_orbit;
}

}
