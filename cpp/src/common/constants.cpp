#include "supermag/constants.h"

extern "C" {

double supermag_const_hbar(void)  { return 1.054571817e-34; }  /* J·s */
double supermag_const_kB(void)    { return 1.380649e-23; }     /* J/K */
double supermag_const_mu_B(void)  { return 9.2740100783e-24; } /* J/T */
double supermag_const_Phi0(void)  { return 2.067833848e-15; }  /* Wb  */
double supermag_const_e(void)     { return 1.602176634e-19; }  /* C   */
double supermag_const_m_e(void)   { return 9.1093837015e-31; } /* kg  */

}
