/*
 * Set up the scalar field initial data at all points on the numerical grid.
 */
void fields_initial_data(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "./set_Cparameters.h"

    #pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        const REAL xx2 = xx[2][i2];
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            const REAL xx1 = xx[1][i1];
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                const REAL xx0 = xx[0][i0];
                   {
                         const double tmp_0 = (1.0/((w)*(w)));
                         const double tmp_1 = (1.0/2.0)*exp(-tmp_0*((r0 + xx0)*(r0 + xx0)));
                         const double tmp_2 = (1.0/2.0)*exp(-tmp_0*((-r0 + xx0)*(-r0 + xx0)));
                         in_gfs[IDX4S(PHIGF,i0,i1,i2)] = A*tmp_1 + A*tmp_2;
                         in_gfs[IDX4S(PIGF,i0,i1,i2)] = -A*(-tmp_0*tmp_1*xx0*(-2*r0 - 2*xx0) - tmp_0*tmp_2*xx0*(-2*r0 + 2*xx0) - tmp_1 + tmp_2)/xx0;
                   }
                
            } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
