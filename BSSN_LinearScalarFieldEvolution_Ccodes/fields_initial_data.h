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
                         const double tmp_1 = (1.0/(SINHW));
                         const double tmp_2 = exp(tmp_1) - exp(-tmp_1);
                         const double tmp_4 = exp(tmp_1*xx0) - exp(-tmp_1*xx0);
                         const double tmp_5 = AMPL*tmp_4/tmp_2;
                         const double tmp_6 = A*exp(-((-r0 + tmp_5)*(-r0 + tmp_5))/((w)*(w)))*sin(xx1)/sqrt(M_PI);
                         const double tmp_10 = 1 + 0.5*tmp_2/(AMPL*tmp_4);
                         const double tmp_11 = ((tmp_10)*(tmp_10)*(tmp_10)*(tmp_10));
                         in_gfs[IDX4S(PHIGF,i0,i1,i2)] = tmp_6*cos(m*xx2);
                         in_gfs[IDX4S(PIGF,i0,i1,i2)] = omega*tmp_6*sin(m*xx2)/sqrt(((tmp_2)*(tmp_2))*(((AMPL)*(AMPL))*tmp_11*((tmp_4)*(tmp_4))/((tmp_2)*(tmp_2)) - 2*((tmp_10)*(tmp_10))*tmp_5)/(((AMPL)*(AMPL))*tmp_11*((tmp_4)*(tmp_4))));
                   }
                
            } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
