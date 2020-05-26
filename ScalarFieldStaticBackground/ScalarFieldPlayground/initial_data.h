/*
Set up the initial data at all points on the numerical grid.
 */
void initial_data(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "set_Cparameters.h"

#pragma omp parallel for
    for(int i2=0; i2<Nxx_plus_2NGHOSTS2; i2++) {
        const REAL xx2 = xx[2][i2];
        for(int i1=0; i1<Nxx_plus_2NGHOSTS1; i1++) {
            const REAL xx1 = xx[1][i1];
            for(int i0=0; i0<Nxx_plus_2NGHOSTS0; i0++) {
                const REAL xx0 = xx[0][i0];
                   {
                         const double tmp0 = (1.0/((w)*(w)));
                         const double tmp1 = r0 + xx0;
                         const double tmp2 = ((tmp1)*(tmp1));
                         const double tmp3 = tmp0*tmp2;
                         const double tmp4 = (1.0/(xx0));
                         const double tmp5 = (1.0/2.0)*tmp4;
                         const double tmp6 = ((-r0 + xx0)*(-r0 + xx0));
                         const double tmp7 = tmp0*tmp6;
                         in_gfs[IDX4S(PHIGF,i0,i1,i2)] = tmp5*exp(-tmp7) + tmp5*exp(-tmp3);
                         in_gfs[IDX4S(PIGF,i0,i1,i2)] = tmp0*tmp4*(tmp1*exp(tmp7) + (r0 - xx0)*exp(tmp3))*exp(-tmp0*(tmp2 + tmp6));
                   }
                
            } // END LOOP: for(int i0=0; i0<Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for(int i1=0; i1<Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for(int i2=0; i2<Nxx_plus_2NGHOSTS2; i2++)
}
