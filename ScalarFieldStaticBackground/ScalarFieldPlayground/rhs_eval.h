/*
Evaluate the RHSs
 */
void rhs_eval(const paramstruct *restrict params, REAL ** xx, const REAL *restrict in_gfs, REAL *restrict rhs_gfs) {
#include "set_Cparameters.h"

#pragma omp parallel for
    for(int i2=NGHOSTS; i2<NGHOSTS+Nxx2; i2++) {
        const REAL xx2 = xx[2][i2];
        for(int i1=NGHOSTS; i1<NGHOSTS+Nxx1; i1++) {
            const REAL xx1 = xx[1][i1];
            for(int i0=NGHOSTS; i0<NGHOSTS+Nxx0; i0++) {
                const REAL xx0 = xx[0][i0];
                {
                   /* 
                    * NRPy+ Finite Difference Code Generation, Step 1 of 2: Read from main memory and compute finite difference stencils:
                    */
                   const double Phi_i0_i1_i2m2 = in_gfs[IDX4S(PHIGF, i0,i1,i2-2)];
                   const double Phi_i0_i1_i2m1 = in_gfs[IDX4S(PHIGF, i0,i1,i2-1)];
                   const double Phi_i0_i1m2_i2 = in_gfs[IDX4S(PHIGF, i0,i1-2,i2)];
                   const double Phi_i0_i1m1_i2 = in_gfs[IDX4S(PHIGF, i0,i1-1,i2)];
                   const double Phi_i0m2_i1_i2 = in_gfs[IDX4S(PHIGF, i0-2,i1,i2)];
                   const double Phi_i0m1_i1_i2 = in_gfs[IDX4S(PHIGF, i0-1,i1,i2)];
                   const double Phi = in_gfs[IDX4S(PHIGF, i0,i1,i2)];
                   const double Phi_i0p1_i1_i2 = in_gfs[IDX4S(PHIGF, i0+1,i1,i2)];
                   const double Phi_i0p2_i1_i2 = in_gfs[IDX4S(PHIGF, i0+2,i1,i2)];
                   const double Phi_i0_i1p1_i2 = in_gfs[IDX4S(PHIGF, i0,i1+1,i2)];
                   const double Phi_i0_i1p2_i2 = in_gfs[IDX4S(PHIGF, i0,i1+2,i2)];
                   const double Phi_i0_i1_i2p1 = in_gfs[IDX4S(PHIGF, i0,i1,i2+1)];
                   const double Phi_i0_i1_i2p2 = in_gfs[IDX4S(PHIGF, i0,i1,i2+2)];
                   const double Pi = in_gfs[IDX4S(PIGF, i0,i1,i2)];
                   const double tmpFD0 = (1.0/12.0)*Phi_i0m2_i1_i2;
                   const double tmpFD1 = -1.0/12.0*Phi_i0p2_i1_i2;
                   const double tmpFD2 = (1.0/12.0)*Phi_i0_i1m2_i2;
                   const double tmpFD3 = -1.0/12.0*Phi_i0_i1p2_i2;
                   const double tmpFD4 = -5.0/2.0*Phi;
                   const double Phi_dD0 = invdx0*(-2.0/3.0*Phi_i0m1_i1_i2 + (2.0/3.0)*Phi_i0p1_i1_i2 + tmpFD0 + tmpFD1);
                   const double Phi_dD1 = invdx1*(-2.0/3.0*Phi_i0_i1m1_i2 + (2.0/3.0)*Phi_i0_i1p1_i2 + tmpFD2 + tmpFD3);
                   const double Phi_dDD00 = ((invdx0)*(invdx0))*((4.0/3.0)*Phi_i0m1_i1_i2 + (4.0/3.0)*Phi_i0p1_i1_i2 - tmpFD0 + tmpFD1 + tmpFD4);
                   const double Phi_dDD11 = ((invdx1)*(invdx1))*((4.0/3.0)*Phi_i0_i1m1_i2 + (4.0/3.0)*Phi_i0_i1p1_i2 - tmpFD2 + tmpFD3 + tmpFD4);
                   const double Phi_dDD22 = ((invdx2)*(invdx2))*((4.0/3.0)*Phi_i0_i1_i2m1 - 1.0/12.0*Phi_i0_i1_i2m2 + (4.0/3.0)*Phi_i0_i1_i2p1 - 1.0/12.0*Phi_i0_i1_i2p2 + tmpFD4);
                   /* 
                    * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                    */
                   const double tmp0 = (1.0/((xx0)*(xx0)));
                   const double tmp1 = sin(xx1);
                   rhs_gfs[IDX4S(PHIGF, i0, i1, i2)] = -Pi;
                   rhs_gfs[IDX4S(PIGF, i0, i1, i2)] = Phi*((mu_s)*(mu_s)) - 2*Phi_dD0/xx0 - Phi_dD1*tmp0*cos(xx1)/tmp1 - Phi_dDD00 - Phi_dDD11*tmp0 - Phi_dDD22*tmp0/((tmp1)*(tmp1));
                }
                
                
            } // END LOOP: for(int i0=NGHOSTS; i0<NGHOSTS+Nxx0; i0++)
        } // END LOOP: for(int i1=NGHOSTS; i1<NGHOSTS+Nxx1; i1++)
    } // END LOOP: for(int i2=NGHOSTS; i2<NGHOSTS+Nxx2; i2++)
}
