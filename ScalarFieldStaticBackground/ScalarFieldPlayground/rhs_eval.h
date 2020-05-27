/*
Evaluate the RHSs
 */
void rhs_eval(const paramstruct *restrict params, REAL ** xx, const REAL *restrict in_gfs, REAL *restrict rhs_gfs) {
#include "./set_Cparameters.h"

    #pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
                const REAL xx2 = xx[2][i2];
        for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
                        const REAL xx1 = xx[1][i1];
            for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
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
                   const double FDPart1_Rational_2_3 = 2.0/3.0;
                   const double FDPart1_Rational_1_12 = 1.0/12.0;
                   const double FDPart1_Rational_5_2 = 5.0/2.0;
                   const double FDPart1_Rational_4_3 = 4.0/3.0;
                   const double FDPart1_1 = -Phi_i0_i1p2_i2;
                   const double FDPart1_2 = -FDPart1_Rational_5_2*Phi;
                   const double Phi_dD0 = invdx0*(FDPart1_Rational_1_12*(Phi_i0m2_i1_i2 - Phi_i0p2_i1_i2) + FDPart1_Rational_2_3*(-Phi_i0m1_i1_i2 + Phi_i0p1_i1_i2));
                   const double Phi_dD1 = invdx1*(FDPart1_Rational_1_12*(FDPart1_1 + Phi_i0_i1m2_i2) + FDPart1_Rational_2_3*(-Phi_i0_i1m1_i2 + Phi_i0_i1p1_i2));
                   const double Phi_dDD00 = ((invdx0)*(invdx0))*(FDPart1_2 + FDPart1_Rational_1_12*(-Phi_i0m2_i1_i2 - Phi_i0p2_i1_i2) + FDPart1_Rational_4_3*(Phi_i0m1_i1_i2 + Phi_i0p1_i1_i2));
                   const double Phi_dDD11 = ((invdx1)*(invdx1))*(FDPart1_2 + FDPart1_Rational_1_12*(FDPart1_1 - Phi_i0_i1m2_i2) + FDPart1_Rational_4_3*(Phi_i0_i1m1_i2 + Phi_i0_i1p1_i2));
                   const double Phi_dDD22 = ((invdx2)*(invdx2))*(FDPart1_2 + FDPart1_Rational_1_12*(-Phi_i0_i1_i2m2 - Phi_i0_i1_i2p2) + FDPart1_Rational_4_3*(Phi_i0_i1_i2m1 + Phi_i0_i1_i2p1));
                   /* 
                    * NRPy+ Finite Difference Code Generation, Step 2 of 2: Evaluate SymPy expressions and write to main memory:
                    */
                   const double FDPart3_0 = (1.0/((xx0)*(xx0)));
                   const double FDPart3_1 = sin(xx1);
                   rhs_gfs[IDX4S(PHIGF, i0, i1, i2)] = -Pi;
                   rhs_gfs[IDX4S(PIGF, i0, i1, i2)] = -FDPart3_0*Phi_dDD11 - FDPart3_0*Phi_dD1*cos(xx1)/FDPart3_1 - FDPart3_0*Phi_dDD22/((FDPart3_1)*(FDPart3_1)) + Phi*((mu_s)*(mu_s)) - 2*Phi_dD0/xx0 - Phi_dDD00;
                }
                
                
            } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
        } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
