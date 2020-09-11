/*
 * Find the CFL-constrained timestep
 */
REAL find_timestep(const paramstruct *restrict params, REAL *restrict xx[3]) {
#include "./set_Cparameters.h"
REAL dsmin = 1e38; // Start with a crazy high value... close to the largest number in single precision.
    for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++) {
        const REAL xx2 = xx[2][i2];
        for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++) {
            const REAL xx1 = xx[1][i1];
            for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++) {
                const REAL xx0 = xx[0][i0];
                REAL ds_dirn0, ds_dirn1, ds_dirn2;
                /*
                 *  Original SymPy expressions:
                 *  "[ds_dirn0 = AMPL*dxx0*(exp(xx0/SINHW)/SINHW + exp(-xx0/SINHW)/SINHW)/(exp(1/SINHW) - exp(-1/SINHW)),
                 *    ds_dirn1 = AMPL*dxx1*(exp(xx0/SINHW) - exp(-xx0/SINHW))/(exp(1/SINHW) - exp(-1/SINHW)),
                 *    ds_dirn2 = AMPL*dxx2*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)/(exp(1/SINHW) - exp(-1/SINHW))]"
                 */
                {
                   const double tmp_0 = (1.0/(SINHW));
                   const double tmp_2 = exp(tmp_0*xx0);
                   const double tmp_3 = exp(-tmp_0*xx0);
                   const double tmp_4 = AMPL/(exp(tmp_0) - exp(-tmp_0));
                   const double tmp_5 = tmp_4*(tmp_2 - tmp_3);
                   ds_dirn0 = dxx0*tmp_4*(tmp_0*tmp_2 + tmp_0*tmp_3);
                   ds_dirn1 = dxx1*tmp_5;
                   ds_dirn2 = dxx2*tmp_5*sin(xx1);
                }
                
                #ifndef MIN
                #define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
                #endif
                        // Set dsmin = MIN(dsmin, ds_dirn0, ds_dirn1, ds_dirn2);
                        dsmin = MIN(dsmin,MIN(ds_dirn0,MIN(ds_dirn1,ds_dirn2)));
                
            } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
        } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
return dsmin*CFL_FACTOR/wavespeed;
}
