

void contraction_term(const paramstruct *restrict params, const int which_gf, const REAL *restrict gfs, REAL *restrict xx[3],
           const int8_t FACEXi[3], const int i0, const int i1, const int i2, REAL *restrict _r, REAL *restrict _partial_i_f) {

#include "RELATIVE_PATH__set_Cparameters.h" /* Header file containing correct #include for set_Cparameters.h;
                                             * accounting for the relative path */

// Initialize derivatives to crazy values, to ensure that
//   we will notice in case they aren't set properly.
REAL fdD0=1e100;
REAL fdD1=1e100;
REAL fdD2=1e100;

REAL xx0 = xx[0][i0];
REAL xx1 = xx[1][i1];
REAL xx2 = xx[2][i2];

int8_t SHIFTSTENCIL0;
int8_t SHIFTSTENCIL1;
int8_t SHIFTSTENCIL2;


// On a +x0 or -x0 face, do up/down winding as appropriate:
if(abs(FACEXi[0])==1 || i0+NGHOSTS >= Nxx_plus_2NGHOSTS0 || i0-NGHOSTS <= 0) {
    int8_t SHIFTSTENCIL0 = FACEXi[0];
    if(i0+NGHOSTS >= Nxx_plus_2NGHOSTS0) SHIFTSTENCIL0 = -1;
    if(i0-NGHOSTS <= 0)                  SHIFTSTENCIL0 = +1;
    SHIFTSTENCIL1 = 0;
    SHIFTSTENCIL2 = 0;

    fdD0
        = SHIFTSTENCIL0*(-1.5*gfs[IDX4S(which_gf,i0+0*SHIFTSTENCIL0,i1+0*SHIFTSTENCIL1,i2+0*SHIFTSTENCIL2)]
                         +2.*gfs[IDX4S(which_gf,i0+1*SHIFTSTENCIL0,i1+1*SHIFTSTENCIL1,i2+1*SHIFTSTENCIL2)]
                         -0.5*gfs[IDX4S(which_gf,i0+2*SHIFTSTENCIL0,i1+2*SHIFTSTENCIL1,i2+2*SHIFTSTENCIL2)]
                        )*invdx0;

// Not on a +x0 or -x0 face, using centered difference:
} else {
    fdD0 = (gfs[IDX4S(which_gf,i0+1,i1,i2)]-gfs[IDX4S(which_gf,i0-1,i1,i2)])*0.5*invdx0;
}

/*
 *  Original SymPy expressions:
 *  "[*_r = xx0,
 *    *_partial_i_f = fdD0]"
 */
*_r = xx0;
*_partial_i_f = fdD0;

} // END contraction_term function
