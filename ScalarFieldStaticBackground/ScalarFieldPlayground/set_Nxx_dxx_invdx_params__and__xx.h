
void set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3], 
                                       paramstruct *restrict params, REAL *restrict xx[3]) {
    // Override parameter defaults with values based on command line arguments and NGHOSTS.
    params->Nxx0 = Nxx[0];
    params->Nxx1 = Nxx[1];
    params->Nxx2 = Nxx[2];
    params->Nxx_plus_2NGHOSTS0 = Nxx[0] + 2*NGHOSTS;
    params->Nxx_plus_2NGHOSTS1 = Nxx[1] + 2*NGHOSTS;
    params->Nxx_plus_2NGHOSTS2 = Nxx[2] + 2*NGHOSTS;
    // Step 0d: Set up space and time coordinates
    // Step 0d.i: Declare \Delta x^i=dxx{0,1,2} and invdxx{0,1,2}, as well as xxmin[3] and xxmax[3]:
#include "set_Cparameters.h"
    REAL xxmin[3],xxmax[3];
    if(EigenCoord == 0) {
        xxmin[0] = 0;
        xxmax[0] = RMAX;
        xxmin[1] = 0;
        xxmax[1] = M_PI;
        xxmin[2] = -M_PI;
        xxmax[2] = M_PI;

    } else if (EigenCoord == 1) {
        xxmin[0] = 0;
        xxmax[0] = RMAX;
        xxmin[1] = 0;
        xxmax[1] = M_PI;
        xxmin[2] = -M_PI;
        xxmax[2] = M_PI;

    }

    params->dxx0 = (xxmax[0] - xxmin[0]) / ((REAL)Nxx[0]);
    params->dxx1 = (xxmax[1] - xxmin[1]) / ((REAL)Nxx[1]);
    params->dxx2 = (xxmax[2] - xxmin[2]) / ((REAL)Nxx[2]);
    params->invdx0 = 1.0/params->dxx0;
    params->invdx1 = 1.0/params->dxx1;
    params->invdx2 = 1.0/params->dxx2;

    // Now that params.dxx{0,1,2} and params.invdxx{0,1,2} have been set,
    // Step 0d.iii: Set up uniform coordinate grids
    xx[0] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS0);
    for(int j=0;j<Nxx_plus_2NGHOSTS0;j++) 
        xx[0][j] = xxmin[0] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx0; // Cell-centered grid.
    xx[1] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS1);
    for(int j=0;j<Nxx_plus_2NGHOSTS1;j++) 
        xx[1][j] = xxmin[1] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx1; // Cell-centered grid.
    xx[2] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS2);
    for(int j=0;j<Nxx_plus_2NGHOSTS2;j++) 
        xx[2][j] = xxmin[2] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx2; // Cell-centered grid.
    //fprintf(stderr,"hey inside setxx: %e %e %e | %e %e\n",xxmin[0],xxmin[1],xxmin[2],xx[0][0],params->dxx0);
}
