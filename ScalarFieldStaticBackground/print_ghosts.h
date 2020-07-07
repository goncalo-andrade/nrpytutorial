void print_ob_ghost_rhs(const paramstruct *restrict params, REAL *restrict xx[3], const bc_struct *restrict bcstruct, REAL *restrict rhs_gfs, FILE *out) {

#include "set_Cparameters.h"

    for (int which_gz = 0; which_gz < NGHOSTS; ++which_gz) {

        for (int pt = 0; pt < bcstruct->num_ob_gz_pts[which_gz]; ++pt) {

            const int i0 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i0;
            const int i1 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i1;
            const int i2 = bcstruct->outer[which_gz][pt].outer_bc_dest_pt.i2;

            const int idx = IDX3S(i0, i1, i2);
            fprintf(out, "%d %e %e %e %e %e\n", which_gz, xx[0][i0], xx[1][i1], xx[2][i2], rhs_gfs[IDX4ptS(PHIGF, idx)], rhs_gfs[IDX4ptS(PIGF, idx)]);
        } // END LOOP for (int pt = 0; pt < bcstruct->num_ob_gz_pts[which_gz]; ++pt)
    }     // END LOOP for (int which_gz = 0; which_gz < NGHOSTS; ++which_gz)
}         // END FUNCTION print_ghost_rhs

void print_ib_ghosts(const paramstruct *restrict params, REAL *restrict xx[3], const bc_struct *restrict bcstruct, REAL *restrict gfs, FILE *out) {

#include "set_Cparameters.h"

    for (int which_gz = 0; which_gz < NGHOSTS; ++which_gz) {
        for (int pt = 0; pt < bcstruct->num_ib_gz_pts[which_gz]; ++pt) {

            const int i0 = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
            const int i1 = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
            const int i2 = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;

            const int idx = IDX3S(i0, i1, i2);
            fprintf(out, "%d %e %e %e %e %e\n", which_gz, xx[0][i0], xx[1][i1], xx[2][i2], gfs[IDX4ptS(PHIGF, idx)], gfs[IDX4ptS(PIGF, idx)]);
        }
    }
}
