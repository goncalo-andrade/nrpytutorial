
void apply_inner_parity_conditions(const paramstruct *restrict params, const bc_struct *restrict bcstruct, const int NUM_GFS, const int8_t *restrict gfs_parity, REAL *restrict gfs) {

#include "set_Cparameters.h"

#pragma omp parallel for
    for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
        for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
            for (int pt = 0; pt < bcstruct->num_ib_gz_pts[which_gz]; pt++)
            {
                const int i0dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i0;
                const int i1dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i1;
                const int i2dest = bcstruct->inner[which_gz][pt].inner_bc_dest_pt.i2;
                const int i0src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i0;
                const int i1src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i1;
                const int i2src = bcstruct->inner[which_gz][pt].inner_bc_src_pt.i2;
                const int8_t *prty = bcstruct->inner[which_gz][pt].parity;
                gfs[IDX4S(which_gf, i0dest, i1dest, i2dest)] = bcstruct->inner[which_gz][pt].parity[gfs_parity[which_gf]] * gfs[IDX4S(which_gf, i0src, i1src, i2src)];
            }
        }
    }
}