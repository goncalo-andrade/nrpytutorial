
// Define REAL and NGHOSTS; and declare CFL_FACTOR. This header is generated in NRPy+.
#include "BSSN_Playground_REAL__NGHOSTS__CFL_FACTOR.h"

#include "rfm_files/rfm_struct__declare.h"

#include "declare_Cparameters_struct.h"

// All SIMD intrinsics used in SIMD-enabled C code loops are defined here:
#include "SIMD/SIMD_intrinsics.h"
#ifdef SIMD_IS_DISABLED
// Algorithm for upwinding, SIMD-disabled version.
// *NOTE*: This upwinding is backwards from
//  usual upwinding algorithms, because the
//  upwinding control vector in BSSN (the shift)
//  acts like a *negative* velocity.
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0
#endif

// Import needed header files
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdint.h> // Needed for Windows GCC 6.x compatibility
#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif
#ifndef M_SQRT1_2
#define M_SQRT1_2 0.707106781186547524400844362104849039L
#endif
#define wavespeed 1.0 // Set CFL-based "wavespeed" to 1.0.

// Declare the IDX4S(gf,i,j,k) macro, which enables us to store 4-dimensions of
// data in a 1D array. In this case, consecutive values of "i"
// (all other indices held to a fixed value) are consecutive in memory, where
// consecutive values of "j" (fixing all other indices) are separated by
// Nxx_plus_2NGHOSTS0 elements in memory. Similarly, consecutive values of
// "k" are separated by Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1 in memory, etc.
#define IDX4S(g, i, j, k) \
    ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k) + Nxx_plus_2NGHOSTS2 * (g))))
#define IDX4ptS(g, idx) ((idx) + (Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2) * (g))
#define IDX3S(i, j, k) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k))))
#define LOOP_REGION(i0min, i0max, i1min, i1max, i2min, i2max) \
    for (int i2 = i2min; i2 < i2max; i2++)                    \
        for (int i1 = i1min; i1 < i1max; i1++)                \
            for (int i0 = i0min; i0 < i0max; i0++)
#define LOOP_ALL_GFS_GPS(ii) _Pragma("omp parallel for") for (int(ii) = 0; (ii) < Nxx_plus_2NGHOSTS_tot * NUM_EVOL_GFS; (ii)++)

// Set *GF macros, as well as xxCart()
#include "boundary_conditions/gridfunction_defines.h"

// Set xxCart(const paramstruct *restrict params,
//            REAL *restrict xx[3],
//            const int i0,const int i1,const int i2,
//            REAL xCart[3]),
// which maps xx->Cartesian via
// {xx[0][i0],xx[1][i1],xx[2][i2]}->{xCart[0],xCart[1],xCart[2]}
#include "xxCart.h"

// Defines set_Nxx_dxx_invdx_params__and__xx(const int EigenCoord, const int Nxx[3], paramstruct *restrict params, REAL *restrict xx[3]),
// which sets params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for
// the chosen Eigen-CoordSystem if EigenCoord==1, or
// CoordSystem if EigenCoord==0.
#include "set_Nxx_dxx_invdx_params__and__xx.h"

// Include basic functions needed to impose curvilinear
// parity and boundary conditions.
#include "boundary_conditions/CurviBC_include_Cfunctions.h"
#include "boundary_conditions/apply_inner_parity_conditions.h"

// Include function for enforcing detgammabar constraint.
#include "enforce_detgammabar_constraint.h"

// Find the CFL-constrained timestep
#include "find_timestep.h"

// Declare function necessary for setting up the initial data.
// Define BSSN_ID() for UIUC Black Hole initial data
// Set the generic driver function for setting up BSSN initial data
#include "initial_data.h"

// Declare function for evaluating Hamiltonian constraint (diagnostic)
#include "Hamiltonian_constraint.h"

// Declare rhs_eval function, which evaluates BSSN RHSs
#include "rhs_eval.h"

// Declare Ricci_eval function, which evaluates Ricci tensor
#include "Ricci_eval.h"

// main() function:
// Step 0: Read command-line input, set up grid structure, allocate memory for gridfunctions, set up coordinates
// Step 1: Set up initial data to an exact solution
// Step 2: Start the timer, for keeping track of how fast the simulation is progressing.
// Step 3: Integrate the initial data forward in time using the chosen RK-like Method of
//         Lines timestepping algorithm, and output periodic simulation diagnostics
// Step 3.a: Output 2D data file periodically, for visualization
// Step 3.b: Step forward one timestep (t -> t+dt) in time using
//           chosen RK-like MoL timestepping algorithm
// Step 3.c: If t=t_final, output conformal factor & Hamiltonian
//           constraint violation to 2D data file
// Step 3.d: Progress indicator printing to stderr
// Step 4: Free all allocated memory
int main(int argc, const char *argv[]) {

    paramstruct params;
#include "set_Cparameters_default.h"

    // Read command-line input, error out if nonconformant
    if ((argc != 4 && argc != 5) || atoi(argv[1]) < NGHOSTS || atoi(argv[2]) < NGHOSTS || atoi(argv[3]) < 2 /* FIXME; allow for axisymmetric sims */)
    {
        fprintf(stderr, "Error: Expected three command-line arguments: ./exec Nx0 Nx1 Nx2,\n");
        fprintf(stderr, "where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
        fprintf(stderr, "Nx[] MUST BE larger than NGHOSTS (= %d)\n", NGHOSTS);
        exit(1);
    }
    if (argc == 5)
    {
        CFL_FACTOR = strtod(argv[4], NULL);
        if (CFL_FACTOR > 0.5 && atoi(argv[3]) != 2)
        {
            fprintf(stderr, "WARNING: CFL_FACTOR was set to %e, which is > 0.5.\n", CFL_FACTOR);
            fprintf(stderr, "         This will generally only be stable if the simulation is purely axisymmetric\n");
            fprintf(stderr, "         However, Nx2 was set to %d>2, which implies a non-axisymmetric simulation\n", atoi(argv[3]));
        }
    }
    // Set up numerical grid structure, first in space...
    const int Nxx[3] = {atoi(argv[1]), atoi(argv[2]), atoi(argv[3])};
    if (Nxx[0] % 2 != 0 || Nxx[1] % 2 != 0 || Nxx[2] % 2 != 0)
    {
        fprintf(stderr, "Error: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
        fprintf(stderr, "       For example, in case of angular directions, proper symmetry zones will not exist.\n");
        exit(1);
    }

    // Set free parameters, overwriting Cparameters defaults
    //          by hand or with command-line input, as desired.
#include "free_parameters.h"

    // Uniform coordinate grids are stored to *xx[3]
    REAL *xx[3];
    // Set bcstruct
    bc_struct bcstruct;
    {
        int EigenCoord = 1;
        // Call set_Nxx_dxx_invdx_params__and__xx(), which sets
        // params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
        // chosen Eigen-CoordSystem.
        set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &params, xx);
        // Set Nxx_plus_2NGHOSTS_tot
#include "set_Cparameters-nopointer.h"
        const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
        // Find ghostzone mappings; set up bcstruct
#include "boundary_conditions/driver_bcstruct.h"
        // Free allocated space for xx[][] array
        for (int i = 0; i < 3; i++)
            free(xx[i]);
    }

    // Call set_Nxx_dxx_invdx_params__and__xx(), which sets
    // params Nxx,Nxx_plus_2NGHOSTS,dxx,invdx, and xx[] for the
    // chosen (non-Eigen) CoordSystem.
    int EigenCoord = 0;
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &params, xx);

    // Set all C parameters "blah" for params.blah, including
    // Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0, etc.
#include "set_Cparameters-nopointer.h"
    const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

    // Time coordinate parameters
    const REAL t_final = domain_size; /* Final time is set so that at t=t_final,
                                        * data at the origin have not been corrupted
                                        * by the approximate outer boundary condition */

    // Set timestep based on smallest proper distance between gridpoints and CFL factor
    REAL dt = find_timestep(&params, xx);
    //fprintf(stderr,"# Timestep set to = %e\n",(double)dt);
    int N_final = (int)(t_final / dt + 0.5); // The number of points in time.
                                             // Add 0.5 to account for C rounding down
                                             // typecasts to integers.
    int output_every_N = (int)((REAL)N_final / 10.0);
    if (output_every_N == 0)
        output_every_N = 1;

    // Error out if the number of auxiliary gridfunctions outnumber evolved gridfunctions.
    // This is a limitation of the RK method. You are always welcome to declare & allocate
    // additional gridfunctions by hand.
    if (NUM_AUX_GFS > NUM_EVOL_GFS)
    {
        fprintf(stderr, "Error: NUM_AUX_GFS > NUM_EVOL_GFS. Either reduce the number of auxiliary gridfunctions,\n");
        fprintf(stderr, "       or allocate (malloc) by hand storage for *diagnostic_output_gfs. \n");
        exit(1);
    }

    // Allocate memory for gridfunctions
#include "MoLtimestepping/RK_Allocate_Memory.h"
    REAL *restrict auxevol_gfs = (REAL *)malloc(sizeof(REAL) * NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot);

    // Set up precomputed reference metric arrays
    // Allocate space for precomputed reference metric arrays.
#include "rfm_files/rfm_struct__malloc.h"

    // Define precomputed reference metric arrays.
    {
#include "set_Cparameters-nopointer.h"
#include "rfm_files/rfm_struct__define.h"
    }

    // Set up initial data to an exact solution
    initial_data(&params, xx, y_n_gfs);

    // Apply inner parity conditions, as initial data
    // are sometimes ill-defined in ghost zones.
    // E.g., spherical initial data might not be
    // properly defined at points where r=-1.
    apply_inner_parity_conditions(&params, &bcstruct, NUM_EVOL_GFS, evol_gf_parity, y_n_gfs);
    enforce_detgammabar_constraint(&rfmstruct, &params, y_n_gfs);

    // Start the timer, for keeping track of how fast the simulation is progressing.
#ifdef __linux__ // Use high-precision timer in Linux.
    struct timespec start, end;
    clock_gettime(CLOCK_REALTIME, &start);
#else // Resort to low-resolution, standards-compliant timer in non-Linux OSs
    // http://www.cplusplus.com/reference/ctime/time/
    time_t start_timer, end_timer;
    time(&start_timer); // Resolution of one second...
#endif

    // Integrate the initial data forward in time using the chosen RK-like Method of
    // Lines timestepping algorithm, and output periodic simulation diagnostics
    for (int n = 0; n <= N_final; n++)
    { // Main loop to progress forward in time.

        // Output 2D data file periodically, for visualization
        if (n % 100 == 0)
        {
            // Evaluate Hamiltonian constraint violation
            Hamiltonian_constraint(&rfmstruct, &params, y_n_gfs, diagnostic_output_gfs);

            // char filename[100];
            // sprintf(filename, "out%d-%08d.txt", Nxx[0], n);
            // FILE *out2D = fopen(filename, "w");
            // LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
            //             NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
            //             NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
            // {
            //     const int idx = IDX3S(i0, i1, i2);
            //     REAL xx0 = xx[0][i0];
            //     REAL xx1 = xx[1][i1];
            //     REAL xx2 = xx[2][i2];
            //     REAL xCart[3];
            //     xxCart(&params, xx, i0, i1, i2, xCart);
            //     fprintf(out2D, "%e %e %e %e\n",
            //             xCart[1], xCart[2],
            //             y_n_gfs[IDX4ptS(CFGF, idx)], log10(fabs(diagnostic_output_gfs[IDX4ptS(HGF, idx)])));
            // }
            // fclose(out2D);

            char filename[100];
            sprintf(filename, "out%d-%08d.txt", Nxx[0], n);
            FILE *out = fopen(filename, "w");
            LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
                        NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
                        NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
            {
                const int idx = IDX3S(i0, i1, i2);
                REAL xx0 = xx[0][i0];
                REAL xx1 = xx[1][i1];
                REAL xx2 = xx[2][i2];
                fprintf(out, "%e %e %e %e %e\n", xx0, xx1, xx2, y_n_gfs[IDX4ptS(CFGF, idx)], diagnostic_output_gfs[IDX4ptS(HGF, idx)]);
            }

            fclose(out);
            }

        // Step forward one timestep (t -> t+dt) in time using
        // chosen RK-like MoL timestepping algorithm
#include "MoLtimestepping/RK_MoL.h"

        // If t=t_final, output conformal factor & Hamiltonian
        // constraint violation to 2D data file
        if (n == N_final - 1)
        {
            // Evaluate Hamiltonian constraint violation
            Hamiltonian_constraint(&rfmstruct, &params, y_n_gfs, diagnostic_output_gfs);

            // char filename[100];
            // sprintf(filename, "out%d.txt", Nxx[0]);
            // FILE *out2D = fopen(filename, "w");
            // const int i0MIN = NGHOSTS; // In spherical, r=Delta r/2.
            // const int i1mid = Nxx_plus_2NGHOSTS1 / 2;
            // const int i2mid = Nxx_plus_2NGHOSTS2 / 2;
            // LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
            //             NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
            //             NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
            // {
            //     REAL xx0 = xx[0][i0];
            //     REAL xx1 = xx[1][i1];
            //     REAL xx2 = xx[2][i2];
            //     REAL xCart[3];
            //     xxCart(&params, xx, i0, i1, i2, xCart);
            //     int idx = IDX3S(i0, i1, i2);
            //     fprintf(out2D, "%e %e %e %e\n", xCart[1], xCart[2], y_n_gfs[IDX4ptS(CFGF, idx)],
            //             log10(fabs(diagnostic_output_gfs[IDX4ptS(HGF, idx)])));
            // }
            // fclose(out2D);

            char filename[100];
            sprintf(filename, "out%d.txt", Nxx[0]);
            FILE *out = fopen(filename, "w");
            LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
                        NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
                        NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
            {
                const int idx = IDX3S(i0, i1, i2);
                REAL xx0 = xx[0][i0];
                REAL xx1 = xx[1][i1];
                REAL xx2 = xx[2][i2];
                fprintf(out, "%e %e %e %e %e\n", xx0, xx1, xx2, y_n_gfs[IDX4ptS(CFGF, idx)], diagnostic_output_gfs[IDX4ptS(HGF, idx)]);
            }

            fclose(out);
        }
        // Progress indicator printing to stderr

        // Measure average time per iteration
#ifdef __linux__ // Use high-precision timer in Linux.
        clock_gettime(CLOCK_REALTIME, &end);
        const long long unsigned int time_in_ns = 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
#else // Resort to low-resolution, standards-compliant timer in non-Linux OSs
        time(&end_timer);                                                 // Resolution of one second...
        REAL time_in_ns = difftime(end_timer, start_timer) * 1.0e9 + 0.5; // Round up to avoid divide-by-zero.
#endif
        const REAL s_per_iteration_avg = ((REAL)time_in_ns / (REAL)n) / 1.0e9;

        const int iterations_remaining = N_final - n;
        const REAL time_remaining_in_mins = s_per_iteration_avg * (REAL)iterations_remaining / 60.0;

        const REAL num_RHS_pt_evals = (REAL)(Nxx[0] * Nxx[1] * Nxx[2]) * 4.0 * (REAL)n; // 4 RHS evals per gridpoint for RK4
        const REAL RHS_pt_evals_per_sec = num_RHS_pt_evals / ((REAL)time_in_ns / 1.0e9);

        // Output simulation progress to stderr
        if (n % 10 == 0)
        {
            fprintf(stderr, "%c[2K", 27);                                                          // Clear the line
            fprintf(stderr, "It: %d t=%.2f dt=%.2e | %.1f%%; ETA %.0f s | t/h %.2f | gp/s %.2e\r", // \r is carriage return, move cursor to the beginning of the line
                    n, n * (double)dt, (double)dt, (double)(100.0 * (REAL)n / (REAL)N_final),
                    (double)time_remaining_in_mins * 60, (double)(dt * 3600.0 / s_per_iteration_avg), (double)RHS_pt_evals_per_sec);
            fflush(stderr); // Flush the stderr buffer
        }                   // End progress indicator if(n % 10 == 0)
    }                       // End main loop to progress forward in time.
    fprintf(stderr, "\n");  // Clear the final line of output from progress indicator.

    // Free all allocated memory
#include "rfm_files/rfm_struct__freemem.h"
#include "boundary_conditions/bcstruct_freemem.h"
#include "MoLtimestepping/RK_Free_Memory.h"
    free(auxevol_gfs);
    for (int i = 0; i < 3; i++)
        free(xx[i]);

    return 0;
}
