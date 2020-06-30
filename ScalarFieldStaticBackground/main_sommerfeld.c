#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fenv.h> // Used to enable floating point exceptions

// Declare important global quantities and macros

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

// This returns the index for gridfunction g, evaulated at the point with spatial coordinates
// corresponding to the indices (i, j, k)
#define IDX4S(g, i, j, k) \
    ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k) + Nxx_plus_2NGHOSTS2 * (g))))
// This returns the index for gridfunction, evaluated at the point corresponding to the FULL
// spatial coordinate index idx
#define IDX4ptS(g, idx) ((idx) + (Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2) * (g))
// This returns the FULL spatial coordinate index for the point with spatial coordinates
// corresponding to indices (i, j, k)
#define IDX3S(i, j, k) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k))))
// This defines a loop through a spatial region of indices (i0min, i1min, i2min) to (i0max, i1max, i2max)
#define LOOP_REGION(i0min, i0max, i1min, i1max, i2min, i2max) \
    for (int i2 = i2min; i2 < i2max; i2++)                    \
        for (int i1 = i1min; i1 < i1max; i1++)                \
            for (int i0 = i0min; i0 < i0max; i0++)
// This defines a loop through all gridfunctions, evaluated at all grid points, in
// a sequential manner corresponding to the order by which they are written in memory
#define LOOP_ALL_GFS_GPS(ii) _Pragma("omp parallel for") for (int(ii) = 0; (ii) < Nxx_plus_2NGHOSTS_tot * NUM_EVOL_GFS; (ii)++)

// Define Upwinding algorithm
#ifndef UPWIND_ALG
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0
#endif

#include "REAL__NGHOSTS__CFL_FACTOR.h"
#include "declare_Cparameters_struct.h"
#include "boundary_conditions/gridfunction_defines.h"
#include "xxCart.h"
#include "set_Nxx_dxx_invdx_params__and__xx.h"
#include "boundary_conditions/CurviBC_include_Cfunctions.h"
#include "find_timestep.h"
#include "initial_data.h"
//#include "SIMD/SIMD_intrinsics.h"
#include "rhs_eval.h"
#include "print_ghosts.h"
#include "apply_inner_parity_conditions.h"

int main(int argc, const char *argv[])
{

    // Enable floating point exceptions
    // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
    // feenableexcept(FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID);

    // Define parameter struct
    paramstruct params;
// Include default parameter values
#include "set_Cparameters_default.h"
    // Set default initial time
    REAL t_initial = 0;

    // Read command line input and print errors
    // Expected arguments: Nxx0 Nxx1 Nxx2 [CFL_FACTOR] [A00/A11] [r0] [w] [mu_s]
    // Arguments in brackets are optional

    // Check if enough were provided to establish the numerical grid
    if (argc < 4 || atoi(argv[1]) < NGHOSTS || atoi(argv[2]) < NGHOSTS || atoi(argv[3]) < 2 /* FIXME; allow for axisymmetric sims */)
    {
        fprintf(stderr, "ERROR: Expected three command-line arguments: ./main Nx0 Nx1 Nx2,\n");
        fprintf(stderr, "       where Nx[0,1,2] is the number of grid points in the 0, 1, and 2 directions.\n");
        fprintf(stderr, "       Nx[] MUST BE larger than NGHOSTS (= %d)\n", NGHOSTS);
        exit(1);
    }

    // Check if the CFL_FACTOR was provided and check its (theoretical) stability
    // given the axisymmetry (or lack thereof) conditions
    // Warn the user they didn't provide custom parameters for the scalar field
    if (argc == 5)
    {
        CFL_FACTOR = strtod(argv[4], NULL);
        if (CFL_FACTOR > 0.5 && atoi(argv[3]) != 2)
        {
            fprintf(stderr, "WARNING: CFL_FACTOR was set to %e, which is > 0.5.\n", CFL_FACTOR);
            fprintf(stderr, "         This will generally only be stable if the simulation is purely axisymmetric.\n");
            fprintf(stderr, "         However, Nx2 was set to %d>2, which implies a non-axisymmetric simulation\n.", atoi(argv[3]));
        }

        fprintf(stderr, "WARNING: values for the scalar field parameters were not provided. Default values will be taken instead.\n");
        fprintf(stderr, "Default values are:\nA = %.2f\tr0 = %.2f\tw = %.2f\tmu_s = %.2f\n", params.A, params.r0, params.w, params.mu_s);
    }

    // Check if all field parameters were provided and if they're > 0
    // except for mu_s, which can also be = 0
    if (argc == 9)
    {

        CFL_FACTOR = strtod(argv[4], NULL);
        if (CFL_FACTOR > 0.5 && atoi(argv[3]) != 2)
        {
            fprintf(stderr, "WARNING: CFL_FACTOR was set to %e, which is > 0.5.\n", CFL_FACTOR);
            fprintf(stderr, "         This will generally only be stable if the simulation is purely axisymmetric.\n");
            fprintf(stderr, "         However, Nx2 was set to %d>2, which implies a non-axisymmetric simulation\n.", atoi(argv[3]));
        }

        REAL A = strtod(argv[5], NULL);
        REAL r0 = strtod(argv[6], NULL);
        REAL w = strtod(argv[7], NULL);
        REAL mu_s = strtod(argv[8], NULL);

        if (A <= 0)
        {
            fprintf(stderr, "WARNING: scalar field amplitude was set to %.2f, which is <= 0. \n", A);
            fprintf(stderr, "         This might produce unexpected results.\n");
        }
        params.A = A;

        if (r0 < 0)
        {
            fprintf(stderr, "ERROR: scalar field central radius r0 was set to %.2f, which is < 0.\n", r0);
            exit(1);
        }
        else
        {
            params.r0 = r0;
        }

        if (w <= 0)
        {
            fprintf(stderr, "ERROR: scalar field distribution width w was set to %.2f, which is <= 0.\n", w);
            exit(1);
        }
        else
        {
            params.w = w;
        }

        if (mu_s < 0)
        {
            fprintf(stderr, "ERROR: scalar field mass parameter mu_s was set to %.2f, which is < 0.\n", mu_s);
            exit(1);
        }
        else
        {
            params.mu_s = mu_s;
        }
    }

    if (argc > 5 && argc < 9)
    {
        fprintf(stderr, "ERROR: not enough scalar field parameters were provided.\n");
        fprintf(stderr, "       Scalar field command line inputs are: A r0 w mu_s\n");
        exit(1);
    }

    // Set up the numerical grid
    // Produces an error if number of gridpoints in at least one direction is not even
    const int Nxx[3] = {atoi(argv[1]), atoi(argv[2]), atoi(argv[3])};
    if (Nxx[0] % 2 != 0 || Nxx[1] % 2 != 0 || Nxx[2] % 2 != 0)
    {
        fprintf(stderr, "ERROR: Cannot guarantee a proper cell-centered grid if number of grid cells not set to even number.\n");
        fprintf(stderr, "       For example, in case of angular directions, proper symmetry zones will not exist.\n");
        exit(1);
    }

// Set free parameters, which overwrite the defaults
#include "free_parameters.h"

    // Coordinate grid
    REAL *xx[3];

    // Set bcstruct
    bc_struct bcstruct;

    // Compute bsstruct quantities in the eigen-coordinate system
    {
        int EigenCoord = 1;

        // Set Nxx, Nxx_plus_2NGHOSTS, dxx, indvx and xx[] for the chosen eigen-coord. system
        set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &params, xx);

// Set Nxx_plus_2NGHOSTS_tot
#include "set_Cparameters-nopointer.h"
        const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

// Find the ghost zone mappings and set up bcstruct
#include "boundary_conditions/driver_bcstruct.h"

        // Free allocated space for xx[][] array
        for (int i = 0; i < 3; i++)
            free(xx[i]);
    }

    // Set EigenCoord to 0 to allow the computation to consider the chosen coordinate system
    int EigenCoord = 0;

    // Set Nxx, Nxx_plus_2NGHOSTS, dxx, indvx and xx[] for the chosen coordinate system
    set_Nxx_dxx_invdx_params__and__xx(EigenCoord, Nxx, &params, xx);

// Set Nxx_plus_2NGHOSTS_tot
#include "set_Cparameters-nopointer.h"
    const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

    // Time parameters

    // Set final time so that the approximate outer BCs don't contaminate the data at the origin
    const REAL t_final = t_initial + 20.0;

    // Timestep based on the CFL condition
    REAL dt = find_timestep(&params, xx);

    // Number of points in time
    int N_final = (int)((t_final - t_initial) / dt + 0.5); // Add 0.5 to account for C rounding down

    // Number of iterations before outputting data
    int output_every_N = (int)((REAL)N_final / 50.0);
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

    // Set up initial data
    initial_data(&params, xx, y_n_gfs);

    // Apply inner parity conditions to inital data
    apply_inner_parity_conditions(&params, &bcstruct, NUM_EVOL_GFS, evol_gf_parity, y_n_gfs);

    // We don't apply boundary conditions to the initial data
    // Because Sommerfeld BCs are applied within the RK evolution

/************************ START TIMER ************************/

// Start timer to keep track of simulation times
#ifdef __linux__ // Use high-precision timer in Linux.
        struct timespec start,
        end;
    clock_gettime(CLOCK_REALTIME, &start);
#else // Resort to low-resolution, standards-compliant timer in non-Linux OSs
        // http://www.cplusplus.com/reference/ctime/time/
        time_t start_timer,
        end_timer;
    time(&start_timer);                                                   // Resolution of one second...
#endif

    // Integrate forward in time
    for (int n = 0; n <= N_final; n++)
    {

        // Output
        if (n % output_every_N == 0)
        {

            // Create variable to determine the current time of the simulation
            double current_t = n * dt + t_initial;

            // File to print scalar field and derivative
            char filename_num_all[100];
            sprintf(filename_num_all, "fields_num_all_%d_%d_%d_%.5f.txt", Nxx[0], Nxx[1], Nxx[2], current_t);
            FILE *out_num_all = fopen(filename_num_all, "w");

            LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
                        NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
                        NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
            {

                const int idx = IDX3S(i0, i1, i2);
                // Print fields
                fprintf(out_num_all, "%e %e %e %e %e\n", xx[0][i0], xx[1][i1], xx[2][i2],
                        y_n_gfs[IDX4ptS(PHIGF, idx)], y_n_gfs[IDX4ptS(PIGF, idx)]);
            }

            fclose(out_num_all);

            // Print interior ghost points
            char filename_ghosts[100];
            sprintf(filename_ghosts, "inner_ghosts_%d_%d_%d_%.5f.txt", Nxx[0], Nxx[1], Nxx[2], current_t);
            FILE *inner_ghosts_out = fopen(filename_ghosts, "w");
            print_ib_ghosts(&params, xx, &bcstruct, y_n_gfs, inner_ghosts_out);
            fclose(inner_ghosts_out);
        }

// Step forward in time
#include "MoLtimestepping/RK_MoL.h"

        // Output final data
        if (n == N_final - 1)
        {

            // Create variable to determine the current time of the simulation
            double current_t = n * dt + t_initial;

            // File to print scalar field and derivative
            char filename_num_all[100];
            sprintf(filename_num_all, "fields_num_all_%d_%d_%d_%.5f.txt", Nxx[0], Nxx[1], Nxx[2], current_t);
            FILE *out_num_all = fopen(filename_num_all, "w");

            LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
                        NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
                        NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS)
            {

                const int idx = IDX3S(i0, i1, i2);
                // Print fields
                fprintf(out_num_all, "%e %e %e %e %e\n", xx[0][i0], xx[1][i1], xx[2][i2],
                        y_n_gfs[IDX4ptS(PHIGF, idx)], y_n_gfs[IDX4ptS(PIGF, idx)]);
            }

            fclose(out_num_all);

            // Print interior ghost points
            char filename_ghosts[100];
            sprintf(filename_ghosts, "inner_ghosts_%d_%d_%d_%.5f.txt", Nxx[0], Nxx[1], Nxx[2], current_t);
            FILE *inner_ghosts_out = fopen(filename_ghosts, "w");
            print_ib_ghosts(&params, xx, &bcstruct, y_n_gfs, inner_ghosts_out);
            fclose(inner_ghosts_out);
        }

// Print progress indicator
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
                    n, n * (double)dt + (double)t_initial, (double)dt, (double)(100.0 * (REAL)n / (REAL)N_final),
                    (double)time_remaining_in_mins * 60, (double)(dt * 3600.0 / s_per_iteration_avg), (double)RHS_pt_evals_per_sec);
            fflush(stderr); // Flush the stderr buffer
        }                   // End progress indicator if(n % 10 == 0)
    }                       // End main loop to progress forward in time.
    fprintf(stderr, "\n");  // Clear the final line of output from progress indicator.

/************************ FREE MEMORY ************************/

// Free all allocated memory
#include "boundary_conditions/bcstruct_freemem.h"
#include "MoLtimestepping/RK_Free_Memory.h"
    free(auxevol_gfs);
    for (int i = 0; i < 3; i++)
        free(xx[i]);

    return 0;
}