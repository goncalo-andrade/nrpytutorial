
// Set simulation time
REAL sim_time = 100;

// Set number of outputs
int N_outputs = 250.0;

// Set variable for the coordinate system
char coord_system[] = "SinhSpherical";

// Set Black Hole mass and spin
params.M = 1;
params.chi = 0.0;
// Set Scalar Field mass
params.mu_s = 0.0;

// Set Scalar Field gaussian profile variables
params.A = 0.15;
params.r0 = 5;
params.w = 0.8;

// Set parameters for pseudo-bound state
params.omega = 0.3929;
params.m = 1;
    
// Set free-parameter values.

const REAL domain_size    = 100;
const REAL sinh_width     = 0.2;
const REAL sinhv2_const_dr= 0.05;
const REAL SymTP_bScale   = 0.5;

params.AMPL = domain_size;
params.SINHW=  sinh_width;

