
// Set simulation time
REAL sim_time = 100;

// Set number of outputs
int N_outputs = 100;

// Set variable for the coordinate system
char coord_system[] = "SinhSpherical";

// Set Black Hole mass
params.M = 1;
// Set parameters for pseudo-bound state
params.omega = 0.3929;
params.m = 1;
    
// Set free-parameter values.

const REAL domain_size    = 500;
const REAL sinh_width     = 0.2;
const REAL sinhv2_const_dr= 0.05;
const REAL SymTP_bScale   = 0.5;

params.AMPL = domain_size;
params.SINHW=  sinh_width;

