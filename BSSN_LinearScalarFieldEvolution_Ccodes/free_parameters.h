
// Set Black Hole mass and spin
params.M = 1;
params.chi = 0.99;

// Set Scalar Field mass
params.mu_s = 0.0;

// Set Scalar Field gaussian profile variables
params.A = 1.0;
params.r0 = 32;
params.w = 8;
    
// Set free-parameter values.

const REAL domain_size    = 128;
const REAL sinh_width     = 0.4;
const REAL sinhv2_const_dr= 0.05;
const REAL SymTP_bScale   = 0.5;

params.RMAX = domain_size;

