
// Part P0.a: Set the number of ghost cells, from NRPy+'s FD_CENTDERIVS_ORDER
#define NGHOSTS 3
// Part P0.b: Set the numerical precision (REAL) to double, ensuring all floating point
//            numbers are stored to at least ~16 significant digits
#define REAL double
// Part P0.c: Set the CFL Factor. Can be overwritten at command line.
REAL CFL_FACTOR = 0.5;