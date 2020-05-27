import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import reference_metric as rfm
import sympy as sp
import sys

thismodule = __name__
# par.set_parval_from_str("reference_metric::CoordSystem", "Spherical")


def declare_gauge_and_field_gridfunctions_if_not_declared_already():

    # Define the global field variables
    global Phi, Pi

    # Check if they were already declared
    # You can't declare a gridfunction twice
    for i in range(len(gri.glb_gridfcs_list)):
        if "Phi" in gri.glb_gridfcs_list[i].name:
            Phi, Pi = sp.symbols('Phi Pi', real=True)
            return Phi, Pi

    # Declare the gridfunctions
    Phi, Pi = gri.register_gridfunctions("EVOL", ["Phi", "Pi"])

    return Phi, Pi


# def betaU_deriv():

#     # Declare the global variables
#     global betaU, BU, betaU_dupD

#     # Set spatial dimension to 3
#     DIM = 3

#     # Get rescaled BSSN variables
#     vetU, betU, _, _, _ = declare_gauge_and_field_gridfunctions_if_not_declared_already()

#     # Un-rescale BSSN variables
#     betaU = ixp.zerorank1()
#     BU = ixp.zerorank1()
#     for i in range(DIM):
#         betaU[i] = vetU[i] * rfm.ReU[i]
#         BU[i] = betU[i] * rfm.ReU[i]

#     # Define derivatives of rescaled BSSN variables
#     vetU_dupD = ixp.declarerank2("vetU_dupD", "nosym")
#     betaU_dupD = ixp.zerorank2()
#     for i in range(DIM):
#         for j in range(DIM):
#             betaU_dupD[i][j] = vetU_dupD[i][j] * \
#                 rfm.ReU[i] + vetU[i] * rfm.ReUdD[i][j]


def flat_metric_quantities():

    # Call declare_gauge_and_field_gridfunctions_if_not_declared_already()
    _, _ = declare_gauge_and_field_gridfunctions_if_not_declared_already()

    # Set spatial dimension to 3
    DIM = 3

    # Call reference metric if it hasn't been called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Check if the coordinate system is spherical
    # This can be extended if desired
    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    if not "Spherical" in CoordSystem:
        print("Error: CoordSystem = " + CoordSystem + " unsupported!")
        sys.exit(1)

    # Define the spherical coordinate variables
    r, th, ph = sp.symbols("r th ph", real=True)
    coords = [r, th, ph]

    # Declare global 3-metric variables
    global gammabarUU, GammabarUDD
    global gammabarDD

    # Define the 3-metric gammabarDD
    gammabarDD = ixp.zerorank2()

    # Compute the diagonal values of gammabarDD
    # representing a flat metric in spherical coordinates
    gammabarDD[0][0] += 1
    gammabarDD[1][1] += r**2
    gammabarDD[2][2] += (r * sp.sin(th))**2

    # Define gammabarUU as the matrix inverse of gammabarDD
    gammabarUU, _ = ixp.symm_matrix_inverter3x3(gammabarDD)

    # Define the derivatives of gammabarDD
    # needed to compute the Christoffel symbols
    gammabarDD_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabarDD_dD[i][j][k] += sp.diff(gammabarDD[i][j], coords[k])

    # Compute the Christoffel symbols
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    GammabarUDD[i][j][k] += sp.Rational(1, 2) * gammabarUU[i][m] * \
                        (gammabarDD_dD[m][j][k] + gammabarDD_dD[m]
                         [k][j] - gammabarDD_dD[j][k][m])

    # This function replaces the spherical coordinates r, th, ph for the numerical grid xx0, xx1, xx2
    # Taken from BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear.py
    def sympify_integers__replace_rthph_or_Cartxyz(obj, rthph_or_xyz, rthph_or_xyz_of_xx):
        if isinstance(obj, int):
            return sp.sympify(obj)
        else:
            return obj.subs(rthph_or_xyz[0], rthph_or_xyz_of_xx[0]).\
                subs(rthph_or_xyz[1], rthph_or_xyz_of_xx[1]).\
                subs(rthph_or_xyz[2], rthph_or_xyz_of_xx[2])

    # Replace the spherical coordinates r, th and ph in gammabarDD, gammabarUU and GammabarUDD
    # by the computational grid coordinates xx0, xx1, xx2
    for i in range(DIM):
        for j in range(DIM):
            gammabarDD[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                gammabarDD[i][j], coords, rfm.xxSph)
            gammabarUU[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                gammabarUU[i][j], coords, rfm.xxSph)
            for k in range(DIM):
                GammabarUDD[i][j][k] = sympify_integers__replace_rthph_or_Cartxyz(
                    GammabarUDD[i][j][k], coords, rfm.xxSph)

    # Declare the global quantities for the trace of the extrisic curvature,
    # for the conformal factors psi and phi, and for related quantities
    global trK, psi, phi, exp_m4phi, phi_dD

    # Set trK to 0 (harmonic slicing)
    trK = sp.sympify(0)

    # Set the conformal factor \psi to 1
    # From the definition, \psi = (\frac{\gamma}{\bar{\gamma}})^{\frac{1}{3}}
    # Since \gamma_{ij} and \bar{\gamma}_{ij} are the same, \psi = 1
    psi = sp.sympify(1)

    # Define the conformal factor phi and related quantities as functions of psi
    phi = sp.log(psi)
    exp_m4phi = sp.exp(-4*phi)

    # Define the derivatives of phi, taken symbolically as this is not a grid function
    phi_dD = ixp.zerorank1()
    for i in range(DIM):
        phi_dD[i] = sp.diff(phi, coords[i])

    # Declare global variables for \bar{\Lambda} and derivatives
    global LambdabarU, LambdabarU_dD

    # Define these quantities as zero tensors
    LambdabarU = ixp.zerorank1()
    LambdabarU_dD = ixp.zerorank2()

    # Compute LambdabarU (should give 0 in this case)
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                LambdabarU[i] += gammabarUU[j][k] * \
                    (GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k])

    # Simplify LambdabarU, otherwise it gives a complicated expression which will just reduce to 0
    for i in range(DIM):
        LambdabarU[i] = sp.simplify(LambdabarU[i])

    # Compute the symbolic derivative of LambdabarU (should also give 0)
    # We differentiate in the computational grid coordinates as LambdabarU is already a function of these
    for i in range(DIM):
        for j in range(DIM):
            LambdabarU_dD[i][j] += sp.diff(LambdabarU[i], coords[j])

    # Convert all quantities to functions of the computational grid coordinates
    # instead of the usual spherical coordinates
    trK = sympify_integers__replace_rthph_or_Cartxyz(trK, coords, rfm.xxSph)
    psi = sympify_integers__replace_rthph_or_Cartxyz(psi, coords, rfm.xxSph)
    phi = sympify_integers__replace_rthph_or_Cartxyz(phi, coords, rfm.xxSph)
    exp_m4phi = sympify_integers__replace_rthph_or_Cartxyz(
        exp_m4phi, coords, rfm.xxSph)
    for i in range(DIM):
        phi_dD[i] = sympify_integers__replace_rthph_or_Cartxyz(
            phi_dD[i], coords, rfm.xxSph)

    # Declare global gauge variables and derivatives
    global alpha, alpha_dD, betaU, betaU_dD

    # Set alpha to be 1/psi**2
    alpha = 1/psi**2

    # Define alpha_dD
    alpha_dD = ixp.zerorank1()

    # Compute alpha_dD as the symbolic derivative of alpha
    for i in range(DIM):
        alpha_dD[i] = sp.diff(alpha, coords[i])

    # Set betaU to 0
    betaU = ixp.zerorank1()

    # Compute betaU_dD as the symbolic derivative of betaU
    betaU_dD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            betaU_dD[i][j] = sp.diff(betaU[i], coords[j])

    # Convert all quantities to functions of the computational grid variables
    alpha = sympify_integers__replace_rthph_or_Cartxyz(
        alpha, coords, rfm.xxSph)
    for i in range(DIM):
        betaU[i] = sympify_integers__replace_rthph_or_Cartxyz(
            betaU[i], coords, rfm.xxSph)
        for j in range(DIM):
            betaU_dD[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                betaU_dD[i][j], coords, rfm.xxSph)
