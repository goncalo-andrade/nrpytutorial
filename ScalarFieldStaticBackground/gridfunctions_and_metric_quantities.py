import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import reference_metric as rfm
import sympy as sp
import sys

thismodule = __name__
# Background Metric
par.initialize_param(par.glb_param(
    "char", thismodule, "BackgroundMetric", "Kerr"))
# Mass (for the Schwarzschild and Kerr solutions) and spin (for the Kerr solution)
M, chi = par.Cparameters("REAL", thismodule, ["M", "chi"], [1.0, 0.99])
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
    global gammabarDD, gammabarUU, GammabarUDD

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
                gammabarDD_dD[i][j][k] = sp.diff(gammabarDD[i][j], coords[k])

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


def schwarzschild_metric_quantities():

    # Based on the BSSN/StaticTrumpet module

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

    # Declare the global variable for the conformal factor psi
    global psi

    # Set the conformal factor psi
    psi = sp.sqrt(1 + M/r)

    # Declare global 3-metric variables
    global gammabarDD, gammabarUU, GammabarUDD

    # Define the physical 3-metric
    gammaDD = ixp.zerorank2()
    gammaDD[0][0] = psi**4
    gammaDD[1][1] = psi**4 * r**2
    gammaDD[2][2] = psi**4 * r**2 * sp.sin(th)**2

    # Define the determinant of gammaDD
    gammaUU, gammaDET = ixp.symm_matrix_inverter3x3(gammaDD)

    # Define the conformal 3-metric gammabarDD
    gammabarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammabarDD[i][j] = (
                rfm.detgammahat / gammaDET)**sp.Rational(1, 3) * gammaDD[i][j]

    # Define gammabarUU as the matrix inverse of gammabarDD
    gammabarUU, _ = ixp.symm_matrix_inverter3x3(gammabarDD)

    # Define the derivatives of gammabarDD
    # needed to compute the Christoffel symbols
    gammabarDD_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabarDD_dD[i][j][k] = sp.diff(gammabarDD[i][j], coords[k])

    # Compute the Christoffel symbols
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    GammabarUDD[i][j][k] += sp.Rational(1, 2) * gammabarUU[i][m] * \
                        (gammabarDD_dD[m][j][k] + gammabarDD_dD[m]
                         [k][j] - gammabarDD_dD[j][k][m])

    # Define the extrinsic curvature
    KDD = ixp.zerorank2()

    # Set the non-zero terms for KDD

    # K_{rr} = M / r^2
    KDD[0][0] = -M / r**2

    # K_{theta theta} = K_{phi phi} / sin^2 theta = M
    KDD[1][1] = M
    KDD[2][2] = M * sp.sin(th)**2

    # Declare global variables for the trace of the extrinsic curvaturature,
    # for the conformal factor phi and for related quantities
    global trK, phi, exp_m4phi, phi_dD

    # Compute trK from the definition
    trK = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trK += gammaUU[i][j] * KDD[i][j]

    # Define quantities related to psi
    phi = sp.log(psi)
    exp_m4phi = sp.exp(-4*phi)

    # Define the derivatives of phi
    phi_dD = ixp.zerorank2()
    for i in range(DIM):
        phi_dD[i] = sp.diff(phi, coords[i])

    # Declare global gauge variables and derivatives
    global alpha, alpha_dD, betaU, betaU_dD

    # Set \alpha = r / (r + M)
    alpha = r / (r + M)

    # Define alpha_dD
    alpha_dD = ixp.zerorank1()

    # Compute alpha_dD symbolically
    for i in range(DIM):
        alpha_dD[i] = sp.diff(alpha, coords[i])

    # Define betaU
    betaU = ixp.zerorank1()

    # Set beta^r = Mr / (r + M)^2
    betaU[0] = M * r / (r + M)**2

    # Define betaU_dD
    betaU_dD = ixp.zerorank2()

    # Compute betaU_dD symbolically
    for i in range(DIM):
        for j in range(DIM):
            betaU_dD[i][j] = sp.diff(betaU[i], coords[j])

    # This function replaces the spherical coordinates r, th, ph for the numerical grid xx0, xx1, xx2
    # Taken from BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear.py
    def sympify_integers__replace_rthph_or_Cartxyz(obj, rthph_or_xyz, rthph_or_xyz_of_xx):
        if isinstance(obj, int):
            return sp.sympify(obj)
        else:
            return obj.subs(rthph_or_xyz[0], rthph_or_xyz_of_xx[0]).\
                subs(rthph_or_xyz[1], rthph_or_xyz_of_xx[1]).\
                subs(rthph_or_xyz[2], rthph_or_xyz_of_xx[2])

    # Replace spherical coordinates for computational grid variables everywhere
    trK = sympify_integers__replace_rthph_or_Cartxyz(trK, coords, rfm.xxSph)
    psi = sympify_integers__replace_rthph_or_Cartxyz(psi, coords, rfm.xxSph)
    phi = sympify_integers__replace_rthph_or_Cartxyz(phi, coords, rfm.xxSph)
    exp_m4phi = sympify_integers__replace_rthph_or_Cartxyz(
        exp_m4phi, coords, rfm.xxSph)
    alpha = sympify_integers__replace_rthph_or_Cartxyz(
        alpha, coords, rfm.xxSph)
    for i in range(DIM):
        alpha_dD[i] = sympify_integers__replace_rthph_or_Cartxyz(
            alpha_dD[i], coords, rfm.xxSph)
        betaU[i] = sympify_integers__replace_rthph_or_Cartxyz(
            betaU[i], coords, rfm.xxSph)
        phi_dD[i] = sympify_integers__replace_rthph_or_Cartxyz(
            phi_dD[i], coords, rfm.xxSph)
        for j in range(DIM):
            gammabarDD[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                gammabarDD[i][j], coords, rfm.xxSph)
            gammabarUU[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                gammabarUU[i][j], coords, rfm.xxSph)
            betaU_dD[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                betaU_dD[i][j], coords, rfm.xxSph)
            for k in range(DIM):
                GammabarUDD[i][j][k] = sympify_integers__replace_rthph_or_Cartxyz(
                    GammabarUDD[i][j][k], coords, rfm.xxSph)


def kerr_metric_quantities():

    # Based on the BSSN/UIUCBlackHole module

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

    # Set up auxiliary quantities

    # Spin per unit mass
    a = M * chi

    # Boyer-Lindquist outer and inner horizons
    rp = M + sp.sqrt(M**2 - a**2)
    rm = M - sp.sqrt(M**2 - a**2)

    # Boyer-Lindquist radius in terms of UIUC radius
    rBL = r * (1 + rp / (4 * r))**2

    # \Sigma = r_{BL}^2 + a^2 \cos^2 \theta
    SIG = rBL**2 + a**2 * sp.cos(th)**2

    # \Delta = r_{BL}^2 - 2 M r_{BL} + a^2
    DEL = rBL**2 - 2 * M * rBL + a**2

    # A = (r_{BL}^2 + a^2)^2 - \Delta a^2 \sin^2 \theta
    AA = (rBL**2 + a**2)**2 - DEL * a**2 * sp.sin(th)**2

    # Declare global 3-metric variables
    global gammabarDD, gammabarUU, GammabarUDD

    # Define the 3-metric
    gammabarDD = ixp.zerorank2()
    gammabarDD[0][0] = (SIG * (r + rp / 4)**2) / (r**3 * (rBL - rm))
    gammabarDD[1][1] = SIG
    gammabarDD[2][2] = AA * sp.sin(th)**2 / SIG

    # Define gammabarUU as the matrix inverse of gammabarDD
    gammabarUU, _ = ixp.symm_matrix_inverter3x3(gammabarDD)

    # Define the derivatives of gammabarDD
    # needed to compute the Christoffel symbols
    gammabarDD_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabarDD_dD[i][j][k] = sp.diff(gammabarDD[i][j], coords[k])

    # Compute the Christoffel symbols
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    GammabarUDD[i][j][k] += sp.Rational(1, 2) * gammabarUU[i][m] * \
                        (gammabarDD_dD[m][j][k] + gammabarDD_dD[m]
                         [k][j] - gammabarDD_dD[j][k][m])

    # We don't need to define the extrinsic curvature
    # All we care about is its trace, which is 0

    global trK, psi, phi, exp_m4phi, phi_dD

    # Set trK to 0
    trK = sp.sympify(0)

    # Define the conformal factor psi
    psi = AA / SIG

    # Define phi and exp_m4phi
    phi = sp.log(psi)
    exp_m4phi = sp.exp(-4*phi)

    # Define the derivatives of phi
    phi_dD = ixp.zerorank1()
    for i in range(DIM):
        phi_dD[i] = sp.diff(phi, coords[i])

    # Declare global gague variables and derivatives
    global alpha, alpha_dD, betaU, betaU_dD

    # Set alpha = 1 / psi^2
    alpha = 1 / psi**2

    # Define alpha_dD
    alpha_dD = ixp.zerorank1()

    # Compute alpha_dD symbolically
    for i in range(DIM):
        alpha_dD[i] = sp.diff(alpha, coords[i])

    # Define betaU = 0
    betaU = ixp.zerorank1()

    # Define betaU_dD = 0
    betaU_dD = ixp.zerorank2()

    # This function replaces the spherical coordinates r, th, ph for the numerical grid xx0, xx1, xx2
    # Taken from BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear.py
    def sympify_integers__replace_rthph_or_Cartxyz(obj, rthph_or_xyz, rthph_or_xyz_of_xx):
        if isinstance(obj, int):
            return sp.sympify(obj)
        else:
            return obj.subs(rthph_or_xyz[0], rthph_or_xyz_of_xx[0]).\
                subs(rthph_or_xyz[1], rthph_or_xyz_of_xx[1]).\
                subs(rthph_or_xyz[2], rthph_or_xyz_of_xx[2])

    # Replace spherical coordinates for computational grid variables everywhere
    trK = sympify_integers__replace_rthph_or_Cartxyz(trK, coords, rfm.xxSph)
    psi = sympify_integers__replace_rthph_or_Cartxyz(psi, coords, rfm.xxSph)
    phi = sympify_integers__replace_rthph_or_Cartxyz(phi, coords, rfm.xxSph)
    exp_m4phi = sympify_integers__replace_rthph_or_Cartxyz(
        exp_m4phi, coords, rfm.xxSph)
    alpha = sympify_integers__replace_rthph_or_Cartxyz(
        alpha, coords, rfm.xxSph)
    for i in range(DIM):
        alpha_dD[i] = sympify_integers__replace_rthph_or_Cartxyz(
            alpha_dD[i], coords, rfm.xxSph)
        betaU[i] = sympify_integers__replace_rthph_or_Cartxyz(
            betaU[i], coords, rfm.xxSph)
        phi_dD[i] = sympify_integers__replace_rthph_or_Cartxyz(
            phi_dD[i], coords, rfm.xxSph)
        for j in range(DIM):
            gammabarDD[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                gammabarDD[i][j], coords, rfm.xxSph)
            gammabarUU[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                gammabarUU[i][j], coords, rfm.xxSph)
            betaU_dD[i][j] = sympify_integers__replace_rthph_or_Cartxyz(
                betaU_dD[i][j], coords, rfm.xxSph)
            for k in range(DIM):
                GammabarUDD[i][j][k] = sympify_integers__replace_rthph_or_Cartxyz(
                    GammabarUDD[i][j][k], coords, rfm.xxSph)


def metric_quantities():

    metric = par.parval_from_str(thismodule + "::BackgroundMetric")

    if metric == "Minkowski":
        flat_metric_quantities()
    elif metric == "Schwarzschild":
        schwarzschild_metric_quantities()
    elif metric == "Kerr":
        kerr_metric_quantities()
    else:
        print(f"Error: Background Metric {metric} unsupported!")
        print(
            f"The suported background metrics by the {thismodule} module are:")
        print("\t - Minkowski (flat spacetime)")
        print("\t - Schwarzschild (spherically symmetric spacetime)")
        # print("\t - Kerr (rotating spacetime)")
        print("Please choose one of the above.")
        sys.exit(1)


def print_gmq_quantities():

    DIM = 3

    print('Printing metric quantities...')
    print(f'trK = {trK}')
    print(f'psi = {psi}')
    print(f'phi = {phi}')
    print(f'exp_m4phi = {exp_m4phi}')
    print(f'alpha = {alpha}')
    print('Printing alpha_dD...')
    for i in range(DIM):
        print(f'alpha_dD[{i}] = {alpha_dD[i]}')
    print('Printing betaU...')
    for i in range(DIM):
        print(f'betaU[{i}] = {betaU[i]}')
    print('Printing gammabarDD...')
    for i in range(DIM):
        for j in range(DIM):
            print(f'gammabarDD[{i}][{j}] = {gammabarDD[i][j]}')
    print('Printing gammabarUU...')
    for i in range(DIM):
        for j in range(DIM):
            print(f'gammabarUU[{i}][{j}] = {gammabarUU[i][j]}')
    print('Printing betaU_dD...')
    for i in range(DIM):
        for j in range(DIM):
            print(f'betaU[{i}][{j}] = {betaU_dD[i][j]}')
    print('Printing GammabarUDD...')
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                print(f'GammabarUDD[{i}][{j}][{k}] = {GammabarUDD[i][j][k]}')


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
