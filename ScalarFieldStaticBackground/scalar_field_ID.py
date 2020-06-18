import indexedexp as ixp
import reference_metric as rfm
import NRPy_param_funcs as par
import ScalarFieldStaticBackground.gridfunctions_and_metric_quantities as gmq
import sympy as sp
import sys

thismodule = __name__
# par.set_parval_from_str("reference_metric::CoordSystem", "Spherical")
par.initialize_param(par.glb_param("char", thismodule,
                                   "IDStateType", "SphericalGaussian"))
A, r0, w = par.Cparameters(
    "REAL", thismodule, ["A", "r0", "w"], [1.0, 12.0, 2.0])


def scalar_field_ID():

    # Define global variables for the initial data
    global Phi, Pi

    # Set spatial dimension to 3
    # DIM = 3

    # Call reference metric if it hasn't been called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Get coordinate system and see if it is of spherical type
    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    if not "Spherical" in CoordSystem:
        print("Error: CoordSystem = " + CoordSystem + " unsupported!")
        sys.exit(1)

    # Define spherical coordinate variables
    r, th, ph = sp.symbols("r th ph", real=True)
    coords = [r, th, ph]

    # Define the conformal factor psi
    gmq.flat_metric_quantities()
    psi = gmq.psi

    # Set up the initial data for the field depending on IDStateType
    IDStateType = par.parval_from_str(thismodule + "::IDStateType")
    if IDStateType == "SphericalGaussian":

        StateTypeCparam = 0

        # Define the radial function F and angular function Z to compute Pi
        F = A * sp.sqrt(r) * sp.exp(-(r - r0)**2 / w**2)
        Z = 1 / sp.sqrt(4 * sp.pi)

    elif IDStateType == "DipoleGaussian":

        StateTypeCparam = 1

        # Define the radial function F and angular function Z to compute Pi
        F = A * r * sp.exp(-(r - r0)**2 / w**2)
        Z = sp.simplify(
            (sp.Ynm(1, -1, th, ph) - sp.Ynm(1, 1, th, ph)).expand(func=True))

    else:
        print("Error: IDStateType = " + IDStateType + " unsupported!")
        sys.exit(1)

    # Define C parameter to identify IDStateType
    _ = par.Cparameters("int", thismodule, "IDStateType", StateTypeCparam)

    # Define the initial data for Phi and Pi from psi, F and Z
    Phi = sp.sympify(0)
    Pi = psi**sp.Rational(-5, 2) / sp.sqrt(r * sp.pi) * F * Z

    # Define the gauge quantities and set initial data
    # alpha = 1 / psi**2
    # betaU = ixp.zerorank1()
    # BU = ixp.zerorank1()

    # This function replaces the spherical coordinates r, th, ph for the numerical grid xx0, xx1, xx2
    # Taken from BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear.py

    def sympify_integers__replace_rthph_or_Cartxyz(obj, rthph_or_xyz, rthph_or_xyz_of_xx):
        if isinstance(obj, int):
            return sp.sympify(obj)
        else:
            return obj.subs(rthph_or_xyz[0], rthph_or_xyz_of_xx[0]).\
                subs(rthph_or_xyz[1], rthph_or_xyz_of_xx[1]).\
                subs(rthph_or_xyz[2], rthph_or_xyz_of_xx[2])

    # Substitute the spherical coordinates by the computational grid coordinates
    Phi = sympify_integers__replace_rthph_or_Cartxyz(Phi, coords, rfm.xxSph)
    Pi = sympify_integers__replace_rthph_or_Cartxyz(Pi, coords, rfm.xxSph)
    # alpha = sympify_integers__replace_rthph_or_Cartxyz(
    #     alpha, coords, rfm.xxSph)
    # for i in range(DIM):
    #     betaU[i] = sympify_integers__replace_rthph_or_Cartxyz(
    #         betaU[i], coords, rfm.xxSph)
    #     BU[i] = sympify_integers__replace_rthph_or_Cartxyz(
    #         BU[i], coords, rfm.xxSph)

    # Rescale the shift variables
    # vetU = ixp.zerorank1()
    # betU = ixp.zerorank1()
    # for i in range(DIM):
    #     vetU[i] += betaU[i] / rfm.ReU[i]
    #     betU[i] += BU[i] / rfm.ReU[i]

    # Add initial data function to C functions dictionary
    import ScalarFieldStaticBackground.scalar_field_ID_function_string as sfIDf
    global returnfunction
    returnfunction = sfIDf.scalar_field_ID_function_string(
        Phi, Pi)


def scalar_field_test_ID():

    global Phi, Pi
    DIM = 3

    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    CoordSystem = par.parval_from_str("reference_metric::CoordSystem")
    if not "Spherical" in CoordSystem:
        print("Error: CoordSystem = " + CoordSystem + " unsupported!")
        sys.exit(1)

    r, th, ph, t = sp.symbols("r th ph t", real=True)
    coords = [r, th, ph]

    gmq.flat_metric_quantities()
    alpha = gmq.alpha
    betaU = gmq.betaU

    f = sp.Rational(1/2) * sp.exp(-(r + t - r0)**2 / w**2)
    g = sp.Rational(1/2) * sp.exp(-(r - t + r0)**2 / w**2)

    Phi_t = 1/r * (f + g)
    Phi = sp.simplify(Phi_t.subs(t, sp.sympify(0)))

    Pi_t = sp.diff(Phi_t, t)
    for i in range(DIM):
        Pi_t -= betaU[i] * sp.diff(Phi_t, coords[i])
    Pi = - 1/alpha * Pi_t.subs(t, sp.sympify(0))
    
    # This function replaces the spherical coordinates r, th, ph for the numerical grid xx0, xx1, xx2
    # Taken from BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear.py
    def sympify_integers__replace_rthph_or_Cartxyz(obj, rthph_or_xyz, rthph_or_xyz_of_xx):
        if isinstance(obj, int):
            return sp.sympify(obj)
        else:
            return obj.subs(rthph_or_xyz[0], rthph_or_xyz_of_xx[0]).\
                subs(rthph_or_xyz[1], rthph_or_xyz_of_xx[1]).\
                subs(rthph_or_xyz[2], rthph_or_xyz_of_xx[2])

    # Substitute the spherical coordinates by the computational grid coordinates
    Phi = sympify_integers__replace_rthph_or_Cartxyz(Phi, coords, rfm.xxSph)
    Pi = sympify_integers__replace_rthph_or_Cartxyz(Pi, coords, rfm.xxSph)
    # alpha = sympify_integers__replace_rthph_or_Cartxyz(
    #     alpha, coords, rfm.xxSph)
    # for i in range(DIM):
    #     betaU[i] = sympify_integers__replace_rthph_or_Cartxyz(
    #         betaU[i], coords, rfm.xxSph)
    #     BU[i] = sympify_integers__replace_rthph_or_Cartxyz(
    #         BU[i], coords, rfm.xxSph)

    # Rescale the shift variables
    # vetU = ixp.zerorank1()
    # betU = ixp.zerorank1()
    # for i in range(DIM):
    #     vetU[i] += betaU[i] / rfm.ReU[i]
    #     betU[i] += BU[i] / rfm.ReU[i]

    # Add initial data function to C functions dictionary
    import ScalarFieldStaticBackground.scalar_field_ID_function_string as sfIDf
    global returnfunction
    returnfunction = sfIDf.scalar_field_ID_function_string(Phi, Pi)
