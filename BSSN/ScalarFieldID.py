import indexedexp as ixp
import reference_metric as rfm
import NRPy_param_funcs as par
import sympy as sp

thismodule = __name__
A, r0, w = par.Cparameters("REAL", thismodule, ["A", "r0", "w"], [0.15, 5.0, 0.5])
omega = par.Cparameters("REAL", thismodule, ["omega"], [0.3929])
m = par.Cparameters("int", thismodule, ["m"], [1]);

# This function replaces the spherical coordinates r, th, ph for the numerical grid xx0, xx1, xx2
# Taken from BSSN.ADM_Exact_Spherical_or_Cartesian_to_BSSNCurvilinear.py
def sympify_integers__replace_rthph_or_Cartxyz(obj, rthph_or_xyz, rthph_or_xyz_of_xx):
    if isinstance(obj, int):
        return sp.sympify(obj)
    else:
        return obj.subs(rthph_or_xyz[0], rthph_or_xyz_of_xx[0]).\
            subs(rthph_or_xyz[1], rthph_or_xyz_of_xx[1]).\
            subs(rthph_or_xyz[2], rthph_or_xyz_of_xx[2])


def ScalarFieldID():

    # Declare global variables for the ID expressions of the fields
    global SphPhi, SphPi

    # Call reference metric if it hasn't been called called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Define the spherical coordinate variables
    r, th, ph , t = sp.symbols("r th ph t", real=True)
    coords = [r, th, ph]

    f = sp.Rational(1/2) * sp.exp(-(r + t - r0)**2 / w**2)*(r+t)
    g = sp.Rational(1/2) * sp.exp(-(r - t + r0)**2 / w**2)*(r-t)

    Phi_t = A / r * (f + g)
    SphPhi = sp.simplify(Phi_t.subs(t, sp.sympify(0)))

    Pi_t = sp.diff(Phi_t, t)
    SphPi = - Pi_t.subs(t, sp.sympify(0))

    SphPhi = sympify_integers__replace_rthph_or_Cartxyz(SphPhi, coords, rfm.xxSph)
    SphPi = sympify_integers__replace_rthph_or_Cartxyz(SphPi, coords, rfm.xxSph)

    import BSSN.ScalarField_ID_function_string as sfIDf
    global returnfunction
    returnfunction = sfIDf.scalar_field_ID_function_string(SphPhi, SphPi)


def ScalarFieldID_Schwarzschild(psi, ID_Type="Gaussian"):

    # Declare global variables for the ID expressions of the fields
    global SphPhi, SphPi

    # Call reference metric if it hasn't been called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Define the spherical coordinate variables
    r, th, ph = sp.symbols("r th ph", real=True)
    coords = [r, th, ph]

    # Set the F(r) and Z(\theta, \varphi) functions depending on ID_Type
    if ID_Type == "Gaussian":
        F = A * sp.sqrt(r) * sp.exp(- (r - r0)**2 / w**2)
        Z = 1 / sp.sqrt(4 * sp.pi)
    elif ID_Type == "Dipole":
        F = A * r * sp.exp(- (r - r0)**2 / w**2)
        Z = sp.sqrt(3 / (2 * sp.pi)) * sp.sin(th) * sp.cos(ph)
    else:
        print(f"ID_Type = {ID_Type} not supported!")
        exit(1)

    # Set the Phi initial data to zero
    SphPhi = 0
    # Set the Pi initial data to its expression
    SphPi = psi**(-sp.Rational(5, 2)) / sp.sqrt(r * sp.pi) * F * Z

    SphPhi = sympify_integers__replace_rthph_or_Cartxyz(
        SphPhi, coords, rfm.xxSph)
    SphPi = sympify_integers__replace_rthph_or_Cartxyz(
        SphPi, coords, rfm.xxSph)

    import BSSN.ScalarField_ID_function_string as sfIDf
    global returnfunction
    returnfunction = sfIDf.scalar_field_ID_function_string(SphPhi, SphPi)


def ScalarFieldID_QuasiBound(alpha, beta, coords):

    # Declare global variables for the ID expressions of the fields
    global SphPhi, SphPi

    # Call reference metric if it hasn't been called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Define the spherical coordinate variables and time variable
    t = sp.symbols("t", real=True)
    r, th, ph = coords

    Phi_t = A / sp.sqrt(sp.pi) * sp.exp(- (r - r0)**2 / w**2) * sp.cos(omega * t + m * ph) * sp.sin(th)
    Pi_t = 1 / alpha * (beta * sp.diff(Phi_t, ph) - sp.diff(Phi_t, t))

    SphPhi = Phi_t.subs(t, 0)
    SphPi  = Pi_t.subs(t, 0)

    SphPhi = sympify_integers__replace_rthph_or_Cartxyz(
        SphPhi, coords, rfm.xxSph)
    SphPi = sympify_integers__replace_rthph_or_Cartxyz(
        SphPi, coords, rfm.xxSph)

    import BSSN.ScalarField_ID_function_string as sfIDf
    global returnfunction
    returnfunction = sfIDf.scalar_field_ID_function_string(SphPhi, SphPi)
