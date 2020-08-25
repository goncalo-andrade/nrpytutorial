import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import BSSN.BSSN_quantities as Bq
import sympy as sp

thismodule = __name__


def declare_field_gridfunctions_if_not_declared_already():

    # Declare global variables for the fields
    global Phi, Pi

    # Check if gridfunctions were declared already
    # You can't declare a gridfunction twice
    for i in range(len(gri.glb_gridfcs_list)):
        if "Phi" in gri.glb_gridfcs_list[i].name:
            Phi, Pi = sp.symbols("Phi Pi", real=True)
            return Phi, Pi

    # Declare the gridfunctions
    Phi, Pi = gri.register_gridfunctions("EVOL", ["Phi", "Pi"])
    return Phi, Pi


def ScalarFieldRHSs():

    # Declare the field gridfunctions
    declare_field_gridfunctions_if_not_declared_already()

    # Set spatial dimension to 3
    DIM = 3

    # Call Bq.phi_and_derivs() to define all the quantities we need
    # It automatically calls Bq.BSSN_basic_tensors() and Bq.gammabar__inverse_and_derivs
    Bq.phi_and_derivs()

    # Define all the needed quantities
    gammabarUU = Bq.gammabarUU
    GammabarUDD = Bq.GammabarUDD
    alpha = Bq.alpha
    betaU = Bq.betaU
    trK = Bq.trK
    exp_m4phi = Bq.exp_m4phi
    phi_dD = Bq.phi_dD

    # Define some needed derivatives of the above quantities
    alpha_dD = ixp.declarerank1("alpha_dD")
    Phi_dD = ixp.declarerank1("Phi_dD")
    Phi_dDD = ixp.declarerank2("Phi_dDD", "sym01")
    # Upstream derivatives for contraction with the shift vector
    Phi_dupD = ixp.declarerank1("Phi_dupD")
    Pi_dupD = ixp.declarerank1("Pi_dupD")

    # Define the scalar field mass as a C parameter
    global mu_s
    mu_s = par.Cparameters("REAL", thismodule, ["mu_s"], 0.0)

    # Declare global variables for the RHS expressions
    global Phi_rhs, Pi_rhs

    # Compute the RHS expressions

    # RHS for the field Phi
    # \partial_t \Phi = -\alpha \Pi + \beta^i \partial_i \Phi
    Phi_rhs = - alpha * Pi
    for i in range(DIM):
        # Shift advection derivatives are upwinded
        Phi_rhs += betaU[i] * Phi_dupD[i]

    # RHS for the conjugate momentum Pi
    # \partial_t \Pi = \alpha \left(- e^{-4 \phi} \bar{\gamma}^{ij} \partial_i \partial_j \Phi +
    # e^{-4 \phi} \bar{\gamma}^{ij} \bar{\Gamma}^k_{ij} \partial_k \Phi -
    # 2 e^{-4 \phi} \bar{\gamma}^{ij} \partial_j \Phi \partial_i \phi + K \Pi + \mu_s^2 \Phi \right) -
    # e^{-4 \phi} \bar{\gamma}^{ij} \partial_j \alpha \partial_i \Phi + \beta^i \partial_i \Pi

    # Scalar part inside braces: K \Pi + \mu_s^2 \Phi
    term_in_braces = trK * Pi + mu_s * mu_s * Phi

    # Sums inside braces

    # First sum: - e^{-4 \phi} \bar{\gamma}^{ij} \partial_i \partial_j \Phi
    for i in range(DIM):
        for j in range(DIM):
            term_in_braces += -exp_m4phi * \
                gammabarUU[i][j] * Phi_dDD[i][j]

    # Second sum: e^{-4 \phi} \bar{\gamma}^{ij} \bar{\Gamma}^k_{ij} \partial_k \Phi
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                term_in_braces += exp_m4phi * \
                    gammabarUU[i][j] * GammabarUDD[k][i][j] * Phi_dD[k]

    # Third sum: - 2 e^{-4 \phi} \bar{\gamma}^{ij} \partial_j \Phi \partial_i \phi
    for i in range(DIM):
        for j in range(DIM):
            term_in_braces += -2 * exp_m4phi * \
                gammabarUU[i][j] * Phi_dD[j] * phi_dD[i]

    # Term in braces complete
    # Now multiply by alpha and assign to the RHS variable
    Pi_rhs = alpha * term_in_braces

    # First term after braces: - e^{-4 \phi} \bar{\gamma}^{ij} \partial_j \alpha \partial_i \Phi
    for i in range(DIM):
        for j in range(DIM):
            Pi_rhs += - exp_m4phi * gammabarUU[i][j] * alpha_dD[j] * Phi_dD[i]

    # Second term after braces: \beta^i \partial_i \Pi
    for i in range(DIM):
        # Shift advection derivatives are upwinded
        Pi_rhs += betaU[i] * Pi_dupD[i]
