import NRPy_param_funcs as par    # NRPy+: Parameter interface
# SymPy: The Python computer algebra package upon which NRPy+ depends
import sympy as sp
# NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import indexedexp as ixp
import reference_metric as rfm    # NRPy+: Reference metric support
# Standard Python modules for multiplatform OS-level functions
import ScalarFieldStaticBackground.gridfunctions_and_metric_quantities as gmq
import sys

thismodule = __name__


def field_RHSs():

    # Set spatial dimension to 3
    DIM = 3

    # Call reference metric function if it has not been called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Call functions to define gridfunctions and establish related metric and gauge quantities
    # gmq.betaU_deriv()
    gmq.flat_metric_quantities()

    # Rename variables
    Phi = gmq.Phi
    Pi = gmq.Pi
    gammabarUU = gmq.gammabarUU
    GammabarUDD = gmq.GammabarUDD
    exp_m4phi = gmq.exp_m4phi
    phi_dD = gmq.phi_dD
    trK = gmq.trK
    alpha = gmq.alpha
    alpha_dD = gmq.alpha_dD
    betaU = gmq.betaU

    # Undefined derivatives
    Phi_dD = ixp.declarerank1("Phi_dD")
    Phi_dDD = ixp.declarerank2("Phi_dDD", "sym01")
    Phi_dupD = ixp.declarerank1("Phi_dupD")
    Pi_dupD = ixp.declarerank1("Pi_dupD")

    # Scalar field mass
    mu_s = par.Cparameters("REAL", thismodule, ["mu_s"], 1.0)

    # Variables for the RHS expressions
    global Phi_rhs, Pi_rhs

    # RHS for the field Phi
    # \partial_t \Phi = -\alpha \Pi + \beta^i \partial_i \Phi
    Phi_rhs = - alpha * Pi
    for i in range(DIM):
        # Shift advection derivatives are upwinded
        Phi_rhs += betaU[i] * Phi_dupD[i]

    # RHS for the conjugate momentum Pi
    # \partial_t \Pi = \alpha (- e^{-4 \phi} \bar{\gamma}^{ij} \partial_i \partial_j \Phi +
    # e^{-4 \phi} \bar{\gamma}^{ij} \bar{\Gamma}^k_{ij} \partial_k \Phi -
    # 2 e^{-4 \phi} \bar{\gamma}^{ij} \partial_j \Phi \partial_i \phi + K \Pi + \mu_s^2 \Phi ) -
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
