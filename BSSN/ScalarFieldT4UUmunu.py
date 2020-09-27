import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import reference_metric as rfm
import BSSN.BSSN_quantities as Bq
import BSSN.ADM_in_terms_of_BSSN as AtoB
import BSSN.ScalarFieldRHSs as sfrhs
import BSSN.ADMBSSN_tofrom_4metric as abtfm
import sympy as sp

thismodule = __name__


def ScalarFieldT4UU():

    # Set spacetime dimension to 4
    DIM = 4
    
    # Call reference metric if it hasn't been called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Define the inverse 4-metric tensor
    abtfm.g4UU_ito_BSSN_or_ADM("BSSN")
    g4UU = abtfm.g4UU

    # Declare stress-energy tensor global variable
    global T4UU
    T4UU = ixp.zerorank2(DIM=DIM)

    # Obtain Phi, its RHS expression and mass
    Phi = sfrhs.Phi
    Phi_rhs = sfrhs.Phi_rhs
    mu_s = sfrhs.mu_s

    # Declare the finite difference derivative of Phi
    Phi_dD = ixp.declarerank1("Phi_dD")

    # Define the stress-energy tensor expression

    for mu in range(DIM):
        for nu in range(DIM):
            for rho in range(DIM):
                for sigma in range(DIM):

                    if rho == 0:
                        rho_deriv = Phi_rhs
                    else:
                        rho_deriv = Phi_dD[rho - 1]

                    if sigma == 0:
                        sig_deriv = Phi_rhs
                    else:
                        sig_deriv = Phi_dD[sigma - 1]

                    contraction = g4UU[rho][sigma] * rho_deriv * sig_deriv

            T4UU[mu][nu] = - sp.Rational(1, 2) * g4UU[mu][nu] * \
                (contraction + mu_s**2 * Phi**2)

    for mu in range(DIM):
        for nu in range(DIM):
            for alpha in range(DIM):
                for beta in range(DIM):

                    if alpha == 0:
                        alpha_deriv = Phi_rhs
                    else:
                        alpha_deriv = Phi_dD[alpha - 1]
                    
                    if beta == 0:
                        beta_deriv = Phi_rhs
                    else:
                        beta_deriv = Phi_dD[beta - 1]

                    T4UU[mu][nu] += sp.Rational(1,2) * g4UU[alpha][mu] * g4UU[beta][nu] * (2 * alpha_deriv * beta_deriv)

def ScalarFieldSourceTerms():

    # Set spatial dimension to 3
    DIM = 3

    # Define PI
    PI = par.Cparameters("REAL", thismodule, [
                         "PI"], "3.14159265358979323846264338327950288")

    # Call reference metric if it hasn't been called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Obtain the BSSN and field variables from other modules
    AtoB.ADM_in_terms_of_BSSN()
    gammaDD = AtoB.gammaDD
    gammaUU = AtoB.gammaUU
    Phi = sfrhs.Phi
    Pi = sfrhs.Pi
    mu_s = sfrhs.mu_s

    # Declare the spatial derivatives of Phi
    Phi_dD = ixp.declarerank1("Phi_dD")

    # Declare the BSSN source terms
    global SDD, SD, S, rho

    SDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            SDD[i][j] += sp.Rational(1,2) * gammaDD[i][j] * (Pi**2 - mu_s**2 * Phi**2) + Phi_dD[i] * Phi_dD[j]
            for k in range(DIM):
                for l in range(DIM):
                    SDD[i][j] += -sp.Rational(1,2) * gammaDD[i][j] * gammaUU[k][l] * Phi_dD[k] * Phi_dD[l]
    
    SD = ixp.zerorank1()
    for i in range(DIM):
        SD[i] += - Pi * Phi_dD[i]

    S = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            S += gammaUU[i][j] * SDD[i][j]

    rho = -sp.Rational(1,2) * (Pi**2 + mu_s**2 * Phi**2)
    for k in range(DIM):
        for l in range(DIM):
            rho += -sp.Rational(1,2) * gammaUU[k][l] * Phi_dD[k] * Phi_dD[l]

    # Declare global variables for source terms
    global sourceterm_trK_rhs, sourceterm_Lambdabar_rhsU, sourceterm_lambda_rhsU, sourceterm_a_rhsDD

    # Obtain the inverse conformal metric as the matrix inverse of the conformal metric
    Bq.BSSN_basic_tensors()
    gammabarUU, _ = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)

    # Define exp_m4phi and the lapse
    Bq.phi_and_derivs()
    exp_m4phi = Bq.exp_m4phi
    alpha = Bq.alpha

    # Compute the trace-free part of SDD
    tracefree_SDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            tracefree_SDD[i][j] = SDD[i][j]
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    tracefree_SDD[i][j] += -sp.Rational(1,3) * Bq.gammabarDD[i][j] * gammabarUU[k][l] * SDD[k][l]

    # Source term for trK
    sourceterm_trK_rhs = 4 * PI * alpha * (rho + S)

    # Source term for AbarDD
    sourceterm_a_rhsDD  = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            Abar_rhsDDij = -8 * PI * alpha * exp_m4phi * tracefree_SDD[i][j]
            sourceterm_a_rhsDD[i][j] = Abar_rhsDDij / rfm.ReDD[i][j]

    # Source term for Lambdabar (needed for betaU rhs)
    sourceterm_Lambdabar_rhsU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            sourceterm_Lambdabar_rhsU[i] += -16 * PI * alpha * gammabarUU[i][j] * SD[j]
    
    # Source term for lambda
    sourceterm_lambda_rhsU = ixp.zerorank1()
    for i in range(DIM):
        sourceterm_lambda_rhsU[i] = sourceterm_lambda_rhsU[i] / rfm.ReU[i]

    # Declare the source term for the hamiltonian constraint
    global sourceterm_H
    sourceterm_H = -16 * PI * rho
