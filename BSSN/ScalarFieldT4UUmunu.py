import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import reference_metric as rfm
import BSSN.BSSN_quantities as Bq
import BSSN.ADM_in_terms_of_BSSN as AtoB
import BSSN.ScalarFieldRHSs as sfrhs
import sympy as sp

thismodule = __name__


def ScalarFieldT4UU():

    # Set spacetime dimension to 4
    DIM = 4

    # Call reference_metric() if it hasn't been called already
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()
    
    # Declare and define all needed BSSN quantities
    Bq.BSSN_basic_tensors()

    # Rename needed BSSN quantities
    alpha = Bq.alpha
    betaU = Bq.betaU

    # Obtain the physical inverse 3-metric from the ADM_in_terms_of_BSSN module
    AtoB.ADM_in_terms_of_BSSN()
    gammaUU = AtoB.gammaUU

    # Declare the inverse 4-metric
    g4UU = ixp.zerorank2(DIM=DIM)

    # Obtain the correct expressions for g4UU
    g4UU[0][0] = 1 / alpha**2
    for i in range(3):
        g4UU[0][i] = g4UU[i][0] = betaU[i] / alpha**2
        for j in range(1, DIM):
            g4UU[i][j] = gammaUU[i-1][j-1] - (betaU[i-1] * betaU[j-1]) / alpha**2

    # Declare stress-energy tensor
    global T4UU
    T4UU = ixp.zerorank2(DIM=DIM)

    # Obtain Phi, the RHS expression for Phi and the scalar field mass from the ScalarFieldRHSs module
    Phi = sfrhs.Phi
    Phi_rhs = sfrhs.Phi_rhs
    mu_s = sfrhs.mu_s

    # Declare the derivative of Phi (finite differencing)
    Phi_dD = ixp.declarerank1("Phi_dD")

    # Define the stress-energy tensor symbolic expression
    # T^{\mu \nu} = - \frac{1}{2} g^{\mu \nu} \left[ g^{00} \left( \partial_t \Phi \right)^2 
    # + 2 g^{0i} \partial_t \Phi \partial_i \Phi + g^{ij} \partial_i \Phi \partial_j \Phi + \mu_s^2 \Phi^2 \right] 
    # + g^{0\mu} g^{0\nu} \left( \partial_t \Phi \right)^2
    # + \left( g^{0\mu} g^{i\nu} + g^{i\mu} g^{0\nu} \right) \partial_t \Phi \partial_i \Phi 
    # + g^{i\mu} g^{j\nu} \partial_i \Phi \partial_j \Phi
    # Where \partial_t \Phi appears, we use the RHS expression for Phi

    # First term inside []: g^{00} \left( \partial_t \Phi \right)^2
    braces = g4UU[0][0] * Phi_rhs**2

    # Second term inside []: 2 g^{0i} \partial_t \Phi \partial_i \Phi
    for i in range(1, DIM):
        braces += 2 * g4UU[0][i] * Phi_rhs * Phi_dD[i-1]

    # Third term inside []: g^{ij} \partial_i \Phi \partial_j \Phi
    for i in range(1, DIM):
        for j in range(1, DIM):
            braces += g4UU[i][j] * Phi_dD[i-1] * Phi_dD[j-1]

    # Forth term inside []: \mu_s^2 \Phi^2
    braces += mu_s**2 * Phi**2

    # Multiply term inside [] by - \frac{1}{2} g^{\mu \nu} and assign to T^{\mu \nu}
    for mu in range(DIM):
        for nu in range(DIM):
            T4UU[mu][nu] = - sp.Rational(1, 2) * g4UU[mu][nu] * braces

    # Add first term after []: g^{0\mu} g^{0\nu} \left( \partial_t \Phi \right)^2
    for mu in range(DIM):
        for nu in range(DIM):
            T4UU[mu][nu] += g4UU[0][mu] * g4UU[0][nu] * Phi_rhs**2

    # Add second term after []: \left( g^{0\mu} g^{i\nu} + g^{i\mu} g^{0\nu} \right) \partial_t \Phi \partial_i \Phi
    for mu in range(DIM):
        for nu in range(DIM):
            for i in range(1, DIM):
                T4UU[mu][nu] += (g4UU[0][mu] * g4UU[i][nu] + g4UU[i][mu] * g4UU[0][nu]) * Phi_rhs * Phi_dD[i-1] 
    
    # Add final term after []: g^{i\mu} g^{j\nu} \partial_i \Phi \partial_j \Phi
    for mu in range(DIM):
        for nu in range(DIM):
            for i in range(1, DIM):
                for j in range(1, DIM):
                    T4UU[mu][nu] += g4UU[i][mu] * g4UU[j][nu] * Phi_dD[i-1] * Phi_dD[j-1]