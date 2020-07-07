import NRPy_param_funcs as par
import indexedexp as ixp
import reference_metric as rfm
import ScalarFieldStaticBackground.gridfunctions_and_metric_quantities as gmq
import sympy as sp
import sys

thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "LapseEvolutionOption", "OnePlusLog"))
par.initialize_param(par.glb_param("char", thismodule, "ShiftEvolutionOption", "GammaDriving2ndOrder_Covariant"))


def gauge_RHSs():

    # Set spatial dimension to 3
    DIM = 3

    # Set up the reference metric and related quantitites
    # needed for the rescaling of the right hand sides
    if not rfm.have_already_called_reference_metric_function:
        rfm.reference_metric()

    # Call functions to define gridfunctions and establish metric and gauge quantities
    gmq.betaU_deriv()
    gmq.flat_metric_quantities()

    # Lapse evolution equation RHS

    # Get lapse evolution option
    LapseEvolutionOption = par.parval_from_str(thismodule + "::LapseEvolutionOption")

    # Get expressions for metric and gauge quantities
    alpha = gmq.alpha
    trK = gmq.trK
    betaU = gmq.betaU

    # Define the global quantity to store the RHS of \alpha
    global alpha_rhs

    # Compute RHS depending on the lapse evolution option
    # Only 1 + log slicing supported, but can be easily extended
    if LapseEvolutionOption == "OnePlusLog":
        
        # Define the upwinded derivative of \alpha
        alpha_dupD = ixp.declarerank1("alpha_dupD")
        
        # Compute RHS
        alpha_rhs = -2 * alpha * trK
        for i in range(DIM):
            alpha_rhs += betaU[i] * alpha_dupD[i]

    else:
        print("Error: LapseEvolutionOption = " + LapseEvolutionOption + " unsupported!")
        sys.exit()

    # Shift evolution equations RHSs

    # Get shift evolution option
    ShiftEvolutionOption = par.parval_from_str(thismodule + "::ShiftEvolutionOption")

    # Get expressions for metric and gauge quantities
    BU = gmq.BU
    betU = gmq.betU
    GammabarUDD = gmq.GammabarUDD
    LambdabarU = gmq.LambdabarU
    LambdabarU_dD = gmq.LambdabarU_dD
    Lambdabar_rhsU = ixp.zerorank1()    # Zero because the spacetime doesn't change
    betaU_dupD = gmq.betaU_dupD

    # Generate expressions for the RHSs and set them to 0
    beta_rhsU = ixp.zerorank1()
    B_rhsU = ixp.zerorank1()

    # Compute expressions for the RHSs depending on the shift evolution option
    # Only covariant 2nd order gamma driver is supported, but it can be easily extended
    if ShiftEvolutionOption == "GammaDriving2ndOrder_Covariant":

        # RHS for B^i

        # Define needed derivatives
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD", "nosym")
        for i in range(DIM):
            for j in range(DIM):
                BU_dupD[i][j] += betU_dupD[i][j] * rfm.ReU[i] + betU[i] * rfm.ReUdD[i][j]

        # Term 1: \beta^j \partial_j B^i
        for i in range(DIM):
            for j in range(DIM):
                B_rhsU[i] += betaU[j] * BU_dupD[i][j]

        # Term 2: \beta^j \bar{\Gamma}^i_{kj} B^k
        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    B_rhsU[i] += betaU[j] * GammabarUDD[i][k][j] * BU[k]
                
        # Term 3: \frac{3}{4} \partial_t \bar{\Lambda}^i
        for i in range(DIM):
            B_rhsU[i] += sp.Rational(3,4) * Lambdabar_rhsU[i]

        # Term 4: - \frac{3}{4} \beta^j \partial_j \bar{\Lambda}^i
        for i in range(DIM):
            for j in range(DIM):
                B_rhsU[i] += - sp.Rational(3,4) * betaU[j] * LambdabarU_dD[i][j]

        # Term 5: - \frac{3}{4} \beta^j \bar{\Gamma}^i_{kj} \bar{\Lambda}^k
        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    B_rhsU[i] += - sp.Rational(3,4) * GammabarUDD[i][k][j] * LambdabarU[k]

        # Term 6: - \eta B^i
        # eta is a free parameter; we declare it here
        eta = par.Cparameters("REAL", thismodule, ["eta"], 2.0)
        for i in range(DIM):
            B_rhsU[i] += eta * BU[i]

        # RHS for \beta^i

        # Term 1: \beta^j \partial_j \beta^i
        for i in range(DIM):
            for j in range(DIM):
                beta_rhsU[i] += betaU[j] * betaU_dupD[i][j]

        # Term 2: \beta^j \bar{\Gamma}^i_{kj} \beta^k
        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    beta_rhsU[i] += betaU[j] * GammabarUDD[i][k][j] * betaU[k]

        # Term 3: B^i
        for i in range(DIM):
            beta_rhsU[i] += BU[i]

    else:
        print("Error: ShiftEvolutionOption = " + ShiftEvolutionOption + " unsupported!")
        sys.exit(1)

    # Rescaling of the RHSs
    global vet_rhsU, bet_rhsU

    vet_rhsU = ixp.zerorank1()
    bet_rhsU = ixp.zerorank1()

    for i in range(DIM):
        vet_rhsU[i] = beta_rhsU[i] / rfm.ReU[i]
        bet_rhsU[i] = B_rhsU[i] / rfm.ReU[i]
    