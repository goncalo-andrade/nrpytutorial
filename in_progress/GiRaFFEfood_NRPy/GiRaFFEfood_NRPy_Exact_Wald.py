#!/usr/bin/env python
# coding: utf-8

# <a id='top'></a>
# # $\texttt{GiRaFFEfood}$: Initial data for $\texttt{GiRaFFE}$
# 
# With the $\texttt{GiRaFFE}$ evolution thorn constructed, we now need to "feed" our giraffe with initial data to evolve. While there are several different choices of initial data we can use here, for the moment, we will only be implementing the "Exact Wald" initial data, given by Table 3 in [the original paper](https://arxiv.org/pdf/1704.00599.pdf):
# \begin{align}
# A_{\phi} &= \frac{C_0}{2} r^2 \sin^2 \theta \\
# E_{\phi} &= 2 M C_0 \left( 1+ \frac {2M}{r} \right)^{-1/2} \sin^2 \theta \\
# \end{align}
# (the unspecified components are set to 0). Here, $C_0$ is a constant that we will set to $1$ in our simulations. Now, to use this initial data scheme, we need to transform the above into the quantities actually tracked by $\texttt{GiRaFFE}$ and HydroBase: $A_i$, $B^i$, $\tilde{S}_i$, $v^i$, and $\Phi$. Of these quantities, $\texttt{GiRaFFEfood}$ will only set $A_i$, $v^i$, and $\Phi=0$; $\texttt{GiRaFFE}$ itself will call functions to set $B^i$ and $\tilde{S}_i$ before the time-evolution begins. This can be done with eqs. 16 and 18, here given in that same order:
# \begin{align}
# v^i &= \alpha \frac{\epsilon^{ijk} E_j B_k}{B^2} -\beta^i \\
# B^i &= \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k \\
# \end{align}
# In the simulations, $B^i$ will be calculated numerically from $A_i$; however, it will be useful to analytically calculate $B^i$ to use calculating the initial $v^i$.
# 
# #### Table of Contents:
# 1. [Steps 0-1:](#preliminaries) Preliminaries
# 1. [Step 2:](#step2) Set the vectors A and E in Spherical coordinates
# 1. [Step 3:](#step3) Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
# 1. [Step 4:](#step4) Calculate $v^i$ from $A_i$ and $E_i$
# 1. [Step 5:](#step5) Build the expression for $\tilde{S}_i$
# 1. [Step 6:](#step6) NRPy+ Module Code Validation

# <a id='preliminaries'></a>
# 
# ### Steps 0-1: Preliminaries
# $$\label{preliminaries}$$
# 
# \[Back to [top](#top)\]
# 
# Here, we will import the NRPy+ core modules and set the reference metric to Cartesian, set commonly used NRPy+ parameters, and set C parameters that will be set from outside the code eventually generated from these expressions. We will also set up a parameter to determine what initial data is set up, although it won't do much yet.


# Step 0: Import the NRPy+ core modules and set the reference metric to Cartesian
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
from outputC import *
import loop

import reference_metric as rfm
par.set_parval_from_str("reference_metric::CoordSystem","Cartesian")
rfm.reference_metric()

# Step 1a: Set commonly used parameters.
thismodule = __name__

def GiRaFFEfood_NRPy_Exact_Wald(gammaDD,M,KerrSchild_radial_shift):

    # <a id='step2'></a>
    # 
    # ### Step 2: Set the vectors A and E in Spherical coordinates
    # $$\label{step2}$$
    # 
    # \[Back to [top](#top)\]
    # 
    # We will first build the fundamental vectors $A_i$ and $E_i$ in spherical coordinates (see [Table 3](https://arxiv.org/pdf/1704.00599.pdf)). Note that we use reference_metric.py to set $r$ and $\theta$ in terms of Cartesian coordinates; this will save us a step later when we convert to Cartesian coordinates. Since $C_0 = 1$,
    # \begin{align}
    # A_{\phi} &= \frac{1}{2} r^2 \sin^2 \theta \\
    # E_{\phi} &= 2 M \left( 1+ \frac {2M}{r} \right)^{-1/2} \sin^2 \theta. \\
    # \end{align}
    # While we have $E_i$ set as a variable in NRPy+, note that the final C code won't store these values.


    # Step 2: Set the vectors A and E in Spherical coordinates

    r     = rfm.xxSph[0] + KerrSchild_radial_shift # We are setting the data up in Shifted Kerr-Schild coordinates
    theta = rfm.xxSph[1]

    # Initialize all components of A and E in the *spherical basis* to zero
    ASphD = ixp.zerorank1()
    ESphD = ixp.zerorank1()
    ASphD[2] = (r * r * sp.sin(theta)**2)/2
    ESphD[2] = 2 * M * sp.sin(theta)**2 / sp.sqrt(1+2*M/r)


    # <a id='step3'></a>
    # 
    # ### Step 3: Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
    # $$\label{step3}$$
    # 
    # \[Back to [top](#top)\]
    # 
    # Now, we will use the coordinate transformation definitions provided by reference_metric.py to build the Jacobian 
    # $$ 
    # \frac{\partial x_{\rm Sph}^j}{\partial x_{\rm Cart}^i},
    # $$ 
    # where $x_{\rm Sph}^j \in \{r,\theta,\phi\}$ and $x_{\rm Cart}^i \in \{x,y,z\}$. We would normally compute its inverse, but since none of the quantities we need to transform have upper indices, it is not necessary. Then, since both $A_i$ and $E_i$ have one lower index, both will need to be multiplied by the Jacobian:
    # 
    # $$
    # A_i^{\rm Cart} = A_j^{\rm Sph} \frac{\partial x_{\rm Sph}^j}{\partial x_{\rm Cart}^i},
    # $$


    # Step 3: Use the Jacobian matrix to transform the vectors to Cartesian coordinates.
    drrefmetric__dx_0UDmatrix = sp.Matrix([[sp.diff(rfm.xxSph[0],rfm.xx[0]), sp.diff(rfm.xxSph[0],rfm.xx[1]), sp.diff(rfm.xxSph[0],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[1],rfm.xx[0]), sp.diff(rfm.xxSph[1],rfm.xx[1]), sp.diff(rfm.xxSph[1],rfm.xx[2])],
                                           [sp.diff(rfm.xxSph[2],rfm.xx[0]), sp.diff(rfm.xxSph[2],rfm.xx[1]), sp.diff(rfm.xxSph[2],rfm.xx[2])]])
    #dx__drrefmetric_0UDmatrix = drrefmetric__dx_0UDmatrix.inv() # We don't actually need this in this case.

    global AD
    AD = ixp.zerorank1(DIM=3)
    ED = ixp.zerorank1(DIM=3)

    for i in range(3):
        for j in range(3):
            AD[i] = drrefmetric__dx_0UDmatrix[(j,i)]*ASphD[j]
            ED[i] = drrefmetric__dx_0UDmatrix[(j,i)]*ESphD[j]

    import GRHD.equations as GRHD
    GRHD.compute_sqrtgammaDET(gammaDD)


    # <a id='step4'></a>
    # 
    # ### Step 4: Calculate $v^i$ from $A_i$ and $E_i$
    # $$\label{step4}$$
    # 
    # \[Back to [top](#top)\]
    # 
    # We will now find the magnetic field using equation 18 in [the original paper](https://arxiv.org/pdf/1704.00599.pdf) $$B^i = \frac{[ijk]}{\sqrt{\gamma}} \partial_j A_k. $$ We will need the metric quantites: the lapse $\alpha$, the shift $\beta^i$, and the three-metric $\gamma_{ij}$. We will also need the Levi-Civita symbol, provided by $\text{WeylScal4NRPy}$. 


    # Step 4: Calculate v^i from A_i and E_i
    # Step 4a: Calculate the magnetic field B^i
    # Here, we build the Levi-Civita tensor from the Levi-Civita symbol.
    # Here, we build the Levi-Civita tensor from the Levi-Civita symbol.
    import WeylScal4NRPy.WeylScalars_Cartesian as weyl
    LeviCivitaSymbolDDD = weyl.define_LeviCivitaSymbol_rank3()
    LeviCivitaTensorUUU = ixp.zerorank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LeviCivitaTensorUUU[i][j][k] = LeviCivitaSymbolDDD[i][j][k] / GRHD.sqrtgammaDET

    # For the initial data, we can analytically take the derivatives of A_i
    ADdD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            ADdD[i][j] = sp.simplify(sp.diff(AD[i],rfm.xxCart[j]))

    BU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                BU[i] += LeviCivitaTensorUUU[i][j][k] * ADdD[k][j]


    # We will now build the initial velocity using equation 152 in [this paper,](https://arxiv.org/pdf/1310.3274v2.pdf) cited in the original $\texttt{GiRaFFE}$ code: $$ v^i = \alpha \frac{\epsilon^{ijk} E_j B_k}{B^2} -\beta^i. $$ 
    # However, our code needs the Valencia 3-velocity while this expression is for the drift velocity. So, we will need to transform it to the Valencia 3-velocity using the rule $\bar{v}^i = \frac{1}{\alpha} \left(v^i +\beta^i \right)$.
    # Thus, $$\bar{v}^i = \frac{\epsilon^{ijk} E_j B_k}{B^2}$$


    # Step 4b: Calculate B^2 and B_i
    # B^2 is an inner product defined in the usual way:
    B2 = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            B2 += gammaDD[i][j] * BU[i] * BU[j]

    # Lower the index on B^i
    BD = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            BD[i] += gammaDD[i][j] * BU[j]

    # Step 4c: Calculate the Valencia 3-velocity 
    global ValenciavU
    ValenciavU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ValenciavU[i] += LeviCivitaTensorUUU[i][j][k]*ED[j]*BD[k]/B2
