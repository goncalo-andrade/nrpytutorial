/*
 * Set up the scalar field initial data at all points on the numerical grid.
 */
void fields_initial_data(const paramstruct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "./set_Cparameters.h"

    #pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        const REAL xx2 = xx[2][i2];
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            const REAL xx1 = xx[1][i1];
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                const REAL xx0 = xx[0][i0];
                   {
                         const double tmp_0 = sin(xx1);
                         const double tmp_1 = sqrt(6);
                         const double tmp_2 = (1.0/(SINHW));
                         const double tmp_3 = exp(tmp_2) - exp(-tmp_2);
                         const double tmp_5 = exp(tmp_2*xx0) - exp(-tmp_2*xx0);
                         const double tmp_6 = AMPL*tmp_5/tmp_3;
                         const double tmp_7 = ((w)*(w));
                         const double tmp_8 = (1.0/(tmp_7));
                         const double tmp_9 = -r0 + tmp_6;
                         const double tmp_10 = tmp_8*((tmp_9)*(tmp_9));
                         const double tmp_12 = ((r0)*(r0));
                         const double tmp_14 = exp(-2*tmp_10);
                         const double tmp_15 = (1.0/sqrt(M_PI));
                         const double tmp_17 = M_SQRT2;
                         const double tmp_18 = tmp_17/w;
                         const double tmp_19 = erf(tmp_18*tmp_9) - 1;
                         const double tmp_22 = tmp_3/(AMPL*tmp_5);
                         const double tmp_24 = ((M)*(M))*((chi)*(chi));
                         const double tmp_25 = sqrt(((M)*(M)) - tmp_24);
                         const double tmp_26 = ((cos(xx1))*(cos(xx1)));
                         const double tmp_28 = (1.0/4.0)*tmp_22*(M + tmp_25) + 1;
                         const double tmp_32 = ((AMPL)*(AMPL))*((tmp_5)*(tmp_5))/((tmp_3)*(tmp_3));
                         const double tmp_33 = ((tmp_28)*(tmp_28)*(tmp_28)*(tmp_28))*tmp_32;
                         const double tmp_34 = ((tmp_28)*(tmp_28))*tmp_6;
                         const double tmp_37 = ((r0)*(r0)*(r0)*(r0));
                         const double tmp_38 = ((w)*(w)*(w)*(w));
                         const double tmp_41 = ((A)*(A))*((tmp_3)*(tmp_3))/(((AMPL)*(AMPL))*((tmp_5)*(tmp_5)));
                         const double tmp_42 = (1.0/1600.0)*sqrt(15)*tmp_41*w;
                         const double tmp_44 = 4*tmp_37 + 2*tmp_38;
                         const double tmp_45 = sqrt(30)*tmp_15;
                         const double tmp_46 = (1.0/800.0)*tmp_41*tmp_45*tmp_7;
                         const double tmp_47 = ((r0)*(r0)*(r0));
                         const double tmp_49 = r0*tmp_42*(erf(r0*tmp_18) + 1)*(40*tmp_12*tmp_7 + 16*tmp_37 + 15*tmp_38) - tmp_14*tmp_46*(4*((AMPL)*(AMPL)*(AMPL)*(AMPL))*((tmp_5)*(tmp_5)*(tmp_5)*(tmp_5))/((tmp_3)*(tmp_3)*(tmp_3)*(tmp_3)) + 4*((AMPL)*(AMPL)*(AMPL))*r0*((tmp_5)*(tmp_5)*(tmp_5))/((tmp_3)*(tmp_3)*(tmp_3)) + 4*tmp_12*tmp_32 + tmp_44 + 4*tmp_47*tmp_6 + tmp_7*(7*r0*tmp_6 + 9*tmp_12 + 4*tmp_32)) + tmp_19*tmp_42*(-16*((AMPL)*(AMPL)*(AMPL)*(AMPL)*(AMPL))*((tmp_5)*(tmp_5)*(tmp_5)*(tmp_5)*(tmp_5))/((tmp_3)*(tmp_3)*(tmp_3)*(tmp_3)*(tmp_3)) + 16*((r0)*(r0)*(r0)*(r0)*(r0)) + 15*r0*tmp_38 + 40*tmp_47*tmp_7) + tmp_46*(9*tmp_12*tmp_7 + tmp_44)*exp(-2*tmp_12*tmp_8);
                         const double tmp_50 = sqrt(5)*tmp_15;
                         in_gfs[IDX4S(PHIGF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(PIGF,i0,i1,i2)] = (1.0/2.0)*A*tmp_0*tmp_1*sqrt(tmp_6)*exp(-tmp_10)*cos(xx2)/(M_PI*pow((1.0/32.0)*((A)*(A))*tmp_15*tmp_22*w*(-2*tmp_14*tmp_15*w*(2*tmp_12 + tmp_7) - tmp_17*tmp_19*(tmp_12*(-4*r0 + 4*tmp_6) + tmp_7*(-3*r0 + tmp_6))) + (1.0/4.0)*((tmp_0)*(tmp_0))*tmp_22*tmp_45*tmp_49*cos(2*xx2) - 1.0/3.0*tmp_1*tmp_22*tmp_49*((3.0/4.0)*tmp_26*tmp_50 - 1.0/4.0*tmp_50) + pow(pow(tmp_3, 7)*(tmp_24*tmp_26 + tmp_33)*(-((tmp_0)*(tmp_0))*tmp_24*(-2*M*tmp_34 + tmp_24 + tmp_33) + ((tmp_24 + tmp_33)*(tmp_24 + tmp_33)))*(((1.0/4.0)*M + (1.0/4.0)*tmp_25 + tmp_6)*((1.0/4.0)*M + (1.0/4.0)*tmp_25 + tmp_6))/(pow(AMPL, 7)*pow(tmp_5, 7)*(-M + tmp_25 + tmp_34)), 1.0/12.0), 5.0/2.0));
                   }
                
            } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
