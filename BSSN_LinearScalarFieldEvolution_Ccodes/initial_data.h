/*
 * Set up the initial data at all points on the numerical grid.
 */
void initial_data(const paramstruct *restrict params,REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "./set_Cparameters.h"

    #pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        const REAL xx2 = xx[2][i2];
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
            const REAL xx1 = xx[1][i1];
            for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
                const REAL xx0 = xx[0][i0];
                   {
                         const double tmp_0 = ((M)*(M));
                         const double tmp_2 = ((chi)*(chi))*tmp_0;
                         const double tmp_4 = sqrt(tmp_0 - tmp_2);
                         const double tmp_6 = (1.0/4.0)*M + (1.0/4.0)*tmp_4 + xx0;
                         const double tmp_7 = ((tmp_6)*(tmp_6));
                         const double tmp_8 = (1.0/(tmp_7));
                         const double tmp_9 = (1.0/(xx0));
                         const double tmp_11 = tmp_9*(M + tmp_4);
                         const double tmp_12 = (1.0/4.0)*tmp_11 + 1;
                         const double tmp_13 = ((tmp_12)*(tmp_12));
                         const double tmp_15 = -M + tmp_13*xx0 + tmp_4;
                         const double tmp_16 = cos(xx1);
                         const double tmp_17 = ((xx0)*(xx0));
                         const double tmp_18 = ((tmp_12)*(tmp_12)*(tmp_12)*(tmp_12));
                         const double tmp_19 = tmp_17*tmp_18;
                         const double tmp_20 = ((tmp_16)*(tmp_16))*tmp_2 + tmp_19;
                         const double tmp_21 = (1.0/(tmp_20));
                         const double tmp_22 = pow(xx0, 7);
                         const double tmp_23 = tmp_19 + tmp_2;
                         const double tmp_25 = -2*M*tmp_13*xx0 + tmp_23;
                         const double tmp_26 = sin(xx1);
                         const double tmp_27 = ((tmp_26)*(tmp_26));
                         const double tmp_28 = tmp_2*tmp_27;
                         const double tmp_29 = ((tmp_23)*(tmp_23)) - tmp_25*tmp_28;
                         const double tmp_30 = (1.0/(tmp_29));
                         const double tmp_31 = tmp_15*tmp_21*tmp_22*tmp_30;
                         const double tmp_32 = tmp_31*tmp_8;
                         const double tmp_33 = cbrt(tmp_32);
                         const double tmp_34 = ((M)*(M)*(M)*(M));
                         const double tmp_35 = ((xx0)*(xx0)*(xx0)*(xx0));
                         const double tmp_36 = 2*tmp_2;
                         const double tmp_38 = tmp_9/sqrt(tmp_20*tmp_29);
                         const double tmp_40 = tmp_21*tmp_27*tmp_33;
                         const double tmp_41 = (1.0/(tmp_15));
                         const double tmp_42 = tmp_41*tmp_7;
                         const double tmp_43 = tmp_42/tmp_22;
                         const double tmp_44 = pow(tmp_20*tmp_29*tmp_43, -1.0/6.0);
                         const double tmp_45 = tmp_20*tmp_33;
                         const double tmp_46 = (1.0/((xx0)*(xx0)*(xx0)));
                         const double tmp_47 = tmp_42*tmp_46;
                         const double tmp_48 = (1.0/(tmp_17));
                         const double tmp_49 = (1.0/(tmp_27));
                         const double tmp_50 = tmp_29*tmp_33;
                         const double tmp_52 = 2*tmp_18*xx0;
                         const double tmp_53 = ((tmp_12)*(tmp_12)*(tmp_12))*(M + tmp_4);
                         const double tmp_54 = tmp_52 - tmp_53;
                         const double tmp_56 = tmp_21*tmp_30*tmp_8;
                         const double tmp_58 = (1.0/2.0)*tmp_11*tmp_12;
                         const double tmp_60 = ((tmp_20)*(tmp_20));
                         const double tmp_61 = (1.0/(tmp_60));
                         const double tmp_62 = tmp_61*(-tmp_52 + tmp_53);
                         const double tmp_64 = (1.0/3.0)*tmp_15*tmp_22*tmp_8;
                         const double tmp_65 = ((tmp_29)*(tmp_29));
                         const double tmp_66 = (1.0/(tmp_65));
                         const double tmp_67 = tmp_23*(4*tmp_18*xx0 - 2*tmp_53);
                         const double tmp_68 = tmp_28*(M*tmp_11*tmp_12 - 2*M*tmp_13 + tmp_54);
                         const double tmp_69 = (7.0/3.0)*tmp_15*tmp_56*pow(xx0, 6) + tmp_21*tmp_64*tmp_66*(-tmp_67 + tmp_68) + (1.0/3.0)*tmp_22*tmp_56*(tmp_13 - tmp_58) + tmp_30*tmp_62*tmp_64 - 2.0/3.0*tmp_31/((tmp_6)*(tmp_6)*(tmp_6));
                         const double tmp_71 = tmp_50*tmp_60*tmp_69;
                         const double tmp_72 = (1.0/(tmp_35));
                         const double tmp_73 = pow(tmp_32, 2.0/3.0);
                         const double tmp_75 = (1.0/2.0)*tmp_29*tmp_73;
                         const double tmp_77 = tmp_27*tmp_33*tmp_43*tmp_65;
                         const double tmp_78 = tmp_43*tmp_49*tmp_60*tmp_73;
                         const double tmp_80 = (1.0/((tmp_15)*(tmp_15)));
                         const double tmp_81 = ((tmp_6)*(tmp_6)*(tmp_6)*(tmp_6))*tmp_80;
                         const double tmp_82 = tmp_81/pow(xx0, 10);
                         const double tmp_83 = (1.0/2.0)*pow(tmp_32, 4.0/3.0)*tmp_65;
                         const double tmp_84 = tmp_16*tmp_26;
                         const double tmp_86 = tmp_33*tmp_36*tmp_84;
                         const double tmp_87 = (2.0/3.0)*tmp_15*tmp_2*tmp_22*tmp_8;
                         const double tmp_89 = tmp_21*tmp_25*tmp_66*tmp_84*tmp_87 + tmp_30*tmp_61*tmp_84*tmp_87;
                         const double tmp_90 = tmp_50*tmp_60*tmp_89;
                         const double tmp_91 = tmp_16*((tmp_26)*(tmp_26)*(tmp_26));
                         in_gfs[IDX4S(ADD00GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(ADD01GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(ADD02GF,i0,i1,i2)] = chi*tmp_0*tmp_12*tmp_21*tmp_26*tmp_33*tmp_38*(-((chi)*(chi)*(chi)*(chi))*tmp_34 + 3*pow(tmp_12, 8)*tmp_35 + tmp_19*tmp_36 - tmp_28*(tmp_19 - tmp_2))/sqrt(tmp_15*xx0);
                         in_gfs[IDX4S(ADD11GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(ADD12GF,i0,i1,i2)] = -2*((chi)*(chi)*(chi))*tmp_13*tmp_16*tmp_34*tmp_38*tmp_40*sqrt(tmp_15*tmp_9)*(-1.0/4.0*M - 1.0/4.0*tmp_4 + xx0);
                         in_gfs[IDX4S(ADD22GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(ALPHAGF,i0,i1,i2)] = tmp_44;
                         in_gfs[IDX4S(BETU0GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(BETU1GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(BETU2GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(CFGF,i0,i1,i2)] = tmp_44;
                         in_gfs[IDX4S(HDD00GF,i0,i1,i2)] = tmp_45*tmp_47 - 1;
                         in_gfs[IDX4S(HDD01GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(HDD02GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(HDD11GF,i0,i1,i2)] = tmp_48*(-tmp_17 + tmp_45);
                         in_gfs[IDX4S(HDD12GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(HDD22GF,i0,i1,i2)] = tmp_48*tmp_49*(-tmp_17*tmp_27 + tmp_21*tmp_27*tmp_50);
                         in_gfs[IDX4S(LAMBDAU0GF,i0,i1,i2)] = tmp_29*tmp_43*tmp_73*(tmp_72*tmp_75*(-tmp_33*tmp_54 - tmp_43*tmp_71) + xx0) + tmp_78*(tmp_27*xx0 + tmp_72*tmp_75*(-tmp_27*tmp_50*tmp_62 - tmp_40*(tmp_67 - tmp_68) - tmp_69*tmp_77)) + tmp_83*(tmp_33*tmp_47*tmp_54 + tmp_41*tmp_45*tmp_46*((1.0/2.0)*M + (1.0/2.0)*tmp_4 + 2*xx0) - 3*tmp_42*tmp_45*tmp_72 + tmp_45*tmp_46*tmp_7*tmp_80*(-tmp_13 + tmp_58) + tmp_71*tmp_82)/pow(xx0, 8);
                         in_gfs[IDX4S(LAMBDAU1GF,i0,i1,i2)] = xx0*(tmp_42*tmp_83*(tmp_47*tmp_86 - tmp_82*tmp_90)/pow(xx0, 11) + tmp_78*(tmp_43*tmp_75*(tmp_21*tmp_25*tmp_33*tmp_36*tmp_91 - 2*tmp_21*tmp_50*tmp_84 - tmp_36*tmp_50*tmp_61*tmp_91 - tmp_77*tmp_89) + (1.0/2.0)*sin(2*xx1)) + tmp_81*tmp_83*(tmp_43*tmp_90 - tmp_86)/pow(xx0, 14));
                         in_gfs[IDX4S(LAMBDAU2GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(TRKGF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(VETU0GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(VETU1GF,i0,i1,i2)] = 0;
                         in_gfs[IDX4S(VETU2GF,i0,i1,i2)] = 0;
                   }
                
            } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
        } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
}
