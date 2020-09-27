
#ifndef __FD_FUNCTIONS_H__
#define __FD_FUNCTIONS_H__
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#define _UNUSED   __attribute__((unused))
#define _NOINLINE __attribute__((noinline))
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_3_20 = 3.0/20.0;
      const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

      const double tmp_Rational_3_4 = 3.0/4.0;
      const REAL_SIMD_ARRAY _Rational_3_4 = ConstSIMD(tmp_Rational_3_4);

      return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_3_20, SubSIMD(f_i0m2_i1_i2, f_i0p2_i1_i2), FusedMulAddSIMD(_Rational_3_4, SubSIMD(f_i0p1_i1_i2, f_i0m1_i1_i2), MulSIMD(_Rational_1_60, SubSIMD(f_i0p3_i1_i2, f_i0m3_i1_i2)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_3_20 = 3.0/20.0;
      const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

      const double tmp_Rational_3_4 = 3.0/4.0;
      const REAL_SIMD_ARRAY _Rational_3_4 = ConstSIMD(tmp_Rational_3_4);

      return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_3_20, SubSIMD(f_i0_i1m2_i2, f_i0_i1p2_i2), FusedMulAddSIMD(_Rational_3_4, SubSIMD(f_i0_i1p1_i2, f_i0_i1m1_i2), MulSIMD(_Rational_1_60, SubSIMD(f_i0_i1p3_i2, f_i0_i1m3_i2)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_3_20 = 3.0/20.0;
      const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

      const double tmp_Rational_3_4 = 3.0/4.0;
      const REAL_SIMD_ARRAY _Rational_3_4 = ConstSIMD(tmp_Rational_3_4);

      return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_3_20, SubSIMD(f_i0_i1_i2m2, f_i0_i1_i2p2), FusedMulAddSIMD(_Rational_3_4, SubSIMD(f_i0_i1_i2p1, f_i0_i1_i2m1), MulSIMD(_Rational_1_60, SubSIMD(f_i0_i1_i2p3, f_i0_i1_i2m3)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 0
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_ddnD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_30 = 1.0/30.0;
      const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_2_15 = 2.0/15.0;
      const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

      const double tmp_Rational_2_5 = 2.0/5.0;
      const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_7_12 = 7.0/12.0;
      const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

      return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_7_12, f, SubSIMD(FusedMulAddSIMD(_Rational_1_60, f_i0m4_i1_i2, FusedMulAddSIMD(_Rational_2_5, f_i0p1_i1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0m2_i1_i2, MulSIMD(_Rational_1_30, f_i0p2_i1_i2)), MulSIMD(_Rational_4_3, f_i0m1_i1_i2)))), MulSIMD(_Rational_2_15, f_i0m3_i1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_ddnD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_30 = 1.0/30.0;
      const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_2_15 = 2.0/15.0;
      const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

      const double tmp_Rational_2_5 = 2.0/5.0;
      const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_7_12 = 7.0/12.0;
      const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

      return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_7_12, f, SubSIMD(FusedMulAddSIMD(_Rational_1_60, f_i0_i1m4_i2, FusedMulAddSIMD(_Rational_2_5, f_i0_i1p1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0_i1m2_i2, MulSIMD(_Rational_1_30, f_i0_i1p2_i2)), MulSIMD(_Rational_4_3, f_i0_i1m1_i2)))), MulSIMD(_Rational_2_15, f_i0_i1m3_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_ddnD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_30 = 1.0/30.0;
      const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_2_15 = 2.0/15.0;
      const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

      const double tmp_Rational_2_5 = 2.0/5.0;
      const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_7_12 = 7.0/12.0;
      const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

      return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_7_12, f, SubSIMD(FusedMulAddSIMD(_Rational_1_60, f_i0_i1_i2m4, FusedMulAddSIMD(_Rational_2_5, f_i0_i1_i2p1, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0_i1_i2m2, MulSIMD(_Rational_1_30, f_i0_i1_i2p2)), MulSIMD(_Rational_4_3, f_i0_i1_i2m1)))), MulSIMD(_Rational_2_15, f_i0_i1_i2m3))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 0
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dupD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_30 = 1.0/30.0;
      const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_2_15 = 2.0/15.0;
      const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

      const double tmp_Rational_2_5 = 2.0/5.0;
      const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_7_12 = 7.0/12.0;
      const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

      return MulSIMD(invdx0, SubSIMD(SubSIMD(FusedMulAddSIMD(_Rational_2_15, f_i0p3_i1_i2, FusedMulAddSIMD(_Rational_4_3, f_i0p1_i1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_30, f_i0m2_i1_i2, MulSIMD(_Rational_1_2, f_i0p2_i1_i2)), MulSIMD(_Rational_7_12, f)))), MulSIMD(_Rational_2_5, f_i0m1_i1_i2)), MulSIMD(_Rational_1_60, f_i0p4_i1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dupD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_30 = 1.0/30.0;
      const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_2_15 = 2.0/15.0;
      const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

      const double tmp_Rational_2_5 = 2.0/5.0;
      const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_7_12 = 7.0/12.0;
      const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

      return MulSIMD(invdx1, SubSIMD(SubSIMD(FusedMulAddSIMD(_Rational_2_15, f_i0_i1p3_i2, FusedMulAddSIMD(_Rational_4_3, f_i0_i1p1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_30, f_i0_i1m2_i2, MulSIMD(_Rational_1_2, f_i0_i1p2_i2)), MulSIMD(_Rational_7_12, f)))), MulSIMD(_Rational_2_5, f_i0_i1m1_i2)), MulSIMD(_Rational_1_60, f_i0_i1p4_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dupD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_30 = 1.0/30.0;
      const REAL_SIMD_ARRAY _Rational_1_30 = ConstSIMD(tmp_Rational_1_30);

      const double tmp_Rational_1_60 = 1.0/60.0;
      const REAL_SIMD_ARRAY _Rational_1_60 = ConstSIMD(tmp_Rational_1_60);

      const double tmp_Rational_2_15 = 2.0/15.0;
      const REAL_SIMD_ARRAY _Rational_2_15 = ConstSIMD(tmp_Rational_2_15);

      const double tmp_Rational_2_5 = 2.0/5.0;
      const REAL_SIMD_ARRAY _Rational_2_5 = ConstSIMD(tmp_Rational_2_5);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_7_12 = 7.0/12.0;
      const REAL_SIMD_ARRAY _Rational_7_12 = ConstSIMD(tmp_Rational_7_12);

      return MulSIMD(invdx2, SubSIMD(SubSIMD(FusedMulAddSIMD(_Rational_2_15, f_i0_i1_i2p3, FusedMulAddSIMD(_Rational_4_3, f_i0_i1_i2p1, SubSIMD(FusedMulSubSIMD(_Rational_1_30, f_i0_i1_i2m2, MulSIMD(_Rational_1_2, f_i0_i1_i2p2)), MulSIMD(_Rational_7_12, f)))), MulSIMD(_Rational_2_5, f_i0_i1_i2m1)), MulSIMD(_Rational_1_60, f_i0_i1_i2p4)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD00(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_90 = 1.0/90.0;
      const REAL_SIMD_ARRAY _Rational_1_90 = ConstSIMD(tmp_Rational_1_90);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_3_20 = 3.0/20.0;
      const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

      const double tmp_Rational_49_18 = 49.0/18.0;
      const REAL_SIMD_ARRAY _Rational_49_18 = ConstSIMD(tmp_Rational_49_18);

      return MulSIMD(MulSIMD(invdx0, invdx0), FusedMulAddSIMD(_Rational_3_2, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), FusedMulAddSIMD(_Rational_3_20, MulSIMD(_NegativeOne_, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2)), FusedMulSubSIMD(_Rational_1_90, AddSIMD(f_i0m3_i1_i2, f_i0p3_i1_i2), MulSIMD(_Rational_49_18, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD01(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0m3_i1m3_i2,const REAL_SIMD_ARRAY f_i0m2_i1m3_i2,const REAL_SIMD_ARRAY f_i0m1_i1m3_i2,const REAL_SIMD_ARRAY f_i0p1_i1m3_i2,const REAL_SIMD_ARRAY f_i0p2_i1m3_i2,const REAL_SIMD_ARRAY f_i0p3_i1m3_i2,const REAL_SIMD_ARRAY f_i0m3_i1m2_i2,const REAL_SIMD_ARRAY f_i0m2_i1m2_i2,const REAL_SIMD_ARRAY f_i0m1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p2_i1m2_i2,const REAL_SIMD_ARRAY f_i0p3_i1m2_i2,const REAL_SIMD_ARRAY f_i0m3_i1m1_i2,const REAL_SIMD_ARRAY f_i0m2_i1m1_i2,const REAL_SIMD_ARRAY f_i0m1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p2_i1m1_i2,const REAL_SIMD_ARRAY f_i0p3_i1m1_i2,const REAL_SIMD_ARRAY f_i0m3_i1p1_i2,const REAL_SIMD_ARRAY f_i0m2_i1p1_i2,const REAL_SIMD_ARRAY f_i0m1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p2_i1p1_i2,const REAL_SIMD_ARRAY f_i0p3_i1p1_i2,const REAL_SIMD_ARRAY f_i0m3_i1p2_i2,const REAL_SIMD_ARRAY f_i0m2_i1p2_i2,const REAL_SIMD_ARRAY f_i0m1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p2_i1p2_i2,const REAL_SIMD_ARRAY f_i0p3_i1p2_i2,const REAL_SIMD_ARRAY f_i0m3_i1p3_i2,const REAL_SIMD_ARRAY f_i0m2_i1p3_i2,const REAL_SIMD_ARRAY f_i0m1_i1p3_i2,const REAL_SIMD_ARRAY f_i0p1_i1p3_i2,const REAL_SIMD_ARRAY f_i0p2_i1p3_i2,const REAL_SIMD_ARRAY f_i0p3_i1p3_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_3600 = 1.0/3600.0;
      const REAL_SIMD_ARRAY _Rational_1_3600 = ConstSIMD(tmp_Rational_1_3600);

      const double tmp_Rational_1_400 = 1.0/400.0;
      const REAL_SIMD_ARRAY _Rational_1_400 = ConstSIMD(tmp_Rational_1_400);

      const double tmp_Rational_1_80 = 1.0/80.0;
      const REAL_SIMD_ARRAY _Rational_1_80 = ConstSIMD(tmp_Rational_1_80);

      const double tmp_Rational_9_16 = 9.0/16.0;
      const REAL_SIMD_ARRAY _Rational_9_16 = ConstSIMD(tmp_Rational_9_16);

      const double tmp_Rational_9_400 = 9.0/400.0;
      const REAL_SIMD_ARRAY _Rational_9_400 = ConstSIMD(tmp_Rational_9_400);

      const double tmp_Rational_9_80 = 9.0/80.0;
      const REAL_SIMD_ARRAY _Rational_9_80 = ConstSIMD(tmp_Rational_9_80);

      return MulSIMD(invdx0, MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_80, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1p1_i2, f_i0m3_i1p1_i2), AddSIMD(AddSIMD(f_i0m3_i1m1_i2, f_i0p1_i1p3_i2), SubSIMD(SubSIMD(f_i0m1_i1m3_i2, f_i0m1_i1p3_i2), f_i0p3_i1m1_i2))), f_i0p1_i1m3_i2), FusedMulAddSIMD(_Rational_9_16, AddSIMD(f_i0p1_i1p1_i2, SubSIMD(SubSIMD(f_i0m1_i1m1_i2, f_i0m1_i1p1_i2), f_i0p1_i1m1_i2)), FusedMulAddSIMD(_Rational_9_400, AddSIMD(f_i0p2_i1p2_i2, SubSIMD(SubSIMD(f_i0m2_i1m2_i2, f_i0m2_i1p2_i2), f_i0p2_i1m2_i2)), FusedMulAddSIMD(_Rational_9_80, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1m1_i2, f_i0m2_i1m1_i2), AddSIMD(AddSIMD(f_i0m2_i1p1_i2, f_i0p1_i1m2_i2), SubSIMD(SubSIMD(f_i0m1_i1p2_i2, f_i0m1_i1m2_i2), f_i0p2_i1p1_i2))), f_i0p1_i1p2_i2), FusedMulAddSIMD(_Rational_1_3600, AddSIMD(f_i0p3_i1p3_i2, SubSIMD(SubSIMD(f_i0m3_i1m3_i2, f_i0m3_i1p3_i2), f_i0p3_i1m3_i2)), MulSIMD(_Rational_1_400, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1m2_i2, f_i0m3_i1m2_i2), AddSIMD(AddSIMD(f_i0m3_i1p2_i2, f_i0p2_i1m3_i2), SubSIMD(SubSIMD(f_i0m2_i1p3_i2, f_i0m2_i1m3_i2), f_i0p3_i1p2_i2))), f_i0p2_i1p3_i2)))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD02(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0m3_i1_i2m3,const REAL_SIMD_ARRAY f_i0m2_i1_i2m3,const REAL_SIMD_ARRAY f_i0m1_i1_i2m3,const REAL_SIMD_ARRAY f_i0p1_i1_i2m3,const REAL_SIMD_ARRAY f_i0p2_i1_i2m3,const REAL_SIMD_ARRAY f_i0p3_i1_i2m3,const REAL_SIMD_ARRAY f_i0m3_i1_i2m2,const REAL_SIMD_ARRAY f_i0m2_i1_i2m2,const REAL_SIMD_ARRAY f_i0m1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p2_i1_i2m2,const REAL_SIMD_ARRAY f_i0p3_i1_i2m2,const REAL_SIMD_ARRAY f_i0m3_i1_i2m1,const REAL_SIMD_ARRAY f_i0m2_i1_i2m1,const REAL_SIMD_ARRAY f_i0m1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p2_i1_i2m1,const REAL_SIMD_ARRAY f_i0p3_i1_i2m1,const REAL_SIMD_ARRAY f_i0m3_i1_i2p1,const REAL_SIMD_ARRAY f_i0m2_i1_i2p1,const REAL_SIMD_ARRAY f_i0m1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p2_i1_i2p1,const REAL_SIMD_ARRAY f_i0p3_i1_i2p1,const REAL_SIMD_ARRAY f_i0m3_i1_i2p2,const REAL_SIMD_ARRAY f_i0m2_i1_i2p2,const REAL_SIMD_ARRAY f_i0m1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p2_i1_i2p2,const REAL_SIMD_ARRAY f_i0p3_i1_i2p2,const REAL_SIMD_ARRAY f_i0m3_i1_i2p3,const REAL_SIMD_ARRAY f_i0m2_i1_i2p3,const REAL_SIMD_ARRAY f_i0m1_i1_i2p3,const REAL_SIMD_ARRAY f_i0p1_i1_i2p3,const REAL_SIMD_ARRAY f_i0p2_i1_i2p3,const REAL_SIMD_ARRAY f_i0p3_i1_i2p3) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_3600 = 1.0/3600.0;
      const REAL_SIMD_ARRAY _Rational_1_3600 = ConstSIMD(tmp_Rational_1_3600);

      const double tmp_Rational_1_400 = 1.0/400.0;
      const REAL_SIMD_ARRAY _Rational_1_400 = ConstSIMD(tmp_Rational_1_400);

      const double tmp_Rational_1_80 = 1.0/80.0;
      const REAL_SIMD_ARRAY _Rational_1_80 = ConstSIMD(tmp_Rational_1_80);

      const double tmp_Rational_9_16 = 9.0/16.0;
      const REAL_SIMD_ARRAY _Rational_9_16 = ConstSIMD(tmp_Rational_9_16);

      const double tmp_Rational_9_400 = 9.0/400.0;
      const REAL_SIMD_ARRAY _Rational_9_400 = ConstSIMD(tmp_Rational_9_400);

      const double tmp_Rational_9_80 = 9.0/80.0;
      const REAL_SIMD_ARRAY _Rational_9_80 = ConstSIMD(tmp_Rational_9_80);

      return MulSIMD(invdx0, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_80, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1_i2p1, f_i0m3_i1_i2p1), AddSIMD(AddSIMD(f_i0m3_i1_i2m1, f_i0p1_i1_i2p3), SubSIMD(SubSIMD(f_i0m1_i1_i2m3, f_i0m1_i1_i2p3), f_i0p3_i1_i2m1))), f_i0p1_i1_i2m3), FusedMulAddSIMD(_Rational_9_16, AddSIMD(f_i0p1_i1_i2p1, SubSIMD(SubSIMD(f_i0m1_i1_i2m1, f_i0m1_i1_i2p1), f_i0p1_i1_i2m1)), FusedMulAddSIMD(_Rational_9_400, AddSIMD(f_i0p2_i1_i2p2, SubSIMD(SubSIMD(f_i0m2_i1_i2m2, f_i0m2_i1_i2p2), f_i0p2_i1_i2m2)), FusedMulAddSIMD(_Rational_9_80, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1_i2m1, f_i0m2_i1_i2m1), AddSIMD(AddSIMD(f_i0m2_i1_i2p1, f_i0p1_i1_i2m2), SubSIMD(SubSIMD(f_i0m1_i1_i2p2, f_i0m1_i1_i2m2), f_i0p2_i1_i2p1))), f_i0p1_i1_i2p2), FusedMulAddSIMD(_Rational_1_3600, AddSIMD(f_i0p3_i1_i2p3, SubSIMD(SubSIMD(f_i0m3_i1_i2m3, f_i0m3_i1_i2p3), f_i0p3_i1_i2m3)), MulSIMD(_Rational_1_400, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1_i2m2, f_i0m3_i1_i2m2), AddSIMD(AddSIMD(f_i0m3_i1_i2p2, f_i0p2_i1_i2m3), SubSIMD(SubSIMD(f_i0m2_i1_i2p3, f_i0m2_i1_i2m3), f_i0p3_i1_i2p2))), f_i0p2_i1_i2p3)))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD11(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_90 = 1.0/90.0;
      const REAL_SIMD_ARRAY _Rational_1_90 = ConstSIMD(tmp_Rational_1_90);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_3_20 = 3.0/20.0;
      const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

      const double tmp_Rational_49_18 = 49.0/18.0;
      const REAL_SIMD_ARRAY _Rational_49_18 = ConstSIMD(tmp_Rational_49_18);

      return MulSIMD(MulSIMD(invdx1, invdx1), FusedMulAddSIMD(_Rational_3_2, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), FusedMulAddSIMD(_Rational_3_20, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2)), FusedMulSubSIMD(_Rational_1_90, AddSIMD(f_i0_i1m3_i2, f_i0_i1p3_i2), MulSIMD(_Rational_49_18, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD12(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1m3_i2m3,const REAL_SIMD_ARRAY f_i0_i1m2_i2m3,const REAL_SIMD_ARRAY f_i0_i1m1_i2m3,const REAL_SIMD_ARRAY f_i0_i1p1_i2m3,const REAL_SIMD_ARRAY f_i0_i1p2_i2m3,const REAL_SIMD_ARRAY f_i0_i1p3_i2m3,const REAL_SIMD_ARRAY f_i0_i1m3_i2m2,const REAL_SIMD_ARRAY f_i0_i1m2_i2m2,const REAL_SIMD_ARRAY f_i0_i1m1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p2_i2m2,const REAL_SIMD_ARRAY f_i0_i1p3_i2m2,const REAL_SIMD_ARRAY f_i0_i1m3_i2m1,const REAL_SIMD_ARRAY f_i0_i1m2_i2m1,const REAL_SIMD_ARRAY f_i0_i1m1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p2_i2m1,const REAL_SIMD_ARRAY f_i0_i1p3_i2m1,const REAL_SIMD_ARRAY f_i0_i1m3_i2p1,const REAL_SIMD_ARRAY f_i0_i1m2_i2p1,const REAL_SIMD_ARRAY f_i0_i1m1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p2_i2p1,const REAL_SIMD_ARRAY f_i0_i1p3_i2p1,const REAL_SIMD_ARRAY f_i0_i1m3_i2p2,const REAL_SIMD_ARRAY f_i0_i1m2_i2p2,const REAL_SIMD_ARRAY f_i0_i1m1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p2_i2p2,const REAL_SIMD_ARRAY f_i0_i1p3_i2p2,const REAL_SIMD_ARRAY f_i0_i1m3_i2p3,const REAL_SIMD_ARRAY f_i0_i1m2_i2p3,const REAL_SIMD_ARRAY f_i0_i1m1_i2p3,const REAL_SIMD_ARRAY f_i0_i1p1_i2p3,const REAL_SIMD_ARRAY f_i0_i1p2_i2p3,const REAL_SIMD_ARRAY f_i0_i1p3_i2p3) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_3600 = 1.0/3600.0;
      const REAL_SIMD_ARRAY _Rational_1_3600 = ConstSIMD(tmp_Rational_1_3600);

      const double tmp_Rational_1_400 = 1.0/400.0;
      const REAL_SIMD_ARRAY _Rational_1_400 = ConstSIMD(tmp_Rational_1_400);

      const double tmp_Rational_1_80 = 1.0/80.0;
      const REAL_SIMD_ARRAY _Rational_1_80 = ConstSIMD(tmp_Rational_1_80);

      const double tmp_Rational_9_16 = 9.0/16.0;
      const REAL_SIMD_ARRAY _Rational_9_16 = ConstSIMD(tmp_Rational_9_16);

      const double tmp_Rational_9_400 = 9.0/400.0;
      const REAL_SIMD_ARRAY _Rational_9_400 = ConstSIMD(tmp_Rational_9_400);

      const double tmp_Rational_9_80 = 9.0/80.0;
      const REAL_SIMD_ARRAY _Rational_9_80 = ConstSIMD(tmp_Rational_9_80);

      return MulSIMD(invdx1, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_80, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p3_i2p1, f_i0_i1m3_i2p1), AddSIMD(AddSIMD(f_i0_i1m3_i2m1, f_i0_i1p1_i2p3), SubSIMD(SubSIMD(f_i0_i1m1_i2m3, f_i0_i1m1_i2p3), f_i0_i1p3_i2m1))), f_i0_i1p1_i2m3), FusedMulAddSIMD(_Rational_9_16, AddSIMD(f_i0_i1p1_i2p1, SubSIMD(SubSIMD(f_i0_i1m1_i2m1, f_i0_i1m1_i2p1), f_i0_i1p1_i2m1)), FusedMulAddSIMD(_Rational_9_400, AddSIMD(f_i0_i1p2_i2p2, SubSIMD(SubSIMD(f_i0_i1m2_i2m2, f_i0_i1m2_i2p2), f_i0_i1p2_i2m2)), FusedMulAddSIMD(_Rational_9_80, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p2_i2m1, f_i0_i1m2_i2m1), AddSIMD(AddSIMD(f_i0_i1m2_i2p1, f_i0_i1p1_i2m2), SubSIMD(SubSIMD(f_i0_i1m1_i2p2, f_i0_i1m1_i2m2), f_i0_i1p2_i2p1))), f_i0_i1p1_i2p2), FusedMulAddSIMD(_Rational_1_3600, AddSIMD(f_i0_i1p3_i2p3, SubSIMD(SubSIMD(f_i0_i1m3_i2m3, f_i0_i1m3_i2p3), f_i0_i1p3_i2m3)), MulSIMD(_Rational_1_400, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p3_i2m2, f_i0_i1m3_i2m2), AddSIMD(AddSIMD(f_i0_i1m3_i2p2, f_i0_i1p2_i2m3), SubSIMD(SubSIMD(f_i0_i1m2_i2p3, f_i0_i1m2_i2m3), f_i0_i1p3_i2p2))), f_i0_i1p2_i2p3)))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_6_f_dDD22(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_90 = 1.0/90.0;
      const REAL_SIMD_ARRAY _Rational_1_90 = ConstSIMD(tmp_Rational_1_90);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_3_20 = 3.0/20.0;
      const REAL_SIMD_ARRAY _Rational_3_20 = ConstSIMD(tmp_Rational_3_20);

      const double tmp_Rational_49_18 = 49.0/18.0;
      const REAL_SIMD_ARRAY _Rational_49_18 = ConstSIMD(tmp_Rational_49_18);

      return MulSIMD(MulSIMD(invdx2, invdx2), FusedMulAddSIMD(_Rational_3_2, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), FusedMulAddSIMD(_Rational_3_20, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2)), FusedMulSubSIMD(_Rational_1_90, AddSIMD(f_i0_i1_i2m3, f_i0_i1_i2p3), MulSIMD(_Rational_49_18, f)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0
 */
static REAL _NOINLINE _UNUSED order_6_f_dD0(const REAL invdx0,const REAL f_i0m3_i1_i2,const REAL f_i0m2_i1_i2,const REAL f_i0m1_i1_i2,const REAL f_i0p1_i1_i2,const REAL f_i0p2_i1_i2,const REAL f_i0p3_i1_i2) {

      const double _Rational_3_4 = 3.0/4.0;
      const double _Rational_3_20 = 3.0/20.0;
      const double _Rational_1_60 = 1.0/60.0;
      return invdx0*(_Rational_1_60*(-f_i0m3_i1_i2 + f_i0p3_i1_i2) + _Rational_3_20*(f_i0m2_i1_i2 - f_i0p2_i1_i2) + _Rational_3_4*(-f_i0m1_i1_i2 + f_i0p1_i1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1
 */
static REAL _NOINLINE _UNUSED order_6_f_dD1(const REAL invdx1,const REAL f_i0_i1m3_i2,const REAL f_i0_i1m2_i2,const REAL f_i0_i1m1_i2,const REAL f_i0_i1p1_i2,const REAL f_i0_i1p2_i2,const REAL f_i0_i1p3_i2) {

      const double _Rational_3_4 = 3.0/4.0;
      const double _Rational_3_20 = 3.0/20.0;
      const double _Rational_1_60 = 1.0/60.0;
      return invdx1*(_Rational_1_60*(-f_i0_i1m3_i2 + f_i0_i1p3_i2) + _Rational_3_20*(f_i0_i1m2_i2 - f_i0_i1p2_i2) + _Rational_3_4*(-f_i0_i1m1_i2 + f_i0_i1p1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2
 */
static REAL _NOINLINE _UNUSED order_6_f_dD2(const REAL invdx2,const REAL f_i0_i1_i2m3,const REAL f_i0_i1_i2m2,const REAL f_i0_i1_i2m1,const REAL f_i0_i1_i2p1,const REAL f_i0_i1_i2p2,const REAL f_i0_i1_i2p3) {

      const double _Rational_3_4 = 3.0/4.0;
      const double _Rational_3_20 = 3.0/20.0;
      const double _Rational_1_60 = 1.0/60.0;
      return invdx2*(_Rational_1_60*(-f_i0_i1_i2m3 + f_i0_i1_i2p3) + _Rational_3_20*(f_i0_i1_i2m2 - f_i0_i1_i2p2) + _Rational_3_4*(-f_i0_i1_i2m1 + f_i0_i1_i2p1));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00
 */
static REAL _NOINLINE _UNUSED order_6_f_dDD00(const REAL invdx0,const REAL f_i0m3_i1_i2,const REAL f_i0m2_i1_i2,const REAL f_i0m1_i1_i2,const REAL f,const REAL f_i0p1_i1_i2,const REAL f_i0p2_i1_i2,const REAL f_i0p3_i1_i2) {

      const double _Rational_49_18 = 49.0/18.0;
      const double _Rational_3_20 = 3.0/20.0;
      const double _Rational_1_90 = 1.0/90.0;
      const double _Rational_3_2 = 3.0/2.0;
      return ((invdx0)*(invdx0))*(_Rational_1_90*(f_i0m3_i1_i2 + f_i0p3_i1_i2) + _Rational_3_2*(f_i0m1_i1_i2 + f_i0p1_i1_i2) + _Rational_3_20*(-f_i0m2_i1_i2 - f_i0p2_i1_i2) - _Rational_49_18*f);
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01
 */
static REAL _NOINLINE _UNUSED order_6_f_dDD01(const REAL invdx0,const REAL invdx1,const REAL f_i0m3_i1m3_i2,const REAL f_i0m2_i1m3_i2,const REAL f_i0m1_i1m3_i2,const REAL f_i0p1_i1m3_i2,const REAL f_i0p2_i1m3_i2,const REAL f_i0p3_i1m3_i2,const REAL f_i0m3_i1m2_i2,const REAL f_i0m2_i1m2_i2,const REAL f_i0m1_i1m2_i2,const REAL f_i0p1_i1m2_i2,const REAL f_i0p2_i1m2_i2,const REAL f_i0p3_i1m2_i2,const REAL f_i0m3_i1m1_i2,const REAL f_i0m2_i1m1_i2,const REAL f_i0m1_i1m1_i2,const REAL f_i0p1_i1m1_i2,const REAL f_i0p2_i1m1_i2,const REAL f_i0p3_i1m1_i2,const REAL f_i0m3_i1p1_i2,const REAL f_i0m2_i1p1_i2,const REAL f_i0m1_i1p1_i2,const REAL f_i0p1_i1p1_i2,const REAL f_i0p2_i1p1_i2,const REAL f_i0p3_i1p1_i2,const REAL f_i0m3_i1p2_i2,const REAL f_i0m2_i1p2_i2,const REAL f_i0m1_i1p2_i2,const REAL f_i0p1_i1p2_i2,const REAL f_i0p2_i1p2_i2,const REAL f_i0p3_i1p2_i2,const REAL f_i0m3_i1p3_i2,const REAL f_i0m2_i1p3_i2,const REAL f_i0m1_i1p3_i2,const REAL f_i0p1_i1p3_i2,const REAL f_i0p2_i1p3_i2,const REAL f_i0p3_i1p3_i2) {

      const double _Rational_9_16 = 9.0/16.0;
      const double _Rational_9_80 = 9.0/80.0;
      const double _Rational_9_400 = 9.0/400.0;
      const double _Rational_1_80 = 1.0/80.0;
      const double _Rational_1_400 = 1.0/400.0;
      const double _Rational_1_3600 = 1.0/3600.0;
      return invdx0*invdx1*(_Rational_1_3600*(f_i0m3_i1m3_i2 - f_i0m3_i1p3_i2 - f_i0p3_i1m3_i2 + f_i0p3_i1p3_i2) + _Rational_1_400*(-f_i0m2_i1m3_i2 + f_i0m2_i1p3_i2 - f_i0m3_i1m2_i2 + f_i0m3_i1p2_i2 + f_i0p2_i1m3_i2 - f_i0p2_i1p3_i2 + f_i0p3_i1m2_i2 - f_i0p3_i1p2_i2) + _Rational_1_80*(f_i0m1_i1m3_i2 - f_i0m1_i1p3_i2 + f_i0m3_i1m1_i2 - f_i0m3_i1p1_i2 - f_i0p1_i1m3_i2 + f_i0p1_i1p3_i2 - f_i0p3_i1m1_i2 + f_i0p3_i1p1_i2) + _Rational_9_16*(f_i0m1_i1m1_i2 - f_i0m1_i1p1_i2 - f_i0p1_i1m1_i2 + f_i0p1_i1p1_i2) + _Rational_9_400*(f_i0m2_i1m2_i2 - f_i0m2_i1p2_i2 - f_i0p2_i1m2_i2 + f_i0p2_i1p2_i2) + _Rational_9_80*(-f_i0m1_i1m2_i2 + f_i0m1_i1p2_i2 - f_i0m2_i1m1_i2 + f_i0m2_i1p1_i2 + f_i0p1_i1m2_i2 - f_i0p1_i1p2_i2 + f_i0p2_i1m1_i2 - f_i0p2_i1p1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02
 */
static REAL _NOINLINE _UNUSED order_6_f_dDD02(const REAL invdx0,const REAL invdx2,const REAL f_i0m3_i1_i2m3,const REAL f_i0m2_i1_i2m3,const REAL f_i0m1_i1_i2m3,const REAL f_i0p1_i1_i2m3,const REAL f_i0p2_i1_i2m3,const REAL f_i0p3_i1_i2m3,const REAL f_i0m3_i1_i2m2,const REAL f_i0m2_i1_i2m2,const REAL f_i0m1_i1_i2m2,const REAL f_i0p1_i1_i2m2,const REAL f_i0p2_i1_i2m2,const REAL f_i0p3_i1_i2m2,const REAL f_i0m3_i1_i2m1,const REAL f_i0m2_i1_i2m1,const REAL f_i0m1_i1_i2m1,const REAL f_i0p1_i1_i2m1,const REAL f_i0p2_i1_i2m1,const REAL f_i0p3_i1_i2m1,const REAL f_i0m3_i1_i2p1,const REAL f_i0m2_i1_i2p1,const REAL f_i0m1_i1_i2p1,const REAL f_i0p1_i1_i2p1,const REAL f_i0p2_i1_i2p1,const REAL f_i0p3_i1_i2p1,const REAL f_i0m3_i1_i2p2,const REAL f_i0m2_i1_i2p2,const REAL f_i0m1_i1_i2p2,const REAL f_i0p1_i1_i2p2,const REAL f_i0p2_i1_i2p2,const REAL f_i0p3_i1_i2p2,const REAL f_i0m3_i1_i2p3,const REAL f_i0m2_i1_i2p3,const REAL f_i0m1_i1_i2p3,const REAL f_i0p1_i1_i2p3,const REAL f_i0p2_i1_i2p3,const REAL f_i0p3_i1_i2p3) {

      const double _Rational_9_16 = 9.0/16.0;
      const double _Rational_9_80 = 9.0/80.0;
      const double _Rational_9_400 = 9.0/400.0;
      const double _Rational_1_80 = 1.0/80.0;
      const double _Rational_1_400 = 1.0/400.0;
      const double _Rational_1_3600 = 1.0/3600.0;
      return invdx0*invdx2*(_Rational_1_3600*(f_i0m3_i1_i2m3 - f_i0m3_i1_i2p3 - f_i0p3_i1_i2m3 + f_i0p3_i1_i2p3) + _Rational_1_400*(-f_i0m2_i1_i2m3 + f_i0m2_i1_i2p3 - f_i0m3_i1_i2m2 + f_i0m3_i1_i2p2 + f_i0p2_i1_i2m3 - f_i0p2_i1_i2p3 + f_i0p3_i1_i2m2 - f_i0p3_i1_i2p2) + _Rational_1_80*(f_i0m1_i1_i2m3 - f_i0m1_i1_i2p3 + f_i0m3_i1_i2m1 - f_i0m3_i1_i2p1 - f_i0p1_i1_i2m3 + f_i0p1_i1_i2p3 - f_i0p3_i1_i2m1 + f_i0p3_i1_i2p1) + _Rational_9_16*(f_i0m1_i1_i2m1 - f_i0m1_i1_i2p1 - f_i0p1_i1_i2m1 + f_i0p1_i1_i2p1) + _Rational_9_400*(f_i0m2_i1_i2m2 - f_i0m2_i1_i2p2 - f_i0p2_i1_i2m2 + f_i0p2_i1_i2p2) + _Rational_9_80*(-f_i0m1_i1_i2m2 + f_i0m1_i1_i2p2 - f_i0m2_i1_i2m1 + f_i0m2_i1_i2p1 + f_i0p1_i1_i2m2 - f_i0p1_i1_i2p2 + f_i0p2_i1_i2m1 - f_i0p2_i1_i2p1));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11
 */
static REAL _NOINLINE _UNUSED order_6_f_dDD11(const REAL invdx1,const REAL f_i0_i1m3_i2,const REAL f_i0_i1m2_i2,const REAL f_i0_i1m1_i2,const REAL f,const REAL f_i0_i1p1_i2,const REAL f_i0_i1p2_i2,const REAL f_i0_i1p3_i2) {

      const double _Rational_49_18 = 49.0/18.0;
      const double _Rational_3_20 = 3.0/20.0;
      const double _Rational_1_90 = 1.0/90.0;
      const double _Rational_3_2 = 3.0/2.0;
      return ((invdx1)*(invdx1))*(_Rational_1_90*(f_i0_i1m3_i2 + f_i0_i1p3_i2) + _Rational_3_2*(f_i0_i1m1_i2 + f_i0_i1p1_i2) + _Rational_3_20*(-f_i0_i1m2_i2 - f_i0_i1p2_i2) - _Rational_49_18*f);
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12
 */
static REAL _NOINLINE _UNUSED order_6_f_dDD12(const REAL invdx1,const REAL invdx2,const REAL f_i0_i1m3_i2m3,const REAL f_i0_i1m2_i2m3,const REAL f_i0_i1m1_i2m3,const REAL f_i0_i1p1_i2m3,const REAL f_i0_i1p2_i2m3,const REAL f_i0_i1p3_i2m3,const REAL f_i0_i1m3_i2m2,const REAL f_i0_i1m2_i2m2,const REAL f_i0_i1m1_i2m2,const REAL f_i0_i1p1_i2m2,const REAL f_i0_i1p2_i2m2,const REAL f_i0_i1p3_i2m2,const REAL f_i0_i1m3_i2m1,const REAL f_i0_i1m2_i2m1,const REAL f_i0_i1m1_i2m1,const REAL f_i0_i1p1_i2m1,const REAL f_i0_i1p2_i2m1,const REAL f_i0_i1p3_i2m1,const REAL f_i0_i1m3_i2p1,const REAL f_i0_i1m2_i2p1,const REAL f_i0_i1m1_i2p1,const REAL f_i0_i1p1_i2p1,const REAL f_i0_i1p2_i2p1,const REAL f_i0_i1p3_i2p1,const REAL f_i0_i1m3_i2p2,const REAL f_i0_i1m2_i2p2,const REAL f_i0_i1m1_i2p2,const REAL f_i0_i1p1_i2p2,const REAL f_i0_i1p2_i2p2,const REAL f_i0_i1p3_i2p2,const REAL f_i0_i1m3_i2p3,const REAL f_i0_i1m2_i2p3,const REAL f_i0_i1m1_i2p3,const REAL f_i0_i1p1_i2p3,const REAL f_i0_i1p2_i2p3,const REAL f_i0_i1p3_i2p3) {

      const double _Rational_9_16 = 9.0/16.0;
      const double _Rational_9_80 = 9.0/80.0;
      const double _Rational_9_400 = 9.0/400.0;
      const double _Rational_1_80 = 1.0/80.0;
      const double _Rational_1_400 = 1.0/400.0;
      const double _Rational_1_3600 = 1.0/3600.0;
      return invdx1*invdx2*(_Rational_1_3600*(f_i0_i1m3_i2m3 - f_i0_i1m3_i2p3 - f_i0_i1p3_i2m3 + f_i0_i1p3_i2p3) + _Rational_1_400*(-f_i0_i1m2_i2m3 + f_i0_i1m2_i2p3 - f_i0_i1m3_i2m2 + f_i0_i1m3_i2p2 + f_i0_i1p2_i2m3 - f_i0_i1p2_i2p3 + f_i0_i1p3_i2m2 - f_i0_i1p3_i2p2) + _Rational_1_80*(f_i0_i1m1_i2m3 - f_i0_i1m1_i2p3 + f_i0_i1m3_i2m1 - f_i0_i1m3_i2p1 - f_i0_i1p1_i2m3 + f_i0_i1p1_i2p3 - f_i0_i1p3_i2m1 + f_i0_i1p3_i2p1) + _Rational_9_16*(f_i0_i1m1_i2m1 - f_i0_i1m1_i2p1 - f_i0_i1p1_i2m1 + f_i0_i1p1_i2p1) + _Rational_9_400*(f_i0_i1m2_i2m2 - f_i0_i1m2_i2p2 - f_i0_i1p2_i2m2 + f_i0_i1p2_i2p2) + _Rational_9_80*(-f_i0_i1m1_i2m2 + f_i0_i1m1_i2p2 - f_i0_i1m2_i2m1 + f_i0_i1m2_i2p1 + f_i0_i1p1_i2m2 - f_i0_i1p1_i2p2 + f_i0_i1p2_i2m1 - f_i0_i1p2_i2p1));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22
 */
static REAL _NOINLINE _UNUSED order_6_f_dDD22(const REAL invdx2,const REAL f_i0_i1_i2m3,const REAL f_i0_i1_i2m2,const REAL f_i0_i1_i2m1,const REAL f,const REAL f_i0_i1_i2p1,const REAL f_i0_i1_i2p2,const REAL f_i0_i1_i2p3) {

      const double _Rational_49_18 = 49.0/18.0;
      const double _Rational_3_20 = 3.0/20.0;
      const double _Rational_1_90 = 1.0/90.0;
      const double _Rational_3_2 = 3.0/2.0;
      return ((invdx2)*(invdx2))*(_Rational_1_90*(f_i0_i1_i2m3 + f_i0_i1_i2p3) + _Rational_3_2*(f_i0_i1_i2m1 + f_i0_i1_i2p1) + _Rational_3_20*(-f_i0_i1_i2m2 - f_i0_i1_i2p2) - _Rational_49_18*f);
}
#endif // #ifndef __FD_FUNCTIONS_H__
