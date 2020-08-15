
#ifndef __FD_FUNCTIONS_H__
#define __FD_FUNCTIONS_H__
#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#define _UNUSED   __attribute__((unused))
#define _NOINLINE __attribute__((noinline))
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 0
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_ddnD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_4 = 1.0/4.0;
      const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_5_6 = 5.0/6.0;
      const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

      return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_1_4, f_i0p1_i1_i2, FusedMulAddSIMD(_Rational_5_6, f, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0m2_i1_i2, MulSIMD(_Rational_1_12, f_i0m3_i1_i2)), MulSIMD(_Rational_3_2, f_i0m1_i1_i2)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_ddnD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_4 = 1.0/4.0;
      const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_5_6 = 5.0/6.0;
      const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

      return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_4, f_i0_i1p1_i2, FusedMulAddSIMD(_Rational_5_6, f, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0_i1m2_i2, MulSIMD(_Rational_1_12, f_i0_i1m3_i2)), MulSIMD(_Rational_3_2, f_i0_i1m1_i2)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_ddnD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_4 = 1.0/4.0;
      const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_5_6 = 5.0/6.0;
      const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

      return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_4, f_i0_i1_i2p1, FusedMulAddSIMD(_Rational_5_6, f, SubSIMD(FusedMulSubSIMD(_Rational_1_2, f_i0_i1_i2m2, MulSIMD(_Rational_1_12, f_i0_i1_i2m3)), MulSIMD(_Rational_3_2, f_i0_i1_i2m1)))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 0
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dupD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_4 = 1.0/4.0;
      const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_5_6 = 5.0/6.0;
      const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

      return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_3_2, f_i0p1_i1_i2, SubSIMD(SubSIMD(FusedMulSubSIMD(_Rational_1_12, f_i0p3_i1_i2, MulSIMD(_Rational_1_2, f_i0p2_i1_i2)), MulSIMD(_Rational_5_6, f)), MulSIMD(_Rational_1_4, f_i0m1_i1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dupD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_4 = 1.0/4.0;
      const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_5_6 = 5.0/6.0;
      const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

      return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_3_2, f_i0_i1p1_i2, SubSIMD(SubSIMD(FusedMulSubSIMD(_Rational_1_12, f_i0_i1p3_i2, MulSIMD(_Rational_1_2, f_i0_i1p2_i2)), MulSIMD(_Rational_5_6, f)), MulSIMD(_Rational_1_4, f_i0_i1m1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dupD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_4 = 1.0/4.0;
      const REAL_SIMD_ARRAY _Rational_1_4 = ConstSIMD(tmp_Rational_1_4);

      const double tmp_Rational_3_2 = 3.0/2.0;
      const REAL_SIMD_ARRAY _Rational_3_2 = ConstSIMD(tmp_Rational_3_2);

      const double tmp_Rational_5_6 = 5.0/6.0;
      const REAL_SIMD_ARRAY _Rational_5_6 = ConstSIMD(tmp_Rational_5_6);

      return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_3_2, f_i0_i1_i2p1, SubSIMD(SubSIMD(FusedMulSubSIMD(_Rational_1_12, f_i0_i1_i2p3, MulSIMD(_Rational_1_2, f_i0_i1_i2p2)), MulSIMD(_Rational_5_6, f)), MulSIMD(_Rational_1_4, f_i0_i1_i2m1))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_2_3 = 2.0/3.0;
      const REAL_SIMD_ARRAY _Rational_2_3 = ConstSIMD(tmp_Rational_2_3);

      return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_1_12, SubSIMD(f_i0m2_i1_i2, f_i0p2_i1_i2), MulSIMD(_Rational_2_3, SubSIMD(f_i0p1_i1_i2, f_i0m1_i1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_2_3 = 2.0/3.0;
      const REAL_SIMD_ARRAY _Rational_2_3 = ConstSIMD(tmp_Rational_2_3);

      return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_12, SubSIMD(f_i0_i1m2_i2, f_i0_i1p2_i2), MulSIMD(_Rational_2_3, SubSIMD(f_i0_i1p1_i2, f_i0_i1m1_i2))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_2_3 = 2.0/3.0;
      const REAL_SIMD_ARRAY _Rational_2_3 = ConstSIMD(tmp_Rational_2_3);

      return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_12, SubSIMD(f_i0_i1_i2m2, f_i0_i1_i2p2), MulSIMD(_Rational_2_3, SubSIMD(f_i0_i1_i2p1, f_i0_i1_i2m1))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD00(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_5_2 = 5.0/2.0;
      const REAL_SIMD_ARRAY _Rational_5_2 = ConstSIMD(tmp_Rational_5_2);

      return MulSIMD(MulSIMD(invdx0, invdx0), FusedMulAddSIMD(_Rational_4_3, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), NegFusedMulSubSIMD(_Rational_1_12, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2), MulSIMD(_Rational_5_2, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD01(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0m2_i1m2_i2,const REAL_SIMD_ARRAY f_i0m1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p2_i1m2_i2,const REAL_SIMD_ARRAY f_i0m2_i1m1_i2,const REAL_SIMD_ARRAY f_i0m1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p2_i1m1_i2,const REAL_SIMD_ARRAY f_i0m2_i1p1_i2,const REAL_SIMD_ARRAY f_i0m1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p2_i1p1_i2,const REAL_SIMD_ARRAY f_i0m2_i1p2_i2,const REAL_SIMD_ARRAY f_i0m1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p2_i1p2_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_144 = 1.0/144.0;
      const REAL_SIMD_ARRAY _Rational_1_144 = ConstSIMD(tmp_Rational_1_144);

      const double tmp_Rational_1_18 = 1.0/18.0;
      const REAL_SIMD_ARRAY _Rational_1_18 = ConstSIMD(tmp_Rational_1_18);

      const double tmp_Rational_4_9 = 4.0/9.0;
      const REAL_SIMD_ARRAY _Rational_4_9 = ConstSIMD(tmp_Rational_4_9);

      return MulSIMD(invdx0, MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_18, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1m1_i2, f_i0m2_i1m1_i2), AddSIMD(AddSIMD(f_i0m2_i1p1_i2, f_i0p1_i1m2_i2), SubSIMD(SubSIMD(f_i0m1_i1p2_i2, f_i0m1_i1m2_i2), f_i0p2_i1p1_i2))), f_i0p1_i1p2_i2), FusedMulAddSIMD(_Rational_4_9, AddSIMD(f_i0p1_i1p1_i2, SubSIMD(SubSIMD(f_i0m1_i1m1_i2, f_i0m1_i1p1_i2), f_i0p1_i1m1_i2)), MulSIMD(_Rational_1_144, AddSIMD(f_i0p2_i1p2_i2, SubSIMD(SubSIMD(f_i0m2_i1m2_i2, f_i0m2_i1p2_i2), f_i0p2_i1m2_i2)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD02(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0m2_i1_i2m2,const REAL_SIMD_ARRAY f_i0m1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p2_i1_i2m2,const REAL_SIMD_ARRAY f_i0m2_i1_i2m1,const REAL_SIMD_ARRAY f_i0m1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p2_i1_i2m1,const REAL_SIMD_ARRAY f_i0m2_i1_i2p1,const REAL_SIMD_ARRAY f_i0m1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p2_i1_i2p1,const REAL_SIMD_ARRAY f_i0m2_i1_i2p2,const REAL_SIMD_ARRAY f_i0m1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p2_i1_i2p2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_144 = 1.0/144.0;
      const REAL_SIMD_ARRAY _Rational_1_144 = ConstSIMD(tmp_Rational_1_144);

      const double tmp_Rational_1_18 = 1.0/18.0;
      const REAL_SIMD_ARRAY _Rational_1_18 = ConstSIMD(tmp_Rational_1_18);

      const double tmp_Rational_4_9 = 4.0/9.0;
      const REAL_SIMD_ARRAY _Rational_4_9 = ConstSIMD(tmp_Rational_4_9);

      return MulSIMD(invdx0, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_18, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1_i2m1, f_i0m2_i1_i2m1), AddSIMD(AddSIMD(f_i0m2_i1_i2p1, f_i0p1_i1_i2m2), SubSIMD(SubSIMD(f_i0m1_i1_i2p2, f_i0m1_i1_i2m2), f_i0p2_i1_i2p1))), f_i0p1_i1_i2p2), FusedMulAddSIMD(_Rational_4_9, AddSIMD(f_i0p1_i1_i2p1, SubSIMD(SubSIMD(f_i0m1_i1_i2m1, f_i0m1_i1_i2p1), f_i0p1_i1_i2m1)), MulSIMD(_Rational_1_144, AddSIMD(f_i0p2_i1_i2p2, SubSIMD(SubSIMD(f_i0m2_i1_i2m2, f_i0m2_i1_i2p2), f_i0p2_i1_i2m2)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD11(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_5_2 = 5.0/2.0;
      const REAL_SIMD_ARRAY _Rational_5_2 = ConstSIMD(tmp_Rational_5_2);

      return MulSIMD(MulSIMD(invdx1, invdx1), FusedMulAddSIMD(_Rational_4_3, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), NegFusedMulSubSIMD(_Rational_1_12, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2), MulSIMD(_Rational_5_2, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD12(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1m2_i2m2,const REAL_SIMD_ARRAY f_i0_i1m1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p2_i2m2,const REAL_SIMD_ARRAY f_i0_i1m2_i2m1,const REAL_SIMD_ARRAY f_i0_i1m1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p2_i2m1,const REAL_SIMD_ARRAY f_i0_i1m2_i2p1,const REAL_SIMD_ARRAY f_i0_i1m1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p2_i2p1,const REAL_SIMD_ARRAY f_i0_i1m2_i2p2,const REAL_SIMD_ARRAY f_i0_i1m1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p2_i2p2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_144 = 1.0/144.0;
      const REAL_SIMD_ARRAY _Rational_1_144 = ConstSIMD(tmp_Rational_1_144);

      const double tmp_Rational_1_18 = 1.0/18.0;
      const REAL_SIMD_ARRAY _Rational_1_18 = ConstSIMD(tmp_Rational_1_18);

      const double tmp_Rational_4_9 = 4.0/9.0;
      const REAL_SIMD_ARRAY _Rational_4_9 = ConstSIMD(tmp_Rational_4_9);

      return MulSIMD(invdx1, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_18, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p2_i2m1, f_i0_i1m2_i2m1), AddSIMD(AddSIMD(f_i0_i1m2_i2p1, f_i0_i1p1_i2m2), SubSIMD(SubSIMD(f_i0_i1m1_i2p2, f_i0_i1m1_i2m2), f_i0_i1p2_i2p1))), f_i0_i1p1_i2p2), FusedMulAddSIMD(_Rational_4_9, AddSIMD(f_i0_i1p1_i2p1, SubSIMD(SubSIMD(f_i0_i1m1_i2m1, f_i0_i1m1_i2p1), f_i0_i1p1_i2m1)), MulSIMD(_Rational_1_144, AddSIMD(f_i0_i1p2_i2p2, SubSIMD(SubSIMD(f_i0_i1m2_i2m2, f_i0_i1m2_i2p2), f_i0_i1p2_i2m2)))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_4_f_dDD22(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_12 = 1.0/12.0;
      const REAL_SIMD_ARRAY _Rational_1_12 = ConstSIMD(tmp_Rational_1_12);

      const double tmp_Rational_4_3 = 4.0/3.0;
      const REAL_SIMD_ARRAY _Rational_4_3 = ConstSIMD(tmp_Rational_4_3);

      const double tmp_Rational_5_2 = 5.0/2.0;
      const REAL_SIMD_ARRAY _Rational_5_2 = ConstSIMD(tmp_Rational_5_2);

      return MulSIMD(MulSIMD(invdx2, invdx2), FusedMulAddSIMD(_Rational_4_3, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), NegFusedMulSubSIMD(_Rational_1_12, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2), MulSIMD(_Rational_5_2, f))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0
 */
static REAL _NOINLINE _UNUSED order_4_f_dD0(const REAL invdx0,const REAL f_i0m2_i1_i2,const REAL f_i0m1_i1_i2,const REAL f_i0p1_i1_i2,const REAL f_i0p2_i1_i2) {

      const double _Rational_2_3 = 2.0/3.0;
      const double _Rational_1_12 = 1.0/12.0;
      return invdx0*(_Rational_1_12*(f_i0m2_i1_i2 - f_i0p2_i1_i2) + _Rational_2_3*(-f_i0m1_i1_i2 + f_i0p1_i1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1
 */
static REAL _NOINLINE _UNUSED order_4_f_dD1(const REAL invdx1,const REAL f_i0_i1m2_i2,const REAL f_i0_i1m1_i2,const REAL f_i0_i1p1_i2,const REAL f_i0_i1p2_i2) {

      const double _Rational_2_3 = 2.0/3.0;
      const double _Rational_1_12 = 1.0/12.0;
      return invdx1*(_Rational_1_12*(f_i0_i1m2_i2 - f_i0_i1p2_i2) + _Rational_2_3*(-f_i0_i1m1_i2 + f_i0_i1p1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2
 */
static REAL _NOINLINE _UNUSED order_4_f_dD2(const REAL invdx2,const REAL f_i0_i1_i2m2,const REAL f_i0_i1_i2m1,const REAL f_i0_i1_i2p1,const REAL f_i0_i1_i2p2) {

      const double _Rational_2_3 = 2.0/3.0;
      const double _Rational_1_12 = 1.0/12.0;
      return invdx2*(_Rational_1_12*(f_i0_i1_i2m2 - f_i0_i1_i2p2) + _Rational_2_3*(-f_i0_i1_i2m1 + f_i0_i1_i2p1));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00
 */
static REAL _NOINLINE _UNUSED order_4_f_dDD00(const REAL invdx0,const REAL f_i0m2_i1_i2,const REAL f_i0m1_i1_i2,const REAL f,const REAL f_i0p1_i1_i2,const REAL f_i0p2_i1_i2) {

      const double _Rational_5_2 = 5.0/2.0;
      const double _Rational_1_12 = 1.0/12.0;
      const double _Rational_4_3 = 4.0/3.0;
      return ((invdx0)*(invdx0))*(_Rational_1_12*(-f_i0m2_i1_i2 - f_i0p2_i1_i2) + _Rational_4_3*(f_i0m1_i1_i2 + f_i0p1_i1_i2) - _Rational_5_2*f);
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01
 */
static REAL _NOINLINE _UNUSED order_4_f_dDD01(const REAL invdx0,const REAL invdx1,const REAL f_i0m2_i1m2_i2,const REAL f_i0m1_i1m2_i2,const REAL f_i0p1_i1m2_i2,const REAL f_i0p2_i1m2_i2,const REAL f_i0m2_i1m1_i2,const REAL f_i0m1_i1m1_i2,const REAL f_i0p1_i1m1_i2,const REAL f_i0p2_i1m1_i2,const REAL f_i0m2_i1p1_i2,const REAL f_i0m1_i1p1_i2,const REAL f_i0p1_i1p1_i2,const REAL f_i0p2_i1p1_i2,const REAL f_i0m2_i1p2_i2,const REAL f_i0m1_i1p2_i2,const REAL f_i0p1_i1p2_i2,const REAL f_i0p2_i1p2_i2) {

      const double _Rational_4_9 = 4.0/9.0;
      const double _Rational_1_18 = 1.0/18.0;
      const double _Rational_1_144 = 1.0/144.0;
      return invdx0*invdx1*(_Rational_1_144*(f_i0m2_i1m2_i2 - f_i0m2_i1p2_i2 - f_i0p2_i1m2_i2 + f_i0p2_i1p2_i2) + _Rational_1_18*(-f_i0m1_i1m2_i2 + f_i0m1_i1p2_i2 - f_i0m2_i1m1_i2 + f_i0m2_i1p1_i2 + f_i0p1_i1m2_i2 - f_i0p1_i1p2_i2 + f_i0p2_i1m1_i2 - f_i0p2_i1p1_i2) + _Rational_4_9*(f_i0m1_i1m1_i2 - f_i0m1_i1p1_i2 - f_i0p1_i1m1_i2 + f_i0p1_i1p1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02
 */
static REAL _NOINLINE _UNUSED order_4_f_dDD02(const REAL invdx0,const REAL invdx2,const REAL f_i0m2_i1_i2m2,const REAL f_i0m1_i1_i2m2,const REAL f_i0p1_i1_i2m2,const REAL f_i0p2_i1_i2m2,const REAL f_i0m2_i1_i2m1,const REAL f_i0m1_i1_i2m1,const REAL f_i0p1_i1_i2m1,const REAL f_i0p2_i1_i2m1,const REAL f_i0m2_i1_i2p1,const REAL f_i0m1_i1_i2p1,const REAL f_i0p1_i1_i2p1,const REAL f_i0p2_i1_i2p1,const REAL f_i0m2_i1_i2p2,const REAL f_i0m1_i1_i2p2,const REAL f_i0p1_i1_i2p2,const REAL f_i0p2_i1_i2p2) {

      const double _Rational_4_9 = 4.0/9.0;
      const double _Rational_1_18 = 1.0/18.0;
      const double _Rational_1_144 = 1.0/144.0;
      return invdx0*invdx2*(_Rational_1_144*(f_i0m2_i1_i2m2 - f_i0m2_i1_i2p2 - f_i0p2_i1_i2m2 + f_i0p2_i1_i2p2) + _Rational_1_18*(-f_i0m1_i1_i2m2 + f_i0m1_i1_i2p2 - f_i0m2_i1_i2m1 + f_i0m2_i1_i2p1 + f_i0p1_i1_i2m2 - f_i0p1_i1_i2p2 + f_i0p2_i1_i2m1 - f_i0p2_i1_i2p1) + _Rational_4_9*(f_i0m1_i1_i2m1 - f_i0m1_i1_i2p1 - f_i0p1_i1_i2m1 + f_i0p1_i1_i2p1));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11
 */
static REAL _NOINLINE _UNUSED order_4_f_dDD11(const REAL invdx1,const REAL f_i0_i1m2_i2,const REAL f_i0_i1m1_i2,const REAL f,const REAL f_i0_i1p1_i2,const REAL f_i0_i1p2_i2) {

      const double _Rational_5_2 = 5.0/2.0;
      const double _Rational_1_12 = 1.0/12.0;
      const double _Rational_4_3 = 4.0/3.0;
      return ((invdx1)*(invdx1))*(_Rational_1_12*(-f_i0_i1m2_i2 - f_i0_i1p2_i2) + _Rational_4_3*(f_i0_i1m1_i2 + f_i0_i1p1_i2) - _Rational_5_2*f);
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12
 */
static REAL _NOINLINE _UNUSED order_4_f_dDD12(const REAL invdx1,const REAL invdx2,const REAL f_i0_i1m2_i2m2,const REAL f_i0_i1m1_i2m2,const REAL f_i0_i1p1_i2m2,const REAL f_i0_i1p2_i2m2,const REAL f_i0_i1m2_i2m1,const REAL f_i0_i1m1_i2m1,const REAL f_i0_i1p1_i2m1,const REAL f_i0_i1p2_i2m1,const REAL f_i0_i1m2_i2p1,const REAL f_i0_i1m1_i2p1,const REAL f_i0_i1p1_i2p1,const REAL f_i0_i1p2_i2p1,const REAL f_i0_i1m2_i2p2,const REAL f_i0_i1m1_i2p2,const REAL f_i0_i1p1_i2p2,const REAL f_i0_i1p2_i2p2) {

      const double _Rational_4_9 = 4.0/9.0;
      const double _Rational_1_18 = 1.0/18.0;
      const double _Rational_1_144 = 1.0/144.0;
      return invdx1*invdx2*(_Rational_1_144*(f_i0_i1m2_i2m2 - f_i0_i1m2_i2p2 - f_i0_i1p2_i2m2 + f_i0_i1p2_i2p2) + _Rational_1_18*(-f_i0_i1m1_i2m2 + f_i0_i1m1_i2p2 - f_i0_i1m2_i2m1 + f_i0_i1m2_i2p1 + f_i0_i1p1_i2m2 - f_i0_i1p1_i2p2 + f_i0_i1p2_i2m1 - f_i0_i1p2_i2p1) + _Rational_4_9*(f_i0_i1m1_i2m1 - f_i0_i1m1_i2p1 - f_i0_i1p1_i2m1 + f_i0_i1p1_i2p1));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22
 */
static REAL _NOINLINE _UNUSED order_4_f_dDD22(const REAL invdx2,const REAL f_i0_i1_i2m2,const REAL f_i0_i1_i2m1,const REAL f,const REAL f_i0_i1_i2p1,const REAL f_i0_i1_i2p2) {

      const double _Rational_5_2 = 5.0/2.0;
      const double _Rational_1_12 = 1.0/12.0;
      const double _Rational_4_3 = 4.0/3.0;
      return ((invdx2)*(invdx2))*(_Rational_1_12*(-f_i0_i1_i2m2 - f_i0_i1_i2p2) + _Rational_4_3*(f_i0_i1_i2m1 + f_i0_i1_i2p1) - _Rational_5_2*f);
}
#endif // #ifndef __FD_FUNCTIONS_H__
