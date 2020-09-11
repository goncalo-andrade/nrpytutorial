
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
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_ddnD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m5_i1_i2,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_14 = 1.0/14.0;
      const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

      const double tmp_Rational_1_168 = 1.0/168.0;
      const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_28 = 1.0/28.0;
      const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_6 = 1.0/6.0;
      const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

      const double tmp_Rational_5_4 = 5.0/4.0;
      const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

      const double tmp_Rational_9_20 = 9.0/20.0;
      const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

      return MulSIMD(invdx0, SubSIMD(FusedMulAddSIMD(_Rational_9_20, f, SubSIMD(FusedMulAddSIMD(_Rational_1_2, AddSIMD(f_i0m2_i1_i2, f_i0p1_i1_i2), FusedMulAddSIMD(_Rational_1_28, f_i0m4_i1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_168, f_i0p3_i1_i2, MulSIMD(_Rational_1_14, f_i0p2_i1_i2)), MulSIMD(_Rational_5_4, f_i0m1_i1_i2)))), MulSIMD(_Rational_1_280, f_i0m5_i1_i2))), MulSIMD(_Rational_1_6, f_i0m3_i1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_ddnD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m5_i2,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_14 = 1.0/14.0;
      const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

      const double tmp_Rational_1_168 = 1.0/168.0;
      const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_28 = 1.0/28.0;
      const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_6 = 1.0/6.0;
      const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

      const double tmp_Rational_5_4 = 5.0/4.0;
      const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

      const double tmp_Rational_9_20 = 9.0/20.0;
      const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

      return MulSIMD(invdx1, SubSIMD(FusedMulAddSIMD(_Rational_9_20, f, SubSIMD(FusedMulAddSIMD(_Rational_1_2, AddSIMD(f_i0_i1m2_i2, f_i0_i1p1_i2), FusedMulAddSIMD(_Rational_1_28, f_i0_i1m4_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_168, f_i0_i1p3_i2, MulSIMD(_Rational_1_14, f_i0_i1p2_i2)), MulSIMD(_Rational_5_4, f_i0_i1m1_i2)))), MulSIMD(_Rational_1_280, f_i0_i1m5_i2))), MulSIMD(_Rational_1_6, f_i0_i1m3_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for downwinded derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_ddnD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m5,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_14 = 1.0/14.0;
      const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

      const double tmp_Rational_1_168 = 1.0/168.0;
      const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_28 = 1.0/28.0;
      const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_6 = 1.0/6.0;
      const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

      const double tmp_Rational_5_4 = 5.0/4.0;
      const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

      const double tmp_Rational_9_20 = 9.0/20.0;
      const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

      return MulSIMD(invdx2, SubSIMD(FusedMulAddSIMD(_Rational_9_20, f, SubSIMD(FusedMulAddSIMD(_Rational_1_2, AddSIMD(f_i0_i1_i2m2, f_i0_i1_i2p1), FusedMulAddSIMD(_Rational_1_28, f_i0_i1_i2m4, SubSIMD(FusedMulSubSIMD(_Rational_1_168, f_i0_i1_i2p3, MulSIMD(_Rational_1_14, f_i0_i1_i2p2)), MulSIMD(_Rational_5_4, f_i0_i1_i2m1)))), MulSIMD(_Rational_1_280, f_i0_i1_i2m5))), MulSIMD(_Rational_1_6, f_i0_i1_i2m3)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 0
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dupD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2,const REAL_SIMD_ARRAY f_i0p5_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_14 = 1.0/14.0;
      const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

      const double tmp_Rational_1_168 = 1.0/168.0;
      const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_28 = 1.0/28.0;
      const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_6 = 1.0/6.0;
      const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

      const double tmp_Rational_5_4 = 5.0/4.0;
      const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

      const double tmp_Rational_9_20 = 9.0/20.0;
      const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

      return MulSIMD(invdx0, SubSIMD(FusedMulAddSIMD(_Rational_1_6, f_i0p3_i1_i2, FusedMulAddSIMD(_Rational_5_4, f_i0p1_i1_i2, FusedMulAddSIMD(_Rational_1_2, MulSIMD(_NegativeOne_, AddSIMD(f_i0p2_i1_i2, f_i0m1_i1_i2)), FusedMulAddSIMD(_Rational_1_280, f_i0p5_i1_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_14, f_i0m2_i1_i2, MulSIMD(_Rational_1_168, f_i0m3_i1_i2)), MulSIMD(_Rational_9_20, f)))))), MulSIMD(_Rational_1_28, f_i0p4_i1_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dupD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2,const REAL_SIMD_ARRAY f_i0_i1p5_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_14 = 1.0/14.0;
      const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

      const double tmp_Rational_1_168 = 1.0/168.0;
      const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_28 = 1.0/28.0;
      const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_6 = 1.0/6.0;
      const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

      const double tmp_Rational_5_4 = 5.0/4.0;
      const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

      const double tmp_Rational_9_20 = 9.0/20.0;
      const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

      return MulSIMD(invdx1, SubSIMD(FusedMulAddSIMD(_Rational_1_6, f_i0_i1p3_i2, FusedMulAddSIMD(_Rational_5_4, f_i0_i1p1_i2, FusedMulAddSIMD(_Rational_1_2, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p2_i2, f_i0_i1m1_i2)), FusedMulAddSIMD(_Rational_1_280, f_i0_i1p5_i2, SubSIMD(FusedMulSubSIMD(_Rational_1_14, f_i0_i1m2_i2, MulSIMD(_Rational_1_168, f_i0_i1m3_i2)), MulSIMD(_Rational_9_20, f)))))), MulSIMD(_Rational_1_28, f_i0_i1p4_i2)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for upwinded derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dupD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4,const REAL_SIMD_ARRAY f_i0_i1_i2p5) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_14 = 1.0/14.0;
      const REAL_SIMD_ARRAY _Rational_1_14 = ConstSIMD(tmp_Rational_1_14);

      const double tmp_Rational_1_168 = 1.0/168.0;
      const REAL_SIMD_ARRAY _Rational_1_168 = ConstSIMD(tmp_Rational_1_168);

      const double tmp_Rational_1_2 = 1.0/2.0;
      const REAL_SIMD_ARRAY _Rational_1_2 = ConstSIMD(tmp_Rational_1_2);

      const double tmp_Rational_1_28 = 1.0/28.0;
      const REAL_SIMD_ARRAY _Rational_1_28 = ConstSIMD(tmp_Rational_1_28);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_6 = 1.0/6.0;
      const REAL_SIMD_ARRAY _Rational_1_6 = ConstSIMD(tmp_Rational_1_6);

      const double tmp_Rational_5_4 = 5.0/4.0;
      const REAL_SIMD_ARRAY _Rational_5_4 = ConstSIMD(tmp_Rational_5_4);

      const double tmp_Rational_9_20 = 9.0/20.0;
      const REAL_SIMD_ARRAY _Rational_9_20 = ConstSIMD(tmp_Rational_9_20);

      return MulSIMD(invdx2, SubSIMD(FusedMulAddSIMD(_Rational_1_6, f_i0_i1_i2p3, FusedMulAddSIMD(_Rational_5_4, f_i0_i1_i2p1, FusedMulAddSIMD(_Rational_1_2, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m1)), FusedMulAddSIMD(_Rational_1_280, f_i0_i1_i2p5, SubSIMD(FusedMulSubSIMD(_Rational_1_14, f_i0_i1_i2m2, MulSIMD(_Rational_1_168, f_i0_i1_i2m3)), MulSIMD(_Rational_9_20, f)))))), MulSIMD(_Rational_1_28, f_i0_i1_i2p4)));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dD0(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_5 = 1.0/5.0;
      const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

      const double tmp_Rational_4_105 = 4.0/105.0;
      const REAL_SIMD_ARRAY _Rational_4_105 = ConstSIMD(tmp_Rational_4_105);

      const double tmp_Rational_4_5 = 4.0/5.0;
      const REAL_SIMD_ARRAY _Rational_4_5 = ConstSIMD(tmp_Rational_4_5);

      return MulSIMD(invdx0, FusedMulAddSIMD(_Rational_4_105, SubSIMD(f_i0p3_i1_i2, f_i0m3_i1_i2), FusedMulAddSIMD(_Rational_4_5, SubSIMD(f_i0p1_i1_i2, f_i0m1_i1_i2), FusedMulAddSIMD(_Rational_1_280, SubSIMD(f_i0m4_i1_i2, f_i0p4_i1_i2), MulSIMD(_Rational_1_5, SubSIMD(f_i0m2_i1_i2, f_i0p2_i1_i2))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dD1(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_5 = 1.0/5.0;
      const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

      const double tmp_Rational_4_105 = 4.0/105.0;
      const REAL_SIMD_ARRAY _Rational_4_105 = ConstSIMD(tmp_Rational_4_105);

      const double tmp_Rational_4_5 = 4.0/5.0;
      const REAL_SIMD_ARRAY _Rational_4_5 = ConstSIMD(tmp_Rational_4_5);

      return MulSIMD(invdx1, FusedMulAddSIMD(_Rational_4_105, SubSIMD(f_i0_i1p3_i2, f_i0_i1m3_i2), FusedMulAddSIMD(_Rational_4_5, SubSIMD(f_i0_i1p1_i2, f_i0_i1m1_i2), FusedMulAddSIMD(_Rational_1_280, SubSIMD(f_i0_i1m4_i2, f_i0_i1p4_i2), MulSIMD(_Rational_1_5, SubSIMD(f_i0_i1m2_i2, f_i0_i1p2_i2))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dD2(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_280 = 1.0/280.0;
      const REAL_SIMD_ARRAY _Rational_1_280 = ConstSIMD(tmp_Rational_1_280);

      const double tmp_Rational_1_5 = 1.0/5.0;
      const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

      const double tmp_Rational_4_105 = 4.0/105.0;
      const REAL_SIMD_ARRAY _Rational_4_105 = ConstSIMD(tmp_Rational_4_105);

      const double tmp_Rational_4_5 = 4.0/5.0;
      const REAL_SIMD_ARRAY _Rational_4_5 = ConstSIMD(tmp_Rational_4_5);

      return MulSIMD(invdx2, FusedMulAddSIMD(_Rational_4_105, SubSIMD(f_i0_i1_i2p3, f_i0_i1_i2m3), FusedMulAddSIMD(_Rational_4_5, SubSIMD(f_i0_i1_i2p1, f_i0_i1_i2m1), FusedMulAddSIMD(_Rational_1_280, SubSIMD(f_i0_i1_i2m4, f_i0_i1_i2p4), MulSIMD(_Rational_1_5, SubSIMD(f_i0_i1_i2m2, f_i0_i1_i2p2))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD00(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY f_i0m4_i1_i2,const REAL_SIMD_ARRAY f_i0m3_i1_i2,const REAL_SIMD_ARRAY f_i0m2_i1_i2,const REAL_SIMD_ARRAY f_i0m1_i1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0p1_i1_i2,const REAL_SIMD_ARRAY f_i0p2_i1_i2,const REAL_SIMD_ARRAY f_i0p3_i1_i2,const REAL_SIMD_ARRAY f_i0p4_i1_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_5 = 1.0/5.0;
      const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

      const double tmp_Rational_1_560 = 1.0/560.0;
      const REAL_SIMD_ARRAY _Rational_1_560 = ConstSIMD(tmp_Rational_1_560);

      const double tmp_Rational_205_72 = 205.0/72.0;
      const REAL_SIMD_ARRAY _Rational_205_72 = ConstSIMD(tmp_Rational_205_72);

      const double tmp_Rational_8_315 = 8.0/315.0;
      const REAL_SIMD_ARRAY _Rational_8_315 = ConstSIMD(tmp_Rational_8_315);

      const double tmp_Rational_8_5 = 8.0/5.0;
      const REAL_SIMD_ARRAY _Rational_8_5 = ConstSIMD(tmp_Rational_8_5);

      return MulSIMD(MulSIMD(invdx0, invdx0), FusedMulAddSIMD(_Rational_1_560, MulSIMD(_NegativeOne_, AddSIMD(f_i0p4_i1_i2, f_i0m4_i1_i2)), FusedMulAddSIMD(_Rational_8_315, AddSIMD(f_i0m3_i1_i2, f_i0p3_i1_i2), FusedMulAddSIMD(_Rational_8_5, AddSIMD(f_i0m1_i1_i2, f_i0p1_i1_i2), NegFusedMulSubSIMD(_Rational_1_5, AddSIMD(f_i0p2_i1_i2, f_i0m2_i1_i2), MulSIMD(_Rational_205_72, f))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD01(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0m4_i1m4_i2,const REAL_SIMD_ARRAY f_i0m3_i1m4_i2,const REAL_SIMD_ARRAY f_i0m2_i1m4_i2,const REAL_SIMD_ARRAY f_i0m1_i1m4_i2,const REAL_SIMD_ARRAY f_i0p1_i1m4_i2,const REAL_SIMD_ARRAY f_i0p2_i1m4_i2,const REAL_SIMD_ARRAY f_i0p3_i1m4_i2,const REAL_SIMD_ARRAY f_i0p4_i1m4_i2,const REAL_SIMD_ARRAY f_i0m4_i1m3_i2,const REAL_SIMD_ARRAY f_i0m3_i1m3_i2,const REAL_SIMD_ARRAY f_i0m2_i1m3_i2,const REAL_SIMD_ARRAY f_i0m1_i1m3_i2,const REAL_SIMD_ARRAY f_i0p1_i1m3_i2,const REAL_SIMD_ARRAY f_i0p2_i1m3_i2,const REAL_SIMD_ARRAY f_i0p3_i1m3_i2,const REAL_SIMD_ARRAY f_i0p4_i1m3_i2,const REAL_SIMD_ARRAY f_i0m4_i1m2_i2,const REAL_SIMD_ARRAY f_i0m3_i1m2_i2,const REAL_SIMD_ARRAY f_i0m2_i1m2_i2,const REAL_SIMD_ARRAY f_i0m1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p1_i1m2_i2,const REAL_SIMD_ARRAY f_i0p2_i1m2_i2,const REAL_SIMD_ARRAY f_i0p3_i1m2_i2,const REAL_SIMD_ARRAY f_i0p4_i1m2_i2,const REAL_SIMD_ARRAY f_i0m4_i1m1_i2,const REAL_SIMD_ARRAY f_i0m3_i1m1_i2,const REAL_SIMD_ARRAY f_i0m2_i1m1_i2,const REAL_SIMD_ARRAY f_i0m1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p1_i1m1_i2,const REAL_SIMD_ARRAY f_i0p2_i1m1_i2,const REAL_SIMD_ARRAY f_i0p3_i1m1_i2,const REAL_SIMD_ARRAY f_i0p4_i1m1_i2,const REAL_SIMD_ARRAY f_i0m4_i1p1_i2,const REAL_SIMD_ARRAY f_i0m3_i1p1_i2,const REAL_SIMD_ARRAY f_i0m2_i1p1_i2,const REAL_SIMD_ARRAY f_i0m1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p1_i1p1_i2,const REAL_SIMD_ARRAY f_i0p2_i1p1_i2,const REAL_SIMD_ARRAY f_i0p3_i1p1_i2,const REAL_SIMD_ARRAY f_i0p4_i1p1_i2,const REAL_SIMD_ARRAY f_i0m4_i1p2_i2,const REAL_SIMD_ARRAY f_i0m3_i1p2_i2,const REAL_SIMD_ARRAY f_i0m2_i1p2_i2,const REAL_SIMD_ARRAY f_i0m1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p1_i1p2_i2,const REAL_SIMD_ARRAY f_i0p2_i1p2_i2,const REAL_SIMD_ARRAY f_i0p3_i1p2_i2,const REAL_SIMD_ARRAY f_i0p4_i1p2_i2,const REAL_SIMD_ARRAY f_i0m4_i1p3_i2,const REAL_SIMD_ARRAY f_i0m3_i1p3_i2,const REAL_SIMD_ARRAY f_i0m2_i1p3_i2,const REAL_SIMD_ARRAY f_i0m1_i1p3_i2,const REAL_SIMD_ARRAY f_i0p1_i1p3_i2,const REAL_SIMD_ARRAY f_i0p2_i1p3_i2,const REAL_SIMD_ARRAY f_i0p3_i1p3_i2,const REAL_SIMD_ARRAY f_i0p4_i1p3_i2,const REAL_SIMD_ARRAY f_i0m4_i1p4_i2,const REAL_SIMD_ARRAY f_i0m3_i1p4_i2,const REAL_SIMD_ARRAY f_i0m2_i1p4_i2,const REAL_SIMD_ARRAY f_i0m1_i1p4_i2,const REAL_SIMD_ARRAY f_i0p1_i1p4_i2,const REAL_SIMD_ARRAY f_i0p2_i1p4_i2,const REAL_SIMD_ARRAY f_i0p3_i1p4_i2,const REAL_SIMD_ARRAY f_i0p4_i1p4_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_16_11025 = 16.0/11025.0;
      const REAL_SIMD_ARRAY _Rational_16_11025 = ConstSIMD(tmp_Rational_16_11025);

      const double tmp_Rational_16_25 = 16.0/25.0;
      const REAL_SIMD_ARRAY _Rational_16_25 = ConstSIMD(tmp_Rational_16_25);

      const double tmp_Rational_16_525 = 16.0/525.0;
      const REAL_SIMD_ARRAY _Rational_16_525 = ConstSIMD(tmp_Rational_16_525);

      const double tmp_Rational_1_1400 = 1.0/1400.0;
      const REAL_SIMD_ARRAY _Rational_1_1400 = ConstSIMD(tmp_Rational_1_1400);

      const double tmp_Rational_1_25 = 1.0/25.0;
      const REAL_SIMD_ARRAY _Rational_1_25 = ConstSIMD(tmp_Rational_1_25);

      const double tmp_Rational_1_350 = 1.0/350.0;
      const REAL_SIMD_ARRAY _Rational_1_350 = ConstSIMD(tmp_Rational_1_350);

      const double tmp_Rational_1_7350 = 1.0/7350.0;
      const REAL_SIMD_ARRAY _Rational_1_7350 = ConstSIMD(tmp_Rational_1_7350);

      const double tmp_Rational_1_78400 = 1.0/78400.0;
      const REAL_SIMD_ARRAY _Rational_1_78400 = ConstSIMD(tmp_Rational_1_78400);

      const double tmp_Rational_4_25 = 4.0/25.0;
      const REAL_SIMD_ARRAY _Rational_4_25 = ConstSIMD(tmp_Rational_4_25);

      const double tmp_Rational_4_525 = 4.0/525.0;
      const REAL_SIMD_ARRAY _Rational_4_525 = ConstSIMD(tmp_Rational_4_525);

      return MulSIMD(invdx0, MulSIMD(invdx1, FusedMulAddSIMD(_Rational_1_7350, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1m3_i2, f_i0m4_i1m3_i2), AddSIMD(AddSIMD(f_i0m4_i1p3_i2, f_i0p3_i1m4_i2), SubSIMD(SubSIMD(f_i0m3_i1p4_i2, f_i0m3_i1m4_i2), f_i0p4_i1p3_i2))), f_i0p3_i1p4_i2), FusedMulAddSIMD(_Rational_1_78400, AddSIMD(f_i0p4_i1p4_i2, SubSIMD(SubSIMD(f_i0m4_i1m4_i2, f_i0m4_i1p4_i2), f_i0p4_i1m4_i2)), FusedMulAddSIMD(_Rational_1_25, AddSIMD(f_i0p2_i1p2_i2, SubSIMD(SubSIMD(f_i0m2_i1m2_i2, f_i0m2_i1p2_i2), f_i0p2_i1m2_i2)), FusedMulAddSIMD(_Rational_1_350, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1m1_i2, f_i0m4_i1m1_i2), AddSIMD(AddSIMD(f_i0m4_i1p1_i2, f_i0p1_i1m4_i2), SubSIMD(SubSIMD(f_i0m1_i1p4_i2, f_i0m1_i1m4_i2), f_i0p4_i1p1_i2))), f_i0p1_i1p4_i2), FusedMulAddSIMD(_Rational_16_525, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1p1_i2, f_i0m3_i1p1_i2), AddSIMD(AddSIMD(f_i0m3_i1m1_i2, f_i0p1_i1p3_i2), SubSIMD(SubSIMD(f_i0m1_i1m3_i2, f_i0m1_i1p3_i2), f_i0p3_i1m1_i2))), f_i0p1_i1m3_i2), FusedMulAddSIMD(_Rational_1_1400, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1p2_i2, f_i0m4_i1p2_i2), AddSIMD(AddSIMD(f_i0m4_i1m2_i2, f_i0p2_i1p4_i2), SubSIMD(SubSIMD(f_i0m2_i1m4_i2, f_i0m2_i1p4_i2), f_i0p4_i1m2_i2))), f_i0p2_i1m4_i2), FusedMulAddSIMD(_Rational_4_25, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1m1_i2, f_i0m2_i1m1_i2), AddSIMD(AddSIMD(f_i0m2_i1p1_i2, f_i0p1_i1m2_i2), SubSIMD(SubSIMD(f_i0m1_i1p2_i2, f_i0m1_i1m2_i2), f_i0p2_i1p1_i2))), f_i0p1_i1p2_i2), FusedMulAddSIMD(_Rational_4_525, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1m2_i2, f_i0m3_i1m2_i2), AddSIMD(AddSIMD(f_i0m3_i1p2_i2, f_i0p2_i1m3_i2), SubSIMD(SubSIMD(f_i0m2_i1p3_i2, f_i0m2_i1m3_i2), f_i0p3_i1p2_i2))), f_i0p2_i1p3_i2), FusedMulAddSIMD(_Rational_16_11025, AddSIMD(f_i0p3_i1p3_i2, SubSIMD(SubSIMD(f_i0m3_i1m3_i2, f_i0m3_i1p3_i2), f_i0p3_i1m3_i2)), MulSIMD(_Rational_16_25, AddSIMD(f_i0p1_i1p1_i2, SubSIMD(SubSIMD(f_i0m1_i1m1_i2, f_i0m1_i1p1_i2), f_i0p1_i1m1_i2))))))))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD02(const REAL_SIMD_ARRAY invdx0,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0m4_i1_i2m4,const REAL_SIMD_ARRAY f_i0m3_i1_i2m4,const REAL_SIMD_ARRAY f_i0m2_i1_i2m4,const REAL_SIMD_ARRAY f_i0m1_i1_i2m4,const REAL_SIMD_ARRAY f_i0p1_i1_i2m4,const REAL_SIMD_ARRAY f_i0p2_i1_i2m4,const REAL_SIMD_ARRAY f_i0p3_i1_i2m4,const REAL_SIMD_ARRAY f_i0p4_i1_i2m4,const REAL_SIMD_ARRAY f_i0m4_i1_i2m3,const REAL_SIMD_ARRAY f_i0m3_i1_i2m3,const REAL_SIMD_ARRAY f_i0m2_i1_i2m3,const REAL_SIMD_ARRAY f_i0m1_i1_i2m3,const REAL_SIMD_ARRAY f_i0p1_i1_i2m3,const REAL_SIMD_ARRAY f_i0p2_i1_i2m3,const REAL_SIMD_ARRAY f_i0p3_i1_i2m3,const REAL_SIMD_ARRAY f_i0p4_i1_i2m3,const REAL_SIMD_ARRAY f_i0m4_i1_i2m2,const REAL_SIMD_ARRAY f_i0m3_i1_i2m2,const REAL_SIMD_ARRAY f_i0m2_i1_i2m2,const REAL_SIMD_ARRAY f_i0m1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p1_i1_i2m2,const REAL_SIMD_ARRAY f_i0p2_i1_i2m2,const REAL_SIMD_ARRAY f_i0p3_i1_i2m2,const REAL_SIMD_ARRAY f_i0p4_i1_i2m2,const REAL_SIMD_ARRAY f_i0m4_i1_i2m1,const REAL_SIMD_ARRAY f_i0m3_i1_i2m1,const REAL_SIMD_ARRAY f_i0m2_i1_i2m1,const REAL_SIMD_ARRAY f_i0m1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p1_i1_i2m1,const REAL_SIMD_ARRAY f_i0p2_i1_i2m1,const REAL_SIMD_ARRAY f_i0p3_i1_i2m1,const REAL_SIMD_ARRAY f_i0p4_i1_i2m1,const REAL_SIMD_ARRAY f_i0m4_i1_i2p1,const REAL_SIMD_ARRAY f_i0m3_i1_i2p1,const REAL_SIMD_ARRAY f_i0m2_i1_i2p1,const REAL_SIMD_ARRAY f_i0m1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p1_i1_i2p1,const REAL_SIMD_ARRAY f_i0p2_i1_i2p1,const REAL_SIMD_ARRAY f_i0p3_i1_i2p1,const REAL_SIMD_ARRAY f_i0p4_i1_i2p1,const REAL_SIMD_ARRAY f_i0m4_i1_i2p2,const REAL_SIMD_ARRAY f_i0m3_i1_i2p2,const REAL_SIMD_ARRAY f_i0m2_i1_i2p2,const REAL_SIMD_ARRAY f_i0m1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p1_i1_i2p2,const REAL_SIMD_ARRAY f_i0p2_i1_i2p2,const REAL_SIMD_ARRAY f_i0p3_i1_i2p2,const REAL_SIMD_ARRAY f_i0p4_i1_i2p2,const REAL_SIMD_ARRAY f_i0m4_i1_i2p3,const REAL_SIMD_ARRAY f_i0m3_i1_i2p3,const REAL_SIMD_ARRAY f_i0m2_i1_i2p3,const REAL_SIMD_ARRAY f_i0m1_i1_i2p3,const REAL_SIMD_ARRAY f_i0p1_i1_i2p3,const REAL_SIMD_ARRAY f_i0p2_i1_i2p3,const REAL_SIMD_ARRAY f_i0p3_i1_i2p3,const REAL_SIMD_ARRAY f_i0p4_i1_i2p3,const REAL_SIMD_ARRAY f_i0m4_i1_i2p4,const REAL_SIMD_ARRAY f_i0m3_i1_i2p4,const REAL_SIMD_ARRAY f_i0m2_i1_i2p4,const REAL_SIMD_ARRAY f_i0m1_i1_i2p4,const REAL_SIMD_ARRAY f_i0p1_i1_i2p4,const REAL_SIMD_ARRAY f_i0p2_i1_i2p4,const REAL_SIMD_ARRAY f_i0p3_i1_i2p4,const REAL_SIMD_ARRAY f_i0p4_i1_i2p4) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_16_11025 = 16.0/11025.0;
      const REAL_SIMD_ARRAY _Rational_16_11025 = ConstSIMD(tmp_Rational_16_11025);

      const double tmp_Rational_16_25 = 16.0/25.0;
      const REAL_SIMD_ARRAY _Rational_16_25 = ConstSIMD(tmp_Rational_16_25);

      const double tmp_Rational_16_525 = 16.0/525.0;
      const REAL_SIMD_ARRAY _Rational_16_525 = ConstSIMD(tmp_Rational_16_525);

      const double tmp_Rational_1_1400 = 1.0/1400.0;
      const REAL_SIMD_ARRAY _Rational_1_1400 = ConstSIMD(tmp_Rational_1_1400);

      const double tmp_Rational_1_25 = 1.0/25.0;
      const REAL_SIMD_ARRAY _Rational_1_25 = ConstSIMD(tmp_Rational_1_25);

      const double tmp_Rational_1_350 = 1.0/350.0;
      const REAL_SIMD_ARRAY _Rational_1_350 = ConstSIMD(tmp_Rational_1_350);

      const double tmp_Rational_1_7350 = 1.0/7350.0;
      const REAL_SIMD_ARRAY _Rational_1_7350 = ConstSIMD(tmp_Rational_1_7350);

      const double tmp_Rational_1_78400 = 1.0/78400.0;
      const REAL_SIMD_ARRAY _Rational_1_78400 = ConstSIMD(tmp_Rational_1_78400);

      const double tmp_Rational_4_25 = 4.0/25.0;
      const REAL_SIMD_ARRAY _Rational_4_25 = ConstSIMD(tmp_Rational_4_25);

      const double tmp_Rational_4_525 = 4.0/525.0;
      const REAL_SIMD_ARRAY _Rational_4_525 = ConstSIMD(tmp_Rational_4_525);

      return MulSIMD(invdx0, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_7350, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1_i2m3, f_i0m4_i1_i2m3), AddSIMD(AddSIMD(f_i0m4_i1_i2p3, f_i0p3_i1_i2m4), SubSIMD(SubSIMD(f_i0m3_i1_i2p4, f_i0m3_i1_i2m4), f_i0p4_i1_i2p3))), f_i0p3_i1_i2p4), FusedMulAddSIMD(_Rational_1_78400, AddSIMD(f_i0p4_i1_i2p4, SubSIMD(SubSIMD(f_i0m4_i1_i2m4, f_i0m4_i1_i2p4), f_i0p4_i1_i2m4)), FusedMulAddSIMD(_Rational_1_25, AddSIMD(f_i0p2_i1_i2p2, SubSIMD(SubSIMD(f_i0m2_i1_i2m2, f_i0m2_i1_i2p2), f_i0p2_i1_i2m2)), FusedMulAddSIMD(_Rational_1_350, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1_i2m1, f_i0m4_i1_i2m1), AddSIMD(AddSIMD(f_i0m4_i1_i2p1, f_i0p1_i1_i2m4), SubSIMD(SubSIMD(f_i0m1_i1_i2p4, f_i0m1_i1_i2m4), f_i0p4_i1_i2p1))), f_i0p1_i1_i2p4), FusedMulAddSIMD(_Rational_16_525, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1_i2p1, f_i0m3_i1_i2p1), AddSIMD(AddSIMD(f_i0m3_i1_i2m1, f_i0p1_i1_i2p3), SubSIMD(SubSIMD(f_i0m1_i1_i2m3, f_i0m1_i1_i2p3), f_i0p3_i1_i2m1))), f_i0p1_i1_i2m3), FusedMulAddSIMD(_Rational_1_1400, SubSIMD(AddSIMD(SubSIMD(f_i0p4_i1_i2p2, f_i0m4_i1_i2p2), AddSIMD(AddSIMD(f_i0m4_i1_i2m2, f_i0p2_i1_i2p4), SubSIMD(SubSIMD(f_i0m2_i1_i2m4, f_i0m2_i1_i2p4), f_i0p4_i1_i2m2))), f_i0p2_i1_i2m4), FusedMulAddSIMD(_Rational_4_25, SubSIMD(AddSIMD(SubSIMD(f_i0p2_i1_i2m1, f_i0m2_i1_i2m1), AddSIMD(AddSIMD(f_i0m2_i1_i2p1, f_i0p1_i1_i2m2), SubSIMD(SubSIMD(f_i0m1_i1_i2p2, f_i0m1_i1_i2m2), f_i0p2_i1_i2p1))), f_i0p1_i1_i2p2), FusedMulAddSIMD(_Rational_4_525, SubSIMD(AddSIMD(SubSIMD(f_i0p3_i1_i2m2, f_i0m3_i1_i2m2), AddSIMD(AddSIMD(f_i0m3_i1_i2p2, f_i0p2_i1_i2m3), SubSIMD(SubSIMD(f_i0m2_i1_i2p3, f_i0m2_i1_i2m3), f_i0p3_i1_i2p2))), f_i0p2_i1_i2p3), FusedMulAddSIMD(_Rational_16_11025, AddSIMD(f_i0p3_i1_i2p3, SubSIMD(SubSIMD(f_i0m3_i1_i2m3, f_i0m3_i1_i2p3), f_i0p3_i1_i2m3)), MulSIMD(_Rational_16_25, AddSIMD(f_i0p1_i1_i2p1, SubSIMD(SubSIMD(f_i0m1_i1_i2m1, f_i0m1_i1_i2p1), f_i0p1_i1_i2m1))))))))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD11(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY f_i0_i1m4_i2,const REAL_SIMD_ARRAY f_i0_i1m3_i2,const REAL_SIMD_ARRAY f_i0_i1m2_i2,const REAL_SIMD_ARRAY f_i0_i1m1_i2,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1p1_i2,const REAL_SIMD_ARRAY f_i0_i1p2_i2,const REAL_SIMD_ARRAY f_i0_i1p3_i2,const REAL_SIMD_ARRAY f_i0_i1p4_i2) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_5 = 1.0/5.0;
      const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

      const double tmp_Rational_1_560 = 1.0/560.0;
      const REAL_SIMD_ARRAY _Rational_1_560 = ConstSIMD(tmp_Rational_1_560);

      const double tmp_Rational_205_72 = 205.0/72.0;
      const REAL_SIMD_ARRAY _Rational_205_72 = ConstSIMD(tmp_Rational_205_72);

      const double tmp_Rational_8_315 = 8.0/315.0;
      const REAL_SIMD_ARRAY _Rational_8_315 = ConstSIMD(tmp_Rational_8_315);

      const double tmp_Rational_8_5 = 8.0/5.0;
      const REAL_SIMD_ARRAY _Rational_8_5 = ConstSIMD(tmp_Rational_8_5);

      return MulSIMD(MulSIMD(invdx1, invdx1), FusedMulAddSIMD(_Rational_1_560, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1p4_i2, f_i0_i1m4_i2)), FusedMulAddSIMD(_Rational_8_315, AddSIMD(f_i0_i1m3_i2, f_i0_i1p3_i2), FusedMulAddSIMD(_Rational_8_5, AddSIMD(f_i0_i1m1_i2, f_i0_i1p1_i2), NegFusedMulSubSIMD(_Rational_1_5, AddSIMD(f_i0_i1p2_i2, f_i0_i1m2_i2), MulSIMD(_Rational_205_72, f))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD12(const REAL_SIMD_ARRAY invdx1,const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1m4_i2m4,const REAL_SIMD_ARRAY f_i0_i1m3_i2m4,const REAL_SIMD_ARRAY f_i0_i1m2_i2m4,const REAL_SIMD_ARRAY f_i0_i1m1_i2m4,const REAL_SIMD_ARRAY f_i0_i1p1_i2m4,const REAL_SIMD_ARRAY f_i0_i1p2_i2m4,const REAL_SIMD_ARRAY f_i0_i1p3_i2m4,const REAL_SIMD_ARRAY f_i0_i1p4_i2m4,const REAL_SIMD_ARRAY f_i0_i1m4_i2m3,const REAL_SIMD_ARRAY f_i0_i1m3_i2m3,const REAL_SIMD_ARRAY f_i0_i1m2_i2m3,const REAL_SIMD_ARRAY f_i0_i1m1_i2m3,const REAL_SIMD_ARRAY f_i0_i1p1_i2m3,const REAL_SIMD_ARRAY f_i0_i1p2_i2m3,const REAL_SIMD_ARRAY f_i0_i1p3_i2m3,const REAL_SIMD_ARRAY f_i0_i1p4_i2m3,const REAL_SIMD_ARRAY f_i0_i1m4_i2m2,const REAL_SIMD_ARRAY f_i0_i1m3_i2m2,const REAL_SIMD_ARRAY f_i0_i1m2_i2m2,const REAL_SIMD_ARRAY f_i0_i1m1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p1_i2m2,const REAL_SIMD_ARRAY f_i0_i1p2_i2m2,const REAL_SIMD_ARRAY f_i0_i1p3_i2m2,const REAL_SIMD_ARRAY f_i0_i1p4_i2m2,const REAL_SIMD_ARRAY f_i0_i1m4_i2m1,const REAL_SIMD_ARRAY f_i0_i1m3_i2m1,const REAL_SIMD_ARRAY f_i0_i1m2_i2m1,const REAL_SIMD_ARRAY f_i0_i1m1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p1_i2m1,const REAL_SIMD_ARRAY f_i0_i1p2_i2m1,const REAL_SIMD_ARRAY f_i0_i1p3_i2m1,const REAL_SIMD_ARRAY f_i0_i1p4_i2m1,const REAL_SIMD_ARRAY f_i0_i1m4_i2p1,const REAL_SIMD_ARRAY f_i0_i1m3_i2p1,const REAL_SIMD_ARRAY f_i0_i1m2_i2p1,const REAL_SIMD_ARRAY f_i0_i1m1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p1_i2p1,const REAL_SIMD_ARRAY f_i0_i1p2_i2p1,const REAL_SIMD_ARRAY f_i0_i1p3_i2p1,const REAL_SIMD_ARRAY f_i0_i1p4_i2p1,const REAL_SIMD_ARRAY f_i0_i1m4_i2p2,const REAL_SIMD_ARRAY f_i0_i1m3_i2p2,const REAL_SIMD_ARRAY f_i0_i1m2_i2p2,const REAL_SIMD_ARRAY f_i0_i1m1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p1_i2p2,const REAL_SIMD_ARRAY f_i0_i1p2_i2p2,const REAL_SIMD_ARRAY f_i0_i1p3_i2p2,const REAL_SIMD_ARRAY f_i0_i1p4_i2p2,const REAL_SIMD_ARRAY f_i0_i1m4_i2p3,const REAL_SIMD_ARRAY f_i0_i1m3_i2p3,const REAL_SIMD_ARRAY f_i0_i1m2_i2p3,const REAL_SIMD_ARRAY f_i0_i1m1_i2p3,const REAL_SIMD_ARRAY f_i0_i1p1_i2p3,const REAL_SIMD_ARRAY f_i0_i1p2_i2p3,const REAL_SIMD_ARRAY f_i0_i1p3_i2p3,const REAL_SIMD_ARRAY f_i0_i1p4_i2p3,const REAL_SIMD_ARRAY f_i0_i1m4_i2p4,const REAL_SIMD_ARRAY f_i0_i1m3_i2p4,const REAL_SIMD_ARRAY f_i0_i1m2_i2p4,const REAL_SIMD_ARRAY f_i0_i1m1_i2p4,const REAL_SIMD_ARRAY f_i0_i1p1_i2p4,const REAL_SIMD_ARRAY f_i0_i1p2_i2p4,const REAL_SIMD_ARRAY f_i0_i1p3_i2p4,const REAL_SIMD_ARRAY f_i0_i1p4_i2p4) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_16_11025 = 16.0/11025.0;
      const REAL_SIMD_ARRAY _Rational_16_11025 = ConstSIMD(tmp_Rational_16_11025);

      const double tmp_Rational_16_25 = 16.0/25.0;
      const REAL_SIMD_ARRAY _Rational_16_25 = ConstSIMD(tmp_Rational_16_25);

      const double tmp_Rational_16_525 = 16.0/525.0;
      const REAL_SIMD_ARRAY _Rational_16_525 = ConstSIMD(tmp_Rational_16_525);

      const double tmp_Rational_1_1400 = 1.0/1400.0;
      const REAL_SIMD_ARRAY _Rational_1_1400 = ConstSIMD(tmp_Rational_1_1400);

      const double tmp_Rational_1_25 = 1.0/25.0;
      const REAL_SIMD_ARRAY _Rational_1_25 = ConstSIMD(tmp_Rational_1_25);

      const double tmp_Rational_1_350 = 1.0/350.0;
      const REAL_SIMD_ARRAY _Rational_1_350 = ConstSIMD(tmp_Rational_1_350);

      const double tmp_Rational_1_7350 = 1.0/7350.0;
      const REAL_SIMD_ARRAY _Rational_1_7350 = ConstSIMD(tmp_Rational_1_7350);

      const double tmp_Rational_1_78400 = 1.0/78400.0;
      const REAL_SIMD_ARRAY _Rational_1_78400 = ConstSIMD(tmp_Rational_1_78400);

      const double tmp_Rational_4_25 = 4.0/25.0;
      const REAL_SIMD_ARRAY _Rational_4_25 = ConstSIMD(tmp_Rational_4_25);

      const double tmp_Rational_4_525 = 4.0/525.0;
      const REAL_SIMD_ARRAY _Rational_4_525 = ConstSIMD(tmp_Rational_4_525);

      return MulSIMD(invdx1, MulSIMD(invdx2, FusedMulAddSIMD(_Rational_1_7350, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p4_i2m3, f_i0_i1m4_i2m3), AddSIMD(AddSIMD(f_i0_i1m4_i2p3, f_i0_i1p3_i2m4), SubSIMD(SubSIMD(f_i0_i1m3_i2p4, f_i0_i1m3_i2m4), f_i0_i1p4_i2p3))), f_i0_i1p3_i2p4), FusedMulAddSIMD(_Rational_1_78400, AddSIMD(f_i0_i1p4_i2p4, SubSIMD(SubSIMD(f_i0_i1m4_i2m4, f_i0_i1m4_i2p4), f_i0_i1p4_i2m4)), FusedMulAddSIMD(_Rational_1_25, AddSIMD(f_i0_i1p2_i2p2, SubSIMD(SubSIMD(f_i0_i1m2_i2m2, f_i0_i1m2_i2p2), f_i0_i1p2_i2m2)), FusedMulAddSIMD(_Rational_1_350, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p4_i2m1, f_i0_i1m4_i2m1), AddSIMD(AddSIMD(f_i0_i1m4_i2p1, f_i0_i1p1_i2m4), SubSIMD(SubSIMD(f_i0_i1m1_i2p4, f_i0_i1m1_i2m4), f_i0_i1p4_i2p1))), f_i0_i1p1_i2p4), FusedMulAddSIMD(_Rational_16_525, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p3_i2p1, f_i0_i1m3_i2p1), AddSIMD(AddSIMD(f_i0_i1m3_i2m1, f_i0_i1p1_i2p3), SubSIMD(SubSIMD(f_i0_i1m1_i2m3, f_i0_i1m1_i2p3), f_i0_i1p3_i2m1))), f_i0_i1p1_i2m3), FusedMulAddSIMD(_Rational_1_1400, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p4_i2p2, f_i0_i1m4_i2p2), AddSIMD(AddSIMD(f_i0_i1m4_i2m2, f_i0_i1p2_i2p4), SubSIMD(SubSIMD(f_i0_i1m2_i2m4, f_i0_i1m2_i2p4), f_i0_i1p4_i2m2))), f_i0_i1p2_i2m4), FusedMulAddSIMD(_Rational_4_25, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p2_i2m1, f_i0_i1m2_i2m1), AddSIMD(AddSIMD(f_i0_i1m2_i2p1, f_i0_i1p1_i2m2), SubSIMD(SubSIMD(f_i0_i1m1_i2p2, f_i0_i1m1_i2m2), f_i0_i1p2_i2p1))), f_i0_i1p1_i2p2), FusedMulAddSIMD(_Rational_4_525, SubSIMD(AddSIMD(SubSIMD(f_i0_i1p3_i2m2, f_i0_i1m3_i2m2), AddSIMD(AddSIMD(f_i0_i1m3_i2p2, f_i0_i1p2_i2m3), SubSIMD(SubSIMD(f_i0_i1m2_i2p3, f_i0_i1m2_i2m3), f_i0_i1p3_i2p2))), f_i0_i1p2_i2p3), FusedMulAddSIMD(_Rational_16_11025, AddSIMD(f_i0_i1p3_i2p3, SubSIMD(SubSIMD(f_i0_i1m3_i2m3, f_i0_i1m3_i2p3), f_i0_i1p3_i2m3)), MulSIMD(_Rational_16_25, AddSIMD(f_i0_i1p1_i2p1, SubSIMD(SubSIMD(f_i0_i1m1_i2m1, f_i0_i1m1_i2p1), f_i0_i1p1_i2m1))))))))))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22
 */
static REAL_SIMD_ARRAY _NOINLINE _UNUSED SIMD_order_8_f_dDD22(const REAL_SIMD_ARRAY invdx2,const REAL_SIMD_ARRAY f_i0_i1_i2m4,const REAL_SIMD_ARRAY f_i0_i1_i2m3,const REAL_SIMD_ARRAY f_i0_i1_i2m2,const REAL_SIMD_ARRAY f_i0_i1_i2m1,const REAL_SIMD_ARRAY f,const REAL_SIMD_ARRAY f_i0_i1_i2p1,const REAL_SIMD_ARRAY f_i0_i1_i2p2,const REAL_SIMD_ARRAY f_i0_i1_i2p3,const REAL_SIMD_ARRAY f_i0_i1_i2p4) {

      const double tmp_NegativeOne_ = -1.0;
      const REAL_SIMD_ARRAY __attribute__((unused)) _NegativeOne_ = ConstSIMD(tmp_NegativeOne_);

      const double tmp_Rational_1_5 = 1.0/5.0;
      const REAL_SIMD_ARRAY _Rational_1_5 = ConstSIMD(tmp_Rational_1_5);

      const double tmp_Rational_1_560 = 1.0/560.0;
      const REAL_SIMD_ARRAY _Rational_1_560 = ConstSIMD(tmp_Rational_1_560);

      const double tmp_Rational_205_72 = 205.0/72.0;
      const REAL_SIMD_ARRAY _Rational_205_72 = ConstSIMD(tmp_Rational_205_72);

      const double tmp_Rational_8_315 = 8.0/315.0;
      const REAL_SIMD_ARRAY _Rational_8_315 = ConstSIMD(tmp_Rational_8_315);

      const double tmp_Rational_8_5 = 8.0/5.0;
      const REAL_SIMD_ARRAY _Rational_8_5 = ConstSIMD(tmp_Rational_8_5);

      return MulSIMD(MulSIMD(invdx2, invdx2), FusedMulAddSIMD(_Rational_1_560, MulSIMD(_NegativeOne_, AddSIMD(f_i0_i1_i2p4, f_i0_i1_i2m4)), FusedMulAddSIMD(_Rational_8_315, AddSIMD(f_i0_i1_i2m3, f_i0_i1_i2p3), FusedMulAddSIMD(_Rational_8_5, AddSIMD(f_i0_i1_i2m1, f_i0_i1_i2p1), NegFusedMulSubSIMD(_Rational_1_5, AddSIMD(f_i0_i1_i2p2, f_i0_i1_i2m2), MulSIMD(_Rational_205_72, f))))));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 0
 */
static REAL _NOINLINE _UNUSED order_8_f_dD0(const REAL invdx0,const REAL f_i0m4_i1_i2,const REAL f_i0m3_i1_i2,const REAL f_i0m2_i1_i2,const REAL f_i0m1_i1_i2,const REAL f_i0p1_i1_i2,const REAL f_i0p2_i1_i2,const REAL f_i0p3_i1_i2,const REAL f_i0p4_i1_i2) {

      const double _Rational_4_5 = 4.0/5.0;
      const double _Rational_4_105 = 4.0/105.0;
      const double _Rational_1_5 = 1.0/5.0;
      const double _Rational_1_280 = 1.0/280.0;
      return invdx0*(_Rational_1_280*(f_i0m4_i1_i2 - f_i0p4_i1_i2) + _Rational_1_5*(f_i0m2_i1_i2 - f_i0p2_i1_i2) + _Rational_4_105*(-f_i0m3_i1_i2 + f_i0p3_i1_i2) + _Rational_4_5*(-f_i0m1_i1_i2 + f_i0p1_i1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 1
 */
static REAL _NOINLINE _UNUSED order_8_f_dD1(const REAL invdx1,const REAL f_i0_i1m4_i2,const REAL f_i0_i1m3_i2,const REAL f_i0_i1m2_i2,const REAL f_i0_i1m1_i2,const REAL f_i0_i1p1_i2,const REAL f_i0_i1p2_i2,const REAL f_i0_i1p3_i2,const REAL f_i0_i1p4_i2) {

      const double _Rational_4_5 = 4.0/5.0;
      const double _Rational_4_105 = 4.0/105.0;
      const double _Rational_1_5 = 1.0/5.0;
      const double _Rational_1_280 = 1.0/280.0;
      return invdx1*(_Rational_1_280*(f_i0_i1m4_i2 - f_i0_i1p4_i2) + _Rational_1_5*(f_i0_i1m2_i2 - f_i0_i1p2_i2) + _Rational_4_105*(-f_i0_i1m3_i2 + f_i0_i1p3_i2) + _Rational_4_5*(-f_i0_i1m1_i2 + f_i0_i1p1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for first derivative: 2
 */
static REAL _NOINLINE _UNUSED order_8_f_dD2(const REAL invdx2,const REAL f_i0_i1_i2m4,const REAL f_i0_i1_i2m3,const REAL f_i0_i1_i2m2,const REAL f_i0_i1_i2m1,const REAL f_i0_i1_i2p1,const REAL f_i0_i1_i2p2,const REAL f_i0_i1_i2p3,const REAL f_i0_i1_i2p4) {

      const double _Rational_4_5 = 4.0/5.0;
      const double _Rational_4_105 = 4.0/105.0;
      const double _Rational_1_5 = 1.0/5.0;
      const double _Rational_1_280 = 1.0/280.0;
      return invdx2*(_Rational_1_280*(f_i0_i1_i2m4 - f_i0_i1_i2p4) + _Rational_1_5*(f_i0_i1_i2m2 - f_i0_i1_i2p2) + _Rational_4_105*(-f_i0_i1_i2m3 + f_i0_i1_i2p3) + _Rational_4_5*(-f_i0_i1_i2m1 + f_i0_i1_i2p1));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 00
 */
static REAL _NOINLINE _UNUSED order_8_f_dDD00(const REAL invdx0,const REAL f_i0m4_i1_i2,const REAL f_i0m3_i1_i2,const REAL f_i0m2_i1_i2,const REAL f_i0m1_i1_i2,const REAL f,const REAL f_i0p1_i1_i2,const REAL f_i0p2_i1_i2,const REAL f_i0p3_i1_i2,const REAL f_i0p4_i1_i2) {

      const double _Rational_205_72 = 205.0/72.0;
      const double _Rational_1_5 = 1.0/5.0;
      const double _Rational_1_560 = 1.0/560.0;
      const double _Rational_8_5 = 8.0/5.0;
      const double _Rational_8_315 = 8.0/315.0;
      return ((invdx0)*(invdx0))*(_Rational_1_5*(-f_i0m2_i1_i2 - f_i0p2_i1_i2) + _Rational_1_560*(-f_i0m4_i1_i2 - f_i0p4_i1_i2) - _Rational_205_72*f + _Rational_8_315*(f_i0m3_i1_i2 + f_i0p3_i1_i2) + _Rational_8_5*(f_i0m1_i1_i2 + f_i0p1_i1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 01
 */
static REAL _NOINLINE _UNUSED order_8_f_dDD01(const REAL invdx0,const REAL invdx1,const REAL f_i0m4_i1m4_i2,const REAL f_i0m3_i1m4_i2,const REAL f_i0m2_i1m4_i2,const REAL f_i0m1_i1m4_i2,const REAL f_i0p1_i1m4_i2,const REAL f_i0p2_i1m4_i2,const REAL f_i0p3_i1m4_i2,const REAL f_i0p4_i1m4_i2,const REAL f_i0m4_i1m3_i2,const REAL f_i0m3_i1m3_i2,const REAL f_i0m2_i1m3_i2,const REAL f_i0m1_i1m3_i2,const REAL f_i0p1_i1m3_i2,const REAL f_i0p2_i1m3_i2,const REAL f_i0p3_i1m3_i2,const REAL f_i0p4_i1m3_i2,const REAL f_i0m4_i1m2_i2,const REAL f_i0m3_i1m2_i2,const REAL f_i0m2_i1m2_i2,const REAL f_i0m1_i1m2_i2,const REAL f_i0p1_i1m2_i2,const REAL f_i0p2_i1m2_i2,const REAL f_i0p3_i1m2_i2,const REAL f_i0p4_i1m2_i2,const REAL f_i0m4_i1m1_i2,const REAL f_i0m3_i1m1_i2,const REAL f_i0m2_i1m1_i2,const REAL f_i0m1_i1m1_i2,const REAL f_i0p1_i1m1_i2,const REAL f_i0p2_i1m1_i2,const REAL f_i0p3_i1m1_i2,const REAL f_i0p4_i1m1_i2,const REAL f_i0m4_i1p1_i2,const REAL f_i0m3_i1p1_i2,const REAL f_i0m2_i1p1_i2,const REAL f_i0m1_i1p1_i2,const REAL f_i0p1_i1p1_i2,const REAL f_i0p2_i1p1_i2,const REAL f_i0p3_i1p1_i2,const REAL f_i0p4_i1p1_i2,const REAL f_i0m4_i1p2_i2,const REAL f_i0m3_i1p2_i2,const REAL f_i0m2_i1p2_i2,const REAL f_i0m1_i1p2_i2,const REAL f_i0p1_i1p2_i2,const REAL f_i0p2_i1p2_i2,const REAL f_i0p3_i1p2_i2,const REAL f_i0p4_i1p2_i2,const REAL f_i0m4_i1p3_i2,const REAL f_i0m3_i1p3_i2,const REAL f_i0m2_i1p3_i2,const REAL f_i0m1_i1p3_i2,const REAL f_i0p1_i1p3_i2,const REAL f_i0p2_i1p3_i2,const REAL f_i0p3_i1p3_i2,const REAL f_i0p4_i1p3_i2,const REAL f_i0m4_i1p4_i2,const REAL f_i0m3_i1p4_i2,const REAL f_i0m2_i1p4_i2,const REAL f_i0m1_i1p4_i2,const REAL f_i0p1_i1p4_i2,const REAL f_i0p2_i1p4_i2,const REAL f_i0p3_i1p4_i2,const REAL f_i0p4_i1p4_i2) {

      const double _Rational_16_25 = 16.0/25.0;
      const double _Rational_16_525 = 16.0/525.0;
      const double _Rational_16_11025 = 16.0/11025.0;
      const double _Rational_4_25 = 4.0/25.0;
      const double _Rational_4_525 = 4.0/525.0;
      const double _Rational_1_25 = 1.0/25.0;
      const double _Rational_1_350 = 1.0/350.0;
      const double _Rational_1_1400 = 1.0/1400.0;
      const double _Rational_1_7350 = 1.0/7350.0;
      const double _Rational_1_78400 = 1.0/78400.0;
      return invdx0*invdx1*(_Rational_16_11025*(f_i0m3_i1m3_i2 - f_i0m3_i1p3_i2 - f_i0p3_i1m3_i2 + f_i0p3_i1p3_i2) + _Rational_16_25*(f_i0m1_i1m1_i2 - f_i0m1_i1p1_i2 - f_i0p1_i1m1_i2 + f_i0p1_i1p1_i2) + _Rational_16_525*(f_i0m1_i1m3_i2 - f_i0m1_i1p3_i2 + f_i0m3_i1m1_i2 - f_i0m3_i1p1_i2 - f_i0p1_i1m3_i2 + f_i0p1_i1p3_i2 - f_i0p3_i1m1_i2 + f_i0p3_i1p1_i2) + _Rational_1_1400*(f_i0m2_i1m4_i2 - f_i0m2_i1p4_i2 + f_i0m4_i1m2_i2 - f_i0m4_i1p2_i2 - f_i0p2_i1m4_i2 + f_i0p2_i1p4_i2 - f_i0p4_i1m2_i2 + f_i0p4_i1p2_i2) + _Rational_1_25*(f_i0m2_i1m2_i2 - f_i0m2_i1p2_i2 - f_i0p2_i1m2_i2 + f_i0p2_i1p2_i2) + _Rational_1_350*(-f_i0m1_i1m4_i2 + f_i0m1_i1p4_i2 - f_i0m4_i1m1_i2 + f_i0m4_i1p1_i2 + f_i0p1_i1m4_i2 - f_i0p1_i1p4_i2 + f_i0p4_i1m1_i2 - f_i0p4_i1p1_i2) + _Rational_1_7350*(-f_i0m3_i1m4_i2 + f_i0m3_i1p4_i2 - f_i0m4_i1m3_i2 + f_i0m4_i1p3_i2 + f_i0p3_i1m4_i2 - f_i0p3_i1p4_i2 + f_i0p4_i1m3_i2 - f_i0p4_i1p3_i2) + _Rational_1_78400*(f_i0m4_i1m4_i2 - f_i0m4_i1p4_i2 - f_i0p4_i1m4_i2 + f_i0p4_i1p4_i2) + _Rational_4_25*(-f_i0m1_i1m2_i2 + f_i0m1_i1p2_i2 - f_i0m2_i1m1_i2 + f_i0m2_i1p1_i2 + f_i0p1_i1m2_i2 - f_i0p1_i1p2_i2 + f_i0p2_i1m1_i2 - f_i0p2_i1p1_i2) + _Rational_4_525*(-f_i0m2_i1m3_i2 + f_i0m2_i1p3_i2 - f_i0m3_i1m2_i2 + f_i0m3_i1p2_i2 + f_i0p2_i1m3_i2 - f_i0p2_i1p3_i2 + f_i0p3_i1m2_i2 - f_i0p3_i1p2_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 02
 */
static REAL _NOINLINE _UNUSED order_8_f_dDD02(const REAL invdx0,const REAL invdx2,const REAL f_i0m4_i1_i2m4,const REAL f_i0m3_i1_i2m4,const REAL f_i0m2_i1_i2m4,const REAL f_i0m1_i1_i2m4,const REAL f_i0p1_i1_i2m4,const REAL f_i0p2_i1_i2m4,const REAL f_i0p3_i1_i2m4,const REAL f_i0p4_i1_i2m4,const REAL f_i0m4_i1_i2m3,const REAL f_i0m3_i1_i2m3,const REAL f_i0m2_i1_i2m3,const REAL f_i0m1_i1_i2m3,const REAL f_i0p1_i1_i2m3,const REAL f_i0p2_i1_i2m3,const REAL f_i0p3_i1_i2m3,const REAL f_i0p4_i1_i2m3,const REAL f_i0m4_i1_i2m2,const REAL f_i0m3_i1_i2m2,const REAL f_i0m2_i1_i2m2,const REAL f_i0m1_i1_i2m2,const REAL f_i0p1_i1_i2m2,const REAL f_i0p2_i1_i2m2,const REAL f_i0p3_i1_i2m2,const REAL f_i0p4_i1_i2m2,const REAL f_i0m4_i1_i2m1,const REAL f_i0m3_i1_i2m1,const REAL f_i0m2_i1_i2m1,const REAL f_i0m1_i1_i2m1,const REAL f_i0p1_i1_i2m1,const REAL f_i0p2_i1_i2m1,const REAL f_i0p3_i1_i2m1,const REAL f_i0p4_i1_i2m1,const REAL f_i0m4_i1_i2p1,const REAL f_i0m3_i1_i2p1,const REAL f_i0m2_i1_i2p1,const REAL f_i0m1_i1_i2p1,const REAL f_i0p1_i1_i2p1,const REAL f_i0p2_i1_i2p1,const REAL f_i0p3_i1_i2p1,const REAL f_i0p4_i1_i2p1,const REAL f_i0m4_i1_i2p2,const REAL f_i0m3_i1_i2p2,const REAL f_i0m2_i1_i2p2,const REAL f_i0m1_i1_i2p2,const REAL f_i0p1_i1_i2p2,const REAL f_i0p2_i1_i2p2,const REAL f_i0p3_i1_i2p2,const REAL f_i0p4_i1_i2p2,const REAL f_i0m4_i1_i2p3,const REAL f_i0m3_i1_i2p3,const REAL f_i0m2_i1_i2p3,const REAL f_i0m1_i1_i2p3,const REAL f_i0p1_i1_i2p3,const REAL f_i0p2_i1_i2p3,const REAL f_i0p3_i1_i2p3,const REAL f_i0p4_i1_i2p3,const REAL f_i0m4_i1_i2p4,const REAL f_i0m3_i1_i2p4,const REAL f_i0m2_i1_i2p4,const REAL f_i0m1_i1_i2p4,const REAL f_i0p1_i1_i2p4,const REAL f_i0p2_i1_i2p4,const REAL f_i0p3_i1_i2p4,const REAL f_i0p4_i1_i2p4) {

      const double _Rational_16_25 = 16.0/25.0;
      const double _Rational_16_525 = 16.0/525.0;
      const double _Rational_16_11025 = 16.0/11025.0;
      const double _Rational_4_25 = 4.0/25.0;
      const double _Rational_4_525 = 4.0/525.0;
      const double _Rational_1_25 = 1.0/25.0;
      const double _Rational_1_350 = 1.0/350.0;
      const double _Rational_1_1400 = 1.0/1400.0;
      const double _Rational_1_7350 = 1.0/7350.0;
      const double _Rational_1_78400 = 1.0/78400.0;
      return invdx0*invdx2*(_Rational_16_11025*(f_i0m3_i1_i2m3 - f_i0m3_i1_i2p3 - f_i0p3_i1_i2m3 + f_i0p3_i1_i2p3) + _Rational_16_25*(f_i0m1_i1_i2m1 - f_i0m1_i1_i2p1 - f_i0p1_i1_i2m1 + f_i0p1_i1_i2p1) + _Rational_16_525*(f_i0m1_i1_i2m3 - f_i0m1_i1_i2p3 + f_i0m3_i1_i2m1 - f_i0m3_i1_i2p1 - f_i0p1_i1_i2m3 + f_i0p1_i1_i2p3 - f_i0p3_i1_i2m1 + f_i0p3_i1_i2p1) + _Rational_1_1400*(f_i0m2_i1_i2m4 - f_i0m2_i1_i2p4 + f_i0m4_i1_i2m2 - f_i0m4_i1_i2p2 - f_i0p2_i1_i2m4 + f_i0p2_i1_i2p4 - f_i0p4_i1_i2m2 + f_i0p4_i1_i2p2) + _Rational_1_25*(f_i0m2_i1_i2m2 - f_i0m2_i1_i2p2 - f_i0p2_i1_i2m2 + f_i0p2_i1_i2p2) + _Rational_1_350*(-f_i0m1_i1_i2m4 + f_i0m1_i1_i2p4 - f_i0m4_i1_i2m1 + f_i0m4_i1_i2p1 + f_i0p1_i1_i2m4 - f_i0p1_i1_i2p4 + f_i0p4_i1_i2m1 - f_i0p4_i1_i2p1) + _Rational_1_7350*(-f_i0m3_i1_i2m4 + f_i0m3_i1_i2p4 - f_i0m4_i1_i2m3 + f_i0m4_i1_i2p3 + f_i0p3_i1_i2m4 - f_i0p3_i1_i2p4 + f_i0p4_i1_i2m3 - f_i0p4_i1_i2p3) + _Rational_1_78400*(f_i0m4_i1_i2m4 - f_i0m4_i1_i2p4 - f_i0p4_i1_i2m4 + f_i0p4_i1_i2p4) + _Rational_4_25*(-f_i0m1_i1_i2m2 + f_i0m1_i1_i2p2 - f_i0m2_i1_i2m1 + f_i0m2_i1_i2p1 + f_i0p1_i1_i2m2 - f_i0p1_i1_i2p2 + f_i0p2_i1_i2m1 - f_i0p2_i1_i2p1) + _Rational_4_525*(-f_i0m2_i1_i2m3 + f_i0m2_i1_i2p3 - f_i0m3_i1_i2m2 + f_i0m3_i1_i2p2 + f_i0p2_i1_i2m3 - f_i0p2_i1_i2p3 + f_i0p3_i1_i2m2 - f_i0p3_i1_i2p2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 11
 */
static REAL _NOINLINE _UNUSED order_8_f_dDD11(const REAL invdx1,const REAL f_i0_i1m4_i2,const REAL f_i0_i1m3_i2,const REAL f_i0_i1m2_i2,const REAL f_i0_i1m1_i2,const REAL f,const REAL f_i0_i1p1_i2,const REAL f_i0_i1p2_i2,const REAL f_i0_i1p3_i2,const REAL f_i0_i1p4_i2) {

      const double _Rational_205_72 = 205.0/72.0;
      const double _Rational_1_5 = 1.0/5.0;
      const double _Rational_1_560 = 1.0/560.0;
      const double _Rational_8_5 = 8.0/5.0;
      const double _Rational_8_315 = 8.0/315.0;
      return ((invdx1)*(invdx1))*(_Rational_1_5*(-f_i0_i1m2_i2 - f_i0_i1p2_i2) + _Rational_1_560*(-f_i0_i1m4_i2 - f_i0_i1p4_i2) - _Rational_205_72*f + _Rational_8_315*(f_i0_i1m3_i2 + f_i0_i1p3_i2) + _Rational_8_5*(f_i0_i1m1_i2 + f_i0_i1p1_i2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 12
 */
static REAL _NOINLINE _UNUSED order_8_f_dDD12(const REAL invdx1,const REAL invdx2,const REAL f_i0_i1m4_i2m4,const REAL f_i0_i1m3_i2m4,const REAL f_i0_i1m2_i2m4,const REAL f_i0_i1m1_i2m4,const REAL f_i0_i1p1_i2m4,const REAL f_i0_i1p2_i2m4,const REAL f_i0_i1p3_i2m4,const REAL f_i0_i1p4_i2m4,const REAL f_i0_i1m4_i2m3,const REAL f_i0_i1m3_i2m3,const REAL f_i0_i1m2_i2m3,const REAL f_i0_i1m1_i2m3,const REAL f_i0_i1p1_i2m3,const REAL f_i0_i1p2_i2m3,const REAL f_i0_i1p3_i2m3,const REAL f_i0_i1p4_i2m3,const REAL f_i0_i1m4_i2m2,const REAL f_i0_i1m3_i2m2,const REAL f_i0_i1m2_i2m2,const REAL f_i0_i1m1_i2m2,const REAL f_i0_i1p1_i2m2,const REAL f_i0_i1p2_i2m2,const REAL f_i0_i1p3_i2m2,const REAL f_i0_i1p4_i2m2,const REAL f_i0_i1m4_i2m1,const REAL f_i0_i1m3_i2m1,const REAL f_i0_i1m2_i2m1,const REAL f_i0_i1m1_i2m1,const REAL f_i0_i1p1_i2m1,const REAL f_i0_i1p2_i2m1,const REAL f_i0_i1p3_i2m1,const REAL f_i0_i1p4_i2m1,const REAL f_i0_i1m4_i2p1,const REAL f_i0_i1m3_i2p1,const REAL f_i0_i1m2_i2p1,const REAL f_i0_i1m1_i2p1,const REAL f_i0_i1p1_i2p1,const REAL f_i0_i1p2_i2p1,const REAL f_i0_i1p3_i2p1,const REAL f_i0_i1p4_i2p1,const REAL f_i0_i1m4_i2p2,const REAL f_i0_i1m3_i2p2,const REAL f_i0_i1m2_i2p2,const REAL f_i0_i1m1_i2p2,const REAL f_i0_i1p1_i2p2,const REAL f_i0_i1p2_i2p2,const REAL f_i0_i1p3_i2p2,const REAL f_i0_i1p4_i2p2,const REAL f_i0_i1m4_i2p3,const REAL f_i0_i1m3_i2p3,const REAL f_i0_i1m2_i2p3,const REAL f_i0_i1m1_i2p3,const REAL f_i0_i1p1_i2p3,const REAL f_i0_i1p2_i2p3,const REAL f_i0_i1p3_i2p3,const REAL f_i0_i1p4_i2p3,const REAL f_i0_i1m4_i2p4,const REAL f_i0_i1m3_i2p4,const REAL f_i0_i1m2_i2p4,const REAL f_i0_i1m1_i2p4,const REAL f_i0_i1p1_i2p4,const REAL f_i0_i1p2_i2p4,const REAL f_i0_i1p3_i2p4,const REAL f_i0_i1p4_i2p4) {

      const double _Rational_16_25 = 16.0/25.0;
      const double _Rational_16_525 = 16.0/525.0;
      const double _Rational_16_11025 = 16.0/11025.0;
      const double _Rational_4_25 = 4.0/25.0;
      const double _Rational_4_525 = 4.0/525.0;
      const double _Rational_1_25 = 1.0/25.0;
      const double _Rational_1_350 = 1.0/350.0;
      const double _Rational_1_1400 = 1.0/1400.0;
      const double _Rational_1_7350 = 1.0/7350.0;
      const double _Rational_1_78400 = 1.0/78400.0;
      return invdx1*invdx2*(_Rational_16_11025*(f_i0_i1m3_i2m3 - f_i0_i1m3_i2p3 - f_i0_i1p3_i2m3 + f_i0_i1p3_i2p3) + _Rational_16_25*(f_i0_i1m1_i2m1 - f_i0_i1m1_i2p1 - f_i0_i1p1_i2m1 + f_i0_i1p1_i2p1) + _Rational_16_525*(f_i0_i1m1_i2m3 - f_i0_i1m1_i2p3 + f_i0_i1m3_i2m1 - f_i0_i1m3_i2p1 - f_i0_i1p1_i2m3 + f_i0_i1p1_i2p3 - f_i0_i1p3_i2m1 + f_i0_i1p3_i2p1) + _Rational_1_1400*(f_i0_i1m2_i2m4 - f_i0_i1m2_i2p4 + f_i0_i1m4_i2m2 - f_i0_i1m4_i2p2 - f_i0_i1p2_i2m4 + f_i0_i1p2_i2p4 - f_i0_i1p4_i2m2 + f_i0_i1p4_i2p2) + _Rational_1_25*(f_i0_i1m2_i2m2 - f_i0_i1m2_i2p2 - f_i0_i1p2_i2m2 + f_i0_i1p2_i2p2) + _Rational_1_350*(-f_i0_i1m1_i2m4 + f_i0_i1m1_i2p4 - f_i0_i1m4_i2m1 + f_i0_i1m4_i2p1 + f_i0_i1p1_i2m4 - f_i0_i1p1_i2p4 + f_i0_i1p4_i2m1 - f_i0_i1p4_i2p1) + _Rational_1_7350*(-f_i0_i1m3_i2m4 + f_i0_i1m3_i2p4 - f_i0_i1m4_i2m3 + f_i0_i1m4_i2p3 + f_i0_i1p3_i2m4 - f_i0_i1p3_i2p4 + f_i0_i1p4_i2m3 - f_i0_i1p4_i2p3) + _Rational_1_78400*(f_i0_i1m4_i2m4 - f_i0_i1m4_i2p4 - f_i0_i1p4_i2m4 + f_i0_i1p4_i2p4) + _Rational_4_25*(-f_i0_i1m1_i2m2 + f_i0_i1m1_i2p2 - f_i0_i1m2_i2m1 + f_i0_i1m2_i2p1 + f_i0_i1p1_i2m2 - f_i0_i1p1_i2p2 + f_i0_i1p2_i2m1 - f_i0_i1p2_i2p1) + _Rational_4_525*(-f_i0_i1m2_i2m3 + f_i0_i1m2_i2p3 - f_i0_i1m3_i2m2 + f_i0_i1m3_i2p2 + f_i0_i1p2_i2m3 - f_i0_i1p2_i2p3 + f_i0_i1p3_i2m2 - f_i0_i1p3_i2p2));
}
/*
 *  * (__FD_OPERATOR_FUNC__) Finite difference operator for second derivative: 22
 */
static REAL _NOINLINE _UNUSED order_8_f_dDD22(const REAL invdx2,const REAL f_i0_i1_i2m4,const REAL f_i0_i1_i2m3,const REAL f_i0_i1_i2m2,const REAL f_i0_i1_i2m1,const REAL f,const REAL f_i0_i1_i2p1,const REAL f_i0_i1_i2p2,const REAL f_i0_i1_i2p3,const REAL f_i0_i1_i2p4) {

      const double _Rational_205_72 = 205.0/72.0;
      const double _Rational_1_5 = 1.0/5.0;
      const double _Rational_1_560 = 1.0/560.0;
      const double _Rational_8_5 = 8.0/5.0;
      const double _Rational_8_315 = 8.0/315.0;
      return ((invdx2)*(invdx2))*(_Rational_1_5*(-f_i0_i1_i2m2 - f_i0_i1_i2p2) + _Rational_1_560*(-f_i0_i1_i2m4 - f_i0_i1_i2p4) - _Rational_205_72*f + _Rational_8_315*(f_i0_i1_i2m3 + f_i0_i1_i2p3) + _Rational_8_5*(f_i0_i1_i2m1 + f_i0_i1_i2p1));
}
#endif // #ifndef __FD_FUNCTIONS_H__
