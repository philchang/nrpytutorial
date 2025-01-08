#include "../BHaH_defines.h"
#include "../intrinsics/simd_intrinsics.h"
/**
 * Finite difference function for operator dD0, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dD0_fdorder4(const REAL_SIMD_ARRAY FDPROTO_i0m1, const REAL_SIMD_ARRAY FDPROTO_i0m2,
                                                               const REAL_SIMD_ARRAY FDPROTO_i0p1, const REAL_SIMD_ARRAY FDPROTO_i0p2,
                                                               const REAL_SIMD_ARRAY invdxx0) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_2_3 = 2.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_2_3 = ConstSIMD(dblFDPart1_Rational_2_3);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(invdxx0, FusedMulAddSIMD(FDPart1_Rational_1_12, SubSIMD(FDPROTO_i0m2, FDPROTO_i0p2),
                                                                     MulSIMD(FDPart1_Rational_2_3, SubSIMD(FDPROTO_i0p1, FDPROTO_i0m1))));

  return FD_result;
}
/**
 * Finite difference function for operator dD1, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dD1_fdorder4(const REAL_SIMD_ARRAY FDPROTO_i1m1, const REAL_SIMD_ARRAY FDPROTO_i1m2,
                                                               const REAL_SIMD_ARRAY FDPROTO_i1p1, const REAL_SIMD_ARRAY FDPROTO_i1p2,
                                                               const REAL_SIMD_ARRAY invdxx1) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_2_3 = 2.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_2_3 = ConstSIMD(dblFDPart1_Rational_2_3);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(invdxx1, FusedMulAddSIMD(FDPart1_Rational_1_12, SubSIMD(FDPROTO_i1m2, FDPROTO_i1p2),
                                                                     MulSIMD(FDPart1_Rational_2_3, SubSIMD(FDPROTO_i1p1, FDPROTO_i1m1))));

  return FD_result;
}
/**
 * Finite difference function for operator dD2, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dD2_fdorder4(const REAL_SIMD_ARRAY FDPROTO_i2m1, const REAL_SIMD_ARRAY FDPROTO_i2m2,
                                                               const REAL_SIMD_ARRAY FDPROTO_i2p1, const REAL_SIMD_ARRAY FDPROTO_i2p2,
                                                               const REAL_SIMD_ARRAY invdxx2) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_2_3 = 2.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_2_3 = ConstSIMD(dblFDPart1_Rational_2_3);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(invdxx2, FusedMulAddSIMD(FDPart1_Rational_1_12, SubSIMD(FDPROTO_i2m2, FDPROTO_i2p2),
                                                                     MulSIMD(FDPart1_Rational_2_3, SubSIMD(FDPROTO_i2p1, FDPROTO_i2m1))));

  return FD_result;
}
/**
 * Finite difference function for operator dDD00, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dDD00_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i0m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0m2, const REAL_SIMD_ARRAY FDPROTO_i0p1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p2, const REAL_SIMD_ARRAY invdxx0) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_3 = ConstSIMD(dblFDPart1_Rational_4_3);

  static const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_2 = ConstSIMD(dblFDPart1_Rational_5_2);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(MulSIMD(invdxx0, invdxx0),
              FusedMulSubSIMD(FDPart1_Rational_4_3, AddSIMD(FDPROTO_i0m1, FDPROTO_i0p1),
                              FusedMulAddSIMD(FDPROTO, FDPart1_Rational_5_2, MulSIMD(FDPart1_Rational_1_12, AddSIMD(FDPROTO_i0m2, FDPROTO_i0p2)))));

  return FD_result;
}
/**
 * Finite difference function for operator dDD01, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dDD01_fdorder4(const REAL_SIMD_ARRAY FDPROTO_i0m1_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0m1_i1m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0m1_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0m1_i1p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0m2_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0m2_i1m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0m2_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0m2_i1p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p1_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0p1_i1m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p1_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0p1_i1p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p2_i1m1, const REAL_SIMD_ARRAY FDPROTO_i0p2_i1m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p2_i1p1, const REAL_SIMD_ARRAY FDPROTO_i0p2_i1p2,
                                                                 const REAL_SIMD_ARRAY invdxx0, const REAL_SIMD_ARRAY invdxx1) {
  static const double dblFDPart1_Rational_1_144 = 1.0 / 144.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_144 = ConstSIMD(dblFDPart1_Rational_1_144);

  static const double dblFDPart1_Rational_1_18 = 1.0 / 18.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_18 = ConstSIMD(dblFDPart1_Rational_1_18);

  static const double dblFDPart1_Rational_4_9 = 4.0 / 9.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_9 = ConstSIMD(dblFDPart1_Rational_4_9);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx0,
      MulSIMD(invdxx1,
              FusedMulAddSIMD(
                  FDPart1_Rational_1_18,
                  AddSIMD(AddSIMD(FDPROTO_i0m2_i1p1, FDPROTO_i0p1_i1m2),
                          AddSIMD(FDPROTO_i0p2_i1m1, SubSIMD(FDPROTO_i0m1_i1p2, AddSIMD(AddSIMD(FDPROTO_i0m1_i1m2, FDPROTO_i0m2_i1m1),
                                                                                        AddSIMD(FDPROTO_i0p1_i1p2, FDPROTO_i0p2_i1p1))))),
                  FusedMulAddSIMD(FDPart1_Rational_4_9,
                                  AddSIMD(FDPROTO_i0p1_i1p1, SubSIMD(FDPROTO_i0m1_i1m1, AddSIMD(FDPROTO_i0m1_i1p1, FDPROTO_i0p1_i1m1))),
                                  MulSIMD(FDPart1_Rational_1_144,
                                          AddSIMD(FDPROTO_i0p2_i1p2, SubSIMD(FDPROTO_i0m2_i1m2, AddSIMD(FDPROTO_i0m2_i1p2, FDPROTO_i0p2_i1m2))))))));

  return FD_result;
}
/**
 * Finite difference function for operator dDD02, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dDD02_fdorder4(const REAL_SIMD_ARRAY FDPROTO_i0m1_i2m1, const REAL_SIMD_ARRAY FDPROTO_i0m1_i2m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0m1_i2p1, const REAL_SIMD_ARRAY FDPROTO_i0m1_i2p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0m2_i2m1, const REAL_SIMD_ARRAY FDPROTO_i0m2_i2m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0m2_i2p1, const REAL_SIMD_ARRAY FDPROTO_i0m2_i2p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p1_i2m1, const REAL_SIMD_ARRAY FDPROTO_i0p1_i2m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p1_i2p1, const REAL_SIMD_ARRAY FDPROTO_i0p1_i2p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p2_i2m1, const REAL_SIMD_ARRAY FDPROTO_i0p2_i2m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p2_i2p1, const REAL_SIMD_ARRAY FDPROTO_i0p2_i2p2,
                                                                 const REAL_SIMD_ARRAY invdxx0, const REAL_SIMD_ARRAY invdxx2) {
  static const double dblFDPart1_Rational_1_144 = 1.0 / 144.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_144 = ConstSIMD(dblFDPart1_Rational_1_144);

  static const double dblFDPart1_Rational_1_18 = 1.0 / 18.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_18 = ConstSIMD(dblFDPart1_Rational_1_18);

  static const double dblFDPart1_Rational_4_9 = 4.0 / 9.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_9 = ConstSIMD(dblFDPart1_Rational_4_9);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx0,
      MulSIMD(invdxx2,
              FusedMulAddSIMD(
                  FDPart1_Rational_1_18,
                  AddSIMD(AddSIMD(FDPROTO_i0m2_i2p1, FDPROTO_i0p1_i2m2),
                          AddSIMD(FDPROTO_i0p2_i2m1, SubSIMD(FDPROTO_i0m1_i2p2, AddSIMD(AddSIMD(FDPROTO_i0m1_i2m2, FDPROTO_i0m2_i2m1),
                                                                                        AddSIMD(FDPROTO_i0p1_i2p2, FDPROTO_i0p2_i2p1))))),
                  FusedMulAddSIMD(FDPart1_Rational_4_9,
                                  AddSIMD(FDPROTO_i0p1_i2p1, SubSIMD(FDPROTO_i0m1_i2m1, AddSIMD(FDPROTO_i0m1_i2p1, FDPROTO_i0p1_i2m1))),
                                  MulSIMD(FDPart1_Rational_1_144,
                                          AddSIMD(FDPROTO_i0p2_i2p2, SubSIMD(FDPROTO_i0m2_i2m2, AddSIMD(FDPROTO_i0m2_i2p2, FDPROTO_i0p2_i2m2))))))));

  return FD_result;
}
/**
 * Finite difference function for operator dDD11, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dDD11_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i1m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1m2, const REAL_SIMD_ARRAY FDPROTO_i1p1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1p2, const REAL_SIMD_ARRAY invdxx1) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_3 = ConstSIMD(dblFDPart1_Rational_4_3);

  static const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_2 = ConstSIMD(dblFDPart1_Rational_5_2);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(MulSIMD(invdxx1, invdxx1),
              FusedMulSubSIMD(FDPart1_Rational_4_3, AddSIMD(FDPROTO_i1m1, FDPROTO_i1p1),
                              FusedMulAddSIMD(FDPROTO, FDPart1_Rational_5_2, MulSIMD(FDPart1_Rational_1_12, AddSIMD(FDPROTO_i1m2, FDPROTO_i1p2)))));

  return FD_result;
}
/**
 * Finite difference function for operator dDD12, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dDD12_fdorder4(const REAL_SIMD_ARRAY FDPROTO_i1m1_i2m1, const REAL_SIMD_ARRAY FDPROTO_i1m1_i2m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1m1_i2p1, const REAL_SIMD_ARRAY FDPROTO_i1m1_i2p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1m2_i2m1, const REAL_SIMD_ARRAY FDPROTO_i1m2_i2m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1m2_i2p1, const REAL_SIMD_ARRAY FDPROTO_i1m2_i2p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1p1_i2m1, const REAL_SIMD_ARRAY FDPROTO_i1p1_i2m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1p1_i2p1, const REAL_SIMD_ARRAY FDPROTO_i1p1_i2p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1p2_i2m1, const REAL_SIMD_ARRAY FDPROTO_i1p2_i2m2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1p2_i2p1, const REAL_SIMD_ARRAY FDPROTO_i1p2_i2p2,
                                                                 const REAL_SIMD_ARRAY invdxx1, const REAL_SIMD_ARRAY invdxx2) {
  static const double dblFDPart1_Rational_1_144 = 1.0 / 144.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_144 = ConstSIMD(dblFDPart1_Rational_1_144);

  static const double dblFDPart1_Rational_1_18 = 1.0 / 18.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_18 = ConstSIMD(dblFDPart1_Rational_1_18);

  static const double dblFDPart1_Rational_4_9 = 4.0 / 9.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_9 = ConstSIMD(dblFDPart1_Rational_4_9);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx1,
      MulSIMD(invdxx2,
              FusedMulAddSIMD(
                  FDPart1_Rational_1_18,
                  AddSIMD(AddSIMD(FDPROTO_i1m2_i2p1, FDPROTO_i1p1_i2m2),
                          AddSIMD(FDPROTO_i1p2_i2m1, SubSIMD(FDPROTO_i1m1_i2p2, AddSIMD(AddSIMD(FDPROTO_i1m1_i2m2, FDPROTO_i1m2_i2m1),
                                                                                        AddSIMD(FDPROTO_i1p1_i2p2, FDPROTO_i1p2_i2p1))))),
                  FusedMulAddSIMD(FDPart1_Rational_4_9,
                                  AddSIMD(FDPROTO_i1p1_i2p1, SubSIMD(FDPROTO_i1m1_i2m1, AddSIMD(FDPROTO_i1m1_i2p1, FDPROTO_i1p1_i2m1))),
                                  MulSIMD(FDPart1_Rational_1_144,
                                          AddSIMD(FDPROTO_i1p2_i2p2, SubSIMD(FDPROTO_i1m2_i2m2, AddSIMD(FDPROTO_i1m2_i2p2, FDPROTO_i1p2_i2m2))))))));

  return FD_result;
}
/**
 * Finite difference function for operator dDD22, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dDD22_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i2m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i2m2, const REAL_SIMD_ARRAY FDPROTO_i2p1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i2p2, const REAL_SIMD_ARRAY invdxx2) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_3 = ConstSIMD(dblFDPart1_Rational_4_3);

  static const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_2 = ConstSIMD(dblFDPart1_Rational_5_2);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(MulSIMD(invdxx2, invdxx2),
              FusedMulSubSIMD(FDPart1_Rational_4_3, AddSIMD(FDPROTO_i2m1, FDPROTO_i2p1),
                              FusedMulAddSIMD(FDPROTO, FDPart1_Rational_5_2, MulSIMD(FDPart1_Rational_1_12, AddSIMD(FDPROTO_i2m2, FDPROTO_i2p2)))));

  return FD_result;
}
/**
 * Finite difference function for operator ddnD0, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_ddnD0_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i0m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0m2, const REAL_SIMD_ARRAY FDPROTO_i0m3,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p1, const REAL_SIMD_ARRAY invdxx0) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  static const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  static const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  static const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_6 = ConstSIMD(dblFDPart1_Rational_5_6);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(invdxx0, FusedMulAddSIMD(FDPROTO_i0m2, FDPart1_Rational_1_2,
                                       FusedMulAddSIMD(FDPROTO_i0p1, FDPart1_Rational_1_4,
                                                       FusedMulSubSIMD(FDPROTO, FDPart1_Rational_5_6,
                                                                       FusedMulAddSIMD(FDPROTO_i0m1, FDPart1_Rational_3_2,
                                                                                       MulSIMD(FDPROTO_i0m3, FDPart1_Rational_1_12))))));

  return FD_result;
}
/**
 * Finite difference function for operator ddnD1, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_ddnD1_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i1m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1m2, const REAL_SIMD_ARRAY FDPROTO_i1m3,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1p1, const REAL_SIMD_ARRAY invdxx1) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  static const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  static const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  static const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_6 = ConstSIMD(dblFDPart1_Rational_5_6);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(invdxx1, FusedMulAddSIMD(FDPROTO_i1m2, FDPart1_Rational_1_2,
                                       FusedMulAddSIMD(FDPROTO_i1p1, FDPart1_Rational_1_4,
                                                       FusedMulSubSIMD(FDPROTO, FDPart1_Rational_5_6,
                                                                       FusedMulAddSIMD(FDPROTO_i1m1, FDPart1_Rational_3_2,
                                                                                       MulSIMD(FDPROTO_i1m3, FDPart1_Rational_1_12))))));

  return FD_result;
}
/**
 * Finite difference function for operator ddnD2, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_ddnD2_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i2m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i2m2, const REAL_SIMD_ARRAY FDPROTO_i2m3,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i2p1, const REAL_SIMD_ARRAY invdxx2) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  static const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  static const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  static const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_6 = ConstSIMD(dblFDPart1_Rational_5_6);

  const REAL_SIMD_ARRAY FD_result =
      MulSIMD(invdxx2, FusedMulAddSIMD(FDPROTO_i2m2, FDPart1_Rational_1_2,
                                       FusedMulAddSIMD(FDPROTO_i2p1, FDPart1_Rational_1_4,
                                                       FusedMulSubSIMD(FDPROTO, FDPart1_Rational_5_6,
                                                                       FusedMulAddSIMD(FDPROTO_i2m1, FDPart1_Rational_3_2,
                                                                                       MulSIMD(FDPROTO_i2m3, FDPart1_Rational_1_12))))));

  return FD_result;
}
/**
 * Finite difference function for operator dupD0, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dupD0_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i0m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p1, const REAL_SIMD_ARRAY FDPROTO_i0p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i0p3, const REAL_SIMD_ARRAY invdxx0) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  static const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  static const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  static const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_6 = ConstSIMD(dblFDPart1_Rational_5_6);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx0,
      FusedMulAddSIMD(FDPROTO_i0p3, FDPart1_Rational_1_12,
                      FusedMulSubSIMD(FDPROTO_i0p1, FDPart1_Rational_3_2,
                                      FusedMulAddSIMD(FDPROTO_i0m1, FDPart1_Rational_1_4,
                                                      FusedMulAddSIMD(FDPROTO_i0p2, FDPart1_Rational_1_2, MulSIMD(FDPROTO, FDPart1_Rational_5_6))))));

  return FD_result;
}
/**
 * Finite difference function for operator dupD1, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dupD1_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i1m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1p1, const REAL_SIMD_ARRAY FDPROTO_i1p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i1p3, const REAL_SIMD_ARRAY invdxx1) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  static const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  static const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  static const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_6 = ConstSIMD(dblFDPart1_Rational_5_6);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx1,
      FusedMulAddSIMD(FDPROTO_i1p3, FDPart1_Rational_1_12,
                      FusedMulSubSIMD(FDPROTO_i1p1, FDPart1_Rational_3_2,
                                      FusedMulAddSIMD(FDPROTO_i1m1, FDPart1_Rational_1_4,
                                                      FusedMulAddSIMD(FDPROTO_i1p2, FDPart1_Rational_1_2, MulSIMD(FDPROTO, FDPart1_Rational_5_6))))));

  return FD_result;
}
/**
 * Finite difference function for operator dupD2, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dupD2_fdorder4(const REAL_SIMD_ARRAY FDPROTO, const REAL_SIMD_ARRAY FDPROTO_i2m1,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i2p1, const REAL_SIMD_ARRAY FDPROTO_i2p2,
                                                                 const REAL_SIMD_ARRAY FDPROTO_i2p3, const REAL_SIMD_ARRAY invdxx2) {
  static const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  static const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  static const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  static const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  static const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_5_6 = ConstSIMD(dblFDPart1_Rational_5_6);

  const REAL_SIMD_ARRAY FD_result = MulSIMD(
      invdxx2,
      FusedMulAddSIMD(FDPROTO_i2p3, FDPart1_Rational_1_12,
                      FusedMulSubSIMD(FDPROTO_i2p1, FDPart1_Rational_3_2,
                                      FusedMulAddSIMD(FDPROTO_i2m1, FDPart1_Rational_1_4,
                                                      FusedMulAddSIMD(FDPROTO_i2p2, FDPart1_Rational_1_2, MulSIMD(FDPROTO, FDPart1_Rational_5_6))))));

  return FD_result;
}

/**
 * Set RHSs for the BSSN evolution equations.
 */
void rhs_eval__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                    const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                                    REAL *restrict rhs_gfs) {
#include "../set_CodeParameters-simd.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    const double NOSIMDf3_of_xx2 = rfmstruct->f3_of_xx2[i2];
    MAYBE_UNUSED const REAL_SIMD_ARRAY f3_of_xx2 = ConstSIMD(NOSIMDf3_of_xx2);
    const double NOSIMDf3_of_xx2__D2 = rfmstruct->f3_of_xx2__D2[i2];
    MAYBE_UNUSED const REAL_SIMD_ARRAY f3_of_xx2__D2 = ConstSIMD(NOSIMDf3_of_xx2__D2);
    const double NOSIMDf3_of_xx2__DD22 = rfmstruct->f3_of_xx2__DD22[i2];
    MAYBE_UNUSED const REAL_SIMD_ARRAY f3_of_xx2__DD22 = ConstSIMD(NOSIMDf3_of_xx2__DD22);

    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0 += simd_width) {
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0 = ReadSIMD(&rfmstruct->f0_of_xx0[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__D0 = ReadSIMD(&rfmstruct->f0_of_xx0__D0[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__DD00 = ReadSIMD(&rfmstruct->f0_of_xx0__DD00[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__DDD000 = ReadSIMD(&rfmstruct->f0_of_xx0__DDD000[i0]);
        /*
         * NRPy+-Generated GF Access/FD Code, Step 1 of 3:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL_SIMD_ARRAY RbarDD00 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD01 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD02 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD11 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD12 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY RbarDD22 = ReadSIMD(&auxevol_gfs[IDX4(RBARDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU00 = ReadSIMD(&auxevol_gfs[IDX4(T4UU00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU01 = ReadSIMD(&auxevol_gfs[IDX4(T4UU01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU02 = ReadSIMD(&auxevol_gfs[IDX4(T4UU02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU03 = ReadSIMD(&auxevol_gfs[IDX4(T4UU03GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU11 = ReadSIMD(&auxevol_gfs[IDX4(T4UU11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU12 = ReadSIMD(&auxevol_gfs[IDX4(T4UU12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU13 = ReadSIMD(&auxevol_gfs[IDX4(T4UU13GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU22 = ReadSIMD(&auxevol_gfs[IDX4(T4UU22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU23 = ReadSIMD(&auxevol_gfs[IDX4(T4UU23GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY T4UU33 = ReadSIMD(&auxevol_gfs[IDX4(T4UU33GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i2m3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY aDD00_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD00_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD00_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD00_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD00_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD00_i2p3 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY aDD01_i2m3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY aDD01_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD01_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD01_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD01_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD01_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD01_i2p3 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY aDD02_i2m3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY aDD02_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD02_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD02_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD02_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD02_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD02_i2p3 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY aDD11_i2m3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY aDD11_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD11_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD11_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD11_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD11_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD11_i2p3 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY aDD12_i2m3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY aDD12_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD12_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD12_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD12_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD12_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD12_i2p3 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY aDD22_i2m3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY aDD22_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD22_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD22_i1m3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY aDD22_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD22_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD22_i2p3 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY alpha_i2m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY alpha_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY alpha_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY alpha_i1m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY alpha = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY alpha_i1p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY alpha_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY alpha_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha_i2p3 = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY betU0_i2m3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY betU0_i2m2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY betU0_i2m1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY betU0_i1m3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY betU0_i1m2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY betU0_i1m1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0m3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0m2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0m1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU0 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0p1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0p2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i0p3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU0_i1p1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY betU0_i1p2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY betU0_i1p3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY betU0_i2p1 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY betU0_i2p2 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY betU0_i2p3 = ReadSIMD(&in_gfs[IDX4(BETU0GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY betU1_i2m3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY betU1_i2m2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY betU1_i2m1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY betU1_i1m3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY betU1_i1m2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY betU1_i1m1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0m3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0m2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0m1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0p1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0p2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i0p3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU1_i1p1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY betU1_i1p2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY betU1_i1p3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY betU1_i2p1 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY betU1_i2p2 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY betU1_i2p3 = ReadSIMD(&in_gfs[IDX4(BETU1GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY betU2_i2m3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY betU2_i2m2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY betU2_i2m1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY betU2_i1m3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY betU2_i1m2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY betU2_i1m1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0m3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0m2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0m1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0p1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0p2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i0p3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY betU2_i1p1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY betU2_i1p2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY betU2_i1p3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY betU2_i2p1 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY betU2_i2p2 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY betU2_i2p3 = ReadSIMD(&in_gfs[IDX4(BETU2GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY cf_i2m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY cf_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY cf_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY cf_i1m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY cf = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY cf_i1p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY cf_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY cf_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY cf_i2p3 = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY hDD00_i2m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY hDD00_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD00_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i2p3 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY hDD01_i2m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY hDD01_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD01_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i2p3 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY hDD02_i2m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY hDD02_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD02_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i2p3 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY hDD11_i2m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY hDD11_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD11_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i2p3 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY hDD12_i2m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY hDD12_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD12_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i2p3 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY hDD22_i2m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY hDD22_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i1m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY hDD22_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i2p3 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY lambdaU0_i2m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY lambdaU0_i2m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY lambdaU0_i2m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i2p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY lambdaU0_i2p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY lambdaU0_i2p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY lambdaU1_i2m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY lambdaU1_i2m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY lambdaU1_i2m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i2p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY lambdaU1_i2p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY lambdaU1_i2p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY lambdaU2_i2m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY lambdaU2_i2m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY lambdaU2_i2m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i2p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY lambdaU2_i2p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY lambdaU2_i2p3 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY trK_i2m3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY trK_i2m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY trK_i2m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY trK_i1m3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY trK_i1m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY trK_i1m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY trK = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i1p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY trK_i1p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY trK_i1p3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY trK_i2p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY trK_i2p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY trK_i2p3 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY vetU0_i2m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY vetU0_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU0_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU0_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU0_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU0_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0_i2p3 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY vetU1_i2m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY vetU1_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU1_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU1_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU1_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU1_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU1_i2p3 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY vetU2_i2m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2 - 3)]);
        const REAL_SIMD_ARRAY vetU2_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY vetU2_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY vetU2_i1m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 3, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 3, i2)]);
        const REAL_SIMD_ARRAY vetU2_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY vetU2_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU2_i2p3 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2 + 3)]);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD000 =
            SIMD_fd_function_ddnD0_fdorder4(aDD00, aDD00_i0m1, aDD00_i0m2, aDD00_i0m3, aDD00_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD001 =
            SIMD_fd_function_ddnD1_fdorder4(aDD00, aDD00_i1m1, aDD00_i1m2, aDD00_i1m3, aDD00_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD002 =
            SIMD_fd_function_ddnD2_fdorder4(aDD00, aDD00_i2m1, aDD00_i2m2, aDD00_i2m3, aDD00_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD010 =
            SIMD_fd_function_ddnD0_fdorder4(aDD01, aDD01_i0m1, aDD01_i0m2, aDD01_i0m3, aDD01_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD011 =
            SIMD_fd_function_ddnD1_fdorder4(aDD01, aDD01_i1m1, aDD01_i1m2, aDD01_i1m3, aDD01_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD012 =
            SIMD_fd_function_ddnD2_fdorder4(aDD01, aDD01_i2m1, aDD01_i2m2, aDD01_i2m3, aDD01_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD020 =
            SIMD_fd_function_ddnD0_fdorder4(aDD02, aDD02_i0m1, aDD02_i0m2, aDD02_i0m3, aDD02_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD021 =
            SIMD_fd_function_ddnD1_fdorder4(aDD02, aDD02_i1m1, aDD02_i1m2, aDD02_i1m3, aDD02_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD022 =
            SIMD_fd_function_ddnD2_fdorder4(aDD02, aDD02_i2m1, aDD02_i2m2, aDD02_i2m3, aDD02_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD110 =
            SIMD_fd_function_ddnD0_fdorder4(aDD11, aDD11_i0m1, aDD11_i0m2, aDD11_i0m3, aDD11_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD111 =
            SIMD_fd_function_ddnD1_fdorder4(aDD11, aDD11_i1m1, aDD11_i1m2, aDD11_i1m3, aDD11_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD112 =
            SIMD_fd_function_ddnD2_fdorder4(aDD11, aDD11_i2m1, aDD11_i2m2, aDD11_i2m3, aDD11_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD120 =
            SIMD_fd_function_ddnD0_fdorder4(aDD12, aDD12_i0m1, aDD12_i0m2, aDD12_i0m3, aDD12_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD121 =
            SIMD_fd_function_ddnD1_fdorder4(aDD12, aDD12_i1m1, aDD12_i1m2, aDD12_i1m3, aDD12_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD122 =
            SIMD_fd_function_ddnD2_fdorder4(aDD12, aDD12_i2m1, aDD12_i2m2, aDD12_i2m3, aDD12_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD220 =
            SIMD_fd_function_ddnD0_fdorder4(aDD22, aDD22_i0m1, aDD22_i0m2, aDD22_i0m3, aDD22_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD221 =
            SIMD_fd_function_ddnD1_fdorder4(aDD22, aDD22_i1m1, aDD22_i1m2, aDD22_i1m3, aDD22_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_ddnD222 =
            SIMD_fd_function_ddnD2_fdorder4(aDD22, aDD22_i2m1, aDD22_i2m2, aDD22_i2m3, aDD22_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD000 =
            SIMD_fd_function_dupD0_fdorder4(aDD00, aDD00_i0m1, aDD00_i0p1, aDD00_i0p2, aDD00_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD001 =
            SIMD_fd_function_dupD1_fdorder4(aDD00, aDD00_i1m1, aDD00_i1p1, aDD00_i1p2, aDD00_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD002 =
            SIMD_fd_function_dupD2_fdorder4(aDD00, aDD00_i2m1, aDD00_i2p1, aDD00_i2p2, aDD00_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD010 =
            SIMD_fd_function_dupD0_fdorder4(aDD01, aDD01_i0m1, aDD01_i0p1, aDD01_i0p2, aDD01_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD011 =
            SIMD_fd_function_dupD1_fdorder4(aDD01, aDD01_i1m1, aDD01_i1p1, aDD01_i1p2, aDD01_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD012 =
            SIMD_fd_function_dupD2_fdorder4(aDD01, aDD01_i2m1, aDD01_i2p1, aDD01_i2p2, aDD01_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD020 =
            SIMD_fd_function_dupD0_fdorder4(aDD02, aDD02_i0m1, aDD02_i0p1, aDD02_i0p2, aDD02_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD021 =
            SIMD_fd_function_dupD1_fdorder4(aDD02, aDD02_i1m1, aDD02_i1p1, aDD02_i1p2, aDD02_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD022 =
            SIMD_fd_function_dupD2_fdorder4(aDD02, aDD02_i2m1, aDD02_i2p1, aDD02_i2p2, aDD02_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD110 =
            SIMD_fd_function_dupD0_fdorder4(aDD11, aDD11_i0m1, aDD11_i0p1, aDD11_i0p2, aDD11_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD111 =
            SIMD_fd_function_dupD1_fdorder4(aDD11, aDD11_i1m1, aDD11_i1p1, aDD11_i1p2, aDD11_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD112 =
            SIMD_fd_function_dupD2_fdorder4(aDD11, aDD11_i2m1, aDD11_i2p1, aDD11_i2p2, aDD11_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD120 =
            SIMD_fd_function_dupD0_fdorder4(aDD12, aDD12_i0m1, aDD12_i0p1, aDD12_i0p2, aDD12_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD121 =
            SIMD_fd_function_dupD1_fdorder4(aDD12, aDD12_i1m1, aDD12_i1p1, aDD12_i1p2, aDD12_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD122 =
            SIMD_fd_function_dupD2_fdorder4(aDD12, aDD12_i2m1, aDD12_i2p1, aDD12_i2p2, aDD12_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD220 =
            SIMD_fd_function_dupD0_fdorder4(aDD22, aDD22_i0m1, aDD22_i0p1, aDD22_i0p2, aDD22_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD221 =
            SIMD_fd_function_dupD1_fdorder4(aDD22, aDD22_i1m1, aDD22_i1p1, aDD22_i1p2, aDD22_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputaDD_dupD222 =
            SIMD_fd_function_dupD2_fdorder4(aDD22, aDD22_i2m1, aDD22_i2p1, aDD22_i2p2, aDD22_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_ddnD0 =
            SIMD_fd_function_ddnD0_fdorder4(alpha, alpha_i0m1, alpha_i0m2, alpha_i0m3, alpha_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_ddnD1 =
            SIMD_fd_function_ddnD1_fdorder4(alpha, alpha_i1m1, alpha_i1m2, alpha_i1m3, alpha_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_ddnD2 =
            SIMD_fd_function_ddnD2_fdorder4(alpha, alpha_i2m1, alpha_i2m2, alpha_i2m3, alpha_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_dupD0 =
            SIMD_fd_function_dupD0_fdorder4(alpha, alpha_i0m1, alpha_i0p1, alpha_i0p2, alpha_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_dupD1 =
            SIMD_fd_function_dupD1_fdorder4(alpha, alpha_i1m1, alpha_i1p1, alpha_i1p2, alpha_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputalpha_dupD2 =
            SIMD_fd_function_dupD2_fdorder4(alpha, alpha_i2m1, alpha_i2p1, alpha_i2p2, alpha_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD00 =
            SIMD_fd_function_ddnD0_fdorder4(betU0, betU0_i0m1, betU0_i0m2, betU0_i0m3, betU0_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD01 =
            SIMD_fd_function_ddnD1_fdorder4(betU0, betU0_i1m1, betU0_i1m2, betU0_i1m3, betU0_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD02 =
            SIMD_fd_function_ddnD2_fdorder4(betU0, betU0_i2m1, betU0_i2m2, betU0_i2m3, betU0_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD10 =
            SIMD_fd_function_ddnD0_fdorder4(betU1, betU1_i0m1, betU1_i0m2, betU1_i0m3, betU1_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD11 =
            SIMD_fd_function_ddnD1_fdorder4(betU1, betU1_i1m1, betU1_i1m2, betU1_i1m3, betU1_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD12 =
            SIMD_fd_function_ddnD2_fdorder4(betU1, betU1_i2m1, betU1_i2m2, betU1_i2m3, betU1_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD20 =
            SIMD_fd_function_ddnD0_fdorder4(betU2, betU2_i0m1, betU2_i0m2, betU2_i0m3, betU2_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD21 =
            SIMD_fd_function_ddnD1_fdorder4(betU2, betU2_i1m1, betU2_i1m2, betU2_i1m3, betU2_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_ddnD22 =
            SIMD_fd_function_ddnD2_fdorder4(betU2, betU2_i2m1, betU2_i2m2, betU2_i2m3, betU2_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD00 =
            SIMD_fd_function_dupD0_fdorder4(betU0, betU0_i0m1, betU0_i0p1, betU0_i0p2, betU0_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD01 =
            SIMD_fd_function_dupD1_fdorder4(betU0, betU0_i1m1, betU0_i1p1, betU0_i1p2, betU0_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD02 =
            SIMD_fd_function_dupD2_fdorder4(betU0, betU0_i2m1, betU0_i2p1, betU0_i2p2, betU0_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD10 =
            SIMD_fd_function_dupD0_fdorder4(betU1, betU1_i0m1, betU1_i0p1, betU1_i0p2, betU1_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD11 =
            SIMD_fd_function_dupD1_fdorder4(betU1, betU1_i1m1, betU1_i1p1, betU1_i1p2, betU1_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD12 =
            SIMD_fd_function_dupD2_fdorder4(betU1, betU1_i2m1, betU1_i2p1, betU1_i2p2, betU1_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD20 =
            SIMD_fd_function_dupD0_fdorder4(betU2, betU2_i0m1, betU2_i0p1, betU2_i0p2, betU2_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD21 =
            SIMD_fd_function_dupD1_fdorder4(betU2, betU2_i1m1, betU2_i1p1, betU2_i1p2, betU2_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputbetU_dupD22 =
            SIMD_fd_function_dupD2_fdorder4(betU2, betU2_i2m1, betU2_i2p1, betU2_i2p2, betU2_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_ddnD0 = SIMD_fd_function_ddnD0_fdorder4(cf, cf_i0m1, cf_i0m2, cf_i0m3, cf_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_ddnD1 = SIMD_fd_function_ddnD1_fdorder4(cf, cf_i1m1, cf_i1m2, cf_i1m3, cf_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_ddnD2 = SIMD_fd_function_ddnD2_fdorder4(cf, cf_i2m1, cf_i2m2, cf_i2m3, cf_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_dupD0 = SIMD_fd_function_dupD0_fdorder4(cf, cf_i0m1, cf_i0p1, cf_i0p2, cf_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_dupD1 = SIMD_fd_function_dupD1_fdorder4(cf, cf_i1m1, cf_i1p1, cf_i1p2, cf_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputcf_dupD2 = SIMD_fd_function_dupD2_fdorder4(cf, cf_i2m1, cf_i2p1, cf_i2p2, cf_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD000 =
            SIMD_fd_function_ddnD0_fdorder4(hDD00, hDD00_i0m1, hDD00_i0m2, hDD00_i0m3, hDD00_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD001 =
            SIMD_fd_function_ddnD1_fdorder4(hDD00, hDD00_i1m1, hDD00_i1m2, hDD00_i1m3, hDD00_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD002 =
            SIMD_fd_function_ddnD2_fdorder4(hDD00, hDD00_i2m1, hDD00_i2m2, hDD00_i2m3, hDD00_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD010 =
            SIMD_fd_function_ddnD0_fdorder4(hDD01, hDD01_i0m1, hDD01_i0m2, hDD01_i0m3, hDD01_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD011 =
            SIMD_fd_function_ddnD1_fdorder4(hDD01, hDD01_i1m1, hDD01_i1m2, hDD01_i1m3, hDD01_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD012 =
            SIMD_fd_function_ddnD2_fdorder4(hDD01, hDD01_i2m1, hDD01_i2m2, hDD01_i2m3, hDD01_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD020 =
            SIMD_fd_function_ddnD0_fdorder4(hDD02, hDD02_i0m1, hDD02_i0m2, hDD02_i0m3, hDD02_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD021 =
            SIMD_fd_function_ddnD1_fdorder4(hDD02, hDD02_i1m1, hDD02_i1m2, hDD02_i1m3, hDD02_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD022 =
            SIMD_fd_function_ddnD2_fdorder4(hDD02, hDD02_i2m1, hDD02_i2m2, hDD02_i2m3, hDD02_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD110 =
            SIMD_fd_function_ddnD0_fdorder4(hDD11, hDD11_i0m1, hDD11_i0m2, hDD11_i0m3, hDD11_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD111 =
            SIMD_fd_function_ddnD1_fdorder4(hDD11, hDD11_i1m1, hDD11_i1m2, hDD11_i1m3, hDD11_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD112 =
            SIMD_fd_function_ddnD2_fdorder4(hDD11, hDD11_i2m1, hDD11_i2m2, hDD11_i2m3, hDD11_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD120 =
            SIMD_fd_function_ddnD0_fdorder4(hDD12, hDD12_i0m1, hDD12_i0m2, hDD12_i0m3, hDD12_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD121 =
            SIMD_fd_function_ddnD1_fdorder4(hDD12, hDD12_i1m1, hDD12_i1m2, hDD12_i1m3, hDD12_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD122 =
            SIMD_fd_function_ddnD2_fdorder4(hDD12, hDD12_i2m1, hDD12_i2m2, hDD12_i2m3, hDD12_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD220 =
            SIMD_fd_function_ddnD0_fdorder4(hDD22, hDD22_i0m1, hDD22_i0m2, hDD22_i0m3, hDD22_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD221 =
            SIMD_fd_function_ddnD1_fdorder4(hDD22, hDD22_i1m1, hDD22_i1m2, hDD22_i1m3, hDD22_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_ddnD222 =
            SIMD_fd_function_ddnD2_fdorder4(hDD22, hDD22_i2m1, hDD22_i2m2, hDD22_i2m3, hDD22_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD000 =
            SIMD_fd_function_dupD0_fdorder4(hDD00, hDD00_i0m1, hDD00_i0p1, hDD00_i0p2, hDD00_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD001 =
            SIMD_fd_function_dupD1_fdorder4(hDD00, hDD00_i1m1, hDD00_i1p1, hDD00_i1p2, hDD00_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD002 =
            SIMD_fd_function_dupD2_fdorder4(hDD00, hDD00_i2m1, hDD00_i2p1, hDD00_i2p2, hDD00_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD010 =
            SIMD_fd_function_dupD0_fdorder4(hDD01, hDD01_i0m1, hDD01_i0p1, hDD01_i0p2, hDD01_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD011 =
            SIMD_fd_function_dupD1_fdorder4(hDD01, hDD01_i1m1, hDD01_i1p1, hDD01_i1p2, hDD01_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD012 =
            SIMD_fd_function_dupD2_fdorder4(hDD01, hDD01_i2m1, hDD01_i2p1, hDD01_i2p2, hDD01_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD020 =
            SIMD_fd_function_dupD0_fdorder4(hDD02, hDD02_i0m1, hDD02_i0p1, hDD02_i0p2, hDD02_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD021 =
            SIMD_fd_function_dupD1_fdorder4(hDD02, hDD02_i1m1, hDD02_i1p1, hDD02_i1p2, hDD02_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD022 =
            SIMD_fd_function_dupD2_fdorder4(hDD02, hDD02_i2m1, hDD02_i2p1, hDD02_i2p2, hDD02_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD110 =
            SIMD_fd_function_dupD0_fdorder4(hDD11, hDD11_i0m1, hDD11_i0p1, hDD11_i0p2, hDD11_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD111 =
            SIMD_fd_function_dupD1_fdorder4(hDD11, hDD11_i1m1, hDD11_i1p1, hDD11_i1p2, hDD11_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD112 =
            SIMD_fd_function_dupD2_fdorder4(hDD11, hDD11_i2m1, hDD11_i2p1, hDD11_i2p2, hDD11_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD120 =
            SIMD_fd_function_dupD0_fdorder4(hDD12, hDD12_i0m1, hDD12_i0p1, hDD12_i0p2, hDD12_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD121 =
            SIMD_fd_function_dupD1_fdorder4(hDD12, hDD12_i1m1, hDD12_i1p1, hDD12_i1p2, hDD12_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD122 =
            SIMD_fd_function_dupD2_fdorder4(hDD12, hDD12_i2m1, hDD12_i2p1, hDD12_i2p2, hDD12_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD220 =
            SIMD_fd_function_dupD0_fdorder4(hDD22, hDD22_i0m1, hDD22_i0p1, hDD22_i0p2, hDD22_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD221 =
            SIMD_fd_function_dupD1_fdorder4(hDD22, hDD22_i1m1, hDD22_i1p1, hDD22_i1p2, hDD22_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputhDD_dupD222 =
            SIMD_fd_function_dupD2_fdorder4(hDD22, hDD22_i2m1, hDD22_i2p1, hDD22_i2p2, hDD22_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD00 =
            SIMD_fd_function_ddnD0_fdorder4(lambdaU0, lambdaU0_i0m1, lambdaU0_i0m2, lambdaU0_i0m3, lambdaU0_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD01 =
            SIMD_fd_function_ddnD1_fdorder4(lambdaU0, lambdaU0_i1m1, lambdaU0_i1m2, lambdaU0_i1m3, lambdaU0_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD02 =
            SIMD_fd_function_ddnD2_fdorder4(lambdaU0, lambdaU0_i2m1, lambdaU0_i2m2, lambdaU0_i2m3, lambdaU0_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD10 =
            SIMD_fd_function_ddnD0_fdorder4(lambdaU1, lambdaU1_i0m1, lambdaU1_i0m2, lambdaU1_i0m3, lambdaU1_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD11 =
            SIMD_fd_function_ddnD1_fdorder4(lambdaU1, lambdaU1_i1m1, lambdaU1_i1m2, lambdaU1_i1m3, lambdaU1_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD12 =
            SIMD_fd_function_ddnD2_fdorder4(lambdaU1, lambdaU1_i2m1, lambdaU1_i2m2, lambdaU1_i2m3, lambdaU1_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD20 =
            SIMD_fd_function_ddnD0_fdorder4(lambdaU2, lambdaU2_i0m1, lambdaU2_i0m2, lambdaU2_i0m3, lambdaU2_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD21 =
            SIMD_fd_function_ddnD1_fdorder4(lambdaU2, lambdaU2_i1m1, lambdaU2_i1m2, lambdaU2_i1m3, lambdaU2_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_ddnD22 =
            SIMD_fd_function_ddnD2_fdorder4(lambdaU2, lambdaU2_i2m1, lambdaU2_i2m2, lambdaU2_i2m3, lambdaU2_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD00 =
            SIMD_fd_function_dupD0_fdorder4(lambdaU0, lambdaU0_i0m1, lambdaU0_i0p1, lambdaU0_i0p2, lambdaU0_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD01 =
            SIMD_fd_function_dupD1_fdorder4(lambdaU0, lambdaU0_i1m1, lambdaU0_i1p1, lambdaU0_i1p2, lambdaU0_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD02 =
            SIMD_fd_function_dupD2_fdorder4(lambdaU0, lambdaU0_i2m1, lambdaU0_i2p1, lambdaU0_i2p2, lambdaU0_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD10 =
            SIMD_fd_function_dupD0_fdorder4(lambdaU1, lambdaU1_i0m1, lambdaU1_i0p1, lambdaU1_i0p2, lambdaU1_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD11 =
            SIMD_fd_function_dupD1_fdorder4(lambdaU1, lambdaU1_i1m1, lambdaU1_i1p1, lambdaU1_i1p2, lambdaU1_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD12 =
            SIMD_fd_function_dupD2_fdorder4(lambdaU1, lambdaU1_i2m1, lambdaU1_i2p1, lambdaU1_i2p2, lambdaU1_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD20 =
            SIMD_fd_function_dupD0_fdorder4(lambdaU2, lambdaU2_i0m1, lambdaU2_i0p1, lambdaU2_i0p2, lambdaU2_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD21 =
            SIMD_fd_function_dupD1_fdorder4(lambdaU2, lambdaU2_i1m1, lambdaU2_i1p1, lambdaU2_i1p2, lambdaU2_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputlambdaU_dupD22 =
            SIMD_fd_function_dupD2_fdorder4(lambdaU2, lambdaU2_i2m1, lambdaU2_i2p1, lambdaU2_i2p2, lambdaU2_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_ddnD0 = SIMD_fd_function_ddnD0_fdorder4(trK, trK_i0m1, trK_i0m2, trK_i0m3, trK_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_ddnD1 = SIMD_fd_function_ddnD1_fdorder4(trK, trK_i1m1, trK_i1m2, trK_i1m3, trK_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_ddnD2 = SIMD_fd_function_ddnD2_fdorder4(trK, trK_i2m1, trK_i2m2, trK_i2m3, trK_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_dupD0 = SIMD_fd_function_dupD0_fdorder4(trK, trK_i0m1, trK_i0p1, trK_i0p2, trK_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_dupD1 = SIMD_fd_function_dupD1_fdorder4(trK, trK_i1m1, trK_i1p1, trK_i1p2, trK_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputtrK_dupD2 = SIMD_fd_function_dupD2_fdorder4(trK, trK_i2m1, trK_i2p1, trK_i2p2, trK_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD00 =
            SIMD_fd_function_ddnD0_fdorder4(vetU0, vetU0_i0m1, vetU0_i0m2, vetU0_i0m3, vetU0_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD01 =
            SIMD_fd_function_ddnD1_fdorder4(vetU0, vetU0_i1m1, vetU0_i1m2, vetU0_i1m3, vetU0_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD02 =
            SIMD_fd_function_ddnD2_fdorder4(vetU0, vetU0_i2m1, vetU0_i2m2, vetU0_i2m3, vetU0_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD10 =
            SIMD_fd_function_ddnD0_fdorder4(vetU1, vetU1_i0m1, vetU1_i0m2, vetU1_i0m3, vetU1_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD11 =
            SIMD_fd_function_ddnD1_fdorder4(vetU1, vetU1_i1m1, vetU1_i1m2, vetU1_i1m3, vetU1_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD12 =
            SIMD_fd_function_ddnD2_fdorder4(vetU1, vetU1_i2m1, vetU1_i2m2, vetU1_i2m3, vetU1_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD20 =
            SIMD_fd_function_ddnD0_fdorder4(vetU2, vetU2_i0m1, vetU2_i0m2, vetU2_i0m3, vetU2_i0p1, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD21 =
            SIMD_fd_function_ddnD1_fdorder4(vetU2, vetU2_i1m1, vetU2_i1m2, vetU2_i1m3, vetU2_i1p1, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_ddnD22 =
            SIMD_fd_function_ddnD2_fdorder4(vetU2, vetU2_i2m1, vetU2_i2m2, vetU2_i2m3, vetU2_i2p1, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD00 =
            SIMD_fd_function_dupD0_fdorder4(vetU0, vetU0_i0m1, vetU0_i0p1, vetU0_i0p2, vetU0_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD01 =
            SIMD_fd_function_dupD1_fdorder4(vetU0, vetU0_i1m1, vetU0_i1p1, vetU0_i1p2, vetU0_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD02 =
            SIMD_fd_function_dupD2_fdorder4(vetU0, vetU0_i2m1, vetU0_i2p1, vetU0_i2p2, vetU0_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD10 =
            SIMD_fd_function_dupD0_fdorder4(vetU1, vetU1_i0m1, vetU1_i0p1, vetU1_i0p2, vetU1_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD11 =
            SIMD_fd_function_dupD1_fdorder4(vetU1, vetU1_i1m1, vetU1_i1p1, vetU1_i1p2, vetU1_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD12 =
            SIMD_fd_function_dupD2_fdorder4(vetU1, vetU1_i2m1, vetU1_i2p1, vetU1_i2p2, vetU1_i2p3, invdxx2);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD20 =
            SIMD_fd_function_dupD0_fdorder4(vetU2, vetU2_i0m1, vetU2_i0p1, vetU2_i0p2, vetU2_i0p3, invdxx0);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD21 =
            SIMD_fd_function_dupD1_fdorder4(vetU2, vetU2_i1m1, vetU2_i1p1, vetU2_i1p2, vetU2_i1p3, invdxx1);
        const REAL_SIMD_ARRAY UpwindAlgInputvetU_dupD22 =
            SIMD_fd_function_dupD2_fdorder4(vetU2, vetU2_i2m1, vetU2_i2p1, vetU2_i2p2, vetU2_i2p3, invdxx2);
        const REAL_SIMD_ARRAY alpha_dD0 = SIMD_fd_function_dD0_fdorder4(alpha_i0m1, alpha_i0m2, alpha_i0p1, alpha_i0p2, invdxx0);
        const REAL_SIMD_ARRAY alpha_dD1 = SIMD_fd_function_dD1_fdorder4(alpha_i1m1, alpha_i1m2, alpha_i1p1, alpha_i1p2, invdxx1);
        const REAL_SIMD_ARRAY alpha_dD2 = SIMD_fd_function_dD2_fdorder4(alpha_i2m1, alpha_i2m2, alpha_i2p1, alpha_i2p2, invdxx2);
        const REAL_SIMD_ARRAY alpha_dDD00 = SIMD_fd_function_dDD00_fdorder4(alpha, alpha_i0m1, alpha_i0m2, alpha_i0p1, alpha_i0p2, invdxx0);
        const REAL_SIMD_ARRAY alpha_dDD01 =
            SIMD_fd_function_dDD01_fdorder4(alpha_i0m1_i1m1, alpha_i0m1_i1m2, alpha_i0m1_i1p1, alpha_i0m1_i1p2, alpha_i0m2_i1m1, alpha_i0m2_i1m2,
                                            alpha_i0m2_i1p1, alpha_i0m2_i1p2, alpha_i0p1_i1m1, alpha_i0p1_i1m2, alpha_i0p1_i1p1, alpha_i0p1_i1p2,
                                            alpha_i0p2_i1m1, alpha_i0p2_i1m2, alpha_i0p2_i1p1, alpha_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY alpha_dDD02 =
            SIMD_fd_function_dDD02_fdorder4(alpha_i0m1_i2m1, alpha_i0m1_i2m2, alpha_i0m1_i2p1, alpha_i0m1_i2p2, alpha_i0m2_i2m1, alpha_i0m2_i2m2,
                                            alpha_i0m2_i2p1, alpha_i0m2_i2p2, alpha_i0p1_i2m1, alpha_i0p1_i2m2, alpha_i0p1_i2p1, alpha_i0p1_i2p2,
                                            alpha_i0p2_i2m1, alpha_i0p2_i2m2, alpha_i0p2_i2p1, alpha_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY alpha_dDD11 = SIMD_fd_function_dDD11_fdorder4(alpha, alpha_i1m1, alpha_i1m2, alpha_i1p1, alpha_i1p2, invdxx1);
        const REAL_SIMD_ARRAY alpha_dDD12 =
            SIMD_fd_function_dDD12_fdorder4(alpha_i1m1_i2m1, alpha_i1m1_i2m2, alpha_i1m1_i2p1, alpha_i1m1_i2p2, alpha_i1m2_i2m1, alpha_i1m2_i2m2,
                                            alpha_i1m2_i2p1, alpha_i1m2_i2p2, alpha_i1p1_i2m1, alpha_i1p1_i2m2, alpha_i1p1_i2p1, alpha_i1p1_i2p2,
                                            alpha_i1p2_i2m1, alpha_i1p2_i2m2, alpha_i1p2_i2p1, alpha_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY alpha_dDD22 = SIMD_fd_function_dDD22_fdorder4(alpha, alpha_i2m1, alpha_i2m2, alpha_i2p1, alpha_i2p2, invdxx2);
        const REAL_SIMD_ARRAY cf_dD0 = SIMD_fd_function_dD0_fdorder4(cf_i0m1, cf_i0m2, cf_i0p1, cf_i0p2, invdxx0);
        const REAL_SIMD_ARRAY cf_dD1 = SIMD_fd_function_dD1_fdorder4(cf_i1m1, cf_i1m2, cf_i1p1, cf_i1p2, invdxx1);
        const REAL_SIMD_ARRAY cf_dD2 = SIMD_fd_function_dD2_fdorder4(cf_i2m1, cf_i2m2, cf_i2p1, cf_i2p2, invdxx2);
        const REAL_SIMD_ARRAY cf_dDD00 = SIMD_fd_function_dDD00_fdorder4(cf, cf_i0m1, cf_i0m2, cf_i0p1, cf_i0p2, invdxx0);
        const REAL_SIMD_ARRAY cf_dDD01 = SIMD_fd_function_dDD01_fdorder4(
            cf_i0m1_i1m1, cf_i0m1_i1m2, cf_i0m1_i1p1, cf_i0m1_i1p2, cf_i0m2_i1m1, cf_i0m2_i1m2, cf_i0m2_i1p1, cf_i0m2_i1p2, cf_i0p1_i1m1,
            cf_i0p1_i1m2, cf_i0p1_i1p1, cf_i0p1_i1p2, cf_i0p2_i1m1, cf_i0p2_i1m2, cf_i0p2_i1p1, cf_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY cf_dDD02 = SIMD_fd_function_dDD02_fdorder4(
            cf_i0m1_i2m1, cf_i0m1_i2m2, cf_i0m1_i2p1, cf_i0m1_i2p2, cf_i0m2_i2m1, cf_i0m2_i2m2, cf_i0m2_i2p1, cf_i0m2_i2p2, cf_i0p1_i2m1,
            cf_i0p1_i2m2, cf_i0p1_i2p1, cf_i0p1_i2p2, cf_i0p2_i2m1, cf_i0p2_i2m2, cf_i0p2_i2p1, cf_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY cf_dDD11 = SIMD_fd_function_dDD11_fdorder4(cf, cf_i1m1, cf_i1m2, cf_i1p1, cf_i1p2, invdxx1);
        const REAL_SIMD_ARRAY cf_dDD12 = SIMD_fd_function_dDD12_fdorder4(
            cf_i1m1_i2m1, cf_i1m1_i2m2, cf_i1m1_i2p1, cf_i1m1_i2p2, cf_i1m2_i2m1, cf_i1m2_i2m2, cf_i1m2_i2p1, cf_i1m2_i2p2, cf_i1p1_i2m1,
            cf_i1p1_i2m2, cf_i1p1_i2p1, cf_i1p1_i2p2, cf_i1p2_i2m1, cf_i1p2_i2m2, cf_i1p2_i2p1, cf_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY cf_dDD22 = SIMD_fd_function_dDD22_fdorder4(cf, cf_i2m1, cf_i2m2, cf_i2p1, cf_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dD000 = SIMD_fd_function_dD0_fdorder4(hDD00_i0m1, hDD00_i0m2, hDD00_i0p1, hDD00_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD001 = SIMD_fd_function_dD1_fdorder4(hDD00_i1m1, hDD00_i1m2, hDD00_i1p1, hDD00_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD002 = SIMD_fd_function_dD2_fdorder4(hDD00_i2m1, hDD00_i2m2, hDD00_i2p1, hDD00_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dD010 = SIMD_fd_function_dD0_fdorder4(hDD01_i0m1, hDD01_i0m2, hDD01_i0p1, hDD01_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD011 = SIMD_fd_function_dD1_fdorder4(hDD01_i1m1, hDD01_i1m2, hDD01_i1p1, hDD01_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD012 = SIMD_fd_function_dD2_fdorder4(hDD01_i2m1, hDD01_i2m2, hDD01_i2p1, hDD01_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dD020 = SIMD_fd_function_dD0_fdorder4(hDD02_i0m1, hDD02_i0m2, hDD02_i0p1, hDD02_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD021 = SIMD_fd_function_dD1_fdorder4(hDD02_i1m1, hDD02_i1m2, hDD02_i1p1, hDD02_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD022 = SIMD_fd_function_dD2_fdorder4(hDD02_i2m1, hDD02_i2m2, hDD02_i2p1, hDD02_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dD110 = SIMD_fd_function_dD0_fdorder4(hDD11_i0m1, hDD11_i0m2, hDD11_i0p1, hDD11_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD111 = SIMD_fd_function_dD1_fdorder4(hDD11_i1m1, hDD11_i1m2, hDD11_i1p1, hDD11_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD112 = SIMD_fd_function_dD2_fdorder4(hDD11_i2m1, hDD11_i2m2, hDD11_i2p1, hDD11_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dD120 = SIMD_fd_function_dD0_fdorder4(hDD12_i0m1, hDD12_i0m2, hDD12_i0p1, hDD12_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD121 = SIMD_fd_function_dD1_fdorder4(hDD12_i1m1, hDD12_i1m2, hDD12_i1p1, hDD12_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD122 = SIMD_fd_function_dD2_fdorder4(hDD12_i2m1, hDD12_i2m2, hDD12_i2p1, hDD12_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dD220 = SIMD_fd_function_dD0_fdorder4(hDD22_i0m1, hDD22_i0m2, hDD22_i0p1, hDD22_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dD221 = SIMD_fd_function_dD1_fdorder4(hDD22_i1m1, hDD22_i1m2, hDD22_i1p1, hDD22_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dD222 = SIMD_fd_function_dD2_fdorder4(hDD22_i2m1, hDD22_i2m2, hDD22_i2p1, hDD22_i2p2, invdxx2);
        const REAL_SIMD_ARRAY trK_dD0 = SIMD_fd_function_dD0_fdorder4(trK_i0m1, trK_i0m2, trK_i0p1, trK_i0p2, invdxx0);
        const REAL_SIMD_ARRAY trK_dD1 = SIMD_fd_function_dD1_fdorder4(trK_i1m1, trK_i1m2, trK_i1p1, trK_i1p2, invdxx1);
        const REAL_SIMD_ARRAY trK_dD2 = SIMD_fd_function_dD2_fdorder4(trK_i2m1, trK_i2m2, trK_i2p1, trK_i2p2, invdxx2);
        const REAL_SIMD_ARRAY vetU_dD00 = SIMD_fd_function_dD0_fdorder4(vetU0_i0m1, vetU0_i0m2, vetU0_i0p1, vetU0_i0p2, invdxx0);
        const REAL_SIMD_ARRAY vetU_dD01 = SIMD_fd_function_dD1_fdorder4(vetU0_i1m1, vetU0_i1m2, vetU0_i1p1, vetU0_i1p2, invdxx1);
        const REAL_SIMD_ARRAY vetU_dD02 = SIMD_fd_function_dD2_fdorder4(vetU0_i2m1, vetU0_i2m2, vetU0_i2p1, vetU0_i2p2, invdxx2);
        const REAL_SIMD_ARRAY vetU_dD10 = SIMD_fd_function_dD0_fdorder4(vetU1_i0m1, vetU1_i0m2, vetU1_i0p1, vetU1_i0p2, invdxx0);
        const REAL_SIMD_ARRAY vetU_dD11 = SIMD_fd_function_dD1_fdorder4(vetU1_i1m1, vetU1_i1m2, vetU1_i1p1, vetU1_i1p2, invdxx1);
        const REAL_SIMD_ARRAY vetU_dD12 = SIMD_fd_function_dD2_fdorder4(vetU1_i2m1, vetU1_i2m2, vetU1_i2p1, vetU1_i2p2, invdxx2);
        const REAL_SIMD_ARRAY vetU_dD20 = SIMD_fd_function_dD0_fdorder4(vetU2_i0m1, vetU2_i0m2, vetU2_i0p1, vetU2_i0p2, invdxx0);
        const REAL_SIMD_ARRAY vetU_dD21 = SIMD_fd_function_dD1_fdorder4(vetU2_i1m1, vetU2_i1m2, vetU2_i1p1, vetU2_i1p2, invdxx1);
        const REAL_SIMD_ARRAY vetU_dD22 = SIMD_fd_function_dD2_fdorder4(vetU2_i2m1, vetU2_i2m2, vetU2_i2p1, vetU2_i2p2, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD000 = SIMD_fd_function_dDD00_fdorder4(vetU0, vetU0_i0m1, vetU0_i0m2, vetU0_i0p1, vetU0_i0p2, invdxx0);
        const REAL_SIMD_ARRAY vetU_dDD001 =
            SIMD_fd_function_dDD01_fdorder4(vetU0_i0m1_i1m1, vetU0_i0m1_i1m2, vetU0_i0m1_i1p1, vetU0_i0m1_i1p2, vetU0_i0m2_i1m1, vetU0_i0m2_i1m2,
                                            vetU0_i0m2_i1p1, vetU0_i0m2_i1p2, vetU0_i0p1_i1m1, vetU0_i0p1_i1m2, vetU0_i0p1_i1p1, vetU0_i0p1_i1p2,
                                            vetU0_i0p2_i1m1, vetU0_i0p2_i1m2, vetU0_i0p2_i1p1, vetU0_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD002 =
            SIMD_fd_function_dDD02_fdorder4(vetU0_i0m1_i2m1, vetU0_i0m1_i2m2, vetU0_i0m1_i2p1, vetU0_i0m1_i2p2, vetU0_i0m2_i2m1, vetU0_i0m2_i2m2,
                                            vetU0_i0m2_i2p1, vetU0_i0m2_i2p2, vetU0_i0p1_i2m1, vetU0_i0p1_i2m2, vetU0_i0p1_i2p1, vetU0_i0p1_i2p2,
                                            vetU0_i0p2_i2m1, vetU0_i0p2_i2m2, vetU0_i0p2_i2p1, vetU0_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD011 = SIMD_fd_function_dDD11_fdorder4(vetU0, vetU0_i1m1, vetU0_i1m2, vetU0_i1p1, vetU0_i1p2, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD012 =
            SIMD_fd_function_dDD12_fdorder4(vetU0_i1m1_i2m1, vetU0_i1m1_i2m2, vetU0_i1m1_i2p1, vetU0_i1m1_i2p2, vetU0_i1m2_i2m1, vetU0_i1m2_i2m2,
                                            vetU0_i1m2_i2p1, vetU0_i1m2_i2p2, vetU0_i1p1_i2m1, vetU0_i1p1_i2m2, vetU0_i1p1_i2p1, vetU0_i1p1_i2p2,
                                            vetU0_i1p2_i2m1, vetU0_i1p2_i2m2, vetU0_i1p2_i2p1, vetU0_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD022 = SIMD_fd_function_dDD22_fdorder4(vetU0, vetU0_i2m1, vetU0_i2m2, vetU0_i2p1, vetU0_i2p2, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD100 = SIMD_fd_function_dDD00_fdorder4(vetU1, vetU1_i0m1, vetU1_i0m2, vetU1_i0p1, vetU1_i0p2, invdxx0);
        const REAL_SIMD_ARRAY vetU_dDD101 =
            SIMD_fd_function_dDD01_fdorder4(vetU1_i0m1_i1m1, vetU1_i0m1_i1m2, vetU1_i0m1_i1p1, vetU1_i0m1_i1p2, vetU1_i0m2_i1m1, vetU1_i0m2_i1m2,
                                            vetU1_i0m2_i1p1, vetU1_i0m2_i1p2, vetU1_i0p1_i1m1, vetU1_i0p1_i1m2, vetU1_i0p1_i1p1, vetU1_i0p1_i1p2,
                                            vetU1_i0p2_i1m1, vetU1_i0p2_i1m2, vetU1_i0p2_i1p1, vetU1_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD102 =
            SIMD_fd_function_dDD02_fdorder4(vetU1_i0m1_i2m1, vetU1_i0m1_i2m2, vetU1_i0m1_i2p1, vetU1_i0m1_i2p2, vetU1_i0m2_i2m1, vetU1_i0m2_i2m2,
                                            vetU1_i0m2_i2p1, vetU1_i0m2_i2p2, vetU1_i0p1_i2m1, vetU1_i0p1_i2m2, vetU1_i0p1_i2p1, vetU1_i0p1_i2p2,
                                            vetU1_i0p2_i2m1, vetU1_i0p2_i2m2, vetU1_i0p2_i2p1, vetU1_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD111 = SIMD_fd_function_dDD11_fdorder4(vetU1, vetU1_i1m1, vetU1_i1m2, vetU1_i1p1, vetU1_i1p2, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD112 =
            SIMD_fd_function_dDD12_fdorder4(vetU1_i1m1_i2m1, vetU1_i1m1_i2m2, vetU1_i1m1_i2p1, vetU1_i1m1_i2p2, vetU1_i1m2_i2m1, vetU1_i1m2_i2m2,
                                            vetU1_i1m2_i2p1, vetU1_i1m2_i2p2, vetU1_i1p1_i2m1, vetU1_i1p1_i2m2, vetU1_i1p1_i2p1, vetU1_i1p1_i2p2,
                                            vetU1_i1p2_i2m1, vetU1_i1p2_i2m2, vetU1_i1p2_i2p1, vetU1_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD122 = SIMD_fd_function_dDD22_fdorder4(vetU1, vetU1_i2m1, vetU1_i2m2, vetU1_i2p1, vetU1_i2p2, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD200 = SIMD_fd_function_dDD00_fdorder4(vetU2, vetU2_i0m1, vetU2_i0m2, vetU2_i0p1, vetU2_i0p2, invdxx0);
        const REAL_SIMD_ARRAY vetU_dDD201 =
            SIMD_fd_function_dDD01_fdorder4(vetU2_i0m1_i1m1, vetU2_i0m1_i1m2, vetU2_i0m1_i1p1, vetU2_i0m1_i1p2, vetU2_i0m2_i1m1, vetU2_i0m2_i1m2,
                                            vetU2_i0m2_i1p1, vetU2_i0m2_i1p2, vetU2_i0p1_i1m1, vetU2_i0p1_i1m2, vetU2_i0p1_i1p1, vetU2_i0p1_i1p2,
                                            vetU2_i0p2_i1m1, vetU2_i0p2_i1m2, vetU2_i0p2_i1p1, vetU2_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD202 =
            SIMD_fd_function_dDD02_fdorder4(vetU2_i0m1_i2m1, vetU2_i0m1_i2m2, vetU2_i0m1_i2p1, vetU2_i0m1_i2p2, vetU2_i0m2_i2m1, vetU2_i0m2_i2m2,
                                            vetU2_i0m2_i2p1, vetU2_i0m2_i2p2, vetU2_i0p1_i2m1, vetU2_i0p1_i2m2, vetU2_i0p1_i2p1, vetU2_i0p1_i2p2,
                                            vetU2_i0p2_i2m1, vetU2_i0p2_i2m2, vetU2_i0p2_i2p1, vetU2_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD211 = SIMD_fd_function_dDD11_fdorder4(vetU2, vetU2_i1m1, vetU2_i1m2, vetU2_i1p1, vetU2_i1p2, invdxx1);
        const REAL_SIMD_ARRAY vetU_dDD212 =
            SIMD_fd_function_dDD12_fdorder4(vetU2_i1m1_i2m1, vetU2_i1m1_i2m2, vetU2_i1m1_i2p1, vetU2_i1m1_i2p2, vetU2_i1m2_i2m1, vetU2_i1m2_i2m2,
                                            vetU2_i1m2_i2p1, vetU2_i1m2_i2p2, vetU2_i1p1_i2m1, vetU2_i1p1_i2m2, vetU2_i1p1_i2p1, vetU2_i1p1_i2p2,
                                            vetU2_i1p2_i2m1, vetU2_i1p2_i2m2, vetU2_i1p2_i2p1, vetU2_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY vetU_dDD222 = SIMD_fd_function_dDD22_fdorder4(vetU2, vetU2_i2m1, vetU2_i2m2, vetU2_i2p1, vetU2_i2p2, invdxx2);
        static const double dblFDPart1_Integer_1 = 1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart1_Integer_1 = ConstSIMD(dblFDPart1_Integer_1);

        static const double dblFDPart1_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart1_NegativeOne_ = ConstSIMD(dblFDPart1_NegativeOne_);

        const REAL_SIMD_ARRAY UpwindControlVectorU0 = DivSIMD(vetU0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY UpwindControlVectorU1 = DivSIMD(vetU1, f0_of_xx0);
        const REAL_SIMD_ARRAY UpwindControlVectorU2 = DivSIMD(vetU2, f3_of_xx2);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 3:
         * Implement upwinding algorithm.
         */
        const double tmp_upwind_Integer_0 = 0.000000000000000000000000000000000;

        const REAL_SIMD_ARRAY upwind_Integer_0 = ConstSIMD(tmp_upwind_Integer_0);
        const double tmp_upwind_Integer_1 = 1.000000000000000000000000000000000;

        const REAL_SIMD_ARRAY upwind_Integer_1 = ConstSIMD(tmp_upwind_Integer_1);
        const REAL_SIMD_ARRAY Upwind0 = UPWIND_ALG(UpwindControlVectorU0);
        const REAL_SIMD_ARRAY Upwind1 = UPWIND_ALG(UpwindControlVectorU1);
        const REAL_SIMD_ARRAY Upwind2 = UPWIND_ALG(UpwindControlVectorU2);
        static const double dblFDPart2_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart2_NegativeOne_ = ConstSIMD(dblFDPart2_NegativeOne_);

        const REAL_SIMD_ARRAY aDD_dupD000 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD000, UpwindAlgInputaDD_ddnD000), UpwindAlgInputaDD_ddnD000);
        const REAL_SIMD_ARRAY aDD_dupD001 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD001, UpwindAlgInputaDD_ddnD001), UpwindAlgInputaDD_ddnD001);
        const REAL_SIMD_ARRAY aDD_dupD002 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputaDD_dupD002, UpwindAlgInputaDD_ddnD002), UpwindAlgInputaDD_ddnD002);
        const REAL_SIMD_ARRAY aDD_dupD010 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD010, UpwindAlgInputaDD_ddnD010), UpwindAlgInputaDD_ddnD010);
        const REAL_SIMD_ARRAY aDD_dupD011 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD011, UpwindAlgInputaDD_ddnD011), UpwindAlgInputaDD_ddnD011);
        const REAL_SIMD_ARRAY aDD_dupD012 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputaDD_dupD012, UpwindAlgInputaDD_ddnD012), UpwindAlgInputaDD_ddnD012);
        const REAL_SIMD_ARRAY aDD_dupD020 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD020, UpwindAlgInputaDD_ddnD020), UpwindAlgInputaDD_ddnD020);
        const REAL_SIMD_ARRAY aDD_dupD021 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD021, UpwindAlgInputaDD_ddnD021), UpwindAlgInputaDD_ddnD021);
        const REAL_SIMD_ARRAY aDD_dupD022 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputaDD_dupD022, UpwindAlgInputaDD_ddnD022), UpwindAlgInputaDD_ddnD022);
        const REAL_SIMD_ARRAY aDD_dupD110 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD110, UpwindAlgInputaDD_ddnD110), UpwindAlgInputaDD_ddnD110);
        const REAL_SIMD_ARRAY aDD_dupD111 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD111, UpwindAlgInputaDD_ddnD111), UpwindAlgInputaDD_ddnD111);
        const REAL_SIMD_ARRAY aDD_dupD112 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputaDD_dupD112, UpwindAlgInputaDD_ddnD112), UpwindAlgInputaDD_ddnD112);
        const REAL_SIMD_ARRAY aDD_dupD120 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD120, UpwindAlgInputaDD_ddnD120), UpwindAlgInputaDD_ddnD120);
        const REAL_SIMD_ARRAY aDD_dupD121 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD121, UpwindAlgInputaDD_ddnD121), UpwindAlgInputaDD_ddnD121);
        const REAL_SIMD_ARRAY aDD_dupD122 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputaDD_dupD122, UpwindAlgInputaDD_ddnD122), UpwindAlgInputaDD_ddnD122);
        const REAL_SIMD_ARRAY aDD_dupD220 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputaDD_dupD220, UpwindAlgInputaDD_ddnD220), UpwindAlgInputaDD_ddnD220);
        const REAL_SIMD_ARRAY aDD_dupD221 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputaDD_dupD221, UpwindAlgInputaDD_ddnD221), UpwindAlgInputaDD_ddnD221);
        const REAL_SIMD_ARRAY aDD_dupD222 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputaDD_dupD222, UpwindAlgInputaDD_ddnD222), UpwindAlgInputaDD_ddnD222);
        const REAL_SIMD_ARRAY alpha_dupD0 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputalpha_dupD0, UpwindAlgInputalpha_ddnD0), UpwindAlgInputalpha_ddnD0);
        const REAL_SIMD_ARRAY alpha_dupD1 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputalpha_dupD1, UpwindAlgInputalpha_ddnD1), UpwindAlgInputalpha_ddnD1);
        const REAL_SIMD_ARRAY alpha_dupD2 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputalpha_dupD2, UpwindAlgInputalpha_ddnD2), UpwindAlgInputalpha_ddnD2);
        const REAL_SIMD_ARRAY betU_dupD00 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputbetU_dupD00, UpwindAlgInputbetU_ddnD00), UpwindAlgInputbetU_ddnD00);
        const REAL_SIMD_ARRAY betU_dupD01 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputbetU_dupD01, UpwindAlgInputbetU_ddnD01), UpwindAlgInputbetU_ddnD01);
        const REAL_SIMD_ARRAY betU_dupD02 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputbetU_dupD02, UpwindAlgInputbetU_ddnD02), UpwindAlgInputbetU_ddnD02);
        const REAL_SIMD_ARRAY betU_dupD10 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputbetU_dupD10, UpwindAlgInputbetU_ddnD10), UpwindAlgInputbetU_ddnD10);
        const REAL_SIMD_ARRAY betU_dupD11 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputbetU_dupD11, UpwindAlgInputbetU_ddnD11), UpwindAlgInputbetU_ddnD11);
        const REAL_SIMD_ARRAY betU_dupD12 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputbetU_dupD12, UpwindAlgInputbetU_ddnD12), UpwindAlgInputbetU_ddnD12);
        const REAL_SIMD_ARRAY betU_dupD20 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputbetU_dupD20, UpwindAlgInputbetU_ddnD20), UpwindAlgInputbetU_ddnD20);
        const REAL_SIMD_ARRAY betU_dupD21 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputbetU_dupD21, UpwindAlgInputbetU_ddnD21), UpwindAlgInputbetU_ddnD21);
        const REAL_SIMD_ARRAY betU_dupD22 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputbetU_dupD22, UpwindAlgInputbetU_ddnD22), UpwindAlgInputbetU_ddnD22);
        const REAL_SIMD_ARRAY cf_dupD0 = FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputcf_dupD0, UpwindAlgInputcf_ddnD0), UpwindAlgInputcf_ddnD0);
        const REAL_SIMD_ARRAY cf_dupD1 = FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputcf_dupD1, UpwindAlgInputcf_ddnD1), UpwindAlgInputcf_ddnD1);
        const REAL_SIMD_ARRAY cf_dupD2 = FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputcf_dupD2, UpwindAlgInputcf_ddnD2), UpwindAlgInputcf_ddnD2);
        const REAL_SIMD_ARRAY hDD_dupD000 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD000, UpwindAlgInputhDD_ddnD000), UpwindAlgInputhDD_ddnD000);
        const REAL_SIMD_ARRAY hDD_dupD001 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD001, UpwindAlgInputhDD_ddnD001), UpwindAlgInputhDD_ddnD001);
        const REAL_SIMD_ARRAY hDD_dupD002 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputhDD_dupD002, UpwindAlgInputhDD_ddnD002), UpwindAlgInputhDD_ddnD002);
        const REAL_SIMD_ARRAY hDD_dupD010 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD010, UpwindAlgInputhDD_ddnD010), UpwindAlgInputhDD_ddnD010);
        const REAL_SIMD_ARRAY hDD_dupD011 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD011, UpwindAlgInputhDD_ddnD011), UpwindAlgInputhDD_ddnD011);
        const REAL_SIMD_ARRAY hDD_dupD012 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputhDD_dupD012, UpwindAlgInputhDD_ddnD012), UpwindAlgInputhDD_ddnD012);
        const REAL_SIMD_ARRAY hDD_dupD020 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD020, UpwindAlgInputhDD_ddnD020), UpwindAlgInputhDD_ddnD020);
        const REAL_SIMD_ARRAY hDD_dupD021 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD021, UpwindAlgInputhDD_ddnD021), UpwindAlgInputhDD_ddnD021);
        const REAL_SIMD_ARRAY hDD_dupD022 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputhDD_dupD022, UpwindAlgInputhDD_ddnD022), UpwindAlgInputhDD_ddnD022);
        const REAL_SIMD_ARRAY hDD_dupD110 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD110, UpwindAlgInputhDD_ddnD110), UpwindAlgInputhDD_ddnD110);
        const REAL_SIMD_ARRAY hDD_dupD111 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD111, UpwindAlgInputhDD_ddnD111), UpwindAlgInputhDD_ddnD111);
        const REAL_SIMD_ARRAY hDD_dupD112 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputhDD_dupD112, UpwindAlgInputhDD_ddnD112), UpwindAlgInputhDD_ddnD112);
        const REAL_SIMD_ARRAY hDD_dupD120 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD120, UpwindAlgInputhDD_ddnD120), UpwindAlgInputhDD_ddnD120);
        const REAL_SIMD_ARRAY hDD_dupD121 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD121, UpwindAlgInputhDD_ddnD121), UpwindAlgInputhDD_ddnD121);
        const REAL_SIMD_ARRAY hDD_dupD122 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputhDD_dupD122, UpwindAlgInputhDD_ddnD122), UpwindAlgInputhDD_ddnD122);
        const REAL_SIMD_ARRAY hDD_dupD220 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputhDD_dupD220, UpwindAlgInputhDD_ddnD220), UpwindAlgInputhDD_ddnD220);
        const REAL_SIMD_ARRAY hDD_dupD221 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputhDD_dupD221, UpwindAlgInputhDD_ddnD221), UpwindAlgInputhDD_ddnD221);
        const REAL_SIMD_ARRAY hDD_dupD222 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputhDD_dupD222, UpwindAlgInputhDD_ddnD222), UpwindAlgInputhDD_ddnD222);
        const REAL_SIMD_ARRAY lambdaU_dupD00 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputlambdaU_dupD00, UpwindAlgInputlambdaU_ddnD00), UpwindAlgInputlambdaU_ddnD00);
        const REAL_SIMD_ARRAY lambdaU_dupD01 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputlambdaU_dupD01, UpwindAlgInputlambdaU_ddnD01), UpwindAlgInputlambdaU_ddnD01);
        const REAL_SIMD_ARRAY lambdaU_dupD02 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputlambdaU_dupD02, UpwindAlgInputlambdaU_ddnD02), UpwindAlgInputlambdaU_ddnD02);
        const REAL_SIMD_ARRAY lambdaU_dupD10 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputlambdaU_dupD10, UpwindAlgInputlambdaU_ddnD10), UpwindAlgInputlambdaU_ddnD10);
        const REAL_SIMD_ARRAY lambdaU_dupD11 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputlambdaU_dupD11, UpwindAlgInputlambdaU_ddnD11), UpwindAlgInputlambdaU_ddnD11);
        const REAL_SIMD_ARRAY lambdaU_dupD12 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputlambdaU_dupD12, UpwindAlgInputlambdaU_ddnD12), UpwindAlgInputlambdaU_ddnD12);
        const REAL_SIMD_ARRAY lambdaU_dupD20 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputlambdaU_dupD20, UpwindAlgInputlambdaU_ddnD20), UpwindAlgInputlambdaU_ddnD20);
        const REAL_SIMD_ARRAY lambdaU_dupD21 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputlambdaU_dupD21, UpwindAlgInputlambdaU_ddnD21), UpwindAlgInputlambdaU_ddnD21);
        const REAL_SIMD_ARRAY lambdaU_dupD22 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputlambdaU_dupD22, UpwindAlgInputlambdaU_ddnD22), UpwindAlgInputlambdaU_ddnD22);
        const REAL_SIMD_ARRAY trK_dupD0 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputtrK_dupD0, UpwindAlgInputtrK_ddnD0), UpwindAlgInputtrK_ddnD0);
        const REAL_SIMD_ARRAY trK_dupD1 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputtrK_dupD1, UpwindAlgInputtrK_ddnD1), UpwindAlgInputtrK_ddnD1);
        const REAL_SIMD_ARRAY trK_dupD2 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputtrK_dupD2, UpwindAlgInputtrK_ddnD2), UpwindAlgInputtrK_ddnD2);
        const REAL_SIMD_ARRAY vetU_dupD00 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputvetU_dupD00, UpwindAlgInputvetU_ddnD00), UpwindAlgInputvetU_ddnD00);
        const REAL_SIMD_ARRAY vetU_dupD01 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputvetU_dupD01, UpwindAlgInputvetU_ddnD01), UpwindAlgInputvetU_ddnD01);
        const REAL_SIMD_ARRAY vetU_dupD02 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputvetU_dupD02, UpwindAlgInputvetU_ddnD02), UpwindAlgInputvetU_ddnD02);
        const REAL_SIMD_ARRAY vetU_dupD10 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputvetU_dupD10, UpwindAlgInputvetU_ddnD10), UpwindAlgInputvetU_ddnD10);
        const REAL_SIMD_ARRAY vetU_dupD11 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputvetU_dupD11, UpwindAlgInputvetU_ddnD11), UpwindAlgInputvetU_ddnD11);
        const REAL_SIMD_ARRAY vetU_dupD12 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputvetU_dupD12, UpwindAlgInputvetU_ddnD12), UpwindAlgInputvetU_ddnD12);
        const REAL_SIMD_ARRAY vetU_dupD20 =
            FusedMulAddSIMD(Upwind0, SubSIMD(UpwindAlgInputvetU_dupD20, UpwindAlgInputvetU_ddnD20), UpwindAlgInputvetU_ddnD20);
        const REAL_SIMD_ARRAY vetU_dupD21 =
            FusedMulAddSIMD(Upwind1, SubSIMD(UpwindAlgInputvetU_dupD21, UpwindAlgInputvetU_ddnD21), UpwindAlgInputvetU_ddnD21);
        const REAL_SIMD_ARRAY vetU_dupD22 =
            FusedMulAddSIMD(Upwind2, SubSIMD(UpwindAlgInputvetU_dupD22, UpwindAlgInputvetU_ddnD22), UpwindAlgInputvetU_ddnD22);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 3 of 3:
         * Evaluate SymPy expressions and write to main memory.
         */
        static const double dblFDPart3_Integer_1 = 1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        static const double dblFDPart3_Integer_10 = 10.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_10 = ConstSIMD(dblFDPart3_Integer_10);

        static const double dblFDPart3_Integer_12 = 12.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_12 = ConstSIMD(dblFDPart3_Integer_12);

        static const double dblFDPart3_Integer_16 = 16.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_16 = ConstSIMD(dblFDPart3_Integer_16);

        static const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        static const double dblFDPart3_Integer_3 = 3.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_3 = ConstSIMD(dblFDPart3_Integer_3);

        static const double dblFDPart3_Integer_4 = 4.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_4 = ConstSIMD(dblFDPart3_Integer_4);

        static const double dblFDPart3_Integer_8 = 8.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_8 = ConstSIMD(dblFDPart3_Integer_8);

        static const double dblFDPart3_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        static const double dblFDPart3_Rational_1_12 = 1.0 / 12.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_12 = ConstSIMD(dblFDPart3_Rational_1_12);

        static const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        static const double dblFDPart3_Rational_1_3 = 1.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_3 = ConstSIMD(dblFDPart3_Rational_1_3);

        static const double dblFDPart3_Rational_1_4 = 1.0 / 4.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_4 = ConstSIMD(dblFDPart3_Rational_1_4);

        static const double dblFDPart3_Rational_1_6 = 1.0 / 6.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_6 = ConstSIMD(dblFDPart3_Rational_1_6);

        static const double dblFDPart3_Rational_2_3 = 2.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_2_3 = ConstSIMD(dblFDPart3_Rational_2_3);

        static const double dblFDPart3_Rational_3_2 = 3.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_3_2 = ConstSIMD(dblFDPart3_Rational_3_2);

        static const double dblFDPart3_Rational_3_4 = 3.0 / 4.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_3_4 = ConstSIMD(dblFDPart3_Rational_3_4);

        static const double dblFDPart3_Rational_4_3 = 4.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_4_3 = ConstSIMD(dblFDPart3_Rational_4_3);

        const REAL_SIMD_ARRAY FDPart3tmp0 = MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp2 = MulSIMD(FDPart3_Integer_2, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp4 = MulSIMD(alpha, trK);
        const REAL_SIMD_ARRAY FDPart3tmp5 = DivSIMD(FDPart3_Integer_1, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp8 = DivSIMD(FDPart3_Integer_1, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp12 = DivSIMD(FDPart3_Integer_1, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp14 = MulSIMD(f0_of_xx0__D0, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp15 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp31 = MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp34 = MulSIMD(f3_of_xx2, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp68 = MulSIMD(f0_of_xx0, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp76 = MulSIMD(f0_of_xx0__D0, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp78 = MulSIMD(f0_of_xx0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp106 = MulSIMD(cf, cf);
        const REAL_SIMD_ARRAY FDPart3tmp110 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf));
        const REAL_SIMD_ARRAY FDPart3tmp118 = MulSIMD(f3_of_xx2, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp122 = MulSIMD(FDPart3_Integer_2, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp123 = MulSIMD(f0_of_xx0, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp129 = MulSIMD(f3_of_xx2, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp146 = MulSIMD(f0_of_xx0__D0, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp192 = MulSIMD(cf, cf);
        const REAL_SIMD_ARRAY FDPart3tmp193 = DivSIMD(FDPart3_Integer_1, cf);
        const REAL_SIMD_ARRAY FDPart3tmp271 = MulSIMD(f0_of_xx0, f0_of_xx0__DD00);
        const REAL_SIMD_ARRAY FDPart3tmp289 = MulSIMD(aDD12, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp302 = MulSIMD(aDD02, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp363 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(T4UU00, alpha));
        const REAL_SIMD_ARRAY FDPart3tmp369 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(f3_of_xx2, f3_of_xx2), f3_of_xx2));
        const REAL_SIMD_ARRAY FDPart3tmp370 = MulSIMD(f3_of_xx2__D2, f3_of_xx2__D2);
        const REAL_SIMD_ARRAY FDPart3tmp381 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(f3_of_xx2, f3_of_xx2), f3_of_xx2), f3_of_xx2));
        const REAL_SIMD_ARRAY FDPart3tmp416 = MulSIMD(FDPart3_Integer_3, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp435 = MulSIMD(FDPart3_Rational_3_2, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp538 = MulSIMD(FDPart3_Integer_3, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp540 = MulSIMD(FDPart3_Integer_2, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp545 = MulSIMD(FDPart3_Integer_4, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp557 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf), cf), cf));
        const REAL_SIMD_ARRAY FDPart3tmp1 = DivSIMD(FDPart3_Integer_1, FDPart3tmp0);
        const REAL_SIMD_ARRAY FDPart3tmp3 = MulSIMD(FDPart3tmp0, aDD00);
        const REAL_SIMD_ARRAY FDPart3tmp6 = MulSIMD(FDPart3tmp5, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp9 = MulSIMD(FDPart3tmp8, vetU2);
        const REAL_SIMD_ARRAY FDPart3tmp11 = MulSIMD(FDPart3tmp2, f0_of_xx0__DD00);
        const REAL_SIMD_ARRAY FDPart3tmp13 = MulSIMD(FDPart3tmp12, vetU0);
        const REAL_SIMD_ARRAY FDPart3tmp16 = DivSIMD(FDPart3_Integer_1, FDPart3tmp15);
        const REAL_SIMD_ARRAY FDPart3tmp20 = MulSIMD(FDPart3tmp2, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp22 = MulSIMD(FDPart3tmp12, vetU_dD00);
        const REAL_SIMD_ARRAY FDPart3tmp29 = MulSIMD(FDPart3tmp8, vetU_dD22);
        const REAL_SIMD_ARRAY FDPart3tmp30 = MulSIMD(FDPart3tmp5, vetU_dD11);
        const REAL_SIMD_ARRAY FDPart3tmp32 = DivSIMD(FDPart3_Integer_1, FDPart3tmp31);
        const REAL_SIMD_ARRAY FDPart3tmp35 = MulSIMD(FDPart3_Integer_2, FDPart3tmp34);
        const REAL_SIMD_ARRAY FDPart3tmp36 = MulSIMD(FDPart3tmp15, FDPart3tmp34);
        const REAL_SIMD_ARRAY FDPart3tmp38 = DivSIMD(FDPart3_Integer_1, FDPart3tmp34);
        const REAL_SIMD_ARRAY FDPart3tmp44 = MulSIMD(FDPart3tmp0, FDPart3tmp15);
        const REAL_SIMD_ARRAY FDPart3tmp47 = MulSIMD(FDPart3tmp34, hDD22);
        const REAL_SIMD_ARRAY FDPart3tmp50 = MulSIMD(FDPart3tmp15, MulSIMD(hDD01, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp53 = MulSIMD(FDPart3tmp0, hDD00);
        const REAL_SIMD_ARRAY FDPart3tmp56 = MulSIMD(FDPart3tmp15, MulSIMD(hDD12, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp59 = MulSIMD(FDPart3tmp15, hDD11);
        const REAL_SIMD_ARRAY FDPart3tmp61 = MulSIMD(FDPart3tmp34, MulSIMD(hDD02, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp71 = MulSIMD(FDPart3tmp68, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp77 = MulSIMD(FDPart3tmp76, aDD02);
        const REAL_SIMD_ARRAY FDPart3tmp79 = MulSIMD(FDPart3tmp34, MulSIMD(hDD02, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp81 = MulSIMD(FDPart3tmp78, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp89 = MulSIMD(FDPart3tmp78, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp90 = MulSIMD(FDPart3tmp76, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp93 = MulSIMD(FDPart3tmp76, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp109 = FusedMulSubSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(FDPart3tmp0, hDD00),
                                                              MulSIMD(FDPart3_Rational_1_3, FDPart3tmp0));
        const REAL_SIMD_ARRAY FDPart3tmp111 = MulSIMD(FDPart3tmp110, T4UU12);
        const REAL_SIMD_ARRAY FDPart3tmp113 = MulSIMD(FDPart3tmp110, T4UU11);
        const REAL_SIMD_ARRAY FDPart3tmp115 = MulSIMD(FDPart3tmp110, T4UU22);
        const REAL_SIMD_ARRAY FDPart3tmp117 = MulSIMD(FDPart3tmp110, T4UU33);
        const REAL_SIMD_ARRAY FDPart3tmp119 = MulSIMD(FDPart3tmp110, T4UU13);
        const REAL_SIMD_ARRAY FDPart3tmp121 = MulSIMD(FDPart3tmp110, T4UU23);
        const REAL_SIMD_ARRAY FDPart3tmp124 = MulSIMD(FDPart3tmp122, FDPart3tmp123);
        const REAL_SIMD_ARRAY FDPart3tmp126 = DivSIMD(FDPart3_Integer_1, FDPart3tmp106);
        const REAL_SIMD_ARRAY FDPart3tmp141 = MulSIMD(FDPart3tmp122, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp156 = MulSIMD(FDPart3tmp2, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp194 = MulSIMD(FDPart3tmp193, cf_dD2);
        const REAL_SIMD_ARRAY FDPart3tmp197 = MulSIMD(FDPart3tmp193, cf_dD0);
        const REAL_SIMD_ARRAY FDPart3tmp200 = DivSIMD(alpha, FDPart3tmp192);
        const REAL_SIMD_ARRAY FDPart3tmp202 = MulSIMD(FDPart3tmp78, hDD_dD012);
        const REAL_SIMD_ARRAY FDPart3tmp203 = MulSIMD(FDPart3tmp76, hDD_dD021);
        const REAL_SIMD_ARRAY FDPart3tmp213 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp193);
        const REAL_SIMD_ARRAY FDPart3tmp217 = MulSIMD(FDPart3tmp193, cf_dD1);
        const REAL_SIMD_ARRAY FDPart3tmp244 = MulSIMD(FDPart3tmp15, hDD_dD112);
        const REAL_SIMD_ARRAY FDPart3tmp255 = MulSIMD(FDPart3tmp122, f3_of_xx2__D2);
        const REAL_SIMD_ARRAY FDPart3tmp292 = MulSIMD(FDPart3tmp15, aDD11);
        const REAL_SIMD_ARRAY FDPart3tmp293 = MulSIMD(FDPart3tmp68, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp299 = MulSIMD(FDPart3tmp12, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp307 = MulSIMD(FDPart3tmp34, aDD22);
        const REAL_SIMD_ARRAY FDPart3tmp312 = MulSIMD(FDPart3tmp12, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp315 = MulSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, vetU_dD01));
        const REAL_SIMD_ARRAY FDPart3tmp319 = FusedMulSubSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(FDPart3tmp15, hDD11),
                                                              MulSIMD(FDPart3_Rational_1_3, FDPart3tmp15));
        const REAL_SIMD_ARRAY FDPart3tmp322 = MulSIMD(FDPart3tmp5, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp326 = FusedMulSubSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(FDPart3tmp34, hDD22),
                                                              MulSIMD(FDPart3_Rational_1_3, FDPart3tmp34));
        const REAL_SIMD_ARRAY FDPart3tmp327 = MulSIMD(FDPart3tmp12, betU0);
        const REAL_SIMD_ARRAY FDPart3tmp335 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(alpha, trK_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp336 = MulSIMD(FDPart3tmp12, vetU_dDD002);
        const REAL_SIMD_ARRAY FDPart3tmp346 = MulSIMD(FDPart3tmp12, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp353 = MulSIMD(FDPart3tmp8, vetU_dDD212);
        const REAL_SIMD_ARRAY FDPart3tmp358 = MulSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp372 =
            FusedMulAddSIMD(MulSIMD(FDPart3_Integer_4, FDPart3tmp31), MulSIMD(FDPart3tmp68, f3_of_xx2__D2),
                            MulSIMD(FDPart3tmp76, MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3tmp15), MulSIMD(f0_of_xx0__DD00, f3_of_xx2__D2))));
        const REAL_SIMD_ARRAY FDPart3tmp380 = MulSIMD(FDPart3tmp8, vetU_dDD202);
        const REAL_SIMD_ARRAY FDPart3tmp393 = MulSIMD(FDPart3tmp8, betU2);
        const REAL_SIMD_ARRAY FDPart3tmp396 = MulSIMD(FDPart3tmp5, betU1);
        const REAL_SIMD_ARRAY FDPart3tmp402 = MulSIMD(FDPart3tmp15, hDD_dD111);
        const REAL_SIMD_ARRAY FDPart3tmp426 = MulSIMD(FDPart3tmp2, FDPart3tmp302);
        const REAL_SIMD_ARRAY FDPart3tmp431 = MulSIMD(FDPart3tmp122, FDPart3tmp289);
        const REAL_SIMD_ARRAY FDPart3tmp477 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(alpha, trK_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp488 = DivSIMD(FDPart3tmp0, MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp33 = MulSIMD(FDPart3tmp32, vetU0);
        const REAL_SIMD_ARRAY FDPart3tmp37 = FusedMulAddSIMD(FDPart3tmp11, FDPart3tmp36, MulSIMD(FDPart3tmp31, MulSIMD(FDPart3tmp35, f0_of_xx0)));
        const REAL_SIMD_ARRAY FDPart3tmp46 = MulSIMD(hDD01, MulSIMD(MulSIMD(FDPart3tmp35, FDPart3tmp44), MulSIMD(hDD02, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp48 = AddSIMD(FDPart3tmp34, FDPart3tmp47);
        const REAL_SIMD_ARRAY FDPart3tmp51 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp0, FDPart3tmp50));
        const REAL_SIMD_ARRAY FDPart3tmp54 = AddSIMD(FDPart3tmp0, FDPart3tmp53);
        const REAL_SIMD_ARRAY FDPart3tmp57 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp34, FDPart3tmp56));
        const REAL_SIMD_ARRAY FDPart3tmp60 = AddSIMD(FDPart3tmp15, FDPart3tmp59);
        const REAL_SIMD_ARRAY FDPart3tmp62 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp0, FDPart3tmp61));
        const REAL_SIMD_ARRAY FDPart3tmp70 = MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp68), MulSIMD(hDD01, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp80 = MulSIMD(FDPart3tmp78, FDPart3tmp79);
        const REAL_SIMD_ARRAY FDPart3tmp92 = MulSIMD(FDPart3tmp15, MulSIMD(FDPart3tmp90, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp114 = MulSIMD(FDPart3tmp0, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp116 = MulSIMD(FDPart3tmp34, FDPart3tmp56);
        const REAL_SIMD_ARRAY FDPart3tmp120 = MulSIMD(FDPart3tmp118, MulSIMD(FDPart3tmp119, FDPart3tmp2));
        const REAL_SIMD_ARRAY FDPart3tmp137 = MulSIMD(FDPart3tmp126, T4UU03);
        const REAL_SIMD_ARRAY FDPart3tmp142 = MulSIMD(FDPart3tmp0, FDPart3tmp50);
        const REAL_SIMD_ARRAY FDPart3tmp144 = MulSIMD(FDPart3tmp126, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp205 = FusedMulAddSIMD(FDPart3tmp68, hDD_dD120, FDPart3tmp90);
        const REAL_SIMD_ARRAY FDPart3tmp216 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp194, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp221 = AddSIMD(FDPart3tmp203, SubSIMD(SubSIMD(FDPart3tmp202, FDPart3tmp90), MulSIMD(FDPart3tmp68, hDD_dD120)));
        const REAL_SIMD_ARRAY FDPart3tmp230 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp197, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp232 = FusedMulAddSIMD(FDPart3tmp20, hDD11, FDPart3tmp20);
        const REAL_SIMD_ARRAY FDPart3tmp245 = FusedMulSubSIMD(FDPart3tmp122, MulSIMD(f0_of_xx0, hDD_dD121), FDPart3tmp244);
        const REAL_SIMD_ARRAY FDPart3tmp246 =
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, f0_of_xx0__D0),
                               FusedMulAddSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD11)),
                                               FusedMulSubSIMD(FDPart3tmp20, hDD_dD011, MulSIMD(FDPart3tmp15, hDD_dD110))));
        const REAL_SIMD_ARRAY FDPart3tmp256 = FusedMulAddSIMD(FDPart3tmp255, hDD22, FDPart3tmp255);
        const REAL_SIMD_ARRAY FDPart3tmp259 =
            FusedMulAddSIMD(FDPart3tmp122, MulSIMD(f0_of_xx0, hDD_dD122),
                            FusedMulSubSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp123, f3_of_xx2__D2), MulSIMD(FDPart3tmp34, hDD_dD221)));
        const REAL_SIMD_ARRAY FDPart3tmp261 =
            FusedMulAddSIMD(FDPart3tmp2, MulSIMD(f3_of_xx2__D2, hDD02), FusedMulSubSIMD(FDPart3tmp156, hDD_dD022, MulSIMD(FDPart3tmp34, hDD_dD220)));
        const REAL_SIMD_ARRAY FDPart3tmp270 =
            FusedMulAddSIMD(FDPart3tmp156, hDD_dD020, FusedMulSubSIMD(FDPart3tmp141, f0_of_xx0__DD00, MulSIMD(FDPart3tmp0, hDD_dD002)));
        const REAL_SIMD_ARRAY FDPart3tmp273 = MulSIMD(hDD01, AddSIMD(FDPart3tmp0, FDPart3tmp271));
        const REAL_SIMD_ARRAY FDPart3tmp275 = FusedMulAddSIMD(FDPart3tmp11, hDD00, FDPart3tmp11);
        const REAL_SIMD_ARRAY FDPart3tmp329 = MulSIMD(FDPart3tmp12, FDPart3tmp6);
        const REAL_SIMD_ARRAY FDPart3tmp330 = MulSIMD(FDPart3tmp12, FDPart3tmp9);
        const REAL_SIMD_ARRAY FDPart3tmp339 = FusedMulSubSIMD(FDPart3tmp12, vetU_dDD012, MulSIMD(FDPart3tmp12, vetU_dD12));
        const REAL_SIMD_ARRAY FDPart3tmp392 = MulSIMD(FDPart3tmp327, FDPart3tmp9);
        const REAL_SIMD_ARRAY FDPart3tmp394 = MulSIMD(FDPart3tmp13, FDPart3tmp393);
        const REAL_SIMD_ARRAY FDPart3tmp397 = MulSIMD(FDPart3tmp396, FDPart3tmp9);
        const REAL_SIMD_ARRAY FDPart3tmp398 = MulSIMD(FDPart3tmp393, FDPart3tmp6);
        const REAL_SIMD_ARRAY FDPart3tmp400 = MulSIMD(FDPart3tmp327, FDPart3tmp6);
        const REAL_SIMD_ARRAY FDPart3tmp401 = MulSIMD(FDPart3tmp13, FDPart3tmp396);
        const REAL_SIMD_ARRAY FDPart3tmp404 = MulSIMD(FDPart3tmp16, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp407 = MulSIMD(FDPart3tmp38, MulSIMD(betU2, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp409 = MulSIMD(FDPart3tmp1, MulSIMD(betU0, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp428 = MulSIMD(FDPart3tmp20, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp442 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp12, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp475 = MulSIMD(FDPart3tmp5, FDPart3tmp9);
        const REAL_SIMD_ARRAY FDPart3tmp478 = FusedMulAddSIMD(FDPart3tmp5, vetU_dD02, MulSIMD(FDPart3tmp5, vetU_dDD112));
        const REAL_SIMD_ARRAY FDPart3tmp485 = FusedMulSubSIMD(FDPart3tmp5, f0_of_xx0__DD00, MulSIMD(FDPart3tmp0, FDPart3tmp16));
        const REAL_SIMD_ARRAY FDPart3tmp505 = MulSIMD(FDPart3tmp13, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp506 = MulSIMD(FDPart3tmp6, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp507 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp322, vetU_dD21));
        const REAL_SIMD_ARRAY FDPart3tmp532 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp12, lambdaU0));
        const REAL_SIMD_ARRAY FDPart3tmp539 = FusedMulAddSIMD(FDPart3tmp197, FDPart3tmp538, alpha_dD0);
        const REAL_SIMD_ARRAY FDPart3tmp543 = FusedMulAddSIMD(FDPart3tmp217, FDPart3tmp538, alpha_dD1);
        const REAL_SIMD_ARRAY FDPart3tmp544 = FusedMulAddSIMD(FDPart3tmp194, FDPart3tmp538, alpha_dD2);
        const REAL_SIMD_ARRAY FDPart3tmp562 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp6, FDPart3tmp9));
        const REAL_SIMD_ARRAY FDPart3tmp564 = MulSIMD(FDPart3tmp16, MulSIMD(vetU1, vetU1));
        const REAL_SIMD_ARRAY FDPart3tmp565 = MulSIMD(FDPart3tmp38, MulSIMD(vetU2, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp566 = MulSIMD(FDPart3tmp1, MulSIMD(vetU0, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp18 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp14, FDPart3tmp16));
        const REAL_SIMD_ARRAY FDPart3tmp25 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp1), MulSIMD(f0_of_xx0__DD00, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp40 = MulSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp16, FDPart3tmp38));
        const REAL_SIMD_ARRAY FDPart3tmp64 = MulSIMD(FDPart3tmp48, FDPart3tmp60);
        const REAL_SIMD_ARRAY FDPart3tmp72 = MulSIMD(FDPart3tmp54, FDPart3tmp71);
        const REAL_SIMD_ARRAY FDPart3tmp82 = MulSIMD(FDPart3tmp48, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp94 = MulSIMD(FDPart3tmp60, FDPart3tmp93);
        const REAL_SIMD_ARRAY FDPart3tmp128 = MulSIMD(FDPart3tmp118, MulSIMD(FDPart3tmp126, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp131 = MulSIMD(FDPart3tmp126, MulSIMD(FDPart3tmp129, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp135 = MulSIMD(FDPart3tmp118, MulSIMD(FDPart3tmp126, FDPart3tmp2));
        const REAL_SIMD_ARRAY FDPart3tmp143 = MulSIMD(FDPart3tmp111, MulSIMD(FDPart3tmp20, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp148 = MulSIMD(FDPart3tmp126, MulSIMD(FDPart3tmp146, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp149 = MulSIMD(FDPart3tmp126, FDPart3tmp54);
        const REAL_SIMD_ARRAY FDPart3tmp157 = MulSIMD(FDPart3tmp144, MulSIMD(f0_of_xx0, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp158 = MulSIMD(FDPart3tmp123, MulSIMD(FDPart3tmp126, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp168 = MulSIMD(FDPart3tmp144, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp171 = MulSIMD(FDPart3tmp126, MulSIMD(FDPart3tmp68, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp196 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp194, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp199 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp197, alpha_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp206 = AddSIMD(FDPart3tmp205, SubSIMD(FDPart3tmp202, FDPart3tmp203));
        const REAL_SIMD_ARRAY FDPart3tmp219 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp217, alpha_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp229 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp217, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp233 = FusedMulAddSIMD(FDPart3tmp15, hDD_dD110, FDPart3tmp232);
        const REAL_SIMD_ARRAY FDPart3tmp234 = AddSIMD(FDPart3tmp205, SubSIMD(FDPart3tmp203, FDPart3tmp202));
        const REAL_SIMD_ARRAY FDPart3tmp257 = FusedMulAddSIMD(FDPart3tmp34, hDD_dD222, FDPart3tmp256);
        const REAL_SIMD_ARRAY FDPart3tmp274 =
            FusedMulAddSIMD(FDPart3tmp20, hDD_dD010, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp273, MulSIMD(FDPart3tmp0, hDD_dD001)));
        const REAL_SIMD_ARRAY FDPart3tmp276 = FusedMulAddSIMD(FDPart3tmp0, hDD_dD000, FDPart3tmp275);
        const REAL_SIMD_ARRAY FDPart3tmp304 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp38, f3_of_xx2__D2));
        const REAL_SIMD_ARRAY FDPart3tmp332 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp1, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp359 = FusedMulAddSIMD(
            vetU0, FusedMulSubSIMD(FDPart3tmp32, FDPart3tmp358, MulSIMD(FDPart3tmp1, f0_of_xx0__DDD000)),
            FusedMulSubSIMD(FDPart3tmp12, vetU_dDD000, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp1), MulSIMD(f0_of_xx0__DD00, vetU_dD00))));
        const REAL_SIMD_ARRAY FDPart3tmp377 = FusedMulAddSIMD(
            vetU2, FusedMulSubSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp369, FDPart3tmp370), MulSIMD(FDPart3tmp38, f3_of_xx2__DD22)),
            FusedMulSubSIMD(FDPart3tmp8, vetU_dDD222, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp38), MulSIMD(f3_of_xx2__D2, vetU_dD22))));
        const REAL_SIMD_ARRAY FDPart3tmp382 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp16, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp405 = MulSIMD(FDPart3tmp404, betU1);
        const REAL_SIMD_ARRAY FDPart3tmp508 = FusedMulAddSIMD(FDPart3tmp507, f0_of_xx0__D0, MulSIMD(FDPart3tmp8, vetU_dDD201));
        const REAL_SIMD_ARRAY FDPart3tmp561 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp13, FDPart3tmp9));
        const REAL_SIMD_ARRAY FDPart3tmp563 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp13, FDPart3tmp6));
        const REAL_SIMD_ARRAY FDPart3tmp19 = FusedMulAddSIMD(FDPart3tmp5, vetU_dD10, FDPart3tmp18);
        const REAL_SIMD_ARRAY FDPart3tmp26 = AddSIMD(FDPart3tmp22, FDPart3tmp25);
        const REAL_SIMD_ARRAY FDPart3tmp87 = FusedMulAddSIMD(FDPart3tmp48, FDPart3tmp54, FDPart3tmp62);
        const REAL_SIMD_ARRAY FDPart3tmp99 = FusedMulAddSIMD(FDPart3tmp54, FDPart3tmp60, FDPart3tmp51);
        const REAL_SIMD_ARRAY FDPart3tmp103 = AddSIMD(FDPart3tmp57, FDPart3tmp64);
        const REAL_SIMD_ARRAY FDPart3tmp133 = MulSIMD(FDPart3tmp126, MulSIMD(FDPart3tmp48, FDPart3tmp9));
        const REAL_SIMD_ARRAY FDPart3tmp150 =
            FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp149, FusedMulAddSIMD(FDPart3tmp14, FDPart3tmp144, FDPart3tmp148));
        const REAL_SIMD_ARRAY FDPart3tmp153 = MulSIMD(FDPart3tmp149, T4UU01);
        const REAL_SIMD_ARRAY FDPart3tmp160 = MulSIMD(FDPart3tmp126, MulSIMD(FDPart3tmp6, FDPart3tmp60));
        const REAL_SIMD_ARRAY FDPart3tmp162 = MulSIMD(FDPart3tmp126, MulSIMD(FDPart3tmp60, T4UU02));
        const REAL_SIMD_ARRAY FDPart3tmp167 = MulSIMD(FDPart3tmp110, MulSIMD(FDPart3tmp54, FDPart3tmp60));
        const REAL_SIMD_ARRAY FDPart3tmp182 = MulSIMD(FDPart3tmp126, MulSIMD(FDPart3tmp48, T4UU03));
        const REAL_SIMD_ARRAY FDPart3tmp187 = MulSIMD(FDPart3tmp110, MulSIMD(FDPart3tmp48, FDPart3tmp54));
        const REAL_SIMD_ARRAY FDPart3tmp374 = MulSIMD(FDPart3tmp1, FDPart3tmp40);
        const REAL_SIMD_ARRAY FDPart3tmp383 = FusedMulAddSIMD(FDPart3tmp382, vetU_dD11, MulSIMD(FDPart3tmp5, vetU_dDD101));
        const REAL_SIMD_ARRAY FDPart3tmp42 =
            AddSIMD(AddSIMD(FDPart3tmp26, FDPart3tmp29), FusedMulAddSIMD(FDPart3tmp33, MulSIMD(FDPart3tmp37, FDPart3tmp40), FDPart3tmp30));
        const REAL_SIMD_ARRAY FDPart3tmp66 =
            FusedMulAddSIMD(FDPart3tmp54, FDPart3tmp57,
                            FusedMulAddSIMD(FDPart3tmp54, FDPart3tmp64,
                                            FusedMulAddSIMD(FDPart3tmp60, FDPart3tmp62, FusedMulAddSIMD(FDPart3tmp48, FDPart3tmp51, FDPart3tmp46))));
        const REAL_SIMD_ARRAY FDPart3tmp74 = SubSIMD(FDPart3tmp70, FDPart3tmp72);
        const REAL_SIMD_ARRAY FDPart3tmp84 = SubSIMD(FDPart3tmp80, FDPart3tmp82);
        const REAL_SIMD_ARRAY FDPart3tmp96 = SubSIMD(FDPart3tmp92, FDPart3tmp94);
        const REAL_SIMD_ARRAY FDPart3tmp134 = AddSIMD(FDPart3tmp128, AddSIMD(FDPart3tmp131, FDPart3tmp133));
        const REAL_SIMD_ARRAY FDPart3tmp161 = AddSIMD(FDPart3tmp157, AddSIMD(FDPart3tmp158, FDPart3tmp160));
        const REAL_SIMD_ARRAY FDPart3tmp306 = FusedMulAddSIMD(FDPart3tmp304, vetU2, FDPart3tmp29);
        const REAL_SIMD_ARRAY FDPart3tmp354 =
            FusedMulAddSIMD(FDPart3tmp332, vetU_dD01,
                            FusedMulAddSIMD(FDPart3tmp5, vetU_dDD111,
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp32, FDPart3tmp37), MulSIMD(FDPart3tmp40, vetU_dD01),
                                                            FusedMulAddSIMD(FDPart3tmp12, vetU_dDD001, FDPart3tmp353))));
        const REAL_SIMD_ARRAY FDPart3tmp364 = FusedMulAddSIMD(
            FDPart3tmp168, MulSIMD(T4UU01, alpha),
            FusedMulAddSIMD(FDPart3tmp171, MulSIMD(T4UU03, alpha),
                            FusedMulAddSIMD(FDPart3tmp162, alpha,
                                            MulSIMD(FDPart3tmp363,
                                                    FusedMulSubSIMD(FDPart3_NegativeOne_, AddSIMD(FDPart3tmp158, FDPart3tmp157), FDPart3tmp160)))));
        const REAL_SIMD_ARRAY FDPart3tmp366 = FusedMulAddSIMD(
            FDPart3tmp171, MulSIMD(T4UU02, alpha),
            FusedMulAddSIMD(MulSIMD(FDPart3tmp126, FDPart3tmp93), MulSIMD(T4UU01, alpha),
                            FusedMulAddSIMD(FDPart3tmp182, alpha,
                                            MulSIMD(FDPart3tmp363,
                                                    FusedMulSubSIMD(FDPart3_NegativeOne_, AddSIMD(FDPart3tmp131, FDPart3tmp128), FDPart3tmp133)))));
        const REAL_SIMD_ARRAY FDPart3tmp367 = FusedMulAddSIMD(
            FDPart3tmp137, MulSIMD(FDPart3tmp93, alpha),
            FusedMulAddSIMD(FDPart3tmp168, MulSIMD(T4UU02, alpha),
                            FusedMulAddSIMD(FDPart3tmp153, alpha,
                                            MulSIMD(FDPart3tmp363, SubSIMD(NegFusedMulSubSIMD(FDPart3tmp13, FDPart3tmp149, FDPart3tmp148),
                                                                           MulSIMD(FDPart3tmp14, FDPart3tmp144))))));
        const REAL_SIMD_ARRAY FDPart3tmp384 = FusedMulAddSIMD(
            FDPart3tmp374,
            FusedMulAddSIMD(FDPart3tmp372, FDPart3tmp9, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp44), MulSIMD(f3_of_xx2__D2, vetU_dD20))),
            FusedMulAddSIMD(
                FDPart3tmp304, vetU_dD20,
                FusedMulAddSIMD(
                    FDPart3tmp374,
                    FusedMulAddSIMD(
                        FDPart3tmp13,
                        FusedMulAddSIMD(
                            FDPart3tmp2, MulSIMD(FDPart3tmp36, f0_of_xx0__DDD000),
                            FusedMulAddSIMD(MulSIMD(FDPart3_Integer_10, FDPart3tmp0), MulSIMD(FDPart3tmp271, FDPart3tmp34),
                                            FusedMulAddSIMD(FDPart3tmp35,
                                                            MulSIMD(MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0), f0_of_xx0__D0),
                                                            MulSIMD(FDPart3tmp358, FDPart3tmp36)))),
                        MulSIMD(FDPart3tmp26, FDPart3tmp37)),
                    AddSIMD(AddSIMD(FDPart3tmp380, FDPart3tmp383),
                            FusedMulAddSIMD(
                                FDPart3tmp381,
                                MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp37),
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp1), MulSIMD(f3_of_xx2__D2, vetU2))),
                                NegFusedMulAddSIMD(
                                    DivSIMD(vetU0, MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0)),
                                    MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp381),
                                            DivSIMD(MulSIMD(FDPart3tmp37, FDPart3tmp37),
                                                    MulSIMD(MulSIMD(MulSIMD(MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0), f0_of_xx0__D0), f0_of_xx0__D0),
                                                            f0_of_xx0__D0))),
                                    FDPart3tmp359))))));
        const REAL_SIMD_ARRAY FDPart3tmp483 = FusedMulAddSIMD(FDPart3tmp26, MulSIMD(FDPart3tmp5, f0_of_xx0__D0), FDPart3tmp383);
        const REAL_SIMD_ARRAY FDPart3tmp558 =
            DivSIMD(FDPart3_Integer_1,
                    FusedMulAddSIMD(FDPart3tmp54, MulSIMD(FDPart3tmp557, FDPart3tmp57),
                                    FusedMulAddSIMD(FDPart3tmp54, MulSIMD(FDPart3tmp557, FDPart3tmp64),
                                                    FusedMulAddSIMD(FDPart3tmp557, MulSIMD(FDPart3tmp60, FDPart3tmp62),
                                                                    FusedMulAddSIMD(FDPart3tmp46, FDPart3tmp557,
                                                                                    MulSIMD(FDPart3tmp48, MulSIMD(FDPart3tmp51, FDPart3tmp557)))))));
        const REAL_SIMD_ARRAY FDPart3tmp67 = DivSIMD(FDPart3_Integer_1, FDPart3tmp66);
        const REAL_SIMD_ARRAY FDPart3tmp154 = FusedMulAddSIMD(
            FDPart3_Integer_2, MulSIMD(FDPart3tmp150, FDPart3tmp153),
            FusedMulAddSIMD(
                FDPart3tmp135, MulSIMD(FDPart3tmp150, T4UU03),
                FusedMulAddSIMD(
                    FDPart3tmp143, FDPart3tmp54,
                    FusedMulAddSIMD(
                        T4UU00, MulSIMD(FDPart3tmp150, FDPart3tmp150),
                        FusedMulAddSIMD(
                            FDPart3tmp115, FDPart3tmp142,
                            FusedMulAddSIMD(FDPart3tmp120, FDPart3tmp54,
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp144, FDPart3tmp150), MulSIMD(FDPart3tmp20, T4UU02),
                                                            FusedMulAddSIMD(FDPart3tmp141,
                                                                            MulSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp121), MulSIMD(f0_of_xx0, hDD01)),
                                                                            FusedMulAddSIMD(FDPart3tmp113, MulSIMD(FDPart3tmp54, FDPart3tmp54),
                                                                                            MulSIMD(FDPart3tmp114, FDPart3tmp117))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp163 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(FDPart3tmp124, FDPart3tmp60),
            FusedMulAddSIMD(
                FDPart3tmp124, MulSIMD(FDPart3tmp137, FDPart3tmp161),
                FusedMulAddSIMD(
                    T4UU00, MulSIMD(FDPart3tmp161, FDPart3tmp161),
                    FusedMulAddSIMD(
                        FDPart3_Integer_2, MulSIMD(FDPart3tmp161, FDPart3tmp162),
                        FusedMulAddSIMD(
                            FDPart3tmp116, FDPart3tmp117,
                            FusedMulAddSIMD(
                                FDPart3tmp143, FDPart3tmp60,
                                FusedMulAddSIMD(MulSIMD(FDPart3tmp144, FDPart3tmp161), MulSIMD(FDPart3tmp20, T4UU01),
                                                FusedMulAddSIMD(FDPart3tmp156, MulSIMD(MulSIMD(FDPart3tmp119, FDPart3tmp15), MulSIMD(hDD01, hDD12)),
                                                                FusedMulAddSIMD(FDPart3tmp113, FDPart3tmp142,
                                                                                MulSIMD(FDPart3tmp115, MulSIMD(FDPart3tmp60, FDPart3tmp60)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp348 =
            SubSIMD(FusedMulSubSIMD(FDPart3tmp12, vetU_dDD001, MulSIMD(FDPart3tmp19, FDPart3tmp346)), MulSIMD(FDPart3tmp5, vetU_dD01));
        const REAL_SIMD_ARRAY FDPart3tmp355 = MulSIMD(FDPart3_Rational_1_4, FDPart3tmp354);
        const REAL_SIMD_ARRAY FDPart3tmp378 = FusedMulAddSIMD(
            FDPart3tmp5, vetU_dDD112,
            FusedMulAddSIMD(
                FDPart3tmp374, FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp372, MulSIMD(FDPart3tmp12, MulSIMD(FDPart3tmp37, vetU_dD02))),
                FusedMulAddSIMD(FDPart3tmp374,
                                FusedMulAddSIMD(FDPart3tmp9,
                                                FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp370, FDPart3tmp44),
                                                                MulSIMD(FDPart3tmp122, MulSIMD(FDPart3tmp44, f3_of_xx2__DD22))),
                                                MulSIMD(FDPart3tmp255, MulSIMD(FDPart3tmp306, FDPart3tmp44))),
                                AddSIMD(FusedMulAddSIMD(FDPart3tmp332, vetU_dD02, FDPart3tmp377),
                                        FusedMulAddSIMD(MulSIMD(FDPart3tmp33, FDPart3tmp369),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp16), MulSIMD(FDPart3tmp37, f3_of_xx2__D2)),
                                                        NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp369), MulSIMD(FDPart3tmp370, vetU2),
                                                                           FDPart3tmp336))))));
        const REAL_SIMD_ARRAY FDPart3tmp385 = MulSIMD(FDPart3_Rational_1_4, FDPart3tmp384);
        const REAL_SIMD_ARRAY FDPart3tmp410 = MulSIMD(FDPart3tmp84, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp411 = DivSIMD(FDPart3_Integer_1, MulSIMD(FDPart3tmp66, FDPart3tmp66));
        const REAL_SIMD_ARRAY FDPart3tmp412 = MulSIMD(FDPart3tmp96, FDPart3tmp96);
        const REAL_SIMD_ARRAY FDPart3tmp421 = MulSIMD(FDPart3tmp74, FDPart3tmp74);
        const REAL_SIMD_ARRAY FDPart3tmp425 = MulSIMD(FDPart3tmp74, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp534 = MulSIMD(FDPart3_Rational_1_3, FDPart3tmp354);
        const REAL_SIMD_ARRAY FDPart3tmp537 = MulSIMD(FDPart3_Rational_1_3, FDPart3tmp384);
        const REAL_SIMD_ARRAY FDPart3tmp559 = MulSIMD(FDPart3_Integer_2, FDPart3tmp558);
        const REAL_SIMD_ARRAY FDPart3tmp75 = MulSIMD(FDPart3tmp67, FDPart3tmp74);
        const REAL_SIMD_ARRAY FDPart3tmp85 = MulSIMD(FDPart3tmp67, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp88 = MulSIMD(FDPart3tmp67, FDPart3tmp87);
        const REAL_SIMD_ARRAY FDPart3tmp97 = MulSIMD(FDPart3tmp67, FDPart3tmp96);
        const REAL_SIMD_ARRAY FDPart3tmp100 = MulSIMD(FDPart3tmp67, FDPart3tmp99);
        const REAL_SIMD_ARRAY FDPart3tmp104 = MulSIMD(FDPart3tmp103, FDPart3tmp67);
        const REAL_SIMD_ARRAY FDPart3tmp139 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(FDPart3tmp124, FDPart3tmp48),
            FusedMulAddSIMD(
                FDPart3tmp134, MulSIMD(FDPart3tmp135, T4UU01),
                FusedMulAddSIMD(
                    T4UU00, MulSIMD(FDPart3tmp134, FDPart3tmp134),
                    FusedMulAddSIMD(
                        FDPart3tmp111, MulSIMD(FDPart3tmp20, FDPart3tmp79),
                        FusedMulAddSIMD(
                            FDPart3tmp117, MulSIMD(FDPart3tmp48, FDPart3tmp48),
                            FusedMulAddSIMD(FDPart3tmp120, FDPart3tmp48,
                                            FusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp134), MulSIMD(FDPart3tmp137, FDPart3tmp48),
                                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp124, FDPart3tmp126), MulSIMD(FDPart3tmp134, T4UU02),
                                                                            FusedMulAddSIMD(FDPart3tmp113, FDPart3tmp114,
                                                                                            MulSIMD(FDPart3tmp115, FDPart3tmp116))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp174 = FusedMulAddSIMD(
            FDPart3tmp150, MulSIMD(FDPart3tmp161, T4UU00),
            FusedMulAddSIMD(
                FDPart3tmp150, MulSIMD(FDPart3tmp168, T4UU01),
                FusedMulAddSIMD(
                    FDPart3tmp115, MulSIMD(FDPart3tmp60, FDPart3tmp81),
                    FusedMulAddSIMD(
                        FDPart3tmp137, MulSIMD(FDPart3tmp161, FDPart3tmp93),
                        FusedMulAddSIMD(
                            FDPart3tmp167, T4UU12,
                            FusedMulAddSIMD(
                                FDPart3tmp113, MulSIMD(FDPart3tmp54, FDPart3tmp81),
                                FusedMulAddSIMD(
                                    FDPart3tmp150, FDPart3tmp162,
                                    FusedMulAddSIMD(
                                        FDPart3tmp153, FDPart3tmp161,
                                        FusedMulAddSIMD(
                                            FDPart3tmp121, FDPart3tmp92,
                                            FusedMulAddSIMD(
                                                FDPart3tmp121, FDPart3tmp94,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp119, FDPart3tmp70,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp119, FDPart3tmp72,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp150, MulSIMD(FDPart3tmp171, T4UU03),
                                                            FusedMulAddSIMD(FDPart3tmp161, MulSIMD(FDPart3tmp168, T4UU02),
                                                                            FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp142,
                                                                                            MulSIMD(FDPart3tmp117, FDPart3tmp80))))))))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp183 = FusedMulAddSIMD(
            FDPart3tmp134, MulSIMD(FDPart3tmp168, T4UU01),
            FusedMulAddSIMD(
                FDPart3tmp134, MulSIMD(FDPart3tmp171, T4UU03),
                FusedMulAddSIMD(
                    FDPart3tmp117, MulSIMD(FDPart3tmp48, FDPart3tmp71),
                    FusedMulAddSIMD(
                        FDPart3tmp134, MulSIMD(FDPart3tmp161, T4UU00),
                        FusedMulAddSIMD(
                            FDPart3tmp110, MulSIMD(FDPart3tmp64, T4UU23),
                            FusedMulAddSIMD(
                                FDPart3tmp115, MulSIMD(FDPart3tmp60, FDPart3tmp71),
                                FusedMulAddSIMD(
                                    FDPart3tmp134, FDPart3tmp162,
                                    FusedMulAddSIMD(
                                        FDPart3tmp161, FDPart3tmp182,
                                        FusedMulAddSIMD(
                                            FDPart3tmp119, FDPart3tmp80,
                                            FusedMulAddSIMD(
                                                FDPart3tmp119, FDPart3tmp82,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp113, FDPart3tmp70,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp116, FDPart3tmp121,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp161, MulSIMD(FDPart3tmp171, T4UU02),
                                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp126, FDPart3tmp161), MulSIMD(FDPart3tmp93, T4UU01),
                                                                            FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp92,
                                                                                            MulSIMD(FDPart3tmp111, FDPart3tmp94))))))))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp189 = FusedMulAddSIMD(
            FDPart3tmp134, MulSIMD(FDPart3tmp150, T4UU00),
            FusedMulAddSIMD(
                FDPart3tmp134, MulSIMD(FDPart3tmp168, T4UU02),
                FusedMulAddSIMD(
                    FDPart3tmp117, MulSIMD(FDPart3tmp48, FDPart3tmp93),
                    FusedMulAddSIMD(
                        FDPart3tmp134, MulSIMD(FDPart3tmp137, FDPart3tmp93),
                        FusedMulAddSIMD(
                            FDPart3tmp187, T4UU13,
                            FusedMulAddSIMD(
                                FDPart3tmp113, MulSIMD(FDPart3tmp54, FDPart3tmp93),
                                FusedMulAddSIMD(
                                    FDPart3tmp134, FDPart3tmp153,
                                    FusedMulAddSIMD(
                                        FDPart3tmp150, FDPart3tmp182,
                                        FusedMulAddSIMD(
                                            FDPart3tmp121, FDPart3tmp80,
                                            FusedMulAddSIMD(
                                                FDPart3tmp121, FDPart3tmp82,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp114, FDPart3tmp119,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp115, FDPart3tmp92,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp150, MulSIMD(FDPart3tmp171, T4UU02),
                                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp126, FDPart3tmp150), MulSIMD(FDPart3tmp93, T4UU01),
                                                                            FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp70,
                                                                                            MulSIMD(FDPart3tmp111, FDPart3tmp72))))))))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp207 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp34, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp96, hDD_dD220)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp103),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, hDD_dD002)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp206), MulSIMD(FDPart3tmp67, FDPart3tmp84))));
        const REAL_SIMD_ARRAY FDPart3tmp209 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp34, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp74, hDD_dD220)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp67),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp84, hDD_dD002)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp206), MulSIMD(FDPart3tmp67, FDPart3tmp87))));
        const REAL_SIMD_ARRAY FDPart3tmp211 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp34, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp99, hDD_dD220)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp0, FDPart3tmp67),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp96, hDD_dD002)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp206), MulSIMD(FDPart3tmp67, FDPart3tmp74))));
        const REAL_SIMD_ARRAY FDPart3tmp222 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp34, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp96, hDD_dD221)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp67),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp84, hDD_dD112)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp103), MulSIMD(FDPart3tmp221, FDPart3tmp67))));
        const REAL_SIMD_ARRAY FDPart3tmp224 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp34, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp74, hDD_dD221)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp67),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp87, hDD_dD112)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp221), MulSIMD(FDPart3tmp67, FDPart3tmp84))));
        const REAL_SIMD_ARRAY FDPart3tmp226 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp34, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp99, hDD_dD221)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp67),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp74, hDD_dD112)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp221), MulSIMD(FDPart3tmp67, FDPart3tmp96))));
        const REAL_SIMD_ARRAY FDPart3tmp235 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp0, FDPart3tmp103), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp234, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp96)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp233), MulSIMD(FDPart3tmp67, FDPart3tmp84))));
        const REAL_SIMD_ARRAY FDPart3tmp237 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp0, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp84, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp234, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp74)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp233), MulSIMD(FDPart3tmp67, FDPart3tmp87))));
        const REAL_SIMD_ARRAY FDPart3tmp239 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp0, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp96, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp234, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp99)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp233), MulSIMD(FDPart3tmp67, FDPart3tmp74))));
        const REAL_SIMD_ARRAY FDPart3tmp247 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp15, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp84, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp245, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp96)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp103), MulSIMD(FDPart3tmp246, FDPart3tmp67))));
        const REAL_SIMD_ARRAY FDPart3tmp249 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp15, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp87, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp246, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp84)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp245), MulSIMD(FDPart3tmp67, FDPart3tmp74))));
        const REAL_SIMD_ARRAY FDPart3tmp251 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp15, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp74, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp246, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp96)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp245), MulSIMD(FDPart3tmp67, FDPart3tmp99))));
        const REAL_SIMD_ARRAY FDPart3tmp262 = FusedMulAddSIMD(
            FDPart3tmp259, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp84)),
            FusedMulSubSIMD(FDPart3tmp257, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp96)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp103), MulSIMD(FDPart3tmp261, FDPart3tmp67))));
        const REAL_SIMD_ARRAY FDPart3tmp264 = FusedMulAddSIMD(
            FDPart3tmp261, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp84)),
            FusedMulSubSIMD(FDPart3tmp259, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp87)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp257), MulSIMD(FDPart3tmp67, FDPart3tmp74))));
        const REAL_SIMD_ARRAY FDPart3tmp266 = FusedMulAddSIMD(
            FDPart3tmp261, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp96)),
            FusedMulSubSIMD(FDPart3tmp259, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp74)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp257), MulSIMD(FDPart3tmp67, FDPart3tmp99))));
        const REAL_SIMD_ARRAY FDPart3tmp277 = FusedMulAddSIMD(
            FDPart3tmp274, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp84)),
            FusedMulSubSIMD(FDPart3tmp270, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp96)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp103), MulSIMD(FDPart3tmp276, FDPart3tmp67))));
        const REAL_SIMD_ARRAY FDPart3tmp279 = FusedMulAddSIMD(
            FDPart3tmp276, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp84)),
            FusedMulSubSIMD(FDPart3tmp274, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp87)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp270), MulSIMD(FDPart3tmp67, FDPart3tmp74))));
        const REAL_SIMD_ARRAY FDPart3tmp281 = FusedMulAddSIMD(
            FDPart3tmp276, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp96)),
            FusedMulSubSIMD(FDPart3tmp274, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp67, FDPart3tmp74)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp270), MulSIMD(FDPart3tmp67, FDPart3tmp99))));
        const REAL_SIMD_ARRAY FDPart3tmp379 = MulSIMD(FDPart3_Rational_1_4, FDPart3tmp378);
        const REAL_SIMD_ARRAY FDPart3tmp415 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp96, aDD02),
            MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp411),
                    MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp34, FDPart3tmp411), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp412, aDD22)),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp84, aDD01),
                    MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp411),
                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                    FusedMulAddSIMD(MulSIMD(FDPart3tmp96, aDD12),
                                    MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp84),
                                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
                                    FusedMulSubSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp410),
                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp411, aDD11)),
                                                    MulSIMD(FDPart3tmp411, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0),
                                                                                   MulSIMD(aDD00, MulSIMD(FDPart3tmp103, FDPart3tmp103)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp420 = MulSIMD(FDPart3tmp3, FDPart3tmp411);
        const REAL_SIMD_ARRAY FDPart3tmp422 = MulSIMD(FDPart3tmp307, FDPart3tmp411);
        const REAL_SIMD_ARRAY FDPart3tmp424 = MulSIMD(FDPart3tmp292, FDPart3tmp411);
        const REAL_SIMD_ARRAY FDPart3tmp427 = MulSIMD(FDPart3tmp411, FDPart3tmp426);
        const REAL_SIMD_ARRAY FDPart3tmp448 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp84, aDD12),
            MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp74), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp96, aDD02),
                MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp84),
                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp74, aDD02),
                    MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp411),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp87, aDD01),
                        MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp411),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                        FusedMulAddSIMD(
                            FDPart3tmp74,
                            MulSIMD(MulSIMD(FDPart3tmp34, FDPart3tmp411),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp96, aDD22))),
                            FusedMulAddSIMD(
                                aDD01,
                                MulSIMD(MulSIMD(FDPart3tmp410, FDPart3tmp411),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp96, aDD12),
                                    MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp87),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
                                    FusedMulSubSIMD(FDPart3tmp84,
                                                    MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp411),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp87, aDD11))),
                                                    MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp411),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0), MulSIMD(FDPart3tmp84, aDD00)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp450 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp96, aDD12),
            MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp74), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp96, aDD01),
                MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp84),
                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp74, aDD01),
                    MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp411),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp99, aDD02),
                        MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp411),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
                        FusedMulAddSIMD(
                            FDPart3tmp96,
                            MulSIMD(MulSIMD(FDPart3tmp34, FDPart3tmp411),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp99, aDD22))),
                            FusedMulAddSIMD(
                                aDD02,
                                MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp412),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp99, aDD12),
                                    MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp84),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
                                    FusedMulSubSIMD(FDPart3tmp74,
                                                    MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp411),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp84, aDD11))),
                                                    MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp411),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0), MulSIMD(FDPart3tmp96, aDD00)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp452 = MulSIMD(FDPart3tmp411, FDPart3tmp77);
        const REAL_SIMD_ARRAY FDPart3tmp453 = MulSIMD(FDPart3tmp411, FDPart3tmp89);
        const REAL_SIMD_ARRAY FDPart3tmp454 = MulSIMD(FDPart3tmp293, FDPart3tmp411);
        const REAL_SIMD_ARRAY FDPart3tmp500 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp87, aDD12),
            MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp74), MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp34, FDPart3tmp411), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp421, aDD22)),
                FusedMulAddSIMD(MulSIMD(FDPart3tmp84, aDD02),
                                MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp74),
                                        MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
                                FusedMulAddSIMD(MulSIMD(FDPart3tmp87, aDD01),
                                                MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp84),
                                                        MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                                                FusedMulSubSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp411),
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                        MulSIMD(aDD11, MulSIMD(FDPart3tmp87, FDPart3tmp87))),
                                                                MulSIMD(FDPart3tmp410, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0),
                                                                                               MulSIMD(FDPart3tmp411, aDD00))))))));
        const REAL_SIMD_ARRAY FDPart3tmp502 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp99, aDD02),
            MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp84),
                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp96, aDD01),
                MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp87),
                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp84, aDD01),
                    MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp74),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp96, aDD02),
                        MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp74),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
                        FusedMulAddSIMD(
                            FDPart3tmp74,
                            MulSIMD(MulSIMD(FDPart3tmp34, FDPart3tmp411),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp99, aDD22))),
                            FusedMulAddSIMD(
                                aDD12,
                                MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp421),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp99, aDD12),
                                    MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp87),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
                                    FusedMulSubSIMD(FDPart3tmp74,
                                                    MulSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp411),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp87, aDD11))),
                                                    MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp84),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0), MulSIMD(FDPart3tmp96, aDD00)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp523 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp99, aDD12),
            MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp74), MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f3_of_xx2))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp34, FDPart3tmp411),
                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD22, MulSIMD(FDPart3tmp99, FDPart3tmp99))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp96, aDD01),
                    MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp74),
                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f0_of_xx0__D0))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp99, aDD02),
                        MulSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp96),
                                MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, f3_of_xx2))),
                        FusedMulSubSIMD(MulSIMD(FDPart3tmp15, FDPart3tmp411),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp421, aDD11)),
                                        MulSIMD(FDPart3tmp411, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0), MulSIMD(FDPart3tmp412, aDD00))))))));
        const REAL_SIMD_ARRAY FDPart3tmp536 = MulSIMD(FDPart3_Rational_1_3, FDPart3tmp378);
        const REAL_SIMD_ARRAY FDPart3tmp140 = MulSIMD(FDPart3tmp100, FDPart3tmp139);
        const REAL_SIMD_ARRAY FDPart3tmp155 = MulSIMD(FDPart3tmp104, FDPart3tmp154);
        const REAL_SIMD_ARRAY FDPart3tmp164 = MulSIMD(FDPart3tmp163, FDPart3tmp88);
        const REAL_SIMD_ARRAY FDPart3tmp175 = MulSIMD(FDPart3_Integer_2, FDPart3tmp85);
        const REAL_SIMD_ARRAY FDPart3tmp185 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp67, FDPart3tmp74));
        const REAL_SIMD_ARRAY FDPart3tmp190 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp67, FDPart3tmp96));
        const REAL_SIMD_ARRAY FDPart3tmp297 =
            FusedMulAddSIMD(FDPart3tmp292, FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp89, FDPart3tmp97, MulSIMD(FDPart3tmp100, FDPart3tmp293)));
        const REAL_SIMD_ARRAY FDPart3tmp298 =
            FusedMulAddSIMD(FDPart3tmp292, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp293, FDPart3tmp97, MulSIMD(FDPart3tmp104, FDPart3tmp89)));
        const REAL_SIMD_ARRAY FDPart3tmp308 =
            FusedMulAddSIMD(FDPart3tmp307, FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp77, FDPart3tmp85, MulSIMD(FDPart3tmp293, FDPart3tmp88)));
        const REAL_SIMD_ARRAY FDPart3tmp311 =
            FusedMulAddSIMD(FDPart3tmp293, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp307, FDPart3tmp97, MulSIMD(FDPart3tmp104, FDPart3tmp77)));
        const REAL_SIMD_ARRAY FDPart3tmp334 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp85, alpha));
        const REAL_SIMD_ARRAY FDPart3tmp337 = MulSIMD(FDPart3_Rational_3_2, FDPart3tmp97);
        const REAL_SIMD_ARRAY FDPart3tmp340 = MulSIMD(FDPart3_Rational_3_2, FDPart3tmp75);
        const REAL_SIMD_ARRAY FDPart3tmp386 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp97);
        const REAL_SIMD_ARRAY FDPart3tmp388 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp104);
        const REAL_SIMD_ARRAY FDPart3tmp390 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp85);
        const REAL_SIMD_ARRAY FDPart3tmp430 = MulSIMD(FDPart3tmp411, MulSIMD(FDPart3tmp428, FDPart3tmp84));
        const REAL_SIMD_ARRAY FDPart3tmp433 = MulSIMD(FDPart3tmp411, MulSIMD(FDPart3tmp431, FDPart3tmp74));
        const REAL_SIMD_ARRAY FDPart3tmp455 = FusedMulAddSIMD(
            FDPart3tmp453, MulSIMD(FDPart3tmp84, FDPart3tmp96),
            FusedMulAddSIMD(
                FDPart3tmp103, MulSIMD(FDPart3tmp453, FDPart3tmp74),
                FusedMulAddSIMD(
                    FDPart3tmp422, MulSIMD(FDPart3tmp96, FDPart3tmp99),
                    FusedMulAddSIMD(FDPart3tmp103, MulSIMD(FDPart3tmp420, FDPart3tmp96),
                                    FusedMulAddSIMD(FDPart3tmp103, MulSIMD(FDPart3tmp452, FDPart3tmp99),
                                                    FusedMulAddSIMD(FDPart3tmp454, MulSIMD(FDPart3tmp74, FDPart3tmp96),
                                                                    FusedMulAddSIMD(FDPart3tmp454, MulSIMD(FDPart3tmp84, FDPart3tmp99),
                                                                                    FusedMulAddSIMD(FDPart3tmp412, FDPart3tmp452,
                                                                                                    MulSIMD(FDPart3tmp424, FDPart3tmp425)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp468 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp97);
        const REAL_SIMD_ARRAY FDPart3tmp469 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp75);
        const REAL_SIMD_ARRAY FDPart3tmp470 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp85);
        const REAL_SIMD_ARRAY FDPart3tmp471 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp88);
        const REAL_SIMD_ARRAY FDPart3tmp472 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp100);
        const REAL_SIMD_ARRAY FDPart3tmp473 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp104);
        const REAL_SIMD_ARRAY FDPart3tmp481 = MulSIMD(
            FDPart3tmp88, FusedMulAddSIMD(FDPart3tmp5, vetU_dDD111,
                                          FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp5, vetU_dD01), MulSIMD(FDPart3tmp19, FDPart3tmp346))));
        const REAL_SIMD_ARRAY FDPart3tmp484 = MulSIMD(
            FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp382, vetU0, FusedMulAddSIMD(FDPart3tmp13, MulSIMD(FDPart3tmp5, f0_of_xx0__DD00), FDPart3tmp483)));
        const REAL_SIMD_ARRAY FDPart3tmp486 = MulSIMD(FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp485, FDPart3tmp483));
        const REAL_SIMD_ARRAY FDPart3tmp491 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp75);
        const REAL_SIMD_ARRAY FDPart3tmp493 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp88);
        const REAL_SIMD_ARRAY FDPart3tmp509 =
            MulSIMD(FDPart3tmp88, FusedMulAddSIMD(FDPart3tmp8, vetU_dDD211, MulSIMD(FDPart3tmp346, MulSIMD(FDPart3tmp8, vetU_dD20))));
        const REAL_SIMD_ARRAY FDPart3tmp515 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp100);
        const REAL_SIMD_ARRAY FDPart3tmp552 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp192, FDPart3tmp97));
        const REAL_SIMD_ARRAY FDPart3tmp554 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp192, FDPart3tmp75));
        const REAL_SIMD_ARRAY FDPart3tmp556 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp192, FDPart3tmp85));
        const REAL_SIMD_ARRAY FDPart3tmp176 = MulSIMD(FDPart3tmp174, FDPart3tmp175);
        const REAL_SIMD_ARRAY FDPart3tmp186 = MulSIMD(FDPart3tmp183, FDPart3tmp185);
        const REAL_SIMD_ARRAY FDPart3tmp191 = MulSIMD(FDPart3tmp189, FDPart3tmp190);
        const REAL_SIMD_ARRAY FDPart3tmp214 = SubSIMD(
            FusedMulAddSIMD(
                FDPart3tmp200, MulSIMD(cf_dD0, cf_dD2),
                SubSIMD(
                    SubSIMD(AddSIMD(FusedMulAddSIMD(RbarDD02, alpha, FDPart3tmp199),
                                    NegFusedMulAddSIMD(
                                        FDPart3_Integer_2,
                                        MulSIMD(alpha,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp193,
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp209, cf_dD1)),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp193,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp211, cf_dD2)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp213, FusedMulSubSIMD(FDPart3tmp197, cf_dD2, cf_dDD02),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp193), MulSIMD(FDPart3tmp207, cf_dD0)))))),
                                        SubSIMD(FDPart3tmp196, alpha_dDD02))),
                            MulSIMD(FDPart3tmp209, alpha_dD1)),
                    MulSIMD(FDPart3tmp207, alpha_dD0))),
            MulSIMD(FDPart3tmp211, alpha_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp228 = SubSIMD(
            FusedMulAddSIMD(
                FDPart3tmp200, MulSIMD(cf_dD1, cf_dD2),
                SubSIMD(
                    SubSIMD(AddSIMD(FusedMulAddSIMD(RbarDD12, alpha, FDPart3tmp219),
                                    NegFusedMulAddSIMD(
                                        FDPart3_Integer_2,
                                        MulSIMD(alpha,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp193,
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp224, cf_dD1)),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp193,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp226, cf_dD2)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp213, FusedMulSubSIMD(FDPart3tmp217, cf_dD2, cf_dDD12),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp193), MulSIMD(FDPart3tmp222, cf_dD0)))))),
                                        SubSIMD(FDPart3tmp216, alpha_dDD12))),
                            MulSIMD(FDPart3tmp224, alpha_dD1)),
                    MulSIMD(FDPart3tmp222, alpha_dD0))),
            MulSIMD(FDPart3tmp226, alpha_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp241 = SubSIMD(
            FusedMulAddSIMD(
                FDPart3tmp200, MulSIMD(cf_dD0, cf_dD1),
                SubSIMD(
                    SubSIMD(AddSIMD(FusedMulAddSIMD(RbarDD01, alpha, FDPart3tmp230),
                                    NegFusedMulAddSIMD(
                                        FDPart3_Integer_2,
                                        MulSIMD(alpha,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp193,
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp237, cf_dD1)),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp193,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp239, cf_dD2)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp213, FusedMulSubSIMD(FDPart3tmp217, cf_dD0, cf_dDD01),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp193), MulSIMD(FDPart3tmp235, cf_dD0)))))),
                                        SubSIMD(FDPart3tmp229, alpha_dDD01))),
                            MulSIMD(FDPart3tmp237, alpha_dD1)),
                    MulSIMD(FDPart3tmp235, alpha_dD0))),
            MulSIMD(FDPart3tmp239, alpha_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp253 = NegFusedMulAddSIMD(
            FDPart3_Integer_2,
            MulSIMD(alpha,
                    FusedMulAddSIMD(
                        FDPart3tmp193, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp249, cf_dD1)),
                        FusedMulAddSIMD(FDPart3tmp193, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp251, cf_dD2)),
                                        FusedMulSubSIMD(FDPart3tmp213, FusedMulSubSIMD(FDPart3tmp193, MulSIMD(cf_dD1, cf_dD1), cf_dDD11),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp193), MulSIMD(FDPart3tmp247, cf_dD0)))))),
            SubSIMD(
                SubSIMD(FusedMulAddSIMD(RbarDD11, alpha,
                                        SubSIMD(FusedMulAddSIMD(FDPart3tmp193,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha_dD1, cf_dD1)),
                                                                FusedMulSubSIMD(FDPart3tmp200, MulSIMD(cf_dD1, cf_dD1), alpha_dDD11)),
                                                MulSIMD(FDPart3tmp247, alpha_dD0))),
                        MulSIMD(FDPart3tmp251, alpha_dD2)),
                MulSIMD(FDPart3tmp249, alpha_dD1)));
        const REAL_SIMD_ARRAY FDPart3tmp268 = NegFusedMulAddSIMD(
            FDPart3_Integer_2,
            MulSIMD(alpha,
                    FusedMulAddSIMD(
                        FDPart3tmp193, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp264, cf_dD1)),
                        FusedMulAddSIMD(FDPart3tmp193, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp266, cf_dD2)),
                                        FusedMulSubSIMD(FDPart3tmp213, FusedMulSubSIMD(FDPart3tmp193, MulSIMD(cf_dD2, cf_dD2), cf_dDD22),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp193), MulSIMD(FDPart3tmp262, cf_dD0)))))),
            SubSIMD(
                SubSIMD(FusedMulAddSIMD(RbarDD22, alpha,
                                        SubSIMD(FusedMulAddSIMD(FDPart3tmp193,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha_dD2, cf_dD2)),
                                                                FusedMulSubSIMD(FDPart3tmp200, MulSIMD(cf_dD2, cf_dD2), alpha_dDD22)),
                                                MulSIMD(FDPart3tmp262, alpha_dD0))),
                        MulSIMD(FDPart3tmp266, alpha_dD2)),
                MulSIMD(FDPart3tmp264, alpha_dD1)));
        const REAL_SIMD_ARRAY FDPart3tmp283 = NegFusedMulAddSIMD(
            FDPart3_Integer_2,
            MulSIMD(alpha,
                    FusedMulAddSIMD(
                        FDPart3tmp193, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp279, cf_dD1)),
                        FusedMulAddSIMD(FDPart3tmp193, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp281, cf_dD2)),
                                        FusedMulSubSIMD(FDPart3tmp213, FusedMulSubSIMD(FDPart3tmp193, MulSIMD(cf_dD0, cf_dD0), cf_dDD00),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp193), MulSIMD(FDPart3tmp277, cf_dD0)))))),
            SubSIMD(
                SubSIMD(FusedMulAddSIMD(RbarDD00, alpha,
                                        SubSIMD(FusedMulAddSIMD(FDPart3tmp193,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha_dD0, cf_dD0)),
                                                                FusedMulSubSIMD(FDPart3tmp200, MulSIMD(cf_dD0, cf_dD0), alpha_dDD00)),
                                                MulSIMD(FDPart3tmp277, alpha_dD0))),
                        MulSIMD(FDPart3tmp281, alpha_dD2)),
                MulSIMD(FDPart3tmp279, alpha_dD1)));
        const REAL_SIMD_ARRAY FDPart3tmp296 =
            FusedMulAddSIMD(FDPart3tmp293, FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp89, MulSIMD(FDPart3tmp292, FDPart3tmp88)));
        const REAL_SIMD_ARRAY FDPart3tmp310 =
            FusedMulAddSIMD(FDPart3tmp293, FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp77, FDPart3tmp97, MulSIMD(FDPart3tmp100, FDPart3tmp307)));
        const REAL_SIMD_ARRAY FDPart3tmp403 =
            FusedMulAddSIMD(FDPart3tmp246, FDPart3tmp388, FusedMulAddSIMD(FDPart3tmp390, FDPart3tmp402, MulSIMD(FDPart3tmp245, FDPart3tmp386)));
        const REAL_SIMD_ARRAY FDPart3tmp406 =
            FusedMulAddSIMD(FDPart3tmp259, FDPart3tmp390, FusedMulAddSIMD(FDPart3tmp261, FDPart3tmp388, MulSIMD(FDPart3tmp257, FDPart3tmp386)));
        const REAL_SIMD_ARRAY FDPart3tmp408 =
            FusedMulAddSIMD(FDPart3tmp274, FDPart3tmp390, FusedMulAddSIMD(FDPart3tmp276, FDPart3tmp388, MulSIMD(FDPart3tmp270, FDPart3tmp386)));
        const REAL_SIMD_ARRAY FDPart3tmp434 = FusedMulAddSIMD(
            FDPart3tmp424, MulSIMD(FDPart3tmp87, FDPart3tmp87),
            FusedMulAddSIMD(FDPart3tmp425, FDPart3tmp427,
                            FusedMulAddSIMD(FDPart3tmp430, FDPart3tmp87,
                                            FusedMulAddSIMD(FDPart3tmp433, FDPart3tmp87,
                                                            FusedMulAddSIMD(FDPart3tmp410, FDPart3tmp420, MulSIMD(FDPart3tmp421, FDPart3tmp422))))));
        const REAL_SIMD_ARRAY FDPart3tmp440 = FusedMulAddSIMD(
            FDPart3tmp422, MulSIMD(FDPart3tmp99, FDPart3tmp99),
            FusedMulAddSIMD(FDPart3tmp433, FDPart3tmp99,
                            FusedMulAddSIMD(FDPart3tmp427, MulSIMD(FDPart3tmp96, FDPart3tmp99),
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp428), MulSIMD(FDPart3tmp74, FDPart3tmp96),
                                                            FusedMulAddSIMD(FDPart3tmp412, FDPart3tmp420, MulSIMD(FDPart3tmp421, FDPart3tmp424))))));
        const REAL_SIMD_ARRAY FDPart3tmp446 = FusedMulAddSIMD(
            FDPart3tmp412, FDPart3tmp422,
            FusedMulAddSIMD(FDPart3tmp420, MulSIMD(FDPart3tmp103, FDPart3tmp103),
                            FusedMulAddSIMD(FDPart3tmp103, MulSIMD(FDPart3tmp427, FDPart3tmp96),
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp411, FDPart3tmp431), MulSIMD(FDPart3tmp84, FDPart3tmp96),
                                                            FusedMulAddSIMD(FDPart3tmp103, FDPart3tmp430, MulSIMD(FDPart3tmp410, FDPart3tmp424))))));
        const REAL_SIMD_ARRAY FDPart3tmp456 = MulSIMD(FDPart3tmp416, FDPart3tmp455);
        const REAL_SIMD_ARRAY FDPart3tmp464 = FusedMulAddSIMD(
            FDPart3tmp452, MulSIMD(FDPart3tmp84, FDPart3tmp99),
            FusedMulAddSIMD(
                FDPart3tmp424, MulSIMD(FDPart3tmp74, FDPart3tmp87),
                FusedMulAddSIMD(
                    FDPart3tmp452, MulSIMD(FDPart3tmp74, FDPart3tmp96),
                    FusedMulAddSIMD(FDPart3tmp420, MulSIMD(FDPart3tmp84, FDPart3tmp96),
                                    FusedMulAddSIMD(FDPart3tmp422, MulSIMD(FDPart3tmp74, FDPart3tmp99),
                                                    FusedMulAddSIMD(FDPart3tmp453, MulSIMD(FDPart3tmp87, FDPart3tmp96),
                                                                    FusedMulAddSIMD(FDPart3tmp454, MulSIMD(FDPart3tmp87, FDPart3tmp99),
                                                                                    FusedMulAddSIMD(FDPart3tmp421, FDPart3tmp454,
                                                                                                    MulSIMD(FDPart3tmp425, FDPart3tmp453)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp466 = FusedMulAddSIMD(
            FDPart3tmp424, MulSIMD(FDPart3tmp84, FDPart3tmp87),
            FusedMulAddSIMD(
                FDPart3tmp103, MulSIMD(FDPart3tmp453, FDPart3tmp87),
                FusedMulAddSIMD(
                    FDPart3tmp422, MulSIMD(FDPart3tmp74, FDPart3tmp96),
                    FusedMulAddSIMD(FDPart3tmp103, MulSIMD(FDPart3tmp420, FDPart3tmp84),
                                    FusedMulAddSIMD(FDPart3tmp103, MulSIMD(FDPart3tmp452, FDPart3tmp74),
                                                    FusedMulAddSIMD(FDPart3tmp452, MulSIMD(FDPart3tmp84, FDPart3tmp96),
                                                                    FusedMulAddSIMD(FDPart3tmp454, MulSIMD(FDPart3tmp87, FDPart3tmp96),
                                                                                    FusedMulAddSIMD(FDPart3tmp410, FDPart3tmp453,
                                                                                                    MulSIMD(FDPart3tmp425, FDPart3tmp454)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp497 =
            FusedMulAddSIMD(FDPart3tmp246, FDPart3tmp390, FusedMulAddSIMD(FDPart3tmp402, FDPart3tmp493, MulSIMD(FDPart3tmp245, FDPart3tmp491)));
        const REAL_SIMD_ARRAY FDPart3tmp498 =
            FusedMulAddSIMD(FDPart3tmp259, FDPart3tmp493, FusedMulAddSIMD(FDPart3tmp261, FDPart3tmp390, MulSIMD(FDPart3tmp257, FDPart3tmp491)));
        const REAL_SIMD_ARRAY FDPart3tmp499 =
            FusedMulAddSIMD(FDPart3tmp274, FDPart3tmp493, FusedMulAddSIMD(FDPart3tmp276, FDPart3tmp390, MulSIMD(FDPart3tmp270, FDPart3tmp491)));
        const REAL_SIMD_ARRAY FDPart3tmp520 =
            FusedMulAddSIMD(FDPart3tmp246, FDPart3tmp386, FusedMulAddSIMD(FDPart3tmp402, FDPart3tmp491, MulSIMD(FDPart3tmp245, FDPart3tmp515)));
        const REAL_SIMD_ARRAY FDPart3tmp521 =
            FusedMulAddSIMD(FDPart3tmp259, FDPart3tmp491, FusedMulAddSIMD(FDPart3tmp261, FDPart3tmp386, MulSIMD(FDPart3tmp257, FDPart3tmp515)));
        const REAL_SIMD_ARRAY FDPart3tmp522 =
            FusedMulAddSIMD(FDPart3tmp274, FDPart3tmp491, FusedMulAddSIMD(FDPart3tmp276, FDPart3tmp386, MulSIMD(FDPart3tmp270, FDPart3tmp515)));
        const REAL_SIMD_ARRAY FDPart3tmp527 = FusedMulAddSIMD(
            FDPart3tmp1, MulSIMD(f0_of_xx0__DD00, vetU0),
            SubSIMD(FusedMulAddSIMD(
                        alpha,
                        FusedMulAddSIMD(FDPart3tmp292, FDPart3tmp88,
                                        FusedMulAddSIMD(FDPart3tmp426, FDPart3tmp97,
                                                        FusedMulAddSIMD(FDPart3tmp428, FDPart3tmp85,
                                                                        FusedMulAddSIMD(FDPart3tmp431, FDPart3tmp75,
                                                                                        FusedMulAddSIMD(FDPart3tmp100, FDPart3tmp307,
                                                                                                        MulSIMD(FDPart3tmp104, FDPart3tmp3)))))),
                        FusedMulAddSIMD(FDPart3tmp37,
                                        MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp32),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp38, vetU0))),
                                        MulSIMD(FDPart3_NegativeOne_, AddSIMD(FDPart3tmp29, FDPart3tmp22)))),
                    FDPart3tmp30));
        const REAL_SIMD_ARRAY FDPart3tmp546 = MulSIMD(FDPart3tmp455, FDPart3tmp545);
        const REAL_SIMD_ARRAY FDPart3tmp550 = FusedMulAddSIMD(
            FDPart3tmp209, alpha_dD1, FusedMulAddSIMD(FDPart3tmp211, alpha_dD2, FusedMulAddSIMD(FDPart3tmp207, alpha_dD0, alpha_dDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp553 = FusedMulAddSIMD(
            FDPart3tmp224, alpha_dD1, FusedMulAddSIMD(FDPart3tmp226, alpha_dD2, FusedMulAddSIMD(FDPart3tmp222, alpha_dD0, alpha_dDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp555 = FusedMulAddSIMD(
            FDPart3tmp237, alpha_dD1, FusedMulAddSIMD(FDPart3tmp239, alpha_dD2, FusedMulAddSIMD(FDPart3tmp235, alpha_dD0, alpha_dDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp284 = FusedMulAddSIMD(
            FDPart3tmp175, FDPart3tmp241,
            FusedMulAddSIMD(FDPart3tmp185, FDPart3tmp228,
                            FusedMulAddSIMD(FDPart3tmp190, FDPart3tmp214,
                                            FusedMulAddSIMD(FDPart3tmp253, FDPart3tmp88,
                                                            FusedMulAddSIMD(FDPart3tmp100, FDPart3tmp268, MulSIMD(FDPart3tmp104, FDPart3tmp283))))));
        const REAL_SIMD_ARRAY FDPart3tmp391 =
            FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp388, hDD_dD002),
                            FusedMulAddSIMD(FDPart3tmp34, MulSIMD(FDPart3tmp386, hDD_dD220), MulSIMD(FDPart3tmp206, FDPart3tmp390)));
        const REAL_SIMD_ARRAY FDPart3tmp395 = FusedMulAddSIMD(
            FDPart3tmp244, FDPart3tmp390, FusedMulAddSIMD(FDPart3tmp34, MulSIMD(FDPart3tmp386, hDD_dD221), MulSIMD(FDPart3tmp221, FDPart3tmp388)));
        const REAL_SIMD_ARRAY FDPart3tmp399 = FusedMulAddSIMD(
            FDPart3tmp234, FDPart3tmp386, FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp388, hDD_dD001), MulSIMD(FDPart3tmp233, FDPart3tmp390)));
        const REAL_SIMD_ARRAY FDPart3tmp419 = AddSIMD(FDPart3tmp346, FDPart3tmp403);
        const REAL_SIMD_ARRAY FDPart3tmp436 = MulSIMD(FDPart3tmp434, FDPart3tmp435);
        const REAL_SIMD_ARRAY FDPart3tmp443 = AddSIMD(FDPart3tmp408, FDPart3tmp442);
        const REAL_SIMD_ARRAY FDPart3tmp447 = MulSIMD(FDPart3tmp435, FDPart3tmp446);
        const REAL_SIMD_ARRAY FDPart3tmp465 = MulSIMD(FDPart3tmp416, FDPart3tmp464);
        const REAL_SIMD_ARRAY FDPart3tmp467 = MulSIMD(FDPart3tmp416, FDPart3tmp466);
        const REAL_SIMD_ARRAY FDPart3tmp494 =
            FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp390, hDD_dD002),
                            FusedMulAddSIMD(FDPart3tmp34, MulSIMD(FDPart3tmp491, hDD_dD220), MulSIMD(FDPart3tmp206, FDPart3tmp493)));
        const REAL_SIMD_ARRAY FDPart3tmp495 = FusedMulAddSIMD(
            FDPart3tmp244, FDPart3tmp493, FusedMulAddSIMD(FDPart3tmp34, MulSIMD(FDPart3tmp491, hDD_dD221), MulSIMD(FDPart3tmp221, FDPart3tmp390)));
        const REAL_SIMD_ARRAY FDPart3tmp496 = FusedMulAddSIMD(
            FDPart3tmp234, FDPart3tmp491, FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp390, hDD_dD001), MulSIMD(FDPart3tmp233, FDPart3tmp493)));
        const REAL_SIMD_ARRAY FDPart3tmp517 =
            FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp386, hDD_dD002),
                            FusedMulAddSIMD(FDPart3tmp34, MulSIMD(FDPart3tmp515, hDD_dD220), MulSIMD(FDPart3tmp206, FDPart3tmp491)));
        const REAL_SIMD_ARRAY FDPart3tmp518 = FusedMulAddSIMD(
            FDPart3tmp244, FDPart3tmp491, FusedMulAddSIMD(FDPart3tmp34, MulSIMD(FDPart3tmp515, hDD_dD221), MulSIMD(FDPart3tmp221, FDPart3tmp386)));
        const REAL_SIMD_ARRAY FDPart3tmp519 = FusedMulAddSIMD(
            FDPart3tmp234, FDPart3tmp515, FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp386, hDD_dD001), MulSIMD(FDPart3tmp233, FDPart3tmp491)));
        const REAL_SIMD_ARRAY FDPart3tmp524 = NegFusedMulAddSIMD(FDPart3tmp8, f3_of_xx2__D2, FDPart3tmp521);
        const REAL_SIMD_ARRAY FDPart3tmp529 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp527);
        const REAL_SIMD_ARRAY FDPart3tmp541 = MulSIMD(FDPart3tmp434, FDPart3tmp540);
        const REAL_SIMD_ARRAY FDPart3tmp542 = MulSIMD(FDPart3tmp446, FDPart3tmp540);
        const REAL_SIMD_ARRAY FDPart3tmp547 = MulSIMD(FDPart3tmp464, FDPart3tmp545);
        const REAL_SIMD_ARRAY FDPart3tmp548 = MulSIMD(FDPart3tmp466, FDPart3tmp545);
        const REAL_SIMD_ARRAY FDPart3tmp474 = MulSIMD(
            FDPart3tmp42, FusedMulAddSIMD(FDPart3tmp399, FDPart3tmp470,
                                          FusedMulAddSIMD(FDPart3tmp406, FDPart3tmp472,
                                                          FusedMulAddSIMD(FDPart3tmp419, FDPart3tmp471,
                                                                          FusedMulAddSIMD(FDPart3tmp443, FDPart3tmp473,
                                                                                          FusedMulAddSIMD(FDPart3tmp391, FDPart3tmp468,
                                                                                                          MulSIMD(FDPart3tmp395, FDPart3tmp469)))))));
        const REAL_SIMD_ARRAY FDPart3tmp503 = NegFusedMulAddSIMD(FDPart3tmp5, f0_of_xx0__D0, FDPart3tmp496);
        const REAL_SIMD_ARRAY FDPart3tmp525 = MulSIMD(
            FDPart3tmp42, FusedMulAddSIMD(FDPart3tmp470, FDPart3tmp519,
                                          FusedMulAddSIMD(FDPart3tmp471, FDPart3tmp520,
                                                          FusedMulAddSIMD(FDPart3tmp472, FDPart3tmp524,
                                                                          FusedMulAddSIMD(FDPart3tmp473, FDPart3tmp522,
                                                                                          FusedMulAddSIMD(FDPart3tmp468, FDPart3tmp517,
                                                                                                          MulSIMD(FDPart3tmp469, FDPart3tmp518)))))));
        const REAL_SIMD_ARRAY FDPart3tmp504 = MulSIMD(
            FDPart3tmp42, FusedMulAddSIMD(FDPart3tmp470, FDPart3tmp503,
                                          FusedMulAddSIMD(FDPart3tmp471, FDPart3tmp497,
                                                          FusedMulAddSIMD(FDPart3tmp472, FDPart3tmp498,
                                                                          FusedMulAddSIMD(FDPart3tmp473, FDPart3tmp499,
                                                                                          FusedMulAddSIMD(FDPart3tmp468, FDPart3tmp494,
                                                                                                          MulSIMD(FDPart3tmp469, FDPart3tmp495)))))));
        const REAL_SIMD_ARRAY __RHS_exp_0 = MulSIMD(
            FDPart3tmp1,
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp106, PI),
                MulSIMD(
                    MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                    MulSIMD(alpha, FusedMulAddSIMD(
                                       FDPart3tmp109, FDPart3tmp176,
                                       FusedMulAddSIMD(FDPart3tmp109, FDPart3tmp155,
                                                       FusedMulAddSIMD(FDPart3tmp109, FDPart3tmp164,
                                                                       FusedMulAddSIMD(FDPart3tmp109, FDPart3tmp186,
                                                                                       FusedMulAddSIMD(FDPart3tmp109, FDPart3tmp191,
                                                                                                       FusedMulAddSIMD(FDPart3tmp109, FDPart3tmp140,
                                                                                                                       FDPart3tmp154)))))))),
                FusedMulAddSIMD(
                    f0_of_xx0__D0,
                    MulSIMD(MulSIMD(aDD01, f0_of_xx0), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                               MulSIMD(alpha, FusedMulAddSIMD(FDPart3tmp75, FDPart3tmp77,
                                                                                              FusedMulAddSIMD(FDPart3tmp88, FDPart3tmp89,
                                                                                                              MulSIMD(FDPart3tmp3, FDPart3tmp85)))))),
                    FusedMulAddSIMD(
                        FDPart3tmp2, MulSIMD(aDD02, vetU_dD20),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp0, aDD00),
                            MulSIMD(
                                MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                MulSIMD(alpha, FusedMulAddSIMD(FDPart3tmp77, FDPart3tmp97,
                                                               FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp89, MulSIMD(FDPart3tmp104, FDPart3tmp3))))),
                            FusedMulAddSIMD(
                                FDPart3tmp0, MulSIMD(FDPart3tmp9, aDD_dupD002),
                                FusedMulAddSIMD(
                                    FDPart3tmp19, MulSIMD(FDPart3tmp20, aDD01),
                                    FusedMulAddSIMD(
                                        FDPart3_Integer_2, MulSIMD(FDPart3tmp26, FDPart3tmp3),
                                        FusedMulAddSIMD(
                                            FDPart3tmp0, MulSIMD(FDPart3tmp6, aDD_dupD001),
                                            FusedMulAddSIMD(
                                                FDPart3tmp192,
                                                NegFusedMulAddSIMD(
                                                    FDPart3tmp284,
                                                    FusedMulAddSIMD(FDPart3_Rational_1_3, FDPart3tmp0, MulSIMD(FDPart3_Rational_1_3, FDPart3tmp53)),
                                                    FDPart3tmp283),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp3, FDPart3tmp4,
                                                    FusedMulAddSIMD(
                                                        f3_of_xx2,
                                                        MulSIMD(MulSIMD(aDD02, f0_of_xx0__D0),
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                        MulSIMD(alpha, FusedMulAddSIMD(
                                                                                           FDPart3tmp3, FDPart3tmp97,
                                                                                           FusedMulAddSIMD(FDPart3tmp75, FDPart3tmp89,
                                                                                                           MulSIMD(FDPart3tmp100, FDPart3tmp77)))))),
                                                        FusedMulSubSIMD(FDPart3tmp13,
                                                                        FusedMulAddSIMD(FDPart3tmp0, aDD_dupD000, MulSIMD(FDPart3tmp11, aDD00)),
                                                                        MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp0),
                                                                                MulSIMD(FDPart3tmp42, aDD00)))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_1 = MulSIMD(
            FDPart3tmp299,
            FusedMulAddSIMD(
                f0_of_xx0,
                MulSIMD(MulSIMD(FDPart3tmp296, aDD01), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, alpha))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp106, PI),
                    MulSIMD(
                        MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                        MulSIMD(alpha,
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp84, f0_of_xx0),
                                    MulSIMD(MulSIMD(FDPart3tmp174, FDPart3tmp67),
                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0__D0, hDD01))),
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3tmp74, f0_of_xx0),
                                        MulSIMD(MulSIMD(FDPart3tmp183, FDPart3tmp67),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0__D0, hDD01))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3tmp99, f0_of_xx0),
                                            MulSIMD(MulSIMD(FDPart3tmp139, FDPart3tmp67),
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f0_of_xx0__D0, hDD01))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp87, f0_of_xx0),
                                                MulSIMD(MulSIMD(FDPart3tmp163, FDPart3tmp67),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f0_of_xx0__D0, hDD01))),
                                                FusedMulAddSIMD(MulSIMD(FDPart3tmp96, f0_of_xx0),
                                                                MulSIMD(MulSIMD(FDPart3tmp189, FDPart3tmp67),
                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3),
                                                                                MulSIMD(f0_of_xx0__D0, hDD01))),
                                                                NegFusedMulAddSIMD(f0_of_xx0,
                                                                                   MulSIMD(MulSIMD(FDPart3tmp154, FDPart3tmp67),
                                                                                           MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp103),
                                                                                                   MulSIMD(f0_of_xx0__D0, hDD01))),
                                                                                   FDPart3tmp174)))))))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp42, aDD01), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0, f0_of_xx0__D0)),
                        FusedMulAddSIMD(
                            aDD01, MulSIMD(f0_of_xx0__D0, vetU_dD11),
                            FusedMulAddSIMD(
                                aDD02, MulSIMD(f0_of_xx0__D0, vetU_dD21),
                                FusedMulAddSIMD(
                                    FDPart3tmp78, MulSIMD(FDPart3tmp9, aDD_dupD012),
                                    FusedMulAddSIMD(
                                        aDD00, MulSIMD(f0_of_xx0__D0, vetU_dD01),
                                        FusedMulAddSIMD(
                                            FDPart3tmp4, FDPart3tmp89,
                                            FusedMulAddSIMD(
                                                FDPart3tmp26, MulSIMD(FDPart3tmp78, aDD01),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp192,
                                                    NegFusedMulAddSIMD(
                                                        f0_of_xx0,
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp284), MulSIMD(f0_of_xx0__D0, hDD01)),
                                                        FDPart3tmp241),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp289, vetU_dD20,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp14, aDD_dupD011,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp19, FDPart3tmp292,
                                                                FusedMulAddSIMD(
                                                                    f0_of_xx0__D0,
                                                                    MulSIMD(MulSIMD(FDPart3tmp297, aDD02),
                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                    MulSIMD(f3_of_xx2, alpha))),
                                                                    FusedMulSubSIMD(
                                                                        FDPart3tmp13,
                                                                        FusedMulAddSIMD(FDPart3tmp78, aDD_dupD010,
                                                                                        MulSIMD(aDD01, AddSIMD(FDPart3tmp0, FDPart3tmp271))),
                                                                        MulSIMD(FDPart3tmp298, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0),
                                                                                                       MulSIMD(aDD00, alpha)))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_2 = MulSIMD(
            FDPart3tmp312,
            FusedMulAddSIMD(
                f0_of_xx0,
                MulSIMD(MulSIMD(FDPart3tmp308, aDD01), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, alpha))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp106, PI),
                    MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                            MulSIMD(alpha,
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3tmp84, f0_of_xx0__D0),
                                        MulSIMD(MulSIMD(FDPart3tmp174, FDPart3tmp67),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f3_of_xx2, hDD02))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3tmp74, f0_of_xx0__D0),
                                            MulSIMD(MulSIMD(FDPart3tmp183, FDPart3tmp67),
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f3_of_xx2, hDD02))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp99, f0_of_xx0__D0),
                                                MulSIMD(MulSIMD(FDPart3tmp139, FDPart3tmp67),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f3_of_xx2, hDD02))),
                                                FusedMulAddSIMD(
                                                    MulSIMD(FDPart3tmp87, f0_of_xx0__D0),
                                                    MulSIMD(MulSIMD(FDPart3tmp163, FDPart3tmp67),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f3_of_xx2, hDD02))),
                                                    FusedMulAddSIMD(MulSIMD(FDPart3tmp96, f0_of_xx0__D0),
                                                                    MulSIMD(MulSIMD(FDPart3tmp189, FDPart3tmp67),
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3),
                                                                                    MulSIMD(f3_of_xx2, hDD02))),
                                                                    NegFusedMulAddSIMD(f0_of_xx0__D0,
                                                                                       MulSIMD(MulSIMD(FDPart3tmp154, FDPart3tmp67),
                                                                                               MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp103),
                                                                                                       MulSIMD(f3_of_xx2, hDD02))),
                                                                                       FDPart3tmp189)))))))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp42, aDD02), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0__D0, f3_of_xx2)),
                        FusedMulAddSIMD(
                            aDD01, MulSIMD(f0_of_xx0__D0, vetU_dD12),
                            FusedMulAddSIMD(
                                aDD22, MulSIMD(f3_of_xx2, vetU_dD20),
                                FusedMulAddSIMD(
                                    FDPart3tmp6, MulSIMD(FDPart3tmp76, aDD_dupD021),
                                    FusedMulAddSIMD(
                                        aDD00, MulSIMD(f0_of_xx0__D0, vetU_dD02),
                                        FusedMulAddSIMD(
                                            FDPart3tmp4, FDPart3tmp77,
                                            FusedMulAddSIMD(
                                                FDPart3tmp9,
                                                FusedMulAddSIMD(FDPart3tmp76, aDD_dupD022, MulSIMD(aDD02, MulSIMD(f0_of_xx0__D0, f3_of_xx2__D2))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp26, FDPart3tmp77,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp306, FDPart3tmp77,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp19, FDPart3tmp293,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp192,
                                                                NegFusedMulAddSIMD(
                                                                    f0_of_xx0__D0,
                                                                    MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp284), MulSIMD(f3_of_xx2, hDD02)),
                                                                    FDPart3tmp214),
                                                                FusedMulAddSIMD(
                                                                    f0_of_xx0__D0,
                                                                    MulSIMD(MulSIMD(FDPart3tmp310, aDD02),
                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                    MulSIMD(f3_of_xx2, alpha))),
                                                                    FusedMulSubSIMD(
                                                                        FDPart3tmp13,
                                                                        FusedMulAddSIMD(FDPart3tmp302, f0_of_xx0__DD00,
                                                                                        MulSIMD(FDPart3tmp76, aDD_dupD020)),
                                                                        MulSIMD(FDPart3tmp311, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0),
                                                                                                       MulSIMD(aDD00, alpha)))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_3 = MulSIMD(
            FDPart3tmp16,
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp106, PI),
                MulSIMD(
                    MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                    MulSIMD(alpha, FusedMulAddSIMD(
                                       FDPart3tmp176, FDPart3tmp319,
                                       FusedMulAddSIMD(FDPart3tmp155, FDPart3tmp319,
                                                       FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp319,
                                                                       FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp319,
                                                                                       FusedMulAddSIMD(FDPart3tmp191, FDPart3tmp319,
                                                                                                       FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp319,
                                                                                                                       FDPart3tmp163)))))))),
                FusedMulAddSIMD(
                    f0_of_xx0,
                    MulSIMD(MulSIMD(FDPart3tmp297, aDD12), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f3_of_xx2, alpha))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3_Integer_2, aDD11), MulSIMD(f0_of_xx0, vetU_dD11),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp15, FDPart3tmp296), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD11, alpha)),
                            FusedMulAddSIMD(
                                FDPart3tmp15, MulSIMD(FDPart3tmp9, aDD_dupD112),
                                FusedMulAddSIMD(
                                    aDD_dupD111, MulSIMD(f0_of_xx0, vetU1),
                                    FusedMulAddSIMD(
                                        FDPart3tmp315, aDD01,
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp289, vetU_dD21),
                                            FusedMulAddSIMD(
                                                FDPart3tmp192,
                                                NegFusedMulAddSIMD(
                                                    FDPart3tmp284,
                                                    FusedMulAddSIMD(FDPart3_Rational_1_3, FDPart3tmp15, MulSIMD(FDPart3_Rational_1_3, FDPart3tmp59)),
                                                    FDPart3tmp253),
                                                FusedMulAddSIMD(FDPart3tmp292, FDPart3tmp4,
                                                                FusedMulAddSIMD(f0_of_xx0,
                                                                                MulSIMD(MulSIMD(FDPart3tmp298, aDD01),
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(f0_of_xx0__D0, alpha))),
                                                                                FusedMulSubSIMD(FDPart3tmp13,
                                                                                                FusedMulAddSIMD(FDPart3tmp15, aDD_dupD110,
                                                                                                                MulSIMD(FDPart3tmp20, aDD11)),
                                                                                                MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp15),
                                                                                                        MulSIMD(FDPart3tmp42, aDD11)))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_4 = MulSIMD(
            FDPart3tmp322,
            FusedMulAddSIMD(
                f0_of_xx0,
                MulSIMD(MulSIMD(FDPart3tmp310, aDD12), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f3_of_xx2, alpha))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp106, PI),
                    MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                            MulSIMD(alpha,
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3tmp84, f0_of_xx0),
                                        MulSIMD(MulSIMD(FDPart3tmp174, FDPart3tmp67),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f3_of_xx2, hDD12))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3tmp74, f0_of_xx0),
                                            MulSIMD(MulSIMD(FDPart3tmp183, FDPart3tmp67),
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f3_of_xx2, hDD12))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp99, f0_of_xx0),
                                                MulSIMD(MulSIMD(FDPart3tmp139, FDPart3tmp67),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f3_of_xx2, hDD12))),
                                                FusedMulAddSIMD(
                                                    MulSIMD(FDPart3tmp87, f0_of_xx0),
                                                    MulSIMD(MulSIMD(FDPart3tmp163, FDPart3tmp67),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f3_of_xx2, hDD12))),
                                                    FusedMulAddSIMD(MulSIMD(FDPart3tmp96, f0_of_xx0),
                                                                    MulSIMD(MulSIMD(FDPart3tmp189, FDPart3tmp67),
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3),
                                                                                    MulSIMD(f3_of_xx2, hDD12))),
                                                                    NegFusedMulAddSIMD(f0_of_xx0,
                                                                                       MulSIMD(MulSIMD(FDPart3tmp154, FDPart3tmp67),
                                                                                               MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp103),
                                                                                                       MulSIMD(f3_of_xx2, hDD12))),
                                                                                       FDPart3tmp183)))))))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp42, aDD12), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0, f3_of_xx2)),
                        FusedMulAddSIMD(
                            aDD12, MulSIMD(f3_of_xx2, vetU_dD11),
                            FusedMulAddSIMD(
                                aDD22, MulSIMD(f3_of_xx2, vetU_dD21),
                                FusedMulAddSIMD(
                                    aDD01, MulSIMD(f0_of_xx0, vetU_dD02),
                                    FusedMulAddSIMD(
                                        aDD11, MulSIMD(f0_of_xx0, vetU_dD12),
                                        FusedMulAddSIMD(
                                            FDPart3tmp302, vetU_dD01,
                                            FusedMulAddSIMD(
                                                FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp289, f3_of_xx2__D2, MulSIMD(FDPart3tmp68, aDD_dupD122)),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp293, FDPart3tmp306,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp293, FDPart3tmp4,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp68, aDD_dupD120, MulSIMD(FDPart3tmp76, aDD12)),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp192,
                                                                NegFusedMulAddSIMD(
                                                                    f0_of_xx0,
                                                                    MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp284), MulSIMD(f3_of_xx2, hDD12)),
                                                                    FDPart3tmp228),
                                                                FusedMulAddSIMD(
                                                                    f0_of_xx0,
                                                                    MulSIMD(MulSIMD(FDPart3tmp311, aDD01),
                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                    MulSIMD(f0_of_xx0__D0, alpha))),
                                                                    FusedMulSubSIMD(
                                                                        FDPart3tmp129, aDD_dupD121,
                                                                        MulSIMD(FDPart3tmp308, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp15),
                                                                                                       MulSIMD(aDD11, alpha)))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_5 = MulSIMD(
            FDPart3tmp38,
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp106, PI),
                MulSIMD(
                    MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                    MulSIMD(alpha, FusedMulAddSIMD(
                                       FDPart3tmp176, FDPart3tmp326,
                                       FusedMulAddSIMD(FDPart3tmp155, FDPart3tmp326,
                                                       FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp326,
                                                                       FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp326,
                                                                                       FusedMulAddSIMD(FDPart3tmp191, FDPart3tmp326,
                                                                                                       FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp326,
                                                                                                                       FDPart3tmp139)))))))),
                FusedMulAddSIMD(
                    f0_of_xx0,
                    MulSIMD(MulSIMD(FDPart3tmp308, aDD12), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f3_of_xx2, alpha))),
                    FusedMulAddSIMD(
                        FDPart3tmp34, MulSIMD(FDPart3tmp6, aDD_dupD221),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp310, FDPart3tmp34), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD22, alpha)),
                            FusedMulAddSIMD(
                                FDPart3tmp13, MulSIMD(FDPart3tmp34, aDD_dupD220),
                                FusedMulAddSIMD(
                                    FDPart3tmp306, MulSIMD(FDPart3tmp35, aDD22),
                                    FusedMulAddSIMD(
                                        FDPart3tmp122, MulSIMD(aDD02, vetU_dD02),
                                        FusedMulAddSIMD(
                                            FDPart3tmp122, MulSIMD(aDD12, vetU_dD12),
                                            FusedMulAddSIMD(
                                                FDPart3tmp307, FDPart3tmp4,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp255, aDD22, MulSIMD(FDPart3tmp34, aDD_dupD222)),
                                                    FusedMulAddSIMD(
                                                        f0_of_xx0__D0,
                                                        MulSIMD(MulSIMD(FDPart3tmp311, aDD02),
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f3_of_xx2, alpha))),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp192,
                                                            NegFusedMulAddSIMD(FDPart3tmp284,
                                                                               FusedMulAddSIMD(FDPart3_Rational_1_3, FDPart3tmp34,
                                                                                               MulSIMD(FDPart3_Rational_1_3, FDPart3tmp47)),
                                                                               FDPart3tmp268),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp34),
                                                                    MulSIMD(FDPart3tmp42, aDD22)))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_6 = FusedMulAddSIMD(
            FDPart3tmp6, alpha_dupD1,
            FusedMulAddSIMD(FDPart3tmp9, alpha_dupD2, FusedMulSubSIMD(FDPart3tmp13, alpha_dupD0, MulSIMD(FDPart3_Integer_2, MulSIMD(alpha, trK)))));
const REAL_SIMD_ARRAY __RHS_exp_7 = MulSIMD(f0_of_xx0__D0, FusedMulAddSIMD(FDPart3tmp5, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp399), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp395, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp391), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))), FusedMulAddSIMD(FDPart3tmp5, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp399), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1))), FusedMulAddSIMD(FDPart3tmp96, MulSIMD(MulSIMD(FDPart3tmp366, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp391), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))), FusedMulAddSIMD(FDPart3tmp67, MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp367), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp84, MulSIMD(MulSIMD(FDPart3tmp364, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp403), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)), FusedMulAddSIMD(MulSIMD(FDPart3tmp38, FDPart3tmp406), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2)), FusedMulAddSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU_dD01)), FusedMulAddSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU_dD02)), FusedMulAddSIMD(FDPart3tmp12, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(FDPart3tmp26, lambdaU0)), FusedMulAddSIMD(MulSIMD(FDPart3tmp1, FDPart3tmp408), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)), FusedMulAddSIMD(FDPart3tmp406, MulSIMD(FDPart3tmp435, FDPart3tmp440), NegFusedMulAddSIMD(FDPart3tmp104, MulSIMD(alpha, trK_dD0), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp332, vetU1, FusedMulAddSIMD(FDPart3tmp6, FusedMulAddSIMD(FDPart3tmp1, FDPart3tmp271, FDPart3_NegativeOne_), AddSIMD(FDPart3tmp348, FDPart3tmp6)))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp88, FusedMulAddSIMD(FDPart3tmp26, FDPart3tmp346, FusedMulAddSIMD(FDPart3tmp271, FDPart3tmp33, NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp12, vetU_dD11), FusedMulSubSIMD(FDPart3tmp12, vetU_dDD011, FDPart3tmp13))))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp448, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD1, cf), alpha_dD1)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp450, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD2, cf), alpha_dD2)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp348, FDPart3tmp85), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp415, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD0, cf), alpha_dD0)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp100, FusedMulSubSIMD(FDPart3tmp12, vetU_dDD022, MulSIMD(FDPart3tmp312, MulSIMD(f3_of_xx2__D2, vetU_dD02)))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp13, FusedMulSubSIMD(FDPart3tmp12, f0_of_xx0__DDD000, MulSIMD(FDPart3tmp1, MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00))), FusedMulAddSIMD(FDPart3tmp12, MulSIMD(FDPart3tmp26, f0_of_xx0__DD00), FDPart3tmp359))), FusedMulAddSIMD(FDPart3tmp419, FDPart3tmp436, FusedMulAddSIMD(FDPart3tmp443, FDPart3tmp447, FusedMulAddSIMD(FDPart3tmp406, FDPart3tmp407, FusedMulAddSIMD(FDPart3tmp408, FDPart3tmp409, FusedMulAddSIMD(FDPart3tmp399, FDPart3tmp467, FusedMulAddSIMD(FDPart3tmp403, FDPart3tmp405, FusedMulAddSIMD(FDPart3tmp399, FDPart3tmp400, FusedMulAddSIMD(FDPart3tmp399, FDPart3tmp401, FusedMulAddSIMD(FDPart3tmp395, FDPart3tmp398, FusedMulAddSIMD(FDPart3tmp395, FDPart3tmp465, FusedMulAddSIMD(FDPart3tmp391, FDPart3tmp456, FusedMulAddSIMD(FDPart3tmp395, FDPart3tmp397, FusedMulAddSIMD(FDPart3tmp391, FDPart3tmp392, FusedMulAddSIMD(FDPart3tmp391, FDPart3tmp394, FusedMulAddSIMD(FDPart3tmp355, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp379, FDPart3tmp97, FusedMulAddSIMD(FDPart3tmp336, FDPart3tmp337, FusedMulAddSIMD(FDPart3tmp339, FDPart3tmp340, FusedMulAddSIMD(FDPart3tmp334, trK_dD1, FusedMulAddSIMD(FDPart3tmp335, FDPart3tmp97, FusedMulAddSIMD(FDPart3tmp329, betU_dupD01, FusedMulAddSIMD(FDPart3tmp330, betU_dupD02, FusedMulAddSIMD(FDPart3tmp104, FDPart3tmp385, FusedMulAddSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp12, betU_dupD00, MulSIMD(FDPart3tmp332, betU0)), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp395, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp474, MulSIMD(FDPart3tmp327, eta))))))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_8 = MulSIMD(f0_of_xx0, FusedMulAddSIMD(FDPart3tmp5, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp496), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp495, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp494), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))), FusedMulAddSIMD(FDPart3tmp5, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp496), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1))), FusedMulAddSIMD(FDPart3tmp84, MulSIMD(MulSIMD(FDPart3tmp367, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp494), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))), FusedMulAddSIMD(FDPart3tmp87, MulSIMD(MulSIMD(FDPart3tmp364, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp74, MulSIMD(MulSIMD(FDPart3tmp366, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp38, FDPart3tmp498), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2)), FusedMulAddSIMD(MulSIMD(FDPart3tmp5, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU_dD12)), FusedMulAddSIMD(MulSIMD(FDPart3tmp1, FDPart3tmp499), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)), FusedMulAddSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp497), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)), FusedMulAddSIMD(FDPart3tmp12, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(FDPart3tmp19, lambdaU0)), FusedMulAddSIMD(FDPart3tmp16, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU_dD11)), FusedMulAddSIMD(FDPart3tmp337, MulSIMD(FDPart3tmp5, vetU_dDD102), FusedMulAddSIMD(FDPart3tmp435, MulSIMD(FDPart3tmp440, FDPart3tmp498), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp500, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD1, cf), alpha_dD1)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp502, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD2, cf), alpha_dD2)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp104, FusedMulAddSIMD(FDPart3tmp19, MulSIMD(FDPart3tmp2, FDPart3tmp5), FusedMulAddSIMD(FDPart3tmp5, vetU_dDD100, FusedMulAddSIMD(vetU1, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp488, MulSIMD(FDPart3tmp16, f0_of_xx0__DD00)), FusedMulAddSIMD(FDPart3tmp485, FDPart3tmp6, FusedMulAddSIMD(FDPart3tmp488, vetU1, FusedMulAddSIMD(FDPart3tmp16, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, vetU_dD10)), FusedMulSubSIMD(FDPart3tmp19, FDPart3tmp442, MulSIMD(FDPart3tmp404, f0_of_xx0__DD00))))))))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp448, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD0, cf), alpha_dD0)), FusedMulAddSIMD(FDPart3tmp477, FDPart3tmp88, FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp100, FusedMulSubSIMD(FDPart3tmp5, vetU_dDD122, MulSIMD(FDPart3tmp322, MulSIMD(f3_of_xx2__D2, vetU_dD12)))), FusedMulAddSIMD(FDPart3tmp467, FDPart3tmp503, FusedMulAddSIMD(FDPart3tmp475, betU_dupD12, FusedMulAddSIMD(FDPart3tmp456, FDPart3tmp494, FusedMulAddSIMD(FDPart3tmp465, FDPart3tmp495, FusedMulAddSIMD(FDPart3tmp436, FDPart3tmp497, FusedMulAddSIMD(FDPart3tmp447, FDPart3tmp499, FusedMulAddSIMD(FDPart3tmp407, FDPart3tmp498, FusedMulAddSIMD(FDPart3tmp409, FDPart3tmp499, FusedMulAddSIMD(FDPart3tmp404, betU_dupD11, FusedMulAddSIMD(FDPart3tmp405, FDPart3tmp497, FusedMulAddSIMD(FDPart3tmp400, FDPart3tmp496, FusedMulAddSIMD(FDPart3tmp401, FDPart3tmp496, FusedMulAddSIMD(FDPart3tmp397, FDPart3tmp495, FusedMulAddSIMD(FDPart3tmp398, FDPart3tmp495, FusedMulAddSIMD(FDPart3tmp392, FDPart3tmp494, FusedMulAddSIMD(FDPart3tmp394, FDPart3tmp494, FusedMulAddSIMD(FDPart3tmp379, FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp385, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp340, FDPart3tmp478, FusedMulAddSIMD(FDPart3tmp355, FDPart3tmp88, FusedMulAddSIMD(FDPart3tmp334, trK_dD0, FusedMulAddSIMD(FDPart3tmp335, FDPart3tmp75, FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp504, FusedMulAddSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp382, betU1, MulSIMD(FDPart3tmp5, betU_dupD10)), FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp484, FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp486, FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp495, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp481, MulSIMD(FDPart3tmp396, eta))))))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_9 = MulSIMD(f3_of_xx2, FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp5, FDPart3tmp518), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp517), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp517), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))), FusedMulAddSIMD(FDPart3tmp519, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1))), FusedMulAddSIMD(FDPart3tmp519, MulSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0))), FusedMulAddSIMD(FDPart3tmp99, MulSIMD(MulSIMD(FDPart3tmp366, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp96, MulSIMD(MulSIMD(FDPart3tmp367, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp5, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU_dD21)), FusedMulAddSIMD(FDPart3tmp74, MulSIMD(MulSIMD(FDPart3tmp364, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp520), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)), FusedMulAddSIMD(MulSIMD(FDPart3tmp38, FDPart3tmp521), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2)), FusedMulAddSIMD(MulSIMD(FDPart3tmp1, FDPart3tmp522), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)), FusedMulAddSIMD(MulSIMD(FDPart3tmp12, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU_dD20)), NegFusedMulAddSIMD(FDPart3tmp97, MulSIMD(alpha, trK_dD0), FusedMulAddSIMD(FDPart3tmp306, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(FDPart3tmp8, lambdaU2)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp523, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD2, cf), alpha_dD2)), FusedMulAddSIMD(FDPart3tmp435, MulSIMD(FDPart3tmp440, FDPart3tmp524), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp450, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD0, cf), alpha_dD0)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp502, FusedMulAddSIMD(FDPart3tmp416, DivSIMD(cf_dD1, cf), alpha_dD1)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp100, FusedMulAddSIMD(FDPart3tmp9, FusedMulSubSIMD(FDPart3tmp8, f3_of_xx2__DD22, MulSIMD(FDPart3tmp370, FDPart3tmp38)), FusedMulAddSIMD(FDPart3tmp306, MulSIMD(FDPart3tmp8, f3_of_xx2__D2), FDPart3tmp377))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp104, FusedMulSubSIMD(FDPart3tmp8, vetU_dDD200, MulSIMD(FDPart3tmp312, MulSIMD(f0_of_xx0__DD00, vetU_dD20)))), FusedMulAddSIMD(FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp304, betU2, MulSIMD(FDPart3tmp8, betU_dupD22)), FusedMulAddSIMD(FDPart3_Rational_3_2, MulSIMD(FDPart3tmp508, FDPart3tmp85), FusedMulAddSIMD(FDPart3tmp505, betU_dupD20, FusedMulAddSIMD(FDPart3tmp506, betU_dupD21, FusedMulAddSIMD(FDPart3tmp467, FDPart3tmp519, FusedMulAddSIMD(FDPart3tmp477, FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp456, FDPart3tmp517, FusedMulAddSIMD(FDPart3tmp465, FDPart3tmp518, FusedMulAddSIMD(FDPart3tmp436, FDPart3tmp520, FusedMulAddSIMD(FDPart3tmp447, FDPart3tmp522, FusedMulAddSIMD(FDPart3tmp407, FDPart3tmp521, FusedMulAddSIMD(FDPart3tmp409, FDPart3tmp522, FusedMulAddSIMD(FDPart3tmp401, FDPart3tmp519, FusedMulAddSIMD(FDPart3tmp405, FDPart3tmp520, FusedMulAddSIMD(FDPart3tmp398, FDPart3tmp518, FusedMulAddSIMD(FDPart3tmp400, FDPart3tmp519, FusedMulAddSIMD(FDPart3tmp394, FDPart3tmp517, FusedMulAddSIMD(FDPart3tmp397, FDPart3tmp518, FusedMulAddSIMD(FDPart3tmp385, FDPart3tmp97, FusedMulAddSIMD(FDPart3tmp392, FDPart3tmp517, FusedMulAddSIMD(FDPart3tmp340, FDPart3tmp353, FusedMulAddSIMD(FDPart3tmp355, FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp100, FDPart3tmp379, FusedMulAddSIMD(FDPart3tmp337, FDPart3tmp380, FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp525, FusedMulAddSIMD(FDPart3tmp100, FDPart3tmp335, FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp5, FDPart3tmp518), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp509, MulSIMD(FDPart3tmp393, eta)))))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_10 = MulSIMD(
    MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
    MulSIMD(
        cf,
        FusedMulAddSIMD(
            MulSIMD(FDPart3tmp12, FDPart3tmp193), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(cf_dupD0, vetU0)),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp193, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(cf_dupD1, vetU1)),
                FusedMulAddSIMD(
                    FDPart3tmp1, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_6), MulSIMD(f0_of_xx0__DD00, vetU0)),
                    FusedMulAddSIMD(
                        FDPart3tmp33, MulSIMD(MulSIMD(FDPart3_Rational_1_12, FDPart3tmp16), MulSIMD(FDPart3tmp37, FDPart3tmp38)),
                        FusedMulAddSIMD(
                            FDPart3_Rational_1_6, FDPart3tmp29,
                            FusedMulAddSIMD(FDPart3_Rational_1_6, FDPart3tmp30,
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp193, FDPart3tmp8),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(cf_dupD2, vetU2)),
                                                            FusedMulSubSIMD(FDPart3_Rational_1_6, FDPart3tmp22,
                                                                            MulSIMD(FDPart3_Rational_1_6, MulSIMD(alpha, trK))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_11 = MulSIMD(
    FDPart3tmp1,
    FusedMulAddSIMD(
        FDPart3tmp19, MulSIMD(FDPart3tmp20, hDD01),
        FusedMulAddSIMD(
            FDPart3tmp0, MulSIMD(FDPart3tmp6, hDD_dupD001),
            FusedMulAddSIMD(
                FDPart3tmp0, MulSIMD(FDPart3tmp9, hDD_dupD002),
                FusedMulAddSIMD(
                    FDPart3tmp527, FusedMulAddSIMD(FDPart3_Rational_2_3, FDPart3tmp0, MulSIMD(FDPart3_Rational_2_3, FDPart3tmp53)),
                    FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp26, FDPart3tmp54),
                                    FusedMulAddSIMD(FDPart3tmp2, MulSIMD(hDD02, vetU_dD20),
                                                    FusedMulSubSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp0, hDD_dupD000, FDPart3tmp275),
                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp0), MulSIMD(aDD00, alpha))))))))));
const REAL_SIMD_ARRAY __RHS_exp_12 = MulSIMD(
    FDPart3tmp299,
    FusedMulAddSIMD(
        FDPart3tmp26, MulSIMD(FDPart3tmp78, hDD01),
        FusedMulAddSIMD(
            FDPart3tmp78, MulSIMD(FDPart3tmp9, hDD_dupD012),
            FusedMulAddSIMD(
                FDPart3tmp529, FDPart3tmp81,
                FusedMulAddSIMD(
                    FDPart3tmp12, MulSIMD(FDPart3tmp54, vetU_dD01),
                    FusedMulAddSIMD(
                        FDPart3tmp146, vetU_dD21,
                        FusedMulAddSIMD(
                            FDPart3tmp19, FDPart3tmp60,
                            FusedMulAddSIMD(
                                FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp78, hDD_dupD010, FDPart3tmp273),
                                FusedMulAddSIMD(FDPart3tmp14, hDD_dupD011,
                                                FusedMulAddSIMD(f0_of_xx0__D0, MulSIMD(hDD01, vetU_dD11),
                                                                FusedMulSubSIMD(FDPart3tmp123, vetU_dD20,
                                                                                MulSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, aDD01),
                                                                                                           MulSIMD(f0_of_xx0__D0, alpha))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_13 = MulSIMD(
    FDPart3tmp312,
    FusedMulAddSIMD(
        FDPart3tmp48, MulSIMD(FDPart3tmp8, vetU_dD20),
        FusedMulAddSIMD(
            FDPart3tmp6, MulSIMD(FDPart3tmp76, hDD_dupD021),
            FusedMulAddSIMD(
                FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp146, f3_of_xx2__D2, MulSIMD(FDPart3tmp76, hDD_dupD022)),
                FusedMulAddSIMD(
                    FDPart3tmp12, MulSIMD(FDPart3tmp54, vetU_dD02),
                    FusedMulAddSIMD(
                        FDPart3tmp306, FDPart3tmp93,
                        FusedMulAddSIMD(
                            FDPart3tmp529, FDPart3tmp93,
                            FusedMulAddSIMD(
                                FDPart3tmp19, FDPart3tmp71,
                                FusedMulAddSIMD(
                                    FDPart3tmp26, FDPart3tmp93,
                                    FusedMulAddSIMD(
                                        f0_of_xx0__D0, MulSIMD(hDD01, vetU_dD12),
                                        FusedMulSubSIMD(
                                            FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp118, f0_of_xx0__DD00, MulSIMD(FDPart3tmp76, hDD_dupD020)),
                                            MulSIMD(f0_of_xx0__D0, MulSIMD(MulSIMD(FDPart3_Integer_2, aDD02), MulSIMD(f3_of_xx2, alpha))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_14 = MulSIMD(
    FDPart3tmp16,
    FusedMulAddSIMD(
        FDPart3tmp15, MulSIMD(FDPart3tmp9, hDD_dupD112),
        FusedMulAddSIMD(
            FDPart3_Integer_2, MulSIMD(FDPart3tmp123, vetU_dD21),
            FusedMulAddSIMD(
                FDPart3_Integer_2, MulSIMD(FDPart3tmp30, FDPart3tmp60),
                FusedMulAddSIMD(
                    FDPart3tmp315, hDD01,
                    FusedMulAddSIMD(FDPart3tmp527, FusedMulAddSIMD(FDPart3_Rational_2_3, FDPart3tmp15, MulSIMD(FDPart3_Rational_2_3, FDPart3tmp59)),
                                    FusedMulAddSIMD(f0_of_xx0, MulSIMD(hDD_dupD111, vetU1),
                                                    FusedMulSubSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp15, hDD_dupD110, FDPart3tmp232),
                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp15), MulSIMD(aDD11, alpha))))))))));
const REAL_SIMD_ARRAY __RHS_exp_15 = MulSIMD(
    FDPart3tmp322,
    FusedMulAddSIMD(
        FDPart3tmp5, MulSIMD(FDPart3tmp60, vetU_dD12),
        FusedMulAddSIMD(
            f0_of_xx0, MulSIMD(hDD01, vetU_dD02),
            FusedMulAddSIMD(
                FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp123, f3_of_xx2__D2, MulSIMD(FDPart3tmp68, hDD_dupD122)),
                FusedMulAddSIMD(
                    FDPart3tmp48, MulSIMD(FDPart3tmp8, vetU_dD21),
                    FusedMulAddSIMD(
                        FDPart3tmp306, FDPart3tmp71,
                        FusedMulAddSIMD(
                            FDPart3tmp529, FDPart3tmp71,
                            FusedMulAddSIMD(
                                FDPart3tmp129, hDD_dupD121,
                                FusedMulAddSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp68, hDD_dupD120, FDPart3tmp90),
                                                FusedMulAddSIMD(f3_of_xx2, MulSIMD(hDD12, vetU_dD11),
                                                                FusedMulSubSIMD(FDPart3tmp118, vetU_dD01,
                                                                                MulSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, aDD12),
                                                                                                           MulSIMD(f3_of_xx2, alpha))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_16 = MulSIMD(
    FDPart3tmp38,
    FusedMulAddSIMD(
        FDPart3tmp13, MulSIMD(FDPart3tmp34, hDD_dupD220),
        FusedMulAddSIMD(
            FDPart3_Integer_2, MulSIMD(FDPart3tmp306, FDPart3tmp48),
            FusedMulAddSIMD(
                FDPart3tmp122, MulSIMD(hDD12, vetU_dD12),
                FusedMulAddSIMD(
                    FDPart3tmp527, FusedMulAddSIMD(FDPart3_Rational_2_3, FDPart3tmp34, MulSIMD(FDPart3_Rational_2_3, FDPart3tmp47)),
                    FusedMulAddSIMD(FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp34, hDD_dupD222, FDPart3tmp256),
                                    FusedMulAddSIMD(FDPart3tmp34, MulSIMD(FDPart3tmp6, hDD_dupD221),
                                                    FusedMulSubSIMD(FDPart3tmp141, vetU_dD02,
                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp34), MulSIMD(aDD22, alpha))))))))));
const REAL_SIMD_ARRAY __RHS_exp_17 = MulSIMD(
    f0_of_xx0__D0,
    FusedMulAddSIMD(
        FDPart3tmp84, MulSIMD(MulSIMD(FDPart3tmp364, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
        FusedMulAddSIMD(
            MulSIMD(FDPart3tmp67, FDPart3tmp96), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD2, alpha)),
            FusedMulAddSIMD(
                FDPart3tmp67,
                MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp367), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp103, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD0, alpha)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp67, FDPart3tmp84), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD1, alpha)),
                        FusedMulAddSIMD(
                            FDPart3tmp406, MulSIMD(FDPart3tmp440, FDPart3tmp540),
                            NegFusedMulAddSIMD(
                                FDPart3tmp312,
                                MulSIMD(lambdaU2, vetU_dD02),
                                FusedMulAddSIMD(
                                    FDPart3tmp85,
                                    FusedMulAddSIMD(FDPart3tmp332, vetU1,
                                                    FusedMulAddSIMD(FDPart3tmp6, FusedMulAddSIMD(FDPart3tmp1, FDPart3tmp271, FDPart3_NegativeOne_),
                                                                    AddSIMD(FDPart3tmp348, FDPart3tmp6))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp88,
                                        FusedMulAddSIMD(
                                            FDPart3tmp26, FDPart3tmp346,
                                            FusedMulAddSIMD(FDPart3tmp271, FDPart3tmp33,
                                                            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp12, vetU_dD11),
                                                                               FusedMulSubSIMD(FDPart3tmp12, vetU_dDD011, FDPart3tmp13)))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp534, FDPart3tmp85,
                                            FusedMulAddSIMD(
                                                FDPart3tmp536, FDPart3tmp97,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp448, FDPart3tmp543,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp450, FDPart3tmp544,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp419, FDPart3tmp541,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp443, FDPart3tmp542,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp399, FDPart3tmp548,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp415, FDPart3tmp539,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp391, FDPart3tmp546,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp395, FDPart3tmp547,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp330, lambdaU_dupD02,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp348, FDPart3tmp85,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp26, FDPart3tmp532,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp329, lambdaU_dupD01,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp185, FDPart3tmp339,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp190, FDPart3tmp336,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp104,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp13,
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3tmp12, f0_of_xx0__DDD000,
                                                                                                                    MulSIMD(
                                                                                                                        FDPart3tmp1,
                                                                                                                        MulSIMD(f0_of_xx0__DD00,
                                                                                                                                f0_of_xx0__DD00))),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp12,
                                                                                                                    MulSIMD(FDPart3tmp26,
                                                                                                                            f0_of_xx0__DD00),
                                                                                                                    FDPart3tmp359)),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp13,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp12, lambdaU_dupD00,
                                                                                                                    MulSIMD(FDPart3tmp332, lambdaU0)),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp100,
                                                                                                                    FusedMulSubSIMD(
                                                                                                                        FDPart3tmp12, vetU_dDD022,
                                                                                                                        MulSIMD(FDPart3tmp312,
                                                                                                                                MulSIMD(f3_of_xx2__D2,
                                                                                                                                        vetU_dD02))),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp104, FDPart3tmp537,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp96,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(FDPart3tmp366,
                                                                                                                                        FDPart3tmp67),
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Integer_16,
                                                                                                                                        FDPart3_NegativeOne_),
                                                                                                                                    MulSIMD(PI,
                                                                                                                                            alpha))),
                                                                                                                            NegFusedMulAddSIMD(
                                                                                                                                FDPart3tmp299,
                                                                                                                                MulSIMD(lambdaU1,
                                                                                                                                        vetU_dD01),
                                                                                                                                FDPart3tmp474))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_18 = MulSIMD(
    f0_of_xx0,
    FusedMulAddSIMD(
        FDPart3tmp74, MulSIMD(MulSIMD(FDPart3tmp366, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
        FusedMulAddSIMD(
            MulSIMD(FDPart3tmp67, FDPart3tmp87), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD1, alpha)),
            FusedMulAddSIMD(
                FDPart3tmp87,
                MulSIMD(MulSIMD(FDPart3tmp364, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp67, FDPart3tmp74), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD2, alpha)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp67, FDPart3tmp84), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD0, alpha)),
                        FusedMulAddSIMD(
                            FDPart3tmp440, MulSIMD(FDPart3tmp498, FDPart3tmp540),
                            NegFusedMulAddSIMD(
                                FDPart3tmp322,
                                MulSIMD(lambdaU2, vetU_dD12),
                                FusedMulAddSIMD(
                                    FDPart3tmp537,
                                    FDPart3tmp85,
                                    FusedMulAddSIMD(
                                        FDPart3tmp190,
                                        MulSIMD(FDPart3tmp5, vetU_dDD102),
                                        FusedMulAddSIMD(
                                            FDPart3tmp534,
                                            FDPart3tmp88,
                                            FusedMulAddSIMD(
                                                FDPart3tmp536,
                                                FDPart3tmp75,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp502,
                                                    FDPart3tmp544,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp503,
                                                        FDPart3tmp548,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp499,
                                                            FDPart3tmp542,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp500,
                                                                FDPart3tmp543,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp495,
                                                                    FDPart3tmp547,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp497,
                                                                        FDPart3tmp541,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp475,
                                                                            lambdaU_dupD12,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp494,
                                                                                FDPart3tmp546,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp404,
                                                                                    lambdaU_dupD11,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp448,
                                                                                        FDPart3tmp539,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp185,
                                                                                            FDPart3tmp478,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp19, FDPart3tmp532,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp104,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp19,
                                                                                                        MulSIMD(FDPart3tmp2, FDPart3tmp5),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp5, vetU_dDD100,
                                                                                                            FusedMulAddSIMD(
                                                                                                                vetU1,
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3_Integer_2, FDPart3tmp488,
                                                                                                                    MulSIMD(FDPart3tmp16,
                                                                                                                            f0_of_xx0__DD00)),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp485, FDPart3tmp6,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp488, vetU1,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp16,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                    FDPart3_NegativeOne_),
                                                                                                                                MulSIMD(f0_of_xx0__D0,
                                                                                                                                        vetU_dD10)),
                                                                                                                            FusedMulSubSIMD(
                                                                                                                                FDPart3tmp19,
                                                                                                                                FDPart3tmp442,
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3tmp404,
                                                                                                                                    f0_of_xx0__DD00)))))))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp13,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp382, lambdaU1,
                                                                                                            MulSIMD(FDPart3tmp5, lambdaU_dupD10)),
                                                                                                        AddSIMD(
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp100,
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3tmp5, vetU_dDD122,
                                                                                                                    MulSIMD(FDPart3tmp322,
                                                                                                                            MulSIMD(f3_of_xx2__D2,
                                                                                                                                    vetU_dD12))),
                                                                                                                FDPart3tmp504),
                                                                                                            AddSIMD(
                                                                                                                AddSIMD(FDPart3tmp484, FDPart3tmp486),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp84,
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp367,
                                                                                                                                FDPart3tmp67),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_16,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                            MulSIMD(PI, alpha))),
                                                                                                                    NegFusedMulAddSIMD(
                                                                                                                        FDPart3tmp16,
                                                                                                                        MulSIMD(lambdaU1, vetU_dD11),
                                                                                                                        FDPart3tmp481))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_19 = MulSIMD(
    f3_of_xx2,
    FusedMulAddSIMD(
        FDPart3tmp74, MulSIMD(MulSIMD(FDPart3tmp364, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
        FusedMulAddSIMD(
            FDPart3tmp99,
            MulSIMD(MulSIMD(FDPart3tmp366, FDPart3tmp67), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp67, FDPart3tmp96), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD0, alpha)),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp67, FDPart3tmp99), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD2, alpha)),
                    NegFusedMulAddSIMD(
                        FDPart3tmp312,
                        MulSIMD(lambdaU0, vetU_dD20),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp67, FDPart3tmp74),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD1, alpha)),
                            FusedMulAddSIMD(
                                FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp304, lambdaU2, MulSIMD(FDPart3tmp8, lambdaU_dupD22)),
                                FusedMulAddSIMD(
                                    FDPart3tmp440,
                                    MulSIMD(FDPart3tmp524, FDPart3tmp540),
                                    FusedMulAddSIMD(
                                        FDPart3tmp534,
                                        FDPart3tmp75,
                                        FusedMulAddSIMD(
                                            FDPart3tmp537,
                                            FDPart3tmp97,
                                            FusedMulAddSIMD(
                                                FDPart3tmp522,
                                                FDPart3tmp542,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp523,
                                                    FDPart3tmp544,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp519,
                                                        FDPart3tmp548,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp520,
                                                            FDPart3tmp541,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp517,
                                                                FDPart3tmp546,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp518,
                                                                    FDPart3tmp547,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp506,
                                                                        lambdaU_dupD21,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp507,
                                                                            lambdaU1,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp502,
                                                                                FDPart3tmp543,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp505,
                                                                                    lambdaU_dupD20,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp190,
                                                                                        FDPart3tmp380,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp450,
                                                                                            FDPart3tmp539,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp175,
                                                                                                FDPart3tmp508,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp185,
                                                                                                    FDPart3tmp353,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp100,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp9,
                                                                                                            FusedMulSubSIMD(
                                                                                                                FDPart3tmp8, f3_of_xx2__DD22,
                                                                                                                MulSIMD(FDPart3tmp370, FDPart3tmp38)),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp306,
                                                                                                                MulSIMD(FDPart3tmp8, f3_of_xx2__D2),
                                                                                                                FDPart3tmp377)),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp104,
                                                                                                            FusedMulSubSIMD(
                                                                                                                FDPart3tmp8, vetU_dDD200,
                                                                                                                MulSIMD(FDPart3tmp312,
                                                                                                                        MulSIMD(f0_of_xx0__DD00,
                                                                                                                                vetU_dD20))),
                                                                                                            AddSIMD(
                                                                                                                FusedMulAddSIMD(FDPart3tmp100,
                                                                                                                                FDPart3tmp536,
                                                                                                                                FDPart3tmp525),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp96,
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp367,
                                                                                                                                FDPart3tmp67),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_16,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                            MulSIMD(PI, alpha))),
                                                                                                                    NegFusedMulAddSIMD(
                                                                                                                        FDPart3tmp306,
                                                                                                                        MulSIMD(FDPart3tmp8,
                                                                                                                                lambdaU2),
                                                                                                                        FDPart3tmp509))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_20 = NegFusedMulAddSIMD(
    FDPart3tmp104,
    MulSIMD(FDPart3tmp192,
            FusedMulAddSIMD(FDPart3tmp277, alpha_dD0,
                            FusedMulAddSIMD(FDPart3tmp279, alpha_dD1,
                                            FusedMulAddSIMD(FDPart3tmp281, alpha_dD2, NegFusedMulAddSIMD(FDPart3tmp197, alpha_dD0, alpha_dDD00))))),
    FusedMulAddSIMD(
        FDPart3tmp431, MulSIMD(FDPart3tmp464, alpha),
        FusedMulAddSIMD(
            FDPart3tmp545,
            MulSIMD(PI,
                    FusedMulAddSIMD(
                        FDPart3tmp174, MulSIMD(FDPart3tmp559, FusedMulSubSIMD(FDPart3tmp110, FDPart3tmp80, MulSIMD(FDPart3tmp110, FDPart3tmp82))),
                        FusedMulAddSIMD(
                            FDPart3tmp154, MulSIMD(FDPart3tmp558, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp57, MulSIMD(FDPart3tmp110, FDPart3tmp64))),
                            FusedMulAddSIMD(
                                FDPart3tmp163, MulSIMD(FDPart3tmp558, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp62, FDPart3tmp187)),
                                FusedMulAddSIMD(
                                    FDPart3tmp183,
                                    MulSIMD(FDPart3tmp559, FusedMulSubSIMD(FDPart3tmp110, FDPart3tmp70, MulSIMD(FDPart3tmp110, FDPart3tmp72))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp189,
                                        MulSIMD(FDPart3tmp559, FusedMulSubSIMD(FDPart3tmp110, FDPart3tmp92, MulSIMD(FDPart3tmp110, FDPart3tmp94))),
                                        FusedMulAddSIMD(T4UU00, MulSIMD(alpha, alpha),
                                                        MulSIMD(FDPart3tmp139, MulSIMD(FDPart3tmp558, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp51,
                                                                                                                      FDPart3tmp167)))))))))),
            FusedMulAddSIMD(
                FDPart3tmp426, MulSIMD(FDPart3tmp455, alpha),
                FusedMulAddSIMD(
                    FDPart3tmp428, MulSIMD(FDPart3tmp466, alpha),
                    FusedMulAddSIMD(
                        FDPart3tmp3, MulSIMD(FDPart3tmp446, alpha),
                        FusedMulAddSIMD(
                            FDPart3tmp307, MulSIMD(FDPart3tmp440, alpha),
                            FusedMulAddSIMD(
                                FDPart3_Rational_1_3, MulSIMD(alpha, MulSIMD(trK, trK)),
                                FusedMulAddSIMD(
                                    FDPart3tmp292, MulSIMD(FDPart3tmp434, alpha),
                                    FusedMulAddSIMD(
                                        FDPart3tmp6, trK_dupD1,
                                        FusedMulAddSIMD(
                                            FDPart3tmp9, trK_dupD2,
                                            FusedMulAddSIMD(
                                                FDPart3tmp556, AddSIMD(FDPart3tmp229, FDPart3tmp555),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp556, AddSIMD(FDPart3tmp230, FDPart3tmp555),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp554, AddSIMD(FDPart3tmp216, FDPart3tmp553),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp554, AddSIMD(FDPart3tmp219, FDPart3tmp553),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp552, AddSIMD(FDPart3tmp196, FDPart3tmp550),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp552, AddSIMD(FDPart3tmp199, FDPart3tmp550),
                                                                    NegFusedMulAddSIMD(
                                                                        FDPart3tmp192,
                                                                        MulSIMD(FDPart3tmp88,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp247, alpha_dD0,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp249, alpha_dD1,
                                                                                        FusedMulAddSIMD(FDPart3tmp251, alpha_dD2,
                                                                                                        NegFusedMulAddSIMD(FDPart3tmp217, alpha_dD1,
                                                                                                                           alpha_dDD11))))),
                                                                        FusedMulSubSIMD(
                                                                            FDPart3tmp13, trK_dupD0,
                                                                            MulSIMD(FDPart3tmp100,
                                                                                    MulSIMD(FDPart3tmp192,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp262, alpha_dD0,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp264, alpha_dD1,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp266, alpha_dD2,
                                                                                                        NegFusedMulAddSIMD(
                                                                                                            FDPart3tmp194,
                                                                                                            alpha_dD2,
                                                                                                            alpha_dDD22)))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_21 = MulSIMD(
    f0_of_xx0__D0,
    FusedMulAddSIMD(
        FDPart3tmp399, FDPart3tmp563,
        FusedMulAddSIMD(
            FDPart3tmp403, FDPart3tmp564,
            FusedMulAddSIMD(
                FDPart3tmp391, FDPart3tmp561,
                FusedMulAddSIMD(
                    FDPart3tmp395, FDPart3tmp562,
                    FusedMulAddSIMD(
                        FDPart3tmp329, vetU_dupD01,
                        FusedMulAddSIMD(
                            FDPart3tmp330, vetU_dupD02,
                            FusedMulAddSIMD(FDPart3tmp406, FDPart3tmp565,
                                            FusedMulAddSIMD(FDPart3tmp408, FDPart3tmp566,
                                                            FusedMulAddSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp12, vetU_dupD00, FDPart3tmp25),
                                                                            FDPart3tmp327))))))))));
const REAL_SIMD_ARRAY __RHS_exp_22 = MulSIMD(
    f0_of_xx0,
    FusedMulAddSIMD(
        FDPart3tmp496, FDPart3tmp563,
        FusedMulAddSIMD(
            FDPart3tmp497, FDPart3tmp564,
            FusedMulAddSIMD(
                FDPart3tmp494, FDPart3tmp561,
                FusedMulAddSIMD(
                    FDPart3tmp495, FDPart3tmp562,
                    FusedMulAddSIMD(
                        FDPart3tmp404, vetU_dupD11,
                        FusedMulAddSIMD(
                            FDPart3tmp475, vetU_dupD12,
                            FusedMulAddSIMD(FDPart3tmp498, FDPart3tmp565,
                                            FusedMulAddSIMD(FDPart3tmp499, FDPart3tmp566,
                                                            FusedMulAddSIMD(FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp5, vetU_dupD10, FDPart3tmp18),
                                                                            FDPart3tmp396))))))))));
const REAL_SIMD_ARRAY __RHS_exp_23 = MulSIMD(
    f3_of_xx2,
    FusedMulAddSIMD(
        FDPart3tmp520, FDPart3tmp564,
        FusedMulAddSIMD(
            FDPart3tmp521, FDPart3tmp565,
            FusedMulAddSIMD(
                FDPart3tmp518, FDPart3tmp562,
                FusedMulAddSIMD(
                    FDPart3tmp519, FDPart3tmp563,
                    FusedMulAddSIMD(
                        FDPart3tmp506, vetU_dupD21,
                        FusedMulAddSIMD(
                            FDPart3tmp517, FDPart3tmp561,
                            FusedMulAddSIMD(FDPart3tmp522, FDPart3tmp566,
                                            FusedMulAddSIMD(FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp304, vetU2, MulSIMD(FDPart3tmp8, vetU_dupD22)),
                                                            FusedMulAddSIMD(FDPart3tmp505, vetU_dupD20, FDPart3tmp393))))))))));

WriteSIMD(&rhs_gfs[IDX4(ADD00GF, i0, i1, i2)], __RHS_exp_0);
WriteSIMD(&rhs_gfs[IDX4(ADD01GF, i0, i1, i2)], __RHS_exp_1);
WriteSIMD(&rhs_gfs[IDX4(ADD02GF, i0, i1, i2)], __RHS_exp_2);
WriteSIMD(&rhs_gfs[IDX4(ADD11GF, i0, i1, i2)], __RHS_exp_3);
WriteSIMD(&rhs_gfs[IDX4(ADD12GF, i0, i1, i2)], __RHS_exp_4);
WriteSIMD(&rhs_gfs[IDX4(ADD22GF, i0, i1, i2)], __RHS_exp_5);
WriteSIMD(&rhs_gfs[IDX4(ALPHAGF, i0, i1, i2)], __RHS_exp_6);
WriteSIMD(&rhs_gfs[IDX4(BETU0GF, i0, i1, i2)], __RHS_exp_7);
WriteSIMD(&rhs_gfs[IDX4(BETU1GF, i0, i1, i2)], __RHS_exp_8);
WriteSIMD(&rhs_gfs[IDX4(BETU2GF, i0, i1, i2)], __RHS_exp_9);
WriteSIMD(&rhs_gfs[IDX4(CFGF, i0, i1, i2)], __RHS_exp_10);
WriteSIMD(&rhs_gfs[IDX4(HDD00GF, i0, i1, i2)], __RHS_exp_11);
WriteSIMD(&rhs_gfs[IDX4(HDD01GF, i0, i1, i2)], __RHS_exp_12);
WriteSIMD(&rhs_gfs[IDX4(HDD02GF, i0, i1, i2)], __RHS_exp_13);
WriteSIMD(&rhs_gfs[IDX4(HDD11GF, i0, i1, i2)], __RHS_exp_14);
WriteSIMD(&rhs_gfs[IDX4(HDD12GF, i0, i1, i2)], __RHS_exp_15);
WriteSIMD(&rhs_gfs[IDX4(HDD22GF, i0, i1, i2)], __RHS_exp_16);
WriteSIMD(&rhs_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)], __RHS_exp_17);
WriteSIMD(&rhs_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)], __RHS_exp_18);
WriteSIMD(&rhs_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)], __RHS_exp_19);
WriteSIMD(&rhs_gfs[IDX4(TRKGF, i0, i1, i2)], __RHS_exp_20);
WriteSIMD(&rhs_gfs[IDX4(VETU0GF, i0, i1, i2)], __RHS_exp_21);
WriteSIMD(&rhs_gfs[IDX4(VETU1GF, i0, i1, i2)], __RHS_exp_22);
WriteSIMD(&rhs_gfs[IDX4(VETU2GF, i0, i1, i2)], __RHS_exp_23);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += simd_width)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
