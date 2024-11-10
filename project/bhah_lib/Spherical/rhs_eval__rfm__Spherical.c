#include "../BHaH_defines.h"
#include "../simd/simd_intrinsics.h"
/**
 * Finite difference function for operator dD0, with FD accuracy order 4.
 */
static NO_INLINE REAL_SIMD_ARRAY SIMD_fd_function_dD0_fdorder4(const REAL_SIMD_ARRAY FDPROTO_i0m1, const REAL_SIMD_ARRAY FDPROTO_i0m2,
                                                               const REAL_SIMD_ARRAY FDPROTO_i0p1, const REAL_SIMD_ARRAY FDPROTO_i0p2,
                                                               const REAL_SIMD_ARRAY invdxx0) {
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_2_3 = 2.0 / 3.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_2_3 = 2.0 / 3.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_2_3 = 2.0 / 3.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_3 = ConstSIMD(dblFDPart1_Rational_4_3);

  const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
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
  const double dblFDPart1_Rational_1_144 = 1.0 / 144.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_144 = ConstSIMD(dblFDPart1_Rational_1_144);

  const double dblFDPart1_Rational_1_18 = 1.0 / 18.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_18 = ConstSIMD(dblFDPart1_Rational_1_18);

  const double dblFDPart1_Rational_4_9 = 4.0 / 9.0;
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
  const double dblFDPart1_Rational_1_144 = 1.0 / 144.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_144 = ConstSIMD(dblFDPart1_Rational_1_144);

  const double dblFDPart1_Rational_1_18 = 1.0 / 18.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_18 = ConstSIMD(dblFDPart1_Rational_1_18);

  const double dblFDPart1_Rational_4_9 = 4.0 / 9.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_3 = ConstSIMD(dblFDPart1_Rational_4_3);

  const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
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
  const double dblFDPart1_Rational_1_144 = 1.0 / 144.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_144 = ConstSIMD(dblFDPart1_Rational_1_144);

  const double dblFDPart1_Rational_1_18 = 1.0 / 18.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_18 = ConstSIMD(dblFDPart1_Rational_1_18);

  const double dblFDPart1_Rational_4_9 = 4.0 / 9.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_4_3 = 4.0 / 3.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_4_3 = ConstSIMD(dblFDPart1_Rational_4_3);

  const double dblFDPart1_Rational_5_2 = 5.0 / 2.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
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
  const double dblFDPart1_Rational_1_12 = 1.0 / 12.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_12 = ConstSIMD(dblFDPart1_Rational_1_12);

  const double dblFDPart1_Rational_1_2 = 1.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_2 = ConstSIMD(dblFDPart1_Rational_1_2);

  const double dblFDPart1_Rational_1_4 = 1.0 / 4.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_1_4 = ConstSIMD(dblFDPart1_Rational_1_4);

  const double dblFDPart1_Rational_3_2 = 3.0 / 2.0;
  const REAL_SIMD_ARRAY FDPart1_Rational_3_2 = ConstSIMD(dblFDPart1_Rational_3_2);

  const double dblFDPart1_Rational_5_6 = 5.0 / 6.0;
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
void rhs_eval__rfm__Spherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                              const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                              REAL *restrict rhs_gfs) {
#include "../set_CodeParameters-simd.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      const double NOSIMDf1_of_xx1 = rfmstruct->f1_of_xx1[i1];
      const REAL_SIMD_ARRAY f1_of_xx1 = ConstSIMD(NOSIMDf1_of_xx1);
      const double NOSIMDf1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
      const REAL_SIMD_ARRAY f1_of_xx1__D1 = ConstSIMD(NOSIMDf1_of_xx1__D1);
      const double NOSIMDf1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
      const REAL_SIMD_ARRAY f1_of_xx1__DD11 = ConstSIMD(NOSIMDf1_of_xx1__DD11);

      for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0 += simd_width) {
        const REAL_SIMD_ARRAY f0_of_xx0 = ReadSIMD(&rfmstruct->f0_of_xx0[i0]);
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
        const double dblFDPart1_Integer_1 = 1.0;
        const REAL_SIMD_ARRAY FDPart1_Integer_1 = ConstSIMD(dblFDPart1_Integer_1);

        const double dblFDPart1_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart1_NegativeOne_ = ConstSIMD(dblFDPart1_NegativeOne_);

        const REAL_SIMD_ARRAY FDPart1tmp0 = DivSIMD(FDPart1_Integer_1, f0_of_xx0);
        const REAL_SIMD_ARRAY UpwindControlVectorU0 = vetU0;
        const REAL_SIMD_ARRAY UpwindControlVectorU1 = MulSIMD(FDPart1tmp0, vetU1);
        const REAL_SIMD_ARRAY UpwindControlVectorU2 = MulSIMD(FDPart1tmp0, DivSIMD(vetU2, f1_of_xx1));

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
        const double dblFDPart2_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart2_NegativeOne_ = ConstSIMD(dblFDPart2_NegativeOne_);

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
        const double dblFDPart3_Integer_1 = 1.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        const double dblFDPart3_Integer_12 = 12.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_12 = ConstSIMD(dblFDPart3_Integer_12);

        const double dblFDPart3_Integer_16 = 16.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_16 = ConstSIMD(dblFDPart3_Integer_16);

        const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        const double dblFDPart3_Integer_3 = 3.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_3 = ConstSIMD(dblFDPart3_Integer_3);

        const double dblFDPart3_Integer_4 = 4.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_4 = ConstSIMD(dblFDPart3_Integer_4);

        const double dblFDPart3_Integer_8 = 8.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_8 = ConstSIMD(dblFDPart3_Integer_8);

        const double dblFDPart3_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        const double dblFDPart3_Rational_1_3 = 1.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_3 = ConstSIMD(dblFDPart3_Rational_1_3);

        const double dblFDPart3_Rational_1_4 = 1.0 / 4.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_4 = ConstSIMD(dblFDPart3_Rational_1_4);

        const double dblFDPart3_Rational_1_6 = 1.0 / 6.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_6 = ConstSIMD(dblFDPart3_Rational_1_6);

        const double dblFDPart3_Rational_2_3 = 2.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_2_3 = ConstSIMD(dblFDPart3_Rational_2_3);

        const double dblFDPart3_Rational_3_2 = 3.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_3_2 = ConstSIMD(dblFDPart3_Rational_3_2);

        const double dblFDPart3_Rational_3_4 = 3.0 / 4.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_3_4 = ConstSIMD(dblFDPart3_Rational_3_4);

        const double dblFDPart3_Rational_4_3 = 4.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_4_3 = ConstSIMD(dblFDPart3_Rational_4_3);

        const REAL_SIMD_ARRAY FDPart3tmp2 = DivSIMD(FDPart3_Integer_1, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp4 = DivSIMD(FDPart3_Integer_1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp7 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp12 = MulSIMD(FDPart3_Integer_2, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp17 = MulSIMD(aDD02, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp28 = MulSIMD(FDPart3_Integer_2, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp29 = MulSIMD(f1_of_xx1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp30 = MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp33 = AddSIMD(FDPart3_Integer_1, hDD00);
        const REAL_SIMD_ARRAY FDPart3tmp53 = MulSIMD(f1_of_xx1, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp60 = MulSIMD(aDD01, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp61 = MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp75 = MulSIMD(f0_of_xx0, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp89 = MulSIMD(cf, cf);
        const REAL_SIMD_ARRAY FDPart3tmp93 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf));
        const REAL_SIMD_ARRAY FDPart3tmp109 = MulSIMD(f0_of_xx0, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp169 = MulSIMD(cf, cf);
        const REAL_SIMD_ARRAY FDPart3tmp172 = DivSIMD(FDPart3_Integer_1, cf);
        const REAL_SIMD_ARRAY FDPart3tmp195 = MulSIMD(f0_of_xx0, hDD_dD012);
        const REAL_SIMD_ARRAY FDPart3tmp242 = MulSIMD(f1_of_xx1, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp267 = MulSIMD(f0_of_xx0, vetU_dD00);
        const REAL_SIMD_ARRAY FDPart3tmp268 = MulSIMD(alpha, trK);
        const REAL_SIMD_ARRAY FDPart3tmp328 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(T4UU00, alpha));
        const REAL_SIMD_ARRAY FDPart3tmp337 = MulSIMD(f1_of_xx1__D1, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp356 = MulSIMD(betU0, vetU0);
        const REAL_SIMD_ARRAY FDPart3tmp382 = MulSIMD(FDPart3_Integer_3, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp398 = MulSIMD(FDPart3_Rational_3_2, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp503 = MulSIMD(FDPart3_Integer_3, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp505 = MulSIMD(FDPart3_Integer_2, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp510 = MulSIMD(FDPart3_Integer_4, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp523 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf), cf), cf));
        const REAL_SIMD_ARRAY FDPart3tmp526 = MulSIMD(vetU0, vetU0);
        const REAL_SIMD_ARRAY FDPart3tmp3 = MulSIMD(FDPart3tmp2, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp5 = MulSIMD(FDPart3tmp2, FDPart3tmp4);
        const REAL_SIMD_ARRAY FDPart3tmp8 = DivSIMD(FDPart3_Integer_1, FDPart3tmp7);
        const REAL_SIMD_ARRAY FDPart3tmp13 = MulSIMD(FDPart3tmp12, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp14 = MulSIMD(FDPart3tmp4, vetU2);
        const REAL_SIMD_ARRAY FDPart3tmp18 = MulSIMD(FDPart3tmp12, FDPart3tmp17);
        const REAL_SIMD_ARRAY FDPart3tmp20 = MulSIMD(FDPart3tmp2, vetU_dD11);
        const REAL_SIMD_ARRAY FDPart3tmp21 = MulSIMD(FDPart3_Integer_2, FDPart3tmp2);
        const REAL_SIMD_ARRAY FDPart3tmp37 = MulSIMD(FDPart3tmp7, hDD11);
        const REAL_SIMD_ARRAY FDPart3tmp39 = MulSIMD(FDPart3tmp29, FDPart3tmp7);
        const REAL_SIMD_ARRAY FDPart3tmp65 = MulSIMD(FDPart3tmp53, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp73 = MulSIMD(FDPart3tmp17, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp92 = NegFusedMulSubSIMD(FDPart3_Rational_1_3, hDD00, FDPart3_Rational_1_3);
        const REAL_SIMD_ARRAY FDPart3tmp94 = MulSIMD(FDPart3tmp93, T4UU11);
        const REAL_SIMD_ARRAY FDPart3tmp95 = MulSIMD(FDPart3tmp93, T4UU22);
        const REAL_SIMD_ARRAY FDPart3tmp96 = MulSIMD(FDPart3tmp93, T4UU12);
        const REAL_SIMD_ARRAY FDPart3tmp97 = MulSIMD(FDPart3tmp93, T4UU33);
        const REAL_SIMD_ARRAY FDPart3tmp98 = MulSIMD(FDPart3tmp93, T4UU13);
        const REAL_SIMD_ARRAY FDPart3tmp99 = MulSIMD(FDPart3tmp12, FDPart3tmp53);
        const REAL_SIMD_ARRAY FDPart3tmp101 = MulSIMD(FDPart3tmp93, T4UU23);
        const REAL_SIMD_ARRAY FDPart3tmp102 = MulSIMD(FDPart3tmp7, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp107 = DivSIMD(FDPart3_Integer_1, FDPart3tmp89);
        const REAL_SIMD_ARRAY FDPart3tmp120 = MulSIMD(FDPart3tmp7, MulSIMD(hDD01, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp121 = MulSIMD(FDPart3tmp28, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp171 = DivSIMD(alpha, FDPart3tmp169);
        const REAL_SIMD_ARRAY FDPart3tmp174 = FusedMulAddSIMD(FDPart3tmp12, hDD_dD010, SubSIMD(FDPart3tmp28, hDD_dD001));
        const REAL_SIMD_ARRAY FDPart3tmp175 = MulSIMD(FDPart3tmp12, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp183 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp172);
        const REAL_SIMD_ARRAY FDPart3tmp185 = MulSIMD(FDPart3tmp172, cf_dD1);
        const REAL_SIMD_ARRAY FDPart3tmp188 = MulSIMD(FDPart3tmp172, cf_dD0);
        const REAL_SIMD_ARRAY FDPart3tmp193 = FusedMulAddSIMD(FDPart3tmp12, hDD11, FDPart3tmp12);
        const REAL_SIMD_ARRAY FDPart3tmp197 = MulSIMD(f0_of_xx0, MulSIMD(f1_of_xx1__D1, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp222 = MulSIMD(FDPart3tmp172, cf_dD2);
        const REAL_SIMD_ARRAY FDPart3tmp229 = MulSIMD(FDPart3tmp12, FDPart3tmp29);
        const REAL_SIMD_ARRAY FDPart3tmp269 = MulSIMD(FDPart3tmp7, aDD11);
        const REAL_SIMD_ARRAY FDPart3tmp272 = DivSIMD(FDPart3_Integer_1, FDPart3tmp29);
        const REAL_SIMD_ARRAY FDPart3tmp290 = FusedMulSubSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(FDPart3tmp7, hDD11),
                                                              MulSIMD(FDPart3_Rational_1_3, FDPart3tmp7));
        const REAL_SIMD_ARRAY FDPart3tmp298 =
            FusedMulSubSIMD(FDPart3tmp29, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(FDPart3tmp7, hDD22)),
                            MulSIMD(FDPart3_Rational_1_3, MulSIMD(FDPart3tmp29, FDPart3tmp7)));
        const REAL_SIMD_ARRAY FDPart3tmp302 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(alpha, trK_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp316 = MulSIMD(FDPart3tmp29, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp323 = MulSIMD(FDPart3tmp4, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp335 = MulSIMD(FDPart3_Integer_2, FDPart3tmp30);
        const REAL_SIMD_ARRAY FDPart3tmp339 = MulSIMD(FDPart3_Integer_4, MulSIMD(FDPart3tmp29, FDPart3tmp61));
        const REAL_SIMD_ARRAY FDPart3tmp441 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(alpha, trK_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp450 = DivSIMD(FDPart3_Integer_2, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp479 = DivSIMD(FDPart3tmp337, MulSIMD(MulSIMD(f1_of_xx1, f1_of_xx1), f1_of_xx1));
        const REAL_SIMD_ARRAY FDPart3tmp6 = MulSIMD(FDPart3tmp5, vetU2);
        const REAL_SIMD_ARRAY FDPart3tmp22 = MulSIMD(FDPart3tmp21, vetU0);
        const REAL_SIMD_ARRAY FDPart3tmp24 = MulSIMD(FDPart3tmp5, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp32 = MulSIMD(FDPart3tmp30, MulSIMD(MulSIMD(FDPart3tmp28, FDPart3tmp29), MulSIMD(hDD02, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp34 = MulSIMD(FDPart3tmp29, MulSIMD(FDPart3tmp30, MulSIMD(hDD12, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp38 = AddSIMD(FDPart3tmp37, FDPart3tmp7);
        const REAL_SIMD_ARRAY FDPart3tmp40 = MulSIMD(FDPart3tmp39, MulSIMD(hDD02, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp43 = MulSIMD(FDPart3tmp39, hDD22);
        const REAL_SIMD_ARRAY FDPart3tmp47 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp7, MulSIMD(hDD01, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp55 = MulSIMD(FDPart3tmp53, MulSIMD(FDPart3tmp7, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp56 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp7, f1_of_xx1));
        const REAL_SIMD_ARRAY FDPart3tmp63 = MulSIMD(FDPart3tmp61, MulSIMD(f1_of_xx1, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp74 = MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp61), MulSIMD(hDD02, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp103 = MulSIMD(FDPart3_Integer_2, FDPart3tmp102);
        const REAL_SIMD_ARRAY FDPart3tmp112 = MulSIMD(FDPart3tmp107, vetU2);
        const REAL_SIMD_ARRAY FDPart3tmp123 = MulSIMD(FDPart3tmp107, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp134 = MulSIMD(FDPart3tmp107, FDPart3tmp33);
        const REAL_SIMD_ARRAY FDPart3tmp142 = MulSIMD(FDPart3tmp102, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp150 = MulSIMD(FDPart3tmp65, T4UU03);
        const REAL_SIMD_ARRAY FDPart3tmp165 = MulSIMD(FDPart3tmp107, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp176 = FusedMulAddSIMD(FDPart3tmp175, hDD_dD020, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp53, hDD_dD002));
        const REAL_SIMD_ARRAY FDPart3tmp190 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp188, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp194 = FusedMulAddSIMD(FDPart3tmp7, hDD_dD110, FDPart3tmp193);
        const REAL_SIMD_ARRAY FDPart3tmp199 = FusedMulAddSIMD(FDPart3tmp109, hDD_dD021, FDPart3tmp197);
        const REAL_SIMD_ARRAY FDPart3tmp200 = MulSIMD(FDPart3tmp175, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp212 = MulSIMD(FDPart3tmp7, MulSIMD(f1_of_xx1__D1, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp214 =
            SubSIMD(NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, hDD11),
                                       FusedMulSubSIMD(FDPart3tmp12, hDD_dD011, MulSIMD(FDPart3_Integer_2, f0_of_xx0))),
                    MulSIMD(FDPart3tmp7, hDD_dD110));
        const REAL_SIMD_ARRAY FDPart3tmp230 = FusedMulAddSIMD(FDPart3tmp229, hDD22, FDPart3tmp229);
        const REAL_SIMD_ARRAY FDPart3tmp239 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp222, alpha_dD1));
        const REAL_SIMD_ARRAY FDPart3tmp244 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp242, FDPart3tmp7));
        const REAL_SIMD_ARRAY FDPart3tmp271 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp2, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp275 = MulSIMD(FDPart3tmp102, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp281 = MulSIMD(FDPart3tmp39, aDD22);
        const REAL_SIMD_ARRAY FDPart3tmp296 = MulSIMD(FDPart3tmp4, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp299 = MulSIMD(FDPart3tmp272, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp309 = FusedMulAddSIMD(FDPart3tmp5, vetU_dDD222, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, vetU_dD02)));
        const REAL_SIMD_ARRAY FDPart3tmp312 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp2, vetU_dD01));
        const REAL_SIMD_ARRAY FDPart3tmp318 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp2, vetU_dD02));
        const REAL_SIMD_ARRAY FDPart3tmp327 = MulSIMD(FDPart3tmp107, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp330 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp2, vetU1));
        const REAL_SIMD_ARRAY FDPart3tmp338 = MulSIMD(FDPart3_Rational_1_2, DivSIMD(FDPart3tmp272, FDPart3tmp30));
        const REAL_SIMD_ARRAY FDPart3tmp341 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp2, FDPart3tmp272));
        const REAL_SIMD_ARRAY FDPart3tmp358 = MulSIMD(FDPart3tmp3, betU0);
        const REAL_SIMD_ARRAY FDPart3tmp360 = MulSIMD(FDPart3tmp2, MulSIMD(betU1, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp363 = MulSIMD(FDPart3tmp8, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp370 = MulSIMD(FDPart3tmp14, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp373 = MulSIMD(FDPart3tmp39, hDD_dD222);
        const REAL_SIMD_ARRAY FDPart3tmp440 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp2, eta));
        const REAL_SIMD_ARRAY FDPart3tmp442 = FusedMulSubSIMD(FDPart3tmp2, vetU_dD00, MulSIMD(FDPart3tmp8, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp474 = FusedMulSubSIMD(FDPart3tmp4, f1_of_xx1__DD11, MulSIMD(FDPart3tmp272, FDPart3tmp337));
        const REAL_SIMD_ARRAY FDPart3tmp504 = FusedMulAddSIMD(FDPart3tmp188, FDPart3tmp503, alpha_dD0);
        const REAL_SIMD_ARRAY FDPart3tmp508 = FusedMulAddSIMD(FDPart3tmp222, FDPart3tmp503, alpha_dD2);
        const REAL_SIMD_ARRAY FDPart3tmp509 = FusedMulAddSIMD(FDPart3tmp185, FDPart3tmp503, alpha_dD1);
        const REAL_SIMD_ARRAY FDPart3tmp528 = MulSIMD(FDPart3tmp8, MulSIMD(vetU1, vetU1));
        const REAL_SIMD_ARRAY FDPart3tmp10 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp8, vetU1));
        const REAL_SIMD_ARRAY FDPart3tmp15 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp14, FDPart3tmp8));
        const REAL_SIMD_ARRAY FDPart3tmp44 = AddSIMD(FDPart3tmp39, FDPart3tmp43);
        const REAL_SIMD_ARRAY FDPart3tmp57 = MulSIMD(FDPart3tmp33, MulSIMD(FDPart3tmp56, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp64 = MulSIMD(FDPart3tmp63, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp66 = MulSIMD(FDPart3tmp38, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp111 = MulSIMD(MulSIMD(FDPart3tmp107, FDPart3tmp109), MulSIMD(hDD12, vetU1));
        const REAL_SIMD_ARRAY FDPart3tmp115 = MulSIMD(FDPart3tmp107, MulSIMD(FDPart3tmp12, FDPart3tmp53));
        const REAL_SIMD_ARRAY FDPart3tmp124 = MulSIMD(FDPart3tmp112, MulSIMD(f0_of_xx0, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp125 = MulSIMD(FDPart3tmp107, FDPart3tmp38);
        const REAL_SIMD_ARRAY FDPart3tmp138 = MulSIMD(FDPart3tmp134, T4UU01);
        const REAL_SIMD_ARRAY FDPart3tmp146 = MulSIMD(FDPart3tmp123, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp148 = MulSIMD(FDPart3tmp142, T4UU03);
        const REAL_SIMD_ARRAY FDPart3tmp187 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp185, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp201 = FusedMulAddSIMD(FDPart3tmp102, hDD_dD120, FDPart3tmp200);
        const REAL_SIMD_ARRAY FDPart3tmp213 =
            FusedMulAddSIMD(FDPart3tmp103, hDD_dD121, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp212, MulSIMD(FDPart3tmp7, hDD_dD112)));
        const REAL_SIMD_ARRAY FDPart3tmp224 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp222, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp226 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp188, alpha_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp231 = FusedMulAddSIMD(FDPart3tmp39, hDD_dD220, FDPart3tmp230);
        const REAL_SIMD_ARRAY FDPart3tmp240 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp185, alpha_dD2));
        const REAL_SIMD_ARRAY FDPart3tmp245 = FusedMulAddSIMD(FDPart3tmp244, hDD22, FDPart3tmp244);
        const REAL_SIMD_ARRAY FDPart3tmp247 = AddSIMD(
            FDPart3tmp199, FusedMulAddSIMD(FDPart3tmp56, hDD_dD120,
                                           NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, f0_of_xx0), MulSIMD(f1_of_xx1, hDD12), FDPart3tmp195)));
        const REAL_SIMD_ARRAY FDPart3tmp256 =
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp29, f0_of_xx0),
                               FusedMulAddSIMD(FDPart3tmp29, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, hDD22)),
                                               FusedMulSubSIMD(FDPart3tmp175, hDD_dD022, MulSIMD(FDPart3tmp39, hDD_dD220))));
        const REAL_SIMD_ARRAY FDPart3tmp257 = FusedMulAddSIMD(
            FDPart3tmp7, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, f1_of_xx1__D1)),
            FusedMulAddSIMD(MulSIMD(FDPart3tmp7, f1_of_xx1), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, hDD22)),
                            FusedMulSubSIMD(FDPart3tmp103, hDD_dD122, MulSIMD(FDPart3tmp39, hDD_dD221))));
        const REAL_SIMD_ARRAY FDPart3tmp273 = MulSIMD(FDPart3tmp271, MulSIMD(FDPart3tmp272, f1_of_xx1__D1));
        const REAL_SIMD_ARRAY FDPart3tmp310 =
            FusedMulAddSIMD(FDPart3tmp2, vetU_dDD112, FusedMulAddSIMD(FDPart3tmp24, vetU_dD12, AddSIMD(FDPart3tmp309, vetU_dDD002)));
        const REAL_SIMD_ARRAY FDPart3tmp342 = MulSIMD(FDPart3tmp341, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp346 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp8, vetU_dD11));
        const REAL_SIMD_ARRAY FDPart3tmp348 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp4, FDPart3tmp8));
        const REAL_SIMD_ARRAY FDPart3tmp364 = MulSIMD(FDPart3tmp363, betU1);
        const REAL_SIMD_ARRAY FDPart3tmp366 = MulSIMD(FDPart3tmp6, betU0);
        const REAL_SIMD_ARRAY FDPart3tmp368 = MulSIMD(FDPart3tmp5, MulSIMD(betU2, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp371 = MulSIMD(FDPart3tmp370, betU1);
        const REAL_SIMD_ARRAY FDPart3tmp372 = MulSIMD(FDPart3tmp296, MulSIMD(betU2, vetU1));
        const REAL_SIMD_ARRAY FDPart3tmp375 = MulSIMD(FDPart3tmp299, vetU2);
        const REAL_SIMD_ARRAY FDPart3tmp393 = MulSIMD(FDPart3tmp103, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp527 = MulSIMD(FDPart3tmp22, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp529 = MulSIMD(FDPart3tmp14, FDPart3tmp22);
        const REAL_SIMD_ARRAY FDPart3tmp530 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp14, FDPart3tmp363));
        const REAL_SIMD_ARRAY FDPart3tmp531 = MulSIMD(FDPart3tmp299, MulSIMD(vetU2, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp11 = FusedMulAddSIMD(FDPart3tmp2, vetU_dD10, FDPart3tmp10);
        const REAL_SIMD_ARRAY FDPart3tmp16 = FusedMulAddSIMD(FDPart3tmp5, vetU_dD20, FDPart3tmp15);
        const REAL_SIMD_ARRAY FDPart3tmp26 =
            AddSIMD(vetU_dD00, FusedMulAddSIMD(FDPart3tmp24, vetU1, FusedMulAddSIMD(FDPart3tmp5, vetU_dD22, AddSIMD(FDPart3tmp20, FDPart3tmp22))));
        const REAL_SIMD_ARRAY FDPart3tmp36 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp33, FDPart3tmp34));
        const REAL_SIMD_ARRAY FDPart3tmp42 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp38, FDPart3tmp40));
        const REAL_SIMD_ARRAY FDPart3tmp49 = MulSIMD(FDPart3tmp38, FDPart3tmp44);
        const REAL_SIMD_ARRAY FDPart3tmp58 = AddSIMD(FDPart3tmp55, FDPart3tmp57);
        const REAL_SIMD_ARRAY FDPart3tmp71 = FusedMulAddSIMD(FDPart3tmp33, FDPart3tmp38, FDPart3tmp47);
        const REAL_SIMD_ARRAY FDPart3tmp76 = MulSIMD(FDPart3tmp44, FDPart3tmp75);
        const REAL_SIMD_ARRAY FDPart3tmp105 = MulSIMD(FDPart3tmp101, MulSIMD(FDPart3tmp103, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp126 =
            FusedMulAddSIMD(FDPart3tmp125, FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp123, MulSIMD(f0_of_xx0, vetU0), FDPart3tmp124));
        const REAL_SIMD_ARRAY FDPart3tmp129 = MulSIMD(FDPart3tmp125, T4UU02);
        const REAL_SIMD_ARRAY FDPart3tmp135 =
            FusedMulAddSIMD(FDPart3tmp123, vetU1, FusedMulAddSIMD(FDPart3tmp134, vetU0, MulSIMD(FDPart3tmp112, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp145 = MulSIMD(FDPart3tmp33, MulSIMD(FDPart3tmp38, FDPart3tmp93));
        const REAL_SIMD_ARRAY FDPart3tmp149 = MulSIMD(FDPart3tmp146, T4UU02);
        const REAL_SIMD_ARRAY FDPart3tmp202 = AddSIMD(FDPart3tmp201, SubSIMD(FDPart3tmp199, FDPart3tmp195));
        const REAL_SIMD_ARRAY FDPart3tmp227 =
            AddSIMD(FDPart3tmp201, SubSIMD(SubSIMD(FDPart3tmp195, FDPart3tmp197), MulSIMD(FDPart3tmp109, hDD_dD021)));
        const REAL_SIMD_ARRAY FDPart3tmp246 = FusedMulAddSIMD(FDPart3tmp39, hDD_dD221, FDPart3tmp245);
        const REAL_SIMD_ARRAY FDPart3tmp274 = FusedMulAddSIMD(FDPart3tmp5, vetU_dD21, FDPart3tmp273);
        const REAL_SIMD_ARRAY FDPart3tmp311 = MulSIMD(FDPart3_Rational_1_4, FDPart3tmp310);
        const REAL_SIMD_ARRAY FDPart3tmp343 = FusedMulAddSIMD(FDPart3tmp342, vetU_dD22, MulSIMD(FDPart3tmp5, vetU_dDD212));
        const REAL_SIMD_ARRAY FDPart3tmp347 = FusedMulAddSIMD(FDPart3tmp2, vetU_dDD101, FDPart3tmp346);
        const REAL_SIMD_ARRAY FDPart3tmp349 = FusedMulAddSIMD(FDPart3tmp348, vetU_dD22, MulSIMD(FDPart3tmp5, vetU_dDD202));
        const REAL_SIMD_ARRAY FDPart3tmp376 = MulSIMD(FDPart3tmp375, betU2);
        const REAL_SIMD_ARRAY FDPart3tmp499 = MulSIMD(FDPart3_Rational_1_3, FDPart3tmp310);
        const REAL_SIMD_ARRAY FDPart3tmp68 = SubSIMD(FDPart3tmp64, FDPart3tmp66);
        const REAL_SIMD_ARRAY FDPart3tmp80 = FusedMulSubSIMD(FDPart3tmp33, FDPart3tmp44, FDPart3tmp40);
        const REAL_SIMD_ARRAY FDPart3tmp86 = SubSIMD(FDPart3tmp49, FDPart3tmp34);
        const REAL_SIMD_ARRAY FDPart3tmp114 =
            FusedMulAddSIMD(FDPart3tmp112, MulSIMD(FDPart3tmp44, FDPart3tmp5),
                            FusedMulAddSIMD(MulSIMD(FDPart3tmp107, FDPart3tmp53), MulSIMD(f0_of_xx0, vetU0), FDPart3tmp111));
        const REAL_SIMD_ARRAY FDPart3tmp117 = MulSIMD(FDPart3tmp107, MulSIMD(FDPart3tmp44, T4UU03));
        const REAL_SIMD_ARRAY FDPart3tmp127 = MulSIMD(FDPart3tmp107, FDPart3tmp126);
        const REAL_SIMD_ARRAY FDPart3tmp164 = MulSIMD(FDPart3tmp33, MulSIMD(FDPart3tmp44, FDPart3tmp93));
        const REAL_SIMD_ARRAY FDPart3tmp313 = MulSIMD(FDPart3tmp11, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp325 =
            SubSIMD(SubSIMD(SubSIMD(vetU_dDD012, vetU_dD12), MulSIMD(FDPart3tmp323, vetU_dD02)), MulSIMD(FDPart3tmp274, FDPart3tmp316));
        const REAL_SIMD_ARRAY FDPart3tmp329 = FusedMulAddSIMD(
            FDPart3tmp150, FDPart3tmp327,
            FusedMulAddSIMD(FDPart3tmp328,
                            SubSIMD(NegFusedMulSubSIMD(FDPart3tmp123, vetU1, MulSIMD(FDPart3tmp112, hDD02)), MulSIMD(FDPart3tmp134, vetU0)),
                            FusedMulAddSIMD(FDPart3tmp138, alpha, MulSIMD(FDPart3tmp149, alpha))));
        const REAL_SIMD_ARRAY FDPart3tmp331 = FusedMulAddSIMD(
            FDPart3tmp328, SubSIMD(FusedMulSubSIMD(FDPart3tmp125, FDPart3tmp330, FDPart3tmp124), MulSIMD(FDPart3tmp146, vetU0)),
            FusedMulAddSIMD(FDPart3tmp146, MulSIMD(T4UU01, alpha), FusedMulAddSIMD(FDPart3tmp129, alpha, MulSIMD(FDPart3tmp148, FDPart3tmp327))));
        const REAL_SIMD_ARRAY FDPart3tmp344 = FusedMulAddSIMD(
            FDPart3tmp338,
            FusedMulAddSIMD(FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp335, FDPart3tmp337, MulSIMD(FDPart3tmp335, MulSIMD(f1_of_xx1, f1_of_xx1__DD11))),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp242), MulSIMD(FDPart3tmp61, vetU_dD11))),
            FusedMulAddSIMD(
                FDPart3tmp338,
                FusedMulAddSIMD(FDPart3tmp339, vetU_dD01, MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3tmp242), MulSIMD(FDPart3tmp61, vetU0))),
                AddSIMD(
                    FusedMulAddSIMD(FDPart3tmp2, vetU_dDD111, vetU_dDD001),
                    FusedMulAddSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp4),
                                    MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, vetU0)),
                                    NegFusedMulAddSIMD(FDPart3tmp272, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp2), MulSIMD(FDPart3tmp337, vetU1)),
                                                       FDPart3tmp343)))));
        const REAL_SIMD_ARRAY FDPart3tmp350 = FusedMulAddSIMD(
            FDPart3tmp338, FusedMulAddSIMD(FDPart3tmp339, vetU_dD00, MulSIMD(FDPart3_Integer_12, MulSIMD(FDPart3tmp39, vetU0))),
            FusedMulAddSIMD(FDPart3tmp338,
                            FusedMulAddSIMD(FDPart3tmp11, MulSIMD(FDPart3tmp242, FDPart3tmp335),
                                            MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3tmp242), MulSIMD(FDPart3tmp7, vetU1))),
                            AddSIMD(AddSIMD(FDPart3tmp349, vetU_dDD000),
                                    FusedMulAddSIMD(MulSIMD(FDPart3tmp4, FDPart3tmp8),
                                                    MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, vetU1)),
                                                    NegFusedMulAddSIMD(FDPart3_Integer_8, MulSIMD(FDPart3tmp8, vetU0), FDPart3tmp347)))));
        const REAL_SIMD_ARRAY FDPart3tmp386 = MulSIMD(FDPart3tmp58, FDPart3tmp58);
        const REAL_SIMD_ARRAY FDPart3tmp444 =
            SubSIMD(FusedMulSubSIMD(FDPart3tmp2, vetU_dDD102, MulSIMD(FDPart3tmp16, FDPart3tmp242)), MulSIMD(FDPart3tmp8, vetU_dD12));
        const REAL_SIMD_ARRAY FDPart3tmp471 = AddSIMD(FDPart3tmp349, FusedMulAddSIMD(FDPart3tmp11, FDPart3tmp323, FDPart3tmp442));
        const REAL_SIMD_ARRAY FDPart3tmp472 = FusedMulAddSIMD(FDPart3tmp2, vetU_dD01, FusedMulAddSIMD(FDPart3tmp20, FDPart3tmp323, FDPart3tmp343));
        const REAL_SIMD_ARRAY FDPart3tmp476 =
            FusedMulAddSIMD(FDPart3tmp348, vetU_dD21,
                            FusedMulAddSIMD(FDPart3tmp375, f1_of_xx1__D1,
                                            FusedMulAddSIMD(FDPart3tmp5, vetU_dDD201,
                                                            FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp323, MulSIMD(FDPart3tmp342, vetU_dD20)))));
        const REAL_SIMD_ARRAY FDPart3tmp51 =
            AddSIMD(FDPart3tmp42,
                    FusedMulAddSIMD(FDPart3tmp33, FDPart3tmp49, FusedMulAddSIMD(FDPart3tmp44, FDPart3tmp47, AddSIMD(FDPart3tmp32, FDPart3tmp36))));
        const REAL_SIMD_ARRAY FDPart3tmp78 = SubSIMD(FDPart3tmp74, FDPart3tmp76);
        const REAL_SIMD_ARRAY FDPart3tmp118 = FusedMulAddSIMD(
            FDPart3_Integer_2, MulSIMD(FDPart3tmp74, FDPart3tmp96),
            FusedMulAddSIMD(
                FDPart3tmp114, MulSIMD(FDPart3tmp115, T4UU01),
                FusedMulAddSIMD(
                    T4UU00, MulSIMD(FDPart3tmp114, FDPart3tmp114),
                    FusedMulAddSIMD(
                        FDPart3_Integer_2, MulSIMD(FDPart3tmp114, FDPart3tmp117),
                        FusedMulAddSIMD(
                            FDPart3tmp40, FDPart3tmp94,
                            FusedMulAddSIMD(
                                FDPart3tmp97, MulSIMD(FDPart3tmp44, FDPart3tmp44),
                                FusedMulAddSIMD(
                                    FDPart3tmp44, MulSIMD(FDPart3tmp98, FDPart3tmp99),
                                    FusedMulAddSIMD(FDPart3tmp114, MulSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp107), MulSIMD(T4UU02, hDD12)),
                                                    FusedMulAddSIMD(FDPart3tmp105, FDPart3tmp44, MulSIMD(FDPart3tmp34, FDPart3tmp95))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp137 = MulSIMD(FDPart3tmp107, MulSIMD(FDPart3tmp135, T4UU02));
        const REAL_SIMD_ARRAY FDPart3tmp152 = FusedMulAddSIMD(
            FDPart3tmp135, MulSIMD(FDPart3tmp146, T4UU01),
            FusedMulAddSIMD(
                FDPart3tmp142, MulSIMD(FDPart3tmp33, FDPart3tmp98),
                FusedMulAddSIMD(
                    FDPart3tmp107, MulSIMD(FDPart3tmp135, FDPart3tmp148),
                    FusedMulAddSIMD(
                        FDPart3tmp126, MulSIMD(FDPart3tmp135, T4UU00),
                        FusedMulAddSIMD(
                            FDPart3tmp55, FDPart3tmp98,
                            FusedMulAddSIMD(
                                FDPart3tmp74, FDPart3tmp97,
                                FusedMulAddSIMD(
                                    FDPart3tmp129, FDPart3tmp135,
                                    FusedMulAddSIMD(
                                        FDPart3tmp145, T4UU12,
                                        FusedMulAddSIMD(
                                            FDPart3tmp126, FDPart3tmp149,
                                            FusedMulAddSIMD(
                                                FDPart3tmp127, FDPart3tmp150,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp120, FDPart3tmp96,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp126, FDPart3tmp138,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp33, MulSIMD(FDPart3tmp75, FDPart3tmp94),
                                                            FusedMulAddSIMD(FDPart3tmp38, MulSIMD(FDPart3tmp75, FDPart3tmp95),
                                                                            FusedMulAddSIMD(FDPart3tmp101, FDPart3tmp64,
                                                                                            MulSIMD(FDPart3tmp101, FDPart3tmp66))))))))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp314 = AddSIMD(vetU_dDD001, SubSIMD(FDPart3tmp312, FDPart3tmp313));
        const REAL_SIMD_ARRAY FDPart3tmp320 = AddSIMD(vetU_dDD002, NegFusedMulAddSIMD(FDPart3tmp16, FDPart3tmp316, FDPart3tmp318));
        const REAL_SIMD_ARRAY FDPart3tmp332 = FusedMulAddSIMD(
            FDPart3tmp142, MulSIMD(FDPart3tmp327, T4UU02),
            FusedMulAddSIMD(FDPart3tmp165, MulSIMD(T4UU01, alpha),
                            FusedMulAddSIMD(FDPart3tmp117, alpha,
                                            MulSIMD(FDPart3tmp328, FusedMulAddSIMD(FDPart3tmp2,
                                                                                   MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp112),
                                                                                           MulSIMD(FDPart3tmp4, FDPart3tmp44)),
                                                                                   NegFusedMulSubSIMD(FDPart3tmp165, vetU0, FDPart3tmp111))))));
        const REAL_SIMD_ARRAY FDPart3tmp345 = MulSIMD(FDPart3_Rational_1_4, FDPart3tmp344);
        const REAL_SIMD_ARRAY FDPart3tmp351 = MulSIMD(FDPart3_Rational_1_4, FDPart3tmp350);
        const REAL_SIMD_ARRAY FDPart3tmp377 = MulSIMD(FDPart3tmp68, FDPart3tmp68);
        const REAL_SIMD_ARRAY FDPart3tmp391 = MulSIMD(FDPart3tmp58, FDPart3tmp68);
        const REAL_SIMD_ARRAY FDPart3tmp448 = FusedMulAddSIMD(
            FDPart3tmp2, vetU_dDD112,
            SubSIMD(FusedMulSubSIMD(FDPart3tmp2, vetU_dD02, MulSIMD(FDPart3tmp24, vetU_dD12)), MulSIMD(FDPart3tmp242, FDPart3tmp274)));
        const REAL_SIMD_ARRAY FDPart3tmp501 = MulSIMD(FDPart3_Rational_1_3, FDPart3tmp344);
        const REAL_SIMD_ARRAY FDPart3tmp502 = MulSIMD(FDPart3_Rational_1_3, FDPart3tmp350);
        const REAL_SIMD_ARRAY FDPart3tmp524 = DivSIMD(
            FDPart3_Integer_1,
            FusedMulAddSIMD(FDPart3tmp42, FDPart3tmp523,
                            FusedMulAddSIMD(FDPart3tmp33, MulSIMD(FDPart3tmp49, FDPart3tmp523),
                                            FusedMulAddSIMD(FDPart3tmp44, MulSIMD(FDPart3tmp47, FDPart3tmp523),
                                                            FusedMulAddSIMD(FDPart3tmp32, FDPart3tmp523, MulSIMD(FDPart3tmp36, FDPart3tmp523))))));
        const REAL_SIMD_ARRAY FDPart3tmp52 = DivSIMD(FDPart3_Integer_1, FDPart3tmp51);
        const REAL_SIMD_ARRAY FDPart3tmp130 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(FDPart3tmp127, T4UU01),
            FusedMulAddSIMD(
                FDPart3tmp121, MulSIMD(FDPart3tmp38, FDPart3tmp96),
                FusedMulAddSIMD(
                    T4UU00, MulSIMD(FDPart3tmp126, FDPart3tmp126),
                    FusedMulAddSIMD(
                        FDPart3_Integer_2, MulSIMD(FDPart3tmp126, FDPart3tmp129),
                        FusedMulAddSIMD(FDPart3tmp34, FDPart3tmp97,
                                        FusedMulAddSIMD(FDPart3tmp95, MulSIMD(FDPart3tmp38, FDPart3tmp38),
                                                        FusedMulAddSIMD(FDPart3tmp28, MulSIMD(FDPart3tmp63, FDPart3tmp98),
                                                                        FusedMulAddSIMD(MulSIMD(FDPart3tmp103, FDPart3tmp127), MulSIMD(T4UU03, hDD12),
                                                                                        FusedMulAddSIMD(FDPart3tmp105, FDPart3tmp38,
                                                                                                        MulSIMD(FDPart3tmp120, FDPart3tmp94))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp139 = FusedMulAddSIMD(
            FDPart3tmp115, MulSIMD(FDPart3tmp135, T4UU03),
            FusedMulAddSIMD(
                FDPart3tmp121, MulSIMD(FDPart3tmp33, FDPart3tmp96),
                FusedMulAddSIMD(
                    T4UU00, MulSIMD(FDPart3tmp135, FDPart3tmp135),
                    FusedMulAddSIMD(
                        FDPart3_Integer_2, MulSIMD(FDPart3tmp135, FDPart3tmp138),
                        FusedMulAddSIMD(
                            FDPart3tmp40, FDPart3tmp97,
                            FusedMulAddSIMD(FDPart3tmp94, MulSIMD(FDPart3tmp33, FDPart3tmp33),
                                            FusedMulAddSIMD(FDPart3tmp33, MulSIMD(FDPart3tmp98, FDPart3tmp99),
                                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp101, FDPart3tmp28), MulSIMD(FDPart3tmp53, FDPart3tmp7),
                                                                            FusedMulAddSIMD(FDPart3tmp120, FDPart3tmp95,
                                                                                            MulSIMD(FDPart3tmp121, FDPart3tmp137))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp161 = FusedMulAddSIMD(
            FDPart3tmp127, MulSIMD(FDPart3tmp65, T4UU01),
            FusedMulAddSIMD(
                FDPart3tmp142, MulSIMD(FDPart3tmp38, FDPart3tmp95),
                FusedMulAddSIMD(
                    FDPart3tmp127, MulSIMD(FDPart3tmp142, T4UU02),
                    FusedMulAddSIMD(
                        FDPart3tmp127, MulSIMD(FDPart3tmp44, T4UU03),
                        FusedMulAddSIMD(
                            FDPart3tmp114, MulSIMD(FDPart3tmp126, T4UU00),
                            FusedMulAddSIMD(
                                FDPart3tmp114, MulSIMD(FDPart3tmp146, T4UU01),
                                FusedMulAddSIMD(
                                    FDPart3tmp76, FDPart3tmp98,
                                    FusedMulAddSIMD(
                                        FDPart3tmp107, MulSIMD(FDPart3tmp114, FDPart3tmp148),
                                        FusedMulAddSIMD(
                                            FDPart3tmp66, FDPart3tmp96,
                                            FusedMulAddSIMD(
                                                FDPart3tmp74, FDPart3tmp98,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp55, FDPart3tmp94,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp64, FDPart3tmp96,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp142, MulSIMD(FDPart3tmp44, FDPart3tmp97),
                                                            FusedMulAddSIMD(FDPart3tmp49, MulSIMD(FDPart3tmp93, T4UU23),
                                                                            FusedMulAddSIMD(FDPart3tmp101, FDPart3tmp34,
                                                                                            MulSIMD(FDPart3tmp114, FDPart3tmp129))))))))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp166 = FusedMulAddSIMD(
            FDPart3tmp135, MulSIMD(FDPart3tmp165, T4UU01),
            FusedMulAddSIMD(
                FDPart3tmp142, MulSIMD(FDPart3tmp33, FDPart3tmp96),
                FusedMulAddSIMD(
                    FDPart3tmp107, MulSIMD(FDPart3tmp114, FDPart3tmp150),
                    FusedMulAddSIMD(
                        FDPart3tmp114, MulSIMD(FDPart3tmp135, T4UU00),
                        FusedMulAddSIMD(
                            FDPart3tmp55, FDPart3tmp96,
                            FusedMulAddSIMD(
                                FDPart3tmp64, FDPart3tmp95,
                                FusedMulAddSIMD(
                                    FDPart3tmp164, T4UU13,
                                    FusedMulAddSIMD(
                                        FDPart3tmp40, FDPart3tmp98,
                                        FusedMulAddSIMD(
                                            FDPart3tmp117, FDPart3tmp135,
                                            FusedMulAddSIMD(
                                                FDPart3tmp137, FDPart3tmp142,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp114, FDPart3tmp138,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp114, FDPart3tmp149,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp33, MulSIMD(FDPart3tmp65, FDPart3tmp94),
                                                            FusedMulAddSIMD(FDPart3tmp44, MulSIMD(FDPart3tmp65, FDPart3tmp97),
                                                                            FusedMulAddSIMD(FDPart3tmp101, FDPart3tmp74,
                                                                                            MulSIMD(FDPart3tmp101, FDPart3tmp76))))))))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp378 = DivSIMD(FDPart3_Integer_1, MulSIMD(FDPart3tmp51, FDPart3tmp51));
        const REAL_SIMD_ARRAY FDPart3tmp379 = MulSIMD(FDPart3tmp78, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp403 = MulSIMD(FDPart3tmp58, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp410 = MulSIMD(FDPart3tmp68, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp525 = MulSIMD(FDPart3_Integer_2, FDPart3tmp524);
        const REAL_SIMD_ARRAY FDPart3tmp59 = MulSIMD(FDPart3tmp52, FDPart3tmp58);
        const REAL_SIMD_ARRAY FDPart3tmp72 = MulSIMD(FDPart3tmp52, FDPart3tmp71);
        const REAL_SIMD_ARRAY FDPart3tmp81 = MulSIMD(FDPart3tmp52, FDPart3tmp80);
        const REAL_SIMD_ARRAY FDPart3tmp82 = MulSIMD(FDPart3tmp52, FDPart3tmp68);
        const REAL_SIMD_ARRAY FDPart3tmp84 = MulSIMD(FDPart3tmp52, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp87 = MulSIMD(FDPart3tmp52, FDPart3tmp86);
        const REAL_SIMD_ARRAY FDPart3tmp177 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp68, hDD_dD000)),
            FusedMulSubSIMD(FDPart3tmp176, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp71)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp174), MulSIMD(FDPart3tmp52, FDPart3tmp58))));
        const REAL_SIMD_ARRAY FDPart3tmp179 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp78, hDD_dD000)),
            FusedMulSubSIMD(FDPart3tmp176, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp58)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp174), MulSIMD(FDPart3tmp52, FDPart3tmp80))));
        const REAL_SIMD_ARRAY FDPart3tmp181 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp86, hDD_dD000)),
            FusedMulSubSIMD(FDPart3tmp176, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp68)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp174), MulSIMD(FDPart3tmp52, FDPart3tmp78))));
        const REAL_SIMD_ARRAY FDPart3tmp203 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp68, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp202, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp71)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp194), MulSIMD(FDPart3tmp52, FDPart3tmp58))));
        const REAL_SIMD_ARRAY FDPart3tmp205 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp78, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp202, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp58)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp194), MulSIMD(FDPart3tmp52, FDPart3tmp80))));
        const REAL_SIMD_ARRAY FDPart3tmp207 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp86, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp202, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp68)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp194), MulSIMD(FDPart3tmp52, FDPart3tmp78))));
        const REAL_SIMD_ARRAY FDPart3tmp215 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp7, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp214, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp68)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp213), MulSIMD(FDPart3tmp52, FDPart3tmp71))));
        const REAL_SIMD_ARRAY FDPart3tmp217 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp7), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp80, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp214, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp78)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp213), MulSIMD(FDPart3tmp52, FDPart3tmp58))));
        const REAL_SIMD_ARRAY FDPart3tmp219 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp7), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp78, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp214, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp86)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp213), MulSIMD(FDPart3tmp52, FDPart3tmp68))));
        const REAL_SIMD_ARRAY FDPart3tmp232 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp68, hDD_dD002)),
            FusedMulSubSIMD(FDPart3tmp231, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp71)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp227), MulSIMD(FDPart3tmp52, FDPart3tmp58))));
        const REAL_SIMD_ARRAY FDPart3tmp234 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp78, hDD_dD002)),
            FusedMulSubSIMD(FDPart3tmp231, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp58)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp227), MulSIMD(FDPart3tmp52, FDPart3tmp80))));
        const REAL_SIMD_ARRAY FDPart3tmp236 = FusedMulAddSIMD(
            FDPart3tmp52, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp86, hDD_dD002)),
            FusedMulSubSIMD(FDPart3tmp231, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp68)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp227), MulSIMD(FDPart3tmp52, FDPart3tmp78))));
        const REAL_SIMD_ARRAY FDPart3tmp248 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp7, hDD_dD112)),
            FusedMulSubSIMD(FDPart3tmp247, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp68)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp246), MulSIMD(FDPart3tmp52, FDPart3tmp71))));
        const REAL_SIMD_ARRAY FDPart3tmp250 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp7), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp80, hDD_dD112)),
            FusedMulSubSIMD(FDPart3tmp247, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp78)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp246), MulSIMD(FDPart3tmp52, FDPart3tmp58))));
        const REAL_SIMD_ARRAY FDPart3tmp252 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp7), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp78, hDD_dD112)),
            FusedMulSubSIMD(FDPart3tmp247, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp86)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp246), MulSIMD(FDPart3tmp52, FDPart3tmp68))));
        const REAL_SIMD_ARRAY FDPart3tmp258 = FusedMulAddSIMD(
            FDPart3tmp7,
            MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp52),
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp71, hDD_dD222))),
            FusedMulSubSIMD(FDPart3tmp257, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp58)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp256), MulSIMD(FDPart3tmp52, FDPart3tmp68))));
        const REAL_SIMD_ARRAY FDPart3tmp260 = FusedMulAddSIMD(
            FDPart3tmp58,
            MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp52),
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp7, hDD_dD222))),
            FusedMulSubSIMD(FDPart3tmp257, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp80)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp256), MulSIMD(FDPart3tmp52, FDPart3tmp78))));
        const REAL_SIMD_ARRAY FDPart3tmp262 = FusedMulAddSIMD(
            FDPart3tmp68,
            MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp52),
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp7, hDD_dD222))),
            FusedMulSubSIMD(FDPart3tmp257, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp52, FDPart3tmp78)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp256), MulSIMD(FDPart3tmp52, FDPart3tmp86))));
        const REAL_SIMD_ARRAY FDPart3tmp381 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp7, FDPart3tmp78),
            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp68), MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
            FusedMulAddSIMD(
                FDPart3tmp378,
                MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp377), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp7, aDD22))),
                FusedMulAddSIMD(
                    FDPart3tmp86,
                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp78),
                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(aDD01, f0_of_xx0))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp86, aDD02),
                        MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp68),
                                MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
                        FusedMulSubSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp379),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp7, aDD11)),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp378), MulSIMD(aDD00, MulSIMD(FDPart3tmp86, FDPart3tmp86))))))));
        const REAL_SIMD_ARRAY FDPart3tmp387 = MulSIMD(FDPart3tmp269, FDPart3tmp378);
        const REAL_SIMD_ARRAY FDPart3tmp388 = MulSIMD(FDPart3tmp378, aDD00);
        const REAL_SIMD_ARRAY FDPart3tmp390 = MulSIMD(FDPart3tmp281, FDPart3tmp378);
        const REAL_SIMD_ARRAY FDPart3tmp394 = MulSIMD(FDPart3tmp378, FDPart3tmp393);
        const REAL_SIMD_ARRAY FDPart3tmp413 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp68, FDPart3tmp7),
            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp71, FDPart3tmp78),
                MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp7), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
                FusedMulAddSIMD(
                    FDPart3tmp78,
                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp68),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD01, f0_of_xx0))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp68, FDPart3tmp7),
                        MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp378),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp71, aDD22))),
                        FusedMulAddSIMD(
                            FDPart3tmp7,
                            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp78, aDD11))),
                            FusedMulAddSIMD(
                                FDPart3tmp86,
                                MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD01, f0_of_xx0))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp86, aDD02),
                                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp71),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
                                    FusedMulSubSIMD(aDD02,
                                                    MulSIMD(MulSIMD(FDPart3tmp377, FDPart3tmp378),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
                                                    MulSIMD(FDPart3tmp68,
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp378), MulSIMD(FDPart3tmp86, aDD00)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp415 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp86, aDD02),
            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp7, FDPart3tmp80),
                MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp68), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp58, FDPart3tmp68),
                    MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp378),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp7, aDD22))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp7, FDPart3tmp78),
                        MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
                        FusedMulAddSIMD(
                            FDPart3tmp78,
                            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp7),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp80, aDD11))),
                            FusedMulAddSIMD(
                                FDPart3tmp86,
                                MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp80),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD01, f0_of_xx0))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp78, aDD02),
                                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp68),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
                                    FusedMulSubSIMD(
                                        MulSIMD(FDPart3tmp378, FDPart3tmp78),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp86, aDD00)),
                                        MulSIMD(FDPart3tmp379, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp378), MulSIMD(aDD01, f0_of_xx0)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp417 = MulSIMD(FDPart3tmp378, FDPart3tmp73);
        const REAL_SIMD_ARRAY FDPart3tmp419 = MulSIMD(FDPart3tmp378, FDPart3tmp60);
        const REAL_SIMD_ARRAY FDPart3tmp422 = MulSIMD(FDPart3tmp275, FDPart3tmp378);
        const REAL_SIMD_ARRAY FDPart3tmp465 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp7, FDPart3tmp80),
            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
            FusedMulAddSIMD(
                FDPart3tmp386,
                MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp378), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp7, aDD22))),
                FusedMulAddSIMD(
                    FDPart3tmp80,
                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp78),
                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(aDD01, f0_of_xx0))),
                    FusedMulAddSIMD(MulSIMD(FDPart3tmp78, aDD02),
                                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58),
                                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
                                    FusedMulSubSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp7),
                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                            MulSIMD(aDD11, MulSIMD(FDPart3tmp80, FDPart3tmp80))),
                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp378), MulSIMD(FDPart3tmp379, aDD00)))))));
        const REAL_SIMD_ARRAY FDPart3tmp468 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp68, aDD02),
            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp71, FDPart3tmp80),
                MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp7), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
                FusedMulAddSIMD(
                    FDPart3tmp80,
                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp68),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD01, f0_of_xx0))),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp58, FDPart3tmp7),
                        MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp378),
                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp71, aDD22))),
                        FusedMulAddSIMD(
                            FDPart3tmp7,
                            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp80, aDD11))),
                            FusedMulAddSIMD(
                                FDPart3tmp78,
                                MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58),
                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD01, f0_of_xx0))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp78, aDD02),
                                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp71),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
                                    FusedMulSubSIMD(FDPart3tmp7,
                                                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp386),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
                                                    MulSIMD(FDPart3tmp68,
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp378), MulSIMD(FDPart3tmp78, aDD00)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp489 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp7, FDPart3tmp71),
            MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(aDD12, f1_of_xx1))),
            FusedMulAddSIMD(
                FDPart3tmp7,
                MulSIMD(MulSIMD(FDPart3tmp29, FDPart3tmp378),
                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD22, MulSIMD(FDPart3tmp71, FDPart3tmp71)))),
                FusedMulAddSIMD(
                    FDPart3tmp68,
                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp58),
                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(aDD01, f0_of_xx0))),
                    FusedMulAddSIMD(MulSIMD(FDPart3tmp71, aDD02),
                                    MulSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp68),
                                            MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, f1_of_xx1))),
                                    FusedMulSubSIMD(MulSIMD(FDPart3tmp378, FDPart3tmp386),
                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp7, aDD11)),
                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp377), MulSIMD(FDPart3tmp378, aDD00)))))));
        const REAL_SIMD_ARRAY FDPart3tmp119 = MulSIMD(FDPart3tmp118, FDPart3tmp72);
        const REAL_SIMD_ARRAY FDPart3tmp131 = MulSIMD(FDPart3tmp130, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp140 = MulSIMD(FDPart3tmp139, FDPart3tmp87);
        const REAL_SIMD_ARRAY FDPart3tmp154 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp52, FDPart3tmp78));
        const REAL_SIMD_ARRAY FDPart3tmp162 = MulSIMD(FDPart3_Integer_2, FDPart3tmp59);
        const REAL_SIMD_ARRAY FDPart3tmp167 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp52, FDPart3tmp68));
        const REAL_SIMD_ARRAY FDPart3tmp276 =
            FusedMulAddSIMD(FDPart3tmp275, FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp60, FDPart3tmp82, MulSIMD(FDPart3tmp269, FDPart3tmp59)));
        const REAL_SIMD_ARRAY FDPart3tmp280 =
            FusedMulAddSIMD(FDPart3tmp275, FDPart3tmp82, FusedMulAddSIMD(FDPart3tmp60, FDPart3tmp87, MulSIMD(FDPart3tmp269, FDPart3tmp84)));
        const REAL_SIMD_ARRAY FDPart3tmp284 =
            FusedMulAddSIMD(FDPart3tmp281, FDPart3tmp59, FusedMulAddSIMD(FDPart3tmp73, FDPart3tmp84, MulSIMD(FDPart3tmp275, FDPart3tmp81)));
        const REAL_SIMD_ARRAY FDPart3tmp285 =
            FusedMulAddSIMD(FDPart3tmp281, FDPart3tmp82, FusedMulAddSIMD(FDPart3tmp73, FDPart3tmp87, MulSIMD(FDPart3tmp275, FDPart3tmp84)));
        const REAL_SIMD_ARRAY FDPart3tmp303 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp84, alpha));
        const REAL_SIMD_ARRAY FDPart3tmp315 = MulSIMD(FDPart3_Rational_3_2, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp321 = MulSIMD(FDPart3_Rational_3_2, FDPart3tmp82);
        const REAL_SIMD_ARRAY FDPart3tmp352 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp82);
        const REAL_SIMD_ARRAY FDPart3tmp353 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp87);
        const REAL_SIMD_ARRAY FDPart3tmp354 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp396 = MulSIMD(FDPart3tmp18, MulSIMD(FDPart3tmp378, FDPart3tmp68));
        const REAL_SIMD_ARRAY FDPart3tmp405 = MulSIMD(FDPart3tmp13, MulSIMD(FDPart3tmp378, FDPart3tmp78));
        const REAL_SIMD_ARRAY FDPart3tmp433 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp87);
        const REAL_SIMD_ARRAY FDPart3tmp434 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp435 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp82);
        const REAL_SIMD_ARRAY FDPart3tmp436 = MulSIMD(FDPart3_Rational_4_3, FDPart3tmp59);
        const REAL_SIMD_ARRAY FDPart3tmp437 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp438 = MulSIMD(FDPart3_Rational_2_3, FDPart3tmp72);
        const REAL_SIMD_ARRAY FDPart3tmp446 =
            MulSIMD(FDPart3tmp81,
                    FusedMulAddSIMD(FDPart3tmp2, vetU_dDD111, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, vetU_dD01), FDPart3tmp313)));
        const REAL_SIMD_ARRAY FDPart3tmp456 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp59);
        const REAL_SIMD_ARRAY FDPart3tmp457 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp473 =
            MulSIMD(FDPart3tmp59, AddSIMD(FDPart3tmp330, FusedMulAddSIMD(FDPart3tmp272, MulSIMD(FDPart3tmp330, FDPart3tmp337), FDPart3tmp472)));
        const REAL_SIMD_ARRAY FDPart3tmp475 = MulSIMD(FDPart3tmp59, FusedMulAddSIMD(FDPart3tmp3, FDPart3tmp474, FDPart3tmp472));
        const REAL_SIMD_ARRAY FDPart3tmp478 =
            MulSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp242, FDPart3tmp274,
                                                  FusedMulAddSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp323, vetU_dD12),
                                                                  FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp316, FDPart3tmp309))));
        const REAL_SIMD_ARRAY FDPart3tmp481 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp72);
        const REAL_SIMD_ARRAY FDPart3tmp518 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp169, FDPart3tmp84));
        const REAL_SIMD_ARRAY FDPart3tmp520 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp169, FDPart3tmp82));
        const REAL_SIMD_ARRAY FDPart3tmp522 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp169, FDPart3tmp59));
        const REAL_SIMD_ARRAY FDPart3tmp155 = MulSIMD(FDPart3tmp152, FDPart3tmp154);
        const REAL_SIMD_ARRAY FDPart3tmp163 = MulSIMD(FDPart3tmp161, FDPart3tmp162);
        const REAL_SIMD_ARRAY FDPart3tmp168 = MulSIMD(FDPart3tmp166, FDPart3tmp167);
        const REAL_SIMD_ARRAY FDPart3tmp184 = NegFusedMulAddSIMD(
            FDPart3_Integer_2,
            MulSIMD(alpha,
                    FusedMulAddSIMD(
                        FDPart3tmp172, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp179, cf_dD1)),
                        FusedMulAddSIMD(FDPart3tmp172, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp181, cf_dD0)),
                                        FusedMulSubSIMD(FDPart3tmp183, FusedMulSubSIMD(FDPart3tmp172, MulSIMD(cf_dD0, cf_dD0), cf_dDD00),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp172), MulSIMD(FDPart3tmp177, cf_dD2)))))),
            SubSIMD(
                SubSIMD(FusedMulAddSIMD(RbarDD00, alpha,
                                        SubSIMD(FusedMulAddSIMD(FDPart3tmp172,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha_dD0, cf_dD0)),
                                                                FusedMulSubSIMD(FDPart3tmp171, MulSIMD(cf_dD0, cf_dD0), alpha_dDD00)),
                                                MulSIMD(FDPart3tmp177, alpha_dD2))),
                        MulSIMD(FDPart3tmp181, alpha_dD0)),
                MulSIMD(FDPart3tmp179, alpha_dD1)));
        const REAL_SIMD_ARRAY FDPart3tmp209 = SubSIMD(
            FusedMulAddSIMD(
                FDPart3tmp171, MulSIMD(cf_dD0, cf_dD1),
                SubSIMD(
                    SubSIMD(AddSIMD(FusedMulAddSIMD(RbarDD01, alpha, FDPart3tmp190),
                                    NegFusedMulAddSIMD(
                                        FDPart3_Integer_2,
                                        MulSIMD(alpha,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp172,
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp205, cf_dD1)),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp172,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp207, cf_dD0)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp183, FusedMulSubSIMD(FDPart3tmp185, cf_dD0, cf_dDD01),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp172), MulSIMD(FDPart3tmp203, cf_dD2)))))),
                                        SubSIMD(FDPart3tmp187, alpha_dDD01))),
                            MulSIMD(FDPart3tmp205, alpha_dD1)),
                    MulSIMD(FDPart3tmp203, alpha_dD2))),
            MulSIMD(FDPart3tmp207, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp221 = NegFusedMulAddSIMD(
            FDPart3_Integer_2,
            MulSIMD(alpha,
                    FusedMulAddSIMD(
                        FDPart3tmp172, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp217, cf_dD1)),
                        FusedMulAddSIMD(FDPart3tmp172, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp219, cf_dD0)),
                                        FusedMulSubSIMD(FDPart3tmp183, FusedMulSubSIMD(FDPart3tmp172, MulSIMD(cf_dD1, cf_dD1), cf_dDD11),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp172), MulSIMD(FDPart3tmp215, cf_dD2)))))),
            SubSIMD(
                SubSIMD(FusedMulAddSIMD(RbarDD11, alpha,
                                        SubSIMD(FusedMulAddSIMD(FDPart3tmp172,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha_dD1, cf_dD1)),
                                                                FusedMulSubSIMD(FDPart3tmp171, MulSIMD(cf_dD1, cf_dD1), alpha_dDD11)),
                                                MulSIMD(FDPart3tmp215, alpha_dD2))),
                        MulSIMD(FDPart3tmp219, alpha_dD0)),
                MulSIMD(FDPart3tmp217, alpha_dD1)));
        const REAL_SIMD_ARRAY FDPart3tmp238 = SubSIMD(
            FusedMulAddSIMD(
                FDPart3tmp171, MulSIMD(cf_dD0, cf_dD2),
                SubSIMD(
                    SubSIMD(AddSIMD(FusedMulAddSIMD(RbarDD02, alpha, FDPart3tmp226),
                                    NegFusedMulAddSIMD(
                                        FDPart3_Integer_2,
                                        MulSIMD(alpha,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp172,
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp234, cf_dD1)),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp172,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp236, cf_dD0)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp183, FusedMulSubSIMD(FDPart3tmp188, cf_dD2, cf_dDD02),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp172), MulSIMD(FDPart3tmp232, cf_dD2)))))),
                                        SubSIMD(FDPart3tmp224, alpha_dDD02))),
                            MulSIMD(FDPart3tmp234, alpha_dD1)),
                    MulSIMD(FDPart3tmp232, alpha_dD2))),
            MulSIMD(FDPart3tmp236, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp254 = SubSIMD(
            FusedMulAddSIMD(
                FDPart3tmp171, MulSIMD(cf_dD1, cf_dD2),
                SubSIMD(
                    SubSIMD(AddSIMD(FusedMulAddSIMD(RbarDD12, alpha, FDPart3tmp240),
                                    NegFusedMulAddSIMD(
                                        FDPart3_Integer_2,
                                        MulSIMD(alpha,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp172,
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp250, cf_dD1)),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp172,
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp252, cf_dD0)),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp183, FusedMulSubSIMD(FDPart3tmp185, cf_dD2, cf_dDD12),
                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp172), MulSIMD(FDPart3tmp248, cf_dD2)))))),
                                        SubSIMD(FDPart3tmp239, alpha_dDD12))),
                            MulSIMD(FDPart3tmp250, alpha_dD1)),
                    MulSIMD(FDPart3tmp248, alpha_dD2))),
            MulSIMD(FDPart3tmp252, alpha_dD0));
        const REAL_SIMD_ARRAY FDPart3tmp264 = NegFusedMulAddSIMD(
            FDPart3_Integer_2,
            MulSIMD(alpha,
                    FusedMulAddSIMD(
                        FDPart3tmp172, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp260, cf_dD1)),
                        FusedMulAddSIMD(FDPart3tmp172, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp262, cf_dD0)),
                                        FusedMulSubSIMD(FDPart3tmp183, FusedMulSubSIMD(FDPart3tmp172, MulSIMD(cf_dD2, cf_dD2), cf_dDD22),
                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp172), MulSIMD(FDPart3tmp258, cf_dD2)))))),
            SubSIMD(
                SubSIMD(FusedMulAddSIMD(RbarDD22, alpha,
                                        SubSIMD(FusedMulAddSIMD(FDPart3tmp172,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(alpha_dD2, cf_dD2)),
                                                                FusedMulSubSIMD(FDPart3tmp171, MulSIMD(cf_dD2, cf_dD2), alpha_dDD22)),
                                                MulSIMD(FDPart3tmp258, alpha_dD2))),
                        MulSIMD(FDPart3tmp262, alpha_dD0)),
                MulSIMD(FDPart3tmp260, alpha_dD1)));
        const REAL_SIMD_ARRAY FDPart3tmp279 =
            FusedMulAddSIMD(FDPart3tmp275, FDPart3tmp59, FusedMulAddSIMD(FDPart3tmp60, FDPart3tmp84, MulSIMD(FDPart3tmp269, FDPart3tmp81)));
        const REAL_SIMD_ARRAY FDPart3tmp283 =
            FusedMulAddSIMD(FDPart3tmp281, FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp73, FDPart3tmp82, MulSIMD(FDPart3tmp275, FDPart3tmp59)));
        const REAL_SIMD_ARRAY FDPart3tmp355 =
            FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp352, FusedMulAddSIMD(FDPart3tmp353, hDD_dD000, MulSIMD(FDPart3tmp174, FDPart3tmp354)));
        const REAL_SIMD_ARRAY FDPart3tmp357 =
            FusedMulAddSIMD(FDPart3tmp202, FDPart3tmp352, FusedMulAddSIMD(FDPart3tmp353, hDD_dD001, MulSIMD(FDPart3tmp194, FDPart3tmp354)));
        const REAL_SIMD_ARRAY FDPart3tmp365 =
            FusedMulAddSIMD(FDPart3tmp231, FDPart3tmp352, FusedMulAddSIMD(FDPart3tmp353, hDD_dD002, MulSIMD(FDPart3tmp227, FDPart3tmp354)));
        const REAL_SIMD_ARRAY FDPart3tmp374 =
            FusedMulAddSIMD(FDPart3tmp257, FDPart3tmp354, FusedMulAddSIMD(FDPart3tmp352, FDPart3tmp373, MulSIMD(FDPart3tmp256, FDPart3tmp353)));
        const REAL_SIMD_ARRAY FDPart3tmp397 = FusedMulAddSIMD(
            FDPart3tmp390, MulSIMD(FDPart3tmp71, FDPart3tmp71),
            FusedMulAddSIMD(FDPart3tmp396, FDPart3tmp71,
                            FusedMulAddSIMD(FDPart3tmp13, MulSIMD(FDPart3tmp378, FDPart3tmp391),
                                            FusedMulAddSIMD(FDPart3tmp394, MulSIMD(FDPart3tmp58, FDPart3tmp71),
                                                            FusedMulAddSIMD(FDPart3tmp377, FDPart3tmp388, MulSIMD(FDPart3tmp386, FDPart3tmp387))))));
        const REAL_SIMD_ARRAY FDPart3tmp406 = FusedMulAddSIMD(
            FDPart3tmp405, FDPart3tmp80,
            FusedMulAddSIMD(FDPart3tmp18, MulSIMD(FDPart3tmp378, FDPart3tmp403),
                            FusedMulAddSIMD(FDPart3tmp281, MulSIMD(FDPart3tmp378, FDPart3tmp386),
                                            FusedMulAddSIMD(FDPart3tmp394, MulSIMD(FDPart3tmp58, FDPart3tmp80),
                                                            FusedMulAddSIMD(FDPart3tmp379, FDPart3tmp388,
                                                                            MulSIMD(FDPart3tmp387, MulSIMD(FDPart3tmp80, FDPart3tmp80)))))));
        const REAL_SIMD_ARRAY FDPart3tmp411 =
            FusedMulAddSIMD(FDPart3tmp396, FDPart3tmp86,
                            FusedMulAddSIMD(FDPart3tmp405, FDPart3tmp86,
                                            FusedMulAddSIMD(FDPart3tmp269, MulSIMD(FDPart3tmp378, FDPart3tmp379),
                                                            FusedMulAddSIMD(FDPart3tmp281, MulSIMD(FDPart3tmp377, FDPart3tmp378),
                                                                            FusedMulAddSIMD(FDPart3tmp388, MulSIMD(FDPart3tmp86, FDPart3tmp86),
                                                                                            MulSIMD(FDPart3tmp394, FDPart3tmp410))))));
        const REAL_SIMD_ARRAY FDPart3tmp423 = FusedMulAddSIMD(
            FDPart3tmp417, MulSIMD(FDPart3tmp71, FDPart3tmp78),
            FusedMulAddSIMD(
                FDPart3tmp387, MulSIMD(FDPart3tmp58, FDPart3tmp80),
                FusedMulAddSIMD(
                    FDPart3tmp390, MulSIMD(FDPart3tmp58, FDPart3tmp71),
                    FusedMulAddSIMD(FDPart3tmp403, FDPart3tmp419,
                                    FusedMulAddSIMD(FDPart3tmp275, MulSIMD(FDPart3tmp378, FDPart3tmp386),
                                                    FusedMulAddSIMD(FDPart3tmp419, MulSIMD(FDPart3tmp68, FDPart3tmp80),
                                                                    FusedMulAddSIMD(FDPart3tmp422, MulSIMD(FDPart3tmp71, FDPart3tmp80),
                                                                                    FusedMulAddSIMD(FDPart3tmp388, FDPart3tmp410,
                                                                                                    MulSIMD(FDPart3tmp391, FDPart3tmp417)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp429 = FusedMulAddSIMD(
            FDPart3tmp417, MulSIMD(FDPart3tmp71, FDPart3tmp86),
            FusedMulAddSIMD(
                FDPart3tmp388, MulSIMD(FDPart3tmp68, FDPart3tmp86),
                FusedMulAddSIMD(
                    FDPart3tmp390, MulSIMD(FDPart3tmp68, FDPart3tmp71),
                    FusedMulAddSIMD(FDPart3tmp377, MulSIMD(FDPart3tmp378, FDPart3tmp73),
                                    FusedMulAddSIMD(FDPart3tmp387, MulSIMD(FDPart3tmp58, FDPart3tmp78),
                                                    FusedMulAddSIMD(FDPart3tmp419, MulSIMD(FDPart3tmp58, FDPart3tmp86),
                                                                    FusedMulAddSIMD(FDPart3tmp422, MulSIMD(FDPart3tmp71, FDPart3tmp78),
                                                                                    FusedMulAddSIMD(FDPart3tmp391, FDPart3tmp422,
                                                                                                    MulSIMD(FDPart3tmp410, FDPart3tmp419)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp431 = FusedMulAddSIMD(
            FDPart3tmp417, MulSIMD(FDPart3tmp58, FDPart3tmp86),
            FusedMulAddSIMD(
                FDPart3tmp387, MulSIMD(FDPart3tmp78, FDPart3tmp80),
                FusedMulAddSIMD(
                    FDPart3tmp388, MulSIMD(FDPart3tmp78, FDPart3tmp86),
                    FusedMulAddSIMD(FDPart3tmp410, FDPart3tmp417,
                                    FusedMulAddSIMD(FDPart3tmp378, MulSIMD(FDPart3tmp379, FDPart3tmp60),
                                                    FusedMulAddSIMD(FDPart3tmp419, MulSIMD(FDPart3tmp80, FDPart3tmp86),
                                                                    FusedMulAddSIMD(FDPart3tmp422, MulSIMD(FDPart3tmp68, FDPart3tmp80),
                                                                                    FusedMulAddSIMD(FDPart3tmp390, FDPart3tmp391,
                                                                                                    MulSIMD(FDPart3tmp403, FDPart3tmp422)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp458 =
            FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp456, FusedMulAddSIMD(FDPart3tmp354, hDD_dD000, MulSIMD(FDPart3tmp174, FDPart3tmp457)));
        const REAL_SIMD_ARRAY FDPart3tmp459 =
            FusedMulAddSIMD(FDPart3tmp202, FDPart3tmp456, FusedMulAddSIMD(FDPart3tmp354, hDD_dD001, MulSIMD(FDPart3tmp194, FDPart3tmp457)));
        const REAL_SIMD_ARRAY FDPart3tmp462 =
            FusedMulAddSIMD(FDPart3tmp231, FDPart3tmp456, FusedMulAddSIMD(FDPart3tmp354, hDD_dD002, MulSIMD(FDPart3tmp227, FDPart3tmp457)));
        const REAL_SIMD_ARRAY FDPart3tmp464 =
            FusedMulAddSIMD(FDPart3tmp257, FDPart3tmp457, FusedMulAddSIMD(FDPart3tmp373, FDPart3tmp456, MulSIMD(FDPart3tmp256, FDPart3tmp354)));
        const REAL_SIMD_ARRAY FDPart3tmp482 =
            FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp481, FusedMulAddSIMD(FDPart3tmp352, hDD_dD000, MulSIMD(FDPart3tmp174, FDPart3tmp456)));
        const REAL_SIMD_ARRAY FDPart3tmp483 =
            FusedMulAddSIMD(FDPart3tmp202, FDPart3tmp481, FusedMulAddSIMD(FDPart3tmp352, hDD_dD001, MulSIMD(FDPart3tmp194, FDPart3tmp456)));
        const REAL_SIMD_ARRAY FDPart3tmp486 =
            FusedMulAddSIMD(FDPart3tmp231, FDPart3tmp481, FusedMulAddSIMD(FDPart3tmp352, hDD_dD002, MulSIMD(FDPart3tmp227, FDPart3tmp456)));
        const REAL_SIMD_ARRAY FDPart3tmp488 =
            FusedMulAddSIMD(FDPart3tmp257, FDPart3tmp456, FusedMulAddSIMD(FDPart3tmp373, FDPart3tmp481, MulSIMD(FDPart3tmp256, FDPart3tmp352)));
        const REAL_SIMD_ARRAY FDPart3tmp495 = NegFusedMulAddSIMD(
            FDPart3_Integer_2, MulSIMD(FDPart3tmp2, vetU0),
            FusedMulAddSIMD(
                FDPart3tmp323, FDPart3tmp330,
                FusedMulAddSIMD(
                    alpha,
                    FusedMulAddSIMD(FDPart3tmp269, FDPart3tmp81,
                                    FusedMulAddSIMD(FDPart3tmp281, FDPart3tmp72,
                                                    FusedMulAddSIMD(FDPart3tmp393, FDPart3tmp59,
                                                                    FusedMulAddSIMD(FDPart3tmp87, aDD00,
                                                                                    FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp84,
                                                                                                    MulSIMD(FDPart3tmp18, FDPart3tmp82)))))),
                    NegFusedMulSubSIMD(FDPart3tmp2, MulSIMD(FDPart3tmp4, vetU_dD22), AddSIMD(vetU_dD00, FDPart3tmp20)))));
        const REAL_SIMD_ARRAY FDPart3tmp517 = FusedMulAddSIMD(
            FDPart3tmp205, alpha_dD1, FusedMulAddSIMD(FDPart3tmp207, alpha_dD0, FusedMulAddSIMD(FDPart3tmp203, alpha_dD2, alpha_dDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp519 = FusedMulAddSIMD(
            FDPart3tmp234, alpha_dD1, FusedMulAddSIMD(FDPart3tmp236, alpha_dD0, FusedMulAddSIMD(FDPart3tmp232, alpha_dD2, alpha_dDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp521 = FusedMulAddSIMD(
            FDPart3tmp250, alpha_dD1, FusedMulAddSIMD(FDPart3tmp252, alpha_dD0, FusedMulAddSIMD(FDPart3tmp248, alpha_dD2, alpha_dDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp265 = FusedMulAddSIMD(
            FDPart3tmp167, FDPart3tmp238,
            FusedMulAddSIMD(FDPart3tmp184, FDPart3tmp87,
                            FusedMulAddSIMD(FDPart3tmp221, FDPart3tmp81,
                                            FusedMulAddSIMD(FDPart3tmp264, FDPart3tmp72,
                                                            FusedMulAddSIMD(FDPart3tmp154, FDPart3tmp209, MulSIMD(FDPart3tmp162, FDPart3tmp254))))));
        const REAL_SIMD_ARRAY FDPart3tmp362 = FusedMulAddSIMD(
            FDPart3tmp214, FDPart3tmp353, FusedMulAddSIMD(FDPart3tmp354, MulSIMD(FDPart3tmp7, hDD_dD111), MulSIMD(FDPart3tmp213, FDPart3tmp352)));
        const REAL_SIMD_ARRAY FDPart3tmp369 = FusedMulAddSIMD(
            FDPart3tmp247, FDPart3tmp353, FusedMulAddSIMD(FDPart3tmp354, MulSIMD(FDPart3tmp7, hDD_dD112), MulSIMD(FDPart3tmp246, FDPart3tmp352)));
        const REAL_SIMD_ARRAY FDPart3tmp385 = AddSIMD(FDPart3tmp316, FDPart3tmp374);
        const REAL_SIMD_ARRAY FDPart3tmp399 = MulSIMD(FDPart3tmp397, FDPart3tmp398);
        const REAL_SIMD_ARRAY FDPart3tmp407 = MulSIMD(FDPart3tmp398, FDPart3tmp406);
        const REAL_SIMD_ARRAY FDPart3tmp424 = MulSIMD(FDPart3tmp382, FDPart3tmp423);
        const REAL_SIMD_ARRAY FDPart3tmp430 = MulSIMD(FDPart3tmp382, FDPart3tmp429);
        const REAL_SIMD_ARRAY FDPart3tmp432 = MulSIMD(FDPart3tmp382, FDPart3tmp431);
        const REAL_SIMD_ARRAY FDPart3tmp461 = FusedMulAddSIMD(
            FDPart3tmp214, FDPart3tmp354, FusedMulAddSIMD(FDPart3tmp457, MulSIMD(FDPart3tmp7, hDD_dD111), MulSIMD(FDPart3tmp213, FDPart3tmp456)));
        const REAL_SIMD_ARRAY FDPart3tmp463 = FusedMulAddSIMD(
            FDPart3tmp247, FDPart3tmp354, FusedMulAddSIMD(FDPart3tmp457, MulSIMD(FDPart3tmp7, hDD_dD112), MulSIMD(FDPart3tmp246, FDPart3tmp456)));
        const REAL_SIMD_ARRAY FDPart3tmp466 = AddSIMD(FDPart3tmp242, FDPart3tmp464);
        const REAL_SIMD_ARRAY FDPart3tmp469 = SubSIMD(FDPart3tmp459, FDPart3tmp2);
        const REAL_SIMD_ARRAY FDPart3tmp485 = FusedMulAddSIMD(
            FDPart3tmp214, FDPart3tmp352, FusedMulAddSIMD(FDPart3tmp456, MulSIMD(FDPart3tmp7, hDD_dD111), MulSIMD(FDPart3tmp213, FDPart3tmp481)));
        const REAL_SIMD_ARRAY FDPart3tmp487 = FusedMulAddSIMD(
            FDPart3tmp247, FDPart3tmp352, FusedMulAddSIMD(FDPart3tmp456, MulSIMD(FDPart3tmp7, hDD_dD112), MulSIMD(FDPart3tmp246, FDPart3tmp481)));
        const REAL_SIMD_ARRAY FDPart3tmp491 = SubSIMD(FDPart3tmp486, FDPart3tmp2);
        const REAL_SIMD_ARRAY FDPart3tmp506 = MulSIMD(FDPart3tmp397, FDPart3tmp505);
        const REAL_SIMD_ARRAY FDPart3tmp507 = MulSIMD(FDPart3tmp406, FDPart3tmp505);
        const REAL_SIMD_ARRAY FDPart3tmp511 = MulSIMD(FDPart3tmp423, FDPart3tmp510);
        const REAL_SIMD_ARRAY FDPart3tmp512 = MulSIMD(FDPart3tmp429, FDPart3tmp510);
        const REAL_SIMD_ARRAY FDPart3tmp513 = MulSIMD(FDPart3tmp431, FDPart3tmp510);
        const REAL_SIMD_ARRAY FDPart3tmp400 = AddSIMD(FDPart3tmp362, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp470 = MulSIMD(
            FDPart3tmp26, FusedMulAddSIMD(FDPart3tmp435, FDPart3tmp462,
                                          FusedMulAddSIMD(FDPart3tmp436, FDPart3tmp463,
                                                          FusedMulAddSIMD(FDPart3tmp437, FDPart3tmp461,
                                                                          FusedMulAddSIMD(FDPart3tmp438, FDPart3tmp466,
                                                                                          FusedMulAddSIMD(FDPart3tmp433, FDPart3tmp458,
                                                                                                          MulSIMD(FDPart3tmp434, FDPart3tmp469)))))));
        const REAL_SIMD_ARRAY FDPart3tmp490 = SubSIMD(FDPart3tmp487, FDPart3tmp323);
        const REAL_SIMD_ARRAY FDPart3tmp439 = MulSIMD(
            FDPart3tmp26, FusedMulAddSIMD(FDPart3tmp365, FDPart3tmp435,
                                          FusedMulAddSIMD(FDPart3tmp369, FDPart3tmp436,
                                                          FusedMulAddSIMD(FDPart3tmp385, FDPart3tmp438,
                                                                          FusedMulAddSIMD(FDPart3tmp400, FDPart3tmp437,
                                                                                          FusedMulAddSIMD(FDPart3tmp355, FDPart3tmp433,
                                                                                                          MulSIMD(FDPart3tmp357, FDPart3tmp434)))))));
        const REAL_SIMD_ARRAY FDPart3tmp492 = MulSIMD(
            FDPart3tmp26, FusedMulAddSIMD(FDPart3tmp435, FDPart3tmp491,
                                          FusedMulAddSIMD(FDPart3tmp436, FDPart3tmp490,
                                                          FusedMulAddSIMD(FDPart3tmp437, FDPart3tmp485,
                                                                          FusedMulAddSIMD(FDPart3tmp438, FDPart3tmp488,
                                                                                          FusedMulAddSIMD(FDPart3tmp433, FDPart3tmp482,
                                                                                                          MulSIMD(FDPart3tmp434, FDPart3tmp483)))))));
        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(
            MulSIMD(aDD01, f0_of_xx0),
            MulSIMD(
                MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                MulSIMD(alpha, FusedMulAddSIMD(FDPart3tmp60, FDPart3tmp81,
                                               FusedMulAddSIMD(FDPart3tmp52, MulSIMD(FDPart3tmp78, aDD00), MulSIMD(FDPart3tmp59, FDPart3tmp73))))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp89, PI),
                MulSIMD(
                    MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                    MulSIMD(alpha, FusedMulAddSIMD(
                                       FDPart3tmp155, FDPart3tmp92,
                                       FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp92,
                                                       FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp92,
                                                                       FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp92,
                                                                                       FusedMulAddSIMD(FDPart3tmp168, FDPart3tmp92,
                                                                                                       FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp92,
                                                                                                                       FDPart3tmp139)))))))),
                FusedMulAddSIMD(
                    aDD00, MulSIMD(alpha, trK),
                    FusedMulAddSIMD(
                        aDD00,
                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                MulSIMD(alpha, FusedMulAddSIMD(FDPart3tmp73, FDPart3tmp82,
                                                               FusedMulAddSIMD(FDPart3tmp87, aDD00, MulSIMD(FDPart3tmp60, FDPart3tmp84))))),
                        FusedMulAddSIMD(
                            aDD_dupD000, vetU0,
                            FusedMulAddSIMD(
                                FDPart3_Integer_2, MulSIMD(aDD00, vetU_dD00),
                                FusedMulAddSIMD(
                                    FDPart3tmp3, aDD_dupD001,
                                    FusedMulAddSIMD(
                                        FDPart3tmp6, aDD_dupD002,
                                        FusedMulAddSIMD(
                                            FDPart3tmp16, FDPart3tmp18,
                                            FusedMulAddSIMD(
                                                FDPart3tmp169,
                                                NegFusedMulAddSIMD(FDPart3tmp265, FusedMulAddSIMD(FDPart3_Rational_1_3, hDD00, FDPart3_Rational_1_3),
                                                                   FDPart3tmp184),
                                                FusedMulAddSIMD(
                                                    f1_of_xx1,
                                                    MulSIMD(MulSIMD(aDD02, f0_of_xx0),
                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                    MulSIMD(alpha, FusedMulAddSIMD(
                                                                                       FDPart3tmp72, FDPart3tmp73,
                                                                                       FusedMulAddSIMD(FDPart3tmp52, MulSIMD(FDPart3tmp68, aDD00),
                                                                                                       MulSIMD(FDPart3tmp59, FDPart3tmp60)))))),
                                                    FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp13,
                                                                    MulSIMD(FDPart3_Rational_2_3, MulSIMD(FDPart3tmp26, aDD00))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_1 = MulSIMD(
            FDPart3tmp2,
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp89, PI),
                MulSIMD(
                    MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                    MulSIMD(alpha,
                            FusedMulAddSIMD(
                                FDPart3tmp78,
                                MulSIMD(MulSIMD(FDPart3tmp152, FDPart3tmp52),
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0, hDD01))),
                                FusedMulAddSIMD(
                                    FDPart3tmp58,
                                    MulSIMD(MulSIMD(FDPart3tmp161, FDPart3tmp52),
                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0, hDD01))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp80,
                                        MulSIMD(MulSIMD(FDPart3tmp130, FDPart3tmp52),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f0_of_xx0, hDD01))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp86,
                                            MulSIMD(MulSIMD(FDPart3tmp139, FDPart3tmp52),
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f0_of_xx0, hDD01))),
                                            FusedMulAddSIMD(
                                                FDPart3tmp68,
                                                MulSIMD(MulSIMD(FDPart3tmp166, FDPart3tmp52),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0, hDD01))),
                                                NegFusedMulAddSIMD(MulSIMD(FDPart3tmp52, FDPart3tmp71),
                                                                   MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp118), MulSIMD(f0_of_xx0, hDD01)),
                                                                   FDPart3tmp152)))))))),
                FusedMulAddSIMD(
                    FDPart3tmp26, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(aDD01, f0_of_xx0)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp279, aDD01), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, alpha)),
                        FusedMulAddSIMD(
                            aDD_dupD011, vetU1,
                            FusedMulAddSIMD(
                                vetU0, FusedMulAddSIMD(aDD_dupD010, f0_of_xx0, aDD01),
                                FusedMulAddSIMD(
                                    aDD00, vetU_dD01,
                                    FusedMulAddSIMD(
                                        aDD01, vetU_dD11,
                                        FusedMulAddSIMD(
                                            FDPart3tmp268, FDPart3tmp60,
                                            FusedMulAddSIMD(
                                                FDPart3tmp274, FDPart3tmp73,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp169,
                                                    NegFusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp265), MulSIMD(f0_of_xx0, hDD01),
                                                                       FDPart3tmp209),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp267, aDD01,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp14, aDD_dupD012,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp16, FDPart3tmp275,
                                                                FusedMulAddSIMD(f0_of_xx0,
                                                                                MulSIMD(MulSIMD(FDPart3tmp276, aDD02),
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(f1_of_xx1, alpha))),
                                                                                FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp269,
                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp280),
                                                                                                        MulSIMD(aDD00, alpha))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_2 = MulSIMD(
            FDPart3tmp5,
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp26, aDD02), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f0_of_xx0, f1_of_xx1)),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp284, aDD01), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, alpha)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp89, PI),
                        MulSIMD(
                            MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                            MulSIMD(alpha,
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3tmp78, f0_of_xx0),
                                        MulSIMD(MulSIMD(FDPart3tmp152, FDPart3tmp52),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f1_of_xx1, hDD02))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3tmp58, f0_of_xx0),
                                            MulSIMD(MulSIMD(FDPart3tmp161, FDPart3tmp52),
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f1_of_xx1, hDD02))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp80, f0_of_xx0),
                                                MulSIMD(MulSIMD(FDPart3tmp130, FDPart3tmp52),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f1_of_xx1, hDD02))),
                                                FusedMulAddSIMD(
                                                    MulSIMD(FDPart3tmp86, f0_of_xx0),
                                                    MulSIMD(MulSIMD(FDPart3tmp139, FDPart3tmp52),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f1_of_xx1, hDD02))),
                                                    FusedMulAddSIMD(MulSIMD(FDPart3tmp68, f0_of_xx0),
                                                                    MulSIMD(MulSIMD(FDPart3tmp166, FDPart3tmp52),
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3),
                                                                                    MulSIMD(f1_of_xx1, hDD02))),
                                                                    NegFusedMulAddSIMD(f0_of_xx0,
                                                                                       MulSIMD(MulSIMD(FDPart3tmp52, FDPart3tmp71),
                                                                                               MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp118),
                                                                                                       MulSIMD(f1_of_xx1, hDD02))),
                                                                                       FDPart3tmp166)))))))),
                        FusedMulAddSIMD(
                            aDD_dupD022, vetU2,
                            FusedMulAddSIMD(
                                vetU0, FusedMulAddSIMD(FDPart3tmp109, aDD_dupD020, FDPart3tmp17),
                                FusedMulAddSIMD(
                                    aDD01, vetU_dD12,
                                    FusedMulAddSIMD(
                                        aDD02, vetU_dD22,
                                        FusedMulAddSIMD(
                                            FDPart3tmp3,
                                            FusedMulAddSIMD(FDPart3tmp109, aDD_dupD021, MulSIMD(aDD02, MulSIMD(f0_of_xx0, f1_of_xx1__D1))),
                                            FusedMulAddSIMD(
                                                aDD00, vetU_dD02,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp17, FDPart3tmp267,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp268, FDPart3tmp73,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp16, FDPart3tmp281,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp169,
                                                                NegFusedMulAddSIMD(
                                                                    f0_of_xx0,
                                                                    MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp265), MulSIMD(f1_of_xx1, hDD02)),
                                                                    FDPart3tmp238),
                                                                FusedMulAddSIMD(f0_of_xx0,
                                                                                MulSIMD(MulSIMD(FDPart3tmp283, aDD02),
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(f1_of_xx1, alpha))),
                                                                                FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp275,
                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp285),
                                                                                                        MulSIMD(aDD00, alpha))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_3 = MulSIMD(
            FDPart3tmp8,
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp280, aDD01), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, alpha)),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp89, PI),
                    MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                            MulSIMD(alpha,
                                    FusedMulAddSIMD(
                                        FDPart3tmp155, FDPart3tmp290,
                                        FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp290,
                                                        FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp290,
                                                                        FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp290,
                                                                                        FusedMulAddSIMD(FDPart3tmp168, FDPart3tmp290,
                                                                                                        FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp290,
                                                                                                                        FDPart3tmp130)))))))),
                    FusedMulAddSIMD(
                        aDD_dupD111, MulSIMD(f0_of_xx0, vetU1),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp279, FDPart3tmp7), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD11, alpha)),
                            FusedMulAddSIMD(
                                FDPart3tmp12, MulSIMD(aDD11, vetU_dD11),
                                FusedMulAddSIMD(
                                    FDPart3tmp14, MulSIMD(aDD_dupD112, f0_of_xx0),
                                    FusedMulAddSIMD(
                                        vetU0, FusedMulAddSIMD(FDPart3tmp12, aDD11, MulSIMD(FDPart3tmp7, aDD_dupD110)),
                                        FusedMulAddSIMD(
                                            FDPart3tmp103, MulSIMD(FDPart3tmp274, aDD12),
                                            FusedMulAddSIMD(
                                                FDPart3tmp169,
                                                NegFusedMulAddSIMD(
                                                    FDPart3tmp265,
                                                    FusedMulAddSIMD(FDPart3_Rational_1_3, FDPart3tmp37, MulSIMD(FDPart3_Rational_1_3, FDPart3tmp7)),
                                                    FDPart3tmp221),
                                                FusedMulAddSIMD(FDPart3tmp268, FDPart3tmp269,
                                                                FusedMulAddSIMD(aDD12,
                                                                                MulSIMD(MulSIMD(FDPart3tmp276, FDPart3tmp7),
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(f1_of_xx1, alpha))),
                                                                                FusedMulSubSIMD(FDPart3tmp13, vetU_dD01,
                                                                                                MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp26),
                                                                                                        MulSIMD(FDPart3tmp7, aDD11)))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_4 = MulSIMD(
            FDPart3tmp296,
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp26, FDPart3tmp7), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(aDD12, f1_of_xx1)),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp285, aDD01), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, alpha)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp89, PI),
                        MulSIMD(
                            MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                            MulSIMD(alpha,
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3tmp7, FDPart3tmp78),
                                        MulSIMD(MulSIMD(FDPart3tmp152, FDPart3tmp52),
                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f1_of_xx1, hDD12))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3tmp58, FDPart3tmp7),
                                            MulSIMD(MulSIMD(FDPart3tmp161, FDPart3tmp52),
                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(f1_of_xx1, hDD12))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp7, FDPart3tmp80),
                                                MulSIMD(MulSIMD(FDPart3tmp130, FDPart3tmp52),
                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f1_of_xx1, hDD12))),
                                                FusedMulAddSIMD(
                                                    MulSIMD(FDPart3tmp7, FDPart3tmp86),
                                                    MulSIMD(MulSIMD(FDPart3tmp139, FDPart3tmp52),
                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_3), MulSIMD(f1_of_xx1, hDD12))),
                                                    FusedMulAddSIMD(MulSIMD(FDPart3tmp68, FDPart3tmp7),
                                                                    MulSIMD(MulSIMD(FDPart3tmp166, FDPart3tmp52),
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3),
                                                                                    MulSIMD(f1_of_xx1, hDD12))),
                                                                    NegFusedMulAddSIMD(FDPart3tmp71,
                                                                                       MulSIMD(MulSIMD(FDPart3tmp52, FDPart3tmp7),
                                                                                               MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp118),
                                                                                                       MulSIMD(f1_of_xx1, hDD12))),
                                                                                       FDPart3tmp161)))))))),
                        FusedMulAddSIMD(
                            aDD12, MulSIMD(f0_of_xx0, vetU_dD22),
                            FusedMulAddSIMD(
                                aDD_dupD122, MulSIMD(f0_of_xx0, vetU2),
                                FusedMulAddSIMD(
                                    aDD01, MulSIMD(f0_of_xx0, vetU_dD02),
                                    FusedMulAddSIMD(
                                        aDD11, MulSIMD(f0_of_xx0, vetU_dD12),
                                        FusedMulAddSIMD(
                                            vetU0, FusedMulAddSIMD(FDPart3tmp102, aDD_dupD120, MulSIMD(FDPart3tmp175, aDD12)),
                                            FusedMulAddSIMD(
                                                FDPart3tmp109, MulSIMD(aDD12, vetU_dD11),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp3,
                                                    FusedMulAddSIMD(FDPart3tmp102, aDD_dupD121, MulSIMD(FDPart3tmp7, MulSIMD(aDD12, f1_of_xx1__D1))),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp73, vetU_dD01,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp268, FDPart3tmp275,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp274, FDPart3tmp281,
                                                                FusedMulAddSIMD(
                                                                    aDD12,
                                                                    MulSIMD(MulSIMD(FDPart3tmp283, FDPart3tmp7),
                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                    MulSIMD(f1_of_xx1, alpha))),
                                                                    FusedMulSubSIMD(
                                                                        FDPart3tmp169,
                                                                        NegFusedMulAddSIMD(FDPart3tmp7,
                                                                                           MulSIMD(MulSIMD(FDPart3_Rational_1_3, FDPart3tmp265),
                                                                                                   MulSIMD(f1_of_xx1, hDD12)),
                                                                                           FDPart3tmp254),
                                                                        MulSIMD(FDPart3tmp7, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp284),
                                                                                                     MulSIMD(aDD11, alpha)))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_5 = MulSIMD(
            FDPart3tmp299,
            FusedMulAddSIMD(
                FDPart3tmp7,
                MulSIMD(MulSIMD(FDPart3tmp283, FDPart3tmp29), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(aDD22, alpha))),
                FusedMulAddSIMD(
                    aDD12,
                    MulSIMD(MulSIMD(FDPart3tmp284, FDPart3tmp7),
                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, alpha))),
                    FusedMulAddSIMD(
                        FDPart3tmp175, MulSIMD(aDD22, vetU_dD22),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp26, FDPart3tmp29),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp7, aDD22)),
                            FusedMulAddSIMD(
                                FDPart3tmp109, MulSIMD(aDD_dupD222, vetU2),
                                FusedMulAddSIMD(
                                    FDPart3tmp175, MulSIMD(aDD12, vetU_dD12),
                                    FusedMulAddSIMD(
                                        FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp244, aDD22, MulSIMD(FDPart3tmp39, aDD_dupD221)),
                                        FusedMulAddSIMD(
                                            vetU0, FusedMulAddSIMD(FDPart3tmp229, aDD22, MulSIMD(FDPart3tmp39, aDD_dupD220)),
                                            FusedMulAddSIMD(
                                                FDPart3tmp18, vetU_dD02,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp268, FDPart3tmp281,
                                                    FusedMulAddSIMD(
                                                        f0_of_xx0,
                                                        MulSIMD(MulSIMD(FDPart3tmp285, aDD02),
                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, alpha))),
                                                        FusedMulSubSIMD(
                                                            FDPart3tmp169,
                                                            NegFusedMulAddSIMD(FDPart3tmp265,
                                                                               FusedMulAddSIMD(FDPart3_Rational_1_3, FDPart3tmp39,
                                                                                               MulSIMD(FDPart3_Rational_1_3, FDPart3tmp43)),
                                                                               FDPart3tmp264),
                                                            MulSIMD(
                                                                PI,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3tmp89),
                                                                        MulSIMD(alpha, FusedMulAddSIMD(
                                                                                           FDPart3tmp155, FDPart3tmp298,
                                                                                           FusedMulAddSIMD(
                                                                                               FDPart3tmp131, FDPart3tmp298,
                                                                                               FusedMulAddSIMD(
                                                                                                   FDPart3tmp140, FDPart3tmp298,
                                                                                                   FusedMulAddSIMD(
                                                                                                       FDPart3tmp163, FDPart3tmp298,
                                                                                                       FusedMulAddSIMD(
                                                                                                           FDPart3tmp168, FDPart3tmp298,
                                                                                                           FusedMulAddSIMD(
                                                                                                               FDPart3tmp119, FDPart3tmp298,
                                                                                                               FDPart3tmp118))))))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_6 = FusedMulAddSIMD(
            FDPart3tmp6, alpha_dupD2,
            FusedMulAddSIMD(alpha_dupD0, vetU0, FusedMulSubSIMD(FDPart3tmp3, alpha_dupD1, MulSIMD(FDPart3_Integer_2, MulSIMD(alpha, trK)))));
const REAL_SIMD_ARRAY __RHS_exp_7 = FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp369, FDPart3tmp4), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))), FusedMulAddSIMD(FDPart3tmp4, MulSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp365), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp272, FDPart3tmp374), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2))), FusedMulAddSIMD(FDPart3tmp68, MulSIMD(MulSIMD(FDPart3tmp332, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp4, MulSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp365), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))), FusedMulAddSIMD(FDPart3tmp86, MulSIMD(MulSIMD(FDPart3tmp329, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp78, MulSIMD(MulSIMD(FDPart3tmp331, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp4), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU_dD02)), FusedMulAddSIMD(MulSIMD(FDPart3tmp362, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)), FusedMulAddSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp357), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1)), FusedMulAddSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp357), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0)), FusedMulAddSIMD(FDPart3tmp2, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU_dD01)), FusedMulAddSIMD(FDPart3tmp355, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)), NegFusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(lambdaU0, vetU_dD00), NegFusedMulAddSIMD(FDPart3tmp87, MulSIMD(alpha, trK_dD0), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp87, vetU_dDD000), FusedMulAddSIMD(FDPart3tmp355, MulSIMD(FDPart3tmp398, FDPart3tmp411), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp72, NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f1_of_xx1, vetU_dD22), FusedMulAddSIMD(FDPart3tmp242, vetU_dD01, FusedMulAddSIMD(FDPart3tmp316, vetU_dD00, FusedMulAddSIMD(f1_of_xx1, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, vetU1)), NegFusedMulAddSIMD(FDPart3tmp29, vetU0, vetU_dDD022)))))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp81, AddSIMD(vetU_dDD011, NegFusedMulAddSIMD(FDPart3_Integer_2, vetU_dD11, SubSIMD(FDPart3tmp267, vetU0)))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp413, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD2, cf), alpha_dD2)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp415, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD1, cf), alpha_dD1)), FusedMulAddSIMD(FDPart3_Rational_3_2, MulSIMD(FDPart3tmp325, FDPart3tmp59), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp381, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD0, cf), alpha_dD0)), FusedMulAddSIMD(FDPart3tmp6, betU_dupD02, FusedMulAddSIMD(betU_dupD00, vetU0, FusedMulAddSIMD(FDPart3tmp385, FDPart3tmp399, FusedMulAddSIMD(FDPart3tmp400, FDPart3tmp407, FusedMulAddSIMD(FDPart3tmp369, FDPart3tmp424, FusedMulAddSIMD(FDPart3tmp374, FDPart3tmp376, FusedMulAddSIMD(FDPart3tmp369, FDPart3tmp371, FusedMulAddSIMD(FDPart3tmp369, FDPart3tmp372, FusedMulAddSIMD(FDPart3tmp365, FDPart3tmp368, FusedMulAddSIMD(FDPart3tmp365, FDPart3tmp430, FusedMulAddSIMD(FDPart3tmp362, FDPart3tmp364, FusedMulAddSIMD(FDPart3tmp365, FDPart3tmp366, FusedMulAddSIMD(FDPart3tmp357, FDPart3tmp360, FusedMulAddSIMD(FDPart3tmp357, FDPart3tmp432, FusedMulAddSIMD(FDPart3tmp355, FDPart3tmp356, FusedMulAddSIMD(FDPart3tmp357, FDPart3tmp358, FusedMulAddSIMD(FDPart3tmp345, FDPart3tmp84, FusedMulAddSIMD(FDPart3tmp351, FDPart3tmp87, FusedMulAddSIMD(FDPart3tmp314, FDPart3tmp315, FusedMulAddSIMD(FDPart3tmp320, FDPart3tmp321, FusedMulAddSIMD(FDPart3tmp303, trK_dD1, FusedMulAddSIMD(FDPart3tmp311, FDPart3tmp82, FusedMulAddSIMD(FDPart3tmp3, betU_dupD01, FusedMulAddSIMD(FDPart3tmp302, FDPart3tmp82, FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp369, FDPart3tmp4), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp439, MulSIMD(betU0, eta))))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_8 = MulSIMD(f0_of_xx0, FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp272, FDPart3tmp464), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp4, FDPart3tmp463), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))), FusedMulAddSIMD(FDPart3tmp462, MulSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp4), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))), FusedMulAddSIMD(FDPart3tmp462, MulSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp4), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))), FusedMulAddSIMD(FDPart3tmp80, MulSIMD(MulSIMD(FDPart3tmp331, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp58, MulSIMD(MulSIMD(FDPart3tmp332, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp461, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)), FusedMulAddSIMD(FDPart3tmp78, MulSIMD(MulSIMD(FDPart3tmp329, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp459), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0)), FusedMulAddSIMD(MulSIMD(FDPart3tmp4, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU_dD12)), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU_dD11)), FusedMulAddSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp459), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1)), FusedMulAddSIMD(FDPart3tmp398, MulSIMD(FDPart3tmp411, FDPart3tmp458), FusedMulAddSIMD(FDPart3tmp458, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp29, FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp2, vetU_dDD122, FusedMulAddSIMD(FDPart3tmp20, FDPart3tmp242, FusedMulAddSIMD(FDPart3tmp330, FDPart3tmp337, FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp316, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp2), MulSIMD(f1_of_xx1__D1, vetU_dD22)))))))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp87, FusedMulAddSIMD(FDPart3tmp2, vetU_dDD100, FusedMulAddSIMD(FDPart3tmp450, vetU1, FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp21, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp8, vetU_dD10)))))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp468, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD2, cf), alpha_dD2)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp59, FusedMulAddSIMD(FDPart3tmp337, FDPart3tmp6, FusedMulAddSIMD(FDPart3tmp6, NegFusedMulSubSIMD(f1_of_xx1, f1_of_xx1__DD11, FDPart3tmp337), FusedMulAddSIMD(FDPart3tmp271, f1_of_xx1, FDPart3tmp448)))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp448, FDPart3tmp59), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp465, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD1, cf), alpha_dD1)), FusedMulAddSIMD(vetU0, FusedMulSubSIMD(FDPart3tmp2, betU_dupD10, MulSIMD(FDPart3tmp8, betU1)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp415, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD0, cf), alpha_dD0)), FusedMulAddSIMD(FDPart3tmp440, betU1, FusedMulAddSIMD(FDPart3tmp441, FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp430, FDPart3tmp462, FusedMulAddSIMD(FDPart3tmp432, FDPart3tmp469, FusedMulAddSIMD(FDPart3tmp407, FDPart3tmp461, FusedMulAddSIMD(FDPart3tmp424, FDPart3tmp463, FusedMulAddSIMD(FDPart3tmp376, FDPart3tmp464, FusedMulAddSIMD(FDPart3tmp399, FDPart3tmp466, FusedMulAddSIMD(FDPart3tmp371, FDPart3tmp463, FusedMulAddSIMD(FDPart3tmp372, FDPart3tmp463, FusedMulAddSIMD(FDPart3tmp368, FDPart3tmp462, FusedMulAddSIMD(FDPart3tmp370, betU_dupD12, FusedMulAddSIMD(FDPart3tmp364, FDPart3tmp461, FusedMulAddSIMD(FDPart3tmp366, FDPart3tmp462, FusedMulAddSIMD(FDPart3tmp360, FDPart3tmp459, FusedMulAddSIMD(FDPart3tmp363, betU_dupD11, FusedMulAddSIMD(FDPart3tmp356, FDPart3tmp458, FusedMulAddSIMD(FDPart3tmp358, FDPart3tmp459, FusedMulAddSIMD(FDPart3tmp345, FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp351, FDPart3tmp84, FusedMulAddSIMD(FDPart3tmp315, AddSIMD(FDPart3tmp347, FDPart3tmp442), FusedMulAddSIMD(FDPart3tmp321, FDPart3tmp444, FusedMulAddSIMD(FDPart3tmp303, trK_dD0, FusedMulAddSIMD(FDPart3tmp311, FDPart3tmp59, FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp470, FusedMulAddSIMD(FDPart3tmp302, FDPart3tmp59, FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp4, FDPart3tmp463), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp446, MulSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp11, lambdaU0)))))))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_9 = MulSIMD(FDPart3tmp109, FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp272, FDPart3tmp488), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU2))), FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp4, FDPart3tmp487), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU2))), FusedMulAddSIMD(FDPart3tmp486, MulSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp4), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU2))), FusedMulAddSIMD(FDPart3tmp486, MulSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp4), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU0))), FusedMulAddSIMD(FDPart3tmp58, MulSIMD(MulSIMD(FDPart3tmp331, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(FDPart3tmp71, MulSIMD(MulSIMD(FDPart3tmp332, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp485, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU1)), FusedMulAddSIMD(FDPart3tmp68, MulSIMD(MulSIMD(FDPart3tmp329, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_12, FDPart3_NegativeOne_), MulSIMD(PI, alpha))), FusedMulAddSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp483), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU1, vetU0)), FusedMulAddSIMD(MulSIMD(FDPart3tmp272, FDPart3tmp8), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU_dD22)), FusedMulAddSIMD(FDPart3tmp482, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU0)), FusedMulAddSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp483), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU0, vetU1)), NegFusedMulAddSIMD(FDPart3tmp82, MulSIMD(alpha, trK_dD0), FusedMulAddSIMD(FDPart3tmp2, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(FDPart3tmp274, lambdaU1)), FusedMulAddSIMD(FDPart3tmp398, MulSIMD(FDPart3tmp411, FDPart3tmp482), FusedMulAddSIMD(FDPart3tmp4, MulSIMD(FDPart3tmp440, betU2), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp81, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp274, FDPart3tmp323), FusedMulAddSIMD(FDPart3tmp5, vetU_dDD211, FusedMulAddSIMD(vetU2, FusedMulAddSIMD(FDPart3tmp21, FDPart3tmp479, MulSIMD(FDPart3tmp341, f1_of_xx1__DD11)), FusedMulAddSIMD(FDPart3tmp16, f0_of_xx0, FusedMulAddSIMD(FDPart3tmp474, FDPart3tmp6, FusedMulAddSIMD(FDPart3tmp2, MulSIMD(FDPart3tmp479, vetU2), NegFusedMulAddSIMD(FDPart3tmp272, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp2), MulSIMD(f1_of_xx1__D1, vetU_dD21)), FDPart3tmp6)))))))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp87, FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp21, FusedMulAddSIMD(FDPart3tmp5, vetU_dDD200, FusedMulSubSIMD(FDPart3tmp14, FDPart3tmp450, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp4), MulSIMD(FDPart3tmp8, vetU_dD20)))))), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp468, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD1, cf), alpha_dD1)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp489, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD2, cf), alpha_dD2)), FusedMulAddSIMD(vetU0, FusedMulAddSIMD(FDPart3tmp348, betU2, MulSIMD(FDPart3tmp5, betU_dupD20)), FusedMulAddSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp413, FusedMulAddSIMD(FDPart3tmp382, DivSIMD(cf_dD0, cf), alpha_dD0)), FusedMulAddSIMD(FDPart3tmp432, FDPart3tmp483, FusedMulAddSIMD(FDPart3tmp441, FDPart3tmp59, FusedMulAddSIMD(FDPart3tmp424, FDPart3tmp490, FusedMulAddSIMD(FDPart3tmp430, FDPart3tmp491, FusedMulAddSIMD(FDPart3tmp399, FDPart3tmp488, FusedMulAddSIMD(FDPart3tmp407, FDPart3tmp485, FusedMulAddSIMD(FDPart3tmp375, betU_dupD22, FusedMulAddSIMD(FDPart3tmp376, FDPart3tmp488, FusedMulAddSIMD(FDPart3tmp371, FDPart3tmp487, FusedMulAddSIMD(FDPart3tmp372, FDPart3tmp487, FusedMulAddSIMD(FDPart3tmp366, FDPart3tmp486, FusedMulAddSIMD(FDPart3tmp368, FDPart3tmp486, FusedMulAddSIMD(FDPart3tmp360, FDPart3tmp483, FusedMulAddSIMD(FDPart3tmp364, FDPart3tmp485, FusedMulAddSIMD(FDPart3tmp356, FDPart3tmp482, FusedMulAddSIMD(FDPart3tmp358, FDPart3tmp483, FusedMulAddSIMD(FDPart3tmp345, FDPart3tmp59, FusedMulAddSIMD(FDPart3tmp351, FDPart3tmp82, FusedMulAddSIMD(FDPart3tmp315, FDPart3tmp476, FusedMulAddSIMD(FDPart3tmp321, FDPart3tmp471, FusedMulAddSIMD(FDPart3tmp302, FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp311, FDPart3tmp72, FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp492, FusedMulAddSIMD(FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp342, betU2, MulSIMD(FDPart3tmp5, betU_dupD21)), FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp475, FusedMulAddSIMD(FDPart3_Rational_3_4, FDPart3tmp478, FusedMulAddSIMD(FDPart3tmp8, MulSIMD(MulSIMD(FDPart3tmp4, FDPart3tmp487), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_3_4), MulSIMD(lambdaU2, vetU1))), FusedMulSubSIMD(FDPart3_Rational_3_4, FDPart3tmp473, MulSIMD(FDPart3_Rational_3_4, MulSIMD(FDPart3tmp16, lambdaU0)))))))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_10 = MulSIMD(
    MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
    MulSIMD(cf, NegFusedMulAddSIMD(
                    FDPart3_Rational_1_6, MulSIMD(alpha, trK),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp172, FDPart3tmp2), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(cf_dupD1, vetU1)),
                        FusedMulAddSIMD(
                            FDPart3_Rational_1_6, MulSIMD(FDPart3tmp24, vetU1),
                            FusedMulAddSIMD(
                                FDPart3_Rational_1_6, MulSIMD(FDPart3tmp5, vetU_dD22),
                                FusedMulAddSIMD(FDPart3_Rational_1_6, vetU_dD00,
                                                FusedMulAddSIMD(FDPart3_Rational_1_3, MulSIMD(FDPart3tmp2, vetU0),
                                                                FusedMulAddSIMD(FDPart3tmp4,
                                                                                MulSIMD(MulSIMD(FDPart3tmp172, FDPart3tmp2),
                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                MulSIMD(cf_dupD2, vetU2))),
                                                                                FusedMulSubSIMD(FDPart3_Rational_1_6, FDPart3tmp20,
                                                                                                MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp172),
                                                                                                        MulSIMD(cf_dupD0, vetU0))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_11 = FusedMulAddSIMD(
    hDD_dupD000, vetU0,
    FusedMulAddSIMD(
        FDPart3tmp495, FusedMulAddSIMD(FDPart3_Rational_2_3, hDD00, FDPart3_Rational_2_3),
        FusedMulAddSIMD(FDPart3tmp6, hDD_dupD002,
                        FusedMulAddSIMD(FDPart3tmp28, FDPart3tmp313,
                                        FusedMulAddSIMD(FDPart3tmp3, hDD_dupD001,
                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp33, vetU_dD00),
                                                                        FusedMulSubSIMD(FDPart3tmp16, FDPart3tmp99,
                                                                                        MulSIMD(FDPart3_Integer_2, MulSIMD(aDD00, alpha)))))))));
const REAL_SIMD_ARRAY __RHS_exp_12 = MulSIMD(
    FDPart3tmp2,
    FusedMulAddSIMD(
        hDD_dupD011, vetU1,
        FusedMulAddSIMD(
            vetU0, FusedMulAddSIMD(f0_of_xx0, hDD_dupD010, hDD01),
            FusedMulAddSIMD(
                FDPart3tmp33, vetU_dD01,
                FusedMulAddSIMD(
                    hDD01, vetU_dD11,
                    FusedMulAddSIMD(
                        FDPart3tmp267, hDD01,
                        FusedMulAddSIMD(FDPart3tmp274, FDPart3tmp65,
                                        FusedMulAddSIMD(FDPart3tmp14, hDD_dupD012,
                                                        FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp16,
                                                                        FusedMulAddSIMD(FDPart3_Rational_2_3, MulSIMD(FDPart3tmp495, FDPart3tmp75),
                                                                                        FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp38,
                                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, aDD01),
                                                                                                                MulSIMD(f0_of_xx0, alpha)))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_13 = MulSIMD(
    FDPart3tmp5,
    FusedMulAddSIMD(
        hDD_dupD022, vetU2,
        FusedMulAddSIMD(
            vetU0, FusedMulAddSIMD(FDPart3tmp109, hDD_dupD020, FDPart3tmp53),
            FusedMulAddSIMD(
                hDD01, vetU_dD12,
                FusedMulAddSIMD(
                    hDD02, vetU_dD22,
                    FusedMulAddSIMD(
                        FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp109, hDD_dupD021, FDPart3tmp197),
                        FusedMulAddSIMD(
                            FDPart3tmp33, vetU_dD02,
                            FusedMulAddSIMD(
                                FDPart3tmp16, FDPart3tmp44,
                                FusedMulAddSIMD(FDPart3tmp267, FDPart3tmp53,
                                                FusedMulAddSIMD(FDPart3_Rational_2_3, MulSIMD(FDPart3tmp495, FDPart3tmp65),
                                                                FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp142,
                                                                                MulSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, aDD02),
                                                                                                           MulSIMD(f1_of_xx1, alpha))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_14 = MulSIMD(
    FDPart3tmp8,
    FusedMulAddSIMD(
        FDPart3tmp14, MulSIMD(f0_of_xx0, hDD_dupD112),
        FusedMulAddSIMD(
            FDPart3_Integer_2, MulSIMD(FDPart3tmp20, FDPart3tmp38),
            FusedMulAddSIMD(
                FDPart3tmp103, MulSIMD(FDPart3tmp274, hDD12),
                FusedMulAddSIMD(
                    FDPart3tmp495, FusedMulAddSIMD(FDPart3_Rational_2_3, FDPart3tmp37, MulSIMD(FDPart3_Rational_2_3, FDPart3tmp7)),
                    FusedMulAddSIMD(vetU0, FusedMulAddSIMD(FDPart3tmp7, hDD_dupD110, FDPart3tmp193),
                                    FusedMulAddSIMD(f0_of_xx0, MulSIMD(hDD_dupD111, vetU1),
                                                    FusedMulSubSIMD(FDPart3tmp121, vetU_dD01,
                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp7), MulSIMD(aDD11, alpha))))))))));
const REAL_SIMD_ARRAY __RHS_exp_15 = MulSIMD(
    FDPart3tmp296,
    FusedMulAddSIMD(
        f0_of_xx0, MulSIMD(hDD12, vetU_dD22),
        FusedMulAddSIMD(
            f0_of_xx0, MulSIMD(hDD_dupD122, vetU2),
            FusedMulAddSIMD(
                FDPart3tmp2, MulSIMD(FDPart3tmp38, vetU_dD12),
                FusedMulAddSIMD(
                    f0_of_xx0, MulSIMD(hDD01, vetU_dD02),
                    FusedMulAddSIMD(
                        vetU0, FusedMulAddSIMD(FDPart3tmp102, hDD_dupD120, FDPart3tmp200),
                        FusedMulAddSIMD(
                            FDPart3tmp109, MulSIMD(hDD12, vetU_dD11),
                            FusedMulAddSIMD(
                                FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp102, hDD_dupD121, FDPart3tmp212),
                                FusedMulAddSIMD(FDPart3tmp65, vetU_dD01,
                                                FusedMulAddSIMD(FDPart3tmp7,
                                                                MulSIMD(MulSIMD(FDPart3_Rational_2_3, FDPart3tmp495), MulSIMD(f1_of_xx1, hDD12)),
                                                                FusedMulSubSIMD(FDPart3tmp274, FDPart3tmp44,
                                                                                MulSIMD(aDD12, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp7),
                                                                                                       MulSIMD(f1_of_xx1, alpha))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_16 = MulSIMD(
    FDPart3tmp299,
    FusedMulAddSIMD(
        FDPart3tmp109, MulSIMD(hDD_dupD222, vetU2),
        FusedMulAddSIMD(
            FDPart3tmp99, vetU_dD02,
            FusedMulAddSIMD(
                vetU0, FusedMulAddSIMD(FDPart3tmp39, hDD_dupD220, FDPart3tmp230),
                FusedMulAddSIMD(FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp39, hDD_dupD221, FDPart3tmp245),
                                FusedMulAddSIMD(FDPart3tmp495,
                                                FusedMulAddSIMD(FDPart3_Rational_2_3, FDPart3tmp39, MulSIMD(FDPart3_Rational_2_3, FDPart3tmp43)),
                                                FusedMulAddSIMD(MulSIMD(FDPart3tmp21, FDPart3tmp4), MulSIMD(FDPart3tmp44, vetU_dD22),
                                                                FusedMulSubSIMD(FDPart3tmp200, vetU_dD12,
                                                                                MulSIMD(FDPart3tmp7, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp29),
                                                                                                             MulSIMD(aDD22, alpha)))))))))));
const REAL_SIMD_ARRAY __RHS_exp_17 = FusedMulAddSIMD(
    FDPart3tmp86, MulSIMD(MulSIMD(FDPart3tmp329, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
    FusedMulAddSIMD(
        FDPart3tmp78, MulSIMD(MulSIMD(FDPart3tmp331, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
        FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp78), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD1, alpha)),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp52, FDPart3tmp86), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD0, alpha)),
                FusedMulAddSIMD(
                    FDPart3tmp355, MulSIMD(FDPart3tmp411, FDPart3tmp505),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp52, FDPart3tmp68), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD2, alpha)),
                        FusedMulAddSIMD(
                            lambdaU_dupD00, vetU0,
                            FusedMulAddSIMD(
                                FDPart3tmp318, MulSIMD(FDPart3tmp4, lambdaU2),
                                FusedMulAddSIMD(
                                    FDPart3tmp81,
                                    AddSIMD(vetU_dDD011, NegFusedMulAddSIMD(FDPart3_Integer_2, vetU_dD11, SubSIMD(FDPart3tmp267, vetU0))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp87, vetU_dDD000,
                                        FusedMulAddSIMD(
                                            FDPart3tmp6, lambdaU_dupD02,
                                            FusedMulAddSIMD(
                                                FDPart3tmp72,
                                                NegFusedMulAddSIMD(
                                                    FDPart3_Integer_2, MulSIMD(f1_of_xx1, vetU_dD22),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp242, vetU_dD01,
                                                        FusedMulAddSIMD(FDPart3tmp316, vetU_dD00,
                                                                        FusedMulAddSIMD(f1_of_xx1,
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(f1_of_xx1__D1, vetU1)),
                                                                                        NegFusedMulAddSIMD(FDPart3tmp29, vetU0, vetU_dDD022))))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp501, FDPart3tmp84,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp502, FDPart3tmp87,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp415, FDPart3tmp509,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp499, FDPart3tmp82,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp400, FDPart3tmp507,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp413, FDPart3tmp508,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp381, FDPart3tmp504,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp385, FDPart3tmp506,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp365, FDPart3tmp512,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp369, FDPart3tmp511,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp312, lambdaU1,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp357, FDPart3tmp513,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp167, FDPart3tmp320,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp3, lambdaU_dupD01,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp154, FDPart3tmp314,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp162, FDPart3tmp325,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp68,
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp332,
                                                                                                                                FDPart3tmp52),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_16,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                            MulSIMD(PI, alpha))),
                                                                                                                    NegFusedMulAddSIMD(
                                                                                                                        lambdaU0, vetU_dD00,
                                                                                                                        FDPart3tmp439))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_18 = MulSIMD(
    f0_of_xx0,
    FusedMulAddSIMD(
        FDPart3tmp80, MulSIMD(MulSIMD(FDPart3tmp331, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
        FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp80), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD1, alpha)),
            FusedMulAddSIMD(
                FDPart3tmp78,
                MulSIMD(MulSIMD(FDPart3tmp329, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp52, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD2, alpha)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp52, FDPart3tmp78), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD0, alpha)),
                        FusedMulAddSIMD(
                            FDPart3tmp411, MulSIMD(FDPart3tmp458, FDPart3tmp505),
                            NegFusedMulAddSIMD(
                                FDPart3tmp296,
                                MulSIMD(lambdaU2, vetU_dD12),
                                FusedMulAddSIMD(
                                    FDPart3tmp87,
                                    FusedMulAddSIMD(FDPart3tmp2, vetU_dDD100,
                                                    FusedMulAddSIMD(FDPart3tmp450, vetU1,
                                                                    FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp21,
                                                                                    MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp8, vetU_dD10))))),
                                    FusedMulAddSIMD(
                                        vetU0, FusedMulSubSIMD(FDPart3tmp2, lambdaU_dupD10, MulSIMD(FDPart3tmp8, lambdaU1)),
                                        FusedMulAddSIMD(
                                            FDPart3tmp59,
                                            FusedMulAddSIMD(FDPart3tmp337, FDPart3tmp6,
                                                            FusedMulAddSIMD(FDPart3tmp6,
                                                                            NegFusedMulSubSIMD(f1_of_xx1, f1_of_xx1__DD11, FDPart3tmp337),
                                                                            FusedMulAddSIMD(FDPart3tmp271, f1_of_xx1, FDPart3tmp448))),
                                            FusedMulAddSIMD(
                                                FDPart3tmp72,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp29, FDPart3tmp3,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp2, vetU_dDD122,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp20, FDPart3tmp242,
                                                            FusedMulAddSIMD(FDPart3tmp330, FDPart3tmp337,
                                                                            FusedMulSubSIMD(FDPart3tmp11, FDPart3tmp316,
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp2),
                                                                                                    MulSIMD(f1_of_xx1__D1, vetU_dD22))))))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp501, FDPart3tmp81,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp502, FDPart3tmp84,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp469, FDPart3tmp513,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp499, FDPart3tmp59,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp466, FDPart3tmp506,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp468, FDPart3tmp508,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp463, FDPart3tmp511,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp465, FDPart3tmp509,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp461, FDPart3tmp507,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp462, FDPart3tmp512,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp415, FDPart3tmp504,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp448, FDPart3tmp59,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp363, lambdaU_dupD11,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp370, lambdaU_dupD12,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp167, FDPart3tmp444,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp346, lambdaU1,
                                                                                                                AddSIMD(
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp154,
                                                                                                                        AddSIMD(FDPart3tmp347,
                                                                                                                                FDPart3tmp442),
                                                                                                                        FDPart3tmp470),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp58,
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(FDPart3tmp332,
                                                                                                                                    FDPart3tmp52),
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_Integer_16,
                                                                                                                                    FDPart3_NegativeOne_),
                                                                                                                                MulSIMD(PI, alpha))),
                                                                                                                        NegFusedMulAddSIMD(
                                                                                                                            FDPart3tmp11, lambdaU0,
                                                                                                                            FDPart3tmp446)))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_19 = MulSIMD(
    FDPart3tmp109,
    FusedMulAddSIMD(
        FDPart3tmp58, MulSIMD(MulSIMD(FDPart3tmp331, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
        FusedMulAddSIMD(
            MulSIMD(FDPart3tmp52, FDPart3tmp71), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD2, alpha)),
            FusedMulAddSIMD(
                FDPart3tmp68,
                MulSIMD(MulSIMD(FDPart3tmp329, FDPart3tmp52), MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, alpha))),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3tmp52, FDPart3tmp58), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD1, alpha)),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3tmp52, FDPart3tmp68), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_4_3), MulSIMD(trK_dD0, alpha)),
                        NegFusedMulAddSIMD(
                            FDPart3tmp2, MulSIMD(FDPart3tmp274, lambdaU1),
                            FusedMulAddSIMD(
                                FDPart3tmp8, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp272), MulSIMD(lambdaU2, vetU_dD22)),
                                FusedMulAddSIMD(
                                    vetU0, FusedMulAddSIMD(FDPart3tmp348, lambdaU2, MulSIMD(FDPart3tmp5, lambdaU_dupD20)),
                                    FusedMulAddSIMD(
                                        FDPart3tmp411, MulSIMD(FDPart3tmp482, FDPart3tmp505),
                                        FusedMulAddSIMD(
                                            FDPart3tmp81,
                                            FusedMulAddSIMD(
                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp274, FDPart3tmp323),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp5, vetU_dDD211,
                                                    FusedMulAddSIMD(
                                                        vetU2, FusedMulAddSIMD(FDPart3tmp21, FDPart3tmp479, MulSIMD(FDPart3tmp341, f1_of_xx1__DD11)),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp16, f0_of_xx0,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp474, FDPart3tmp6,
                                                                FusedMulAddSIMD(FDPart3tmp2, MulSIMD(FDPart3tmp479, vetU2),
                                                                                NegFusedMulAddSIMD(FDPart3tmp272,
                                                                                                   MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp2),
                                                                                                           MulSIMD(f1_of_xx1__D1, vetU_dD21)),
                                                                                                   FDPart3tmp6))))))),
                                            FusedMulAddSIMD(
                                                FDPart3tmp87,
                                                FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp21,
                                                                FusedMulAddSIMD(FDPart3tmp5, vetU_dDD200,
                                                                                FusedMulSubSIMD(FDPart3tmp14, FDPart3tmp450,
                                                                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp4),
                                                                                                        MulSIMD(FDPart3tmp8, vetU_dD20))))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp501, FDPart3tmp59,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp502, FDPart3tmp82,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp491, FDPart3tmp512,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp499, FDPart3tmp72,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp489, FDPart3tmp508,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp490, FDPart3tmp511,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp485, FDPart3tmp507,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp488, FDPart3tmp506,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp468, FDPart3tmp509,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp483, FDPart3tmp513,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp375, lambdaU_dupD22,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp413, FDPart3tmp504,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp167, FDPart3tmp471,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp3,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp342, lambdaU2,
                                                                                                            MulSIMD(FDPart3tmp5, lambdaU_dupD21)),
                                                                                                        AddSIMD(
                                                                                                            FusedMulAddSIMD(FDPart3tmp154,
                                                                                                                            FDPart3tmp476,
                                                                                                                            FDPart3tmp492),
                                                                                                            AddSIMD(
                                                                                                                AddSIMD(FDPart3tmp475, FDPart3tmp478),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp71,
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3tmp332,
                                                                                                                                FDPart3tmp52),
                                                                                                                        MulSIMD(
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_16,
                                                                                                                                FDPart3_NegativeOne_),
                                                                                                                            MulSIMD(PI, alpha))),
                                                                                                                    NegFusedMulAddSIMD(
                                                                                                                        FDPart3tmp16, lambdaU0,
                                                                                                                        FDPart3tmp473))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_20 = NegFusedMulAddSIMD(
    FDPart3tmp169,
    MulSIMD(FDPart3tmp81,
            FusedMulAddSIMD(FDPart3tmp215, alpha_dD2,
                            FusedMulAddSIMD(FDPart3tmp217, alpha_dD1,
                                            FusedMulAddSIMD(FDPart3tmp219, alpha_dD0, NegFusedMulAddSIMD(FDPart3tmp185, alpha_dD1, alpha_dDD11))))),
    FusedMulAddSIMD(
        FDPart3tmp411, MulSIMD(aDD00, alpha),
        FusedMulAddSIMD(
            FDPart3tmp510,
            MulSIMD(PI, FusedMulAddSIMD(
                            FDPart3tmp152, MulSIMD(FDPart3tmp525, FusedMulSubSIMD(FDPart3tmp74, FDPart3tmp93, MulSIMD(FDPart3tmp76, FDPart3tmp93))),
                            FusedMulAddSIMD(
                                FDPart3tmp130, MulSIMD(FDPart3tmp524, NegFusedMulAddSIMD(FDPart3tmp40, FDPart3tmp93, FDPart3tmp164)),
                                FusedMulAddSIMD(
                                    FDPart3tmp139,
                                    MulSIMD(FDPart3tmp524, FusedMulSubSIMD(FDPart3tmp49, FDPart3tmp93, MulSIMD(FDPart3tmp34, FDPart3tmp93))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp161,
                                        MulSIMD(FDPart3tmp525, FusedMulAddSIMD(FDPart3tmp55, FDPart3tmp93, MulSIMD(FDPart3tmp57, FDPart3tmp93))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp166,
                                            MulSIMD(FDPart3tmp525, FusedMulSubSIMD(FDPart3tmp64, FDPart3tmp93, MulSIMD(FDPart3tmp66, FDPart3tmp93))),
                                            FusedMulAddSIMD(T4UU00, MulSIMD(alpha, alpha),
                                                            MulSIMD(FDPart3tmp118, MulSIMD(FDPart3tmp524, FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp93,
                                                                                                                          FDPart3tmp145)))))))))),
            FusedMulAddSIMD(
                FDPart3tmp281, MulSIMD(FDPart3tmp397, alpha),
                FusedMulAddSIMD(
                    FDPart3tmp393, MulSIMD(FDPart3tmp423, alpha),
                    FusedMulAddSIMD(
                        FDPart3tmp18, MulSIMD(FDPart3tmp429, alpha),
                        FusedMulAddSIMD(
                            FDPart3tmp269, MulSIMD(FDPart3tmp406, alpha),
                            FusedMulAddSIMD(
                                FDPart3_Rational_1_3, MulSIMD(alpha, MulSIMD(trK, trK)),
                                FusedMulAddSIMD(
                                    FDPart3tmp13, MulSIMD(FDPart3tmp431, alpha),
                                    FusedMulAddSIMD(
                                        FDPart3tmp6, trK_dupD2,
                                        FusedMulAddSIMD(
                                            trK_dupD0, vetU0,
                                            FusedMulAddSIMD(
                                                FDPart3tmp522, AddSIMD(FDPart3tmp239, FDPart3tmp521),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp522, AddSIMD(FDPart3tmp240, FDPart3tmp521),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp520, AddSIMD(FDPart3tmp224, FDPart3tmp519),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp520, AddSIMD(FDPart3tmp226, FDPart3tmp519),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp518, AddSIMD(FDPart3tmp187, FDPart3tmp517),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp518, AddSIMD(FDPart3tmp190, FDPart3tmp517),
                                                                    NegFusedMulAddSIMD(
                                                                        FDPart3tmp169,
                                                                        MulSIMD(FDPart3tmp87,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp177, alpha_dD2,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp179, alpha_dD1,
                                                                                        FusedMulAddSIMD(FDPart3tmp181, alpha_dD0,
                                                                                                        NegFusedMulAddSIMD(FDPart3tmp188, alpha_dD0,
                                                                                                                           alpha_dDD00))))),
                                                                        FusedMulSubSIMD(
                                                                            FDPart3tmp3, trK_dupD1,
                                                                            MulSIMD(FDPart3tmp169,
                                                                                    MulSIMD(FDPart3tmp72,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp258, alpha_dD2,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp260, alpha_dD1,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp262, alpha_dD0,
                                                                                                        NegFusedMulAddSIMD(
                                                                                                            FDPart3tmp222,
                                                                                                            alpha_dD2,
                                                                                                            alpha_dDD22)))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_21 = FusedMulAddSIMD(
    FDPart3tmp369, FDPart3tmp530,
    FusedMulAddSIMD(
        FDPart3tmp374, FDPart3tmp531,
        FusedMulAddSIMD(
            FDPart3tmp362, FDPart3tmp528,
            FusedMulAddSIMD(FDPart3tmp365, FDPart3tmp529,
                            FusedMulAddSIMD(FDPart3tmp355, FDPart3tmp526,
                                            FusedMulAddSIMD(FDPart3tmp357, FDPart3tmp527,
                                                            FusedMulAddSIMD(FDPart3tmp6, vetU_dupD02,
                                                                            FusedMulAddSIMD(vetU0, vetU_dupD00,
                                                                                            FusedMulAddSIMD(FDPart3tmp3, vetU_dupD01, betU0)))))))));
const REAL_SIMD_ARRAY __RHS_exp_22 = MulSIMD(
    f0_of_xx0,
    FusedMulAddSIMD(
        FDPart3tmp462, FDPart3tmp529,
        FusedMulAddSIMD(
            FDPart3tmp463, FDPart3tmp530,
            FusedMulAddSIMD(
                FDPart3tmp459, FDPart3tmp527,
                FusedMulAddSIMD(
                    FDPart3tmp461, FDPart3tmp528,
                    FusedMulAddSIMD(FDPart3tmp370, vetU_dupD12,
                                    FusedMulAddSIMD(FDPart3tmp458, FDPart3tmp526,
                                                    FusedMulAddSIMD(FDPart3tmp464, FDPart3tmp531,
                                                                    FusedMulAddSIMD(vetU0, FusedMulAddSIMD(FDPart3tmp2, vetU_dupD10, FDPart3tmp10),
                                                                                    FusedMulAddSIMD(FDPart3tmp2, betU1,
                                                                                                    MulSIMD(FDPart3tmp363, vetU_dupD11)))))))))));
const REAL_SIMD_ARRAY __RHS_exp_23 = MulSIMD(
    FDPart3tmp109,
    FusedMulAddSIMD(
        FDPart3tmp487, FDPart3tmp530,
        FusedMulAddSIMD(
            FDPart3tmp488, FDPart3tmp531,
            FusedMulAddSIMD(
                FDPart3tmp485, FDPart3tmp528,
                FusedMulAddSIMD(
                    FDPart3tmp486, FDPart3tmp529,
                    FusedMulAddSIMD(
                        FDPart3tmp482, FDPart3tmp526,
                        FusedMulAddSIMD(
                            FDPart3tmp483, FDPart3tmp527,
                            FusedMulAddSIMD(FDPart3tmp5, betU2,
                                            FusedMulAddSIMD(vetU0, FusedMulAddSIMD(FDPart3tmp5, vetU_dupD20, FDPart3tmp15),
                                                            FusedMulAddSIMD(FDPart3tmp3, FusedMulAddSIMD(FDPart3tmp5, vetU_dupD21, FDPart3tmp273),
                                                                            MulSIMD(FDPart3tmp375, vetU_dupD22)))))))))));

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
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
