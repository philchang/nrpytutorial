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
 * Set Ricci tensor.
 */
void Ricci_eval__rfm__Spherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                const rfm_struct *restrict rfmstruct, const REAL *restrict in_gfs, REAL *restrict auxevol_gfs) {
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
         * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
         * Read gridfunction(s) from main memory and compute FD stencils as needed.
         */
        const REAL_SIMD_ARRAY hDD00_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD00_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i1m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i1m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i1p1_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i1p2_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i1m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i1m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i1p1_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i1p2_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i1m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i1p1_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i1p2_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i1m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i1m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i1p1_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i1p2_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2 + 2)]);
        const REAL_SIMD_ARRAY lambdaU0_i2m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY lambdaU0_i2m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU0_i2p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY lambdaU0_i2p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY lambdaU1_i2m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY lambdaU1_i2m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU1_i2p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY lambdaU1_i2p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY lambdaU2_i2m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY lambdaU2_i2m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0m1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i0p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i1p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY lambdaU2_i2p1 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY lambdaU2_i2p2 = ReadSIMD(&in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2 + 2)]);
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
        const REAL_SIMD_ARRAY hDD_dDD0000 = SIMD_fd_function_dDD00_fdorder4(hDD00, hDD00_i0m1, hDD00_i0m2, hDD00_i0p1, hDD00_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD0001 =
            SIMD_fd_function_dDD01_fdorder4(hDD00_i0m1_i1m1, hDD00_i0m1_i1m2, hDD00_i0m1_i1p1, hDD00_i0m1_i1p2, hDD00_i0m2_i1m1, hDD00_i0m2_i1m2,
                                            hDD00_i0m2_i1p1, hDD00_i0m2_i1p2, hDD00_i0p1_i1m1, hDD00_i0p1_i1m2, hDD00_i0p1_i1p1, hDD00_i0p1_i1p2,
                                            hDD00_i0p2_i1m1, hDD00_i0p2_i1m2, hDD00_i0p2_i1p1, hDD00_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0002 =
            SIMD_fd_function_dDD02_fdorder4(hDD00_i0m1_i2m1, hDD00_i0m1_i2m2, hDD00_i0m1_i2p1, hDD00_i0m1_i2p2, hDD00_i0m2_i2m1, hDD00_i0m2_i2m2,
                                            hDD00_i0m2_i2p1, hDD00_i0m2_i2p2, hDD00_i0p1_i2m1, hDD00_i0p1_i2m2, hDD00_i0p1_i2p1, hDD00_i0p1_i2p2,
                                            hDD00_i0p2_i2m1, hDD00_i0p2_i2m2, hDD00_i0p2_i2p1, hDD00_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD0011 = SIMD_fd_function_dDD11_fdorder4(hDD00, hDD00_i1m1, hDD00_i1m2, hDD00_i1p1, hDD00_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0012 =
            SIMD_fd_function_dDD12_fdorder4(hDD00_i1m1_i2m1, hDD00_i1m1_i2m2, hDD00_i1m1_i2p1, hDD00_i1m1_i2p2, hDD00_i1m2_i2m1, hDD00_i1m2_i2m2,
                                            hDD00_i1m2_i2p1, hDD00_i1m2_i2p2, hDD00_i1p1_i2m1, hDD00_i1p1_i2m2, hDD00_i1p1_i2p1, hDD00_i1p1_i2p2,
                                            hDD00_i1p2_i2m1, hDD00_i1p2_i2m2, hDD00_i1p2_i2p1, hDD00_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD0022 = SIMD_fd_function_dDD22_fdorder4(hDD00, hDD00_i2m1, hDD00_i2m2, hDD00_i2p1, hDD00_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD0100 = SIMD_fd_function_dDD00_fdorder4(hDD01, hDD01_i0m1, hDD01_i0m2, hDD01_i0p1, hDD01_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD0101 =
            SIMD_fd_function_dDD01_fdorder4(hDD01_i0m1_i1m1, hDD01_i0m1_i1m2, hDD01_i0m1_i1p1, hDD01_i0m1_i1p2, hDD01_i0m2_i1m1, hDD01_i0m2_i1m2,
                                            hDD01_i0m2_i1p1, hDD01_i0m2_i1p2, hDD01_i0p1_i1m1, hDD01_i0p1_i1m2, hDD01_i0p1_i1p1, hDD01_i0p1_i1p2,
                                            hDD01_i0p2_i1m1, hDD01_i0p2_i1m2, hDD01_i0p2_i1p1, hDD01_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0102 =
            SIMD_fd_function_dDD02_fdorder4(hDD01_i0m1_i2m1, hDD01_i0m1_i2m2, hDD01_i0m1_i2p1, hDD01_i0m1_i2p2, hDD01_i0m2_i2m1, hDD01_i0m2_i2m2,
                                            hDD01_i0m2_i2p1, hDD01_i0m2_i2p2, hDD01_i0p1_i2m1, hDD01_i0p1_i2m2, hDD01_i0p1_i2p1, hDD01_i0p1_i2p2,
                                            hDD01_i0p2_i2m1, hDD01_i0p2_i2m2, hDD01_i0p2_i2p1, hDD01_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD0111 = SIMD_fd_function_dDD11_fdorder4(hDD01, hDD01_i1m1, hDD01_i1m2, hDD01_i1p1, hDD01_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0112 =
            SIMD_fd_function_dDD12_fdorder4(hDD01_i1m1_i2m1, hDD01_i1m1_i2m2, hDD01_i1m1_i2p1, hDD01_i1m1_i2p2, hDD01_i1m2_i2m1, hDD01_i1m2_i2m2,
                                            hDD01_i1m2_i2p1, hDD01_i1m2_i2p2, hDD01_i1p1_i2m1, hDD01_i1p1_i2m2, hDD01_i1p1_i2p1, hDD01_i1p1_i2p2,
                                            hDD01_i1p2_i2m1, hDD01_i1p2_i2m2, hDD01_i1p2_i2p1, hDD01_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD0122 = SIMD_fd_function_dDD22_fdorder4(hDD01, hDD01_i2m1, hDD01_i2m2, hDD01_i2p1, hDD01_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD0200 = SIMD_fd_function_dDD00_fdorder4(hDD02, hDD02_i0m1, hDD02_i0m2, hDD02_i0p1, hDD02_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD0201 =
            SIMD_fd_function_dDD01_fdorder4(hDD02_i0m1_i1m1, hDD02_i0m1_i1m2, hDD02_i0m1_i1p1, hDD02_i0m1_i1p2, hDD02_i0m2_i1m1, hDD02_i0m2_i1m2,
                                            hDD02_i0m2_i1p1, hDD02_i0m2_i1p2, hDD02_i0p1_i1m1, hDD02_i0p1_i1m2, hDD02_i0p1_i1p1, hDD02_i0p1_i1p2,
                                            hDD02_i0p2_i1m1, hDD02_i0p2_i1m2, hDD02_i0p2_i1p1, hDD02_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0202 =
            SIMD_fd_function_dDD02_fdorder4(hDD02_i0m1_i2m1, hDD02_i0m1_i2m2, hDD02_i0m1_i2p1, hDD02_i0m1_i2p2, hDD02_i0m2_i2m1, hDD02_i0m2_i2m2,
                                            hDD02_i0m2_i2p1, hDD02_i0m2_i2p2, hDD02_i0p1_i2m1, hDD02_i0p1_i2m2, hDD02_i0p1_i2p1, hDD02_i0p1_i2p2,
                                            hDD02_i0p2_i2m1, hDD02_i0p2_i2m2, hDD02_i0p2_i2p1, hDD02_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD0211 = SIMD_fd_function_dDD11_fdorder4(hDD02, hDD02_i1m1, hDD02_i1m2, hDD02_i1p1, hDD02_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD0212 =
            SIMD_fd_function_dDD12_fdorder4(hDD02_i1m1_i2m1, hDD02_i1m1_i2m2, hDD02_i1m1_i2p1, hDD02_i1m1_i2p2, hDD02_i1m2_i2m1, hDD02_i1m2_i2m2,
                                            hDD02_i1m2_i2p1, hDD02_i1m2_i2p2, hDD02_i1p1_i2m1, hDD02_i1p1_i2m2, hDD02_i1p1_i2p1, hDD02_i1p1_i2p2,
                                            hDD02_i1p2_i2m1, hDD02_i1p2_i2m2, hDD02_i1p2_i2p1, hDD02_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD0222 = SIMD_fd_function_dDD22_fdorder4(hDD02, hDD02_i2m1, hDD02_i2m2, hDD02_i2p1, hDD02_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD1100 = SIMD_fd_function_dDD00_fdorder4(hDD11, hDD11_i0m1, hDD11_i0m2, hDD11_i0p1, hDD11_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD1101 =
            SIMD_fd_function_dDD01_fdorder4(hDD11_i0m1_i1m1, hDD11_i0m1_i1m2, hDD11_i0m1_i1p1, hDD11_i0m1_i1p2, hDD11_i0m2_i1m1, hDD11_i0m2_i1m2,
                                            hDD11_i0m2_i1p1, hDD11_i0m2_i1p2, hDD11_i0p1_i1m1, hDD11_i0p1_i1m2, hDD11_i0p1_i1p1, hDD11_i0p1_i1p2,
                                            hDD11_i0p2_i1m1, hDD11_i0p2_i1m2, hDD11_i0p2_i1p1, hDD11_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD1102 =
            SIMD_fd_function_dDD02_fdorder4(hDD11_i0m1_i2m1, hDD11_i0m1_i2m2, hDD11_i0m1_i2p1, hDD11_i0m1_i2p2, hDD11_i0m2_i2m1, hDD11_i0m2_i2m2,
                                            hDD11_i0m2_i2p1, hDD11_i0m2_i2p2, hDD11_i0p1_i2m1, hDD11_i0p1_i2m2, hDD11_i0p1_i2p1, hDD11_i0p1_i2p2,
                                            hDD11_i0p2_i2m1, hDD11_i0p2_i2m2, hDD11_i0p2_i2p1, hDD11_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD1111 = SIMD_fd_function_dDD11_fdorder4(hDD11, hDD11_i1m1, hDD11_i1m2, hDD11_i1p1, hDD11_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD1112 =
            SIMD_fd_function_dDD12_fdorder4(hDD11_i1m1_i2m1, hDD11_i1m1_i2m2, hDD11_i1m1_i2p1, hDD11_i1m1_i2p2, hDD11_i1m2_i2m1, hDD11_i1m2_i2m2,
                                            hDD11_i1m2_i2p1, hDD11_i1m2_i2p2, hDD11_i1p1_i2m1, hDD11_i1p1_i2m2, hDD11_i1p1_i2p1, hDD11_i1p1_i2p2,
                                            hDD11_i1p2_i2m1, hDD11_i1p2_i2m2, hDD11_i1p2_i2p1, hDD11_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD1122 = SIMD_fd_function_dDD22_fdorder4(hDD11, hDD11_i2m1, hDD11_i2m2, hDD11_i2p1, hDD11_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD1200 = SIMD_fd_function_dDD00_fdorder4(hDD12, hDD12_i0m1, hDD12_i0m2, hDD12_i0p1, hDD12_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD1201 =
            SIMD_fd_function_dDD01_fdorder4(hDD12_i0m1_i1m1, hDD12_i0m1_i1m2, hDD12_i0m1_i1p1, hDD12_i0m1_i1p2, hDD12_i0m2_i1m1, hDD12_i0m2_i1m2,
                                            hDD12_i0m2_i1p1, hDD12_i0m2_i1p2, hDD12_i0p1_i1m1, hDD12_i0p1_i1m2, hDD12_i0p1_i1p1, hDD12_i0p1_i1p2,
                                            hDD12_i0p2_i1m1, hDD12_i0p2_i1m2, hDD12_i0p2_i1p1, hDD12_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD1202 =
            SIMD_fd_function_dDD02_fdorder4(hDD12_i0m1_i2m1, hDD12_i0m1_i2m2, hDD12_i0m1_i2p1, hDD12_i0m1_i2p2, hDD12_i0m2_i2m1, hDD12_i0m2_i2m2,
                                            hDD12_i0m2_i2p1, hDD12_i0m2_i2p2, hDD12_i0p1_i2m1, hDD12_i0p1_i2m2, hDD12_i0p1_i2p1, hDD12_i0p1_i2p2,
                                            hDD12_i0p2_i2m1, hDD12_i0p2_i2m2, hDD12_i0p2_i2p1, hDD12_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD1211 = SIMD_fd_function_dDD11_fdorder4(hDD12, hDD12_i1m1, hDD12_i1m2, hDD12_i1p1, hDD12_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD1212 =
            SIMD_fd_function_dDD12_fdorder4(hDD12_i1m1_i2m1, hDD12_i1m1_i2m2, hDD12_i1m1_i2p1, hDD12_i1m1_i2p2, hDD12_i1m2_i2m1, hDD12_i1m2_i2m2,
                                            hDD12_i1m2_i2p1, hDD12_i1m2_i2p2, hDD12_i1p1_i2m1, hDD12_i1p1_i2m2, hDD12_i1p1_i2p1, hDD12_i1p1_i2p2,
                                            hDD12_i1p2_i2m1, hDD12_i1p2_i2m2, hDD12_i1p2_i2p1, hDD12_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD1222 = SIMD_fd_function_dDD22_fdorder4(hDD12, hDD12_i2m1, hDD12_i2m2, hDD12_i2p1, hDD12_i2p2, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD2200 = SIMD_fd_function_dDD00_fdorder4(hDD22, hDD22_i0m1, hDD22_i0m2, hDD22_i0p1, hDD22_i0p2, invdxx0);
        const REAL_SIMD_ARRAY hDD_dDD2201 =
            SIMD_fd_function_dDD01_fdorder4(hDD22_i0m1_i1m1, hDD22_i0m1_i1m2, hDD22_i0m1_i1p1, hDD22_i0m1_i1p2, hDD22_i0m2_i1m1, hDD22_i0m2_i1m2,
                                            hDD22_i0m2_i1p1, hDD22_i0m2_i1p2, hDD22_i0p1_i1m1, hDD22_i0p1_i1m2, hDD22_i0p1_i1p1, hDD22_i0p1_i1p2,
                                            hDD22_i0p2_i1m1, hDD22_i0p2_i1m2, hDD22_i0p2_i1p1, hDD22_i0p2_i1p2, invdxx0, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD2202 =
            SIMD_fd_function_dDD02_fdorder4(hDD22_i0m1_i2m1, hDD22_i0m1_i2m2, hDD22_i0m1_i2p1, hDD22_i0m1_i2p2, hDD22_i0m2_i2m1, hDD22_i0m2_i2m2,
                                            hDD22_i0m2_i2p1, hDD22_i0m2_i2p2, hDD22_i0p1_i2m1, hDD22_i0p1_i2m2, hDD22_i0p1_i2p1, hDD22_i0p1_i2p2,
                                            hDD22_i0p2_i2m1, hDD22_i0p2_i2m2, hDD22_i0p2_i2p1, hDD22_i0p2_i2p2, invdxx0, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD2211 = SIMD_fd_function_dDD11_fdorder4(hDD22, hDD22_i1m1, hDD22_i1m2, hDD22_i1p1, hDD22_i1p2, invdxx1);
        const REAL_SIMD_ARRAY hDD_dDD2212 =
            SIMD_fd_function_dDD12_fdorder4(hDD22_i1m1_i2m1, hDD22_i1m1_i2m2, hDD22_i1m1_i2p1, hDD22_i1m1_i2p2, hDD22_i1m2_i2m1, hDD22_i1m2_i2m2,
                                            hDD22_i1m2_i2p1, hDD22_i1m2_i2p2, hDD22_i1p1_i2m1, hDD22_i1p1_i2m2, hDD22_i1p1_i2p1, hDD22_i1p1_i2p2,
                                            hDD22_i1p2_i2m1, hDD22_i1p2_i2m2, hDD22_i1p2_i2p1, hDD22_i1p2_i2p2, invdxx1, invdxx2);
        const REAL_SIMD_ARRAY hDD_dDD2222 = SIMD_fd_function_dDD22_fdorder4(hDD22, hDD22_i2m1, hDD22_i2m2, hDD22_i2p1, hDD22_i2p2, invdxx2);
        const REAL_SIMD_ARRAY lambdaU_dD00 = SIMD_fd_function_dD0_fdorder4(lambdaU0_i0m1, lambdaU0_i0m2, lambdaU0_i0p1, lambdaU0_i0p2, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dD01 = SIMD_fd_function_dD1_fdorder4(lambdaU0_i1m1, lambdaU0_i1m2, lambdaU0_i1p1, lambdaU0_i1p2, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dD02 = SIMD_fd_function_dD2_fdorder4(lambdaU0_i2m1, lambdaU0_i2m2, lambdaU0_i2p1, lambdaU0_i2p2, invdxx2);
        const REAL_SIMD_ARRAY lambdaU_dD10 = SIMD_fd_function_dD0_fdorder4(lambdaU1_i0m1, lambdaU1_i0m2, lambdaU1_i0p1, lambdaU1_i0p2, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dD11 = SIMD_fd_function_dD1_fdorder4(lambdaU1_i1m1, lambdaU1_i1m2, lambdaU1_i1p1, lambdaU1_i1p2, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dD12 = SIMD_fd_function_dD2_fdorder4(lambdaU1_i2m1, lambdaU1_i2m2, lambdaU1_i2p1, lambdaU1_i2p2, invdxx2);
        const REAL_SIMD_ARRAY lambdaU_dD20 = SIMD_fd_function_dD0_fdorder4(lambdaU2_i0m1, lambdaU2_i0m2, lambdaU2_i0p1, lambdaU2_i0p2, invdxx0);
        const REAL_SIMD_ARRAY lambdaU_dD21 = SIMD_fd_function_dD1_fdorder4(lambdaU2_i1m1, lambdaU2_i1m2, lambdaU2_i1p1, lambdaU2_i1p2, invdxx1);
        const REAL_SIMD_ARRAY lambdaU_dD22 = SIMD_fd_function_dD2_fdorder4(lambdaU2_i2m1, lambdaU2_i2m2, lambdaU2_i2p1, lambdaU2_i2p2, invdxx2);

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const double dblFDPart3_Integer_1 = 1.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        const double dblFDPart3_Integer_3 = 3.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_3 = ConstSIMD(dblFDPart3_Integer_3);

        const double dblFDPart3_Integer_4 = 4.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_4 = ConstSIMD(dblFDPart3_Integer_4);

        const double dblFDPart3_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        const REAL_SIMD_ARRAY FDPart3tmp0 = AddSIMD(FDPart3_Integer_1, hDD00);
        const REAL_SIMD_ARRAY FDPart3tmp2 = MulSIMD(f1_of_xx1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp5 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp13 = MulSIMD(FDPart3_Integer_2, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp21 = DivSIMD(FDPart3_Integer_1, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp24 = MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp25 = MulSIMD(f1_of_xx1, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp26 = MulSIMD(f1_of_xx1, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp32 = MulSIMD(f0_of_xx0, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp43 = DivSIMD(FDPart3_Integer_1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp47 = MulSIMD(f1_of_xx1, hDD_dD020);
        const REAL_SIMD_ARRAY FDPart3tmp57 = MulSIMD(f1_of_xx1__D1, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp60 = MulSIMD(f0_of_xx0, hDD_dD012);
        const REAL_SIMD_ARRAY FDPart3tmp65 = MulSIMD(f1_of_xx1, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp66 = MulSIMD(f1_of_xx1, hDD_dD022);
        const REAL_SIMD_ARRAY FDPart3tmp73 = MulSIMD(FDPart3_Integer_2, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp76 = MulSIMD(FDPart3_Rational_1_2, hDD_dD000);
        const REAL_SIMD_ARRAY FDPart3tmp147 = MulSIMD(f1_of_xx1, hDD_dD120);
        const REAL_SIMD_ARRAY FDPart3tmp207 = SubSIMD(lambdaU_dD01, lambdaU1);
        const REAL_SIMD_ARRAY FDPart3tmp236 = MulSIMD(f1_of_xx1__D1, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp284 = NegFusedMulAddSIMD(f1_of_xx1, lambdaU2, lambdaU_dD02);
        const REAL_SIMD_ARRAY FDPart3tmp289 = MulSIMD(f0_of_xx0, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp296 = MulSIMD(MulSIMD(f1_of_xx1, f1_of_xx1), f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp305 = MulSIMD(f1_of_xx1, f1_of_xx1__DD11);
        const REAL_SIMD_ARRAY FDPart3tmp3 = MulSIMD(FDPart3tmp2, MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp6 = MulSIMD(FDPart3tmp5, hDD11);
        const REAL_SIMD_ARRAY FDPart3tmp8 = MulSIMD(FDPart3tmp2, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp20 = NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f1_of_xx1, hDD02), hDD_dD002);
        const REAL_SIMD_ARRAY FDPart3tmp28 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp26, f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp30 = NegFusedMulAddSIMD(FDPart3_Integer_2, hDD01, hDD_dD001);
        const REAL_SIMD_ARRAY FDPart3tmp39 = MulSIMD(f0_of_xx0, MulSIMD(f1_of_xx1, hDD_dD021));
        const REAL_SIMD_ARRAY FDPart3tmp40 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp25, f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp45 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp43, f1_of_xx1__D1));
        const REAL_SIMD_ARRAY FDPart3tmp48 = FusedMulAddSIMD(FDPart3tmp47, f0_of_xx0, FDPart3tmp26);
        const REAL_SIMD_ARRAY FDPart3tmp52 = FusedMulAddSIMD(f0_of_xx0, hDD_dD011, FusedMulSubSIMD(f0_of_xx0, hDD00, MulSIMD(f0_of_xx0, hDD11)));
        const REAL_SIMD_ARRAY FDPart3tmp56 = FusedMulAddSIMD(f0_of_xx0, hDD_dD010, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp67 = MulSIMD(FDPart3tmp2, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp83 = MulSIMD(FDPart3tmp26, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp90 = MulSIMD(FDPart3tmp5, hDD_dD110);
        const REAL_SIMD_ARRAY FDPart3tmp94 = MulSIMD(FDPart3tmp5, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp101 = MulSIMD(FDPart3tmp2, FDPart3tmp73);
        const REAL_SIMD_ARRAY FDPart3tmp102 = MulSIMD(FDPart3tmp2, hDD22);
        const REAL_SIMD_ARRAY FDPart3tmp130 = MulSIMD(FDPart3tmp25, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp141 = MulSIMD(FDPart3tmp5, hDD_dD112);
        const REAL_SIMD_ARRAY FDPart3tmp142 = MulSIMD(FDPart3_Integer_2, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp156 = MulSIMD(FDPart3tmp5, hDD_dD111);
        const REAL_SIMD_ARRAY FDPart3tmp195 = MulSIMD(FDPart3tmp26, FDPart3tmp73);
        const REAL_SIMD_ARRAY FDPart3tmp196 = MulSIMD(FDPart3tmp13, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp204 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp32);
        const REAL_SIMD_ARRAY FDPart3tmp213 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp147, f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp214 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(f0_of_xx0, f1_of_xx1__D1));
        const REAL_SIMD_ARRAY FDPart3tmp235 = DivSIMD(FDPart3_Integer_1, FDPart3tmp2);
        const REAL_SIMD_ARRAY FDPart3tmp286 = MulSIMD(FDPart3tmp21, FDPart3tmp43);
        const REAL_SIMD_ARRAY FDPart3tmp291 = FusedMulAddSIMD(FDPart3tmp289, hDD_dDD0201, FusedMulAddSIMD(f1_of_xx1, hDD_dD021, FDPart3tmp57));
        const REAL_SIMD_ARRAY FDPart3tmp297 = MulSIMD(FDPart3tmp296, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp377 = MulSIMD(FDPart3_Integer_2, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp4 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp3, MulSIMD(hDD12, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp7 = AddSIMD(FDPart3tmp5, FDPart3tmp6);
        const REAL_SIMD_ARRAY FDPart3tmp9 = MulSIMD(FDPart3tmp8, hDD22);
        const REAL_SIMD_ARRAY FDPart3tmp15 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp8, MulSIMD(hDD02, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp17 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp5, MulSIMD(hDD01, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp23 = NegFusedMulAddSIMD(FDPart3tmp20, FDPart3tmp21, hDD_dDD0002);
        const REAL_SIMD_ARRAY FDPart3tmp31 = NegFusedMulAddSIMD(FDPart3tmp21, FDPart3tmp30, hDD_dDD0001);
        const REAL_SIMD_ARRAY FDPart3tmp36 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp25, FDPart3tmp5));
        const REAL_SIMD_ARRAY FDPart3tmp41 = AddSIMD(FDPart3tmp39, FDPart3tmp40);
        const REAL_SIMD_ARRAY FDPart3tmp46 = FusedMulAddSIMD(FDPart3tmp20, FDPart3tmp45, hDD_dDD0012);
        const REAL_SIMD_ARRAY FDPart3tmp53 = FusedMulSubSIMD(f0_of_xx0, hDD_dD000, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp52)));
        const REAL_SIMD_ARRAY FDPart3tmp59 = FusedMulAddSIMD(FDPart3tmp57, f0_of_xx0, FDPart3tmp39);
        const REAL_SIMD_ARRAY FDPart3tmp61 = NegFusedMulAddSIMD(FDPart3tmp57, f0_of_xx0, FDPart3tmp60);
        const REAL_SIMD_ARRAY FDPart3tmp91 = FusedMulAddSIMD(FDPart3tmp73, hDD11, FDPart3tmp90);
        const REAL_SIMD_ARRAY FDPart3tmp95 = FusedMulAddSIMD(FDPart3tmp25, FDPart3tmp73, MulSIMD(FDPart3tmp94, hDD_dD120));
        const REAL_SIMD_ARRAY FDPart3tmp143 = MulSIMD(FDPart3tmp142, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp157 = MulSIMD(FDPart3tmp94, hDD_dD121);
        const REAL_SIMD_ARRAY FDPart3tmp159 = MulSIMD(FDPart3tmp5, MulSIMD(f1_of_xx1__D1, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp161 = NegFusedMulAddSIMD(
            FDPart3_Integer_2, f0_of_xx0,
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, hDD11), FusedMulSubSIMD(FDPart3tmp73, hDD_dD011, FDPart3tmp90)));
        const REAL_SIMD_ARRAY FDPart3tmp181 = MulSIMD(FDPart3tmp8, hDD_dD222);
        const REAL_SIMD_ARRAY FDPart3tmp210 = FusedMulAddSIMD(FDPart3tmp21, lambdaU0, MulSIMD(FDPart3tmp21, lambdaU_dD11));
        const REAL_SIMD_ARRAY FDPart3tmp218 = FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp5, FDPart3tmp156);
        const REAL_SIMD_ARRAY FDPart3tmp225 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp5, FDPart3tmp65));
        const REAL_SIMD_ARRAY FDPart3tmp231 =
            NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp5), MulSIMD(f1_of_xx1__D1, hDD12), FDPart3tmp141);
        const REAL_SIMD_ARRAY FDPart3tmp238 = FusedMulSubSIMD(FDPart3tmp43, f1_of_xx1__DD11, MulSIMD(FDPart3tmp235, FDPart3tmp236));
        const REAL_SIMD_ARRAY FDPart3tmp283 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp83);
        const REAL_SIMD_ARRAY FDPart3tmp285 = FusedMulSubSIMD(FDPart3tmp21, lambdaU_dD12, MulSIMD(FDPart3tmp21, MulSIMD(f1_of_xx1__D1, lambdaU2)));
        const REAL_SIMD_ARRAY FDPart3tmp288 = FusedMulAddSIMD(
            FDPart3tmp286, lambdaU_dD22, FusedMulAddSIMD(FDPart3tmp286, MulSIMD(f1_of_xx1__D1, lambdaU1), MulSIMD(FDPart3tmp21, lambdaU0)));
        const REAL_SIMD_ARRAY FDPart3tmp300 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp8, f1_of_xx1__D1));
        const REAL_SIMD_ARRAY FDPart3tmp306 = MulSIMD(FDPart3_NegativeOne_, AddSIMD(FDPart3tmp305, FDPart3tmp236));
        const REAL_SIMD_ARRAY FDPart3tmp335 = FusedMulSubSIMD(FDPart3tmp235, FDPart3tmp236, MulSIMD(FDPart3tmp43, f1_of_xx1__DD11));
        const REAL_SIMD_ARRAY FDPart3tmp353 = MulSIMD(FDPart3tmp142, FDPart3tmp25);
        const REAL_SIMD_ARRAY FDPart3tmp356 =
            FusedMulAddSIMD(FDPart3tmp73, MulSIMD(f1_of_xx1, hDD_dD121),
                            FusedMulAddSIMD(FDPart3tmp73, MulSIMD(f1_of_xx1__D1, hDD12), MulSIMD(FDPart3tmp94, hDD_dDD1201)));
        const REAL_SIMD_ARRAY FDPart3tmp374 =
            FusedMulAddSIMD(FDPart3tmp8, hDD_dDD2201, MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3tmp65), MulSIMD(f0_of_xx0, hDD22)));
        const REAL_SIMD_ARRAY FDPart3tmp10 = AddSIMD(FDPart3tmp8, FDPart3tmp9);
        const REAL_SIMD_ARRAY FDPart3tmp29 = FusedMulAddSIMD(FDPart3tmp28, FDPart3tmp7, MulSIMD(FDPart3tmp24, MulSIMD(FDPart3tmp25, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp37 = FusedMulAddSIMD(FDPart3tmp0, FDPart3tmp36, MulSIMD(FDPart3tmp26, MulSIMD(FDPart3tmp5, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp42 = MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp21, FDPart3tmp41));
        const REAL_SIMD_ARRAY FDPart3tmp62 = AddSIMD(FDPart3tmp40, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp64 = FusedMulAddSIMD(FDPart3tmp0, FDPart3tmp7, FDPart3tmp17);
        const REAL_SIMD_ARRAY FDPart3tmp69 = FusedMulAddSIMD(
            FDPart3tmp66, f0_of_xx0, FusedMulAddSIMD(FDPart3tmp67, hDD00, FusedMulSubSIMD(FDPart3tmp32, FDPart3tmp65, MulSIMD(FDPart3tmp67, hDD22))));
        const REAL_SIMD_ARRAY FDPart3tmp92 = AddSIMD(FDPart3tmp73, FDPart3tmp91);
        const REAL_SIMD_ARRAY FDPart3tmp96 = AddSIMD(FDPart3tmp95, SubSIMD(FDPart3tmp59, FDPart3tmp60));
        const REAL_SIMD_ARRAY FDPart3tmp104 = FusedMulAddSIMD(FDPart3tmp102, FDPart3tmp73, MulSIMD(FDPart3tmp8, hDD_dD220));
        const REAL_SIMD_ARRAY FDPart3tmp106 = AddSIMD(FDPart3tmp95, SubSIMD(FDPart3tmp61, FDPart3tmp39));
        const REAL_SIMD_ARRAY FDPart3tmp145 = FusedMulAddSIMD(FDPart3tmp143, hDD22, MulSIMD(FDPart3tmp8, hDD_dD221));
        const REAL_SIMD_ARRAY FDPart3tmp148 =
            AddSIMD(FDPart3tmp60, FusedMulAddSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, hDD12)),
                                                  NegFusedMulAddSIMD(FDPart3tmp147, FDPart3tmp5, FDPart3tmp59)));
        const REAL_SIMD_ARRAY FDPart3tmp160 =
            FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp159, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp157, FDPart3tmp141));
        const REAL_SIMD_ARRAY FDPart3tmp182 =
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, f0_of_xx0),
                               FusedMulAddSIMD(FDPart3tmp2, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, hDD22)),
                                               FusedMulSubSIMD(FDPart3tmp66, FDPart3tmp73, MulSIMD(FDPart3tmp8, hDD_dD220))));
        const REAL_SIMD_ARRAY FDPart3tmp184 = FusedMulAddSIMD(
            FDPart3tmp5, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, f1_of_xx1__D1)),
            FusedMulAddSIMD(MulSIMD(FDPart3tmp5, f1_of_xx1), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, hDD22)),
                            FusedMulSubSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp94, hDD_dD122), MulSIMD(FDPart3tmp8, hDD_dD221))));
        const REAL_SIMD_ARRAY FDPart3tmp219 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp21, FDPart3tmp218));
        const REAL_SIMD_ARRAY FDPart3tmp220 =
            FusedMulAddSIMD(FDPart3tmp45, FDPart3tmp48, NegFusedMulAddSIMD(FDPart3tmp21, FDPart3tmp95, FDPart3tmp25));
        const REAL_SIMD_ARRAY FDPart3tmp222 = FusedMulAddSIMD(FDPart3tmp26, FDPart3tmp5, FDPart3tmp157);
        const REAL_SIMD_ARRAY FDPart3tmp226 = FusedMulAddSIMD(
            FDPart3tmp8, hDD01, FusedMulAddSIMD(FDPart3tmp94, hDD_dD122, FusedMulAddSIMD(FDPart3tmp225, hDD22, MulSIMD(FDPart3tmp6, FDPart3tmp65))));
        const REAL_SIMD_ARRAY FDPart3tmp232 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp21, FDPart3tmp231));
        const REAL_SIMD_ARRAY FDPart3tmp233 = AddSIMD(FDPart3tmp157, FDPart3tmp159);
        const REAL_SIMD_ARRAY FDPart3tmp293 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp67, hDD_dD221));
        const REAL_SIMD_ARRAY FDPart3tmp301 =
            FusedMulAddSIMD(FDPart3tmp300, hDD12, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp297, hDD02), FDPart3tmp181));
        const REAL_SIMD_ARRAY FDPart3tmp332 = FusedMulAddSIMD(FDPart3tmp231, FDPart3tmp45, MulSIMD(FDPart3tmp5, hDD_dDD1112));
        const REAL_SIMD_ARRAY FDPart3tmp34 =
            FusedMulSubSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp24), MulSIMD(hDD02, hDD12), MulSIMD(FDPart3tmp10, FDPart3tmp32));
        const REAL_SIMD_ARRAY FDPart3tmp54 = FusedMulAddSIMD(FDPart3tmp0, FDPart3tmp10, FDPart3tmp15);
        const REAL_SIMD_ARRAY FDPart3tmp63 = MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp21, FDPart3tmp62));
        const REAL_SIMD_ARRAY FDPart3tmp71 =
            FusedMulSubSIMD(FDPart3tmp67, hDD_dD000, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp69)));
        const REAL_SIMD_ARRAY FDPart3tmp105 = AddSIMD(FDPart3tmp101, FDPart3tmp104);
        const REAL_SIMD_ARRAY FDPart3tmp146 = AddSIMD(FDPart3tmp143, FDPart3tmp145);
        const REAL_SIMD_ARRAY FDPart3tmp221 = FusedMulAddSIMD(FDPart3tmp45, FDPart3tmp62, MulSIMD(f0_of_xx0, hDD_dDD0112));
        const REAL_SIMD_ARRAY FDPart3tmp223 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp21, FDPart3tmp222));
        const REAL_SIMD_ARRAY FDPart3tmp227 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp21, FDPart3tmp226));
        const REAL_SIMD_ARRAY FDPart3tmp230 = FusedMulAddSIMD(
            FDPart3tmp52, FDPart3tmp65, FusedMulSubSIMD(FDPart3tmp45, FDPart3tmp69, MulSIMD(f0_of_xx0, MulSIMD(f1_of_xx1, hDD_dD122))));
        const REAL_SIMD_ARRAY FDPart3tmp239 =
            FusedMulAddSIMD(FDPart3tmp45, FDPart3tmp59, FusedMulSubSIMD(FDPart3tmp238, FDPart3tmp28, MulSIMD(FDPart3tmp21, FDPart3tmp233)));
        const REAL_SIMD_ARRAY FDPart3tmp302 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp21, FDPart3tmp301));
        const REAL_SIMD_ARRAY FDPart3tmp334 =
            MulSIMD(FDPart3tmp226, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp43, f1_of_xx1__D1)));
        const REAL_SIMD_ARRAY FDPart3tmp12 = FusedMulAddSIMD(FDPart3tmp10, FDPart3tmp7, FDPart3tmp4);
        const REAL_SIMD_ARRAY FDPart3tmp18 = DivSIMD(
            FDPart3_Integer_1,
            FusedMulAddSIMD(FDPart3tmp15, FDPart3tmp7,
                            FusedMulAddSIMD(FDPart3tmp0, MulSIMD(FDPart3tmp10, FDPart3tmp7),
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp13, FDPart3tmp3), MulSIMD(hDD02, hDD12),
                                                            FusedMulAddSIMD(FDPart3tmp0, FDPart3tmp4, MulSIMD(FDPart3tmp10, FDPart3tmp17))))));
        const REAL_SIMD_ARRAY FDPart3tmp224 = FusedMulAddSIMD(FDPart3tmp41, FDPart3tmp45, FDPart3tmp223);
        const REAL_SIMD_ARRAY FDPart3tmp228 = FusedMulAddSIMD(FDPart3tmp8, hDD_dD010, FDPart3tmp227);
        const REAL_SIMD_ARRAY FDPart3tmp72 = MulSIMD(FDPart3tmp18, FDPart3tmp29);
        const REAL_SIMD_ARRAY FDPart3tmp75 = MulSIMD(FDPart3tmp12, FDPart3tmp18);
        const REAL_SIMD_ARRAY FDPart3tmp77 = MulSIMD(FDPart3tmp18, FDPart3tmp34);
        const REAL_SIMD_ARRAY FDPart3tmp80 = MulSIMD(FDPart3tmp18, FDPart3tmp37);
        const REAL_SIMD_ARRAY FDPart3tmp81 = MulSIMD(FDPart3tmp18, FDPart3tmp64);
        const REAL_SIMD_ARRAY FDPart3tmp84 = MulSIMD(FDPart3tmp18, FDPart3tmp54);
        const REAL_SIMD_ARRAY FDPart3tmp79 = FusedMulAddSIMD(
            FDPart3_Rational_1_2,
            MulSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp73, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp26, hDD_dD002))),
            FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp73, hDD_dD010, SubSIMD(FDPart3tmp13, hDD_dD001))),
                            MulSIMD(FDPart3tmp75, FDPart3tmp76)));
        const REAL_SIMD_ARRAY FDPart3tmp82 =
            FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp73, hDD_dD010, SubSIMD(FDPart3tmp13, hDD_dD001))),
                            FusedMulAddSIMD(FDPart3_Rational_1_2,
                                            MulSIMD(FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp73,
                                                                                  FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp26, hDD_dD002))),
                                            MulSIMD(FDPart3tmp72, FDPart3tmp76)));
        const REAL_SIMD_ARRAY FDPart3tmp85 = FusedMulAddSIMD(
            FDPart3_Rational_1_2,
            MulSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp73, FusedMulSubSIMD(FDPart3_Integer_2, FDPart3tmp26, hDD_dD002))),
            FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp84, FusedMulAddSIMD(FDPart3tmp73, hDD_dD010, SubSIMD(FDPart3tmp13, hDD_dD001))),
                            MulSIMD(FDPart3tmp76, FDPart3tmp77)));
        const REAL_SIMD_ARRAY FDPart3tmp88 = MulSIMD(FDPart3_Integer_3, FDPart3tmp75);
        const REAL_SIMD_ARRAY FDPart3tmp89 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp75);
        const REAL_SIMD_ARRAY FDPart3tmp93 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp77);
        const REAL_SIMD_ARRAY FDPart3tmp97 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp72);
        const REAL_SIMD_ARRAY FDPart3tmp100 = MulSIMD(FDPart3_Integer_3, FDPart3tmp77);
        const REAL_SIMD_ARRAY FDPart3tmp109 = MulSIMD(FDPart3_Integer_3, FDPart3tmp72);
        const REAL_SIMD_ARRAY FDPart3tmp110 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp80);
        const REAL_SIMD_ARRAY FDPart3tmp111 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp113 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp119 = MulSIMD(FDPart3_Integer_3, FDPart3tmp80);
        const REAL_SIMD_ARRAY FDPart3tmp121 = MulSIMD(FDPart3_Integer_3, FDPart3tmp84);
        const REAL_SIMD_ARRAY FDPart3tmp129 = MulSIMD(FDPart3_Integer_3, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp98 =
            FusedMulAddSIMD(FDPart3tmp92, FDPart3tmp93, FusedMulAddSIMD(FDPart3tmp96, FDPart3tmp97, MulSIMD(FDPart3tmp89, hDD_dD001)));
        const REAL_SIMD_ARRAY FDPart3tmp107 =
            FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp93, FusedMulAddSIMD(FDPart3tmp89, hDD_dD002, MulSIMD(FDPart3tmp105, FDPart3tmp97)));
        const REAL_SIMD_ARRAY FDPart3tmp112 =
            FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp96, FusedMulAddSIMD(FDPart3tmp97, hDD_dD001, MulSIMD(FDPart3tmp110, FDPart3tmp92)));
        const REAL_SIMD_ARRAY FDPart3tmp114 = FusedMulAddSIMD(
            FDPart3tmp113, FDPart3tmp92, FusedMulAddSIMD(FDPart3tmp93, hDD_dD001, FusedMulSubSIMD(FDPart3tmp110, FDPart3tmp96, FDPart3tmp21)));
        const REAL_SIMD_ARRAY FDPart3tmp122 = FusedMulAddSIMD(
            FDPart3tmp106, FDPart3tmp110, FusedMulAddSIMD(FDPart3tmp97, hDD_dD002, FusedMulSubSIMD(FDPart3tmp105, FDPart3tmp111, FDPart3tmp21)));
        const REAL_SIMD_ARRAY FDPart3tmp123 =
            FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp113, FusedMulAddSIMD(FDPart3tmp93, hDD_dD002, MulSIMD(FDPart3tmp105, FDPart3tmp110)));
        const REAL_SIMD_ARRAY FDPart3tmp131 =
            FusedMulAddSIMD(FDPart3tmp32, FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp7, FDPart3tmp85, MulSIMD(FDPart3tmp130, FDPart3tmp82)));
        const REAL_SIMD_ARRAY FDPart3tmp132 = MulSIMD(FDPart3_Integer_2, FDPart3tmp85);
        const REAL_SIMD_ARRAY FDPart3tmp133 =
            FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp79, FDPart3tmp83, MulSIMD(FDPart3tmp10, FDPart3tmp82)));
        const REAL_SIMD_ARRAY FDPart3tmp134 = MulSIMD(FDPart3_Integer_2, FDPart3tmp82);
        const REAL_SIMD_ARRAY FDPart3tmp149 = FusedMulAddSIMD(
            FDPart3tmp111, FDPart3tmp146, FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp97, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp141, FDPart3tmp45)));
        const REAL_SIMD_ARRAY FDPart3tmp150 =
            FusedMulAddSIMD(FDPart3tmp113, FDPart3tmp141, FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp93, MulSIMD(FDPart3tmp110, FDPart3tmp146)));
        const REAL_SIMD_ARRAY FDPart3tmp151 =
            FusedMulAddSIMD(FDPart3tmp146, FDPart3tmp97, FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp89, MulSIMD(FDPart3tmp141, FDPart3tmp93)));
        const REAL_SIMD_ARRAY FDPart3tmp162 =
            FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp160, FusedMulAddSIMD(FDPart3tmp161, FDPart3tmp97, MulSIMD(FDPart3tmp110, FDPart3tmp156)));
        const REAL_SIMD_ARRAY FDPart3tmp163 =
            FusedMulAddSIMD(FDPart3tmp113, FDPart3tmp156, FusedMulAddSIMD(FDPart3tmp161, FDPart3tmp93, MulSIMD(FDPart3tmp110, FDPart3tmp160)));
        const REAL_SIMD_ARRAY FDPart3tmp164 = FusedMulAddSIMD(
            FDPart3tmp160, FDPart3tmp97, FusedMulAddSIMD(FDPart3tmp161, FDPart3tmp89, FusedMulAddSIMD(FDPart3tmp156, FDPart3tmp93, f0_of_xx0)));
        const REAL_SIMD_ARRAY FDPart3tmp185 =
            FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp181, FusedMulAddSIMD(FDPart3tmp182, FDPart3tmp97, MulSIMD(FDPart3tmp110, FDPart3tmp184)));
        const REAL_SIMD_ARRAY FDPart3tmp186 = FusedMulAddSIMD(
            FDPart3tmp113, FDPart3tmp184, FusedMulAddSIMD(FDPart3tmp182, FDPart3tmp93, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp181, FDPart3tmp65)));
        const REAL_SIMD_ARRAY FDPart3tmp187 = FusedMulAddSIMD(
            FDPart3tmp182, FDPart3tmp89, FusedMulAddSIMD(FDPart3tmp184, FDPart3tmp93, FusedMulAddSIMD(FDPart3tmp181, FDPart3tmp97, FDPart3tmp67)));
        const REAL_SIMD_ARRAY FDPart3tmp87 =
            FusedMulAddSIMD(FDPart3tmp32, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp82, FDPart3tmp83, MulSIMD(FDPart3tmp0, FDPart3tmp79)));
        const REAL_SIMD_ARRAY FDPart3tmp139 =
            FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp130, FusedMulAddSIMD(FDPart3tmp83, FDPart3tmp98, MulSIMD(FDPart3tmp10, FDPart3tmp112)));
        const REAL_SIMD_ARRAY FDPart3tmp152 =
            FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp32, MulSIMD(FDPart3tmp0, FDPart3tmp151)));
        const REAL_SIMD_ARRAY FDPart3tmp165 =
            FusedMulAddSIMD(FDPart3tmp162, FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp32, MulSIMD(FDPart3tmp0, FDPart3tmp164)));
        const REAL_SIMD_ARRAY FDPart3tmp174 =
            FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp130, FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp7, MulSIMD(FDPart3tmp107, FDPart3tmp32)));
        const REAL_SIMD_ARRAY FDPart3tmp188 =
            FusedMulAddSIMD(FDPart3tmp185, FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp32, MulSIMD(FDPart3tmp0, FDPart3tmp187)));
        const REAL_SIMD_ARRAY FDPart3tmp202 = FusedMulAddSIMD(
            FDPart3tmp113, FDPart3tmp164,
            FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp80,
                            FusedMulAddSIMD(FDPart3tmp77, FDPart3tmp98,
                                            FusedMulAddSIMD(FDPart3tmp79, FDPart3tmp89,
                                                            FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp72, MulSIMD(FDPart3tmp111, FDPart3tmp187))))));
        const REAL_SIMD_ARRAY FDPart3tmp255 =
            FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp163, FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp83, MulSIMD(FDPart3tmp10, FDPart3tmp162)));
        const REAL_SIMD_ARRAY FDPart3tmp277 =
            FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp7, FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp32, MulSIMD(FDPart3tmp130, FDPart3tmp185)));
        const REAL_SIMD_ARRAY FDPart3tmp341 = MulSIMD(FDPart3_Integer_2, FDPart3tmp164);
        const REAL_SIMD_ARRAY FDPart3tmp343 = MulSIMD(FDPart3_Integer_2, FDPart3tmp162);
        const REAL_SIMD_ARRAY FDPart3tmp381 = MulSIMD(FDPart3_Integer_2, FDPart3tmp186);
        const REAL_SIMD_ARRAY FDPart3tmp116 =
            FusedMulAddSIMD(FDPart3tmp112, FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp32, MulSIMD(FDPart3tmp0, FDPart3tmp98)));
        const REAL_SIMD_ARRAY FDPart3tmp125 =
            FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp32, MulSIMD(FDPart3tmp0, FDPart3tmp107)));
        const REAL_SIMD_ARRAY FDPart3tmp154 =
            FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp7, FusedMulAddSIMD(FDPart3tmp32, FDPart3tmp98, MulSIMD(FDPart3tmp112, FDPart3tmp130)));
        const REAL_SIMD_ARRAY FDPart3tmp166 = MulSIMD(FDPart3tmp112, FDPart3tmp139);
        const REAL_SIMD_ARRAY FDPart3tmp179 =
            FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp83, FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp130, MulSIMD(FDPart3tmp10, FDPart3tmp122)));
        const REAL_SIMD_ARRAY FDPart3tmp189 = MulSIMD(FDPart3tmp123, FDPart3tmp174);
        const REAL_SIMD_ARRAY FDPart3tmp199 = FusedMulAddSIMD(
            FDPart3tmp113, FDPart3tmp162,
            FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp72,
                            FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp80,
                                            FusedMulAddSIMD(FDPart3tmp82, FDPart3tmp89,
                                                            FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp185, MulSIMD(FDPart3tmp112, FDPart3tmp77))))));
        const REAL_SIMD_ARRAY FDPart3tmp201 = FusedMulAddSIMD(
            FDPart3tmp114, FDPart3tmp77,
            FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp72,
                            FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp80,
                                            FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp89,
                                                            FusedMulAddSIMD(FDPart3tmp111, FDPart3tmp186, MulSIMD(FDPart3tmp113, FDPart3tmp163))))));
        const REAL_SIMD_ARRAY FDPart3tmp242 =
            FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp7, FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp32, MulSIMD(FDPart3tmp130, FDPart3tmp162)));
        const REAL_SIMD_ARRAY FDPart3tmp246 =
            FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp7, FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp32, MulSIMD(FDPart3tmp130, FDPart3tmp149)));
        const REAL_SIMD_ARRAY FDPart3tmp270 = MulSIMD(FDPart3tmp152, FDPart3tmp98);
        const REAL_SIMD_ARRAY FDPart3tmp273 = MulSIMD(FDPart3tmp107, FDPart3tmp152);
        const REAL_SIMD_ARRAY FDPart3tmp275 =
            FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp150, FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp83, MulSIMD(FDPart3tmp10, FDPart3tmp149)));
        const REAL_SIMD_ARRAY FDPart3tmp310 =
            FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp186, FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp83, MulSIMD(FDPart3tmp10, FDPart3tmp185)));
        const REAL_SIMD_ARRAY FDPart3tmp349 = MulSIMD(FDPart3tmp151, FDPart3tmp152);
        const REAL_SIMD_ARRAY FDPart3tmp118 = MulSIMD(FDPart3tmp107, FDPart3tmp116);
        const REAL_SIMD_ARRAY FDPart3tmp120 = MulSIMD(FDPart3tmp116, FDPart3tmp98);
        const REAL_SIMD_ARRAY FDPart3tmp127 = MulSIMD(FDPart3tmp125, FDPart3tmp98);
        const REAL_SIMD_ARRAY FDPart3tmp128 = MulSIMD(FDPart3tmp107, FDPart3tmp125);
        const REAL_SIMD_ARRAY FDPart3tmp168 = MulSIMD(FDPart3tmp123, FDPart3tmp154);
        const REAL_SIMD_ARRAY FDPart3tmp170 = MulSIMD(FDPart3tmp114, FDPart3tmp154);
        const REAL_SIMD_ARRAY FDPart3tmp193 = MulSIMD(FDPart3tmp122, FDPart3tmp179);
        const REAL_SIMD_ARRAY FDPart3tmp248 = MulSIMD(FDPart3tmp114, FDPart3tmp246);
        const REAL_SIMD_ARRAY FDPart3tmp250 = MulSIMD(FDPart3tmp123, FDPart3tmp246);
        const REAL_SIMD_ARRAY FDPart3tmp261 = MulSIMD(FDPart3tmp116, FDPart3tmp151);
        const REAL_SIMD_ARRAY FDPart3tmp271 = MulSIMD(FDPart3tmp125, FDPart3tmp151);
        const REAL_SIMD_ARRAY FDPart3tmp278 = MulSIMD(FDPart3tmp112, FDPart3tmp275);
        const REAL_SIMD_ARRAY FDPart3tmp280 = MulSIMD(FDPart3tmp122, FDPart3tmp275);
        const REAL_SIMD_ARRAY FDPart3tmp311 = MulSIMD(FDPart3_Integer_2, FDPart3tmp310);
        const REAL_SIMD_ARRAY FDPart3tmp315 = FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp139, MulSIMD(FDPart3tmp114, FDPart3tmp174));
        const REAL_SIMD_ARRAY FDPart3tmp317 = FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp98, FDPart3tmp270);
        const REAL_SIMD_ARRAY FDPart3tmp321 = FusedMulAddSIMD(FDPart3tmp179, FDPart3tmp98, MulSIMD(FDPart3tmp188, FDPart3tmp98));
        const REAL_SIMD_ARRAY FDPart3tmp322 = FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp275, MulSIMD(FDPart3tmp114, FDPart3tmp277));
        const REAL_SIMD_ARRAY FDPart3tmp340 = MulSIMD(FDPart3tmp150, FDPart3tmp246);
        const REAL_SIMD_ARRAY FDPart3tmp351 = MulSIMD(FDPart3tmp149, FDPart3tmp275);
        const REAL_SIMD_ARRAY FDPart3tmp307 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp112, FDPart3tmp179));
        const REAL_SIMD_ARRAY FDPart3tmp313 = FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp98, FDPart3tmp127);
        const REAL_SIMD_ARRAY FDPart3tmp319 = FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp255, FDPart3tmp248);
        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(
            FDPart3tmp18,
            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp54, AddSIMD(hDD_dDD0011, NegFusedMulAddSIMD(FDPart3_Integer_2, hDD_dD011, FDPart3tmp53)))),
            FusedMulAddSIMD(
                FDPart3tmp18,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp37, AddSIMD(FDPart3tmp46, NegFusedMulAddSIMD(FDPart3_Integer_2, hDD_dD012, FDPart3tmp42)))),
                FusedMulAddSIMD(
                    FDPart3tmp18,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp37,
                                    AddSIMD(FDPart3tmp63, NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp59), FDPart3tmp46)))),
                    FusedMulAddSIMD(
                        FDPart3tmp18,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp34, NegFusedMulAddSIMD(FDPart3_Integer_2, hDD_dD010, FDPart3tmp31))),
                        FusedMulAddSIMD(
                            FDPart3tmp18,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp34, FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp21,
                                                                          NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp56),
                                                                                             FDPart3tmp31)))),
                            FusedMulAddSIMD(
                                FDPart3tmp18,
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                        MulSIMD(FDPart3tmp29, NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f1_of_xx1, hDD_dD020), FDPart3tmp23))),
                                FusedMulAddSIMD(
                                    FDPart3tmp18,
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                            MulSIMD(FDPart3tmp29,
                                                    FusedMulAddSIMD(
                                                        FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp26),
                                                        NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp48), FDPart3tmp23)))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp109, MulSIMD(FDPart3tmp125, FDPart3tmp79),
                                        FusedMulAddSIMD(
                                            FDPart3tmp79, MulSIMD(FDPart3tmp87, FDPart3tmp88),
                                            FusedMulAddSIMD(
                                                FDPart3tmp100, MulSIMD(FDPart3tmp87, FDPart3tmp98),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp107, MulSIMD(FDPart3tmp109, FDPart3tmp87),
                                                    FusedMulAddSIMD(
                                                        hDD02, lambdaU_dD20,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp100, MulSIMD(FDPart3tmp116, FDPart3tmp79),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp84,
                                                                FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp114, FDPart3tmp165),
                                                                                FDPart3tmp170),
                                                                FusedMulAddSIMD(
                                                                    hDD01, lambdaU_dD10,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp81,
                                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp122, FDPart3tmp188),
                                                                                        FDPart3tmp193),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp84,
                                                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp112, FDPart3tmp152),
                                                                                            FDPart3tmp166),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp80,
                                                                                FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp139,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp122, FDPart3tmp152))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp81,
                                                                                    FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                    MulSIMD(FDPart3tmp123, FDPart3tmp152),
                                                                                                    FDPart3tmp189),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp80,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp112, FDPart3tmp179,
                                                                                            MulSIMD(FDPart3_Integer_2,
                                                                                                    MulSIMD(FDPart3tmp112, FDPart3tmp188))),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp80,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp114, FDPart3tmp174,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp114, FDPart3tmp152))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp77,
                                                                                                FusedMulAddSIMD(FDPart3tmp134, FDPart3tmp152,
                                                                                                                MulSIMD(FDPart3tmp139, FDPart3tmp82)),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp80,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp123, FDPart3tmp165),
                                                                                                        FDPart3tmp168),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp77,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp114, FDPart3tmp131,
                                                                                                            MulSIMD(FDPart3_Integer_2,
                                                                                                                    MulSIMD(FDPart3tmp114,
                                                                                                                            FDPart3tmp116))),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp77,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp132, FDPart3tmp165,
                                                                                                                MulSIMD(FDPart3tmp154, FDPart3tmp85)),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp75,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp125, FDPart3tmp134,
                                                                                                                    MulSIMD(FDPart3tmp133,
                                                                                                                            FDPart3tmp82)),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp77,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp112, FDPart3tmp133,
                                                                                                                        MulSIMD(
                                                                                                                            FDPart3_Integer_2,
                                                                                                                            MulSIMD(FDPart3tmp112,
                                                                                                                                    FDPart3tmp125))),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp72,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp134,
                                                                                                                            FDPart3tmp188,
                                                                                                                            MulSIMD(FDPart3tmp179,
                                                                                                                                    FDPart3tmp82)),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp75,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp116,
                                                                                                                                FDPart3tmp132,
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3tmp131,
                                                                                                                                    FDPart3tmp85)),
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp72,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp123,
                                                                                                                                    FDPart3tmp131,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp116,
                                                                                                                                            FDPart3tmp123))),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp72,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp132,
                                                                                                                                        FDPart3tmp152,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp174,
                                                                                                                                            FDPart3tmp85)),
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp202,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp196,
                                                                                                                                            FDPart3tmp85,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3_Integer_2,
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp0,
                                                                                                                                                    FDPart3tmp79),
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp195,
                                                                                                                                                    FDPart3tmp82))),
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp72,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp122,
                                                                                                                                                FDPart3tmp133,
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                                    MulSIMD(
                                                                                                                                                        FDPart3tmp122,
                                                                                                                                                        FDPart3tmp125))),
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp199,
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3tmp123,
                                                                                                                                                    FDPart3tmp196,
                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp0,
                                                                                                                                                            FDPart3tmp107),
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp122,
                                                                                                                                                            FDPart3tmp195))),
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3tmp201,
                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                        FDPart3tmp114,
                                                                                                                                                        FDPart3tmp196,
                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                            FDPart3_Integer_2,
                                                                                                                                                            MulSIMD(
                                                                                                                                                                FDPart3tmp0,
                                                                                                                                                                FDPart3tmp98),
                                                                                                                                                            MulSIMD(
                                                                                                                                                                FDPart3tmp112,
                                                                                                                                                                FDPart3tmp195))),
                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                        FDPart3tmp120,
                                                                                                                                                        FDPart3tmp121,
                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp128,
                                                                                                                                                                        FDPart3tmp129,
                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp118,
                                                                                                                                                                                        FDPart3tmp119,
                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp119,
                                                                                                                                                                                                        FDPart3tmp127,
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp18,
                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                            MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp64, AddSIMD(
                                                                                                                                                                                                                                                                                                           hDD_dDD0022,
                                                                                                                                                                                                                                                                                                           FusedMulAddSIMD(FDPart3tmp30, FDPart3tmp65, NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f1_of_xx1, hDD_dD022), FDPart3tmp71))))),
                                                                                                                                                                                                                        FusedMulSubSIMD(FDPart3tmp0,
                                                                                                                                                                                                                                        lambdaU_dD00,
                                                                                                                                                                                                                                        MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp12), MulSIMD(
                                                                                                                                                                                                                                                                                                 FDPart3tmp18, hDD_dDD0000)))))))))))))))))))))))))))))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(
            FDPart3tmp18,
            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp37, AddSIMD(AddSIMD(FDPart3tmp221, FDPart3tmp232), FusedMulAddSIMD(FDPart3tmp20, f0_of_xx0, FDPart3tmp239)))),
            FusedMulAddSIMD(
                FDPart3tmp18,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp54, FusedMulAddSIMD(f0_of_xx0, hDD_dD001,
                                                              FusedMulAddSIMD(FDPart3tmp30, f0_of_xx0,
                                                                              FusedMulAddSIMD(FDPart3tmp5, hDD_dD010,
                                                                                              FusedMulAddSIMD(f0_of_xx0, hDD_dDD0111,
                                                                                                              NegFusedMulAddSIMD(f0_of_xx0, hDD_dD111,
                                                                                                                                 FDPart3tmp219))))))),
                FusedMulAddSIMD(
                    FDPart3tmp18,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp34,
                                    AddSIMD(hDD_dD011, AddSIMD(AddSIMD(hDD00, hDD11),
                                                               FusedMulAddSIMD(f0_of_xx0, hDD_dDD0101,
                                                                               NegFusedMulAddSIMD(FDPart3tmp21, FDPart3tmp91, FDPart3tmp53)))))),
                    FusedMulAddSIMD(
                        FDPart3tmp18,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp37,
                                        AddSIMD(FDPart3tmp224,
                                                FusedMulAddSIMD(f0_of_xx0, hDD_dD002, NegFusedMulAddSIMD(f0_of_xx0, hDD_dD112, FDPart3tmp221))))),
                        FusedMulAddSIMD(
                            FDPart3tmp18,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp29,
                                            AddSIMD(AddSIMD(FDPart3tmp220, FDPart3tmp63), FusedMulAddSIMD(f0_of_xx0, hDD_dDD0102, hDD_dD012)))),
                            FusedMulAddSIMD(
                                FDPart3tmp18,
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                        MulSIMD(FDPart3tmp34,
                                                FusedMulAddSIMD(f0_of_xx0, hDD_dDD0101,
                                                                SubSIMD(FusedMulSubSIMD(f0_of_xx0, hDD_dD000, MulSIMD(FDPart3tmp21, FDPart3tmp52)),
                                                                        MulSIMD(f0_of_xx0, hDD_dD110))))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3_Rational_1_2, f0_of_xx0), MulSIMD(hDD12, lambdaU_dD20),
                                    FusedMulAddSIMD(
                                        FDPart3tmp18,
                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                MulSIMD(FDPart3tmp29, FusedMulAddSIMD(FDPart3tmp214, hDD_dD020,
                                                                                      FusedMulAddSIMD(f0_of_xx0, hDD_dDD0102,
                                                                                                      NegFusedMulAddSIMD(FDPart3tmp21, FDPart3tmp62,
                                                                                                                         FDPart3tmp213))))),
                                        FusedMulAddSIMD(
                                            FDPart3_Rational_1_2, MulSIMD(hDD02, lambdaU_dD21),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3_Rational_1_2, FDPart3tmp21), MulSIMD(FDPart3tmp7, lambdaU_dD10),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp84,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp154, FDPart3tmp98,
                                                        FusedMulAddSIMD(FDPart3tmp165, FDPart3tmp98, MulSIMD(FDPart3tmp116, FDPart3tmp164))),
                                                    FusedMulAddSIMD(
                                                        FDPart3_Rational_1_2, MulSIMD(FDPart3tmp0, FDPart3tmp207),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp84,
                                                            FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp165,
                                                                            MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp114, FDPart3tmp242))),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp84,
                                                                FusedMulAddSIMD(FDPart3tmp112, FDPart3tmp255,
                                                                                FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp162,
                                                                                                MulSIMD(FDPart3tmp112, FDPart3tmp246))),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp81,
                                                                    AddSIMD(FDPart3tmp271,
                                                                            FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp174, FDPart3tmp273)),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp81,
                                                                        FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp277,
                                                                                        FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp188, FDPart3tmp280)),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp80,
                                                                            FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp255,
                                                                                            FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp152,
                                                                                                            MulSIMD(FDPart3tmp122, FDPart3tmp246))),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp81,
                                                                                FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp250,
                                                                                                MulSIMD(FDPart3tmp150, FDPart3tmp152)),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp80,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp125, FDPart3tmp164,
                                                                                        FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp98, FDPart3tmp270)),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp80,
                                                                                        FusedMulAddSIMD(FDPart3tmp112, FDPart3tmp277,
                                                                                                        FusedMulAddSIMD(FDPart3tmp162, FDPart3tmp188,
                                                                                                                        FDPart3tmp278)),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp80,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp150, FDPart3tmp165,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp123, FDPart3tmp242))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp80,
                                                                                                FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp154,
                                                                                                                FusedMulAddSIMD(FDPart3tmp107,
                                                                                                                                FDPart3tmp165,
                                                                                                                                FDPart3tmp261)),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp77,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp246, FDPart3tmp82,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp255, FDPart3tmp82,
                                                                                                            MulSIMD(FDPart3tmp112, FDPart3tmp152))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp80,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3_Integer_2, FDPart3tmp248,
                                                                                                            MulSIMD(FDPart3tmp152, FDPart3tmp163)),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp77,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp154, FDPart3tmp79,
                                                                                                                FusedMulAddSIMD(FDPart3tmp165,
                                                                                                                                FDPart3tmp79,
                                                                                                                                FDPart3tmp120)),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp77,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp112, FDPart3tmp174,
                                                                                                                    FusedMulAddSIMD(FDPart3tmp125,
                                                                                                                                    FDPart3tmp162,
                                                                                                                                    FDPart3tmp166)),
                                                                                                                FusedMulAddSIMD(FDPart3tmp77,
                                                                                                                                FusedMulAddSIMD(FDPart3tmp114,
                                                                                                                                                FDPart3tmp165,
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp132,
                                                                                                                                                    FDPart3tmp242)),
                                                                                                                                FusedMulAddSIMD(FDPart3tmp77,
                                                                                                                                                FusedMulAddSIMD(FDPart3tmp131,
                                                                                                                                                                FDPart3tmp98,
                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp164,
                                                                                                                                                                                FDPart3tmp87,
                                                                                                                                                                                FDPart3tmp120)),
                                                                                                                                                FusedMulAddSIMD(FDPart3tmp75,
                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp131,
                                                                                                                                                                                FDPart3tmp79,
                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp87,
                                                                                                                                                                                                FDPart3tmp98,
                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                    FDPart3tmp116,
                                                                                                                                                                                                    FDPart3tmp79))),
                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp77,
                                                                                                                                                                                FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                FDPart3tmp170,
                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                    FDPart3tmp116,
                                                                                                                                                                                                    FDPart3tmp163)),
                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp75,
                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp116, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp154, FDPart3tmp85))),
                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp75,
                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp139,
                                                                                                                                                                                                                                FDPart3tmp82,
                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                    FDPart3tmp174, FDPart3tmp82,
                                                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                                                        FDPart3tmp112,
                                                                                                                                                                                                                                        FDPart3tmp125))),
                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                    FDPart3tmp72,
                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp275, FDPart3tmp82, FusedMulAddSIMD(FDPart3tmp277, FDPart3tmp82, MulSIMD(FDPart3tmp112, FDPart3tmp188))),
                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp72,
                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp122,
                                                                                                                                                                                                                                                    FDPart3tmp174,
                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                        FDPart3tmp125, FDPart3tmp149, MulSIMD(FDPart3tmp122, FDPart3tmp139))),
                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp72,
                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp107,
                                                                                                                                                                                                                                                                    FDPart3tmp131,
                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp151,
                                                                                                                                                                                                                                                                                    FDPart3tmp87,
                                                                                                                                                                                                                                                                                    FDPart3tmp118)),
                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp72,
                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp79,
                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                        FDPart3tmp174,
                                                                                                                                                                                                                                                                                        FDPart3tmp79,
                                                                                                                                                                                                                                                                                        FDPart3tmp127)),
                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp72,
                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp168, MulSIMD(FDPart3tmp116, FDPart3tmp150)),
                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp152, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp246, FDPart3tmp85))),
                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp204,
                                                                                                                                                                                                                                                                                                                    FDPart3tmp210,
                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp204,
                                                                                                                                                                                                                                                                                                                                    lambdaU_dD00,
                                                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                        FDPart3tmp201, AddSIMD(FDPart3tmp154, FDPart3tmp165),
                                                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                            FDPart3tmp202,
                                                                                                                                                                                                                                                                                                                                            AddSIMD(
                                                                                                                                                                                                                                                                                                                                                FDPart3tmp116,
                                                                                                                                                                                                                                                                                                                                                FDPart3tmp131),
                                                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp18,
                                                                                                                                                                                                                                                                                                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp64, FusedMulAddSIMD(FDPart3tmp214, hDD_dD022, FusedMulAddSIMD(f0_of_xx0, hDD_dDD0122, AddSIMD(FDPart3tmp228, FDPart3tmp230))))),
                                                                                                                                                                                                                                                                                                                                                            FusedMulSubSIMD(FDPart3tmp199, AddSIMD(FDPart3tmp152, FDPart3tmp174),
                                                                                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                MulSIMD(FDPart3_Rational_1_2, FDPart3tmp12), MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                 FDPart3tmp18,
                                                                                                                                                                                                                                                                                                                                                                                                                                 FusedMulAddSIMD(FDPart3tmp21,
                                                                                                                                                                                                                                                                                                                                                                                                                                                 hDD01,
                                                                                                                                                                                                                                                                                                                                                                                                                                                 FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                                                     f0_of_xx0, hDD_dDD0100, NegFusedMulAddSIMD(FDPart3tmp21, FDPart3tmp56, hDD_dD010)))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_2 = FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp37, FusedMulAddSIMD(FDPart3tmp65, MulSIMD(f0_of_xx0, hDD_dD011), FusedMulAddSIMD(f0_of_xx0, MulSIMD(f1_of_xx1__D1, hDD_dD022), SubSIMD(FusedMulAddSIMD(FDPart3tmp65, MulSIMD(FDPart3tmp73, hDD00), FusedMulAddSIMD(FDPart3tmp289, hDD_dDD0212, FusedMulAddSIMD(FDPart3tmp67, hDD_dD001, FusedMulAddSIMD(FDPart3tmp43, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp69, f1_of_xx1__D1)), NegFusedMulAddSIMD(FDPart3tmp145, FDPart3tmp21, FDPart3tmp227))))), MulSIMD(FDPart3tmp306, FDPart3tmp32)))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp54, FusedMulAddSIMD(FDPart3tmp289, hDD_dDD0211, FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp5, FusedMulAddSIMD(FDPart3tmp73, MulSIMD(f1_of_xx1__D1, hDD_dD021), FusedMulAddSIMD(f0_of_xx0, MulSIMD(f1_of_xx1__DD11, hDD02), AddSIMD(FDPart3tmp224, FDPart3tmp239))))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp34, AddSIMD(AddSIMD(FDPart3tmp220, FDPart3tmp291), FusedMulAddSIMD(f0_of_xx0, MulSIMD(f1_of_xx1__D1, hDD_dD020), FDPart3tmp42)))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp37, FusedMulAddSIMD(FDPart3tmp289, hDD_dDD0212, FusedMulAddSIMD(FDPart3tmp30, FDPart3tmp67, AddSIMD(FDPart3tmp230, FDPart3tmp293))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp29, FusedMulAddSIMD(FDPart3tmp2, hDD00, FusedMulAddSIMD(FDPart3tmp289, hDD_dDD0202, AddSIMD(AddSIMD(FDPart3tmp66, FDPart3tmp71), FusedMulAddSIMD(FDPart3tmp56, FDPart3tmp65, NegFusedMulAddSIMD(FDPart3tmp104, FDPart3tmp21, FDPart3tmp102))))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp34, AddSIMD(FDPart3tmp291, SubSIMD(NegFusedMulAddSIMD(FDPart3tmp21, FDPart3tmp41, FDPart3tmp213), MulSIMD(FDPart3tmp21, FDPart3tmp59))))), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp25), MulSIMD(f0_of_xx0, lambdaU_dD10), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp29, FusedMulAddSIMD(FDPart3tmp67, hDD_dD000, SubSIMD(FusedMulAddSIMD(FDPart3tmp65, MulSIMD(f0_of_xx0, hDD_dD010), FusedMulSubSIMD(FDPart3tmp289, hDD_dDD0202, MulSIMD(FDPart3tmp21, FDPart3tmp69))), MulSIMD(FDPart3tmp67, hDD_dD220))))), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp0, FDPart3tmp284), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp10), MulSIMD(FDPart3tmp286, lambdaU_dD20), FusedMulAddSIMD(FDPart3tmp84, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp165, FDPart3tmp319), FusedMulAddSIMD(FDPart3tmp84, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp278, MulSIMD(FDPart3tmp149, FDPart3tmp152)), FusedMulAddSIMD(FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp277, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp186, MulSIMD(FDPart3tmp123, FDPart3tmp275))), FusedMulAddSIMD(FDPart3tmp84, AddSIMD(FDPart3tmp261, FDPart3tmp317), FusedMulAddSIMD(FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp311, MulSIMD(FDPart3tmp185, FDPart3tmp188)), FusedMulAddSIMD(FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp188, FusedMulAddSIMD(FDPart3tmp125, FDPart3tmp187, MulSIMD(FDPart3tmp107, FDPart3tmp179))), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp255, FusedMulAddSIMD(FDPart3tmp165, FDPart3tmp186, FDPart3tmp250)), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp139, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp187, FDPart3tmp273)), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp280, MulSIMD(FDPart3tmp152, FDPart3tmp185)), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp112, FDPart3tmp311, MulSIMD(FDPart3tmp149, FDPart3tmp188)), FusedMulAddSIMD(FDPart3tmp80, AddSIMD(FDPart3tmp271, FDPart3tmp321), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp152, FDPart3tmp322), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp79, FDPart3tmp118)), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp246, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp255, FDPart3tmp85, MulSIMD(FDPart3tmp123, FDPart3tmp165))), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp150, FDPart3tmp315), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp152, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp275, FDPart3tmp82))), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp125, FDPart3tmp149, FDPart3tmp307), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp87, FDPart3tmp313), FusedMulAddSIMD(FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp125, FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp79, MulSIMD(FDPart3tmp107, FDPart3tmp87))), FusedMulAddSIMD(FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp85, MulSIMD(FDPart3tmp116, FDPart3tmp123))), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp275, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp277, FDPart3tmp85, MulSIMD(FDPart3tmp123, FDPart3tmp152))), FusedMulAddSIMD(FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp125, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp179, FDPart3tmp82))), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp179, FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp188, FDPart3tmp79, FDPart3tmp128)), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp186, FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp139, FDPart3tmp189)), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp188, MulSIMD(FDPart3tmp134, FDPart3tmp310)), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp133, FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp87, FDPart3tmp128)), FusedMulAddSIMD(FDPart3tmp283, lambdaU_dD00, FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp193, MulSIMD(FDPart3tmp125, FDPart3tmp185)), FusedMulAddSIMD(FDPart3tmp204, FDPart3tmp285, FusedMulAddSIMD(FDPart3tmp283, FDPart3tmp288, FusedMulAddSIMD(FDPart3tmp201, AddSIMD(FDPart3tmp139, FDPart3tmp152), FusedMulAddSIMD(FDPart3tmp202, AddSIMD(FDPart3tmp125, FDPart3tmp133), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp64, FusedMulAddSIMD(FDPart3tmp60, FDPart3tmp65, FusedMulAddSIMD(FDPart3tmp62, FDPart3tmp65, FusedMulAddSIMD(FDPart3tmp297, hDD_dD020, FusedMulAddSIMD(FDPart3tmp41, FDPart3tmp65, FusedMulAddSIMD(FDPart3tmp20, FDPart3tmp67, FusedMulAddSIMD(FDPart3tmp289, hDD_dDD0222, FusedMulAddSIMD(FDPart3tmp67, hDD_dD002, NegFusedMulAddSIMD(FDPart3tmp67, hDD_dD222, FDPart3tmp302)))))))))), FusedMulSubSIMD(FDPart3tmp199, AddSIMD(FDPart3tmp179, FDPart3tmp188), MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp12), MulSIMD(FDPart3tmp18, FusedMulAddSIMD(FDPart3tmp21, FDPart3tmp26, FusedMulAddSIMD(FDPart3tmp289, hDD_dDD0200, NegFusedMulAddSIMD(FDPart3tmp21, FDPart3tmp48, FDPart3tmp47)))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_3 = FusedMulAddSIMD(
    FDPart3tmp18,
    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
            MulSIMD(FDPart3tmp54, FusedMulAddSIMD(FDPart3tmp5, hDD_dDD1111,
                                                  FusedMulAddSIMD(FDPart3tmp52, FDPart3tmp73,
                                                                  FusedMulAddSIMD(FDPart3tmp142, hDD_dD011, MulSIMD(FDPart3tmp24, hDD_dD110)))))),
    FusedMulAddSIMD(
        FDPart3tmp18,
        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp37, FusedMulAddSIMD(FDPart3tmp142, hDD_dD012,
                                                      NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp222),
                                                                         MulSIMD(FDPart3tmp43, f1_of_xx1__D1), FDPart3tmp332)))),
        FusedMulAddSIMD(
            FDPart3tmp18,
            MulSIMD(
                MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp37,
                        FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp335,
                                        FusedMulAddSIMD(FDPart3tmp238, FDPart3tmp36,
                                                        FusedMulAddSIMD(FDPart3tmp62, FDPart3tmp73,
                                                                        NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp233),
                                                                                           MulSIMD(FDPart3tmp43, f1_of_xx1__D1), FDPart3tmp332)))))),
            FusedMulAddSIMD(
                FDPart3tmp18,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp34, FusedMulAddSIMD(FDPart3tmp5, hDD_dDD1101,
                                                              FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp5, hDD_dD010), FDPart3tmp219)))),
                FusedMulAddSIMD(
                    FDPart3tmp18,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp34,
                                    FusedMulAddSIMD(
                                        FDPart3tmp5, hDD_dDD1101,
                                        FusedMulAddSIMD(FDPart3tmp56, FDPart3tmp73,
                                                        FusedMulAddSIMD(FDPart3tmp73, hDD_dD111,
                                                                        NegFusedMulAddSIMD(FDPart3_Integer_3, MulSIMD(FDPart3tmp21, FDPart3tmp218),
                                                                                           FDPart3tmp196)))))),
                    FusedMulAddSIMD(
                        FDPart3tmp18,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp29, FusedMulAddSIMD(FDPart3tmp5, hDD_dDD1102,
                                                                      NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp5),
                                                                                         MulSIMD(f1_of_xx1__D1, hDD_dD120), FDPart3tmp232)))),
                        FusedMulAddSIMD(
                            FDPart3tmp18,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp29,
                                            FusedMulAddSIMD(
                                                FDPart3tmp73, hDD_dD112,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp43,
                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp95, f1_of_xx1__D1)),
                                                    FusedMulSubSIMD(FDPart3tmp5, hDD_dDD1102,
                                                                    MulSIMD(FDPart3_Integer_3, MulSIMD(FDPart3tmp21, FDPart3tmp231))))))),
                            FusedMulAddSIMD(
                                f0_of_xx0, MulSIMD(hDD12, lambdaU_dD21),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3_Integer_3, FDPart3tmp114), MulSIMD(FDPart3tmp242, FDPart3tmp77),
                                    FusedMulAddSIMD(
                                        FDPart3tmp119, MulSIMD(FDPart3tmp163, FDPart3tmp246),
                                        FusedMulAddSIMD(
                                            FDPart3tmp121, MulSIMD(FDPart3tmp163, FDPart3tmp242),
                                            FusedMulAddSIMD(
                                                FDPart3tmp109, MulSIMD(FDPart3tmp150, FDPart3tmp154),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp119, MulSIMD(FDPart3tmp150, FDPart3tmp242),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp84,
                                                        FusedMulAddSIMD(FDPart3tmp162, FDPart3tmp255, MulSIMD(FDPart3tmp246, FDPart3tmp343)),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp100, MulSIMD(FDPart3tmp154, FDPart3tmp163),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp81,
                                                                FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp149, FDPart3tmp277),
                                                                                FDPart3tmp351),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp84,
                                                                    FusedMulAddSIMD(FDPart3tmp154, FDPart3tmp341,
                                                                                    MulSIMD(FDPart3tmp164, FDPart3tmp165)),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp80,
                                                                        FusedMulAddSIMD(FDPart3tmp162, FDPart3tmp275,
                                                                                        MulSIMD(FDPart3tmp277, FDPart3tmp343)),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp81,
                                                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp151, FDPart3tmp174),
                                                                                            FDPart3tmp349),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp80,
                                                                                FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp165,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp151, FDPart3tmp154))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp80,
                                                                                    FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp164,
                                                                                                    MulSIMD(FDPart3tmp174, FDPart3tmp341)),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp77,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp165, FDPart3tmp98,
                                                                                            MulSIMD(FDPart3_Integer_2,
                                                                                                    MulSIMD(FDPart3tmp154, FDPart3tmp98))),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp80,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp149, FDPart3tmp255,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp149, FDPart3tmp246))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp77,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp116, FDPart3tmp164,
                                                                                                    MulSIMD(FDPart3tmp131, FDPart3tmp341)),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp77,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp139, FDPart3tmp162,
                                                                                                        MulSIMD(FDPart3tmp174, FDPart3tmp343)),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp75,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp112, FDPart3tmp174),
                                                                                                            FDPart3tmp166),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp77,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp112, FDPart3tmp255,
                                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                                        MulSIMD(FDPart3tmp112,
                                                                                                                                FDPart3tmp246))),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp72,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp139, FDPart3tmp149,
                                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                                            MulSIMD(FDPart3tmp149,
                                                                                                                                    FDPart3tmp174))),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp75,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3_Integer_2,
                                                                                                                        MulSIMD(FDPart3tmp131,
                                                                                                                                FDPart3tmp98),
                                                                                                                        FDPart3tmp120),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp72,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3_Integer_2,
                                                                                                                            MulSIMD(FDPart3tmp174,
                                                                                                                                    FDPart3tmp98),
                                                                                                                            FDPart3tmp270),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp72,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3tmp112,
                                                                                                                                    FDPart3tmp277),
                                                                                                                                FDPart3tmp278),
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp210,
                                                                                                                                FDPart3tmp7,
                                                                                                                                FusedMulAddSIMD(FDPart3tmp72,
                                                                                                                                                FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                MulSIMD(
                                                                                                                                                                    FDPart3tmp131,
                                                                                                                                                                    FDPart3tmp151),
                                                                                                                                                                FDPart3tmp261),
                                                                                                                                                FusedMulAddSIMD(FDPart3tmp202,
                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp196,
                                                                                                                                                                                FDPart3tmp98,
                                                                                                                                                                                FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                    FDPart3tmp114,
                                                                                                                                                                                                    FDPart3tmp7),
                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                    FDPart3tmp112,
                                                                                                                                                                                                    FDPart3tmp353))),
                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp207,
                                                                                                                                                                                FDPart3tmp32,
                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp199,
                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp151,
                                                                                                                                                                                                                FDPart3tmp196, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp150, FDPart3tmp7), MulSIMD(FDPart3tmp149, FDPart3tmp353))),
                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                    FDPart3tmp201, FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp196, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp163, FDPart3tmp7), MulSIMD(FDPart3tmp130, FDPart3tmp343))),
                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                        FDPart3tmp129,
                                                                                                                                                                                                        FDPart3tmp340,
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp170,
                                                                                                                                                                                                                        FDPart3tmp88,
                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp18,
                                                                                                                                                                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(
                                                                                                                                                                                                                                                                                                         FDPart3tmp64, FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                           FDPart3tmp218,
                                                                                                                                                                                                                                                                                                                           FDPart3tmp65,
                                                                                                                                                                                                                                                                                                                           FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                               FDPart3tmp5,
                                                                                                                                                                                                                                                                                                                               hDD_dDD1122,
                                                                                                                                                                                                                                                                                                                               FusedMulAddSIMD(FDPart3tmp2, MulSIMD(FDPart3tmp24, hDD_dD110),
                                                                                                                                                                                                                                                                                                                                               NegFusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                   MulSIMD(FDPart3_Integer_2, FDPart3tmp5), MulSIMD(f1_of_xx1__D1, hDD_dD122), FDPart3tmp334)))))),
                                                                                                                                                                                                                                        FusedMulSubSIMD(FDPart3tmp109,
                                                                                                                                                                                                                                                        FDPart3tmp248, MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp12), MulSIMD(FDPart3tmp18, FusedMulAddSIMD(FDPart3tmp5, hDD_dDD1100, FusedMulAddSIMD(FDPart3tmp73, hDD_dD110, FusedMulSubSIMD(FDPart3_Integer_4, hDD11, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp91))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_4 = FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp37, FusedMulAddSIMD(FDPart3tmp13, MulSIMD(FDPart3tmp5, FDPart3tmp65), FusedMulAddSIMD(FDPart3tmp8, hDD_dD011, FusedMulAddSIMD(FDPart3tmp94, hDD_dDD1212, FusedMulAddSIMD(FDPart3tmp335, FDPart3tmp9, FusedMulAddSIMD(FDPart3tmp69, f0_of_xx0, FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp45, FusedMulAddSIMD(FDPart3tmp156, FDPart3tmp65, FusedMulAddSIMD(FDPart3tmp5, MulSIMD(f1_of_xx1__D1, hDD_dD122), NegFusedMulAddSIMD(FDPart3tmp306, FDPart3tmp6, FDPart3tmp334))))))))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp54, FusedMulAddSIMD(FDPart3tmp142, MulSIMD(f1_of_xx1__D1, hDD_dD121), FusedMulAddSIMD(FDPart3tmp59, f0_of_xx0, FusedMulAddSIMD(FDPart3tmp94, hDD_dDD1211, FusedMulAddSIMD(FDPart3tmp238, FDPart3tmp36, FusedMulAddSIMD(FDPart3tmp41, f0_of_xx0, FusedMulAddSIMD(FDPart3tmp24, MulSIMD(f1_of_xx1, hDD_dD120), FusedMulAddSIMD(FDPart3tmp5, MulSIMD(f1_of_xx1__DD11, hDD12), FusedMulAddSIMD(FDPart3tmp222, FDPart3tmp45, MulSIMD(FDPart3tmp233, FDPart3tmp45))))))))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp34, FusedMulAddSIMD(FDPart3tmp48, f0_of_xx0, AddSIMD(FusedMulAddSIMD(FDPart3tmp45, FDPart3tmp95, FDPart3tmp83), FusedMulAddSIMD(FDPart3tmp5, MulSIMD(f1_of_xx1__D1, hDD_dD120), NegFusedMulAddSIMD(FDPart3_Integer_3, MulSIMD(FDPart3tmp21, FDPart3tmp222), FDPart3tmp356)))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp37, FusedMulAddSIMD(FDPart3tmp226, FDPart3tmp45, FusedMulAddSIMD(FDPart3tmp5, FDPart3tmp66, FusedMulAddSIMD(FDPart3tmp52, FDPart3tmp67, FusedMulAddSIMD(FDPart3tmp94, hDD_dDD1212, FusedMulAddSIMD(FDPart3tmp218, FDPart3tmp65, MulSIMD(FDPart3tmp225, hDD_dD221)))))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp29, FusedMulAddSIMD(FDPart3tmp67, hDD01, FusedMulAddSIMD(FDPart3tmp94, hDD_dDD1202, FusedMulAddSIMD(FDPart3tmp56, FDPart3tmp67, FusedMulAddSIMD(FDPart3tmp65, FDPart3tmp91, FusedMulAddSIMD(FDPart3tmp73, MulSIMD(f1_of_xx1, hDD_dD122), FusedMulSubSIMD(FDPart3tmp104, FDPart3tmp45, MulSIMD(FDPart3_Integer_3, MulSIMD(FDPart3tmp21, FDPart3tmp226)))))))))), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp34, AddSIMD(FDPart3tmp356, FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp5, NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp233), FDPart3tmp223))))), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp10), MulSIMD(FDPart3tmp286, lambdaU_dD21), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp29, FusedMulAddSIMD(FDPart3tmp65, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp94, hDD_dDD1202, FusedMulAddSIMD(FDPart3tmp225, hDD_dD220, FDPart3tmp228))))), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp130, FDPart3tmp288), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp285, FDPart3tmp7), FusedMulAddSIMD(FDPart3tmp84, FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp246, FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp255, MulSIMD(FDPart3tmp150, FDPart3tmp242))), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp130, FDPart3tmp210), FusedMulAddSIMD(FDPart3tmp84, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp246, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp162, FDPart3tmp275))), FusedMulAddSIMD(FDPart3tmp84, FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp154, FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp164, MulSIMD(FDPart3tmp139, FDPart3tmp164))), FusedMulAddSIMD(FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp277, FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp246, MulSIMD(FDPart3tmp150, FDPart3tmp275))), FusedMulAddSIMD(FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp188, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp187, MulSIMD(FDPart3tmp151, FDPart3tmp179))), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp179, FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp188, MulSIMD(FDPart3tmp151, FDPart3tmp174))), FusedMulAddSIMD(FDPart3tmp81, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp311, MulSIMD(FDPart3tmp185, FDPart3tmp277)), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp275, FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp277, FDPart3tmp340)), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp151, FusedMulAddSIMD(FDPart3tmp154, FDPart3tmp187, FDPart3tmp349)), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp277, MulSIMD(FDPart3tmp162, FDPart3tmp311)), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp255, FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp242, FDPart3tmp340)), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp154, FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp174, MulSIMD(FDPart3tmp139, FDPart3tmp163))), FusedMulAddSIMD(FDPart3tmp80, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp351, MulSIMD(FDPart3tmp185, FDPart3tmp246)), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp174, MulSIMD(FDPart3tmp179, FDPart3tmp343)), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp151, FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp164, MulSIMD(FDPart3tmp125, FDPart3tmp164))), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp242, FDPart3tmp319), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp278, MulSIMD(FDPart3tmp122, FDPart3tmp246)), FusedMulAddSIMD(FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp131, FDPart3tmp313), FusedMulAddSIMD(FDPart3tmp77, FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp154, FDPart3tmp317), FusedMulAddSIMD(FDPart3tmp75, AddSIMD(FDPart3tmp168, FDPart3tmp315), FusedMulAddSIMD(FDPart3tmp75, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp174, FDPart3tmp307), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp131, FDPart3tmp187, FusedMulAddSIMD(FDPart3tmp133, FDPart3tmp151, FDPart3tmp271)), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp174, FusedMulAddSIMD(FDPart3tmp154, FDPart3tmp186, MulSIMD(FDPart3tmp139, FDPart3tmp150))), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp112, FDPart3tmp311, MulSIMD(FDPart3tmp122, FDPart3tmp277)), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp174, FDPart3tmp185, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp149, FDPart3tmp179))), FusedMulAddSIMD(FDPart3tmp72, AddSIMD(FDPart3tmp250, FDPart3tmp322), FusedMulAddSIMD(FDPart3tmp72, FusedMulAddSIMD(FDPart3tmp107, FDPart3tmp174, FDPart3tmp321), FusedMulAddSIMD(FDPart3tmp204, FDPart3tmp284, FusedMulAddSIMD(FDPart3tmp207, FDPart3tmp283, FusedMulAddSIMD(FDPart3tmp201, AddSIMD(FDPart3tmp246, FDPart3tmp255), FusedMulAddSIMD(FDPart3tmp202, AddSIMD(FDPart3tmp139, FDPart3tmp174), FusedMulAddSIMD(FDPart3tmp18, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp64, FusedMulAddSIMD(FDPart3tmp8, hDD_dD012, FusedMulAddSIMD(FDPart3tmp301, FDPart3tmp45, FusedMulAddSIMD(FDPart3tmp62, FDPart3tmp67, FusedMulAddSIMD(FDPart3tmp225, hDD_dD222, FusedMulAddSIMD(FDPart3tmp231, FDPart3tmp65, FusedMulAddSIMD(FDPart3tmp94, hDD_dDD1222, FusedMulAddSIMD(FDPart3tmp24, MulSIMD(FDPart3tmp296, hDD_dD120), FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp65, MulSIMD(FDPart3tmp222, FDPart3tmp65))))))))))), FusedMulSubSIMD(FDPart3tmp199, AddSIMD(FDPart3tmp275, FDPart3tmp277), MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp12), MulSIMD(FDPart3tmp18, FusedMulAddSIMD(FDPart3tmp147, FDPart3tmp73, FusedMulAddSIMD(FDPart3tmp94, hDD_dDD1200, FusedMulSubSIMD(FDPart3_Integer_4, FDPart3tmp25, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp21, FDPart3tmp95)))))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_5 = FusedMulAddSIMD(
    FDPart3tmp18,
    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
            MulSIMD(FDPart3tmp54,
                    FusedMulAddSIMD(
                        hDD22, FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp236, MulSIMD(FDPart3tmp142, FDPart3tmp305)),
                        FusedMulAddSIMD(FDPart3tmp2, MulSIMD(FDPart3tmp24, hDD_dD220),
                                        FusedMulAddSIMD(FDPart3tmp335, FDPart3tmp9,
                                                        FusedMulAddSIMD(FDPart3tmp8, hDD_dDD2211,
                                                                        FusedMulAddSIMD(FDPart3tmp145,
                                                                                        MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                                                MulSIMD(FDPart3tmp43, f1_of_xx1__D1)),
                                                                                        FusedMulSubSIMD(FDPart3tmp143, hDD_dD221,
                                                                                                        MulSIMD(FDPart3tmp238, FDPart3tmp9))))))))),
    FusedMulAddSIMD(
        FDPart3tmp18,
        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp37,
                        FusedMulAddSIMD(FDPart3tmp8, hDD_dDD2212,
                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp222, FDPart3tmp65),
                                                        FusedMulAddSIMD(FDPart3tmp101, FDPart3tmp41, MulSIMD(FDPart3tmp301, FDPart3tmp45)))))),
        FusedMulAddSIMD(
            FDPart3tmp18,
            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp37,
                            FusedMulAddSIMD(
                                FDPart3tmp8, hDD_dDD2212,
                                FusedMulAddSIMD(
                                    FDPart3tmp233, FDPart3tmp377,
                                    FusedMulAddSIMD(
                                        FDPart3tmp306, FDPart3tmp36,
                                        FusedMulAddSIMD(
                                            FDPart3tmp130, AddSIMD(FDPart3tmp236, FDPart3tmp305),
                                            FusedMulAddSIMD(FDPart3tmp143, hDD_dD222,
                                                            FusedMulAddSIMD(FDPart3_Integer_4, MulSIMD(FDPart3tmp57, FDPart3tmp8),
                                                                            FusedMulSubSIMD(FDPart3tmp101, FDPart3tmp59,
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_3, FDPart3tmp301),
                                                                                                    MulSIMD(FDPart3tmp43, f1_of_xx1__D1))))))))))),
            FusedMulAddSIMD(
                FDPart3tmp18,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp34,
                                FusedMulAddSIMD(FDPart3tmp67, hDD_dD221,
                                                NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp145, FDPart3tmp21), FDPart3tmp374)))),
                FusedMulAddSIMD(
                    FDPart3tmp18,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp34, AddSIMD(FDPart3tmp374,
                                                          FusedMulAddSIMD(FDPart3tmp143, hDD_dD220,
                                                                          NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp104),
                                                                                             MulSIMD(FDPart3tmp43, f1_of_xx1__D1), FDPart3tmp293))))),
                    FusedMulAddSIMD(
                        FDPart3tmp18,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp29, FusedMulAddSIMD(FDPart3tmp8, hDD_dDD2202,
                                                                      FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp297, hDD_dD020),
                                                                                      FusedMulAddSIMD(FDPart3tmp300, hDD_dD120, FDPart3tmp302))))),
                        FusedMulAddSIMD(
                            FDPart3tmp18,
                            MulSIMD(
                                MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp29,
                                        FusedMulAddSIMD(
                                            FDPart3tmp8, hDD_dDD2202,
                                            FusedMulAddSIMD(
                                                FDPart3tmp101, hDD_dD222,
                                                FusedMulAddSIMD(FDPart3tmp377, FDPart3tmp95,
                                                                FusedMulAddSIMD(FDPart3tmp296, MulSIMD(FDPart3tmp73, hDD02),
                                                                                FusedMulSubSIMD(FDPart3tmp101, FDPart3tmp48,
                                                                                                MulSIMD(FDPart3_Integer_3,
                                                                                                        MulSIMD(FDPart3tmp21, FDPart3tmp301))))))))),
                            FusedMulAddSIMD(
                                MulSIMD(FDPart3_Integer_3, FDPart3tmp122), MulSIMD(FDPart3tmp310, FDPart3tmp72),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3_Integer_3, FDPart3tmp149), MulSIMD(FDPart3tmp310, FDPart3tmp80),
                                    FusedMulAddSIMD(
                                        FDPart3tmp119, MulSIMD(FDPart3tmp185, FDPart3tmp275),
                                        FusedMulAddSIMD(
                                            FDPart3tmp129, MulSIMD(FDPart3tmp185, FDPart3tmp310),
                                            FusedMulAddSIMD(
                                                FDPart3tmp100, MulSIMD(FDPart3tmp149, FDPart3tmp179),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp109, MulSIMD(FDPart3tmp179, FDPart3tmp185),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp84,
                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp150, FDPart3tmp255), FDPart3tmp340),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp84,
                                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp139, FDPart3tmp151), FDPart3tmp349),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp81,
                                                                FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp277, MulSIMD(FDPart3tmp275, FDPart3tmp381)),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp81,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp187, FDPart3tmp188,
                                                                        MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp179, FDPart3tmp187))),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp80,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp152, FDPart3tmp187,
                                                                            MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp139, FDPart3tmp187))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp80,
                                                                            FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp246,
                                                                                            MulSIMD(FDPart3tmp255, FDPart3tmp381)),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp80,
                                                                                FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp277,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp150, FDPart3tmp275))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp80,
                                                                                    FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp188,
                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp151, FDPart3tmp179))),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp77,
                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp107, FDPart3tmp139),
                                                                                                        FDPart3tmp273),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp77,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp150, FDPart3tmp174,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp139, FDPart3tmp150))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp77,
                                                                                                FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                MulSIMD(FDPart3tmp123, FDPart3tmp255),
                                                                                                                FDPart3tmp250),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp77,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp133, FDPart3tmp151),
                                                                                                        FDPart3tmp271),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp75,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp107, FDPart3tmp133),
                                                                                                            FDPart3tmp128),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp75,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3_Integer_2,
                                                                                                                MulSIMD(FDPart3tmp123, FDPart3tmp139),
                                                                                                                FDPart3tmp189),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp72,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp125, FDPart3tmp187,
                                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                                            MulSIMD(FDPart3tmp133,
                                                                                                                                    FDPart3tmp187))),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp72,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp139, FDPart3tmp381,
                                                                                                                        MulSIMD(FDPart3tmp174,
                                                                                                                                FDPart3tmp186)),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp72,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp107,
                                                                                                                            FDPart3tmp188,
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3tmp107,
                                                                                                                                    FDPart3tmp179))),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp72,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp123,
                                                                                                                                FDPart3tmp277,
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp123,
                                                                                                                                        FDPart3tmp275))),
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp202,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp123,
                                                                                                                                    FDPart3tmp353,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp10,
                                                                                                                                            FDPart3tmp122),
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp107,
                                                                                                                                            FDPart3tmp195))),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp284,
                                                                                                                                    FDPart3tmp83,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp199,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp187,
                                                                                                                                            FDPart3tmp195,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3_Integer_2,
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp10,
                                                                                                                                                    FDPart3tmp185),
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp186,
                                                                                                                                                    FDPart3tmp353))),
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp201,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp151,
                                                                                                                                                FDPart3tmp195,
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                                    MulSIMD(
                                                                                                                                                        FDPart3tmp10,
                                                                                                                                                        FDPart3tmp149),
                                                                                                                                                    MulSIMD(
                                                                                                                                                        FDPart3tmp150,
                                                                                                                                                        FDPart3tmp353))),
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp130,
                                                                                                                                                FDPart3tmp285,
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3tmp193,
                                                                                                                                                    FDPart3tmp88,
                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp100,
                                                                                                                                                                    FDPart3tmp280,
                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp121,
                                                                                                                                                                                    FDPart3tmp351,
                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp18,
                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                        MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp64,
                                                                                                                                                                                                                                                                     FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp297, hDD_dD022), FusedMulAddSIMD(FDPart3tmp300, hDD_dD122, FusedMulAddSIMD(FDPart3tmp8, hDD_dDD2222, FusedMulAddSIMD(FDPart3tmp24, MulSIMD(hDD_dD220, MulSIMD(MulSIMD(MulSIMD(f1_of_xx1, f1_of_xx1), f1_of_xx1), f1_of_xx1)), FusedMulAddSIMD(FDPart3tmp297, MulSIMD(f1_of_xx1__D1, hDD_dD221), FusedMulAddSIMD(FDPart3tmp101, FDPart3tmp69, MulSIMD(FDPart3tmp226, FDPart3tmp377))))))))),
                                                                                                                                                                                                    FusedMulSubSIMD(
                                                                                                                                                                                                        FDPart3tmp10,
                                                                                                                                                                                                        FDPart3tmp288,
                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                            MulSIMD(FDPart3_Rational_1_2, FDPart3tmp12), MulSIMD(FDPart3tmp18, FusedMulAddSIMD(FDPart3tmp101, hDD_dD220, FusedMulAddSIMD(FDPart3tmp8, hDD_dDD2200, FusedMulSubSIMD(FDPart3_Integer_4, FDPart3tmp102, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp104, FDPart3tmp21))))))))))))))))))))))))))))))))))))))))))))))));

WriteSIMD(&auxevol_gfs[IDX4(RBARDD00GF, i0, i1, i2)], __RHS_exp_0);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD01GF, i0, i1, i2)], __RHS_exp_1);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD02GF, i0, i1, i2)], __RHS_exp_2);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD11GF, i0, i1, i2)], __RHS_exp_3);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD12GF, i0, i1, i2)], __RHS_exp_4);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD22GF, i0, i1, i2)], __RHS_exp_5);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += simd_width)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
