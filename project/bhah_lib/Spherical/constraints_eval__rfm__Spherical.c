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
 * Evaluate BSSN constraints.
 */
void constraints_eval__rfm__Spherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                      const rfm_struct *restrict rfmstruct, const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs,
                                      REAL *restrict diagnostic_output_gfs) {
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
        const REAL_SIMD_ARRAY aDD00_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD00_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD00_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD00_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD00GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD01_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD01_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD01_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD01_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD01GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD02_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD02_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD02_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD02_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD02GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD11_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD11_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD11_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD11_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD11GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD12_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD12_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD12_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD12_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD12GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY aDD22_i2m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY aDD22_i2m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY aDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY aDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY aDD22_i2p1 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY aDD22_i2p2 = ReadSIMD(&in_gfs[IDX4(ADD22GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY alpha = ReadSIMD(&in_gfs[IDX4(ALPHAGF, i0, i1, i2)]);
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
        const REAL_SIMD_ARRAY cf_i0m2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0m1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY cf = ReadSIMD(&in_gfs[IDX4(CFGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p1 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY cf_i0p2 = ReadSIMD(&in_gfs[IDX4(CFGF, i0 + 2, i1, i2)]);
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
        const REAL_SIMD_ARRAY hDD00_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD00_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD00_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD00_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD00_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD00_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD00GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD01_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD01_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD01_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD01_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD01_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD01_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD01GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD02_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD02_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD02_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD02_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD02_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD02_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD02GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD11_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD11_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD11_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD11_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD11_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD11_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD11GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD12_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD12_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD12_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD12_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD12_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD12_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD12GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY hDD22_i2m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY hDD22_i2m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY hDD22_i1m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0m1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i0p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY hDD22_i1p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY hDD22_i2p1 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY hDD22_i2p2 = ReadSIMD(&in_gfs[IDX4(HDD22GF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY trK_i2m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 - 2)]);
        const REAL_SIMD_ARRAY trK_i2m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 - 1)]);
        const REAL_SIMD_ARRAY trK_i1m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 2, i2)]);
        const REAL_SIMD_ARRAY trK_i1m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 - 1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 2, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0m1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 - 1, i1, i2)]);
        const REAL_SIMD_ARRAY trK = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 1, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i0p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0 + 2, i1, i2)]);
        const REAL_SIMD_ARRAY trK_i1p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 1, i2)]);
        const REAL_SIMD_ARRAY trK_i1p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1 + 2, i2)]);
        const REAL_SIMD_ARRAY trK_i2p1 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 + 1)]);
        const REAL_SIMD_ARRAY trK_i2p2 = ReadSIMD(&in_gfs[IDX4(TRKGF, i0, i1, i2 + 2)]);
        const REAL_SIMD_ARRAY vetU0 = ReadSIMD(&in_gfs[IDX4(VETU0GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY vetU1 = ReadSIMD(&in_gfs[IDX4(VETU1GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY vetU2 = ReadSIMD(&in_gfs[IDX4(VETU2GF, i0, i1, i2)]);
        const REAL_SIMD_ARRAY aDD_dD000 = SIMD_fd_function_dD0_fdorder4(aDD00_i0m1, aDD00_i0m2, aDD00_i0p1, aDD00_i0p2, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD001 = SIMD_fd_function_dD1_fdorder4(aDD00_i1m1, aDD00_i1m2, aDD00_i1p1, aDD00_i1p2, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD002 = SIMD_fd_function_dD2_fdorder4(aDD00_i2m1, aDD00_i2m2, aDD00_i2p1, aDD00_i2p2, invdxx2);
        const REAL_SIMD_ARRAY aDD_dD010 = SIMD_fd_function_dD0_fdorder4(aDD01_i0m1, aDD01_i0m2, aDD01_i0p1, aDD01_i0p2, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD011 = SIMD_fd_function_dD1_fdorder4(aDD01_i1m1, aDD01_i1m2, aDD01_i1p1, aDD01_i1p2, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD012 = SIMD_fd_function_dD2_fdorder4(aDD01_i2m1, aDD01_i2m2, aDD01_i2p1, aDD01_i2p2, invdxx2);
        const REAL_SIMD_ARRAY aDD_dD020 = SIMD_fd_function_dD0_fdorder4(aDD02_i0m1, aDD02_i0m2, aDD02_i0p1, aDD02_i0p2, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD021 = SIMD_fd_function_dD1_fdorder4(aDD02_i1m1, aDD02_i1m2, aDD02_i1p1, aDD02_i1p2, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD022 = SIMD_fd_function_dD2_fdorder4(aDD02_i2m1, aDD02_i2m2, aDD02_i2p1, aDD02_i2p2, invdxx2);
        const REAL_SIMD_ARRAY aDD_dD110 = SIMD_fd_function_dD0_fdorder4(aDD11_i0m1, aDD11_i0m2, aDD11_i0p1, aDD11_i0p2, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD111 = SIMD_fd_function_dD1_fdorder4(aDD11_i1m1, aDD11_i1m2, aDD11_i1p1, aDD11_i1p2, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD112 = SIMD_fd_function_dD2_fdorder4(aDD11_i2m1, aDD11_i2m2, aDD11_i2p1, aDD11_i2p2, invdxx2);
        const REAL_SIMD_ARRAY aDD_dD120 = SIMD_fd_function_dD0_fdorder4(aDD12_i0m1, aDD12_i0m2, aDD12_i0p1, aDD12_i0p2, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD121 = SIMD_fd_function_dD1_fdorder4(aDD12_i1m1, aDD12_i1m2, aDD12_i1p1, aDD12_i1p2, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD122 = SIMD_fd_function_dD2_fdorder4(aDD12_i2m1, aDD12_i2m2, aDD12_i2p1, aDD12_i2p2, invdxx2);
        const REAL_SIMD_ARRAY aDD_dD220 = SIMD_fd_function_dD0_fdorder4(aDD22_i0m1, aDD22_i0m2, aDD22_i0p1, aDD22_i0p2, invdxx0);
        const REAL_SIMD_ARRAY aDD_dD221 = SIMD_fd_function_dD1_fdorder4(aDD22_i1m1, aDD22_i1m2, aDD22_i1p1, aDD22_i1p2, invdxx1);
        const REAL_SIMD_ARRAY aDD_dD222 = SIMD_fd_function_dD2_fdorder4(aDD22_i2m1, aDD22_i2m2, aDD22_i2p1, aDD22_i2p2, invdxx2);
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

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
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

        const double dblFDPart3_Integer_4 = 4.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_4 = ConstSIMD(dblFDPart3_Integer_4);

        const double dblFDPart3_Integer_6 = 6.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_6 = ConstSIMD(dblFDPart3_Integer_6);

        const double dblFDPart3_Integer_8 = 8.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_8 = ConstSIMD(dblFDPart3_Integer_8);

        const double dblFDPart3_NegativeOne_ = -1.0;
        const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        const double dblFDPart3_Rational_2_3 = 2.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_2_3 = ConstSIMD(dblFDPart3_Rational_2_3);

        const REAL_SIMD_ARRAY FDPart3tmp1 = MulSIMD(FDPart3_Integer_2, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp2 = MulSIMD(f1_of_xx1, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp5 = AddSIMD(FDPart3_Integer_1, hDD00);
        const REAL_SIMD_ARRAY FDPart3tmp8 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp30 = MulSIMD(hDD12, MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp33 = MulSIMD(f0_of_xx0, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp45 = MulSIMD(FDPart3_Integer_2, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp53 = MulSIMD(aDD02, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp90 = MulSIMD(aDD01, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp122 = MulSIMD(cf, cf);
        const REAL_SIMD_ARRAY FDPart3tmp131 = DivSIMD(FDPart3_Integer_1, cf);
        const REAL_SIMD_ARRAY FDPart3tmp144 = MulSIMD(f0_of_xx0, hDD_dD012);
        const REAL_SIMD_ARRAY FDPart3tmp179 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf));
        const REAL_SIMD_ARRAY FDPart3tmp181 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf), cf), cf));
        const REAL_SIMD_ARRAY FDPart3tmp183 = DivSIMD(FDPart3_Integer_1, MulSIMD(cf, cf));
        const REAL_SIMD_ARRAY FDPart3tmp193 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(T4UU00, alpha));
        const REAL_SIMD_ARRAY FDPart3tmp197 = DivSIMD(FDPart3_Integer_1, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp3 = MulSIMD(FDPart3tmp2, MulSIMD(MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0), f0_of_xx0));
        const REAL_SIMD_ARRAY FDPart3tmp9 = FusedMulAddSIMD(FDPart3tmp8, hDD11, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp10 = MulSIMD(FDPart3tmp2, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp22 = MulSIMD(FDPart3tmp8, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp28 = MulSIMD(FDPart3tmp8, aDD11);
        const REAL_SIMD_ARRAY FDPart3tmp31 = MulSIMD(FDPart3tmp30, MulSIMD(f1_of_xx1, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp34 = MulSIMD(FDPart3tmp33, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp46 = MulSIMD(FDPart3tmp45, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp54 = MulSIMD(FDPart3tmp45, FDPart3tmp53);
        const REAL_SIMD_ARRAY FDPart3tmp57 = MulSIMD(FDPart3tmp2, MulSIMD(FDPart3tmp30, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp59 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(f0_of_xx0, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp87 = MulSIMD(FDPart3tmp53, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp125 = DivSIMD(FDPart3_Integer_1, FDPart3tmp122);
        const REAL_SIMD_ARRAY FDPart3tmp132 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp131);
        const REAL_SIMD_ARRAY FDPart3tmp134 = FusedMulAddSIMD(FDPart3tmp45, hDD_dD010, SubSIMD(FDPart3tmp1, hDD_dD001));
        const REAL_SIMD_ARRAY FDPart3tmp136 = MulSIMD(FDPart3tmp45, f1_of_xx1);
        const REAL_SIMD_ARRAY FDPart3tmp147 = MulSIMD(f0_of_xx0, MulSIMD(f1_of_xx1__D1, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp155 = MulSIMD(FDPart3tmp2, FDPart3tmp45);
        const REAL_SIMD_ARRAY FDPart3tmp184 = MulSIMD(FDPart3tmp183, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp189 = MulSIMD(FDPart3tmp183, vetU1);
        const REAL_SIMD_ARRAY FDPart3tmp201 = MulSIMD(FDPart3_Integer_12, FDPart3tmp90);
        const REAL_SIMD_ARRAY FDPart3tmp4 = MulSIMD(MulSIMD(FDPart3tmp1, FDPart3tmp3), MulSIMD(hDD02, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp6 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp3, MulSIMD(hDD12, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp14 = FusedMulAddSIMD(FDPart3tmp10, hDD22, FDPart3tmp10);
        const REAL_SIMD_ARRAY FDPart3tmp16 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp8, MulSIMD(hDD01, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp23 = MulSIMD(FDPart3tmp22, MulSIMD(hDD01, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp24 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp8, f1_of_xx1));
        const REAL_SIMD_ARRAY FDPart3tmp43 = MulSIMD(FDPart3tmp10, aDD22);
        const REAL_SIMD_ARRAY FDPart3tmp49 = MulSIMD(FDPart3_Integer_2, FDPart3tmp22);
        const REAL_SIMD_ARRAY FDPart3tmp85 = MulSIMD(FDPart3tmp22, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp137 =
            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f1_of_xx1, hDD02), FusedMulSubSIMD(FDPart3tmp136, hDD_dD020, hDD_dD002));
        const REAL_SIMD_ARRAY FDPart3tmp143 = FusedMulAddSIMD(FDPart3tmp45, hDD11, FusedMulAddSIMD(FDPart3tmp8, hDD_dD110, FDPart3tmp45));
        const REAL_SIMD_ARRAY FDPart3tmp148 = FusedMulAddSIMD(FDPart3tmp33, hDD_dD021, FDPart3tmp147);
        const REAL_SIMD_ARRAY FDPart3tmp149 = FusedMulAddSIMD(FDPart3tmp136, hDD12, MulSIMD(FDPart3tmp22, hDD_dD120));
        const REAL_SIMD_ARRAY FDPart3tmp163 =
            SubSIMD(NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, hDD11),
                                       FusedMulSubSIMD(FDPart3tmp45, hDD_dD011, MulSIMD(FDPart3_Integer_2, f0_of_xx0))),
                    MulSIMD(FDPart3tmp8, hDD_dD110));
        const REAL_SIMD_ARRAY FDPart3tmp185 = MulSIMD(FDPart3tmp184, T4UU02);
        const REAL_SIMD_ARRAY FDPart3tmp187 = MulSIMD(FDPart3tmp184, T4UU01);
        const REAL_SIMD_ARRAY FDPart3tmp188 = MulSIMD(FDPart3tmp184, T4UU03);
        const REAL_SIMD_ARRAY FDPart3tmp191 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp183, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp192 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp183, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp12 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp10, MulSIMD(hDD02, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp18 = MulSIMD(FDPart3tmp14, FDPart3tmp9);
        const REAL_SIMD_ARRAY FDPart3tmp25 = MulSIMD(FDPart3tmp24, MulSIMD(FDPart3tmp5, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp35 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp34, FDPart3tmp9));
        const REAL_SIMD_ARRAY FDPart3tmp41 = FusedMulAddSIMD(FDPart3tmp5, FDPart3tmp9, FDPart3tmp16);
        const REAL_SIMD_ARRAY FDPart3tmp50 = MulSIMD(FDPart3tmp49, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp150 = AddSIMD(FDPart3tmp149, SubSIMD(FDPart3tmp148, FDPart3tmp144));
        const REAL_SIMD_ARRAY FDPart3tmp154 =
            AddSIMD(FDPart3tmp149, SubSIMD(SubSIMD(FDPart3tmp144, FDPart3tmp147), MulSIMD(FDPart3tmp33, hDD_dD021)));
        const REAL_SIMD_ARRAY FDPart3tmp157 = FusedMulAddSIMD(FDPart3tmp10, hDD_dD220, FusedMulAddSIMD(FDPart3tmp155, hDD22, FDPart3tmp155));
        const REAL_SIMD_ARRAY FDPart3tmp162 = FusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp8), MulSIMD(f1_of_xx1__D1, hDD12),
                                                              FusedMulSubSIMD(FDPart3tmp49, hDD_dD121, MulSIMD(FDPart3tmp8, hDD_dD112)));
        const REAL_SIMD_ARRAY FDPart3tmp167 = MulSIMD(FDPart3tmp49, f1_of_xx1__D1);
        const REAL_SIMD_ARRAY FDPart3tmp170 = AddSIMD(
            FDPart3tmp148, FusedMulAddSIMD(FDPart3tmp24, hDD_dD120,
                                           NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, f0_of_xx0), MulSIMD(f1_of_xx1, hDD12), FDPart3tmp144)));
        const REAL_SIMD_ARRAY FDPart3tmp174 =
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, f0_of_xx0),
                               FusedMulAddSIMD(FDPart3tmp2, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0, hDD22)),
                                               FusedMulSubSIMD(FDPart3tmp136, hDD_dD022, MulSIMD(FDPart3tmp10, hDD_dD220))));
        const REAL_SIMD_ARRAY FDPart3tmp175 = FusedMulAddSIMD(
            FDPart3tmp8, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1, f1_of_xx1__D1)),
            FusedMulAddSIMD(MulSIMD(FDPart3tmp8, f1_of_xx1), MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f1_of_xx1__D1, hDD22)),
                            FusedMulSubSIMD(FDPart3tmp49, hDD_dD122, MulSIMD(FDPart3tmp10, hDD_dD221))));
        const REAL_SIMD_ARRAY FDPart3tmp202 = MulSIMD(FDPart3_Integer_12, FDPart3tmp85);
        const REAL_SIMD_ARRAY FDPart3tmp26 = AddSIMD(FDPart3tmp23, FDPart3tmp25);
        const REAL_SIMD_ARRAY FDPart3tmp36 = AddSIMD(FDPart3tmp31, FDPart3tmp35);
        const REAL_SIMD_ARRAY FDPart3tmp61 = FusedMulAddSIMD(FDPart3tmp14, FDPart3tmp59, FDPart3tmp57);
        const REAL_SIMD_ARRAY FDPart3tmp65 = FusedMulAddSIMD(FDPart3tmp14, FDPart3tmp5, FDPart3tmp12);
        const REAL_SIMD_ARRAY FDPart3tmp76 = AddSIMD(FDPart3tmp18, FDPart3tmp6);
        const REAL_SIMD_ARRAY FDPart3tmp169 = FusedMulAddSIMD(FDPart3tmp10, hDD_dD221, FusedMulAddSIMD(FDPart3tmp167, hDD22, FDPart3tmp167));
        const REAL_SIMD_ARRAY FDPart3tmp180 = FusedMulAddSIMD(FDPart3tmp179, FDPart3tmp31, MulSIMD(FDPart3tmp179, FDPart3tmp35));
        const REAL_SIMD_ARRAY FDPart3tmp195 = FusedMulAddSIMD(FDPart3tmp179, FDPart3tmp23, MulSIMD(FDPart3tmp179, FDPart3tmp25));
        const REAL_SIMD_ARRAY FDPart3tmp239 =
            FusedMulAddSIMD(FDPart3tmp179, FDPart3tmp57, MulSIMD(FDPart3tmp14, MulSIMD(FDPart3tmp179, FDPart3tmp59)));
        const REAL_SIMD_ARRAY FDPart3tmp20 =
            FusedMulAddSIMD(FDPart3tmp14, FDPart3tmp16,
                            FusedMulAddSIMD(FDPart3tmp18, FDPart3tmp5,
                                            FusedMulAddSIMD(FDPart3tmp5, FDPart3tmp6, FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp9, FDPart3tmp4))));
        const REAL_SIMD_ARRAY FDPart3tmp182 =
            DivSIMD(FDPart3_Integer_1,
                    FusedMulAddSIMD(FDPart3tmp14, MulSIMD(FDPart3tmp16, FDPart3tmp181),
                                    FusedMulAddSIMD(FDPart3tmp18, MulSIMD(FDPart3tmp181, FDPart3tmp5),
                                                    FusedMulAddSIMD(FDPart3tmp181, MulSIMD(FDPart3tmp5, FDPart3tmp6),
                                                                    FusedMulAddSIMD(FDPart3tmp181, FDPart3tmp4,
                                                                                    MulSIMD(FDPart3tmp12, MulSIMD(FDPart3tmp181, FDPart3tmp9)))))));
        const REAL_SIMD_ARRAY FDPart3tmp21 = DivSIMD(FDPart3_Integer_1, MulSIMD(FDPart3tmp20, FDPart3tmp20));
        const REAL_SIMD_ARRAY FDPart3tmp123 = DivSIMD(FDPart3_Integer_1, FDPart3tmp20);
        const REAL_SIMD_ARRAY FDPart3tmp27 = MulSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp26, FDPart3tmp26));
        const REAL_SIMD_ARRAY FDPart3tmp38 = MulSIMD(FDPart3tmp21, aDD00);
        const REAL_SIMD_ARRAY FDPart3tmp44 = MulSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp43, MulSIMD(FDPart3tmp41, FDPart3tmp41)));
        const REAL_SIMD_ARRAY FDPart3tmp47 = MulSIMD(FDPart3tmp21, FDPart3tmp26);
        const REAL_SIMD_ARRAY FDPart3tmp51 = MulSIMD(FDPart3tmp21, FDPart3tmp41);
        const REAL_SIMD_ARRAY FDPart3tmp66 = MulSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp65, FDPart3tmp65));
        const REAL_SIMD_ARRAY FDPart3tmp69 = MulSIMD(FDPart3tmp21, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp72 = MulSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp36, FDPart3tmp36));
        const REAL_SIMD_ARRAY FDPart3tmp74 = MulSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp61, FDPart3tmp61));
        const REAL_SIMD_ARRAY FDPart3tmp80 = MulSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp36, FDPart3tmp61));
        const REAL_SIMD_ARRAY FDPart3tmp124 = MulSIMD(FDPart3_Integer_2, FDPart3tmp123);
        const REAL_SIMD_ARRAY FDPart3tmp138 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp36, hDD_dD000)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp137, FDPart3tmp41)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp134, FDPart3tmp26))));
        const REAL_SIMD_ARRAY FDPart3tmp139 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp61, hDD_dD000)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp137, FDPart3tmp26)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp134, FDPart3tmp65))));
        const REAL_SIMD_ARRAY FDPart3tmp140 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp76, hDD_dD000)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp137, FDPart3tmp36)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp134, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp151 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp36, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp150, FDPart3tmp41)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp143, FDPart3tmp26))));
        const REAL_SIMD_ARRAY FDPart3tmp152 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp61, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp150, FDPart3tmp26)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp143, FDPart3tmp65))));
        const REAL_SIMD_ARRAY FDPart3tmp153 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp76, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp150, FDPart3tmp36)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp143, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp158 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp36, hDD_dD002)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp157, FDPart3tmp41)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp154, FDPart3tmp26))));
        const REAL_SIMD_ARRAY FDPart3tmp159 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp61, hDD_dD002)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp157, FDPart3tmp26)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp154, FDPart3tmp65))));
        const REAL_SIMD_ARRAY FDPart3tmp160 = FusedMulAddSIMD(
            FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp76, hDD_dD002)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp157, FDPart3tmp36)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp154, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp164 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp123, FDPart3tmp26), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp163, FDPart3tmp36)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp162, FDPart3tmp41))));
        const REAL_SIMD_ARRAY FDPart3tmp165 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp123, FDPart3tmp65), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp163, FDPart3tmp61)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp162, FDPart3tmp26))));
        const REAL_SIMD_ARRAY FDPart3tmp166 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp123, FDPart3tmp61), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp163, FDPart3tmp76)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp162, FDPart3tmp36))));
        const REAL_SIMD_ARRAY FDPart3tmp171 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp123, FDPart3tmp26), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD112)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp170, FDPart3tmp36)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp169, FDPart3tmp41))));
        const REAL_SIMD_ARRAY FDPart3tmp172 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp123, FDPart3tmp65), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD112)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp170, FDPart3tmp61)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp169, FDPart3tmp26))));
        const REAL_SIMD_ARRAY FDPart3tmp173 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp123, FDPart3tmp61), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD112)),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp170, FDPart3tmp76)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp169, FDPart3tmp36))));
        const REAL_SIMD_ARRAY FDPart3tmp176 = FusedMulAddSIMD(
            FDPart3tmp41,
            MulSIMD(MulSIMD(FDPart3tmp123, FDPart3tmp2),
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD222))),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp175, FDPart3tmp26)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp174, FDPart3tmp36))));
        const REAL_SIMD_ARRAY FDPart3tmp177 = FusedMulAddSIMD(
            FDPart3tmp26,
            MulSIMD(MulSIMD(FDPart3tmp123, FDPart3tmp2),
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD222))),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp175, FDPart3tmp65)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp174, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp178 = FusedMulAddSIMD(
            FDPart3tmp36,
            MulSIMD(MulSIMD(FDPart3tmp123, FDPart3tmp2),
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp8, hDD_dD222))),
            FusedMulSubSIMD(FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp175, FDPart3tmp61)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp123), MulSIMD(FDPart3tmp174, FDPart3tmp76))));
        const REAL_SIMD_ARRAY FDPart3tmp39 = MulSIMD(FDPart3tmp38, MulSIMD(FDPart3tmp36, FDPart3tmp36));
        const REAL_SIMD_ARRAY FDPart3tmp48 = MulSIMD(FDPart3tmp36, FDPart3tmp47);
        const REAL_SIMD_ARRAY FDPart3tmp52 = MulSIMD(FDPart3tmp26, FDPart3tmp51);
        const REAL_SIMD_ARRAY FDPart3tmp63 = MulSIMD(FDPart3tmp38, MulSIMD(FDPart3tmp61, FDPart3tmp61));
        const REAL_SIMD_ARRAY FDPart3tmp68 = MulSIMD(FDPart3tmp47, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp70 = MulSIMD(FDPart3tmp26, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp71 = MulSIMD(FDPart3tmp61, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp78 = MulSIMD(FDPart3tmp38, MulSIMD(FDPart3tmp76, FDPart3tmp76));
        const REAL_SIMD_ARRAY FDPart3tmp82 = MulSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp36, FDPart3tmp76));
        const REAL_SIMD_ARRAY FDPart3tmp83 = MulSIMD(FDPart3tmp21, MulSIMD(FDPart3tmp61, FDPart3tmp76));
        const REAL_SIMD_ARRAY FDPart3tmp92 = MulSIMD(FDPart3tmp36, MulSIMD(FDPart3tmp38, FDPart3tmp61));
        const REAL_SIMD_ARRAY FDPart3tmp94 = MulSIMD(FDPart3tmp51, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp96 = MulSIMD(FDPart3tmp36, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp98 = MulSIMD(FDPart3tmp51, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp106 = MulSIMD(FDPart3tmp47, FDPart3tmp76);
        const REAL_SIMD_ARRAY FDPart3tmp205 =
            MulSIMD(FDPart3tmp21, FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp46,
                                                  FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp140, aDD00),
                                                                  FusedMulAddSIMD(FDPart3tmp138, FDPart3tmp54, aDD_dD000))));
        const REAL_SIMD_ARRAY FDPart3tmp215 =
            FusedMulAddSIMD(FDPart3tmp8, aDD_dD111,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp165, FDPart3tmp28),
                                            FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp50, MulSIMD(FDPart3tmp166, FDPart3tmp46))));
        const REAL_SIMD_ARRAY FDPart3tmp222 =
            MulSIMD(FDPart3tmp21, FusedMulAddSIMD(FDPart3tmp178, FDPart3tmp54,
                                                  FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp176, FDPart3tmp43),
                                                                  FusedMulAddSIMD(FDPart3tmp10, aDD_dD222, MulSIMD(FDPart3tmp177, FDPart3tmp50)))));
        const REAL_SIMD_ARRAY FDPart3tmp228 =
            FusedMulAddSIMD(FDPart3tmp159, FDPart3tmp28, FusedMulAddSIMD(FDPart3tmp160, FDPart3tmp90, MulSIMD(FDPart3tmp158, FDPart3tmp85)));
        const REAL_SIMD_ARRAY FDPart3tmp229 =
            FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp87, MulSIMD(FDPart3tmp151, FDPart3tmp43)));
        const REAL_SIMD_ARRAY FDPart3tmp231 =
            FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp173, aDD00, MulSIMD(FDPart3tmp171, FDPart3tmp87)));
        const REAL_SIMD_ARRAY FDPart3tmp102 = MulSIMD(FDPart3tmp36, MulSIMD(FDPart3tmp43, FDPart3tmp51));
        const REAL_SIMD_ARRAY FDPart3tmp109 = MulSIMD(FDPart3tmp36, MulSIMD(FDPart3tmp38, FDPart3tmp76));
        const REAL_SIMD_ARRAY FDPart3tmp111 = MulSIMD(FDPart3tmp51, MulSIMD(FDPart3tmp76, FDPart3tmp87));
        const REAL_SIMD_ARRAY FDPart3tmp119 = MulSIMD(FDPart3tmp38, MulSIMD(FDPart3tmp61, FDPart3tmp76));
        const REAL_SIMD_ARRAY FDPart3tmp121 = MulSIMD(FDPart3tmp69, MulSIMD(FDPart3tmp76, FDPart3tmp90));
        const REAL_SIMD_ARRAY FDPart3tmp208 = FusedMulAddSIMD(
            FDPart3tmp152, FDPart3tmp46,
            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp153, aDD00), FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp54, aDD_dD001)));
        const REAL_SIMD_ARRAY FDPart3tmp210 = FusedMulAddSIMD(
            FDPart3tmp159, FDPart3tmp46,
            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp160, aDD00), FusedMulAddSIMD(FDPart3tmp158, FDPart3tmp54, aDD_dD002)));
        const REAL_SIMD_ARRAY FDPart3tmp212 =
            FusedMulAddSIMD(FDPart3tmp45, aDD11,
                            FusedMulAddSIMD(FDPart3tmp8, aDD_dD110,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp152, FDPart3tmp28),
                                                            FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp50, MulSIMD(FDPart3tmp153, FDPart3tmp46)))));
        const REAL_SIMD_ARRAY FDPart3tmp214 =
            FusedMulAddSIMD(FDPart3tmp8, aDD_dD112,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp172, FDPart3tmp28),
                                            FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp50, MulSIMD(FDPart3tmp173, FDPart3tmp46))));
        const REAL_SIMD_ARRAY FDPart3tmp217 =
            FusedMulAddSIMD(FDPart3tmp159, FDPart3tmp50,
                            FusedMulAddSIMD(FDPart3tmp160, FDPart3tmp54,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp158, FDPart3tmp43),
                                                            FusedMulAddSIMD(FDPart3tmp10, aDD_dD220, MulSIMD(FDPart3tmp155, aDD22)))));
        const REAL_SIMD_ARRAY FDPart3tmp220 =
            FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp50,
                            FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp54,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp171, FDPart3tmp43),
                                                            FusedMulAddSIMD(FDPart3tmp10, aDD_dD221, MulSIMD(FDPart3tmp167, aDD22)))));
        const REAL_SIMD_ARRAY FDPart3tmp223 = FusedMulAddSIMD(
            FDPart3tmp151, FDPart3tmp87,
            FusedMulAddSIMD(FDPart3tmp152, FDPart3tmp90,
                            FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp28,
                                            FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp90,
                                                            FusedMulAddSIMD(FDPart3tmp153, aDD00,
                                                                            FusedMulAddSIMD(aDD_dD010, f0_of_xx0,
                                                                                            FusedMulAddSIMD(FDPart3tmp138, FDPart3tmp85, aDD01)))))));
        const REAL_SIMD_ARRAY FDPart3tmp225 = FusedMulAddSIMD(
            FDPart3tmp158, FDPart3tmp87,
            FusedMulAddSIMD(
                FDPart3tmp159, FDPart3tmp90,
                FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp85,
                                FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp87,
                                                FusedMulAddSIMD(FDPart3tmp160, aDD00,
                                                                FusedMulAddSIMD(FDPart3tmp33, aDD_dD020,
                                                                                FusedMulAddSIMD(FDPart3tmp138, FDPart3tmp43, FDPart3tmp53)))))));
        const REAL_SIMD_ARRAY FDPart3tmp227 = FusedMulAddSIMD(
            FDPart3tmp165, FDPart3tmp90,
            FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp90,
                            FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp87,
                                            FusedMulAddSIMD(FDPart3tmp166, aDD00,
                                                            FusedMulAddSIMD(aDD_dD011, f0_of_xx0,
                                                                            FusedMulAddSIMD(FDPart3tmp151, FDPart3tmp85,
                                                                                            MulSIMD(FDPart3tmp152, FDPart3tmp28)))))));
        const REAL_SIMD_ARRAY FDPart3tmp230 =
            FusedMulAddSIMD(FDPart3tmp136, aDD12, FusedMulAddSIMD(FDPart3tmp22, aDD_dD120, AddSIMD(FDPart3tmp228, FDPart3tmp229)));
        const REAL_SIMD_ARRAY FDPart3tmp232 = FusedMulAddSIMD(
            FDPart3tmp33, aDD_dD021, FusedMulAddSIMD(aDD02, MulSIMD(f0_of_xx0, f1_of_xx1__D1), AddSIMD(FDPart3tmp229, FDPart3tmp231)));
        const REAL_SIMD_ARRAY FDPart3tmp233 = FusedMulAddSIMD(
            FDPart3tmp172, FDPart3tmp28,
            FusedMulAddSIMD(
                FDPart3tmp173, FDPart3tmp90,
                FusedMulAddSIMD(FDPart3tmp166, FDPart3tmp87,
                                FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp85,
                                                FusedMulAddSIMD(FDPart3tmp22, aDD_dD121,
                                                                FusedMulAddSIMD(FDPart3tmp8, MulSIMD(aDD12, f1_of_xx1__D1),
                                                                                FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp43,
                                                                                                MulSIMD(FDPart3tmp165, FDPart3tmp85))))))));
        const REAL_SIMD_ARRAY FDPart3tmp234 = AddSIMD(FDPart3tmp228, FusedMulAddSIMD(aDD_dD012, f0_of_xx0, FDPart3tmp231));
        const REAL_SIMD_ARRAY FDPart3tmp235 = FusedMulAddSIMD(
            FDPart3tmp177, FDPart3tmp90,
            FusedMulAddSIMD(FDPart3tmp160, FDPart3tmp87,
                            FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp87,
                                            FusedMulAddSIMD(FDPart3tmp178, aDD00,
                                                            FusedMulAddSIMD(FDPart3tmp33, aDD_dD022,
                                                                            FusedMulAddSIMD(FDPart3tmp158, FDPart3tmp43,
                                                                                            MulSIMD(FDPart3tmp159, FDPart3tmp85)))))));
        const REAL_SIMD_ARRAY FDPart3tmp237 = FusedMulAddSIMD(
            FDPart3tmp177, FDPart3tmp28,
            FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp87,
                            FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp85,
                                            FusedMulAddSIMD(FDPart3tmp178, FDPart3tmp90,
                                                            FusedMulAddSIMD(FDPart3tmp22, aDD_dD122,
                                                                            FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp43,
                                                                                            MulSIMD(FDPart3tmp172, FDPart3tmp85)))))));
        const REAL_SIMD_ARRAY FDPart3tmp203 = FusedMulAddSIMD(
            FDPart3_Integer_6, MulSIMD(FDPart3tmp85, FDPart3tmp98),
            FusedMulAddSIMD(
                FDPart3_Integer_6, MulSIMD(FDPart3tmp48, FDPart3tmp87),
                FusedMulAddSIMD(
                    FDPart3_Integer_6, MulSIMD(FDPart3tmp68, FDPart3tmp90),
                    FusedMulAddSIMD(
                        FDPart3_Integer_6, MulSIMD(FDPart3tmp28, FDPart3tmp70),
                        FusedMulAddSIMD(
                            FDPart3_Integer_6, MulSIMD(FDPart3tmp43, FDPart3tmp52),
                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp87, FDPart3tmp94),
                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp90, FDPart3tmp96),
                                                            FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp92,
                                                                            MulSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp27, FDPart3tmp85))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp204 = FusedMulAddSIMD(
            FDPart3_Integer_6, MulSIMD(FDPart3tmp72, FDPart3tmp87),
            FusedMulAddSIMD(
                FDPart3_Integer_6, MulSIMD(FDPart3tmp28, FDPart3tmp68),
                FusedMulAddSIMD(
                    FDPart3_Integer_6, MulSIMD(FDPart3tmp48, FDPart3tmp85),
                    FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp111,
                                    FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp106, FDPart3tmp90),
                                                    FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp80, FDPart3tmp90),
                                                                    FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp85, FDPart3tmp94),
                                                                                    FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp102,
                                                                                                    MulSIMD(FDPart3_Integer_6, FDPart3tmp109)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp224 = MulSIMD(FDPart3tmp223, FDPart3tmp76);
        const REAL_SIMD_ARRAY FDPart3tmp226 = MulSIMD(FDPart3tmp225, FDPart3tmp76);
        const REAL_SIMD_ARRAY FDPart3tmp241 = FusedMulAddSIMD(
            FDPart3_Integer_6, MulSIMD(FDPart3tmp74, FDPart3tmp90),
            FusedMulAddSIMD(
                FDPart3_Integer_6, MulSIMD(FDPart3tmp43, FDPart3tmp48),
                FusedMulAddSIMD(
                    FDPart3_Integer_6, MulSIMD(FDPart3tmp68, FDPart3tmp85),
                    FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp106, FDPart3tmp87),
                                    FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp28, FDPart3tmp71),
                                                    FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp80, FDPart3tmp87),
                                                                    FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp85, FDPart3tmp96),
                                                                                    FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp119,
                                                                                                    MulSIMD(FDPart3_Integer_6, FDPart3tmp121)))))))));
        const REAL_SIMD_ARRAY FDPart3tmp238 = FusedMulSubSIMD(
            FDPart3tmp122,
            FusedMulAddSIMD(
                FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp36, trK_dD0)),
                FusedMulAddSIMD(
                    FDPart3tmp131,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(cf_dD2, FusedMulAddSIMD(
                                                FDPart3tmp201, FDPart3tmp48,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp202, FDPart3tmp52,
                                                    FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp27, FDPart3tmp28),
                                                                    FusedMulAddSIMD(MulSIMD(FDPart3_Integer_12, FDPart3tmp36),
                                                                                    MulSIMD(FDPart3tmp51, FDPart3tmp87),
                                                                                    FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp39,
                                                                                                    MulSIMD(FDPart3_Integer_6, FDPart3tmp44)))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp26, trK_dD1)),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3_Integer_2, FDPart3tmp235), MulSIMD(FDPart3tmp36, FDPart3tmp51),
                            FusedMulAddSIMD(
                                FDPart3tmp131, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp204, cf_dD0)),
                                FusedMulAddSIMD(
                                    FDPart3tmp217, MulSIMD(FDPart3tmp36, FDPart3tmp51),
                                    FusedMulAddSIMD(
                                        FDPart3tmp220, MulSIMD(FDPart3tmp26, FDPart3tmp51),
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp237, FDPart3tmp52),
                                            FusedMulAddSIMD(
                                                FDPart3tmp205, MulSIMD(FDPart3tmp36, FDPart3tmp76),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp233, FDPart3tmp98,
                                                    FusedMulAddSIMD(
                                                        FDPart3_Integer_2, MulSIMD(FDPart3tmp234, FDPart3tmp48),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp232, FDPart3tmp94,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp233, FDPart3tmp27,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp230, FDPart3tmp94,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp232, FDPart3tmp48,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp227, FDPart3tmp96,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp230, FDPart3tmp48,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp226, FDPart3tmp51,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp227, FDPart3tmp68,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp224, FDPart3tmp47,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp225, FDPart3tmp72,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp222,
                                                                                                    MulSIMD(FDPart3tmp41, FDPart3tmp41),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp223, FDPart3tmp80,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp214, FDPart3tmp27,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp215, FDPart3tmp70,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp210, FDPart3tmp72,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp212, FDPart3tmp68,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp123,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_NegativeOne_,
                                                                                                                                    FDPart3_Rational_2_3),
                                                                                                                                MulSIMD(FDPart3tmp41,
                                                                                                                                        trK_dD2)),
                                                                                                                            FusedMulSubSIMD(
                                                                                                                                FDPart3tmp208,
                                                                                                                                FDPart3tmp80,
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Rational_1_2,
                                                                                                                                        FDPart3tmp131),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp203,
                                                                                                                                        cf_dD1))))))))))))))))))))))))))))))),
            MulSIMD(
                FDPart3_Integer_8,
                MulSIMD(PI,
                        FusedMulAddSIMD(
                            FDPart3tmp182,
                            MulSIMD(FDPart3tmp195,
                                    FusedMulAddSIMD(
                                        FDPart3tmp187, MulSIMD(f0_of_xx0, hDD01),
                                        FusedMulAddSIMD(
                                            FDPart3tmp188, MulSIMD(FDPart3tmp22, hDD12),
                                            FusedMulAddSIMD(
                                                FDPart3tmp185, FDPart3tmp9,
                                                MulSIMD(FDPart3tmp193,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp191, MulSIMD(f0_of_xx0, hDD12),
                                                            FusedMulSubSIMD(FDPart3tmp183, MulSIMD(FDPart3tmp59, vetU0),
                                                                            MulSIMD(FDPart3tmp189, MulSIMD(FDPart3tmp197, FDPart3tmp9))))))))),
                            FusedMulAddSIMD(
                                FDPart3tmp182,
                                MulSIMD(FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp179, MulSIMD(FDPart3tmp179, MulSIMD(FDPart3tmp5, FDPart3tmp9))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp193,
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp14, FDPart3tmp191), DivSIMD(FDPart3tmp197, f1_of_xx1),
                                                FusedMulSubSIMD(FDPart3tmp192, FDPart3tmp34, MulSIMD(FDPart3tmp189, MulSIMD(FDPart3tmp33, hDD12)))),
                                            FusedMulAddSIMD(FDPart3tmp185, MulSIMD(FDPart3tmp22, hDD12),
                                                            FusedMulAddSIMD(FDPart3tmp14, FDPart3tmp188, MulSIMD(FDPart3tmp187, FDPart3tmp34))))),
                                MulSIMD(FDPart3tmp180,
                                        MulSIMD(FDPart3tmp182,
                                                FusedMulAddSIMD(FDPart3tmp193,
                                                                FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp5,
                                                                                FusedMulSubSIMD(FDPart3tmp191, hDD02, MulSIMD(FDPart3tmp189, hDD01))),
                                                                FusedMulAddSIMD(FDPart3tmp185, MulSIMD(f0_of_xx0, hDD01),
                                                                                FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp5,
                                                                                                MulSIMD(FDPart3tmp188, FDPart3tmp34)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp243 = FusedMulSubSIMD(
            FDPart3tmp122,
            FusedMulAddSIMD(
                FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp61, trK_dD0)),
                FusedMulAddSIMD(
                    FDPart3tmp131,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(cf_dD1,
                                    FusedMulAddSIMD(
                                        FDPart3tmp202, FDPart3tmp70,
                                        FusedMulAddSIMD(FDPart3_Integer_12, MulSIMD(FDPart3tmp68, FDPart3tmp87),
                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp27, FDPart3tmp43),
                                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp28, FDPart3tmp66),
                                                                                        FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp63,
                                                                                                        MulSIMD(FDPart3tmp201, FDPart3tmp71)))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp26, trK_dD2)),
                        FusedMulAddSIMD(
                            FDPart3tmp222, MulSIMD(FDPart3tmp26, FDPart3tmp41),
                            FusedMulAddSIMD(
                                FDPart3tmp131, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp241, cf_dD0)),
                                FusedMulAddSIMD(
                                    FDPart3tmp205, MulSIMD(FDPart3tmp61, FDPart3tmp76),
                                    FusedMulAddSIMD(
                                        FDPart3tmp217, MulSIMD(FDPart3tmp36, FDPart3tmp47),
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp232, FDPart3tmp68),
                                            FusedMulAddSIMD(
                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp233, FDPart3tmp70),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp237, FDPart3tmp98,
                                                    FusedMulAddSIMD(
                                                        FDPart3_Integer_2, MulSIMD(FDPart3tmp227, FDPart3tmp71),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp235, FDPart3tmp94,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp237, FDPart3tmp27,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp234, FDPart3tmp96,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp235, FDPart3tmp48,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp230, FDPart3tmp96,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp234, FDPart3tmp68,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp226, FDPart3tmp47,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp230, FDPart3tmp68,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp224, FDPart3tmp69,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp225, FDPart3tmp80,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp220, FDPart3tmp27,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp223, FDPart3tmp74,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp214, FDPart3tmp70,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp215, FDPart3tmp66,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp210, FDPart3tmp80,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp212, FDPart3tmp71,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp123,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_NegativeOne_,
                                                                                                                                    FDPart3_Rational_2_3),
                                                                                                                                MulSIMD(FDPart3tmp65,
                                                                                                                                        trK_dD1)),
                                                                                                                            FusedMulSubSIMD(
                                                                                                                                FDPart3tmp208,
                                                                                                                                FDPart3tmp74,
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Rational_1_2,
                                                                                                                                        FDPart3tmp131),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp203,
                                                                                                                                        cf_dD2))))))))))))))))))))))))))))))),
            MulSIMD(
                FDPart3_Integer_8,
                MulSIMD(
                    PI,
                    FusedMulAddSIMD(
                        FDPart3tmp182,
                        MulSIMD(FDPart3tmp239,
                                FusedMulAddSIMD(
                                    FDPart3tmp193,
                                    FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp5, FusedMulSubSIMD(FDPart3tmp191, hDD02, MulSIMD(FDPart3tmp189, hDD01))),
                                    FusedMulAddSIMD(FDPart3tmp185, MulSIMD(f0_of_xx0, hDD01),
                                                    FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp5, MulSIMD(FDPart3tmp188, FDPart3tmp34))))),
                        FusedMulAddSIMD(
                            FDPart3tmp182,
                            MulSIMD(
                                FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp179, MulSIMD(FDPart3tmp14, MulSIMD(FDPart3tmp179, FDPart3tmp5))),
                                FusedMulAddSIMD(
                                    FDPart3tmp187, MulSIMD(f0_of_xx0, hDD01),
                                    FusedMulAddSIMD(
                                        FDPart3tmp188, MulSIMD(FDPart3tmp22, hDD12),
                                        FusedMulAddSIMD(FDPart3tmp185, FDPart3tmp9,
                                                        MulSIMD(FDPart3tmp193,
                                                                FusedMulAddSIMD(FDPart3tmp191, MulSIMD(f0_of_xx0, hDD12),
                                                                                FusedMulSubSIMD(FDPart3tmp183, MulSIMD(FDPart3tmp59, vetU0),
                                                                                                MulSIMD(FDPart3tmp189,
                                                                                                        MulSIMD(FDPart3tmp197, FDPart3tmp9))))))))),
                            MulSIMD(FDPart3tmp182,
                                    MulSIMD(FDPart3tmp195,
                                            FusedMulAddSIMD(FDPart3tmp193,
                                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp14, FDPart3tmp191), DivSIMD(FDPart3tmp197, f1_of_xx1),
                                                                            FusedMulSubSIMD(FDPart3tmp192, FDPart3tmp34,
                                                                                            MulSIMD(FDPart3tmp189, MulSIMD(FDPart3tmp33, hDD12)))),
                                                            FusedMulAddSIMD(FDPart3tmp185, MulSIMD(FDPart3tmp22, hDD12),
                                                                            FusedMulAddSIMD(FDPart3tmp14, FDPart3tmp188,
                                                                                            MulSIMD(FDPart3tmp187, FDPart3tmp34)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp245 = FusedMulSubSIMD(
            FDPart3tmp122,
            FusedMulAddSIMD(
                FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp61, trK_dD1)),
                FusedMulAddSIMD(
                    FDPart3tmp131,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(cf_dD0,
                                    FusedMulAddSIMD(
                                        FDPart3tmp202, FDPart3tmp80,
                                        FusedMulAddSIMD(FDPart3_Integer_12, MulSIMD(FDPart3tmp82, FDPart3tmp87),
                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp28, FDPart3tmp74),
                                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp43, FDPart3tmp72),
                                                                                        FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp78,
                                                                                                        MulSIMD(FDPart3tmp201, FDPart3tmp83)))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp123, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp36, trK_dD2)),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3tmp21, FDPart3tmp220), MulSIMD(FDPart3tmp26, FDPart3tmp36),
                            FusedMulAddSIMD(
                                FDPart3tmp131, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp241, cf_dD1)),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3_Integer_2, FDPart3tmp21), MulSIMD(FDPart3tmp224, FDPart3tmp61),
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3_Integer_2, FDPart3tmp21), MulSIMD(FDPart3tmp226, FDPart3tmp36),
                                        FusedMulAddSIMD(
                                            FDPart3tmp227, MulSIMD(FDPart3tmp69, FDPart3tmp76),
                                            FusedMulAddSIMD(
                                                FDPart3tmp235, MulSIMD(FDPart3tmp51, FDPart3tmp76),
                                                FusedMulAddSIMD(
                                                    FDPart3_Integer_2, MulSIMD(FDPart3tmp230, FDPart3tmp80),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp222, MulSIMD(FDPart3tmp36, FDPart3tmp41),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp237, FDPart3tmp48,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp237, FDPart3tmp94,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp234, FDPart3tmp80,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp235, FDPart3tmp72,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp233, FDPart3tmp68,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp233, FDPart3tmp96,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp227, FDPart3tmp74,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp232, FDPart3tmp80,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp215, FDPart3tmp71,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp217, FDPart3tmp72,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp212, FDPart3tmp74,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp214, FDPart3tmp68,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp208, FDPart3tmp83,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp210, FDPart3tmp82,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp106, FDPart3tmp234,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp205,
                                                                                                                        MulSIMD(FDPart3tmp76,
                                                                                                                                FDPart3tmp76),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp123,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_NegativeOne_,
                                                                                                                                    FDPart3_Rational_2_3),
                                                                                                                                MulSIMD(FDPart3tmp76,
                                                                                                                                        trK_dD0)),
                                                                                                                            FusedMulSubSIMD(
                                                                                                                                FDPart3tmp106,
                                                                                                                                FDPart3tmp232,
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Rational_1_2,
                                                                                                                                        FDPart3tmp131),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp204,
                                                                                                                                        cf_dD2))))))))))))))))))))))))))))))),
            MulSIMD(FDPart3_Integer_8,
                    MulSIMD(PI, FusedMulAddSIMD(
                                    FDPart3tmp182,
                                    MulSIMD(FDPart3tmp239,
                                            FusedMulAddSIMD(
                                                FDPart3tmp187, MulSIMD(f0_of_xx0, hDD01),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp188, MulSIMD(FDPart3tmp22, hDD12),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp185, FDPart3tmp9,
                                                        MulSIMD(FDPart3tmp193,
                                                                FusedMulAddSIMD(FDPart3tmp191, MulSIMD(f0_of_xx0, hDD12),
                                                                                FusedMulSubSIMD(FDPart3tmp183, MulSIMD(FDPart3tmp59, vetU0),
                                                                                                MulSIMD(FDPart3tmp189,
                                                                                                        MulSIMD(FDPart3tmp197, FDPart3tmp9))))))))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp182,
                                        MulSIMD(FusedMulAddSIMD(FDPart3tmp179, FDPart3tmp18, MulSIMD(FDPart3tmp179, FDPart3tmp6)),
                                                FusedMulAddSIMD(FDPart3tmp193,
                                                                FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp5,
                                                                                FusedMulSubSIMD(FDPart3tmp191, hDD02, MulSIMD(FDPart3tmp189, hDD01))),
                                                                FusedMulAddSIMD(FDPart3tmp185, MulSIMD(f0_of_xx0, hDD01),
                                                                                FusedMulAddSIMD(FDPart3tmp187, FDPart3tmp5,
                                                                                                MulSIMD(FDPart3tmp188, FDPart3tmp34))))),
                                        MulSIMD(FDPart3tmp180,
                                                MulSIMD(FDPart3tmp182,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp193,
                                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp14, FDPart3tmp191), DivSIMD(FDPart3tmp197, f1_of_xx1),
                                                                            FusedMulSubSIMD(FDPart3tmp192, FDPart3tmp34,
                                                                                            MulSIMD(FDPart3tmp189, MulSIMD(FDPart3tmp33, hDD12)))),
                                                            FusedMulAddSIMD(FDPart3tmp185, MulSIMD(FDPart3tmp22, hDD12),
                                                                            FusedMulAddSIMD(FDPart3tmp14, FDPart3tmp188,
                                                                                            MulSIMD(FDPart3tmp187, FDPart3tmp34)))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(
            aDD01,
            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                    MulSIMD(f0_of_xx0,
                            FusedMulAddSIMD(
                                FDPart3tmp74, FDPart3tmp90,
                                FusedMulAddSIMD(
                                    FDPart3tmp43, FDPart3tmp48,
                                    FusedMulAddSIMD(FDPart3tmp68, FDPart3tmp85,
                                                    FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp87,
                                                                    FusedMulAddSIMD(FDPart3tmp28, FDPart3tmp71,
                                                                                    FusedMulAddSIMD(FDPart3tmp80, FDPart3tmp87,
                                                                                                    FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp96,
                                                                                                                    AddSIMD(FDPart3tmp119,
                                                                                                                            FDPart3tmp121)))))))))),
            FusedMulAddSIMD(
                MulSIMD(FDPart3tmp8, aDD12),
                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                        MulSIMD(f1_of_xx1, FusedMulAddSIMD(
                                               FDPart3tmp85, FDPart3tmp98,
                                               FusedMulAddSIMD(
                                                   FDPart3tmp48, FDPart3tmp87,
                                                   FusedMulAddSIMD(
                                                       FDPart3tmp68, FDPart3tmp90,
                                                       FusedMulAddSIMD(
                                                           FDPart3tmp28, FDPart3tmp70,
                                                           FusedMulAddSIMD(FDPart3tmp43, FDPart3tmp52,
                                                                           FusedMulAddSIMD(FDPart3tmp87, FDPart3tmp94,
                                                                                           FusedMulAddSIMD(FDPart3tmp90, FDPart3tmp96,
                                                                                                           FusedMulAddSIMD(FDPart3tmp27, FDPart3tmp85,
                                                                                                                           FDPart3tmp92)))))))))),
                NegFusedMulAddSIMD(
                    FDPart3tmp8,
                    MulSIMD(aDD11, FusedMulAddSIMD(
                                       FDPart3tmp28, FDPart3tmp66,
                                       FusedMulAddSIMD(FDPart3tmp46, FDPart3tmp71,
                                                       FusedMulAddSIMD(FDPart3tmp50, FDPart3tmp70,
                                                                       FusedMulAddSIMD(FDPart3tmp54, FDPart3tmp68,
                                                                                       FusedMulAddSIMD(FDPart3tmp27, FDPart3tmp43, FDPart3tmp63)))))),
                    FusedMulAddSIMD(
                        T4UU00, MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, MulSIMD(alpha, alpha))),
                        FusedMulAddSIMD(
                            FDPart3tmp122,
                            FusedMulAddSIMD(
                                FDPart3tmp36,
                                MulSIMD(MulSIMD(FDPart3tmp123, FDPart3tmp125),
                                        MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(cf_dD0, cf_dD2))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp123, FDPart3tmp125),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp76, MulSIMD(cf_dD0, cf_dD0))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp26,
                                        MulSIMD(MulSIMD(FDPart3tmp123, FDPart3tmp125),
                                                MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(cf_dD1, cf_dD2))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3tmp123, FDPart3tmp125),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp41, MulSIMD(cf_dD2, cf_dD2))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp123, FDPart3tmp125),
                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                        MulSIMD(FDPart3tmp65, MulSIMD(cf_dD1, cf_dD1))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp123,
                                                    MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                                            MulSIMD(FDPart3tmp65,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp131,
                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                MulSIMD(FDPart3tmp165, cf_dD1)),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp131,
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                    MulSIMD(FDPart3tmp166, cf_dD0)),
                                                                            FusedMulSubSIMD(
                                                                                FDPart3tmp132,
                                                                                FusedMulSubSIMD(FDPart3tmp131, MulSIMD(cf_dD1, cf_dD1), cf_dDD11),
                                                                                MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp131),
                                                                                        MulSIMD(FDPart3tmp164, cf_dD2))))))),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp123,
                                                        MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                                                MulSIMD(FDPart3tmp76,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp131,
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                    MulSIMD(FDPart3tmp139, cf_dD1)),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp131,
                                                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                        MulSIMD(FDPart3tmp140, cf_dD0)),
                                                                                FusedMulSubSIMD(
                                                                                    FDPart3tmp132,
                                                                                    FusedMulSubSIMD(FDPart3tmp131, MulSIMD(cf_dD0, cf_dD0), cf_dDD00),
                                                                                    MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp131),
                                                                                            MulSIMD(FDPart3tmp138, cf_dD2))))))),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp123,
                                                            MulSIMD(
                                                                MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_),
                                                                MulSIMD(FDPart3tmp61,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp131,
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                    MulSIMD(FDPart3tmp152, cf_dD1)),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp131,
                                                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                        MulSIMD(FDPart3tmp153, cf_dD0)),
                                                                                FusedMulSubSIMD(
                                                                                    FDPart3tmp132,
                                                                                    FusedMulSubSIMD(FDPart3tmp131, MulSIMD(cf_dD0, cf_dD1), cf_dDD01),
                                                                                    MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp131),
                                                                                            MulSIMD(FDPart3tmp151, cf_dD2))))))),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp123,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                                                        MulSIMD(FDPart3tmp41,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp131,
                                                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                            MulSIMD(FDPart3tmp177, cf_dD1)),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp131,
                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                MulSIMD(FDPart3tmp178, cf_dD0)),
                                                                                        FusedMulSubSIMD(
                                                                                            FDPart3tmp132,
                                                                                            FusedMulSubSIMD(FDPart3tmp131, MulSIMD(cf_dD2, cf_dD2),
                                                                                                            cf_dDD22),
                                                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp131),
                                                                                                    MulSIMD(FDPart3tmp176, cf_dD2))))))),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp124, MulSIMD(FDPart3tmp61, RbarDD01),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp123,
                                                                        MulSIMD(
                                                                            MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_),
                                                                            MulSIMD(
                                                                                FDPart3tmp36,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp131,
                                                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                            MulSIMD(FDPart3tmp159, cf_dD1)),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp131,
                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                MulSIMD(FDPart3tmp160, cf_dD0)),
                                                                                        FusedMulSubSIMD(
                                                                                            FDPart3tmp132,
                                                                                            FusedMulSubSIMD(FDPart3tmp131, MulSIMD(cf_dD0, cf_dD2),
                                                                                                            cf_dDD02),
                                                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp131),
                                                                                                    MulSIMD(FDPart3tmp158, cf_dD2))))))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp124, MulSIMD(FDPart3tmp26, RbarDD12),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp124, MulSIMD(FDPart3tmp36, RbarDD02),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp123, MulSIMD(FDPart3tmp65, RbarDD11),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp123, MulSIMD(FDPart3tmp76, RbarDD00),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp61,
                                                                                            MulSIMD(MulSIMD(FDPart3tmp123, FDPart3tmp125),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_4,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(cf_dD0, cf_dD1))),
                                                                                            FusedMulSubSIMD(
                                                                                                FDPart3tmp123, MulSIMD(FDPart3tmp41, RbarDD22),
                                                                                                MulSIMD(
                                                                                                    MulSIMD(FDPart3_Integer_16, FDPart3tmp123),
                                                                                                    MulSIMD(
                                                                                                        FDPart3tmp26,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp131,
                                                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_,
                                                                                                                            FDPart3_Rational_1_2),
                                                                                                                    MulSIMD(FDPart3tmp172, cf_dD1)),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp131,
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3_NegativeOne_,
                                                                                                                            FDPart3_Rational_1_2),
                                                                                                                    MulSIMD(FDPart3tmp173, cf_dD0)),
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3tmp132,
                                                                                                                    FusedMulSubSIMD(
                                                                                                                        FDPart3tmp131,
                                                                                                                        MulSIMD(cf_dD1, cf_dD2),
                                                                                                                        cf_dDD12),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Rational_1_2,
                                                                                                                                FDPart3tmp131),
                                                                                                                        MulSIMD(
                                                                                                                            FDPart3tmp171,
                                                                                                                            cf_dD2)))))))))))))))))))))))),
                            NegFusedMulAddSIMD(
                                FDPart3tmp10,
                                MulSIMD(aDD22, FusedMulAddSIMD(
                                                   FDPart3tmp27, FDPart3tmp28,
                                                   FusedMulAddSIMD(FDPart3tmp46, FDPart3tmp48,
                                                                   FusedMulAddSIMD(FDPart3tmp50, FDPart3tmp52,
                                                                                   FusedMulAddSIMD(FDPart3tmp36, MulSIMD(FDPart3tmp51, FDPart3tmp54),
                                                                                                   AddSIMD(FDPart3tmp39, FDPart3tmp44)))))),
                                FusedMulAddSIMD(
                                    MulSIMD(aDD02, f0_of_xx0),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                            MulSIMD(f1_of_xx1,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp72, FDPart3tmp87,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp28, FDPart3tmp68,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp48, FDPart3tmp85,
                                                                AddSIMD(FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp90, FDPart3tmp111),
                                                                        FusedMulAddSIMD(FDPart3tmp80, FDPart3tmp90,
                                                                                        FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp94,
                                                                                                        AddSIMD(FDPart3tmp102, FDPart3tmp109))))))))),
                                    FusedMulSubSIMD(
                                        FDPart3_Rational_2_3, MulSIMD(trK, trK),
                                        MulSIMD(aDD00, FusedMulAddSIMD(
                                                           FDPart3tmp43, FDPart3tmp72,
                                                           FusedMulAddSIMD(FDPart3tmp46, FDPart3tmp83,
                                                                           FusedMulAddSIMD(FDPart3tmp50, FDPart3tmp80,
                                                                                           FusedMulAddSIMD(FDPart3tmp54, FDPart3tmp82,
                                                                                                           FusedMulAddSIMD(FDPart3tmp28, FDPart3tmp74,
                                                                                                                           FDPart3tmp78))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(
            FDPart3tmp9, MulSIMD(FDPart3tmp243, FDPart3tmp243),
            FusedMulAddSIMD(MulSIMD(FDPart3tmp238, FDPart3tmp243), MulSIMD(FDPart3tmp49, hDD12),
                            FusedMulAddSIMD(MulSIMD(FDPart3tmp243, FDPart3tmp245), MulSIMD(FDPart3tmp45, hDD01),
                                            FusedMulAddSIMD(FDPart3tmp45, MulSIMD(MulSIMD(FDPart3tmp238, FDPart3tmp245), MulSIMD(f1_of_xx1, hDD02)),
                                                            FusedMulAddSIMD(FDPart3tmp14, MulSIMD(FDPart3tmp238, FDPart3tmp238),
                                                                            MulSIMD(FDPart3tmp5, MulSIMD(FDPart3tmp245, FDPart3tmp245)))))));

        WriteSIMD(&diagnostic_output_gfs[IDX4(HGF, i0, i1, i2)], __RHS_exp_0);
        WriteSIMD(&diagnostic_output_gfs[IDX4(MSQUAREDGF, i0, i1, i2)], __RHS_exp_1);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += simd_width)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
