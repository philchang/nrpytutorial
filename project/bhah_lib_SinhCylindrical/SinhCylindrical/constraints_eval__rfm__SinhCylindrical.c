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
 * Evaluate BSSN constraints.
 */
void constraints_eval__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                            const rfm_struct *restrict rfmstruct, const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs,
                                            REAL *restrict diagnostic_output_gfs) {
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
        static const double dblFDPart3_Integer_1 = 1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        static const double dblFDPart3_Integer_12 = 12.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_12 = ConstSIMD(dblFDPart3_Integer_12);

        static const double dblFDPart3_Integer_16 = 16.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_16 = ConstSIMD(dblFDPart3_Integer_16);

        static const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        static const double dblFDPart3_Integer_4 = 4.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_4 = ConstSIMD(dblFDPart3_Integer_4);

        static const double dblFDPart3_Integer_6 = 6.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_6 = ConstSIMD(dblFDPart3_Integer_6);

        static const double dblFDPart3_Integer_8 = 8.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_8 = ConstSIMD(dblFDPart3_Integer_8);

        static const double dblFDPart3_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        static const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        static const double dblFDPart3_Rational_2_3 = 2.0 / 3.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_2_3 = ConstSIMD(dblFDPart3_Rational_2_3);

        const REAL_SIMD_ARRAY FDPart3tmp1 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp3 = MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp5 = MulSIMD(f3_of_xx2, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp23 = MulSIMD(f0_of_xx0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp32 = MulSIMD(f0_of_xx0__D0, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp47 = MulSIMD(f0_of_xx0, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp120 = MulSIMD(cf, cf);
        const REAL_SIMD_ARRAY FDPart3tmp128 = DivSIMD(FDPart3_Integer_1, cf);
        const REAL_SIMD_ARRAY FDPart3tmp156 = MulSIMD(FDPart3_Integer_2, f3_of_xx2__D2);
        const REAL_SIMD_ARRAY FDPart3tmp176 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf));
        const REAL_SIMD_ARRAY FDPart3tmp178 = DivSIMD(FDPart3_Integer_1, MulSIMD(MulSIMD(MulSIMD(MulSIMD(MulSIMD(cf, cf), cf), cf), cf), cf));
        const REAL_SIMD_ARRAY FDPart3tmp180 = DivSIMD(FDPart3_Integer_1, MulSIMD(cf, cf));
        const REAL_SIMD_ARRAY FDPart3tmp190 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(T4UU00, alpha));
        const REAL_SIMD_ARRAY FDPart3tmp2 = MulSIMD(FDPart3tmp1, aDD11);
        const REAL_SIMD_ARRAY FDPart3tmp6 = MulSIMD(FDPart3tmp5, MulSIMD(hDD02, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp8 = FusedMulAddSIMD(FDPart3tmp5, hDD22, FDPart3tmp5);
        const REAL_SIMD_ARRAY FDPart3tmp12 = FusedMulAddSIMD(FDPart3tmp3, hDD00, FDPart3tmp3);
        const REAL_SIMD_ARRAY FDPart3tmp16 = FusedMulAddSIMD(FDPart3tmp1, hDD11, FDPart3tmp1);
        const REAL_SIMD_ARRAY FDPart3tmp26 = MulSIMD(FDPart3tmp23, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp31 = MulSIMD(FDPart3tmp5, aDD22);
        const REAL_SIMD_ARRAY FDPart3tmp33 = MulSIMD(FDPart3tmp32, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp36 = MulSIMD(FDPart3tmp32, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp41 = MulSIMD(FDPart3tmp3, aDD00);
        const REAL_SIMD_ARRAY FDPart3tmp48 = MulSIMD(FDPart3_Integer_2, FDPart3tmp47);
        const REAL_SIMD_ARRAY FDPart3tmp52 = MulSIMD(FDPart3_Integer_2, FDPart3tmp23);
        const REAL_SIMD_ARRAY FDPart3tmp55 = MulSIMD(FDPart3_Integer_2, FDPart3tmp32);
        const REAL_SIMD_ARRAY FDPart3tmp58 = MulSIMD(MulSIMD(FDPart3tmp3, FDPart3tmp47), MulSIMD(hDD01, hDD02));
        const REAL_SIMD_ARRAY FDPart3tmp85 = MulSIMD(FDPart3tmp23, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp87 = MulSIMD(FDPart3tmp32, aDD02);
        const REAL_SIMD_ARRAY FDPart3tmp89 = MulSIMD(FDPart3tmp47, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp123 = DivSIMD(FDPart3_Integer_1, FDPart3tmp120);
        const REAL_SIMD_ARRAY FDPart3tmp129 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp128);
        const REAL_SIMD_ARRAY FDPart3tmp131 = MulSIMD(FDPart3tmp23, hDD_dD012);
        const REAL_SIMD_ARRAY FDPart3tmp132 = MulSIMD(FDPart3tmp32, hDD_dD021);
        const REAL_SIMD_ARRAY FDPart3tmp157 = MulSIMD(FDPart3tmp156, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp169 = FusedMulAddSIMD(f0_of_xx0, f0_of_xx0__DD00, FDPart3tmp3);
        const REAL_SIMD_ARRAY FDPart3tmp171 = MulSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0__D0, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp181 = MulSIMD(FDPart3tmp180, alpha);
        const REAL_SIMD_ARRAY FDPart3tmp7 = MulSIMD(FDPart3tmp3, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp1), MulSIMD(FDPart3tmp6, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp10 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp1), MulSIMD(FDPart3tmp3, MulSIMD(hDD01, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp14 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp1), MulSIMD(FDPart3tmp5, MulSIMD(hDD12, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp17 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp3), MulSIMD(FDPart3tmp5, MulSIMD(hDD02, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp19 = MulSIMD(FDPart3tmp16, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp34 = MulSIMD(FDPart3tmp1, MulSIMD(FDPart3tmp33, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp49 = MulSIMD(FDPart3tmp48, aDD12);
        const REAL_SIMD_ARRAY FDPart3tmp53 = MulSIMD(FDPart3tmp52, aDD01);
        const REAL_SIMD_ARRAY FDPart3tmp56 = MulSIMD(FDPart3tmp55, aDD02);
        const REAL_SIMD_ARRAY FDPart3tmp60 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp12), MulSIMD(FDPart3tmp47, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp134 = FusedMulAddSIMD(FDPart3tmp47, hDD_dD120, FDPart3tmp33);
        const REAL_SIMD_ARRAY FDPart3tmp140 = AddSIMD(FDPart3tmp132, SubSIMD(SubSIMD(FDPart3tmp131, FDPart3tmp33), MulSIMD(FDPart3tmp47, hDD_dD120)));
        const REAL_SIMD_ARRAY FDPart3tmp145 = FusedMulAddSIMD(FDPart3tmp1, hDD_dD110, FusedMulAddSIMD(FDPart3tmp52, hDD11, FDPart3tmp52));
        const REAL_SIMD_ARRAY FDPart3tmp151 = FusedMulSubSIMD(FDPart3tmp48, hDD_dD121, MulSIMD(FDPart3tmp1, hDD_dD112));
        const REAL_SIMD_ARRAY FDPart3tmp152 =
            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, f0_of_xx0__D0),
                               FusedMulAddSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD11)),
                                               FusedMulSubSIMD(FDPart3tmp52, hDD_dD011, MulSIMD(FDPart3tmp1, hDD_dD110))));
        const REAL_SIMD_ARRAY FDPart3tmp158 = FusedMulAddSIMD(FDPart3tmp157, hDD22, FusedMulAddSIMD(FDPart3tmp5, hDD_dD222, FDPart3tmp157));
        const REAL_SIMD_ARRAY FDPart3tmp160 =
            FusedMulAddSIMD(FDPart3tmp156, MulSIMD(f0_of_xx0, hDD12), FusedMulSubSIMD(FDPart3tmp48, hDD_dD122, MulSIMD(FDPart3tmp5, hDD_dD221)));
        const REAL_SIMD_ARRAY FDPart3tmp162 =
            FusedMulAddSIMD(FDPart3tmp156, MulSIMD(f0_of_xx0__D0, hDD02), FusedMulSubSIMD(FDPart3tmp55, hDD_dD022, MulSIMD(FDPart3tmp5, hDD_dD220)));
        const REAL_SIMD_ARRAY FDPart3tmp168 = FusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, f0_of_xx0__DD00), MulSIMD(f3_of_xx2, hDD02),
                                                              FusedMulSubSIMD(FDPart3tmp55, hDD_dD020, MulSIMD(FDPart3tmp3, hDD_dD002)));
        const REAL_SIMD_ARRAY FDPart3tmp170 = FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp169, hDD01),
                                                              FusedMulSubSIMD(FDPart3tmp52, hDD_dD010, MulSIMD(FDPart3tmp3, hDD_dD001)));
        const REAL_SIMD_ARRAY FDPart3tmp172 = FusedMulAddSIMD(FDPart3tmp171, hDD00, FusedMulAddSIMD(FDPart3tmp3, hDD_dD000, FDPart3tmp171));
        const REAL_SIMD_ARRAY FDPart3tmp182 = MulSIMD(FDPart3tmp181, T4UU01);
        const REAL_SIMD_ARRAY FDPart3tmp183 = MulSIMD(FDPart3tmp181, T4UU03);
        const REAL_SIMD_ARRAY FDPart3tmp185 = MulSIMD(FDPart3tmp181, T4UU02);
        const REAL_SIMD_ARRAY FDPart3tmp187 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp180, vetU0));
        const REAL_SIMD_ARRAY FDPart3tmp188 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp180, vetU2));
        const REAL_SIMD_ARRAY FDPart3tmp196 = MulSIMD(FDPart3_Integer_12, FDPart3tmp89);
        const REAL_SIMD_ARRAY FDPart3tmp197 = MulSIMD(FDPart3_Integer_12, FDPart3tmp85);
        const REAL_SIMD_ARRAY FDPart3tmp198 = MulSIMD(FDPart3_Integer_12, FDPart3tmp87);
        const REAL_SIMD_ARRAY FDPart3tmp27 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp26, FDPart3tmp8));
        const REAL_SIMD_ARRAY FDPart3tmp37 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp16, FDPart3tmp36));
        const REAL_SIMD_ARRAY FDPart3tmp42 = AddSIMD(FDPart3tmp14, FDPart3tmp19);
        const REAL_SIMD_ARRAY FDPart3tmp61 = AddSIMD(FDPart3tmp58, FDPart3tmp60);
        const REAL_SIMD_ARRAY FDPart3tmp65 = FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp8, FDPart3tmp17);
        const REAL_SIMD_ARRAY FDPart3tmp76 = FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp16, FDPart3tmp10);
        const REAL_SIMD_ARRAY FDPart3tmp135 = AddSIMD(FDPart3tmp134, SubSIMD(FDPart3tmp131, FDPart3tmp132));
        const REAL_SIMD_ARRAY FDPart3tmp146 = AddSIMD(FDPart3tmp134, SubSIMD(FDPart3tmp132, FDPart3tmp131));
        const REAL_SIMD_ARRAY FDPart3tmp232 = FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp58, MulSIMD(FDPart3tmp176, FDPart3tmp60));
        const REAL_SIMD_ARRAY FDPart3tmp21 =
            FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp14,
                            FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp19,
                                            FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp17, FusedMulAddSIMD(FDPart3tmp10, FDPart3tmp8, FDPart3tmp7))));
        const REAL_SIMD_ARRAY FDPart3tmp28 = FusedMulAddSIMD(FDPart3tmp23, FDPart3tmp6, FDPart3tmp27);
        const REAL_SIMD_ARRAY FDPart3tmp38 = AddSIMD(FDPart3tmp34, FDPart3tmp37);
        const REAL_SIMD_ARRAY FDPart3tmp177 =
            FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp27, MulSIMD(FDPart3tmp176, MulSIMD(FDPart3tmp23, FDPart3tmp6)));
        const REAL_SIMD_ARRAY FDPart3tmp179 =
            DivSIMD(FDPart3_Integer_1,
                    FusedMulAddSIMD(FDPart3tmp12, MulSIMD(FDPart3tmp14, FDPart3tmp178),
                                    FusedMulAddSIMD(FDPart3tmp12, MulSIMD(FDPart3tmp178, FDPart3tmp19),
                                                    FusedMulAddSIMD(FDPart3tmp16, MulSIMD(FDPart3tmp17, FDPart3tmp178),
                                                                    FusedMulAddSIMD(FDPart3tmp178, FDPart3tmp7,
                                                                                    MulSIMD(FDPart3tmp10, MulSIMD(FDPart3tmp178, FDPart3tmp8)))))));
        const REAL_SIMD_ARRAY FDPart3tmp192 = FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp34, MulSIMD(FDPart3tmp176, FDPart3tmp37));
        const REAL_SIMD_ARRAY FDPart3tmp22 = DivSIMD(FDPart3_Integer_1, MulSIMD(FDPart3tmp21, FDPart3tmp21));
        const REAL_SIMD_ARRAY FDPart3tmp121 = DivSIMD(FDPart3_Integer_1, FDPart3tmp21);
        const REAL_SIMD_ARRAY FDPart3tmp29 = MulSIMD(FDPart3tmp22, MulSIMD(FDPart3tmp28, FDPart3tmp28));
        const REAL_SIMD_ARRAY FDPart3tmp39 = MulSIMD(FDPart3tmp22, MulSIMD(FDPart3tmp38, FDPart3tmp38));
        const REAL_SIMD_ARRAY FDPart3tmp43 = MulSIMD(FDPart3tmp22, MulSIMD(FDPart3tmp42, FDPart3tmp42));
        const REAL_SIMD_ARRAY FDPart3tmp50 = MulSIMD(FDPart3tmp22, FDPart3tmp42);
        const REAL_SIMD_ARRAY FDPart3tmp62 = MulSIMD(FDPart3tmp22, MulSIMD(FDPart3tmp61, FDPart3tmp61));
        const REAL_SIMD_ARRAY FDPart3tmp66 = MulSIMD(FDPart3tmp22, MulSIMD(FDPart3tmp65, FDPart3tmp65));
        const REAL_SIMD_ARRAY FDPart3tmp69 = MulSIMD(FDPart3tmp22, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp77 = MulSIMD(FDPart3tmp22, MulSIMD(FDPart3tmp76, FDPart3tmp76));
        const REAL_SIMD_ARRAY FDPart3tmp79 = MulSIMD(FDPart3tmp22, FDPart3tmp38);
        const REAL_SIMD_ARRAY FDPart3tmp81 = MulSIMD(FDPart3tmp22, FDPart3tmp76);
        const REAL_SIMD_ARRAY FDPart3tmp122 = MulSIMD(FDPart3_Integer_2, FDPart3tmp121);
        const REAL_SIMD_ARRAY FDPart3tmp136 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp38), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp5, hDD_dD220)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp121, FDPart3tmp3),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp42, hDD_dD002)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp135, FDPart3tmp28))));
        const REAL_SIMD_ARRAY FDPart3tmp137 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp61, hDD_dD220)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp121, FDPart3tmp28),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp3, hDD_dD002)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp135, FDPart3tmp65))));
        const REAL_SIMD_ARRAY FDPart3tmp138 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp76, hDD_dD220)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp121, FDPart3tmp3),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp38, hDD_dD002)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp135, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp141 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp38), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp5, hDD_dD221)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp1, FDPart3tmp121),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp28, hDD_dD112)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp140, FDPart3tmp42))));
        const REAL_SIMD_ARRAY FDPart3tmp142 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp61, hDD_dD221)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp1, FDPart3tmp121),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp65, hDD_dD112)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp140, FDPart3tmp28))));
        const REAL_SIMD_ARRAY FDPart3tmp143 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp5), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp76, hDD_dD221)),
            FusedMulSubSIMD(MulSIMD(FDPart3tmp1, FDPart3tmp121),
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp61, hDD_dD112)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp140, FDPart3tmp38))));
        const REAL_SIMD_ARRAY FDPart3tmp147 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp3), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp42, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp146, FDPart3tmp38)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp145, FDPart3tmp28))));
        const REAL_SIMD_ARRAY FDPart3tmp148 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp28), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp3, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp146, FDPart3tmp61)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp145, FDPart3tmp65))));
        const REAL_SIMD_ARRAY FDPart3tmp149 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp121, FDPart3tmp3), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp38, hDD_dD001)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp146, FDPart3tmp76)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp145, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp153 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp1, FDPart3tmp121), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp28, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp152, FDPart3tmp42)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp151, FDPart3tmp38))));
        const REAL_SIMD_ARRAY FDPart3tmp154 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp1, FDPart3tmp121), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp65, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp152, FDPart3tmp28)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp151, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp155 = FusedMulAddSIMD(
            MulSIMD(FDPart3tmp1, FDPart3tmp121), MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp61, hDD_dD111)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp152, FDPart3tmp38)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp151, FDPart3tmp76))));
        const REAL_SIMD_ARRAY FDPart3tmp163 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp162, FDPart3tmp42)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp160, FDPart3tmp28)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp158, FDPart3tmp38))));
        const REAL_SIMD_ARRAY FDPart3tmp164 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp162, FDPart3tmp28)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp160, FDPart3tmp65)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp158, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp165 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp162, FDPart3tmp38)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp160, FDPart3tmp61)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp158, FDPart3tmp76))));
        const REAL_SIMD_ARRAY FDPart3tmp173 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp172, FDPart3tmp42)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp170, FDPart3tmp28)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp168, FDPart3tmp38))));
        const REAL_SIMD_ARRAY FDPart3tmp174 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp172, FDPart3tmp28)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp170, FDPart3tmp65)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp168, FDPart3tmp61))));
        const REAL_SIMD_ARRAY FDPart3tmp175 = FusedMulAddSIMD(
            FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp172, FDPart3tmp38)),
            FusedMulSubSIMD(FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp170, FDPart3tmp61)),
                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp121), MulSIMD(FDPart3tmp168, FDPart3tmp76))));
        const REAL_SIMD_ARRAY FDPart3tmp46 = MulSIMD(FDPart3tmp22, MulSIMD(FDPart3tmp28, FDPart3tmp38));
        const REAL_SIMD_ARRAY FDPart3tmp51 = MulSIMD(FDPart3tmp28, FDPart3tmp50);
        const REAL_SIMD_ARRAY FDPart3tmp54 = MulSIMD(FDPart3tmp38, FDPart3tmp50);
        const REAL_SIMD_ARRAY FDPart3tmp68 = MulSIMD(FDPart3tmp22, MulSIMD(FDPart3tmp28, FDPart3tmp61));
        const REAL_SIMD_ARRAY FDPart3tmp70 = MulSIMD(FDPart3tmp28, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp71 = MulSIMD(FDPart3tmp61, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp80 = MulSIMD(FDPart3tmp61, FDPart3tmp79);
        const REAL_SIMD_ARRAY FDPart3tmp82 = MulSIMD(FDPart3tmp38, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp83 = MulSIMD(FDPart3tmp61, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp92 = MulSIMD(FDPart3tmp50, FDPart3tmp61);
        const REAL_SIMD_ARRAY FDPart3tmp94 = MulSIMD(FDPart3tmp38, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp98 = MulSIMD(FDPart3tmp50, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp105 = MulSIMD(FDPart3tmp28, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp109 = MulSIMD(FDPart3tmp42, FDPart3tmp81);
        const REAL_SIMD_ARRAY FDPart3tmp119 = MulSIMD(FDPart3tmp65, MulSIMD(FDPart3tmp81, FDPart3tmp89));
        const REAL_SIMD_ARRAY FDPart3tmp214 =
            FusedMulAddSIMD(FDPart3tmp155, FDPart3tmp49,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp154, FDPart3tmp2),
                                            FusedMulAddSIMD(FDPart3tmp1, aDD_dD111, MulSIMD(FDPart3tmp153, FDPart3tmp53))));
        const REAL_SIMD_ARRAY FDPart3tmp216 =
            FusedMulAddSIMD(FDPart3tmp175, FDPart3tmp56,
                            FusedMulAddSIMD(FDPart3tmp3, aDD_dD000,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp173, FDPart3tmp41),
                                                            FusedMulAddSIMD(FDPart3tmp171, aDD00, MulSIMD(FDPart3tmp174, FDPart3tmp53)))));
        const REAL_SIMD_ARRAY FDPart3tmp217 =
            FusedMulAddSIMD(FDPart3tmp164, FDPart3tmp49,
                            FusedMulAddSIMD(FDPart3tmp5, aDD_dD222,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp165, FDPart3tmp31),
                                                            FusedMulAddSIMD(FDPart3tmp157, aDD22, MulSIMD(FDPart3tmp163, FDPart3tmp56)))));
        const REAL_SIMD_ARRAY FDPart3tmp218 =
            FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp85, FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp87, MulSIMD(FDPart3tmp141, FDPart3tmp41)));
        const REAL_SIMD_ARRAY FDPart3tmp219 =
            FusedMulAddSIMD(FDPart3tmp137, FDPart3tmp2, FusedMulAddSIMD(FDPart3tmp138, FDPart3tmp89, MulSIMD(FDPart3tmp136, FDPart3tmp85)));
        const REAL_SIMD_ARRAY FDPart3tmp221 =
            FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp89, FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp31, MulSIMD(FDPart3tmp147, FDPart3tmp87)));
        const REAL_SIMD_ARRAY FDPart3tmp202 =
            FusedMulAddSIMD(FDPart3tmp5, aDD_dD220,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp138, FDPart3tmp31),
                                            FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp56, MulSIMD(FDPart3tmp137, FDPart3tmp49))));
        const REAL_SIMD_ARRAY FDPart3tmp204 =
            FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp49,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp142, FDPart3tmp2),
                                            FusedMulAddSIMD(FDPart3tmp1, aDD_dD112, MulSIMD(FDPart3tmp141, FDPart3tmp53))));
        const REAL_SIMD_ARRAY FDPart3tmp206 =
            FusedMulAddSIMD(FDPart3tmp5, aDD_dD221,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp143, FDPart3tmp31),
                                            FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp56, MulSIMD(FDPart3tmp142, FDPart3tmp49))));
        const REAL_SIMD_ARRAY FDPart3tmp209 =
            FusedMulAddSIMD(FDPart3tmp3, aDD_dD002,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp136, FDPart3tmp41),
                                            FusedMulAddSIMD(FDPart3tmp137, FDPart3tmp53, MulSIMD(FDPart3tmp138, FDPart3tmp56))));
        const REAL_SIMD_ARRAY FDPart3tmp211 =
            FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp49,
                            FusedMulAddSIMD(FDPart3tmp52, aDD11,
                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp148, FDPart3tmp2),
                                                            FusedMulAddSIMD(FDPart3tmp1, aDD_dD110, MulSIMD(FDPart3tmp147, FDPart3tmp53)))));
        const REAL_SIMD_ARRAY FDPart3tmp213 =
            FusedMulAddSIMD(FDPart3tmp3, aDD_dD001,
                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp147, FDPart3tmp41),
                                            FusedMulAddSIMD(FDPart3tmp148, FDPart3tmp53, MulSIMD(FDPart3tmp149, FDPart3tmp56))));
        const REAL_SIMD_ARRAY FDPart3tmp220 = AddSIMD(FDPart3tmp218, FusedMulAddSIMD(FDPart3tmp23, aDD_dD012, FDPart3tmp219));
        const REAL_SIMD_ARRAY FDPart3tmp222 = AddSIMD(FDPart3tmp218, FusedMulAddSIMD(FDPart3tmp32, aDD_dD021, FDPart3tmp221));
        const REAL_SIMD_ARRAY FDPart3tmp223 =
            FusedMulAddSIMD(FDPart3tmp32, aDD12, FusedMulAddSIMD(FDPart3tmp47, aDD_dD120, AddSIMD(FDPart3tmp219, FDPart3tmp221)));
        const REAL_SIMD_ARRAY FDPart3tmp224 =
            FusedMulAddSIMD(FDPart3tmp154, FDPart3tmp89,
                            FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp89,
                                            FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp87,
                                                            FusedMulAddSIMD(FDPart3tmp155, FDPart3tmp31,
                                                                            FusedMulAddSIMD(FDPart3tmp47, aDD_dD121,
                                                                                            FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp85,
                                                                                                            MulSIMD(FDPart3tmp142, FDPart3tmp2)))))));
        const REAL_SIMD_ARRAY FDPart3tmp225 = FusedMulAddSIMD(
            FDPart3tmp164, FDPart3tmp85,
            FusedMulAddSIMD(
                FDPart3tmp165, FDPart3tmp87,
                FusedMulAddSIMD(FDPart3tmp138, FDPart3tmp31,
                                FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp41,
                                                FusedMulAddSIMD(FDPart3tmp32, aDD_dD022,
                                                                FusedMulAddSIMD(aDD02, MulSIMD(f0_of_xx0__D0, f3_of_xx2__D2),
                                                                                FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp87,
                                                                                                MulSIMD(FDPart3tmp137, FDPart3tmp89))))))));
        const REAL_SIMD_ARRAY FDPart3tmp227 =
            FusedMulAddSIMD(FDPart3tmp154, FDPart3tmp85,
                            FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp89,
                                            FusedMulAddSIMD(FDPart3tmp153, FDPart3tmp41,
                                                            FusedMulAddSIMD(FDPart3tmp155, FDPart3tmp87,
                                                                            FusedMulAddSIMD(FDPart3tmp23, aDD_dD011,
                                                                                            FusedMulAddSIMD(FDPart3tmp147, FDPart3tmp85,
                                                                                                            MulSIMD(FDPart3tmp148, FDPart3tmp2)))))));
        const REAL_SIMD_ARRAY FDPart3tmp228 = FusedMulAddSIMD(
            FDPart3tmp164, FDPart3tmp2,
            FusedMulAddSIMD(
                FDPart3tmp165, FDPart3tmp89,
                FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp31,
                                FusedMulAddSIMD(FDPart3tmp163, FDPart3tmp85,
                                                FusedMulAddSIMD(FDPart3tmp47, aDD_dD122,
                                                                FusedMulAddSIMD(aDD12, MulSIMD(f0_of_xx0, f3_of_xx2__D2),
                                                                                FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp87,
                                                                                                MulSIMD(FDPart3tmp142, FDPart3tmp89))))))));
        const REAL_SIMD_ARRAY FDPart3tmp229 = FusedMulAddSIMD(
            FDPart3tmp174, FDPart3tmp89,
            FusedMulAddSIMD(
                FDPart3tmp175, FDPart3tmp31,
                FusedMulAddSIMD(FDPart3tmp138, FDPart3tmp87,
                                FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp87,
                                                FusedMulAddSIMD(FDPart3tmp32, aDD_dD020,
                                                                FusedMulAddSIMD(aDD02, MulSIMD(f0_of_xx0__DD00, f3_of_xx2),
                                                                                FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp41,
                                                                                                MulSIMD(FDPart3tmp137, FDPart3tmp85))))))));
        const REAL_SIMD_ARRAY FDPart3tmp230 = FusedMulAddSIMD(
            FDPart3tmp173, FDPart3tmp85,
            FusedMulAddSIMD(
                FDPart3tmp174, FDPart3tmp2,
                FusedMulAddSIMD(FDPart3tmp149, FDPart3tmp87,
                                FusedMulAddSIMD(FDPart3tmp169, aDD01,
                                                FusedMulAddSIMD(FDPart3tmp175, FDPart3tmp89,
                                                                FusedMulAddSIMD(FDPart3tmp23, aDD_dD010,
                                                                                FusedMulAddSIMD(FDPart3tmp147, FDPart3tmp41,
                                                                                                MulSIMD(FDPart3tmp148, FDPart3tmp85))))))));
        const REAL_SIMD_ARRAY FDPart3tmp199 = FusedMulAddSIMD(
            FDPart3_Integer_6, MulSIMD(FDPart3tmp85, FDPart3tmp98),
            FusedMulAddSIMD(
                FDPart3_Integer_6, MulSIMD(FDPart3tmp46, FDPart3tmp87),
                FusedMulAddSIMD(
                    FDPart3_Integer_6, MulSIMD(FDPart3tmp68, FDPart3tmp89),
                    FusedMulAddSIMD(
                        FDPart3_Integer_6, MulSIMD(FDPart3tmp31, FDPart3tmp80),
                        FusedMulAddSIMD(
                            FDPart3_Integer_6, MulSIMD(FDPart3tmp41, FDPart3tmp51),
                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp87, FDPart3tmp92),
                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp89, FDPart3tmp94),
                                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp2, FDPart3tmp70),
                                                                            MulSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp29, FDPart3tmp85))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp200 = FusedMulAddSIMD(
            FDPart3_Integer_6, MulSIMD(FDPart3tmp46, FDPart3tmp85),
            FusedMulAddSIMD(
                FDPart3_Integer_6, MulSIMD(FDPart3tmp39, FDPart3tmp87),
                FusedMulAddSIMD(
                    FDPart3_Integer_6, MulSIMD(FDPart3tmp41, FDPart3tmp54),
                    FusedMulAddSIMD(
                        FDPart3_Integer_6, MulSIMD(FDPart3tmp2, FDPart3tmp68),
                        FusedMulAddSIMD(
                            FDPart3_Integer_6, MulSIMD(FDPart3tmp31, FDPart3tmp82),
                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp80, FDPart3tmp89),
                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp85, FDPart3tmp92),
                                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp105, FDPart3tmp89),
                                                                            MulSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp109, FDPart3tmp87))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp226 = MulSIMD(FDPart3tmp224, FDPart3tmp65);
        const REAL_SIMD_ARRAY FDPart3tmp233 = FusedMulAddSIMD(
            FDPart3_Integer_6, MulSIMD(FDPart3tmp68, FDPart3tmp85),
            FusedMulAddSIMD(
                FDPart3_Integer_6, MulSIMD(FDPart3tmp41, FDPart3tmp46),
                FusedMulAddSIMD(
                    FDPart3_Integer_6, MulSIMD(FDPart3tmp62, FDPart3tmp89),
                    FusedMulAddSIMD(
                        FDPart3_Integer_6, MulSIMD(FDPart3tmp2, FDPart3tmp71),
                        FusedMulAddSIMD(
                            FDPart3_Integer_6, MulSIMD(FDPart3tmp31, FDPart3tmp83),
                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp80, FDPart3tmp87),
                                            FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp85, FDPart3tmp94),
                                                            FusedMulAddSIMD(FDPart3_Integer_6, FDPart3tmp119,
                                                                            MulSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp105, FDPart3tmp87))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp231 = FusedMulSubSIMD(
            FDPart3tmp120,
            FusedMulAddSIMD(
                FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp38, trK_dD2)),
                FusedMulAddSIMD(
                    FDPart3tmp128,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(cf_dD0,
                                    FusedMulAddSIMD(
                                        FDPart3tmp198, FDPart3tmp54,
                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp2, FDPart3tmp29),
                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp31, FDPart3tmp39),
                                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp41, FDPart3tmp43),
                                                                                        FusedMulAddSIMD(FDPart3tmp196, FDPart3tmp46,
                                                                                                        MulSIMD(FDPart3tmp197, FDPart3tmp51)))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp28, trK_dD1)),
                        FusedMulAddSIMD(
                            FDPart3tmp214, MulSIMD(FDPart3tmp28, FDPart3tmp69),
                            FusedMulAddSIMD(
                                FDPart3tmp128, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp200, cf_dD2)),
                                FusedMulAddSIMD(
                                    FDPart3_Integer_2, MulSIMD(FDPart3tmp230, FDPart3tmp51),
                                    FusedMulAddSIMD(
                                        FDPart3tmp206, MulSIMD(FDPart3tmp61, FDPart3tmp79),
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp223, FDPart3tmp46),
                                            FusedMulAddSIMD(
                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp229, FDPart3tmp54),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp227, FDPart3tmp98,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp228, FDPart3tmp80,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp226, FDPart3tmp79,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp227, FDPart3tmp29,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp224, FDPart3tmp68,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp225, FDPart3tmp39,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp222, FDPart3tmp46,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp222, FDPart3tmp92,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp220, FDPart3tmp46,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp220, FDPart3tmp92,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp216, FDPart3tmp43,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp217, FDPart3tmp82,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp211, FDPart3tmp29,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp213, FDPart3tmp51,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp204, FDPart3tmp68,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp209, FDPart3tmp54,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp109, FDPart3tmp225,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp202, FDPart3tmp39,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp121,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_NegativeOne_,
                                                                                                                                    FDPart3_Rational_2_3),
                                                                                                                                MulSIMD(FDPart3tmp42,
                                                                                                                                        trK_dD0)),
                                                                                                                            FusedMulSubSIMD(
                                                                                                                                FDPart3tmp105,
                                                                                                                                FDPart3tmp228,
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Rational_1_2,
                                                                                                                                        FDPart3tmp128),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp199,
                                                                                                                                        cf_dD1))))))))))))))))))))))))))))))),
            MulSIMD(
                FDPart3_Integer_8,
                MulSIMD(PI, FusedMulAddSIMD(
                                FDPart3tmp179,
                                MulSIMD(FDPart3tmp192,
                                        FusedMulAddSIMD(
                                            FDPart3tmp190,
                                            FusedMulAddSIMD(f3_of_xx2, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp180), MulSIMD(hDD12, vetU1)),
                                                            FusedMulSubSIMD(FDPart3tmp187, MulSIMD(f3_of_xx2, hDD02),
                                                                            MulSIMD(MulSIMD(FDPart3tmp180, FDPart3tmp8), DivSIMD(vetU2, f3_of_xx2)))),
                                            FusedMulAddSIMD(FDPart3tmp185, MulSIMD(FDPart3tmp47, hDD12),
                                                            FusedMulAddSIMD(FDPart3tmp182, FDPart3tmp36, MulSIMD(FDPart3tmp183, FDPart3tmp8))))),
                                FusedMulAddSIMD(
                                    FDPart3tmp179,
                                    MulSIMD(FusedMulAddSIMD(FDPart3tmp14, FDPart3tmp176, MulSIMD(FDPart3tmp176, FDPart3tmp19)),
                                            FusedMulAddSIMD(
                                                FDPart3tmp185, FDPart3tmp26,
                                                FusedMulAddSIMD(FDPart3tmp190,
                                                                FusedMulAddSIMD(FDPart3tmp188, MulSIMD(f0_of_xx0__D0, hDD02),
                                                                                FusedMulSubSIMD(FDPart3tmp12, DivSIMD(FDPart3tmp187, f0_of_xx0__D0),
                                                                                                MulSIMD(MulSIMD(FDPart3tmp180, f0_of_xx0__D0),
                                                                                                        MulSIMD(hDD01, vetU1)))),
                                                                FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp182, MulSIMD(FDPart3tmp183, FDPart3tmp36))))),
                                    MulSIMD(FDPart3tmp177,
                                            MulSIMD(FDPart3tmp179,
                                                    FusedMulAddSIMD(FDPart3tmp190,
                                                                    FusedMulAddSIMD(FDPart3tmp188, MulSIMD(f0_of_xx0, hDD12),
                                                                                    FusedMulSubSIMD(FDPart3tmp187, MulSIMD(f0_of_xx0, hDD01),
                                                                                                    MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp180),
                                                                                                            DivSIMD(vetU1, f0_of_xx0)))),
                                                                    FusedMulAddSIMD(FDPart3tmp183, MulSIMD(FDPart3tmp47, hDD12),
                                                                                    FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp185,
                                                                                                    MulSIMD(FDPart3tmp182, FDPart3tmp26)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp236 = FusedMulSubSIMD(
            FDPart3tmp120,
            FusedMulAddSIMD(
                FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp61, trK_dD2)),
                FusedMulAddSIMD(
                    FDPart3tmp128,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(cf_dD1,
                                    FusedMulAddSIMD(
                                        FDPart3tmp198, FDPart3tmp68,
                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp2, FDPart3tmp66),
                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp29, FDPart3tmp41),
                                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp31, FDPart3tmp62),
                                                                                        FusedMulAddSIMD(FDPart3tmp196, FDPart3tmp71,
                                                                                                        MulSIMD(FDPart3tmp197, FDPart3tmp70)))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp28, trK_dD0)),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3_Integer_2, FDPart3tmp22), MulSIMD(FDPart3tmp226, FDPart3tmp61),
                            FusedMulAddSIMD(
                                FDPart3tmp128, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp233, cf_dD2)),
                                FusedMulAddSIMD(
                                    FDPart3_Integer_2, MulSIMD(FDPart3tmp227, FDPart3tmp70),
                                    FusedMulAddSIMD(
                                        FDPart3tmp228, MulSIMD(FDPart3tmp65, FDPart3tmp81),
                                        FusedMulAddSIMD(
                                            FDPart3tmp230, FDPart3tmp98,
                                            FusedMulAddSIMD(
                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp222, FDPart3tmp68),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp229, FDPart3tmp92,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp230, FDPart3tmp29,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp228, FDPart3tmp62,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp229, FDPart3tmp46,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp223, FDPart3tmp94,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp225, FDPart3tmp80,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp220, FDPart3tmp94,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp223, FDPart3tmp68,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp217, FDPart3tmp83,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp220, FDPart3tmp68,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp214, FDPart3tmp66,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp216, FDPart3tmp51,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp211, FDPart3tmp70,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp213, FDPart3tmp29,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp206, FDPart3tmp62,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp209, FDPart3tmp46,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp202, FDPart3tmp80,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp204, FDPart3tmp71,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp121,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_NegativeOne_,
                                                                                                                                    FDPart3_Rational_2_3),
                                                                                                                                MulSIMD(FDPart3tmp65,
                                                                                                                                        trK_dD1)),
                                                                                                                            FusedMulSubSIMD(
                                                                                                                                FDPart3tmp105,
                                                                                                                                FDPart3tmp225,
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Rational_1_2,
                                                                                                                                        FDPart3tmp128),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp199,
                                                                                                                                        cf_dD0))))))))))))))))))))))))))))))),
            MulSIMD(
                FDPart3_Integer_8,
                MulSIMD(PI,
                        FusedMulAddSIMD(
                            FDPart3tmp179,
                            MulSIMD(FDPart3tmp232,
                                    FusedMulAddSIMD(
                                        FDPart3tmp190,
                                        FusedMulAddSIMD(f3_of_xx2, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp180), MulSIMD(hDD12, vetU1)),
                                                        FusedMulSubSIMD(FDPart3tmp187, MulSIMD(f3_of_xx2, hDD02),
                                                                        MulSIMD(MulSIMD(FDPart3tmp180, FDPart3tmp8), DivSIMD(vetU2, f3_of_xx2)))),
                                        FusedMulAddSIMD(FDPart3tmp185, MulSIMD(FDPart3tmp47, hDD12),
                                                        FusedMulAddSIMD(FDPart3tmp182, FDPart3tmp36, MulSIMD(FDPart3tmp183, FDPart3tmp8))))),
                            FusedMulAddSIMD(
                                FDPart3tmp179,
                                MulSIMD(
                                    FusedMulAddSIMD(FDPart3tmp17, FDPart3tmp176, MulSIMD(FDPart3tmp12, MulSIMD(FDPart3tmp176, FDPart3tmp8))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp190,
                                        FusedMulAddSIMD(FDPart3tmp188, MulSIMD(f0_of_xx0, hDD12),
                                                        FusedMulSubSIMD(FDPart3tmp187, MulSIMD(f0_of_xx0, hDD01),
                                                                        MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp180), DivSIMD(vetU1, f0_of_xx0)))),
                                        FusedMulAddSIMD(FDPart3tmp183, MulSIMD(FDPart3tmp47, hDD12),
                                                        FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp185, MulSIMD(FDPart3tmp182, FDPart3tmp26))))),
                                MulSIMD(FDPart3tmp177,
                                        MulSIMD(FDPart3tmp179,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp185, FDPart3tmp26,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp190,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp188, MulSIMD(f0_of_xx0__D0, hDD02),
                                                            FusedMulSubSIMD(FDPart3tmp12, DivSIMD(FDPart3tmp187, f0_of_xx0__D0),
                                                                            MulSIMD(MulSIMD(FDPart3tmp180, f0_of_xx0__D0), MulSIMD(hDD01, vetU1)))),
                                                        FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp182, MulSIMD(FDPart3tmp183, FDPart3tmp36)))))))))));
        const REAL_SIMD_ARRAY FDPart3tmp237 = FusedMulSubSIMD(
            FDPart3tmp120,
            FusedMulAddSIMD(
                FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp61, trK_dD1)),
                FusedMulAddSIMD(
                    FDPart3tmp128,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(cf_dD2,
                                    FusedMulAddSIMD(
                                        FDPart3tmp198, FDPart3tmp82,
                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp2, FDPart3tmp62),
                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp31, FDPart3tmp77),
                                                                        FusedMulAddSIMD(FDPart3_Integer_6, MulSIMD(FDPart3tmp39, FDPart3tmp41),
                                                                                        FusedMulAddSIMD(FDPart3tmp196, FDPart3tmp83,
                                                                                                        MulSIMD(FDPart3tmp197, FDPart3tmp80)))))))),
                    FusedMulAddSIMD(
                        FDPart3tmp121, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_2_3), MulSIMD(FDPart3tmp38, trK_dD0)),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3_Integer_2, FDPart3tmp228), MulSIMD(FDPart3tmp61, FDPart3tmp81),
                            FusedMulAddSIMD(
                                FDPart3tmp128, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp233, cf_dD1)),
                                FusedMulAddSIMD(
                                    FDPart3tmp206, MulSIMD(FDPart3tmp61, FDPart3tmp81),
                                    FusedMulAddSIMD(
                                        FDPart3tmp214, MulSIMD(FDPart3tmp61, FDPart3tmp69),
                                        FusedMulAddSIMD(
                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp220, FDPart3tmp80),
                                            FusedMulAddSIMD(
                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp225, FDPart3tmp82),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp230, FDPart3tmp46,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp230, FDPart3tmp92,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp227, FDPart3tmp94,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp229, FDPart3tmp39,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp226, FDPart3tmp81,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp227, FDPart3tmp68,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp223, FDPart3tmp80,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp224, FDPart3tmp62,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp217, FDPart3tmp77,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp222, FDPart3tmp80,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp213, FDPart3tmp46,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp216, FDPart3tmp54,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp209, FDPart3tmp39,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp211, FDPart3tmp68,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp202, FDPart3tmp82,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp204, FDPart3tmp62,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp105, FDPart3tmp223,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp109, FDPart3tmp229,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp121,
                                                                                                                            MulSIMD(
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_NegativeOne_,
                                                                                                                                    FDPart3_Rational_2_3),
                                                                                                                                MulSIMD(FDPart3tmp76,
                                                                                                                                        trK_dD2)),
                                                                                                                            FusedMulSubSIMD(
                                                                                                                                FDPart3tmp105,
                                                                                                                                FDPart3tmp222,
                                                                                                                                MulSIMD(
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Rational_1_2,
                                                                                                                                        FDPart3tmp128),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp200,
                                                                                                                                        cf_dD0))))))))))))))))))))))))))))))),
            MulSIMD(
                FDPart3_Integer_8,
                MulSIMD(PI,
                        FusedMulAddSIMD(
                            FDPart3tmp179,
                            MulSIMD(FDPart3tmp232,
                                    FusedMulAddSIMD(
                                        FDPart3tmp190,
                                        FusedMulAddSIMD(FDPart3tmp188, MulSIMD(f0_of_xx0, hDD12),
                                                        FusedMulSubSIMD(FDPart3tmp187, MulSIMD(f0_of_xx0, hDD01),
                                                                        MulSIMD(MulSIMD(FDPart3tmp16, FDPart3tmp180), DivSIMD(vetU1, f0_of_xx0)))),
                                        FusedMulAddSIMD(FDPart3tmp183, MulSIMD(FDPart3tmp47, hDD12),
                                                        FusedMulAddSIMD(FDPart3tmp16, FDPart3tmp185, MulSIMD(FDPart3tmp182, FDPart3tmp26))))),
                            FusedMulAddSIMD(
                                FDPart3tmp179,
                                MulSIMD(FusedMulAddSIMD(FDPart3tmp10, FDPart3tmp176, MulSIMD(FDPart3tmp12, MulSIMD(FDPart3tmp16, FDPart3tmp176))),
                                        FusedMulAddSIMD(
                                            FDPart3tmp190,
                                            FusedMulAddSIMD(f3_of_xx2, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp180), MulSIMD(hDD12, vetU1)),
                                                            FusedMulSubSIMD(FDPart3tmp187, MulSIMD(f3_of_xx2, hDD02),
                                                                            MulSIMD(MulSIMD(FDPart3tmp180, FDPart3tmp8), DivSIMD(vetU2, f3_of_xx2)))),
                                            FusedMulAddSIMD(FDPart3tmp185, MulSIMD(FDPart3tmp47, hDD12),
                                                            FusedMulAddSIMD(FDPart3tmp182, FDPart3tmp36, MulSIMD(FDPart3tmp183, FDPart3tmp8))))),
                                MulSIMD(FDPart3tmp179,
                                        MulSIMD(FDPart3tmp192,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp185, FDPart3tmp26,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp190,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp188, MulSIMD(f0_of_xx0__D0, hDD02),
                                                            FusedMulSubSIMD(FDPart3tmp12, DivSIMD(FDPart3tmp187, f0_of_xx0__D0),
                                                                            MulSIMD(MulSIMD(FDPart3tmp180, f0_of_xx0__D0), MulSIMD(hDD01, vetU1)))),
                                                        FusedMulAddSIMD(FDPart3tmp12, FDPart3tmp182, MulSIMD(FDPart3tmp183, FDPart3tmp36)))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(
            MulSIMD(aDD01, f0_of_xx0),
            MulSIMD(
                MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                MulSIMD(f0_of_xx0__D0,
                        FusedMulAddSIMD(
                            FDPart3tmp85, FDPart3tmp98,
                            FusedMulAddSIMD(
                                FDPart3tmp46, FDPart3tmp87,
                                FusedMulAddSIMD(
                                    FDPart3tmp68, FDPart3tmp89,
                                    FusedMulAddSIMD(
                                        FDPart3tmp31, FDPart3tmp80,
                                        FusedMulAddSIMD(FDPart3tmp41, FDPart3tmp51,
                                                        FusedMulAddSIMD(FDPart3tmp87, FDPart3tmp92,
                                                                        FusedMulAddSIMD(FDPart3tmp89, FDPart3tmp94,
                                                                                        FusedMulAddSIMD(FDPart3tmp2, FDPart3tmp70,
                                                                                                        MulSIMD(FDPart3tmp29, FDPart3tmp85))))))))))),
            FusedMulAddSIMD(
                MulSIMD(aDD02, f0_of_xx0__D0),
                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                        MulSIMD(f3_of_xx2,
                                FusedMulAddSIMD(
                                    FDPart3tmp46, FDPart3tmp85,
                                    FusedMulAddSIMD(
                                        FDPart3tmp39, FDPart3tmp87,
                                        FusedMulAddSIMD(
                                            FDPart3tmp41, FDPart3tmp54,
                                            FusedMulAddSIMD(
                                                FDPart3tmp2, FDPart3tmp68,
                                                FusedMulAddSIMD(
                                                    FDPart3tmp31, FDPart3tmp82,
                                                    FusedMulAddSIMD(FDPart3tmp80, FDPart3tmp89,
                                                                    FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp92,
                                                                                    FusedMulAddSIMD(FDPart3tmp105, FDPart3tmp89,
                                                                                                    MulSIMD(FDPart3tmp109, FDPart3tmp87))))))))))),
                NegFusedMulAddSIMD(
                    FDPart3tmp5,
                    MulSIMD(aDD22,
                            FusedMulAddSIMD(FDPart3tmp39, FDPart3tmp41,
                                            FusedMulAddSIMD(FDPart3tmp49, FDPart3tmp83,
                                                            FusedMulAddSIMD(FDPart3tmp53, FDPart3tmp80,
                                                                            FusedMulAddSIMD(FDPart3tmp56, FDPart3tmp82,
                                                                                            FusedMulAddSIMD(FDPart3tmp2, FDPart3tmp62,
                                                                                                            MulSIMD(FDPart3tmp31, FDPart3tmp77))))))),
                    FusedMulAddSIMD(
                        T4UU00, MulSIMD(MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_), MulSIMD(PI, MulSIMD(alpha, alpha))),
                        FusedMulAddSIMD(
                            FDPart3tmp120,
                            FusedMulAddSIMD(
                                FDPart3tmp38,
                                MulSIMD(MulSIMD(FDPart3tmp121, FDPart3tmp123),
                                        MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(cf_dD0, cf_dD2))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3tmp121, FDPart3tmp123),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp76, MulSIMD(cf_dD2, cf_dD2))),
                                    FusedMulAddSIMD(
                                        FDPart3tmp28,
                                        MulSIMD(MulSIMD(FDPart3tmp121, FDPart3tmp123),
                                                MulSIMD(MulSIMD(FDPart3_Integer_4, FDPart3_NegativeOne_), MulSIMD(cf_dD0, cf_dD1))),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3tmp121, FDPart3tmp123),
                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp42, MulSIMD(cf_dD0, cf_dD0))),
                                            FusedMulAddSIMD(
                                                MulSIMD(FDPart3tmp121, FDPart3tmp123),
                                                MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                        MulSIMD(FDPart3tmp65, MulSIMD(cf_dD1, cf_dD1))),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp121,
                                                    MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                                            MulSIMD(FDPart3tmp65,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp128,
                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                MulSIMD(FDPart3tmp154, cf_dD1)),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp128,
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                    MulSIMD(FDPart3tmp155, cf_dD2)),
                                                                            FusedMulSubSIMD(
                                                                                FDPart3tmp129,
                                                                                FusedMulSubSIMD(FDPart3tmp128, MulSIMD(cf_dD1, cf_dD1), cf_dDD11),
                                                                                MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp128),
                                                                                        MulSIMD(FDPart3tmp153, cf_dD0))))))),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp121,
                                                        MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                                                MulSIMD(FDPart3tmp76,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp128,
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                    MulSIMD(FDPart3tmp164, cf_dD1)),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp128,
                                                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                        MulSIMD(FDPart3tmp165, cf_dD2)),
                                                                                FusedMulSubSIMD(
                                                                                    FDPart3tmp129,
                                                                                    FusedMulSubSIMD(FDPart3tmp128, MulSIMD(cf_dD2, cf_dD2), cf_dDD22),
                                                                                    MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp128),
                                                                                            MulSIMD(FDPart3tmp163, cf_dD0))))))),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp121,
                                                            MulSIMD(
                                                                MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_),
                                                                MulSIMD(FDPart3tmp61,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp128,
                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                    MulSIMD(FDPart3tmp142, cf_dD1)),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp128,
                                                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                        MulSIMD(FDPart3tmp143, cf_dD2)),
                                                                                FusedMulSubSIMD(
                                                                                    FDPart3tmp129,
                                                                                    FusedMulSubSIMD(FDPart3tmp128, MulSIMD(cf_dD1, cf_dD2), cf_dDD12),
                                                                                    MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp128),
                                                                                            MulSIMD(FDPart3tmp141, cf_dD0))))))),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp121,
                                                                MulSIMD(MulSIMD(FDPart3_Integer_8, FDPart3_NegativeOne_),
                                                                        MulSIMD(FDPart3tmp42,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp128,
                                                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                            MulSIMD(FDPart3tmp174, cf_dD1)),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp128,
                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                MulSIMD(FDPart3tmp175, cf_dD2)),
                                                                                        FusedMulSubSIMD(
                                                                                            FDPart3tmp129,
                                                                                            FusedMulSubSIMD(FDPart3tmp128, MulSIMD(cf_dD0, cf_dD0),
                                                                                                            cf_dDD00),
                                                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp128),
                                                                                                    MulSIMD(FDPart3tmp173, cf_dD0))))))),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp122, MulSIMD(FDPart3tmp61, RbarDD12),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp121,
                                                                        MulSIMD(
                                                                            MulSIMD(FDPart3_Integer_16, FDPart3_NegativeOne_),
                                                                            MulSIMD(
                                                                                FDPart3tmp38,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp128,
                                                                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                            MulSIMD(FDPart3tmp137, cf_dD1)),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp128,
                                                                                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                MulSIMD(FDPart3tmp138, cf_dD2)),
                                                                                        FusedMulSubSIMD(
                                                                                            FDPart3tmp129,
                                                                                            FusedMulSubSIMD(FDPart3tmp128, MulSIMD(cf_dD0, cf_dD2),
                                                                                                            cf_dDD02),
                                                                                            MulSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp128),
                                                                                                    MulSIMD(FDPart3tmp136, cf_dD0))))))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp122, MulSIMD(FDPart3tmp28, RbarDD01),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp122, MulSIMD(FDPart3tmp38, RbarDD02),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp121, MulSIMD(FDPart3tmp65, RbarDD11),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp121, MulSIMD(FDPart3tmp76, RbarDD22),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp61,
                                                                                            MulSIMD(MulSIMD(FDPart3tmp121, FDPart3tmp123),
                                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_4,
                                                                                                                    FDPart3_NegativeOne_),
                                                                                                            MulSIMD(cf_dD1, cf_dD2))),
                                                                                            FusedMulSubSIMD(
                                                                                                FDPart3tmp121, MulSIMD(FDPart3tmp42, RbarDD00),
                                                                                                MulSIMD(
                                                                                                    MulSIMD(FDPart3_Integer_16, FDPart3tmp121),
                                                                                                    MulSIMD(
                                                                                                        FDPart3tmp28,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp128,
                                                                                                            MulSIMD(MulSIMD(FDPart3_NegativeOne_,
                                                                                                                            FDPart3_Rational_1_2),
                                                                                                                    MulSIMD(FDPart3tmp148, cf_dD1)),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp128,
                                                                                                                MulSIMD(
                                                                                                                    MulSIMD(FDPart3_NegativeOne_,
                                                                                                                            FDPart3_Rational_1_2),
                                                                                                                    MulSIMD(FDPart3tmp149, cf_dD2)),
                                                                                                                FusedMulSubSIMD(
                                                                                                                    FDPart3tmp129,
                                                                                                                    FusedMulSubSIMD(
                                                                                                                        FDPart3tmp128,
                                                                                                                        MulSIMD(cf_dD0, cf_dD1),
                                                                                                                        cf_dDD01),
                                                                                                                    MulSIMD(
                                                                                                                        MulSIMD(FDPart3_Rational_1_2,
                                                                                                                                FDPart3tmp128),
                                                                                                                        MulSIMD(
                                                                                                                            FDPart3tmp147,
                                                                                                                            cf_dD0)))))))))))))))))))))))),
                            NegFusedMulAddSIMD(
                                FDPart3tmp3,
                                MulSIMD(aDD00,
                                        FusedMulAddSIMD(
                                            FDPart3tmp41, FDPart3tmp43,
                                            FusedMulAddSIMD(FDPart3tmp46, FDPart3tmp49,
                                                            FusedMulAddSIMD(FDPart3tmp51, FDPart3tmp53,
                                                                            FusedMulAddSIMD(FDPart3tmp54, FDPart3tmp56,
                                                                                            FusedMulAddSIMD(FDPart3tmp2, FDPart3tmp29,
                                                                                                            MulSIMD(FDPart3tmp31, FDPart3tmp39))))))),
                                FusedMulAddSIMD(
                                    MulSIMD(aDD12, f0_of_xx0),
                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                            MulSIMD(f3_of_xx2,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp68, FDPart3tmp85,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp41, FDPart3tmp46,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp62, FDPart3tmp89,
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp2, FDPart3tmp71,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp31, FDPart3tmp83,
                                                                        FusedMulAddSIMD(FDPart3tmp80, FDPart3tmp87,
                                                                                        FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp94,
                                                                                                        FusedMulAddSIMD(FDPart3tmp105, FDPart3tmp87,
                                                                                                                        FDPart3tmp119)))))))))),
                                    FusedMulSubSIMD(
                                        FDPart3_Rational_2_3, MulSIMD(trK, trK),
                                        MulSIMD(FDPart3tmp1,
                                                MulSIMD(aDD11,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp31, FDPart3tmp62,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp49, FDPart3tmp71,
                                                                FusedMulAddSIMD(FDPart3tmp53, FDPart3tmp70,
                                                                                FusedMulAddSIMD(FDPart3tmp56, FDPart3tmp68,
                                                                                                FusedMulAddSIMD(FDPart3tmp2, FDPart3tmp66,
                                                                                                                MulSIMD(FDPart3tmp29,
                                                                                                                        FDPart3tmp41))))))))))))))));
        const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(
            FDPart3tmp8, MulSIMD(FDPart3tmp237, FDPart3tmp237),
            FusedMulAddSIMD(MulSIMD(FDPart3tmp231, FDPart3tmp236), MulSIMD(FDPart3tmp52, hDD01),
                            FusedMulAddSIMD(MulSIMD(FDPart3tmp231, FDPart3tmp237), MulSIMD(FDPart3tmp55, hDD02),
                                            FusedMulAddSIMD(MulSIMD(FDPart3tmp236, FDPart3tmp237), MulSIMD(FDPart3tmp48, hDD12),
                                                            FusedMulAddSIMD(FDPart3tmp12, MulSIMD(FDPart3tmp231, FDPart3tmp231),
                                                                            MulSIMD(FDPart3tmp16, MulSIMD(FDPart3tmp236, FDPart3tmp236)))))));

        WriteSIMD(&diagnostic_output_gfs[IDX4(HGF, i0, i1, i2)], __RHS_exp_0);
        WriteSIMD(&diagnostic_output_gfs[IDX4(MSQUAREDGF, i0, i1, i2)], __RHS_exp_1);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0 += simd_width)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}
