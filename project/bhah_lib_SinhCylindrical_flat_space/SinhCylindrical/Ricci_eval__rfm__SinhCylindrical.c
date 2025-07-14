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
} // END FUNCTION SIMD_fd_function_dD0_fdorder4
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
} // END FUNCTION SIMD_fd_function_dD1_fdorder4
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
} // END FUNCTION SIMD_fd_function_dD2_fdorder4
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
} // END FUNCTION SIMD_fd_function_dDD00_fdorder4
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
} // END FUNCTION SIMD_fd_function_dDD01_fdorder4
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
} // END FUNCTION SIMD_fd_function_dDD02_fdorder4
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
} // END FUNCTION SIMD_fd_function_dDD11_fdorder4
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
} // END FUNCTION SIMD_fd_function_dDD12_fdorder4
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
} // END FUNCTION SIMD_fd_function_dDD22_fdorder4

/**
 * Set Ricci tensor.
 */
void Ricci_eval__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                      const rfm_struct *restrict rfmstruct, const REAL *restrict in_gfs, REAL *restrict auxevol_gfs) {
#include "../set_CodeParameters-simd.h"
#pragma omp parallel for collapse(2)
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      const double NOSIMDf3_of_xx2 = rfmstruct->f3_of_xx2[i2];
      MAYBE_UNUSED const REAL_SIMD_ARRAY f3_of_xx2 = ConstSIMD(NOSIMDf3_of_xx2);
      const double NOSIMDf3_of_xx2__D2 = rfmstruct->f3_of_xx2__D2[i2];
      MAYBE_UNUSED const REAL_SIMD_ARRAY f3_of_xx2__D2 = ConstSIMD(NOSIMDf3_of_xx2__D2);
      const double NOSIMDf3_of_xx2__DD22 = rfmstruct->f3_of_xx2__DD22[i2];
      MAYBE_UNUSED const REAL_SIMD_ARRAY f3_of_xx2__DD22 = ConstSIMD(NOSIMDf3_of_xx2__DD22);

      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += simd_width) {
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0 = ReadSIMD(&rfmstruct->f0_of_xx0[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__D0 = ReadSIMD(&rfmstruct->f0_of_xx0__D0[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__DD00 = ReadSIMD(&rfmstruct->f0_of_xx0__DD00[i0]);
        MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0__DDD000 = ReadSIMD(&rfmstruct->f0_of_xx0__DDD000[i0]);

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
        static const double dblFDPart3_Integer_1 = 1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_Integer_1 = ConstSIMD(dblFDPart3_Integer_1);

        static const double dblFDPart3_Integer_2 = 2.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_2 = ConstSIMD(dblFDPart3_Integer_2);

        static const double dblFDPart3_Integer_3 = 3.0;
        const REAL_SIMD_ARRAY FDPart3_Integer_3 = ConstSIMD(dblFDPart3_Integer_3);

        static const double dblFDPart3_NegativeOne_ = -1.0;
        MAYBE_UNUSED const REAL_SIMD_ARRAY FDPart3_NegativeOne_ = ConstSIMD(dblFDPart3_NegativeOne_);

        static const double dblFDPart3_Rational_1_2 = 1.0 / 2.0;
        const REAL_SIMD_ARRAY FDPart3_Rational_1_2 = ConstSIMD(dblFDPart3_Rational_1_2);

        const REAL_SIMD_ARRAY FDPart3tmp1 = MulSIMD(f0_of_xx0__D0, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp2 = MulSIMD(f0_of_xx0__D0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp5 = DivSIMD(FDPart3_Integer_1, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp8 = MulSIMD(f3_of_xx2, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp10 = MulSIMD(f0_of_xx0, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp23 = MulSIMD(f0_of_xx0__D0, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp29 = MulSIMD(f0_of_xx0, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp37 = DivSIMD(FDPart3_Integer_1, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp41 = MulSIMD(f0_of_xx0, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp43 = DivSIMD(FDPart3_Integer_1, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp50 = MulSIMD(f0_of_xx0, f0_of_xx0__DD00);
        const REAL_SIMD_ARRAY FDPart3tmp62 = MulSIMD(FDPart3_Integer_2, f0_of_xx0__D0);
        const REAL_SIMD_ARRAY FDPart3tmp112 = MulSIMD(FDPart3_Integer_2, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp169 = MulSIMD(f0_of_xx0, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp282 = MulSIMD(f3_of_xx2, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp290 = MulSIMD(f0_of_xx0__DD00, f3_of_xx2);
        const REAL_SIMD_ARRAY FDPart3tmp293 = MulSIMD(f0_of_xx0__D0, f3_of_xx2__D2);
        const REAL_SIMD_ARRAY FDPart3tmp3 = MulSIMD(FDPart3tmp2, hDD00);
        const REAL_SIMD_ARRAY FDPart3tmp9 = MulSIMD(FDPart3tmp8, MulSIMD(hDD02, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp12 = MulSIMD(FDPart3tmp8, hDD22);
        const REAL_SIMD_ARRAY FDPart3tmp18 = MulSIMD(FDPart3tmp10, hDD11);
        const REAL_SIMD_ARRAY FDPart3tmp24 = MulSIMD(FDPart3tmp23, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp25 = MulSIMD(FDPart3tmp23, hDD02);
        const REAL_SIMD_ARRAY FDPart3tmp31 = MulSIMD(FDPart3tmp29, hDD12);
        const REAL_SIMD_ARRAY FDPart3tmp38 = MulSIMD(FDPart3tmp37, f3_of_xx2__D2);
        const REAL_SIMD_ARRAY FDPart3tmp47 = MulSIMD(FDPart3tmp41, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp52 = AddSIMD(FDPart3tmp2, FDPart3tmp50);
        const REAL_SIMD_ARRAY FDPart3tmp56 = MulSIMD(FDPart3tmp2, hDD_dD001);
        const REAL_SIMD_ARRAY FDPart3tmp58 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp43, f0_of_xx0__D0));
        const REAL_SIMD_ARRAY FDPart3tmp65 = DivSIMD(FDPart3_Integer_1, FDPart3tmp2);
        const REAL_SIMD_ARRAY FDPart3tmp68 = MulSIMD(FDPart3tmp62, f0_of_xx0__DD00);
        const REAL_SIMD_ARRAY FDPart3tmp72 = DivSIMD(FDPart3tmp2, FDPart3tmp10);
        const REAL_SIMD_ARRAY FDPart3tmp77 = MulSIMD(FDPart3tmp8, hDD_dD220);
        const REAL_SIMD_ARRAY FDPart3tmp80 = MulSIMD(FDPart3tmp2, hDD_dD002);
        const REAL_SIMD_ARRAY FDPart3tmp81 = MulSIMD(FDPart3tmp41, hDD_dD012);
        const REAL_SIMD_ARRAY FDPart3tmp82 = MulSIMD(FDPart3tmp23, hDD_dD021);
        const REAL_SIMD_ARRAY FDPart3tmp101 = MulSIMD(FDPart3_Integer_2, FDPart3tmp41);
        const REAL_SIMD_ARRAY FDPart3tmp138 = MulSIMD(FDPart3tmp8, hDD_dD221);
        const REAL_SIMD_ARRAY FDPart3tmp139 = MulSIMD(FDPart3tmp10, hDD_dD112);
        const REAL_SIMD_ARRAY FDPart3tmp184 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp10, hDD_dD112));
        const REAL_SIMD_ARRAY FDPart3tmp186 = MulSIMD(FDPart3tmp10, hDD_dD111);
        const REAL_SIMD_ARRAY FDPart3tmp191 = MulSIMD(FDPart3tmp5, f0_of_xx0);
        const REAL_SIMD_ARRAY FDPart3tmp222 = FusedMulAddSIMD(FDPart3tmp43, lambdaU0, MulSIMD(FDPart3tmp43, lambdaU_dD11));
        const REAL_SIMD_ARRAY FDPart3tmp228 = MulSIMD(FDPart3_Integer_2, FDPart3tmp10);
        const REAL_SIMD_ARRAY FDPart3tmp292 = FusedMulAddSIMD(FDPart3tmp23, hDD_dDD0202, MulSIMD(f0_of_xx0__DD00, MulSIMD(f3_of_xx2__D2, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp302 = DivSIMD(MulSIMD(f3_of_xx2__D2, f3_of_xx2__D2), FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp348 = FusedMulAddSIMD(FDPart3tmp29, hDD_dDD1202, MulSIMD(FDPart3tmp293, hDD12));
        const REAL_SIMD_ARRAY FDPart3tmp4 = AddSIMD(FDPart3tmp2, FDPart3tmp3);
        const REAL_SIMD_ARRAY FDPart3tmp13 = AddSIMD(FDPart3tmp12, FDPart3tmp8);
        const REAL_SIMD_ARRAY FDPart3tmp15 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp10), MulSIMD(FDPart3tmp2, MulSIMD(hDD01, hDD01)));
        const REAL_SIMD_ARRAY FDPart3tmp17 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp10), MulSIMD(FDPart3tmp8, MulSIMD(hDD12, hDD12)));
        const REAL_SIMD_ARRAY FDPart3tmp19 = AddSIMD(FDPart3tmp10, FDPart3tmp18);
        const REAL_SIMD_ARRAY FDPart3tmp20 = MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3tmp2), MulSIMD(FDPart3tmp8, MulSIMD(hDD02, hDD02)));
        const REAL_SIMD_ARRAY FDPart3tmp36 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp2, hDD_dD002));
        const REAL_SIMD_ARRAY FDPart3tmp45 =
            FusedMulAddSIMD(FDPart3tmp41, hDD_dD011, FusedMulSubSIMD(FDPart3tmp41, hDD00, MulSIMD(FDPart3tmp41, hDD11)));
        const REAL_SIMD_ARRAY FDPart3tmp57 = NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, hDD01), FDPart3tmp56);
        const REAL_SIMD_ARRAY FDPart3tmp66 = MulSIMD(FDPart3tmp65, MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp69 = FusedMulAddSIMD(FDPart3tmp2, hDD_dD000, MulSIMD(FDPart3tmp68, hDD00));
        const REAL_SIMD_ARRAY FDPart3tmp73 = NegFusedMulAddSIMD(FDPart3tmp43, f0_of_xx0__DD00, FDPart3tmp72);
        const REAL_SIMD_ARRAY FDPart3tmp84 = FusedMulAddSIMD(FDPart3tmp29, hDD_dD120, FDPart3tmp24);
        const REAL_SIMD_ARRAY FDPart3tmp103 = FusedMulAddSIMD(FDPart3tmp10, hDD_dD110, MulSIMD(FDPart3tmp101, hDD11));
        const REAL_SIMD_ARRAY FDPart3tmp110 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp5, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp140 = SubSIMD(FDPart3tmp82, FDPart3tmp24);
        const REAL_SIMD_ARRAY FDPart3tmp162 = FusedMulAddSIMD(FDPart3tmp8, hDD_dD222, MulSIMD(FDPart3tmp112, MulSIMD(f3_of_xx2__D2, hDD22)));
        const REAL_SIMD_ARRAY FDPart3tmp166 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp8, hDD_dD220));
        const REAL_SIMD_ARRAY FDPart3tmp171 =
            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp29, hDD_dD122),
                            FusedMulSubSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp169, f3_of_xx2__D2), MulSIMD(FDPart3tmp8, hDD_dD221)));
        const REAL_SIMD_ARRAY FDPart3tmp185 = FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp29, hDD_dD121), FDPart3tmp184);
        const REAL_SIMD_ARRAY FDPart3tmp187 = MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp10, hDD_dD110));
        const REAL_SIMD_ARRAY FDPart3tmp209 = MulSIMD(FDPart3tmp101, hDD01);
        const REAL_SIMD_ARRAY FDPart3tmp224 = FusedMulSubSIMD(FDPart3tmp5, lambdaU_dD01, MulSIMD(FDPart3tmp5, lambdaU1));
        const REAL_SIMD_ARRAY FDPart3tmp227 = FusedMulAddSIMD(FDPart3tmp41, hDD_dDD0101, MulSIMD(FDPart3tmp52, hDD_dD011));
        const REAL_SIMD_ARRAY FDPart3tmp229 = FusedMulAddSIMD(FDPart3tmp228, hDD01, FDPart3tmp186);
        const REAL_SIMD_ARRAY FDPart3tmp291 = FusedMulAddSIMD(FDPart3tmp1, f3_of_xx2__D2, MulSIMD(FDPart3tmp23, hDD_dD022));
        const REAL_SIMD_ARRAY FDPart3tmp294 = FusedMulAddSIMD(FDPart3tmp23, hDD_dD020, MulSIMD(FDPart3tmp282, f0_of_xx0__DD00));
        const REAL_SIMD_ARRAY FDPart3tmp296 = FusedMulAddSIMD(FDPart3tmp29, hDD02, MulSIMD(FDPart3tmp29, hDD_dD121));
        const REAL_SIMD_ARRAY FDPart3tmp298 = FusedMulAddSIMD(FDPart3tmp169, f3_of_xx2__D2, MulSIMD(FDPart3tmp29, hDD_dD122));
        const REAL_SIMD_ARRAY FDPart3tmp26 = FusedMulSubSIMD(FDPart3tmp10, MulSIMD(FDPart3tmp24, hDD01), MulSIMD(FDPart3tmp19, FDPart3tmp25));
        const REAL_SIMD_ARRAY FDPart3tmp32 =
            FusedMulSubSIMD(MulSIMD(FDPart3tmp2, FDPart3tmp29), MulSIMD(hDD01, hDD02), MulSIMD(FDPart3tmp31, FDPart3tmp4));
        const REAL_SIMD_ARRAY FDPart3tmp39 = FusedMulAddSIMD(FDPart3tmp19, FDPart3tmp4, FDPart3tmp15);
        const REAL_SIMD_ARRAY FDPart3tmp40 = FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp4, FDPart3tmp20);
        const REAL_SIMD_ARRAY FDPart3tmp46 =
            MulSIMD(FDPart3tmp43, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp45, f0_of_xx0__D0)));
        const REAL_SIMD_ARRAY FDPart3tmp49 = FusedMulSubSIMD(FDPart3tmp41, FDPart3tmp9, MulSIMD(FDPart3tmp13, FDPart3tmp47));
        const REAL_SIMD_ARRAY FDPart3tmp54 = FusedMulAddSIMD(FDPart3tmp41, hDD_dD010, MulSIMD(FDPart3tmp52, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp59 = FusedMulAddSIMD(FDPart3tmp2, hDD_dDD0001, MulSIMD(FDPart3tmp57, FDPart3tmp58));
        const REAL_SIMD_ARRAY FDPart3tmp67 = NegFusedMulAddSIMD(FDPart3tmp5, f0_of_xx0__DDD000, FDPart3tmp66);
        const REAL_SIMD_ARRAY FDPart3tmp75 =
            MulSIMD(FDPart3_NegativeOne_, MulSIMD(FDPart3tmp47, FusedMulSubSIMD(FDPart3tmp43, f0_of_xx0__DD00, FDPart3tmp72)));
        const REAL_SIMD_ARRAY FDPart3tmp85 = AddSIMD(FDPart3tmp84, SubSIMD(FDPart3tmp81, FDPart3tmp82));
        const REAL_SIMD_ARRAY FDPart3tmp104 = MulSIMD(FDPart3_Rational_1_2, AddSIMD(FDPart3tmp101, FDPart3tmp103));
        const REAL_SIMD_ARRAY FDPart3tmp113 = FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp23, hDD_dD020),
                                                              FusedMulAddSIMD(FDPart3tmp112, MulSIMD(f0_of_xx0__DD00, hDD02), FDPart3tmp36));
        const REAL_SIMD_ARRAY FDPart3tmp114 =
            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp52, hDD01), FusedMulSubSIMD(FDPart3tmp101, hDD_dD010, FDPart3tmp56));
        const REAL_SIMD_ARRAY FDPart3tmp115 = AddSIMD(FDPart3tmp68, FDPart3tmp69);
        const REAL_SIMD_ARRAY FDPart3tmp141 = AddSIMD(FDPart3tmp81, NegFusedMulAddSIMD(FDPart3tmp29, hDD_dD120, FDPart3tmp140));
        const REAL_SIMD_ARRAY FDPart3tmp163 = FusedMulAddSIMD(FDPart3tmp112, f3_of_xx2__D2, FDPart3tmp162);
        const REAL_SIMD_ARRAY FDPart3tmp167 = FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp23, hDD_dD022),
                                                              FusedMulAddSIMD(FDPart3tmp62, MulSIMD(f3_of_xx2__D2, hDD02), FDPart3tmp166));
        const REAL_SIMD_ARRAY FDPart3tmp188 =
            FusedMulAddSIMD(FDPart3tmp101, hDD_dD011,
                            FusedMulAddSIMD(f0_of_xx0, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(f0_of_xx0__D0, hDD11)),
                                            NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0, f0_of_xx0__D0), FDPart3tmp187)));
        const REAL_SIMD_ARRAY FDPart3tmp234 = FusedMulAddSIMD(FDPart3tmp50, FDPart3tmp65, FDPart3_NegativeOne_);
        const REAL_SIMD_ARRAY FDPart3tmp285 = MulSIMD(FDPart3tmp13, FDPart3tmp37);
        const REAL_SIMD_ARRAY FDPart3tmp289 = FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp58, MulSIMD(FDPart3tmp23, hDD_dDD0201));
        const REAL_SIMD_ARRAY FDPart3tmp297 = FusedMulAddSIMD(FDPart3tmp29, hDD_dD020, MulSIMD(FDPart3tmp296, FDPart3tmp58));
        const REAL_SIMD_ARRAY FDPart3tmp304 = FusedMulAddSIMD(FDPart3tmp31, FDPart3tmp73, MulSIMD(FDPart3tmp58, FDPart3tmp84));
        const REAL_SIMD_ARRAY FDPart3tmp328 = NegFusedMulAddSIMD(FDPart3tmp50, FDPart3tmp65, FDPart3_Integer_1);
        const REAL_SIMD_ARRAY FDPart3tmp22 =
            DivSIMD(FDPart3_Integer_1,
                    FusedMulAddSIMD(
                        FDPart3tmp19, FDPart3tmp20,
                        FusedMulAddSIMD(FDPart3tmp13, MulSIMD(FDPart3tmp19, FDPart3tmp4),
                                        FusedMulAddSIMD(FDPart3tmp2, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp10), MulSIMD(FDPart3tmp9, hDD01)),
                                                        FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp15, MulSIMD(FDPart3tmp17, FDPart3tmp4))))));
        const REAL_SIMD_ARRAY FDPart3tmp55 = SubSIMD(NegFusedMulAddSIMD(FDPart3tmp2, hDD01, FDPart3tmp54), MulSIMD(FDPart3tmp50, hDD01));
        const REAL_SIMD_ARRAY FDPart3tmp60 = FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp19, FDPart3tmp17);
        const REAL_SIMD_ARRAY FDPart3tmp27 = MulSIMD(FDPart3tmp22, FDPart3tmp26);
        const REAL_SIMD_ARRAY FDPart3tmp33 = MulSIMD(FDPart3tmp22, FDPart3tmp32);
        const REAL_SIMD_ARRAY FDPart3tmp78 = MulSIMD(FDPart3tmp22, FDPart3tmp60);
        const REAL_SIMD_ARRAY FDPart3tmp86 = MulSIMD(FDPart3tmp22, FDPart3tmp49);
        const REAL_SIMD_ARRAY FDPart3tmp90 = MulSIMD(FDPart3tmp22, FDPart3tmp40);
        const REAL_SIMD_ARRAY FDPart3tmp93 = MulSIMD(FDPart3tmp22, FDPart3tmp39);
        const REAL_SIMD_ARRAY FDPart3tmp76 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp27);
        const REAL_SIMD_ARRAY FDPart3tmp79 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp87 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp86);
        const REAL_SIMD_ARRAY FDPart3tmp89 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp33);
        const REAL_SIMD_ARRAY FDPart3tmp91 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp90);
        const REAL_SIMD_ARRAY FDPart3tmp94 = MulSIMD(FDPart3_Rational_1_2, FDPart3tmp93);
        const REAL_SIMD_ARRAY FDPart3tmp100 = MulSIMD(FDPart3_Integer_3, FDPart3tmp93);
        const REAL_SIMD_ARRAY FDPart3tmp108 = MulSIMD(FDPart3_Integer_3, FDPart3tmp33);
        const REAL_SIMD_ARRAY FDPart3tmp118 = MulSIMD(FDPart3_Integer_3, FDPart3tmp27);
        const REAL_SIMD_ARRAY FDPart3tmp126 = MulSIMD(FDPart3_Integer_3, FDPart3tmp90);
        const REAL_SIMD_ARRAY FDPart3tmp128 = MulSIMD(FDPart3_Integer_3, FDPart3tmp86);
        const REAL_SIMD_ARRAY FDPart3tmp135 = MulSIMD(FDPart3_Integer_3, FDPart3tmp78);
        const REAL_SIMD_ARRAY FDPart3tmp88 =
            FusedMulAddSIMD(FDPart3tmp79, FDPart3tmp80, FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp87, MulSIMD(FDPart3tmp76, FDPart3tmp77)));
        const REAL_SIMD_ARRAY FDPart3tmp92 =
            FusedMulAddSIMD(FDPart3tmp80, FDPart3tmp87, FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp91, MulSIMD(FDPart3tmp77, FDPart3tmp89)));
        const REAL_SIMD_ARRAY FDPart3tmp95 =
            FusedMulAddSIMD(FDPart3tmp77, FDPart3tmp94, FusedMulAddSIMD(FDPart3tmp85, FDPart3tmp89, MulSIMD(FDPart3tmp76, FDPart3tmp80)));
        const REAL_SIMD_ARRAY FDPart3tmp106 =
            FusedMulAddSIMD(FDPart3tmp56, FDPart3tmp79,
                            FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp27, AddSIMD(FDPart3tmp84, SubSIMD(FDPart3tmp82, FDPart3tmp81))),
                                            MulSIMD(FDPart3tmp104, FDPart3tmp86)));
        const REAL_SIMD_ARRAY FDPart3tmp116 = FusedMulAddSIMD(
            FDPart3tmp114, FDPart3tmp87, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp113, FDPart3tmp76, FDPart3tmp110)));
        const REAL_SIMD_ARRAY FDPart3tmp119 =
            FusedMulAddSIMD(FDPart3tmp56, FDPart3tmp76,
                            FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp93, AddSIMD(FDPart3tmp84, SubSIMD(FDPart3tmp82, FDPart3tmp81))),
                                            MulSIMD(FDPart3tmp104, FDPart3tmp33)));
        const REAL_SIMD_ARRAY FDPart3tmp122 =
            FusedMulAddSIMD(FDPart3tmp56, FDPart3tmp87,
                            FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp33, AddSIMD(FDPart3tmp84, SubSIMD(FDPart3tmp82, FDPart3tmp81))),
                                            FusedMulAddSIMD(FDPart3tmp104, FDPart3tmp90, FDPart3tmp58)));
        const REAL_SIMD_ARRAY FDPart3tmp129 =
            FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp91, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp87, MulSIMD(FDPart3tmp113, FDPart3tmp89)));
        const REAL_SIMD_ARRAY FDPart3tmp130 =
            FusedMulAddSIMD(FDPart3tmp114, FDPart3tmp89, FusedMulAddSIMD(FDPart3tmp115, FDPart3tmp76, MulSIMD(FDPart3tmp113, FDPart3tmp94)));
        const REAL_SIMD_ARRAY FDPart3tmp142 =
            FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp91, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp87, MulSIMD(FDPart3tmp138, FDPart3tmp89)));
        const REAL_SIMD_ARRAY FDPart3tmp143 =
            FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp89, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp76, MulSIMD(FDPart3tmp138, FDPart3tmp94)));
        const REAL_SIMD_ARRAY FDPart3tmp144 =
            FusedMulAddSIMD(FDPart3tmp139, FDPart3tmp87, FusedMulAddSIMD(FDPart3tmp141, FDPart3tmp79, MulSIMD(FDPart3tmp138, FDPart3tmp76)));
        const REAL_SIMD_ARRAY FDPart3tmp172 =
            FusedMulAddSIMD(FDPart3tmp167, FDPart3tmp87, FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp91, MulSIMD(FDPart3tmp163, FDPart3tmp89)));
        const REAL_SIMD_ARRAY FDPart3tmp173 =
            FusedMulAddSIMD(FDPart3tmp167, FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp87, MulSIMD(FDPart3tmp163, FDPart3tmp76)));
        const REAL_SIMD_ARRAY FDPart3tmp175 = FusedMulAddSIMD(
            FDPart3tmp167, FDPart3tmp76, FusedMulAddSIMD(FDPart3tmp171, FDPart3tmp89, FusedMulSubSIMD(FDPart3tmp163, FDPart3tmp94, FDPart3tmp38)));
        const REAL_SIMD_ARRAY FDPart3tmp189 =
            FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp91, FusedMulAddSIMD(FDPart3tmp188, FDPart3tmp87, MulSIMD(FDPart3tmp185, FDPart3tmp89)));
        const REAL_SIMD_ARRAY FDPart3tmp190 =
            FusedMulAddSIMD(FDPart3tmp186, FDPart3tmp89, FusedMulAddSIMD(FDPart3tmp188, FDPart3tmp76, MulSIMD(FDPart3tmp185, FDPart3tmp94)));
        const REAL_SIMD_ARRAY FDPart3tmp192 = FusedMulAddSIMD(
            FDPart3tmp186, FDPart3tmp87, FusedMulAddSIMD(FDPart3tmp188, FDPart3tmp79, FusedMulAddSIMD(FDPart3tmp185, FDPart3tmp76, FDPart3tmp191)));
        const REAL_SIMD_ARRAY FDPart3tmp136 =
            FusedMulAddSIMD(FDPart3tmp31, FDPart3tmp95, FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp88, MulSIMD(FDPart3tmp19, FDPart3tmp92)));
        const REAL_SIMD_ARRAY FDPart3tmp145 =
            FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp25, FusedMulAddSIMD(FDPart3tmp144, FDPart3tmp4, MulSIMD(FDPart3tmp142, FDPart3tmp47)));
        const REAL_SIMD_ARRAY FDPart3tmp150 =
            FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp13, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp31, MulSIMD(FDPart3tmp106, FDPart3tmp25)));
        const REAL_SIMD_ARRAY FDPart3tmp153 = MulSIMD(FDPart3_Integer_2, FDPart3tmp129);
        const REAL_SIMD_ARRAY FDPart3tmp176 =
            FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp4, FusedMulAddSIMD(FDPart3tmp175, FDPart3tmp25, MulSIMD(FDPart3tmp172, FDPart3tmp47)));
        const REAL_SIMD_ARRAY FDPart3tmp193 =
            FusedMulAddSIMD(FDPart3tmp190, FDPart3tmp25, FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp4, MulSIMD(FDPart3tmp189, FDPart3tmp47)));
        const REAL_SIMD_ARRAY FDPart3tmp196 =
            FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp31, FusedMulAddSIMD(FDPart3tmp13, FDPart3tmp130, MulSIMD(FDPart3tmp116, FDPart3tmp25)));
        const REAL_SIMD_ARRAY FDPart3tmp199 = MulSIMD(FDPart3_Integer_2, FDPart3tmp130);
        const REAL_SIMD_ARRAY FDPart3tmp203 =
            FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp19, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp31, MulSIMD(FDPart3tmp116, FDPart3tmp47)));
        const REAL_SIMD_ARRAY FDPart3tmp215 = FusedMulAddSIMD(
            FDPart3tmp144, FDPart3tmp33,
            FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp94,
                            FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp91,
                                            FusedMulAddSIMD(FDPart3tmp27, FDPart3tmp88,
                                                            FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp86, MulSIMD(FDPart3tmp116, FDPart3tmp79))))));
        const REAL_SIMD_ARRAY FDPart3tmp257 =
            FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp31, FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp25, MulSIMD(FDPart3tmp13, FDPart3tmp190)));
        const REAL_SIMD_ARRAY FDPart3tmp269 =
            FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp47, FusedMulAddSIMD(FDPart3tmp175, FDPart3tmp31, MulSIMD(FDPart3tmp172, FDPart3tmp19)));
        const REAL_SIMD_ARRAY FDPart3tmp337 = MulSIMD(FDPart3_Integer_2, FDPart3tmp192);
        const REAL_SIMD_ARRAY FDPart3tmp340 = MulSIMD(FDPart3_Integer_2, FDPart3tmp190);
        const REAL_SIMD_ARRAY FDPart3tmp363 = MulSIMD(FDPart3_Integer_2, FDPart3tmp172);
        const REAL_SIMD_ARRAY FDPart3tmp364 = MulSIMD(FDPart3_Integer_2, FDPart3tmp173);
        const REAL_SIMD_ARRAY FDPart3tmp98 =
            FusedMulAddSIMD(FDPart3tmp4, FDPart3tmp88, FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp92, MulSIMD(FDPart3tmp25, FDPart3tmp95)));
        const REAL_SIMD_ARRAY FDPart3tmp123 =
            FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp25, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp47, MulSIMD(FDPart3tmp106, FDPart3tmp4)));
        const REAL_SIMD_ARRAY FDPart3tmp132 =
            FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp47, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp25, MulSIMD(FDPart3tmp116, FDPart3tmp4)));
        const REAL_SIMD_ARRAY FDPart3tmp137 = MulSIMD(FDPart3tmp136, FDPart3tmp92);
        const REAL_SIMD_ARRAY FDPart3tmp155 = MulSIMD(FDPart3tmp119, FDPart3tmp150);
        const REAL_SIMD_ARRAY FDPart3tmp159 =
            FusedMulAddSIMD(FDPart3tmp25, FDPart3tmp88, FusedMulAddSIMD(FDPart3tmp31, FDPart3tmp92, MulSIMD(FDPart3tmp13, FDPart3tmp95)));
        const REAL_SIMD_ARRAY FDPart3tmp180 =
            FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp31, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp19, MulSIMD(FDPart3tmp106, FDPart3tmp47)));
        const REAL_SIMD_ARRAY FDPart3tmp212 = FusedMulAddSIMD(
            FDPart3tmp143, FDPart3tmp33,
            FusedMulAddSIMD(FDPart3tmp175, FDPart3tmp94,
                            FusedMulAddSIMD(FDPart3tmp190, FDPart3tmp91,
                                            FusedMulAddSIMD(FDPart3tmp27, FDPart3tmp95,
                                                            FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp86, MulSIMD(FDPart3tmp130, FDPart3tmp79))))));
        const REAL_SIMD_ARRAY FDPart3tmp214 = FusedMulAddSIMD(
            FDPart3tmp142, FDPart3tmp33,
            FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp94,
                            FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp91,
                                            FusedMulAddSIMD(FDPart3tmp27, FDPart3tmp92,
                                                            FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp86, MulSIMD(FDPart3tmp129, FDPart3tmp79))))));
        const REAL_SIMD_ARRAY FDPart3tmp238 =
            FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp31, FusedMulAddSIMD(FDPart3tmp144, FDPart3tmp47, MulSIMD(FDPart3tmp142, FDPart3tmp19)));
        const REAL_SIMD_ARRAY FDPart3tmp245 =
            FusedMulAddSIMD(FDPart3tmp190, FDPart3tmp31, FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp47, MulSIMD(FDPart3tmp189, FDPart3tmp19)));
        const REAL_SIMD_ARRAY FDPart3tmp249 = MulSIMD(FDPart3tmp145, FDPart3tmp88);
        const REAL_SIMD_ARRAY FDPart3tmp253 = MulSIMD(FDPart3tmp106, FDPart3tmp145);
        const REAL_SIMD_ARRAY FDPart3tmp267 =
            FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp31, FusedMulAddSIMD(FDPart3tmp144, FDPart3tmp25, MulSIMD(FDPart3tmp13, FDPart3tmp143)));
        const REAL_SIMD_ARRAY FDPart3tmp308 =
            FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp31, FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp25, MulSIMD(FDPart3tmp13, FDPart3tmp175)));
        const REAL_SIMD_ARRAY FDPart3tmp335 = MulSIMD(FDPart3tmp144, FDPart3tmp145);
        const REAL_SIMD_ARRAY FDPart3tmp99 = MulSIMD(FDPart3tmp88, FDPart3tmp98);
        const REAL_SIMD_ARRAY FDPart3tmp107 = MulSIMD(FDPart3tmp106, FDPart3tmp98);
        const REAL_SIMD_ARRAY FDPart3tmp124 = MulSIMD(FDPart3tmp123, FDPart3tmp88);
        const REAL_SIMD_ARRAY FDPart3tmp125 = MulSIMD(FDPart3tmp106, FDPart3tmp123);
        const REAL_SIMD_ARRAY FDPart3tmp160 = MulSIMD(FDPart3tmp159, FDPart3tmp95);
        const REAL_SIMD_ARRAY FDPart3tmp181 = MulSIMD(FDPart3tmp180, FDPart3tmp92);
        const REAL_SIMD_ARRAY FDPart3tmp205 = MulSIMD(FDPart3tmp122, FDPart3tmp180);
        const REAL_SIMD_ARRAY FDPart3tmp239 = MulSIMD(FDPart3tmp238, FDPart3tmp92);
        const REAL_SIMD_ARRAY FDPart3tmp240 = MulSIMD(FDPart3tmp122, FDPart3tmp238);
        const REAL_SIMD_ARRAY FDPart3tmp248 = MulSIMD(FDPart3tmp144, FDPart3tmp98);
        const REAL_SIMD_ARRAY FDPart3tmp263 = MulSIMD(FDPart3tmp123, FDPart3tmp144);
        const REAL_SIMD_ARRAY FDPart3tmp268 = MulSIMD(FDPart3tmp267, FDPart3tmp95);
        const REAL_SIMD_ARRAY FDPart3tmp275 = MulSIMD(FDPart3tmp119, FDPart3tmp267);
        const REAL_SIMD_ARRAY FDPart3tmp309 = MulSIMD(FDPart3_Integer_2, FDPart3tmp308);
        const REAL_SIMD_ARRAY FDPart3tmp311 = FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp159, MulSIMD(FDPart3tmp106, FDPart3tmp176));
        const REAL_SIMD_ARRAY FDPart3tmp319 = FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp150, FDPart3tmp253);
        const REAL_SIMD_ARRAY FDPart3tmp320 = FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp136, MulSIMD(FDPart3tmp122, FDPart3tmp150));
        const REAL_SIMD_ARRAY FDPart3tmp321 = FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp267, MulSIMD(FDPart3tmp122, FDPart3tmp269));
        const REAL_SIMD_ARRAY FDPart3tmp330 = MulSIMD(FDPart3tmp142, FDPart3tmp238);
        const REAL_SIMD_ARRAY FDPart3tmp341 = MulSIMD(FDPart3tmp143, FDPart3tmp267);
        const REAL_SIMD_ARRAY FDPart3tmp305 = MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp119, FDPart3tmp159));
        const REAL_SIMD_ARRAY FDPart3tmp323 = FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp257, FDPart3tmp240);
        const REAL_SIMD_ARRAY FDPart3tmp324 = FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp196, FDPart3tmp107);
        const REAL_SIMD_ARRAY __RHS_exp_0 = FusedMulAddSIMD(
            FDPart3tmp22,
            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp49,
                            FusedMulAddSIMD(FDPart3tmp68, hDD_dD001,
                                            AddSIMD(FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp73, FDPart3tmp75),
                                                    FusedMulAddSIMD(FDPart3tmp5,
                                                                    MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_),
                                                                            MulSIMD(FDPart3tmp57, f0_of_xx0__DD00)),
                                                                    NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp43),
                                                                                       MulSIMD(FDPart3tmp54, f0_of_xx0__D0), FDPart3tmp59)))))),
            FusedMulAddSIMD(
                FDPart3tmp22,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp40, FusedMulAddSIMD(FDPart3tmp2, hDD_dDD0011,
                                                              FusedMulAddSIMD(FDPart3tmp41, hDD_dD000,
                                                                              NegFusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, hDD_dD011),
                                                                                                 FDPart3tmp46))))),
                FusedMulAddSIMD(
                    FDPart3tmp22,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp49, NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp43), MulSIMD(FDPart3tmp55, f0_of_xx0__D0),
                                                                     FDPart3tmp59))),
                    NegFusedMulAddSIMD(
                        FDPart3tmp2, MulSIMD(FDPart3tmp27, hDD_dDD0002),
                        FusedMulAddSIMD(
                            FDPart3tmp22,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp39, FusedMulAddSIMD(FDPart3tmp2, hDD_dDD0022, MulSIMD(FDPart3tmp36, FDPart3tmp38)))),
                            FusedMulAddSIMD(
                                FDPart3tmp4, MulSIMD(FDPart3tmp5, lambdaU_dD00),
                                FusedMulAddSIMD(
                                    f0_of_xx0__D0, MulSIMD(hDD01, lambdaU_dD10),
                                    FusedMulAddSIMD(
                                        FDPart3tmp116, MulSIMD(FDPart3tmp132, FDPart3tmp135),
                                        FusedMulAddSIMD(
                                            FDPart3tmp118, MulSIMD(FDPart3tmp132, FDPart3tmp88),
                                            FusedMulAddSIMD(
                                                FDPart3tmp116, MulSIMD(FDPart3tmp118, FDPart3tmp98),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp116, MulSIMD(FDPart3tmp123, FDPart3tmp128),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp93,
                                                        FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp176, FDPart3tmp95), FDPart3tmp160),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp106, MulSIMD(FDPart3tmp128, FDPart3tmp132),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp90,
                                                                FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp122, FDPart3tmp193),
                                                                                FDPart3tmp205),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp93,
                                                                    FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp145, FDPart3tmp92),
                                                                                    FDPart3tmp137),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp86,
                                                                        FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp150,
                                                                                        MulSIMD(FDPart3tmp145, FDPart3tmp199)),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp90,
                                                                            FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp119, FDPart3tmp145),
                                                                                            FDPart3tmp155),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp86,
                                                                                FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp203,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp122, FDPart3tmp123))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp86,
                                                                                    FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp180,
                                                                                                    MulSIMD(FDPart3tmp153, FDPart3tmp193)),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp78,
                                                                                        FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp196,
                                                                                                        MulSIMD(FDPart3tmp199, FDPart3tmp98)),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp86,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp119, FDPart3tmp196,
                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                        MulSIMD(FDPart3tmp119, FDPart3tmp98))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp33,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp150, FDPart3tmp95,
                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp145, FDPart3tmp95))),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp78,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp123, FDPart3tmp153,
                                                                                                        MulSIMD(FDPart3tmp129, FDPart3tmp203)),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp33,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp119, FDPart3tmp159,
                                                                                                            MulSIMD(FDPart3_Integer_2,
                                                                                                                    MulSIMD(FDPart3tmp119,
                                                                                                                            FDPart3tmp176))),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp33,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp122, FDPart3tmp136,
                                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                                        MulSIMD(FDPart3tmp122,
                                                                                                                                FDPart3tmp145))),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp27,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp203, FDPart3tmp92,
                                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                                            MulSIMD(FDPart3tmp123,
                                                                                                                                    FDPart3tmp92))),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp33,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3_Integer_2,
                                                                                                                        MulSIMD(FDPart3tmp193,
                                                                                                                                FDPart3tmp92),
                                                                                                                        FDPart3tmp181),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp27,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp130,
                                                                                                                            FDPart3tmp159,
                                                                                                                            MulSIMD(FDPart3tmp176,
                                                                                                                                    FDPart3tmp199)),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp27,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp196,
                                                                                                                                FDPart3tmp95,
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp95,
                                                                                                                                        FDPart3tmp98))),
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp215,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp199,
                                                                                                                                    FDPart3tmp25,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp116,
                                                                                                                                            FDPart3tmp4),
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp129,
                                                                                                                                            FDPart3tmp209))),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp27,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp129,
                                                                                                                                        FDPart3tmp136,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp145,
                                                                                                                                            FDPart3tmp153)),
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp212,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3_Integer_2,
                                                                                                                                            MulSIMD(
                                                                                                                                                FDPart3tmp25,
                                                                                                                                                FDPart3tmp95),
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3_Integer_2,
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp4,
                                                                                                                                                    FDPart3tmp88),
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp209,
                                                                                                                                                    FDPart3tmp92))),
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp214,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3_Integer_2,
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp106,
                                                                                                                                                    FDPart3tmp4),
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                                    MulSIMD(
                                                                                                                                                        FDPart3tmp119,
                                                                                                                                                        FDPart3tmp25),
                                                                                                                                                    MulSIMD(
                                                                                                                                                        FDPart3tmp122,
                                                                                                                                                        FDPart3tmp209))),
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp108,
                                                                                                                                                FDPart3tmp124,
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3tmp125,
                                                                                                                                                    FDPart3tmp126,
                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                        FDPart3tmp100,
                                                                                                                                                        FDPart3tmp99,
                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                            FDPart3tmp107,
                                                                                                                                                            FDPart3tmp108,
                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp22,
                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp60, FusedMulAddSIMD(f0_of_xx0__D0,
                                                                                                                                                                                                                                                                           MulSIMD(f0_of_xx0__DD00, hDD_dD000), FusedMulAddSIMD(FDPart3tmp3, FDPart3tmp67, FusedMulAddSIMD(hDD00, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f0_of_xx0__DD00, f0_of_xx0__DD00), MulSIMD(FDPart3tmp62, f0_of_xx0__DDD000)), FusedMulAddSIMD(FDPart3tmp5, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp69, f0_of_xx0__DD00)), FusedMulSubSIMD(FDPart3tmp2, hDD_dDD0000, MulSIMD(FDPart3tmp3, FusedMulSubSIMD(FDPart3tmp5, f0_of_xx0__DDD000, FDPart3tmp66))))))))),
                                                                                                                                                                            FusedMulSubSIMD(FDPart3tmp1, lambdaU_dD20, MulSIMD(FDPart3tmp33, FusedMulSubSIMD(FDPart3tmp2, hDD_dDD0012, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp2, hDD_dD012)))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_1 = FusedMulAddSIMD(FDPart3tmp22, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp49, SubSIMD(FusedMulAddSIMD(FDPart3tmp41, hDD_dD000, FusedMulAddSIMD(FDPart3tmp45, FDPart3tmp58, SubSIMD(NegFusedMulAddSIMD(FDPart3tmp2, hDD_dD011, FDPart3tmp227), MulSIMD(FDPart3tmp50, hDD_dD011)))), MulSIMD(FDPart3tmp41, hDD_dD110)))), FusedMulAddSIMD(FDPart3tmp22, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp49, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp45, FusedMulAddSIMD(FDPart3tmp18, FDPart3tmp73, AddSIMD(FusedMulAddSIMD(FDPart3tmp103, FDPart3tmp58, FDPart3tmp46), FusedMulAddSIMD(FDPart3tmp191, FDPart3tmp69, NegFusedMulAddSIMD(FDPart3tmp234, FDPart3tmp3, FDPart3tmp227))))))), FusedMulAddSIMD(FDPart3tmp22, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp39, FusedMulSubSIMD(FDPart3tmp41, hDD_dDD0122, MulSIMD(FDPart3tmp38, FDPart3tmp81)))), FusedMulAddSIMD(FDPart3tmp22, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp40, FusedMulAddSIMD(FDPart3tmp41, hDD_dD001, FusedMulAddSIMD(FDPart3tmp191, FDPart3tmp57, FusedMulAddSIMD(FDPart3tmp229, FDPart3tmp58, FusedMulAddSIMD(FDPart3tmp41, hDD_dDD0111, FusedMulSubSIMD(FDPart3tmp191, FDPart3tmp55, MulSIMD(FDPart3tmp41, hDD_dD111)))))))), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, FDPart3tmp19), MulSIMD(FDPart3tmp43, lambdaU_dD10), FusedMulAddSIMD(MulSIMD(FDPart3_Rational_1_2, f0_of_xx0), MulSIMD(hDD01, lambdaU_dD00), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp222, FDPart3tmp47), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp224, FDPart3tmp4), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp1, lambdaU_dD21), FusedMulAddSIMD(FDPart3_Rational_1_2, MulSIMD(FDPart3tmp169, lambdaU_dD20), FusedMulAddSIMD(FDPart3tmp93, FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp176, FusedMulAddSIMD(FDPart3tmp269, FDPart3tmp95, FDPart3tmp268)), SubSIMD(FusedMulAddSIMD(FDPart3tmp93, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp239, MulSIMD(FDPart3tmp142, FDPart3tmp145)), FusedMulAddSIMD(FDPart3tmp93, AddSIMD(FDPart3tmp248, FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp88, FDPart3tmp249)), FusedMulAddSIMD(FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp193, FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp192, MulSIMD(FDPart3tmp106, FDPart3tmp180))), FusedMulAddSIMD(FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp257, FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp190, MulSIMD(FDPart3tmp119, FDPart3tmp238))), FusedMulAddSIMD(FDPart3tmp86, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp238, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp257, MulSIMD(FDPart3tmp119, FDPart3tmp145))), FusedMulAddSIMD(FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp193, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp122, FDPart3tmp245))), FusedMulAddSIMD(FDPart3tmp86, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp180, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp193, FDPart3tmp125)), FusedMulAddSIMD(FDPart3tmp86, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp136, FusedMulAddSIMD(FDPart3tmp190, FDPart3tmp98, FDPart3tmp155)), FusedMulAddSIMD(FDPart3tmp86, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp193, MulSIMD(FDPart3tmp153, FDPart3tmp245)), FusedMulAddSIMD(FDPart3tmp86, FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp203, FusedMulAddSIMD(FDPart3tmp132, FDPart3tmp192, FDPart3tmp125)), FusedMulAddSIMD(FDPart3tmp78, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp136, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp150, MulSIMD(FDPart3tmp119, FDPart3tmp98))), FusedMulAddSIMD(FDPart3tmp86, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp205, MulSIMD(FDPart3tmp123, FDPart3tmp189)), FusedMulAddSIMD(FDPart3tmp78, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp123, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp129, FDPart3tmp180))), FusedMulAddSIMD(FDPart3tmp78, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp123, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp203, MulSIMD(FDPart3tmp106, FDPart3tmp132))), FusedMulAddSIMD(FDPart3tmp33, FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp269, FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp190, FDPart3tmp275)), FusedMulAddSIMD(FDPart3tmp33, FusedMulAddSIMD(FDPart3tmp238, FDPart3tmp95, FusedMulAddSIMD(FDPart3tmp257, FDPart3tmp95, MulSIMD(FDPart3tmp143, FDPart3tmp145))), FusedMulAddSIMD(FDPart3tmp33, FusedMulAddSIMD(FDPart3tmp106, FDPart3tmp136, FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp98, FDPart3tmp253)), FusedMulAddSIMD(FDPart3tmp33, FusedMulAddSIMD(FDPart3tmp180, FDPart3tmp88, FusedMulAddSIMD(FDPart3tmp193, FDPart3tmp88, FDPart3tmp263)), FusedMulAddSIMD(FDPart3tmp33, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp240, MulSIMD(FDPart3tmp145, FDPart3tmp189)), FusedMulAddSIMD(FDPart3tmp33, FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp193, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp245, FDPart3tmp92))), FusedMulAddSIMD(FDPart3tmp27, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp267, FusedMulAddSIMD(FDPart3tmp130, FDPart3tmp269, MulSIMD(FDPart3tmp119, FDPart3tmp176))), FusedMulAddSIMD(FDPart3tmp27, FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp98, FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp95, MulSIMD(FDPart3tmp136, FDPart3tmp95))), FusedMulAddSIMD(FDPart3tmp27, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp136, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp145, FDPart3tmp107)), FusedMulAddSIMD(FDPart3tmp27, FusedMulAddSIMD(FDPart3tmp132, FDPart3tmp144, FusedMulAddSIMD(FDPart3tmp203, FDPart3tmp88, FDPart3tmp124)), FusedMulAddSIMD(FDPart3tmp27, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp181, MulSIMD(FDPart3tmp123, FDPart3tmp142)), FusedMulAddSIMD(FDPart3tmp27, FusedMulAddSIMD(FDPart3tmp122, FDPart3tmp145, MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp129, FDPart3tmp238))), FusedMulAddSIMD(FDPart3tmp214, AddSIMD(FDPart3tmp180, FDPart3tmp193), FusedMulAddSIMD(FDPart3tmp215, AddSIMD(FDPart3tmp123, FDPart3tmp203), FusedMulAddSIMD(FDPart3tmp22, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp60, FusedMulAddSIMD(FDPart3tmp55, FDPart3tmp58, FusedMulAddSIMD(hDD01, FusedMulAddSIMD(f0_of_xx0, f0_of_xx0__DDD000, MulSIMD(FDPart3_Integer_3, MulSIMD(f0_of_xx0__D0, f0_of_xx0__DD00))), FusedMulAddSIMD(FDPart3tmp47, FDPart3tmp67, FusedMulAddSIMD(FDPart3tmp54, FDPart3tmp58, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp54, FusedMulAddSIMD(FDPart3tmp41, hDD_dDD0100, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp52, hDD_dD010), NegFusedMulAddSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp5), MulSIMD(FDPart3tmp55, f0_of_xx0__DD00), FDPart3tmp75)))))))))), FusedMulSubSIMD(FDPart3tmp212, AddSIMD(FDPart3tmp136, FDPart3tmp145), MulSIMD(FDPart3tmp27, FusedMulAddSIMD(FDPart3tmp52, hDD_dD012, SubSIMD(FusedMulSubSIMD(FDPart3tmp41, hDD_dDD0102, MulSIMD(FDPart3tmp2, hDD_dD012)), MulSIMD(FDPart3tmp50, hDD_dD012)))))))))))))))))))))))))))))))))), MulSIMD(FDPart3tmp33, FusedMulAddSIMD(FDPart3tmp41, hDD_dDD0112, FusedMulSubSIMD(FDPart3tmp41, hDD_dD002, MulSIMD(FDPart3tmp41, hDD_dD112))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_2 = FusedMulAddSIMD(
    FDPart3tmp22,
    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp49, NegFusedMulAddSIMD(FDPart3tmp23, hDD_dD120, FDPart3tmp289))),
    FusedMulAddSIMD(
        FDPart3tmp22,
        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp49, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp140,
                                                      FusedMulAddSIMD(FDPart3tmp290, hDD_dD021, AddSIMD(FDPart3tmp289, FDPart3tmp304))))),
        FusedMulAddSIMD(
            FDPart3tmp22,
            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp39, FusedMulAddSIMD(FDPart3tmp23, hDD_dDD0222,
                                                          SubSIMD(FusedMulSubSIMD(FDPart3tmp1, f3_of_xx2__DD22,
                                                                                  MulSIMD(FDPart3tmp25, FusedMulSubSIMD(FDPart3tmp37, f3_of_xx2__DD22,
                                                                                                                        FDPart3tmp302))),
                                                                  MulSIMD(FDPart3tmp291, FDPart3tmp38))))),
            FusedMulAddSIMD(
                FDPart3tmp22,
                MulSIMD(
                    MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp40, FusedMulAddSIMD(FDPart3tmp23, hDD_dDD0211, NegFusedMulAddSIMD(FDPart3tmp23, hDD_dD121, FDPart3tmp297)))),
                FusedMulAddSIMD(
                    FDPart3tmp22,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp32, FusedMulSubSIMD(FDPart3tmp23, hDD_dDD0212, MulSIMD(FDPart3tmp23, hDD_dD122)))),
                    FusedMulAddSIMD(
                        FDPart3tmp22,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp32, FusedMulAddSIMD(FDPart3tmp293, hDD_dD021,
                                                                      FusedMulAddSIMD(FDPart3tmp298, FDPart3tmp58,
                                                                                      FusedMulSubSIMD(FDPart3tmp23, hDD_dDD0212,
                                                                                                      MulSIMD(FDPart3tmp140, FDPart3tmp38)))))),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3_Rational_1_2, f3_of_xx2), MulSIMD(hDD12, lambdaU_dD10),
                            FusedMulAddSIMD(
                                FDPart3tmp22,
                                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                        MulSIMD(FDPart3tmp26, FusedMulAddSIMD(FDPart3tmp293, hDD_dD020,
                                                                              NegFusedMulAddSIMD(FDPart3tmp294, FDPart3tmp38, FDPart3tmp292)))),
                                FusedMulAddSIMD(
                                    MulSIMD(FDPart3_Rational_1_2, FDPart3tmp4), MulSIMD(FDPart3tmp5, lambdaU_dD02),
                                    FusedMulAddSIMD(
                                        MulSIMD(FDPart3_Rational_1_2, f0_of_xx0__D0), MulSIMD(hDD01, lambdaU_dD12),
                                        FusedMulAddSIMD(
                                            FDPart3_Rational_1_2, MulSIMD(FDPart3tmp282, lambdaU_dD00),
                                            FusedMulAddSIMD(
                                                FDPart3_Rational_1_2, MulSIMD(FDPart3tmp285, lambdaU_dD20),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp93,
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp173, FDPart3tmp98,
                                                        FusedMulAddSIMD(FDPart3tmp176, FDPart3tmp88, MulSIMD(FDPart3tmp159, FDPart3tmp88))),
                                                    FusedMulAddSIMD(
                                                        FDPart3_Rational_1_2, MulSIMD(FDPart3tmp1, lambdaU_dD22),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp93,
                                                            FusedMulAddSIMD(FDPart3tmp175, FDPart3tmp176, MulSIMD(FDPart3tmp309, FDPart3tmp95)),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp93,
                                                                FusedMulAddSIMD(FDPart3tmp267, FDPart3tmp92,
                                                                                FusedMulAddSIMD(FDPart3tmp269, FDPart3tmp92,
                                                                                                MulSIMD(FDPart3tmp145, FDPart3tmp172))),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp90, FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp193, FDPart3tmp323),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp90,
                                                                        FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp275,
                                                                                        MulSIMD(FDPart3tmp143, FDPart3tmp145)),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp86,
                                                                            FusedMulAddSIMD(FDPart3tmp129, FDPart3tmp257,
                                                                                            FusedMulAddSIMD(FDPart3tmp193, FDPart3tmp92,
                                                                                                            MulSIMD(FDPart3tmp129, FDPart3tmp238))),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp90, AddSIMD(FDPart3tmp263, FDPart3tmp319),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp86,
                                                                                    FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp95,
                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp130, FDPart3tmp267))),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp86,
                                                                                        FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp145,
                                                                                                        FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp150,
                                                                                                                        FDPart3tmp124)),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp86,
                                                                                            FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp142,
                                                                                                            FDPart3tmp320),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp86,
                                                                                                FusedMulAddSIMD(FDPart3tmp132, FDPart3tmp144,
                                                                                                                FDPart3tmp324),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp78,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp129, FDPart3tmp136,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp129, FDPart3tmp150,
                                                                                                            MulSIMD(FDPart3tmp123, FDPart3tmp92))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp86,
                                                                                                        FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp98,
                                                                                                                        FDPart3tmp305),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp78,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp95, FDPart3tmp98,
                                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                                        MulSIMD(FDPart3tmp130,
                                                                                                                                FDPart3tmp159))),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp78,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp116, FDPart3tmp98,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp132, FDPart3tmp88,
                                                                                                                        MulSIMD(FDPart3tmp116,
                                                                                                                                FDPart3tmp196))),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp33,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp172, FDPart3tmp193,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp257,
                                                                                                                            FDPart3tmp92,
                                                                                                                            FDPart3tmp239)),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp33,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp123,
                                                                                                                            FDPart3tmp173,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp150,
                                                                                                                                FDPart3tmp88,
                                                                                                                                FDPart3tmp249)),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp33,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                FDPart3tmp268,
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3tmp145,
                                                                                                                                    FDPart3tmp175)),
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp33,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp119,
                                                                                                                                    FDPart3tmp309,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp143,
                                                                                                                                        FDPart3tmp176)),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp33,
                                                                                                                                    AddSIMD(
                                                                                                                                        FDPart3tmp248,
                                                                                                                                        FDPart3tmp311),
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp33,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp142,
                                                                                                                                            FDPart3tmp145,
                                                                                                                                            FDPart3tmp321),
                                                                                                                                        FusedMulAddSIMD(FDPart3tmp27,
                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp132,
                                                                                                                                                                        FDPart3tmp173,
                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp196,
                                                                                                                                                                                        FDPart3tmp88,
                                                                                                                                                                                        FDPart3tmp99)),
                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp27,
                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp129,
                                                                                                                                                                                        FDPart3tmp269,
                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp145,
                                                                                                                                                                                                        FDPart3tmp92,
                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                            FDPart3tmp129,
                                                                                                                                                                                                            FDPart3tmp267))),
                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp27,
                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp123,
                                                                                                                                                                                                        FDPart3tmp172,
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp150,
                                                                                                                                                                                                                        FDPart3tmp92,
                                                                                                                                                                                                                        FDPart3tmp137)),
                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp27,
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp159, FusedMulAddSIMD(FDPart3tmp116, FDPart3tmp176, FDPart3tmp99)),
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp27,
                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                            FDPart3_Integer_2,
                                                                                                                                                                                                                            FDPart3tmp160,
                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                FDPart3tmp175,
                                                                                                                                                                                                                                FDPart3tmp98)),
                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp27,
                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                            FDPart3tmp176, FDPart3tmp95, MulSIMD(FDPart3tmp199, FDPart3tmp308)),
                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                            FDPart3tmp214,
                                                                                                                                                                                                                                            AddSIMD(
                                                                                                                                                                                                                                                FDPart3tmp145,
                                                                                                                                                                                                                                                FDPart3tmp150),
                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp215,
                                                                                                                                                                                                                                                            AddSIMD(
                                                                                                                                                                                                                                                                FDPart3tmp196,
                                                                                                                                                                                                                                                                FDPart3tmp98),
                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp22,
                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                                                                                                    FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                                                                                                                                                                                                                                                                MulSIMD(FDPart3tmp60,
                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp25,
                                                                                                                                                                                                                                                                                                        FDPart3tmp67,
                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp282, f0_of_xx0__DDD000, FusedMulAddSIMD(FDPart3tmp110, FDPart3tmp294, MulSIMD(FDPart3tmp23, hDD_dDD0200)))))),
                                                                                                                                                                                                                                                                            FusedMulSubSIMD(
                                                                                                                                                                                                                                                                                FDPart3tmp212, AddSIMD(FDPart3tmp159, FDPart3tmp176),
                                                                                                                                                                                                                                                                                MulSIMD(
                                                                                                                                                                                                                                                                                    MulSIMD(FDPart3_Rational_1_2, FDPart3tmp22), MulSIMD(FDPart3tmp26,
                                                                                                                                                                                                                                                                                                                                         FusedMulAddSIMD(FDPart3tmp110,
                                                                                                                                                                                                                                                                                                                                                         FDPart3tmp291,
                                                                                                                                                                                                                                                                                                                                                         FusedMulAddSIMD(FDPart3tmp290, hDD_dD022, FDPart3tmp292))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_3 = FusedMulAddSIMD(
    FDPart3tmp22,
    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
            MulSIMD(FDPart3tmp49,
                    FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp191, FDPart3tmp54),
                                    FusedMulAddSIMD(FDPart3tmp101, hDD_dD111,
                                                    FusedMulAddSIMD(FDPart3tmp328, FDPart3tmp47,
                                                                    FusedMulAddSIMD(FDPart3tmp229,
                                                                                    MulSIMD(MulSIMD(FDPart3_Integer_3, FDPart3_NegativeOne_),
                                                                                            MulSIMD(FDPart3tmp43, f0_of_xx0__D0)),
                                                                                    FusedMulSubSIMD(FDPart3tmp10, hDD_dDD1101,
                                                                                                    MulSIMD(FDPart3tmp234, FDPart3tmp47)))))))),
    FusedMulAddSIMD(
        FDPart3tmp22,
        MulSIMD(
            MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
            MulSIMD(FDPart3tmp40, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp191, FDPart3tmp45),
                                                  FusedMulAddSIMD(FDPart3tmp5, MulSIMD(hDD_dD110, MulSIMD(MulSIMD(f0_of_xx0, f0_of_xx0), f0_of_xx0)),
                                                                  FusedMulAddSIMD(FDPart3tmp10, hDD_dDD1111, MulSIMD(FDPart3tmp228, hDD_dD011)))))),
        FusedMulAddSIMD(
            FDPart3tmp22,
            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp49, FusedMulAddSIMD(FDPart3tmp229, FDPart3tmp58,
                                                          FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp191, FDPart3tmp55),
                                                                          MulSIMD(FDPart3tmp10, hDD_dDD1101))))),
            NegFusedMulAddSIMD(
                FDPart3tmp10, MulSIMD(FDPart3tmp27, hDD_dDD1102),
                FusedMulAddSIMD(
                    FDPart3tmp22,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp39, FusedMulAddSIMD(FDPart3tmp10, hDD_dDD1122, MulSIMD(FDPart3tmp184, FDPart3tmp38)))),
                    FusedMulAddSIMD(
                        FDPart3tmp128, MulSIMD(FDPart3tmp180, FDPart3tmp189),
                        FusedMulAddSIMD(
                            MulSIMD(FDPart3_Integer_3, FDPart3tmp122),
                            MulSIMD(FDPart3tmp245, FDPart3tmp86),
                            FusedMulAddSIMD(
                                FDPart3tmp118,
                                MulSIMD(FDPart3tmp142, FDPart3tmp180),
                                FusedMulAddSIMD(
                                    FDPart3tmp126,
                                    MulSIMD(FDPart3tmp189, FDPart3tmp245),
                                    FusedMulAddSIMD(
                                        FDPart3tmp108,
                                        MulSIMD(FDPart3tmp142, FDPart3tmp245),
                                        FusedMulAddSIMD(
                                            FDPart3tmp108,
                                            MulSIMD(FDPart3tmp189, FDPart3tmp238),
                                            FusedMulAddSIMD(
                                                FDPart3tmp93,
                                                FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp136, FDPart3tmp144), FDPart3tmp335),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp93,
                                                    FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp143, FDPart3tmp269), FDPart3tmp341),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp90,
                                                        FusedMulAddSIMD(FDPart3tmp180, FDPart3tmp337, MulSIMD(FDPart3tmp192, FDPart3tmp193)),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp90,
                                                            FusedMulAddSIMD(FDPart3tmp190, FDPart3tmp257, MulSIMD(FDPart3tmp238, FDPart3tmp340)),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp86,
                                                                FusedMulAddSIMD(FDPart3tmp123, FDPart3tmp192, MulSIMD(FDPart3tmp203, FDPart3tmp337)),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp86,
                                                                    FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp340,
                                                                                    MulSIMD(FDPart3tmp150, FDPart3tmp190)),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp86,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp106, FDPart3tmp193,
                                                                            MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp106, FDPart3tmp180))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp86,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp119, FDPart3tmp257,
                                                                                MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp119, FDPart3tmp238))),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp78,
                                                                                FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                MulSIMD(FDPart3tmp106, FDPart3tmp203), FDPart3tmp125),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp78,
                                                                                    FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                    MulSIMD(FDPart3tmp119, FDPart3tmp136),
                                                                                                    FDPart3tmp155),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp33,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp144, FDPart3tmp193,
                                                                                            MulSIMD(FDPart3_Integer_2,
                                                                                                    MulSIMD(FDPart3tmp144, FDPart3tmp180))),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp33,
                                                                                            FusedMulAddSIMD(FDPart3tmp190, FDPart3tmp267,
                                                                                                            MulSIMD(FDPart3tmp269, FDPart3tmp340)),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp33,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp136, FDPart3tmp337,
                                                                                                    MulSIMD(FDPart3tmp145, FDPart3tmp192)),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp33,
                                                                                                    FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp257,
                                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                                            MulSIMD(FDPart3tmp143,
                                                                                                                                    FDPart3tmp238))),
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp27,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp119, FDPart3tmp269),
                                                                                                            FDPart3tmp275),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp27,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp143, FDPart3tmp150,
                                                                                                                MulSIMD(FDPart3_Integer_2,
                                                                                                                        MulSIMD(FDPart3tmp136,
                                                                                                                                FDPart3tmp143))),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp27,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3_Integer_2,
                                                                                                                    MulSIMD(FDPart3tmp106,
                                                                                                                            FDPart3tmp136),
                                                                                                                    FDPart3tmp253),
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp27,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3_Integer_2,
                                                                                                                        MulSIMD(FDPart3tmp144,
                                                                                                                                FDPart3tmp203),
                                                                                                                        FDPart3tmp263),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp215,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3_Integer_2,
                                                                                                                            MulSIMD(FDPart3tmp119,
                                                                                                                                    FDPart3tmp31),
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                MulSIMD(FDPart3tmp122,
                                                                                                                                        FDPart3tmp19),
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3tmp106,
                                                                                                                                    FDPart3tmp209))),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp224,
                                                                                                                            FDPart3tmp47,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp212,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp142,
                                                                                                                                        FDPart3tmp19),
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp143,
                                                                                                                                            FDPart3tmp31),
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp144,
                                                                                                                                            FDPart3tmp209))),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp214,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp31,
                                                                                                                                        FDPart3tmp340,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3_Integer_2,
                                                                                                                                            MulSIMD(
                                                                                                                                                FDPart3tmp189,
                                                                                                                                                FDPart3tmp19),
                                                                                                                                            MulSIMD(
                                                                                                                                                FDPart3tmp192,
                                                                                                                                                FDPart3tmp209))),
                                                                                                                                    FusedMulAddSIMD(FDPart3tmp169,
                                                                                                                                                    lambdaU_dD21,
                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp19,
                                                                                                                                                                    FDPart3tmp222,
                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp118,
                                                                                                                                                                                    FDPart3tmp240,
                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp135,
                                                                                                                                                                                                    FDPart3tmp205,
                                                                                                                                                                                                    FusedMulAddSIMD(FDPart3tmp22, MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2), MulSIMD(FDPart3tmp60, FusedMulAddSIMD(hDD11, FusedMulAddSIMD(FDPart3_Integer_2, FDPart3tmp2, MulSIMD(FDPart3_Integer_2, FDPart3tmp50)), FusedMulAddSIMD(FDPart3tmp187, MulSIMD(FDPart3tmp5, f0_of_xx0__DD00), FusedMulAddSIMD(FDPart3tmp101, hDD_dD110, FusedMulAddSIMD(FDPart3tmp18, FDPart3tmp73, FusedMulAddSIMD(FDPart3tmp103, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp43, f0_of_xx0__D0)), FusedMulSubSIMD(FDPart3tmp10, hDD_dDD1100, MulSIMD(FDPart3tmp18, FusedMulSubSIMD(FDPart3tmp43, f0_of_xx0__DD00, FDPart3tmp72)))))))))),
                                                                                                                                                                                                                    FusedMulSubSIMD(
                                                                                                                                                                                                                        FDPart3tmp100,
                                                                                                                                                                                                                        FDPart3tmp330, MulSIMD(FDPart3tmp33, FusedMulAddSIMD(FDPart3tmp10, hDD_dDD1112, MulSIMD(FDPart3tmp228, hDD_dD012))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_4 =
    FusedMulAddSIMD(
        FDPart3tmp22,
        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp49, FusedMulAddSIMD(FDPart3tmp29, hDD_dDD1201, FDPart3tmp297))),
        FusedMulAddSIMD(
            FDPart3tmp22,
            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                    MulSIMD(FDPart3tmp49,
                            FusedMulAddSIMD(FDPart3tmp23, hDD_dD121,
                                            FusedMulAddSIMD(FDPart3tmp25, FDPart3tmp328,
                                                            FusedMulAddSIMD(FDPart3tmp29, hDD_dDD1201,
                                                                            FusedMulSubSIMD(FDPart3tmp191, FDPart3tmp294,
                                                                                            MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3tmp296),
                                                                                                    MulSIMD(FDPart3tmp43, f0_of_xx0__D0)))))))),
            FusedMulAddSIMD(
                FDPart3tmp22,
                MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                        MulSIMD(FDPart3tmp39,
                                FusedMulAddSIMD(FDPart3tmp29, hDD_dDD1222,
                                                SubSIMD(FusedMulSubSIMD(FDPart3tmp169, f3_of_xx2__DD22, MulSIMD(FDPart3tmp298, FDPart3tmp38)),
                                                        MulSIMD(FDPart3tmp31, FusedMulSubSIMD(FDPart3tmp37, f3_of_xx2__DD22, FDPart3tmp302)))))),
                FusedMulAddSIMD(
                    FDPart3tmp22,
                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                            MulSIMD(FDPart3tmp40, FusedMulAddSIMD(FDPart3tmp29, hDD_dDD1211,
                                                                  FusedMulAddSIMD(MulSIMD(FDPart3tmp10, FDPart3tmp5), MulSIMD(f3_of_xx2, hDD_dD120),
                                                                                  FusedMulAddSIMD(FDPart3tmp140, FDPart3tmp191,
                                                                                                  MulSIMD(FDPart3tmp29, hDD_dD021)))))),
                    FusedMulAddSIMD(
                        FDPart3tmp22,
                        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                MulSIMD(FDPart3tmp32, FusedMulAddSIMD(FDPart3tmp29, hDD_dD022, MulSIMD(FDPart3tmp29, hDD_dDD1212)))),
                        FusedMulAddSIMD(
                            FDPart3tmp22,
                            MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                    MulSIMD(FDPart3tmp32, FusedMulAddSIMD(FDPart3tmp29, hDD_dDD1212,
                                                                          FusedMulAddSIMD(f0_of_xx0, MulSIMD(f3_of_xx2__D2, hDD_dD121),
                                                                                          FusedMulSubSIMD(FDPart3tmp191, FDPart3tmp291,
                                                                                                          MulSIMD(FDPart3tmp296, FDPart3tmp38)))))),
                            FusedMulAddSIMD(
                                MulSIMD(FDPart3_Rational_1_2, f0_of_xx0), MulSIMD(hDD01, lambdaU_dD02),
                                FusedMulAddSIMD(
                                    FDPart3tmp22,
                                    MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                                            MulSIMD(FDPart3tmp26, FusedMulAddSIMD(f0_of_xx0, MulSIMD(f3_of_xx2__D2, hDD_dD120),
                                                                                  NegFusedMulAddSIMD(FDPart3tmp38, FDPart3tmp84, FDPart3tmp348)))),
                                    FusedMulAddSIMD(
                                        FDPart3_Rational_1_2, MulSIMD(FDPart3tmp285, lambdaU_dD21),
                                        FusedMulAddSIMD(
                                            MulSIMD(FDPart3_Rational_1_2, FDPart3tmp19), MulSIMD(FDPart3tmp43, lambdaU_dD12),
                                            FusedMulAddSIMD(
                                                FDPart3_Rational_1_2, MulSIMD(FDPart3tmp222, FDPart3tmp31),
                                                FusedMulAddSIMD(
                                                    FDPart3_Rational_1_2, MulSIMD(FDPart3tmp224, FDPart3tmp25),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp93,
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp142, FDPart3tmp269,
                                                            FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp238, MulSIMD(FDPart3tmp142, FDPart3tmp267))),
                                                        FusedMulAddSIMD(
                                                            FDPart3_Rational_1_2, MulSIMD(FDPart3tmp169, lambdaU_dD22),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp93,
                                                                FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp309, MulSIMD(FDPart3tmp175, FDPart3tmp269)),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp93,
                                                                    FusedMulAddSIMD(FDPart3tmp144, FDPart3tmp159,
                                                                                    FusedMulAddSIMD(FDPart3tmp144, FDPart3tmp176,
                                                                                                    MulSIMD(FDPart3tmp136, FDPart3tmp173))),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp90,
                                                                        FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp238,
                                                                                        FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp257,
                                                                                                        MulSIMD(FDPart3tmp142, FDPart3tmp245))),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp90,
                                                                            FusedMulAddSIMD(FDPart3tmp145, FDPart3tmp192,
                                                                                            FusedMulAddSIMD(FDPart3tmp150, FDPart3tmp192,
                                                                                                            MulSIMD(FDPart3tmp144, FDPart3tmp180))),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp86,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp192, FDPart3tmp196,
                                                                                    FusedMulAddSIMD(FDPart3tmp192, FDPart3tmp98,
                                                                                                    MulSIMD(FDPart3tmp144, FDPart3tmp203))),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp90,
                                                                                    FusedMulAddSIMD(FDPart3tmp143, FDPart3tmp238,
                                                                                                    MulSIMD(FDPart3_Integer_2,
                                                                                                            MulSIMD(FDPart3tmp190, FDPart3tmp267))),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp86,
                                                                                        FusedMulAddSIMD(FDPart3tmp136, FDPart3tmp143,
                                                                                                        MulSIMD(FDPart3tmp159, FDPart3tmp340)),
                                                                                        FusedMulAddSIMD(FDPart3tmp86,
                                                                                                        FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp180,
                                                                                                                        FusedMulAddSIMD(FDPart3tmp150,
                                                                                                                                        FDPart3tmp189,
                                                                                                                                        MulSIMD(FDPart3tmp136,
                                                                                                                                                FDPart3tmp189))),
                                                                                                        FusedMulAddSIMD(FDPart3tmp86,
                                                                                                                        FusedMulAddSIMD(FDPart3tmp245,
                                                                                                                                        FDPart3tmp92,
                                                                                                                                        FDPart3tmp323),
                                                                                                                        FusedMulAddSIMD(FDPart3tmp86,
                                                                                                                                        FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                        FDPart3tmp275,
                                                                                                                                                        MulSIMD(FDPart3tmp238,
                                                                                                                                                                FDPart3tmp95)),
                                                                                                                                        FusedMulAddSIMD(FDPart3tmp78,
                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp203,
                                                                                                                                                                        FDPart3tmp88,
                                                                                                                                                                        FDPart3tmp324),
                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp86,
                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp180,
                                                                                                                                                                                        FDPart3tmp88,
                                                                                                                                                                                        FDPart3tmp319),
                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp78,
                                                                                                                                                                                        AddSIMD(FDPart3tmp181,
                                                                                                                                                                                                FDPart3tmp320),
                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp78,
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp136,
                                                                                                                                                                                                                        FDPart3tmp95,
                                                                                                                                                                                                                        FDPart3tmp305),
                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp33,
                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp144,
                                                                                                                                                                                                                                        FDPart3tmp150,
                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp180, FDPart3tmp335)),
                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp33,
                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                            FDPart3tmp159,
                                                                                                                                                                                                                                            FDPart3tmp192,
                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp176,
                                                                                                                                                                                                                                                            FDPart3tmp192,
                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                FDPart3tmp136,
                                                                                                                                                                                                                                                                FDPart3tmp144))),
                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp33,
                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp142, FDPart3tmp257, FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp245, FDPart3tmp330)),
                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                            FDPart3tmp33,
                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                FDPart3tmp189, FDPart3tmp267, FusedMulAddSIMD(FDPart3tmp189, FDPart3tmp269, FDPart3tmp330)),
                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp33,
                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                            FDPart3tmp341,
                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                FDPart3tmp175,
                                                                                                                                                                                                                                                                                                FDPart3tmp238)),
                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp33,
                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp143,
                                                                                                                                                                                                                                                                                                            FDPart3tmp269,
                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                FDPart3tmp190,
                                                                                                                                                                                                                                                                                                                FDPart3tmp309)),
                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                FDPart3tmp27, FusedMulAddSIMD(FDPart3tmp144, FDPart3tmp196, FusedMulAddSIMD(FDPart3tmp173, FDPart3tmp203, FDPart3tmp248)),
                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                    FDPart3tmp27,
                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                        FDPart3tmp142,
                                                                                                                                                                                                                                                                                                        FDPart3tmp150,
                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp180, MulSIMD(FDPart3tmp136, FDPart3tmp142))),
                                                                                                                                                                                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                        FDPart3tmp27,
                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(FDPart3tmp119, FDPart3tmp309, MulSIMD(FDPart3tmp269, FDPart3tmp95)),
                                                                                                                                                                                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                            FDPart3tmp27,
                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(FDPart3tmp136,
                                                                                                                                                                                                                                                                                                                            FDPart3tmp175,
                                                                                                                                                                                                                                                                                                                            MulSIMD(
                                                                                                                                                                                                                                                                                                                                FDPart3_Integer_2,
                                                                                                                                                                                                                                                                                                                                MulSIMD(FDPart3tmp143, FDPart3tmp159))),
                                                                                                                                                                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                FDPart3tmp27,
                                                                                                                                                                                                                                                                                                                AddSIMD(
                                                                                                                                                                                                                                                                                                                    FDPart3tmp239,
                                                                                                                                                                                                                                                                                                                    FDPart3tmp321),
                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp27,
                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                    FDPart3tmp136, FDPart3tmp88,
                                                                                                                                                                                                                                                                                                                                    FDPart3tmp311),
                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp214,
                                                                                                                                                                                                                                                                                                                                                AddSIMD(
                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp238,
                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp257),
                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp215, AddSIMD(FDPart3tmp136, FDPart3tmp150),
                                                                                                                                                                                                                                                                                                                                                                FusedMulAddSIMD(FDPart3tmp22,
                                                                                                                                                                                                                                                                                                                                                                                MulSIMD(MulSIMD(FDPart3_NegativeOne_,
                                                                                                                                                                                                                                                                                                                                                                                                FDPart3_Rational_1_2),
                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp60, FusedMulAddSIMD(FDPart3tmp23,
                                                                                                                                                                                                                                                                                                                                                                                                                          hDD_dD120,
                                                                                                                                                                                                                                                                                                                                                                                                                          FusedMulAddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                                                              FDPart3tmp29, hDD_dDD1200,
                                                                                                                                                                                                                                                                                                                                                                                                                              FusedMulAddSIMD(FDPart3tmp290,
                                                                                                                                                                                                                                                                                                                                                                                                                                              hDD12,
                                                                                                                                                                                                                                                                                                                                                                                                                                              NegFusedMulAddSIMD(MulSIMD(FDPart3tmp5, FDPart3tmp50), MulSIMD(f3_of_xx2, hDD_dD120), FDPart3tmp304)))))),
                                                                                                                                                                                                                                                                                                                                                                                FusedMulSubSIMD(
                                                                                                                                                                                                                                                                                                                                                                                    FDPart3tmp212,
                                                                                                                                                                                                                                                                                                                                                                                    AddSIMD(
                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp267,
                                                                                                                                                                                                                                                                                                                                                                                        FDPart3tmp269),
                                                                                                                                                                                                                                                                                                                                                                                    MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(
                                                                                                                                                                                                                                                                                                                                                                                            FDPart3_Rational_1_2,
                                                                                                                                                                                                                                                                                                                                                                                            FDPart3tmp22),
                                                                                                                                                                                                                                                                                                                                                                                        MulSIMD(FDPart3tmp26, FusedMulAddSIMD(FDPart3tmp23,
                                                                                                                                                                                                                                                                                                                                                                                                                              hDD_dD122,
                                                                                                                                                                                                                                                                                                                                                                                                                              FusedMulAddSIMD(FDPart3tmp298, FDPart3tmp58, FDPart3tmp348))))))))))))))))))))))))))))))))))))))))))))))));
const REAL_SIMD_ARRAY __RHS_exp_5 = FusedMulAddSIMD(
    FDPart3tmp22,
    MulSIMD(
        MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
        MulSIMD(FDPart3tmp39,
                FusedMulAddSIMD(
                    f3_of_xx2, MulSIMD(f3_of_xx2__D2, hDD_dD222),
                    FusedMulAddSIMD(
                        FDPart3tmp8, hDD_dDD2222,
                        FusedMulAddSIMD(
                            hDD22, FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(f3_of_xx2__D2, f3_of_xx2__D2), MulSIMD(FDPart3tmp112, f3_of_xx2__DD22)),
                            FusedMulAddSIMD(
                                FDPart3tmp162, MulSIMD(MulSIMD(FDPart3_Integer_2, FDPart3_NegativeOne_), MulSIMD(FDPart3tmp37, f3_of_xx2__D2)),
                                FusedMulSubSIMD(FDPart3tmp12, NegFusedMulAddSIMD(FDPart3tmp37, f3_of_xx2__DD22, FDPart3tmp302),
                                                MulSIMD(FDPart3tmp12, FusedMulSubSIMD(FDPart3tmp37, f3_of_xx2__DD22, FDPart3tmp302))))))))),
    FusedMulAddSIMD(
        FDPart3tmp22,
        MulSIMD(MulSIMD(FDPart3_NegativeOne_, FDPart3_Rational_1_2),
                MulSIMD(FDPart3tmp40, FusedMulAddSIMD(FDPart3tmp191, FDPart3tmp77, MulSIMD(FDPart3tmp8, hDD_dDD2211)))),
        NegFusedMulAddSIMD(
            FDPart3tmp27, MulSIMD(FDPart3tmp8, hDD_dDD2202),
            NegFusedMulAddSIMD(
                FDPart3tmp33, MulSIMD(FDPart3tmp8, hDD_dDD2212),
                FusedMulAddSIMD(
                    MulSIMD(FDPart3_Integer_3, FDPart3tmp143), MulSIMD(FDPart3tmp308, FDPart3tmp33),
                    FusedMulAddSIMD(
                        MulSIMD(FDPart3_Integer_3, FDPart3tmp27),
                        MulSIMD(FDPart3tmp308, FDPart3tmp95),
                        FusedMulAddSIMD(
                            FDPart3tmp128,
                            MulSIMD(FDPart3tmp143, FDPart3tmp159),
                            FusedMulAddSIMD(
                                f3_of_xx2, MulSIMD(hDD12, lambdaU_dD12),
                                FusedMulAddSIMD(
                                    FDPart3tmp108,
                                    MulSIMD(FDPart3tmp175, FDPart3tmp267),
                                    FusedMulAddSIMD(
                                        FDPart3tmp118,
                                        MulSIMD(FDPart3tmp159, FDPart3tmp175),
                                        FusedMulAddSIMD(
                                            FDPart3tmp93,
                                            FusedMulAddSIMD(
                                                FDPart3tmp173,
                                                FDPart3tmp176,
                                                MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp159, FDPart3tmp173))),
                                            FusedMulAddSIMD(
                                                FDPart3tmp100,
                                                MulSIMD(FDPart3tmp175, FDPart3tmp308),
                                                FusedMulAddSIMD(
                                                    FDPart3tmp90,
                                                    FusedMulAddSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp144, FDPart3tmp150), FDPart3tmp335),
                                                    FusedMulAddSIMD(
                                                        FDPart3tmp93,
                                                        FusedMulAddSIMD(FDPart3tmp172, FDPart3tmp269, MulSIMD(FDPart3tmp267, FDPart3tmp363)),
                                                        FusedMulAddSIMD(
                                                            FDPart3tmp86,
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp136,
                                                                FDPart3tmp142,
                                                                MulSIMD(FDPart3_Integer_2, MulSIMD(FDPart3tmp142, FDPart3tmp150))),
                                                            FusedMulAddSIMD(
                                                                FDPart3tmp90,
                                                                FusedMulAddSIMD(
                                                                    FDPart3_Integer_2, MulSIMD(FDPart3tmp142, FDPart3tmp257), FDPart3tmp330),
                                                                FusedMulAddSIMD(
                                                                    FDPart3tmp86,
                                                                    FusedMulAddSIMD(
                                                                        FDPart3_Integer_2, MulSIMD(FDPart3tmp144, FDPart3tmp196), FDPart3tmp248),
                                                                    FusedMulAddSIMD(
                                                                        FDPart3tmp86,
                                                                        FusedMulAddSIMD(
                                                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp150, FDPart3tmp88), FDPart3tmp249),
                                                                        FusedMulAddSIMD(
                                                                            FDPart3tmp78,
                                                                            FusedMulAddSIMD(
                                                                                FDPart3_Integer_2, MulSIMD(FDPart3tmp196, FDPart3tmp88),
                                                                                FDPart3tmp99),
                                                                            FusedMulAddSIMD(
                                                                                FDPart3tmp86,
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3_Integer_2, MulSIMD(FDPart3tmp257, FDPart3tmp92),
                                                                                    FDPart3tmp239),
                                                                                FusedMulAddSIMD(
                                                                                    FDPart3tmp33,
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp172,
                                                                                        FDPart3tmp238,
                                                                                        MulSIMD(FDPart3tmp257, FDPart3tmp363)),
                                                                                    FusedMulAddSIMD(
                                                                                        FDPart3tmp78,
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3_Integer_2, MulSIMD(FDPart3tmp150, FDPart3tmp92),
                                                                                            FDPart3tmp137),
                                                                                        FusedMulAddSIMD(
                                                                                            FDPart3tmp33,
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp144,
                                                                                                FDPart3tmp176,
                                                                                                MulSIMD(
                                                                                                    FDPart3_Integer_2,
                                                                                                    MulSIMD(FDPart3tmp144, FDPart3tmp159))),
                                                                                            FusedMulAddSIMD(
                                                                                                FDPart3tmp33,
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp145, FDPart3tmp173,
                                                                                                    MulSIMD(FDPart3tmp150, FDPart3tmp364)),
                                                                                                FusedMulAddSIMD(
                                                                                                    FDPart3tmp285, lambdaU_dD22,
                                                                                                    FusedMulAddSIMD(
                                                                                                        FDPart3tmp33,
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp142, FDPart3tmp269,
                                                                                                            MulSIMD(
                                                                                                                FDPart3_Integer_2,
                                                                                                                MulSIMD(
                                                                                                                    FDPart3tmp142, FDPart3tmp267))),
                                                                                                        FusedMulAddSIMD(
                                                                                                            FDPart3tmp27,
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp269, FDPart3tmp92,
                                                                                                                MulSIMD(
                                                                                                                    FDPart3_Integer_2,
                                                                                                                    MulSIMD(
                                                                                                                        FDPart3tmp267,
                                                                                                                        FDPart3tmp92))),
                                                                                                            FusedMulAddSIMD(
                                                                                                                FDPart3tmp282, lambdaU_dD02,
                                                                                                                FusedMulAddSIMD(
                                                                                                                    FDPart3tmp27,
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp173, FDPart3tmp98,
                                                                                                                        MulSIMD(
                                                                                                                            FDPart3tmp196,
                                                                                                                            FDPart3tmp364)),
                                                                                                                    FusedMulAddSIMD(
                                                                                                                        FDPart3tmp27,
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp176,
                                                                                                                            FDPart3tmp88,
                                                                                                                            MulSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                MulSIMD(
                                                                                                                                    FDPart3tmp159,
                                                                                                                                    FDPart3tmp88))),
                                                                                                                        FusedMulAddSIMD(
                                                                                                                            FDPart3tmp215,
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3_Integer_2,
                                                                                                                                MulSIMD(FDPart3tmp25,
                                                                                                                                        FDPart3tmp88),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp31,
                                                                                                                                        FDPart3tmp92),
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3_Integer_2,
                                                                                                                                        MulSIMD(
                                                                                                                                            FDPart3tmp13,
                                                                                                                                            FDPart3tmp95)))),
                                                                                                                            FusedMulAddSIMD(
                                                                                                                                FDPart3tmp27,
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp136,
                                                                                                                                    FDPart3tmp172,
                                                                                                                                    MulSIMD(
                                                                                                                                        FDPart3tmp150,
                                                                                                                                        FDPart3tmp363)),
                                                                                                                                FusedMulAddSIMD(
                                                                                                                                    FDPart3tmp212,
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp31,
                                                                                                                                        FDPart3tmp363,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3_Integer_2,
                                                                                                                                            MulSIMD(
                                                                                                                                                FDPart3tmp13,
                                                                                                                                                FDPart3tmp175),
                                                                                                                                            MulSIMD(
                                                                                                                                                FDPart3tmp25,
                                                                                                                                                FDPart3tmp364))),
                                                                                                                                    FusedMulAddSIMD(
                                                                                                                                        FDPart3tmp214,
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3_Integer_2,
                                                                                                                                            MulSIMD(
                                                                                                                                                FDPart3tmp142,
                                                                                                                                                FDPart3tmp31),
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3_Integer_2,
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3tmp144,
                                                                                                                                                    FDPart3tmp25),
                                                                                                                                                MulSIMD(
                                                                                                                                                    FDPart3_Integer_2,
                                                                                                                                                    MulSIMD(
                                                                                                                                                        FDPart3tmp13,
                                                                                                                                                        FDPart3tmp143)))),
                                                                                                                                        FusedMulAddSIMD(
                                                                                                                                            FDPart3tmp128,
                                                                                                                                            FDPart3tmp268,
                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                FDPart3tmp135,
                                                                                                                                                FDPart3tmp160,
                                                                                                                                                FusedMulAddSIMD(
                                                                                                                                                    FDPart3tmp22,
                                                                                                                                                    MulSIMD(
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3_NegativeOne_,
                                                                                                                                                            FDPart3_Rational_1_2),
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp60,
                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                FDPart3tmp8,
                                                                                                                                                                hDD_dDD2200,
                                                                                                                                                                MulSIMD(
                                                                                                                                                                    FDPart3tmp166,
                                                                                                                                                                    MulSIMD(
                                                                                                                                                                        FDPart3tmp5,
                                                                                                                                                                        f0_of_xx0__DD00))))),
                                                                                                                                                    FusedMulSubSIMD(
                                                                                                                                                        FDPart3tmp126,
                                                                                                                                                        FDPart3tmp341,
                                                                                                                                                        MulSIMD(
                                                                                                                                                            FDPart3tmp86,
                                                                                                                                                            FusedMulAddSIMD(
                                                                                                                                                                FDPart3tmp138,
                                                                                                                                                                FDPart3tmp58,
                                                                                                                                                                MulSIMD(
                                                                                                                                                                    FDPart3tmp8,
                                                                                                                                                                    hDD_dDD2201)))))))))))))))))))))))))))))))))))))))));

WriteSIMD(&auxevol_gfs[IDX4(RBARDD00GF, i0, i1, i2)], __RHS_exp_0);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD01GF, i0, i1, i2)], __RHS_exp_1);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD02GF, i0, i1, i2)], __RHS_exp_2);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD11GF, i0, i1, i2)], __RHS_exp_3);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD12GF, i0, i1, i2)], __RHS_exp_4);
WriteSIMD(&auxevol_gfs[IDX4(RBARDD22GF, i0, i1, i2)], __RHS_exp_5);

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 += simd_width)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
} // END FUNCTION Ricci_eval__rfm__SinhCylindrical
