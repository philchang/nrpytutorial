#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 *
 * compute time component of four velocity, via
 *
 * // Derivation of first equation:
 * // \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2
 * //   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
 * //   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
 * //   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
 * //   = 1 - 1/(u^0 \alpha)^2 <= 1
 *
 */
void compute_u4Ut__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                        const REAL max_Lorentz_factor, const int i0, const int i1, const int i2, REAL *restrict in_gfs,
                                        REAL *restrict rescaledvU0, REAL *restrict rescaledvU1, REAL *restrict rescaledvU2, REAL *restrict u4Ut) {
#include "../set_CodeParameters.h"

  // 1 - W_max^{-2}
  const REAL inv_sq_max_Lorentz_factor = 1.0 / SQR(max_Lorentz_factor);
  const REAL one_minus_one_over_W_max_squared = 1.0 - inv_sq_max_Lorentz_factor;

  const REAL rescaledvU_old0 = *rescaledvU0;
  const REAL rescaledvU_old1 = *rescaledvU1;
  const REAL rescaledvU_old2 = *rescaledvU2;
  /*
   * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
   * Read gridfunction(s) from main memory and compute FD stencils as needed.
   */
  const REAL alpha = in_gfs[IDX4(ALPHAGF, i0, i1, i2)];
  const REAL cf = in_gfs[IDX4(CFGF, i0, i1, i2)];
  const REAL hDD00 = in_gfs[IDX4(HDD00GF, i0, i1, i2)];
  const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
  const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
  const REAL hDD11 = in_gfs[IDX4(HDD11GF, i0, i1, i2)];
  const REAL hDD12 = in_gfs[IDX4(HDD12GF, i0, i1, i2)];
  const REAL hDD22 = in_gfs[IDX4(HDD22GF, i0, i1, i2)];
  const REAL vetU0 = in_gfs[IDX4(VETU0GF, i0, i1, i2)];
  const REAL vetU1 = in_gfs[IDX4(VETU1GF, i0, i1, i2)];
  const REAL vetU2 = in_gfs[IDX4(VETU2GF, i0, i1, i2)];

  /*
   * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
   * Evaluate SymPy expressions and write to main memory.
   */
  const REAL FDPart3tmp0 = (1.0 / 2.0) * one_minus_one_over_W_max_squared;
  const REAL FDPart3tmp7 = 2 * rescaledvU_old0;
  const REAL FDPart3tmp9 = 2 * rescaledvU_old1;
  const REAL FDPart3tmp11 = 2 * rescaledvU_old2;
  const REAL FDPart3tmp19 =
      (FDPart3tmp11 * hDD02 * vetU0 + FDPart3tmp11 * hDD12 * vetU1 + FDPart3tmp11 * hDD22 * vetU2 + FDPart3tmp11 * vetU2 +
       FDPart3tmp7 * hDD00 * vetU0 + FDPart3tmp7 * hDD01 * rescaledvU_old1 + FDPart3tmp7 * hDD01 * vetU1 + FDPart3tmp7 * hDD02 * rescaledvU_old2 +
       FDPart3tmp7 * hDD02 * vetU2 + FDPart3tmp7 * vetU0 + FDPart3tmp9 * hDD01 * vetU0 + FDPart3tmp9 * hDD11 * vetU1 +
       FDPart3tmp9 * hDD12 * rescaledvU_old2 + FDPart3tmp9 * hDD12 * vetU2 + FDPart3tmp9 * vetU1 + hDD00 * ((rescaledvU_old0) * (rescaledvU_old0)) +
       hDD00 * ((vetU0) * (vetU0)) + 2 * hDD01 * vetU0 * vetU1 + 2 * hDD02 * vetU0 * vetU2 + hDD11 * ((rescaledvU_old1) * (rescaledvU_old1)) +
       hDD11 * ((vetU1) * (vetU1)) + 2 * hDD12 * vetU1 * vetU2 + hDD22 * ((rescaledvU_old2) * (rescaledvU_old2)) + hDD22 * ((vetU2) * (vetU2)) +
       ((rescaledvU_old0) * (rescaledvU_old0)) + ((rescaledvU_old1) * (rescaledvU_old1)) + ((rescaledvU_old2) * (rescaledvU_old2)) +
       ((vetU0) * (vetU0)) + ((vetU1) * (vetU1)) + ((vetU2) * (vetU2))) /
      (((alpha) * (alpha)) * ((cf) * (cf)));
  const REAL FDPart3tmp20 = (1.0 / 2.0) * FDPart3tmp19;
  const REAL FDPart3tmp21 = -FDPart3tmp19 + one_minus_one_over_W_max_squared;
  const REAL FDPart3tmp23 = sqrt(one_minus_one_over_W_max_squared / (FDPart3tmp19 + TINYDOUBLE));
  const REAL FDPart3tmp24 = FDPart3tmp19 + TINYDOUBLE - one_minus_one_over_W_max_squared;
  const REAL FDPart3tmp22 = (-FDPart3tmp0 + FDPart3tmp20 - 1.0 / 2.0 * fabs(-FDPart3tmp21)) / (-FDPart3tmp21 - TINYDOUBLE);
  const REAL FDPart3tmp25 = (-FDPart3tmp0 + FDPart3tmp20 + (1.0 / 2.0) * TINYDOUBLE + (1.0 / 2.0) * fabs(FDPart3tmp24)) / FDPart3tmp24;
  *u4Ut = 1.0 / (alpha * sqrt(-FDPart3tmp0 - FDPart3tmp20 + (1.0 / 2.0) * fabs(FDPart3tmp21) + 1.0));
  *rescaledvU0 = FDPart3tmp22 * rescaledvU_old0 + FDPart3tmp25 * (FDPart3tmp23 * (rescaledvU_old0 + vetU0) - vetU0);
  *rescaledvU1 = FDPart3tmp22 * rescaledvU_old1 + FDPart3tmp25 * (FDPart3tmp23 * (rescaledvU_old1 + vetU1) - vetU1);
  *rescaledvU2 = FDPart3tmp22 * rescaledvU_old2 + FDPart3tmp25 * (FDPart3tmp23 * (rescaledvU_old2 + vetU2) - vetU2);
} // END FUNCTION compute_u4Ut__rfm__SinhCylindrical
