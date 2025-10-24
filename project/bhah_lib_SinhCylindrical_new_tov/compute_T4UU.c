#include "BHaH_defines.h"
/**
 * Compute T4UU in curvilinear coordinates in a point-wise fashion.
 */
void compute_T4UU(const commondata_struct *restrict commondata, const params_struct *restrict params,
                  const int i0, const int i1, const int i2, REAL *restrict xx[3], const REAL rhob,
                  const REAL P, const REAL h, const REAL u4Ut, const REAL rescaledvU0, const REAL rescaledvU1,
                  const REAL rescaledvU2, const REAL *restrict in_gfs, REAL *restrict auxevol_gfs) {
#include "set_CodeParameters.h"
  // Set up coordinates at (i0, i1, i2)
  MAYBE_UNUSED REAL xx0 = xx[0][i0];
  MAYBE_UNUSED REAL xx1 = xx[1][i1];
  MAYBE_UNUSED REAL xx2 = xx[2][i2];
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
  const REAL FDPart3tmp0 = ((alpha) * (alpha));
  const REAL FDPart3tmp5 = (1.0 / (SINHWRHO));
  const REAL FDPart3tmp17 = (1.0 / (SINHWZ));
  const REAL FDPart3tmp24 = hDD22 + 1;
  const REAL FDPart3tmp26 = hDD11 + 1;
  const REAL FDPart3tmp28 = hDD00 + 1;
  const REAL FDPart3tmp1 = (1.0 / (FDPart3tmp0));
  const REAL FDPart3tmp3 = FDPart3tmp0 * h * rhob * ((u4Ut) * (u4Ut));
  const REAL FDPart3tmp6 = 2 * FDPart3tmp5;
  const REAL FDPart3tmp18 = 2 * FDPart3tmp17;
  const REAL FDPart3tmp29 = FDPart3tmp24 * FDPart3tmp26 * FDPart3tmp28 - FDPart3tmp24 * ((hDD01) * (hDD01)) - FDPart3tmp26 * ((hDD02) * (hDD02)) -
                            FDPart3tmp28 * ((hDD12) * (hDD12)) + 2 * hDD01 * hDD02 * hDD12;
  const REAL FDPart3tmp31 = FDPart3tmp0 * ((cf) * (cf));
  const REAL FDPart3tmp7 = exp(FDPart3tmp6) - 1;
  const REAL FDPart3tmp9 = exp(FDPart3tmp6 * xx0);
  const REAL FDPart3tmp19 = exp(FDPart3tmp18 * xx2);
  const REAL FDPart3tmp21 = exp(FDPart3tmp18) - 1;
  const REAL FDPart3tmp30 = FDPart3tmp29 * FDPart3tmp3;
  const REAL FDPart3tmp32 = (1.0 / (FDPart3tmp29));
  const REAL FDPart3tmp8 = FDPart3tmp7 * exp(-FDPart3tmp5) * exp(FDPart3tmp5 * xx0) / AMPLRHO;
  const REAL FDPart3tmp34 = ((FDPart3tmp7) * (FDPart3tmp7)) * FDPart3tmp9 * exp(-FDPart3tmp6) / ((AMPLRHO) * (AMPLRHO));
  const REAL FDPart3tmp36 = FDPart3tmp29 * FDPart3tmp3 * rescaledvU0;
  const REAL FDPart3tmp11 = SINHWRHO / (FDPart3tmp9 + 1);
  const REAL FDPart3tmp15 = FDPart3tmp1 / (FDPart3tmp9 - 1);
  const REAL FDPart3tmp22 = FDPart3tmp21 * SINHWZ * exp(-FDPart3tmp17) * exp(FDPart3tmp17 * xx2) / (AMPLZ * (FDPart3tmp19 + 1));
  const REAL FDPart3tmp35 = FDPart3tmp1 * FDPart3tmp32 * FDPart3tmp34;
  const REAL FDPart3tmp12 = FDPart3tmp1 * FDPart3tmp11 * FDPart3tmp8;
  auxevol_gfs[IDX4(T4UU00GF, i0, i1, i2)] = FDPart3tmp1 * (FDPart3tmp0 * h * rhob * ((u4Ut) * (u4Ut)) - P);
  auxevol_gfs[IDX4(T4UU01GF, i0, i1, i2)] = FDPart3tmp12 * (FDPart3tmp3 * rescaledvU0 + P * vetU0);
  auxevol_gfs[IDX4(T4UU02GF, i0, i1, i2)] = FDPart3tmp15 * FDPart3tmp8 * (FDPart3tmp3 * rescaledvU1 + P * vetU1);
  auxevol_gfs[IDX4(T4UU03GF, i0, i1, i2)] = FDPart3tmp1 * FDPart3tmp22 * (FDPart3tmp3 * rescaledvU2 + P * vetU2);
  auxevol_gfs[IDX4(T4UU11GF, i0, i1, i2)] =
      FDPart3tmp35 * ((SINHWRHO) * (SINHWRHO)) *
      (FDPart3tmp30 * ((rescaledvU0) * (rescaledvU0)) +
       P * (-FDPart3tmp29 * ((vetU0) * (vetU0)) + FDPart3tmp31 * (FDPart3tmp24 * FDPart3tmp26 - ((hDD12) * (hDD12))))) /
      ((FDPart3tmp9 + 1) * (FDPart3tmp9 + 1));
  auxevol_gfs[IDX4(T4UU12GF, i0, i1, i2)] =
      FDPart3tmp11 * FDPart3tmp15 * FDPart3tmp32 * FDPart3tmp34 *
      (FDPart3tmp36 * rescaledvU1 + P * (-FDPart3tmp29 * vetU0 * vetU1 + FDPart3tmp31 * (-FDPart3tmp24 * hDD01 + hDD02 * hDD12)));
  auxevol_gfs[IDX4(T4UU13GF, i0, i1, i2)] =
      FDPart3tmp12 * FDPart3tmp22 * FDPart3tmp32 *
      (FDPart3tmp36 * rescaledvU2 + P * (-FDPart3tmp29 * vetU0 * vetU2 + FDPart3tmp31 * (-FDPart3tmp26 * hDD02 + hDD01 * hDD12)));
  auxevol_gfs[IDX4(T4UU22GF, i0, i1, i2)] =
      FDPart3tmp35 *
      (FDPart3tmp30 * ((rescaledvU1) * (rescaledvU1)) +
       P * (-FDPart3tmp29 * ((vetU1) * (vetU1)) + FDPart3tmp31 * (FDPart3tmp24 * FDPart3tmp28 - ((hDD02) * (hDD02))))) /
      ((FDPart3tmp9 - 1) * (FDPart3tmp9 - 1));
  auxevol_gfs[IDX4(T4UU23GF, i0, i1, i2)] = FDPart3tmp15 * FDPart3tmp22 * FDPart3tmp32 * FDPart3tmp8 *
                                            (FDPart3tmp29 * FDPart3tmp3 * rescaledvU1 * rescaledvU2 +
                                             P * (-FDPart3tmp29 * vetU1 * vetU2 + FDPart3tmp31 * (-FDPart3tmp28 * hDD12 + hDD01 * hDD02)));
  auxevol_gfs[IDX4(T4UU33GF, i0, i1, i2)] =
      FDPart3tmp1 * FDPart3tmp19 * ((FDPart3tmp21) * (FDPart3tmp21)) * FDPart3tmp32 * ((SINHWZ) * (SINHWZ)) *
      (FDPart3tmp30 * ((rescaledvU2) * (rescaledvU2)) +
       P * (-FDPart3tmp29 * ((vetU2) * (vetU2)) + FDPart3tmp31 * (FDPart3tmp26 * FDPart3tmp28 - ((hDD01) * (hDD01))))) *
      exp(-FDPart3tmp18) / (((AMPLZ) * (AMPLZ)) * ((FDPart3tmp19 + 1) * (FDPart3tmp19 + 1)));
} // END FUNCTION compute_T4UU
