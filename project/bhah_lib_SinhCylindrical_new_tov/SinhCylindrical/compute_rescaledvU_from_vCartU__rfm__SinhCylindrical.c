#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Compute rescaled 3-vector in curvilinar coordinates from 3-vector in Cartesian coordinates.
 */
void compute_rescaledvU_from_vCartU__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                          const REAL vCartU[3], const REAL xx0, const REAL xx1, const REAL xx2,
                                                          REAL *restrict rescaledvU0, REAL *restrict rescaledvU1, REAL *restrict rescaledvU2) {
#include "../set_CodeParameters.h"
  const REAL tmp0 = cos(xx1);
  const REAL tmp1 = (1.0 / (SINHWRHO));
  const REAL tmp10 = sin(xx1);
  const REAL tmp14 = (1.0 / (SINHWZ));
  const REAL tmp6 = exp(tmp1) - exp(-tmp1);
  const REAL tmp15 = (1.0 / (exp(tmp14) - exp(-tmp14)));
  const REAL tmp3 = exp(tmp1 * xx0);
  const REAL tmp4 = exp(-tmp1 * xx0);
  const REAL tmp7 = (1.0 / (tmp6));
  const REAL tmp17 = tmp14 * exp(tmp14 * xx2) + tmp14 * exp(-tmp14 * xx2);
  const REAL tmp5 = tmp3 - tmp4;
  const REAL tmp11 = tmp1 * tmp3 + tmp1 * tmp4;
  const REAL tmp18 = AMPLZ * tmp15 * tmp17;
  const REAL tmp9 = AMPLRHO * tmp5 * tmp7;
  const REAL tmp12 = ((AMPLRHO) * (AMPLRHO)) * tmp11 * tmp5 / ((tmp6) * (tmp6));
  const REAL tmp23 = AMPLRHO * tmp11 * tmp7;
  const REAL tmp13 = ((tmp10) * (tmp10)) * tmp12;
  const REAL tmp19 = ((tmp0) * (tmp0)) * tmp12;
  const REAL tmp20 = (1.0 / (tmp13 * tmp18 + tmp18 * tmp19));
  const REAL tmp21 = tmp18 * tmp20;
  *rescaledvU0 = tmp23 * (tmp0 * tmp21 * tmp9 * vCartU[0] + tmp10 * tmp21 * tmp9 * vCartU[1]);
  *rescaledvU1 = tmp9 * (AMPLRHO * AMPLZ * tmp0 * tmp11 * tmp15 * tmp17 * tmp20 * tmp7 * vCartU[1] - tmp10 * tmp21 * tmp23 * vCartU[0]);
  *rescaledvU2 = tmp21 * vCartU[2] * (tmp13 + tmp19);
} // END FUNCTION compute_rescaledvU_from_vCartU__rfm__SinhCylindrical
