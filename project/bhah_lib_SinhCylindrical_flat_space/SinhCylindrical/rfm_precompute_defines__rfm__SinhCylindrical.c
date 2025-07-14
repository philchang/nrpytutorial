#include "../BHaH_defines.h"
/**
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                  rfm_struct *restrict rfmstruct, REAL *restrict xx[3]) {
#include "../set_CodeParameters.h"
  for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f0_of_xx0[i0] = AMPLRHO * (exp(xx0 / SINHWRHO) - exp(-xx0 / SINHWRHO)) / (exp(1.0 / SINHWRHO) - exp(-1 / SINHWRHO));
  }

  for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f0_of_xx0__D0[i0] =
        AMPLRHO * (exp(xx0 / SINHWRHO) / SINHWRHO + exp(-xx0 / SINHWRHO) / SINHWRHO) / (exp(1.0 / SINHWRHO) - exp(-1 / SINHWRHO));
  }

  for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f0_of_xx0__DD00[i0] =
        AMPLRHO * (exp(xx0 / SINHWRHO) / pow(SINHWRHO, 2) - exp(-xx0 / SINHWRHO) / pow(SINHWRHO, 2)) / (exp(1.0 / SINHWRHO) - exp(-1 / SINHWRHO));
  }

  for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
    const REAL xx0 = xx[0][i0];
    rfmstruct->f0_of_xx0__DDD000[i0] =
        AMPLRHO * (exp(xx0 / SINHWRHO) / pow(SINHWRHO, 3) + exp(-xx0 / SINHWRHO) / pow(SINHWRHO, 3)) / (exp(1.0 / SINHWRHO) - exp(-1 / SINHWRHO));
  }

  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = xx[2][i2];
    rfmstruct->f3_of_xx2[i2] = AMPLZ * (exp(xx2 / SINHWZ) / SINHWZ + exp(-xx2 / SINHWZ) / SINHWZ) / (exp(1.0 / SINHWZ) - exp(-1 / SINHWZ));
  }

  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = xx[2][i2];
    rfmstruct->f3_of_xx2__D2[i2] =
        AMPLZ * (exp(xx2 / SINHWZ) / pow(SINHWZ, 2) - exp(-xx2 / SINHWZ) / pow(SINHWZ, 2)) / (exp(1.0 / SINHWZ) - exp(-1 / SINHWZ));
  }

  for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
    const REAL xx2 = xx[2][i2];
    rfmstruct->f3_of_xx2__DD22[i2] =
        AMPLZ * (exp(xx2 / SINHWZ) / pow(SINHWZ, 3) + exp(-xx2 / SINHWZ) / pow(SINHWZ, 3)) / (exp(1.0 / SINHWZ) - exp(-1 / SINHWZ));
  }
}
