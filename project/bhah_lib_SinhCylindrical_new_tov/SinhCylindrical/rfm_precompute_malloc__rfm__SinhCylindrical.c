#include "../BHaH_defines.h"
/**
 * rfm_precompute_malloc: reference metric precomputed lookup arrays: malloc
 */
void rfm_precompute_malloc__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                 rfm_struct *restrict rfmstruct) {
#include "../set_CodeParameters.h"
  rfmstruct->f0_of_xx0 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__D0 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DD00 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DDD000 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  rfmstruct->f3_of_xx2 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS2);
  rfmstruct->f3_of_xx2__D2 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS2);
  rfmstruct->f3_of_xx2__DD22 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS2);
}
