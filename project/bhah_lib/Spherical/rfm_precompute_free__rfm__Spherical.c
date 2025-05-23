#include "../BHaH_defines.h"
/**
 * rfm_precompute_free: reference metric precomputed lookup arrays: free
 */
void rfm_precompute_free__rfm__Spherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                         rfm_struct *restrict rfmstruct) {
#include "../set_CodeParameters.h"
  free(rfmstruct->f0_of_xx0);
  free(rfmstruct->f1_of_xx1);
  free(rfmstruct->f1_of_xx1__D1);
  free(rfmstruct->f1_of_xx1__DD11);
}
