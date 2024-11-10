#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Compute minimum timestep dt = CFL_FACTOR * ds_min on a Spherical numerical grid.
 */
void cfl_limited_timestep__rfm__Spherical(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]) {
#include "../set_CodeParameters.h"

  REAL ds_min = 1e38;
#pragma omp parallel for reduction(min : ds_min)
  LOOP_NOOMP(i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    const REAL xx2 = xx[2][i2];
    REAL dsmin0, dsmin1, dsmin2;
    /*
     *  Original SymPy expressions:
     *  "[dsmin0 = Abs(dxx0)]"
     *  "[dsmin1 = Abs(dxx1*xx0)]"
     *  "[dsmin2 = Abs(dxx2*xx0*sin(xx1))]"
     */
    dsmin0 = fabs(dxx0);
    dsmin1 = fabs(dxx1 * xx0);
    dsmin2 = fabs(dxx2 * xx0 * sin(xx1));
    ds_min = MIN(ds_min, MIN(dsmin0, MIN(dsmin1, dsmin2)));
  }
  commondata->dt = MIN(commondata->dt, ds_min * commondata->CFL_FACTOR);
}
