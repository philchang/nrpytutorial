#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Compute minimum timestep dt = CFL_FACTOR * ds_min on a SinhCylindrical numerical grid.
 */
void cfl_limited_timestep__rfm__SinhCylindrical(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]) {
#include "../set_CodeParameters.h"

  REAL ds_min = 1e38;
#pragma omp parallel for reduction(min : ds_min)
  LOOP_NOOMP(i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
    MAYBE_UNUSED const REAL xx0 = xx[0][i0];
    MAYBE_UNUSED const REAL xx1 = xx[1][i1];
    MAYBE_UNUSED const REAL xx2 = xx[2][i2];
    REAL dsmin0, dsmin1, dsmin2;
    /*
     *  Original SymPy expressions:
     *  "[dsmin0 = Abs(AMPLRHO*dxx0*(exp(xx0/SINHWRHO)/SINHWRHO + exp(-xx0/SINHWRHO)/SINHWRHO)/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)))]"
     *  "[dsmin1 = Abs(AMPLRHO*dxx1*(exp(xx0/SINHWRHO) - exp(-xx0/SINHWRHO))/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)))]"
     *  "[dsmin2 = Abs(AMPLZ*dxx2*(exp(xx2/SINHWZ)/SINHWZ + exp(-xx2/SINHWZ)/SINHWZ)/(exp(1/SINHWZ) - exp(-1/SINHWZ)))]"
     */
    const REAL tmp0 = (1.0 / (SINHWRHO));
    const REAL tmp5 = (1.0 / (SINHWZ));
    const REAL tmp4 = AMPLRHO / (exp(tmp0) - exp(-tmp0));
    const REAL tmp2 = exp(tmp0 * xx0);
    const REAL tmp3 = exp(-tmp0 * xx0);
    dsmin0 = fabs(dxx0 * tmp4 * (tmp0 * tmp2 + tmp0 * tmp3));
    dsmin1 = fabs(dxx1 * tmp4 * (tmp2 - tmp3));
    dsmin2 = fabs(AMPLZ * dxx2 * (tmp5 * exp(tmp5 * xx2) + tmp5 * exp(-tmp5 * xx2)) / (exp(tmp5) - exp(-tmp5)));
    ds_min = MIN(ds_min, MIN(dsmin0, MIN(dsmin1, dsmin2)));
  }
  commondata->dt = MIN(commondata->dt, ds_min * commondata->CFL_FACTOR);
}
