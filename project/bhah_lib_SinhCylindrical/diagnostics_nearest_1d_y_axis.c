#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Output diagnostic quantities at gridpoints closest to y axis.
 */
void diagnostics_nearest_1d_y_axis(commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                   MoL_gridfunctions_struct *restrict gridfuncs) {
  switch (params->CoordSystem_hash) {
  case SINHCYLINDRICAL:
    diagnostics_nearest_1d_y_axis__rfm__SinhCylindrical(commondata, params, xx, gridfuncs);
    break;
  default:
    fprintf(stderr, "ERROR in diagnostics_nearest_1d_y_axis(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}
