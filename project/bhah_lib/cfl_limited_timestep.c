#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Compute minimum timestep dt = CFL_FACTOR * ds_min on a Spherical numerical grid.
 */
void cfl_limited_timestep(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]) {
  switch (params->CoordSystem_hash) {
  case SPHERICAL:
    cfl_limited_timestep__rfm__Spherical(commondata, params, xx);
    break;
  default:
    fprintf(stderr, "ERROR in cfl_limited_timestep(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}
