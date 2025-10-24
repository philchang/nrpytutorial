#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Perform spacetime evolution in BlackHoles@Home
 */
void bhah_evolve(BHaH_struct *bhah_struct) {

  commondata_struct *commondata = bhah_struct->commondata;
  griddata_struct *griddata = bhah_struct->griddata;
  while (commondata->time < commondata->t_final) {
    bhah_diagnostics(bhah_struct);
    MoL_step_forward_in_time(commondata, griddata);
  }
}
