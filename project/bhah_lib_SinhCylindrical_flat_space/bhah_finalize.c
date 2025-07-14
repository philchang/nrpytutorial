#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Finalize the BlackHoles@Home library
 */
void bhah_finalize(BHaH_struct *bhah_struct) {

  commondata_struct *commondata = bhah_struct->commondata;
  griddata_struct *griddata = bhah_struct->griddata;
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    MoL_free_memory_y_n_gfs(&griddata[grid].gridfuncs);
    MoL_free_memory_non_y_n_gfs(&griddata[grid].gridfuncs);
    rfm_precompute_free(commondata, &griddata[grid].params, &griddata[grid].rfmstruct);
  }
  free(bhah_struct->commondata);
  free(bhah_struct->griddata);
  free(bhah_struct);
}
