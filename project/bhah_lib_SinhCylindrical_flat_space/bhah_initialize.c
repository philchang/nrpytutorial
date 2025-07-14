#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * BlackHoles@Home library setup
 */
BHaH_struct *bhah_initialize() {

  // Step 1.a: Allocate memory for the BH@H struct
  BHaH_struct *bhah_struct = (BHaH_struct *)malloc(sizeof(BHaH_struct));
  commondata_struct *commondata = (commondata_struct *)malloc(sizeof(commondata_struct));

  // Step 1.b: Set each commondata CodeParameter to default.
  commondata_struct_set_to_default(commondata);

  // Step 1.c: Allocate NUMGRIDS griddata arrays, each containing
  //           data specific to an individual grid.
  const int n_grids = commondata->NUMGRIDS;
  griddata_struct *griddata = (griddata_struct *)malloc(sizeof(griddata_struct) * n_grids);

  // Step 1.d: Set each CodeParameter in griddata.params to default.
  params_struct_set_to_default(commondata, griddata);

  // Step 1.e: Set up numerical grids: xx[3], masks, Nxx, dxx, invdxx, bcstruct,
  //           rfm_precompute, timestep, etc. With calling_for_first_time = true,
  //           commondata time is set to zero (as is the iteration, etc).
  const bool calling_for_first_time = true;
  numerical_grids_and_timestep(commondata, griddata, calling_for_first_time);

  // Step 1.f: Allocate memory for all gridfunctions
  for (int grid = 0; grid < n_grids; grid++) {
    MoL_malloc_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
    MoL_malloc_non_y_n_gfs(commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
  }

  // Step 1.g: Set up initial data
  initial_data(commondata, griddata);

  // Step 1.h: Set pointers to griddata and commondata
  bhah_struct->griddata = griddata;
  bhah_struct->commondata = commondata;

  return bhah_struct;
}
