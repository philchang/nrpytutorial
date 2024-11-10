#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Set up numerical grids and timestep.
 */
void numerical_grids_and_timestep(commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time) {

  // Thiago says: I am commenting out the following line as it would set Nxx0, Nxx1, and Nxx2, RMAX and, grid_physical_size to the default values
  // // Step 1.a: Set each CodeParameter in griddata.params to default, for MAXNUMGRIDS grids.
  // params_struct_set_to_default(commondata, griddata);

  // Thiago says: BHaH tries really hard to set Nxx0, Nxx1 and Nxx2 to default values. 
  //              The numerical_grid_params_Nxx_dxx_xx() function overwrites 
  //              (griddata->params).Nxx0, (griddata->params).Nxx1, and (griddata->params).Nxx2
  //              with the values of int Nx[3]. So here we set them to the current values, which will overwrite the default values.
  // // Independent grids
  // int Nx[3] = {-1, -1, -1};
  int Nx[3] = {(griddata->params).Nxx0, (griddata->params).Nxx1, (griddata->params).Nxx2};

  // Step 1.b: Set commondata->NUMGRIDS to number of CoordSystems we have
  commondata->NUMGRIDS = 1;

  {
    // Step 1.c: For each grid, set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
    const bool set_xxmin_xxmax_to_defaults = true;
    int grid = 0;
    griddata[grid].params.CoordSystem_hash = SPHERICAL;
    // griddata[grid].params.grid_physical_size = 1.6; // Thiago says: commenting this out to keep the value of grid_physical_size to rmax (set in MaNGa)
    numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx, Nx, set_xxmin_xxmax_to_defaults);
    grid++;
  }

  // Step 1.d: Allocate memory for and define reference-metric precomputation lookup tables
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    rfm_precompute_malloc(commondata, &griddata[grid].params, &griddata[grid].rfmstruct);
    rfm_precompute_defines(commondata, &griddata[grid].params, &griddata[grid].rfmstruct, griddata[grid].xx);
  }

  // Step 1.e: Set up curvilinear boundary condition struct (bcstruct)
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    bcstruct_set_up(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
  }

  // Step 1.f: Set timestep based on minimum spacing between neighboring gridpoints.
  commondata->dt = 1e30;
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx);
  }
  // Step 1.g: Initialize timestepping parameters to zero if this is the first time this function is called.
  if (calling_for_first_time) {
    commondata->nn = 0;
    commondata->nn_0 = 0;
    commondata->t_0 = 0.0;
    commondata->time = 0.0;
  }
}
