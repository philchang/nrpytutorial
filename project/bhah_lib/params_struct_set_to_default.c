#include "BHaH_defines.h"
/**
 * Set params_struct to default values specified within NRPy+.
 */
void params_struct_set_to_default(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  // Loop over params structs:
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    params_struct *restrict params = &griddata[grid].params;
    // Set params_struct variables to default
    params->Cart_originx = 0.0;                         // nrpy.grid::Cart_originx
    params->Cart_originy = 0.0;                         // nrpy.grid::Cart_originy
    params->Cart_originz = 0.0;                         // nrpy.grid::Cart_originz
    params->Nxx0 = 96;                                  // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx0
    params->Nxx1 = 16;                                  // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx1
    params->Nxx2 = 2;                                   // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx2
    params->PI = 3.14159265358979323846264338327950288; // nrpy.reference_metric::PI
    params->RMAX = 10.0;                                // nrpy.reference_metric_Spherical::RMAX
    params->grid_hole_radius = 2.0;                     // nrpy.reference_metric::grid_hole_radius
    params->grid_physical_size = 10.0;                  // nrpy.reference_metric::grid_physical_size
    params->grid_rotates = false;                       // nrpy.grid::grid_rotates
    params->xxmax0 = 10.0;                              // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax0
    params->xxmax1 = 10.0;                              // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax1
    params->xxmax2 = 10.0;                              // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax2
    params->xxmin0 = -10.0;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin0
    params->xxmin1 = -10.0;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin1
    params->xxmin2 = -10.0;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin2
    snprintf(params->CoordSystemName, 50, "Spherical"); // nrpy.reference_metric::CoordSystemName
  }
}
