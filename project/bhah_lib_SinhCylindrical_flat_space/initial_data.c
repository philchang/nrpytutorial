#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Set initial data.
 */
void initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  // Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
  if (read_checkpoint(commondata, griddata))
    return;
  ID_persist_struct ID_persist;

  TOVola_solve(commondata, &ID_persist);

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    // Unpack griddata struct:
    params_struct *restrict params = &griddata[grid].params;
    initial_data_reader__convert_ADM_Spherical_to_BSSN(commondata, params, griddata[grid].xx, &griddata[grid].bcstruct, &griddata[grid].gridfuncs,
                                                       &ID_persist, TOVola_interp);
    apply_bcs_outerextrap_and_inner(commondata, params, &griddata[grid].bcstruct, griddata[grid].gridfuncs.y_n_gfs);
  }

  {
    free(ID_persist.r_Schw_arr);
    free(ID_persist.rho_energy_arr);
    free(ID_persist.rho_baryon_arr);
    free(ID_persist.P_arr);
    free(ID_persist.M_arr);
    free(ID_persist.expnu_arr);
    free(ID_persist.exp4phi_arr);
    free(ID_persist.r_iso_arr);
  }
}

/**
 * Set radial initial data quantities for MaNGa
 */
void manga_radial_initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                               const int num_radial_pts, REAL *restrict r_axis, REAL *restrict rho_baryon, REAL *restrict pressure) {

  // Set params struct for a single grid
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;

  // Declare ID struct and populate it with TOV solution
  ID_persist_struct ID_persist;
  TOVola_solve(commondata, &ID_persist);

  // Declare quantities for sampling local radial axis using SinhSpherical coordinates (this has nothing to do with the coordinates used by BHaH for evolution)
  const REAL domain_size = params->AMPLRHO;
  const REAL local_SINHW = 0.2;
  const REAL dx = 1.0 / ((REAL)num_radial_pts);

  // Populate arrays by interpolating rho_baryon and pressure onto r_axis

  const int R_idx = ID_persist.numpoints_arr - 1;
  const REAL r_iso_max_inside_star = ID_persist.r_iso_arr[R_idx];

  for (int i = 0; i < num_radial_pts; i++) {
    const REAL x = dx / 2.0 + dx * i;
    const REAL r_iso = domain_size * sinh(x / local_SINHW) / sinh(1.0 / local_SINHW);
    r_axis[i] = r_iso;
    if (r_iso <= (r_iso_max_inside_star / 5.0)) {
      // rho_baryon[i] = commondata->initial_central_density;
      rho_baryon[i] = 1.0e-3;
      pressure[i] = 0.0;
    }
    else {
      rho_baryon[i] = 1.0e-10;
      pressure[i] = 0.0;
    }

    // TOVola_radial_only_interp(commondata, params, r_axis[i], &ID_persist, &rho_baryon[i], &pressure[i]);
  }

  /** START DEBUGGING **/
  // printf("Inside manga_radial_initial_data():\n");
  // printf("%15s %15s %15s\n", "r_axis", "rho_baryon", "pressure");
  // printf("-------------------------------------------------------------\n");
  // for (int i = 0; i < num_radial_pts; i++) {
  //   printf("%15.6e %15.6e %15.6e\n", r_axis[i], rho_baryon[i], pressure[i]);
  // }
  /** END DEBUGGING **/

  // Free memory allocated for ID struct
  {
    free(ID_persist.r_Schw_arr);
    free(ID_persist.rho_energy_arr);
    free(ID_persist.rho_baryon_arr);
    free(ID_persist.P_arr);
    free(ID_persist.M_arr);
    free(ID_persist.expnu_arr);
    free(ID_persist.exp4phi_arr);
    free(ID_persist.r_iso_arr);
  }
}