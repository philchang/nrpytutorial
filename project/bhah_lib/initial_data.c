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
void manga_radial_initial_data(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {

  // Set single grid index
  const int grid =0;
  params_struct *restrict params = &griddata[grid].params;


  // Manually sample radial grid for MaNGa
  const int num_radial_pts = 100;
  const REAL RMAX = params->RMAX;
  const REAL local_SINHW = 0.2;
  REAL *restrict r_axis = (REAL *restrict)malloc(sizeof(REAL) * num_radial_pts);
  const REAL dx = 1.0 / ((REAL) num_radial_pts);
  for (int i=0; i < num_radial_pts; i++){
    const x = dx / 2.0 + dx * i;
    r_axis[i] = RMAX * sinh(x / local_SINHW) / sinh(1.0 / local_SINHW);
  }
  REAL *restrict rho_baryon = (REAL *restrict)malloc(sizeof(REAL) * num_radial_pts);
  REAL *restrict pressure = (REAL *restrict)malloc(sizeof(REAL) * num_radial_pts);


  ID_persist_struct ID_persist;


  for (int i=0; i < num_radial_pts; i++){
    const x = dx / 2.0 + dx * i;
    r_axis[i] = RMAX * sinh(x / local_SINHW) / sinh(1.0 / local_SINHW);
    TOVola_radial_only_interp(commondata, params, r_axis[i], &ID_persist,  &rho_baryon[i], &pressure[i]);
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