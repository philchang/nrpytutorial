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
