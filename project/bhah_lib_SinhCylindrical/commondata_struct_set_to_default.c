#include "BHaH_defines.h"
/**
 * Set commondata_struct to default values specified within NRPy+.
 */
void commondata_struct_set_to_default(commondata_struct *restrict commondata) {

  // Set commondata_struct variables to default
  commondata->CFL_FACTOR = 0.5;                         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::CFL_FACTOR
  commondata->NUMGRIDS = 1;                             // nrpy.grid::NUMGRIDS
  commondata->checkpoint_every = 2.0;                   // nrpy.infrastructures.BHaH.checkpointing::checkpoint_every
  commondata->convergence_factor = 1.0;                 // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
  commondata->diagnostics_output_every = 100;           // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::diagnostics_output_every
  commondata->eta = 1.0;                                // nrpy.equations.general_relativity.BSSN_gauge_RHSs::eta
  commondata->initial_central_density = 0.129285;       // TOVola::initial_central_density
  commondata->initial_ode_step_size = 1e-20;            // TOVola::initial_ode_step_size
  commondata->interpolation_stencil_size = 11;          // TOVola::interpolation_stencil_size
  commondata->max_interpolation_stencil_size = 13;      // TOVola::max_interpolation_stencil_size
  commondata->max_step_size = 1.0;                      // TOVola::max_step_size
  commondata->min_step_size = 1e-20;                    // TOVola::min_step_size
  commondata->ode_error_limit = 1e-06;                  // TOVola::ode_error_limit
  commondata->ode_max_steps = 5000000;                  // TOVola::ode_max_steps
  commondata->output_progress_every = 1;                // nrpy.infrastructures.BHaH.diagnostics.progress_indicator::output_progress_every
  commondata->poly_eos_Gamma = 2.0;                     // TOVola::poly_eos_Gamma
  commondata->poly_eos_K = 1.0;                         // TOVola::poly_eos_K
  commondata->t_final = 10.0;                           // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_final
  snprintf(commondata->outer_bc_type, 50, "radiation"); // nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions::outer_bc_type
}
