const REAL Cart_originx = params.Cart_originx;                 // nrpy.grid::Cart_originx
const REAL Cart_originy = params.Cart_originy;                 // nrpy.grid::Cart_originy
const REAL Cart_originz = params.Cart_originz;                 // nrpy.grid::Cart_originz
const REAL CFL_FACTOR = commondata.CFL_FACTOR;                 // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::CFL_FACTOR
const REAL checkpoint_every = commondata.checkpoint_every;     // nrpy.infrastructures.BHaH.checkpointing::checkpoint_every
const REAL convergence_factor = commondata.convergence_factor; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
const int CoordSystem_hash = params.CoordSystem_hash;          // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::CoordSystem_hash
char CoordSystemName[50];                                      // nrpy.reference_metric::CoordSystemName
{
  // Copy up to 49 characters from params.CoordSystemName to CoordSystemName
  strncpy(CoordSystemName, params.CoordSystemName, 50 - 1);
  // Explicitly null-terminate CoordSystemName to ensure it is a valid C-string
  CoordSystemName[50 - 1] = '\0'; // Properly null terminate char array.
}
const REAL diagnostics_output_every =
    commondata.diagnostics_output_every; // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::diagnostics_output_every
const REAL dt = commondata.dt;           // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::dt
const REAL dxx0 = params.dxx0;           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx0
const REAL dxx1 = params.dxx1;           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx1
const REAL dxx2 = params.dxx2;           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx2
const REAL eta = commondata.eta;         // nrpy.equations.general_relativity.BSSN_gauge_RHSs::eta
const REAL f0_of_xx0 = params.f0_of_xx0; // nrpy.reference_metric_Spherical::f0_of_xx0
const REAL f1_of_xx1 = params.f1_of_xx1; // nrpy.reference_metric_Spherical::f1_of_xx1
const REAL f2_of_xx0 = params.f2_of_xx0; // nrpy.reference_metric_Spherical::f2_of_xx0
const REAL f3_of_xx2 = params.f3_of_xx2; // nrpy.reference_metric_Spherical::f3_of_xx2
const REAL f4_of_xx1 = params.f4_of_xx1; // nrpy.reference_metric_Spherical::f4_of_xx1
const REAL grid_hole_radius = params.grid_hole_radius;                        // nrpy.reference_metric::grid_hole_radius
const REAL grid_physical_size = params.grid_physical_size;                    // nrpy.reference_metric::grid_physical_size
const bool grid_rotates = params.grid_rotates;                                // nrpy.grid::grid_rotates
const REAL initial_central_density = commondata.initial_central_density;      // TOVola::initial_central_density
const REAL initial_ode_step_size = commondata.initial_ode_step_size;          // TOVola::initial_ode_step_size
const int interpolation_stencil_size = commondata.interpolation_stencil_size; // TOVola::interpolation_stencil_size
const REAL invdxx0 = params.invdxx0;                                          // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx0
const REAL invdxx1 = params.invdxx1;                                          // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx1
const REAL invdxx2 = params.invdxx2;                                          // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx2
const int max_interpolation_stencil_size = commondata.max_interpolation_stencil_size; // TOVola::max_interpolation_stencil_size
const REAL max_step_size = commondata.max_step_size;                                  // TOVola::max_step_size
const REAL min_step_size = commondata.min_step_size;                                  // TOVola::min_step_size
const int nn = commondata.nn;                                                         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::nn
const int nn_0 = commondata.nn_0;                                                     // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::nn_0
const int NUMGRIDS = commondata.NUMGRIDS;                                             // nrpy.grid::NUMGRIDS
const int Nxx0 = params.Nxx0;                                                         // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx0
const int Nxx1 = params.Nxx1;                                                         // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx1
const int Nxx2 = params.Nxx2;                                                         // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx2
const int Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS0
const int Nxx_plus_2NGHOSTS1 = params.Nxx_plus_2NGHOSTS1; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS1
const int Nxx_plus_2NGHOSTS2 = params.Nxx_plus_2NGHOSTS2; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS2
const REAL ode_error_limit = commondata.ode_error_limit;  // TOVola::ode_error_limit
const int ode_max_steps = commondata.ode_max_steps;       // TOVola::ode_max_steps
char outer_bc_type[50];                                   // nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions::outer_bc_type
{
  // Copy up to 49 characters from commondata.outer_bc_type to outer_bc_type
  strncpy(outer_bc_type, commondata.outer_bc_type, 50 - 1);
  // Explicitly null-terminate outer_bc_type to ensure it is a valid C-string
  outer_bc_type[50 - 1] = '\0'; // Properly null terminate char array.
}
const REAL PI = params.PI;                             // nrpy.reference_metric::PI
const REAL poly_eos_Gamma = commondata.poly_eos_Gamma; // TOVola::poly_eos_Gamma
const REAL poly_eos_K = commondata.poly_eos_K;         // TOVola::poly_eos_K
const REAL RMAX = params.RMAX;                         // nrpy.reference_metric_Spherical::RMAX
const REAL t_0 = commondata.t_0;                       // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_0
const REAL t_final = commondata.t_final;               // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_final
const REAL time = commondata.time;                     // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::time
const REAL xxmax0 = params.xxmax0;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax0
const REAL xxmax1 = params.xxmax1;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax1
const REAL xxmax2 = params.xxmax2;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax2
const REAL xxmin0 = params.xxmin0;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin0
const REAL xxmin1 = params.xxmin1;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin1
const REAL xxmin2 = params.xxmin2;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin2
