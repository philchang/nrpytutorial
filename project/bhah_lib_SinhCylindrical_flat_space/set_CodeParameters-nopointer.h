MAYBE_UNUSED const REAL AMPLRHO = params.AMPLRHO;                       // nrpy.reference_metric_SinhCylindrical::AMPLRHO
MAYBE_UNUSED const REAL AMPLZ = params.AMPLZ;                           // nrpy.reference_metric_SinhCylindrical::AMPLZ
MAYBE_UNUSED const REAL Cart_originx = params.Cart_originx;             // nrpy.grid::Cart_originx
MAYBE_UNUSED const REAL Cart_originy = params.Cart_originy;             // nrpy.grid::Cart_originy
MAYBE_UNUSED const REAL Cart_originz = params.Cart_originz;             // nrpy.grid::Cart_originz
MAYBE_UNUSED const REAL CFL_FACTOR = commondata.CFL_FACTOR;             // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::CFL_FACTOR
MAYBE_UNUSED const REAL checkpoint_every = commondata.checkpoint_every; // nrpy.infrastructures.BHaH.checkpointing::checkpoint_every
MAYBE_UNUSED const REAL convergence_factor =
    commondata.convergence_factor;                                 // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
MAYBE_UNUSED const int CoordSystem_hash = params.CoordSystem_hash; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::CoordSystem_hash
char CoordSystemName[50];                                          // nrpy.reference_metric::CoordSystemName
{
  // Copy up to 49 characters from params.CoordSystemName to CoordSystemName
  strncpy(CoordSystemName, params.CoordSystemName, 50 - 1);
  // Explicitly null-terminate CoordSystemName to ensure it is a valid C-string
  CoordSystemName[50 - 1] = '\0'; // Properly null terminate char array.
}
MAYBE_UNUSED const REAL diagnostics_output_every =
    commondata.diagnostics_output_every;              // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::diagnostics_output_every
MAYBE_UNUSED const REAL dt = commondata.dt;           // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::dt
MAYBE_UNUSED const REAL dxx0 = params.dxx0;           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx0
MAYBE_UNUSED const REAL dxx1 = params.dxx1;           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx1
MAYBE_UNUSED const REAL dxx2 = params.dxx2;           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx2
MAYBE_UNUSED const REAL eta = commondata.eta;         // nrpy.equations.general_relativity.BSSN_gauge_RHSs::eta
MAYBE_UNUSED const REAL f0_of_xx0 = params.f0_of_xx0; // nrpy.reference_metric_SinhCylindrical::f0_of_xx0
MAYBE_UNUSED const REAL f1_of_xx1 = params.f1_of_xx1; // nrpy.reference_metric_SinhCylindrical::f1_of_xx1
MAYBE_UNUSED const REAL f2_of_xx0 = params.f2_of_xx0; // nrpy.reference_metric_SinhCylindrical::f2_of_xx0
MAYBE_UNUSED const REAL f3_of_xx2 = params.f3_of_xx2; // nrpy.reference_metric_SinhCylindrical::f3_of_xx2
MAYBE_UNUSED const REAL f4_of_xx1 = params.f4_of_xx1; // nrpy.reference_metric_SinhCylindrical::f4_of_xx1
MAYBE_UNUSED const REAL grid_hole_radius = params.grid_hole_radius;                        // nrpy.reference_metric::grid_hole_radius
MAYBE_UNUSED const REAL grid_physical_size = params.grid_physical_size;                    // nrpy.reference_metric::grid_physical_size
MAYBE_UNUSED const bool grid_rotates = params.grid_rotates;                                // nrpy.grid::grid_rotates
MAYBE_UNUSED const REAL initial_central_density = commondata.initial_central_density;      // TOVola::initial_central_density
MAYBE_UNUSED const REAL initial_ode_step_size = commondata.initial_ode_step_size;          // TOVola::initial_ode_step_size
MAYBE_UNUSED const int interpolation_stencil_size = commondata.interpolation_stencil_size; // TOVola::interpolation_stencil_size
MAYBE_UNUSED const REAL invdxx0 = params.invdxx0; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx0
MAYBE_UNUSED const REAL invdxx1 = params.invdxx1; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx1
MAYBE_UNUSED const REAL invdxx2 = params.invdxx2; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx2
MAYBE_UNUSED const int max_interpolation_stencil_size = commondata.max_interpolation_stencil_size; // TOVola::max_interpolation_stencil_size
MAYBE_UNUSED const REAL max_step_size = commondata.max_step_size;                                  // TOVola::max_step_size
MAYBE_UNUSED const REAL min_step_size = commondata.min_step_size;                                  // TOVola::min_step_size
MAYBE_UNUSED const int nn = commondata.nn;                             // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::nn
MAYBE_UNUSED const int nn_0 = commondata.nn_0;                         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::nn_0
MAYBE_UNUSED const int NUMGRIDS = commondata.NUMGRIDS;                 // nrpy.grid::NUMGRIDS
MAYBE_UNUSED const int Nxx0 = params.Nxx0;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx0
MAYBE_UNUSED const int Nxx1 = params.Nxx1;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx1
MAYBE_UNUSED const int Nxx2 = params.Nxx2;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx2
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params.Nxx_plus_2NGHOSTS0; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS0
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params.Nxx_plus_2NGHOSTS1; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS1
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params.Nxx_plus_2NGHOSTS2; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS2
MAYBE_UNUSED const REAL ode_error_limit = commondata.ode_error_limit;  // TOVola::ode_error_limit
MAYBE_UNUSED const int ode_max_steps = commondata.ode_max_steps;       // TOVola::ode_max_steps
char outer_bc_type[50]; // nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions::outer_bc_type
{
  // Copy up to 49 characters from commondata.outer_bc_type to outer_bc_type
  strncpy(outer_bc_type, commondata.outer_bc_type, 50 - 1);
  // Explicitly null-terminate outer_bc_type to ensure it is a valid C-string
  outer_bc_type[50 - 1] = '\0'; // Properly null terminate char array.
}
MAYBE_UNUSED const REAL PI = params.PI;                             // nrpy.reference_metric::PI
MAYBE_UNUSED const REAL poly_eos_Gamma = commondata.poly_eos_Gamma; // TOVola::poly_eos_Gamma
MAYBE_UNUSED const REAL poly_eos_K = commondata.poly_eos_K;         // TOVola::poly_eos_K
MAYBE_UNUSED const REAL SINHWRHO = params.SINHWRHO;                 // nrpy.reference_metric_SinhCylindrical::SINHWRHO
MAYBE_UNUSED const REAL SINHWZ = params.SINHWZ;                     // nrpy.reference_metric_SinhCylindrical::SINHWZ
MAYBE_UNUSED const REAL t_0 = commondata.t_0;                       // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_0
MAYBE_UNUSED const REAL t_final = commondata.t_final;               // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_final
MAYBE_UNUSED const REAL time = commondata.time;                     // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::time
MAYBE_UNUSED const REAL xxmax0 = params.xxmax0;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax0
MAYBE_UNUSED const REAL xxmax1 = params.xxmax1;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax1
MAYBE_UNUSED const REAL xxmax2 = params.xxmax2;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax2
MAYBE_UNUSED const REAL xxmin0 = params.xxmin0;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin0
MAYBE_UNUSED const REAL xxmin1 = params.xxmin1;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin1
MAYBE_UNUSED const REAL xxmin2 = params.xxmin2;                     // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin2
