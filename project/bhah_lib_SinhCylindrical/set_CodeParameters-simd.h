const REAL NOSIMDAMPLRHO = params->AMPLRHO;                                              // nrpy.reference_metric_SinhCylindrical::AMPLRHO
MAYBE_UNUSED const REAL_SIMD_ARRAY AMPLRHO = ConstSIMD(NOSIMDAMPLRHO);                   // nrpy.reference_metric_SinhCylindrical::AMPLRHO
const REAL NOSIMDAMPLZ = params->AMPLZ;                                                  // nrpy.reference_metric_SinhCylindrical::AMPLZ
MAYBE_UNUSED const REAL_SIMD_ARRAY AMPLZ = ConstSIMD(NOSIMDAMPLZ);                       // nrpy.reference_metric_SinhCylindrical::AMPLZ
const REAL NOSIMDCart_originx = params->Cart_originx;                                    // nrpy.grid::Cart_originx
MAYBE_UNUSED const REAL_SIMD_ARRAY Cart_originx = ConstSIMD(NOSIMDCart_originx);         // nrpy.grid::Cart_originx
const REAL NOSIMDCart_originy = params->Cart_originy;                                    // nrpy.grid::Cart_originy
MAYBE_UNUSED const REAL_SIMD_ARRAY Cart_originy = ConstSIMD(NOSIMDCart_originy);         // nrpy.grid::Cart_originy
const REAL NOSIMDCart_originz = params->Cart_originz;                                    // nrpy.grid::Cart_originz
MAYBE_UNUSED const REAL_SIMD_ARRAY Cart_originz = ConstSIMD(NOSIMDCart_originz);         // nrpy.grid::Cart_originz
const REAL NOSIMDCFL_FACTOR = commondata->CFL_FACTOR;                                    // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::CFL_FACTOR
MAYBE_UNUSED const REAL_SIMD_ARRAY CFL_FACTOR = ConstSIMD(NOSIMDCFL_FACTOR);             // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::CFL_FACTOR
const REAL NOSIMDcheckpoint_every = commondata->checkpoint_every;                        // nrpy.infrastructures.BHaH.checkpointing::checkpoint_every
MAYBE_UNUSED const REAL_SIMD_ARRAY checkpoint_every = ConstSIMD(NOSIMDcheckpoint_every); // nrpy.infrastructures.BHaH.checkpointing::checkpoint_every
const REAL NOSIMDconvergence_factor = commondata->convergence_factor; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
MAYBE_UNUSED const REAL_SIMD_ARRAY convergence_factor =
    ConstSIMD(NOSIMDconvergence_factor);                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::convergence_factor
MAYBE_UNUSED const int CoordSystem_hash = params->CoordSystem_hash; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::CoordSystem_hash
const REAL NOSIMDdiagnostics_output_every =
    commondata->diagnostics_output_every; // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::diagnostics_output_every
MAYBE_UNUSED const REAL_SIMD_ARRAY diagnostics_output_every =
    ConstSIMD(NOSIMDdiagnostics_output_every); // nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library::diagnostics_output_every
const REAL NOSIMDdt = commondata->dt;          // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::dt
MAYBE_UNUSED const REAL_SIMD_ARRAY dt = ConstSIMD(NOSIMDdt);               // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::dt
const REAL NOSIMDdxx0 = params->dxx0;                                      // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx0
MAYBE_UNUSED const REAL_SIMD_ARRAY dxx0 = ConstSIMD(NOSIMDdxx0);           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx0
const REAL NOSIMDdxx1 = params->dxx1;                                      // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx1
MAYBE_UNUSED const REAL_SIMD_ARRAY dxx1 = ConstSIMD(NOSIMDdxx1);           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx1
const REAL NOSIMDdxx2 = params->dxx2;                                      // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx2
MAYBE_UNUSED const REAL_SIMD_ARRAY dxx2 = ConstSIMD(NOSIMDdxx2);           // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::dxx2
const REAL NOSIMDeta = commondata->eta;                                    // nrpy.equations.general_relativity.BSSN_gauge_RHSs::eta
MAYBE_UNUSED const REAL_SIMD_ARRAY eta = ConstSIMD(NOSIMDeta);             // nrpy.equations.general_relativity.BSSN_gauge_RHSs::eta
const REAL NOSIMDf0_of_xx0 = params->f0_of_xx0;                            // nrpy.reference_metric_SinhCylindrical::f0_of_xx0
MAYBE_UNUSED const REAL_SIMD_ARRAY f0_of_xx0 = ConstSIMD(NOSIMDf0_of_xx0); // nrpy.reference_metric_SinhCylindrical::f0_of_xx0
const REAL NOSIMDf1_of_xx1 = params->f1_of_xx1;                            // nrpy.reference_metric_SinhCylindrical::f1_of_xx1
MAYBE_UNUSED const REAL_SIMD_ARRAY f1_of_xx1 = ConstSIMD(NOSIMDf1_of_xx1); // nrpy.reference_metric_SinhCylindrical::f1_of_xx1
const REAL NOSIMDf2_of_xx0 = params->f2_of_xx0;                            // nrpy.reference_metric_SinhCylindrical::f2_of_xx0
MAYBE_UNUSED const REAL_SIMD_ARRAY f2_of_xx0 = ConstSIMD(NOSIMDf2_of_xx0); // nrpy.reference_metric_SinhCylindrical::f2_of_xx0
const REAL NOSIMDf3_of_xx2 = params->f3_of_xx2;                            // nrpy.reference_metric_SinhCylindrical::f3_of_xx2
MAYBE_UNUSED const REAL_SIMD_ARRAY f3_of_xx2 = ConstSIMD(NOSIMDf3_of_xx2); // nrpy.reference_metric_SinhCylindrical::f3_of_xx2
const REAL NOSIMDf4_of_xx1 = params->f4_of_xx1;                            // nrpy.reference_metric_SinhCylindrical::f4_of_xx1
MAYBE_UNUSED const REAL_SIMD_ARRAY f4_of_xx1 = ConstSIMD(NOSIMDf4_of_xx1); // nrpy.reference_metric_SinhCylindrical::f4_of_xx1
const REAL NOSIMDgrid_hole_radius = params->grid_hole_radius;              // nrpy.reference_metric::grid_hole_radius
MAYBE_UNUSED const REAL_SIMD_ARRAY grid_hole_radius = ConstSIMD(NOSIMDgrid_hole_radius);               // nrpy.reference_metric::grid_hole_radius
const REAL NOSIMDgrid_physical_size = params->grid_physical_size;                                      // nrpy.reference_metric::grid_physical_size
MAYBE_UNUSED const REAL_SIMD_ARRAY grid_physical_size = ConstSIMD(NOSIMDgrid_physical_size);           // nrpy.reference_metric::grid_physical_size
MAYBE_UNUSED const bool grid_rotates = params->grid_rotates;                                           // nrpy.grid::grid_rotates
const REAL NOSIMDinitial_central_density = commondata->initial_central_density;                        // TOVola::initial_central_density
MAYBE_UNUSED const REAL_SIMD_ARRAY initial_central_density = ConstSIMD(NOSIMDinitial_central_density); // TOVola::initial_central_density
const REAL NOSIMDinitial_ode_step_size = commondata->initial_ode_step_size;                            // TOVola::initial_ode_step_size
MAYBE_UNUSED const REAL_SIMD_ARRAY initial_ode_step_size = ConstSIMD(NOSIMDinitial_ode_step_size);     // TOVola::initial_ode_step_size
MAYBE_UNUSED const int interpolation_stencil_size = commondata->interpolation_stencil_size;            // TOVola::interpolation_stencil_size
const REAL NOSIMDinvdxx0 = params->invdxx0;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx0
MAYBE_UNUSED const REAL_SIMD_ARRAY invdxx0 = ConstSIMD(NOSIMDinvdxx0); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx0
const REAL NOSIMDinvdxx1 = params->invdxx1;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx1
MAYBE_UNUSED const REAL_SIMD_ARRAY invdxx1 = ConstSIMD(NOSIMDinvdxx1); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx1
const REAL NOSIMDinvdxx2 = params->invdxx2;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx2
MAYBE_UNUSED const REAL_SIMD_ARRAY invdxx2 = ConstSIMD(NOSIMDinvdxx2); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::invdxx2
MAYBE_UNUSED const int max_interpolation_stencil_size = commondata->max_interpolation_stencil_size; // TOVola::max_interpolation_stencil_size
const REAL NOSIMDmax_step_size = commondata->max_step_size;                                         // TOVola::max_step_size
MAYBE_UNUSED const REAL_SIMD_ARRAY max_step_size = ConstSIMD(NOSIMDmax_step_size);                  // TOVola::max_step_size
const REAL NOSIMDmin_step_size = commondata->min_step_size;                                         // TOVola::min_step_size
MAYBE_UNUSED const REAL_SIMD_ARRAY min_step_size = ConstSIMD(NOSIMDmin_step_size);                  // TOVola::min_step_size
MAYBE_UNUSED const int nn = commondata->nn;                             // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::nn
MAYBE_UNUSED const int nn_0 = commondata->nn_0;                         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::nn_0
MAYBE_UNUSED const int NUMGRIDS = commondata->NUMGRIDS;                 // nrpy.grid::NUMGRIDS
MAYBE_UNUSED const int Nxx0 = params->Nxx0;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx0
MAYBE_UNUSED const int Nxx1 = params->Nxx1;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx1
MAYBE_UNUSED const int Nxx2 = params->Nxx2;                             // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx2
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS0
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS1
MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2; // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::Nxx_plus_2NGHOSTS2
const REAL NOSIMDode_error_limit = commondata->ode_error_limit;         // TOVola::ode_error_limit
MAYBE_UNUSED const REAL_SIMD_ARRAY ode_error_limit = ConstSIMD(NOSIMDode_error_limit); // TOVola::ode_error_limit
MAYBE_UNUSED const int ode_max_steps = commondata->ode_max_steps;                      // TOVola::ode_max_steps
const REAL NOSIMDPI = params->PI;                                                      // nrpy.reference_metric::PI
MAYBE_UNUSED const REAL_SIMD_ARRAY PI = ConstSIMD(NOSIMDPI);                           // nrpy.reference_metric::PI
const REAL NOSIMDpoly_eos_Gamma = commondata->poly_eos_Gamma;                          // TOVola::poly_eos_Gamma
MAYBE_UNUSED const REAL_SIMD_ARRAY poly_eos_Gamma = ConstSIMD(NOSIMDpoly_eos_Gamma);   // TOVola::poly_eos_Gamma
const REAL NOSIMDpoly_eos_K = commondata->poly_eos_K;                                  // TOVola::poly_eos_K
MAYBE_UNUSED const REAL_SIMD_ARRAY poly_eos_K = ConstSIMD(NOSIMDpoly_eos_K);           // TOVola::poly_eos_K
const REAL NOSIMDSINHWRHO = params->SINHWRHO;                                          // nrpy.reference_metric_SinhCylindrical::SINHWRHO
MAYBE_UNUSED const REAL_SIMD_ARRAY SINHWRHO = ConstSIMD(NOSIMDSINHWRHO);               // nrpy.reference_metric_SinhCylindrical::SINHWRHO
const REAL NOSIMDSINHWZ = params->SINHWZ;                                              // nrpy.reference_metric_SinhCylindrical::SINHWZ
MAYBE_UNUSED const REAL_SIMD_ARRAY SINHWZ = ConstSIMD(NOSIMDSINHWZ);                   // nrpy.reference_metric_SinhCylindrical::SINHWZ
const REAL NOSIMDt_0 = commondata->t_0;                                                // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_0
MAYBE_UNUSED const REAL_SIMD_ARRAY t_0 = ConstSIMD(NOSIMDt_0);                         // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_0
const REAL NOSIMDt_final = commondata->t_final;                                        // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_final
MAYBE_UNUSED const REAL_SIMD_ARRAY t_final = ConstSIMD(NOSIMDt_final);                 // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::t_final
const REAL NOSIMDtime = commondata->time;                                              // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::time
MAYBE_UNUSED const REAL_SIMD_ARRAY time = ConstSIMD(NOSIMDtime);                       // nrpy.infrastructures.BHaH.MoLtimestepping.MoL::time
const REAL NOSIMDxxmax0 = params->xxmax0;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax0
MAYBE_UNUSED const REAL_SIMD_ARRAY xxmax0 = ConstSIMD(NOSIMDxxmax0); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax0
const REAL NOSIMDxxmax1 = params->xxmax1;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax1
MAYBE_UNUSED const REAL_SIMD_ARRAY xxmax1 = ConstSIMD(NOSIMDxxmax1); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax1
const REAL NOSIMDxxmax2 = params->xxmax2;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax2
MAYBE_UNUSED const REAL_SIMD_ARRAY xxmax2 = ConstSIMD(NOSIMDxxmax2); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmax2
const REAL NOSIMDxxmin0 = params->xxmin0;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin0
MAYBE_UNUSED const REAL_SIMD_ARRAY xxmin0 = ConstSIMD(NOSIMDxxmin0); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin0
const REAL NOSIMDxxmin1 = params->xxmin1;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin1
MAYBE_UNUSED const REAL_SIMD_ARRAY xxmin1 = ConstSIMD(NOSIMDxxmin1); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin1
const REAL NOSIMDxxmin2 = params->xxmin2;                            // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin2
MAYBE_UNUSED const REAL_SIMD_ARRAY xxmin2 = ConstSIMD(NOSIMDxxmin2); // nrpy.infrastructures.BHaH.numerical_grids_and_timestep::xxmin2
