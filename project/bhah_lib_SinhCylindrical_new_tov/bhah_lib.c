#include "bhah_lib.h"
#include "uniform_Lagrange_interp_3D.h"
#include <limits.h>
#include <time.h>

#define IDX3_TO_ijk(idx, i, j, k)                            \
  do {                                                       \
    (k) = (idx) / (Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1); \
    (j) = ((idx) / Nxx_plus_2NGHOSTS0) % Nxx_plus_2NGHOSTS1; \
    (i) = (idx) % Nxx_plus_2NGHOSTS0;                        \
  } while (0)

/**
 * List of functions that need to be defined.
 * Some of these functions will be merged and/or deprecated.
 *
 * void BHaH_setup(const int nx1, const int nx2, const int nx3,
 *                 const REAL cfl, const REAL rmax, BHaH_struct *bhahstruct);
 *
 * void BHaH_initialize( BHaH_struct *bhahstruct);
 *
 * double BHaH_get_timestep( BHaH_struct *bhahstruct); <- this one should be deprecated
 *
 * void BHaH_advance_timestep( const REAL t_final, BHaH_struct *bhahstruct);
 *
 * void BHaH_get_metric_extrinsic_curvature(REAL Cartx, REAL Carty, REAL Cartz,
 *                                          BHaH_struct *bhahstruct, REAL *alpha_out,
 *                                          REAL *beta0, REAL *beta1, REAL *beta2,
 *                                          REAL (*gammaDD)[3], REAL (*KDD)[3]);
 *
 * void BHaH_set_Tmunu(const int rgrid, const REAL rmax, REAL *in_Tmunu, BHaH_struct *bhahstruct);
 *
 *     NOTE: This function is no longer used and will not be included in the new library.
 *
 * BHaH_output_file(char *filename, BHaH_struct *bhahstruct, int dim);
 *
 * int BHaH_get_gridpoints(int *indices, REAL *xCartGrid, REAL xCartMax[3], BHaH_struct *bhahstruct);
 *
 * void BHaH_set_Tmunu_gridpoints(const int nCartGrid, int *indices, REAL *TmunuGrid, BHaH_struct *bhahstruct);
 **/

void BHaH_setup(const int nxx0, const int nxx1, const int nxx2, const REAL cfl,
                const REAL rmax, const REAL sinhwrho, const REAL sinhwz,
                BHaH_struct *bhahstruct) {

  // Step 1.a: Allocate memory for the commondata struct
  commondata_struct *commondata = (commondata_struct *)malloc(sizeof(commondata_struct));

  // Step 1.b: Set each commondata CodeParameter to default.
  commondata_struct_set_to_default(commondata);

  // Step 1.c: Allocate NUMGRIDS griddata arrays, each containing
  //           data specific to an individual grid.
  const int n_grids = commondata->NUMGRIDS;
  griddata_struct *griddata = (griddata_struct *)malloc(sizeof(griddata_struct) * n_grids);

  // Step 1.d: Set each CodeParameter in griddata.params to default.
  params_struct_set_to_default(commondata, griddata);

  // Thiago says: Overwrite parameters for single grid
  const int grid = 0;

  // Thiago says: Overwrite grid size based on rmax parameter (to be set by MaNGa)
  (griddata[grid].params).AMPLRHO = rmax;
  (griddata[grid].params).AMPLZ = rmax;
  (griddata[grid].params).grid_physical_size = rmax;

  // Thiago says: Overwrite SINH* parameters (to be set by MaNGa)
  (griddata[grid].params).SINHWRHO = sinhwrho;
  (griddata[grid].params).SINHWZ = sinhwz;

  // Thiago says: Overwrite Nxx0, Nxx1, and Nxx2 for the single grid
  (griddata[grid].params).Nxx0 = nxx0;
  (griddata[grid].params).Nxx1 = nxx1;
  (griddata[grid].params).Nxx2 = nxx2;

  // Thiago says: Overwrite CFL value:
  commondata->CFL_FACTOR = cfl;

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
  bhahstruct->griddata = griddata;
  bhahstruct->commondata = commondata;
}

void BHaH_get_metric_extrinsic_curvature(
    BHaH_struct *bhahstruct,
    REAL Cartx, REAL Carty, REAL Cartz,
    REAL *alpha_out, REAL *beta0, REAL *beta1, REAL *beta2,
    REAL gammaDD[3][3], REAL KDD[3][3]) {
  // Unpack commondata & griddata
  commondata_struct *commondata = bhahstruct->commondata;
  griddata_struct *restrict griddata = bhahstruct->griddata;
  const int grid = 0; // single grid
  params_struct *restrict params = &griddata[grid].params;

  // Set grid coordinates
  REAL *restrict xx[3] = {griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2]};

  // Set evolution grid functions
  REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;

#include "set_CodeParameters.h"

  // Map Cartesian -> reference-metric
  REAL xCart[3] = {Cartx, Carty, Cartz};
  REAL xx012[3];
  int i012[3];
  Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart, xx012, i012);

  // interpolation setup
  const int num_interp_gfs = 18;
  static const int list_of_interp_gfs[18] = {
      ALPHAGF, CFGF, TRKGF, VETU0GF, VETU1GF, VETU2GF,
      HDD00GF, HDD01GF, HDD02GF, HDD11GF, HDD12GF, HDD22GF,
      ADD00GF, ADD01GF, ADD02GF, ADD11GF, ADD12GF, ADD22GF};
  const int num_interp_pts = 1;
  const int N_interp_GHOSTS = 3; // interpolation order = 2 * N_interp_GHOSTS + 1
  const REAL src_dxx0 = dxx0;
  const REAL src_dxx1 = dxx1;
  const REAL src_dxx2 = dxx2;

  // Single destination point array
  REAL dst_x0x1x2[1][3] = {{xx012[0], xx012[1], xx012[2]}};

  // Allocate interpolation buffers
  REAL interp_output[num_interp_gfs * num_interp_pts];
  const REAL *restrict src_gf_ptrs[num_interp_gfs];
  REAL *restrict dst_data[num_interp_gfs];

  // Compute total pts per GF
  int Nxx_plus_2NGHOSTS_TOTAL = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  // Build pointers into source and destination buffers
  for (int gf = 0; gf < num_interp_gfs; gf++) {
    int which = list_of_interp_gfs[gf];
    src_gf_ptrs[gf] = &y_n_gfs[which * Nxx_plus_2NGHOSTS_TOTAL];
    dst_data[gf] = &interp_output[gf * num_interp_pts];
  }

  // Call the new interpolator
  int err = interpolation_3d_general__uniform_src_grid(
      N_interp_GHOSTS,
      src_dxx0, src_dxx1, src_dxx2,
      Nxx_plus_2NGHOSTS0,
      Nxx_plus_2NGHOSTS1,
      Nxx_plus_2NGHOSTS2,
      num_interp_gfs,
      xx,
      src_gf_ptrs,
      num_interp_pts,
      dst_x0x1x2,
      dst_data);
  if (err != 0) {
    fprintf(stderr, "BHaH_get_metric_extrinsic_curvature: interpolation error %d\n", err);
    exit(1);
  }

  // Unpack interpolated values at a single point
  REAL alpha = dst_data[0][0];
  REAL cf = dst_data[1][0];
  REAL trK = dst_data[2][0];
  REAL vetU0 = dst_data[3][0];
  REAL vetU1 = dst_data[4][0];
  REAL vetU2 = dst_data[5][0];
  REAL hDD00 = dst_data[6][0];
  REAL hDD01 = dst_data[7][0];
  REAL hDD02 = dst_data[8][0];
  REAL hDD11 = dst_data[9][0];
  REAL hDD12 = dst_data[10][0];
  REAL hDD22 = dst_data[11][0];
  REAL aDD00 = dst_data[12][0];
  REAL aDD01 = dst_data[13][0];
  REAL aDD02 = dst_data[14][0];
  REAL aDD11 = dst_data[15][0];
  REAL aDD12 = dst_data[16][0];
  REAL aDD22 = dst_data[17][0];

  // Update lapse
  *alpha_out = alpha;

  // Transform to Cartesian and compute extrinsic curvature
#include "NRPY+unrescale+basis_transform_to_Cartesian_metric.h"
}

void BHaH_output_file(char *filename, BHaH_struct *bhahstruct, int dim) {

  // Define local structures
  commondata_struct *commondata = bhahstruct->commondata;

  // Define param struct for a single grid
  griddata_struct *griddata = bhahstruct->griddata;
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;
  // Define rfm_struct for a single grid
  const rfm_struct *restrict rfmstruct = &griddata[grid].rfmstruct;

  // Define grid points
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  // Define evolution grid functions
  REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
  // Define auxiliary grid functions
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  // Define diagnostics grid functions
  REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;

#include "set_CodeParameters.h"

  // Evaluate Hamiltonian and momentum constraints
  constraints_eval(commondata, params, rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);

  FILE *outfile = fopen(filename, "w");

  /**
   * Thiago says: FIXME: TOV_Mass macro is not defined.
   *
   * Thiago says: TODO: Output momentum constraint as well.
   *
   */
  const REAL TOV_Mass = 1.0; /** FIXME: set it to 1 for now **/

  LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
              NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
              NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {
    const int idx = IDX3(i0, i1, i2);
    REAL xCart[3];
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
    if (dim == 2) {
      fprintf(outfile, "%e %e %e %e\n",
              xCart[1] / TOV_Mass, xCart[2] / TOV_Mass,
              y_n_gfs[IDX4pt(CFGF, idx)], log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx)])));
    } else {
      fprintf(outfile, "%e %e %e %e %e\n",
              xCart[0] / TOV_Mass, xCart[1] / TOV_Mass, xCart[2] / TOV_Mass,
              y_n_gfs[IDX4pt(CFGF, idx)], log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx)])));
    }
  }
  fclose(outfile);
}

int BHaH_get_gridpoints(int *indices, REAL *xCartGrid, REAL xCartMax[3], BHaH_struct *bhahstruct) {

  // Define local structures
  commondata_struct *commondata = bhahstruct->commondata;

  // Define param struct for a single grid
  griddata_struct *griddata = bhahstruct->griddata;
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;

  // Define grid points
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

#include "set_CodeParameters.h"

  int nCartGrid = 0;
  // THIAGO says: the line below was being incorrectly set
  // LOOP_REGION(0, Nxx_plus_2NGHOSTS0, 0, Nxx_plus_2NGHOSTS0, 0, Nxx_plus_2NGHOSTS2) {

  // THIAGO SAYS: Unrelated to the issue above, let us try only sending points from the interior grid
  // LOOP_REGION(0, Nxx_plus_2NGHOSTS0, 0, Nxx_plus_2NGHOSTS1, 0, Nxx_plus_2NGHOSTS2) {
  LOOP_REGION(NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS, NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS, NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {
    const int idx = IDX3(i0, i1, i2);
    REAL xCart[3];
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

    if (fabs(xCart[0]) > xCartMax[0] || fabs(xCart[1]) > xCartMax[1] || fabs(xCart[2]) > xCartMax[2])
      continue;

    xCartGrid[3 * nCartGrid + 0] = xCart[0];
    xCartGrid[3 * nCartGrid + 1] = xCart[1];
    xCartGrid[3 * nCartGrid + 2] = xCart[2];
    indices[nCartGrid++] = idx;
  }
  return nCartGrid;
}

/**
 * Thiago says: NOTE: TMUNU_AVG_COMP is also defined in GR_Utils.h
 */
#define TMUNU_COMP 10

void BHaH_set_Tmunu_gridpoints(const int nCartGrid, int *indices, REAL *TmunuGrid, BHaH_struct *bhahstruct) {

  // Define local structures
  commondata_struct *commondata = bhahstruct->commondata;

  // Define param struct for a single grid
  griddata_struct *griddata = bhahstruct->griddata;
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;

  // Define grid points
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Define auxiliary grid functions
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;

#include "set_CodeParameters.h"

  for (int i = 0; i < nCartGrid; i++) {
    const int idx = indices[i];
    int i0 = -10, i1 = -10, i2 = -10;
    IDX3_TO_ijk(idx, i0, i1, i2);

    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    const REAL xx2 = xx[2][i2];

    // Set components of Tmunu in Cartesian coordinates
    REAL T4CartUU00 = TmunuGrid[i * TMUNU_COMP + 0];
    REAL T4CartUU01 = TmunuGrid[i * TMUNU_COMP + 1];
    REAL T4CartUU02 = TmunuGrid[i * TMUNU_COMP + 2];
    REAL T4CartUU03 = TmunuGrid[i * TMUNU_COMP + 3];
    REAL T4CartUU11 = TmunuGrid[i * TMUNU_COMP + 4];
    REAL T4CartUU12 = TmunuGrid[i * TMUNU_COMP + 5];
    REAL T4CartUU13 = TmunuGrid[i * TMUNU_COMP + 6];
    REAL T4CartUU22 = TmunuGrid[i * TMUNU_COMP + 7];
    REAL T4CartUU23 = TmunuGrid[i * TMUNU_COMP + 8];
    REAL T4CartUU33 = TmunuGrid[i * TMUNU_COMP + 9];

    // Compute components of Tmunu in reference-metric coordinates
    REAL T4UU00, T4UU01, T4UU02, T4UU03, T4UU11, T4UU12, T4UU13, T4UU22, T4UU23, T4UU33;
#include "transform_T4UU_from_Cart_to_SinhCylindrical.h"

    auxevol_gfs[IDX4(T4UU00GF, i0, i1, i2)] = T4UU00;
    auxevol_gfs[IDX4(T4UU01GF, i0, i1, i2)] = T4UU01;
    auxevol_gfs[IDX4(T4UU02GF, i0, i1, i2)] = T4UU02;
    auxevol_gfs[IDX4(T4UU03GF, i0, i1, i2)] = T4UU03;
    auxevol_gfs[IDX4(T4UU11GF, i0, i1, i2)] = T4UU11;
    auxevol_gfs[IDX4(T4UU12GF, i0, i1, i2)] = T4UU12;
    auxevol_gfs[IDX4(T4UU13GF, i0, i1, i2)] = T4UU13;
    auxevol_gfs[IDX4(T4UU22GF, i0, i1, i2)] = T4UU22;
    auxevol_gfs[IDX4(T4UU23GF, i0, i1, i2)] = T4UU23;
    auxevol_gfs[IDX4(T4UU33GF, i0, i1, i2)] = T4UU33;
  }
}

void BHaH_set_TOV(BHaH_struct *bhahstruct, const int num_radial_pts, REAL *restrict r_axis, REAL *restrict rho_baryon, REAL *restrict pressure) {

  // Define local structures
  commondata_struct *commondata = bhahstruct->commondata;

  // Define griddata struct
  griddata_struct *griddata = bhahstruct->griddata;

  manga_radial_initial_data(commondata, griddata, num_radial_pts, r_axis, rho_baryon, pressure);
}

void BHaH_write_checkpoint(BHaH_struct *bhahstruct, char *filename) {
  // Define local structures
  commondata_struct *commondata = bhahstruct->commondata;
  griddata_struct *griddata = bhahstruct->griddata;

  // Write checkpoint using default function, which creates a file named 'checkpoint-conv_factor1.00.dat'
  write_checkpoint(commondata, griddata);

  // Rename file to the parameter 'filename'
  const char *default_filename = "checkpoint-conv_factor1.00.dat";
  if (rename(default_filename, filename) != 0) {
    perror("BHaH_write_checkpoint: ERROR renaming checkpoint file");
    exit(EXIT_FAILURE);
  }
}

// The value of the number of components is hard coded at this point, and must match the same in BHaH_interface.h in MaNGa
// TODO: consolidate into one definition
#ifndef BHAH_DATA_COMPONENTS
#define BHAH_DATA_COMPONENTS 23
#endif

// void BHaH_set_PrimsAndPos_particles(int nParticles, const double *prims_and_pos, BHaH_struct *bhahstruct) {

//   // Define local structures
//   commondata_struct *commondata = bhahstruct->commondata;

//   // Define param struct for a single grid
//   griddata_struct *griddata = bhahstruct->griddata;
//   const int grid = 0;
//   params_struct *restrict params = &griddata[grid].params;

//   // Define grid points
//   REAL *restrict xx[3];
//   for (int ww = 0; ww < 3; ww++)
//     xx[ww] = griddata[grid].xx[ww];

//   // Set up grid functions
//   REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
//   REAL *restrict y_n = griddata[grid].gridfuncs.y_n_gfs;

// #include "set_CodeParameters.h"

//   // Loop through particles once to identify hydro grid bounds and atmosphere values
//   REAL MAXHYDROSIZE = -1.0;
//   REAL rho_min = 1.0e100;
//   REAL ie_min  = 1.0e100;
//   for (int p = 0; p < nParticles; p++) {
//     int base_idx = p * BHAH_DATA_COMPONENTS;
//     const REAL _rho = prims_and_pos[base_idx + 0];
//     const REAL _ie = prims_and_pos[base_idx + 4];
//     const REAL xHydro = prims_and_pos[base_idx + 5];
//     // Since hydro grid is a square box centered at the origin, we only need to check one dimension
//     if (xHydro > MAXHYDROSIZE)
//       MAXHYDROSIZE = xHydro;
//     if (_rho < rho_min)
//       rho_min = _rho;
//     if (_ie < ie_min)
//       ie_min = _ie;
//   }  // END for (int p = 0; p < nParticles; p++)

//   // Adjust hydro size to 95%
//   MAXHYDROSIZE = 0.95 * MAXHYDROSIZE;

//   // BHaH loop to populate Tmunu
//   LOOP_OMP("omp parallel for",
//           i0, NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
//           i1, NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
//           i2, NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {

//     // Set curvilinear coordinates
//     const REAL xx0 = xx[0][i0];
//     const REAL xx1 = xx[1][i1];
//     const REAL xx2 = xx[2][i2];

//     // Compute Cartesian coordinate of the grid point
//     REAL xCart[3] = {0.0, 0.0, 0.0};
//     xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

//     // Declare variables for hydro particles (initialize to NANs to catch any errors)
//     REAL rho = NAN;
//     REAL vx = NAN;
//     REAL vy = NAN;
//     REAL vz = NAN;
//     REAL ie = NAN;

//     // Set primitive values to atmosphere if outside hydro grid
//     if ((fabs(xCart[0]) > MAXHYDROSIZE) || (fabs(xCart[1]) > MAXHYDROSIZE) || (fabs(xCart[2]) > MAXHYDROSIZE)) {
//       rho = rho_min;
//       vx = 0.0;
//       vy = 0.0;
//       vz = 0.0;
//       ie = ie_min;
//     }
//     else {
//       // Initialize distance square to a very high value
//       REAL dist2 = 1.0e300;

//       // Loop all particles to perform a nearest neighbor search
//       for (int p = 0; p < nParticles; p++) {
//         int base_idx = p * BHAH_DATA_COMPONENTS;
//         const REAL _rho = prims_and_pos[base_idx + 0];
//         const REAL _vx = prims_and_pos[base_idx + 1];
//         const REAL _vy = prims_and_pos[base_idx + 2];
//         const REAL _vz = prims_and_pos[base_idx + 3];
//         const REAL _ie = prims_and_pos[base_idx + 4];
//         const REAL xHydro = prims_and_pos[base_idx + 5];
//         const REAL yHydro = prims_and_pos[base_idx + 6];
//         const REAL zHydro = prims_and_pos[base_idx + 7];

//         const REAL d2 = (
//           (xCart[0] - xHydro) * (xCart[0] - xHydro) +
//           (xCart[1] - yHydro) * (xCart[1] - yHydro) +
//           (xCart[2] - zHydro) * (xCart[2] - zHydro)
//         );

//         // Assign primitive variables to the values of the nearest hydro particle
//         if (d2 < dist2) {
//           dist2 = d2;
//           rho = _rho;
//           vx = _vx;
//           vy = _vy;
//           vz = _vz;
//           ie = _ie;
//         }  // END if (d2 < dist2)
//       }  // END for (int p = 0; p < nParticles; p++)
//     }  // END else

//     // Compute corresponding pressure
//     const REAL press = (commondata->poly_eos_Gamma - 1.0) * rho * ie;

//     // Compute rescaled 3-velocity in curvilinear coordinates
//     const REAL vCartU[3] = {vx, vy, vz};
//     REAL rescaledvU0, rescaledvU1, rescaledvU2;
//     compute_rescaledvU_from_vCartU__rfm__SinhCylindrical(commondata, params, vCartU, xx0, xx1, xx2,
//                                                          &rescaledvU0, &rescaledvU1, &rescaledvU2);

//     // Compute time component of 4-velocity and restrict 3-velocity
//     REAL u4Ut;
//     const REAL max_lorentz_factor = 10.0;
//     compute_u4Ut__rfm__SinhCylindrical(commondata, params, max_lorentz_factor, i0, i1, i2, y_n,
//                                        &rescaledvU0, &rescaledvU1, &rescaledvU2, &u4Ut);

//     // Compute specific enthalpy
//     const REAL h = 1.0 + ie + press / rho;

//     // Compute TUU[4][4] in SinhCylindrical coordinates at point (i0, i1, i2)
//     compute_T4UU(commondata, params, i0, i1, i2, xx, rho, press, h, u4Ut,
//                  rescaledvU0, rescaledvU1, rescaledvU2, y_n, auxevol_gfs);

//   } // END LOOP_OMP

//   // // Load primitives and positions of the hydro particles
//   // FILE *fp = fopen("bhah_prims_and_pos_dump.dat", "w");
//   // for (int p = 0; p < nParticles; p++) {
//   //   int base_idx = p * BHAH_DATA_COMPONENTS;

//   //   const double rho = prims_and_pos[base_idx + 0];
//   //   const double vx = prims_and_pos[base_idx + 1];
//   //   const double vy = prims_and_pos[base_idx + 2];
//   //   const double vz = prims_and_pos[base_idx + 3];
//   //   const double ie = prims_and_pos[base_idx + 4];
//   //   const double x = prims_and_pos[base_idx + 5];
//   //   const double y = prims_and_pos[base_idx + 6];
//   //   const double z = prims_and_pos[base_idx + 7];

//   //   fprintf(fp, "%g %g %g %g %g %g %g %g\n", rho, vx, vy, vz, ie, x, y, z);
//   // }
//   // fclose(fp);
//   // printf("[BHaH_set_PrimsAndPos_particles] wrote %d entries to bhah_prims_and_pos_dump.dat\n", nParticles);
// }  // END void BHaH_set_PrimsAndPos_particles()

/*
 * =====================================================================
 * Function: BHaH_build_particle_bins
 * ---------------------------------------------------------------------
 * Purpose:
 *   Organize hydro particles into a 3D spatial grid of bins for efficient
 *   nearest-neighbor queries.
 *
 * Description:
 *   This function implements a simple "binned spatial search" structure:
 *   the simulation domain [xmin,xmax) × [ymin,ymax) × [zmin,zmax) is
 *   subdivided into Nx × Ny × Nz uniform bins. Each bin stores the list
 *   of particle indices whose positions fall inside its bounds.
 *
 *   Internally, the bin data are stored in Compressed Sparse Row (CSR)
 *   format:
 *     - `offsets[b]` and `offsets[b+1]` mark the range of particles
 *       belonging to bin b inside the `indices[]` array.
 *     - `indices[pos]` contains the particle index (0 ≤ index < nParticles).
 *
 *   The build process occurs in two passes:
 *     1. Pass 1: Count how many particles fall into each bin (histogram).
 *     2. Pass 2: Compute prefix sums (offsets) and fill the indices array.
 *
 *   Once built, this structure allows fast neighborhood lookups by
 *   scanning only the bin containing a query point and its immediate
 *   neighbors (typically 27 bins total), instead of all particles.
 *
 * Inputs:
 *   nParticles     - total number of particles.
 *   prims_and_pos  - array of primitive variables + positions, where
 *                    positions are stored at offsets +5,+6,+7.
 *   xmin,xmax,...  - spatial domain limits.
 *   Nx,Ny,Nz       - number of bins per dimension.
 *
 * Outputs:
 *   out_indices    - pointer to allocated array of particle indices (CSR data).
 *   out_offsets    - pointer to allocated array of offsets (CSR row boundaries).
 *
 * Returns:
 *   0 on success, nonzero on error.
 * =====================================================================
 */

int BHaH_build_particle_bins(const int nParticles,
                             const REAL *restrict prims_and_pos,
                             const REAL xmin, const REAL xmax,
                             const REAL ymin, const REAL ymax,
                             const REAL zmin, const REAL zmax,
                             const int Nx, const int Ny, const int Nz,
                             int **out_indices,
                             int **out_offsets) {
  if (prims_and_pos == NULL || out_indices == NULL || out_offsets == NULL)
    return 1;
  if (nParticles < 0 || Nx < 1 || Ny < 1 || Nz < 1)
    return 2;

  const REAL Lx = xmax - xmin;
  const REAL Ly = ymax - ymin;
  const REAL Lz = zmax - zmin;
  if (!(Lx > 0 && Ly > 0 && Lz > 0))
    return 3;

  const int Nbins = Nx * Ny * Nz;
  const REAL dx = Lx / (REAL)Nx;
  const REAL dy = Ly / (REAL)Ny;
  const REAL dz = Lz / (REAL)Nz;

  int *counts = (int *)calloc((size_t)Nbins, sizeof(int));
  if (counts == NULL)
    return 4;

  printf("[BHaH_build_particle_bins]: starting Pass 1: histogram per bin...\n");

  /* Pass 1: histogram per bin */
  for (int p = 0; p < nParticles; p++) {
    const int base_idx = p * BHAH_DATA_COMPONENTS;
    const REAL x = prims_and_pos[base_idx + 5];
    const REAL y = prims_and_pos[base_idx + 6];
    const REAL z = prims_and_pos[base_idx + 7];

    int ix = (int)floor((x - xmin) / dx);
    int iy = (int)floor((y - ymin) / dy);
    int iz = (int)floor((z - zmin) / dz);

    if (ix < 0)
      ix = 0;
    else if (ix >= Nx)
      ix = Nx - 1;
    if (iy < 0)
      iy = 0;
    else if (iy >= Ny)
      iy = Ny - 1;
    if (iz < 0)
      iz = 0;
    else if (iz >= Nz)
      iz = Nz - 1;

    const int b = ix + Nx * (iy + Ny * iz);
    counts[b] += 1;
  }

  /* Build CSR offsets (exclusive prefix sum) */
  int *offsets = (int *)malloc((size_t)(Nbins + 1) * sizeof(int));
  if (offsets == NULL) {
    free(counts);
    return 5;
  }

  offsets[0] = 0;
  for (int b = 0; b < Nbins; b++) {
    offsets[b + 1] = offsets[b] + counts[b];
  }
  if (offsets[Nbins] != nParticles) {
    free(counts);
    free(offsets);
    return 6;
  }

  /* Allocate indices; reuse counts as write cursors */
  int *indices = (int *)malloc((size_t)nParticles * sizeof(int));
  if (indices == NULL) {
    free(counts);
    free(offsets);
    return 7;
  }
  for (int b = 0; b < Nbins; b++)
    counts[b] = offsets[b];

  printf("[BHaH_build_particle_bins]: Pass 2: scatter particle indices...\n");

  /* Pass 2: scatter particle indices */
  for (int p = 0; p < nParticles; p++) {
    const int base_idx = p * BHAH_DATA_COMPONENTS;
    const REAL x = prims_and_pos[base_idx + 5];
    const REAL y = prims_and_pos[base_idx + 6];
    const REAL z = prims_and_pos[base_idx + 7];

    int ix = (int)floor((x - xmin) / dx);
    int iy = (int)floor((y - ymin) / dy);
    int iz = (int)floor((z - zmin) / dz);

    if (ix < 0)
      ix = 0;
    else if (ix >= Nx)
      ix = Nx - 1;
    if (iy < 0)
      iy = 0;
    else if (iy >= Ny)
      iy = Ny - 1;
    if (iz < 0)
      iz = 0;
    else if (iz >= Nz)
      iz = Nz - 1;

    const int b = ix + Nx * (iy + Ny * iz);
    const int pos = counts[b]++;
    indices[pos] = p;
  }

  free(counts);
  *out_indices = indices;
  *out_offsets = offsets;

  printf("[BHaH_build_particle_bins]: finished building bins\n");
  return 0;
} // END BHaH_build_particle_bins

/*
 * =====================================================================
 * Function: BHaH_find_nearest_particle
 * ---------------------------------------------------------------------
 * Purpose:
 *   Given a Cartesian query position xCart = {x,y,z}, find the index of
 *   the nearest hydro particle using the CSR bins produced by
 *   BHaH_build_particle_bins().
 *
 * Description:
 *   - The domain [xmin,xmax)×[ymin,ymax)×[zmin,zmax) is partitioned into
 *     a uniform grid of Nx×Ny×Nz bins. Particles are stored in CSR:
 *       * offsets[b] .. offsets[b+1]-1 are indices into `indices[]`
 *       * indices[...] holds particle IDs p (0..nParticles-1)
 *       * particle p position is prims_and_pos[p*BHAH_DATA_COMPONENTS + {5,6,7}]
 *   - Starting from the query’s own bin, an expanding “cube” search over
 *     neighboring bins is performed. After each expansion, an exact
 *     stopping test compares the best point distance found so far against
 *     a lower bound given by the next shell’s bin AABBs. If no unseen bin
 *     can contain a closer particle, the search terminates with the exact
 *     nearest neighbor.
 *
 * Returns:
 *   0 on success; nonzero error code otherwise.
 * =====================================================================
 */

/* Branchless helper: squared distance from xCart to the axis-aligned bounding box of bin (i,j,k). */
static inline REAL BHaH_bin_aabb_dist2(const REAL xCart[3],
                                       const REAL xmin, const REAL ymin, const REAL zmin,
                                       const REAL dx, const REAL dy, const REAL dz,
                                       const int i, const int j, const int k) {
  const REAL xc = xCart[0], yc = xCart[1], zc = xCart[2];

  const REAL xlo = xmin + (REAL)i * dx, xhi = xlo + dx;
  const REAL ylo = ymin + (REAL)j * dy, yhi = ylo + dy;
  const REAL zlo = zmin + (REAL)k * dz, zhi = zlo + dz;

  /* Branchless distance to interval [lo, hi] */
  const REAL dx1 = fmax(xlo - xc, (REAL)0.0);
  const REAL dx2 = fmax(dx1, xc - xhi);
  const REAL dy1 = fmax(ylo - yc, (REAL)0.0);
  const REAL dy2 = fmax(dy1, yc - yhi);
  const REAL dz1 = fmax(zlo - zc, (REAL)0.0);
  const REAL dz2 = fmax(dz1, zc - zhi);

  return dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
}

int BHaH_find_nearest_particle(const REAL xCart[3],
                               const REAL xmin, const REAL xmax,
                               const REAL ymin, const REAL ymax,
                               const REAL zmin, const REAL zmax,
                               const int Nx, const int Ny, const int Nz,
                               const int nParticles,
                               const int *restrict offsets, /* length: Nx*Ny*Nz + 1 */
                               const int *restrict indices, /* length: nParticles */
                               const REAL *restrict prims_and_pos,
                               int *out_pid,
                               REAL *out_dist2) {
  if (xCart == NULL || offsets == NULL || indices == NULL ||
      prims_and_pos == NULL || out_pid == NULL)
    return 1;
  if (nParticles <= 0 || Nx < 1 || Ny < 1 || Nz < 1)
    return 2;

  const REAL Lx = xmax - xmin;
  const REAL Ly = ymax - ymin;
  const REAL Lz = zmax - zmin;
  if (!(Lx > 0 && Ly > 0 && Lz > 0))
    return 3;

  const REAL dx = Lx / (REAL)Nx;
  const REAL dy = Ly / (REAL)Ny;
  const REAL dz = Lz / (REAL)Nz;

  /* Map query to its bin; clamp to domain. */
  int iq = (int)floor((xCart[0] - xmin) / dx);
  int jq = (int)floor((xCart[1] - ymin) / dy);
  int kq = (int)floor((xCart[2] - zmin) / dz);
  if (iq < 0)
    iq = 0;
  else if (iq >= Nx)
    iq = Nx - 1;
  if (jq < 0)
    jq = 0;
  else if (jq >= Ny)
    jq = Ny - 1;
  if (kq < 0)
    kq = 0;
  else if (kq >= Nz)
    kq = Nz - 1;

  /* Fixed 3x3x3 neighborhood (clamped). */
  const int i_min = (iq - 1) < 0 ? 0 : (iq - 1);
  const int i_max = (iq + 1) > Nx - 1 ? Nx - 1 : (iq + 1);
  const int j_min = (jq - 1) < 0 ? 0 : (jq - 1);
  const int j_max = (jq + 1) > Ny - 1 ? Ny - 1 : (jq + 1);
  const int k_min = (kq - 1) < 0 ? 0 : (kq - 1);
  const int k_max = (kq + 1) > Nz - 1 ? Nz - 1 : (kq + 1);

  const REAL xc = xCart[0], yc = xCart[1], zc = xCart[2];

  REAL best_r2 = 1.0e300;
  int best_pid = -1;

  for (int k = k_min; k <= k_max; k++) {
    for (int j = j_min; j <= j_max; j++) {
      for (int i = i_min; i <= i_max; i++) {
        const int b = i + Nx * (j + Ny * k);
        const int begin = offsets[b];
        const int end = offsets[b + 1];
        if (begin == end)
          continue;

        for (int pos = begin; pos < end; ++pos) {
          const int pid = indices[pos];
          const int base_idx = pid * BHAH_DATA_COMPONENTS;
          const REAL px = prims_and_pos[base_idx + 5];
          const REAL py = prims_and_pos[base_idx + 6];
          const REAL pz = prims_and_pos[base_idx + 7];

          const REAL dxp = xc - px;
          const REAL dyp = yc - py;
          const REAL dzp = zc - pz;
          const REAL r2 = dxp * dxp + dyp * dyp + dzp * dzp;

          if (r2 < best_r2) {
            best_r2 = r2;
            best_pid = pid;
            if (best_r2 == 0.0)
              goto BH_NEAREST_DONE; /* exact match */
          }
        }
      }
    }
  }

BH_NEAREST_DONE:
  if (best_pid < 0)
    return 4;

  *out_pid = best_pid;
  if (out_dist2)
    *out_dist2 = best_r2;
  return 0;
}
// END BHaH_find_nearest_particle

// void BHaH_set_PrimsAndPos_particles(int nParticles, const double *prims_and_pos, BHaH_struct *bhahstruct) {

//   printf("[bhah_lib]: Skipping BHaH_set_PrimsAndPos_particles\n");
// }

void BHaH_set_PrimsAndPos_particles(int nParticles, const double *prims_and_pos, BHaH_struct *bhahstruct) {

  // Define local structures
  commondata_struct *commondata = bhahstruct->commondata;

  // Define param struct for a single grid
  griddata_struct *griddata = bhahstruct->griddata;
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;

  // Define grid points
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Set up grid functions
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  REAL *restrict y_n = griddata[grid].gridfuncs.y_n_gfs;

#include "set_CodeParameters.h"

  // Loop through particles once to identify hydro grid bounds and atmosphere values
  REAL MAXHYDROSIZE = -1.0;
  REAL rho_min = 1.0e100;
  REAL ie_min = 1.0e100;
  for (int p = 0; p < nParticles; p++) {
    int base_idx = p * BHAH_DATA_COMPONENTS;
    const REAL _rho = prims_and_pos[base_idx + BHAH_DATA_RHO];
    const REAL _ie = prims_and_pos[base_idx + BHAH_DATA_IE];
    const REAL xHydro = prims_and_pos[base_idx + BHAH_DATA_X];
    // Since hydro grid is a square box centered at the origin, we only need to check one dimension
    if (xHydro > MAXHYDROSIZE)
      MAXHYDROSIZE = xHydro;
    if (_rho < rho_min)
      rho_min = _rho;
    if (_ie < ie_min)
      ie_min = _ie;
  } // END for (int p = 0; p < nParticles; p++)

  // // Adjust hydro size to 95%
  // MAXHYDROSIZE = 0.95 * MAXHYDROSIZE;
  // Adjust hydro size to 3.0
  MAXHYDROSIZE = 3.0;

  /* === Insert after computing MAXHYDROSIZE, rho_min, ie_min === */
  int *bin_indices = NULL, *bin_offsets = NULL;
  /* Use cubic hydro domain centered at origin, sized by MAXHYDROSIZE */
  const REAL xmin_bins = -MAXHYDROSIZE, xmax_bins = MAXHYDROSIZE;
  const REAL ymin_bins = -MAXHYDROSIZE, ymax_bins = MAXHYDROSIZE;
  const REAL zmin_bins = -MAXHYDROSIZE, zmax_bins = MAXHYDROSIZE;
  /* Tunable bin resolution */
  const int Nx_bins = 50, Ny_bins = 50, Nz_bins = 50;
  {
    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    const int rc = BHaH_build_particle_bins(nParticles, prims_and_pos,
                                            xmin_bins, xmax_bins,
                                            ymin_bins, ymax_bins,
                                            zmin_bins, zmax_bins,
                                            Nx_bins, Ny_bins, Nz_bins,
                                            &bin_indices, &bin_offsets);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double elapsed_sec = (t1.tv_sec - t0.tv_sec) + 1.0e-9 * (t1.tv_nsec - t0.tv_nsec);
    printf("[Timing] BHaH_build_particle_bins: %.6f seconds\n", elapsed_sec);
    if (rc != 0) {
      printf("[BHaH_set_PrimsAndPos_particles] ERROR: failed to build particle bins (rc=%d)\n", rc);
      exit(EXIT_FAILURE);
    }
  }

  // BHaH loop to populate Tmunu
  struct timespec t0, t1;
  clock_gettime(CLOCK_MONOTONIC, &t0);
  LOOP_OMP("omp parallel for",
           i0, NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
           i1, NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
           i2, NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {

    // Set curvilinear coordinates
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    const REAL xx2 = xx[2][i2];

    // Compute Cartesian coordinate of the grid point
    REAL xCart[3] = {0.0, 0.0, 0.0};
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

    // Declare variables for hydro particles (initialize to NANs to catch any errors)
    REAL rho = NAN;
    REAL vx = NAN;
    REAL vy = NAN;
    REAL vz = NAN;
    REAL ie = NAN;

    // Set primitive values to atmosphere if outside hydro grid
    if ((fabs(xCart[0]) > MAXHYDROSIZE) || (fabs(xCart[1]) > MAXHYDROSIZE) || (fabs(xCart[2]) > MAXHYDROSIZE)) {
      rho = rho_min;
      vx = 0.0;
      vy = 0.0;
      vz = 0.0;
      ie = ie_min;
    } else {
      int pid;
      REAL dist2;

      const int rc_nn = BHaH_find_nearest_particle(
          xCart,
          xmin_bins, xmax_bins,
          ymin_bins, ymax_bins,
          zmin_bins, zmax_bins,
          Nx_bins, Ny_bins, Nz_bins,
          nParticles,
          bin_offsets, bin_indices,
          prims_and_pos,
          &pid, &dist2);

      if (rc_nn != 0) {
        printf("[BHaH_set_PrimsAndPos_particles] WARNING: nearest particle not found (rc_nn=%d)\n", rc_nn);
        printf("[BHaH_set_PrimsAndPos_particles] setting primitives to atmosphere at (x, y, z) = (%.2e, %.2e, %.2e)\n",
               xCart[0], xCart[1], xCart[2]);
        rho = rho_min;
        vx = 0.0;
        vy = 0.0;
        vz = 0.0;
        ie = ie_min;
      } else {
        const int base_idx = pid * BHAH_DATA_COMPONENTS;
        /* Primitive values at particle centre */
        const REAL rho0 = prims_and_pos[base_idx + BHAH_DATA_RHO];
        const REAL vx0 = prims_and_pos[base_idx + BHAH_DATA_VX];
        const REAL vy0 = prims_and_pos[base_idx + BHAH_DATA_VY];
        const REAL vz0 = prims_and_pos[base_idx + BHAH_DATA_VZ];
        const REAL ie0 = prims_and_pos[base_idx + BHAH_DATA_IE];

        /* Particle position */
        const REAL px = prims_and_pos[base_idx + BHAH_DATA_X];
        const REAL py = prims_and_pos[base_idx + BHAH_DATA_Y];
        const REAL pz = prims_and_pos[base_idx + BHAH_DATA_Z];

        /* Offset from particle to gridpoint (in Cartesian coords) */
        const REAL dx = xCart[0] - px;
        const REAL dy = xCart[1] - py;
        const REAL dz = xCart[2] - pz;

        /* Gradients d(q)/dx^i at the particle */
        const REAL gradrho_x = prims_and_pos[base_idx + BHAH_DATA_GRADRHO_X];
        const REAL gradrho_y = prims_and_pos[base_idx + BHAH_DATA_GRADRHO_Y];
        const REAL gradrho_z = prims_and_pos[base_idx + BHAH_DATA_GRADRHO_Z];

        const REAL gradvx_x = prims_and_pos[base_idx + BHAH_DATA_GRADVX_X];
        const REAL gradvx_y = prims_and_pos[base_idx + BHAH_DATA_GRADVX_Y];
        const REAL gradvx_z = prims_and_pos[base_idx + BHAH_DATA_GRADVX_Z];

        const REAL gradvy_x = prims_and_pos[base_idx + BHAH_DATA_GRADVY_X];
        const REAL gradvy_y = prims_and_pos[base_idx + BHAH_DATA_GRADVY_Y];
        const REAL gradvy_z = prims_and_pos[base_idx + BHAH_DATA_GRADVY_Z];

        const REAL gradvz_x = prims_and_pos[base_idx + BHAH_DATA_GRADVZ_X];
        const REAL gradvz_y = prims_and_pos[base_idx + BHAH_DATA_GRADVZ_Y];
        const REAL gradvz_z = prims_and_pos[base_idx + BHAH_DATA_GRADVZ_Z];

        const REAL gradie_x = prims_and_pos[base_idx + BHAH_DATA_GRADIE_X];
        const REAL gradie_y = prims_and_pos[base_idx + BHAH_DATA_GRADIE_Y];
        const REAL gradie_z = prims_and_pos[base_idx + BHAH_DATA_GRADIE_Z];

        /* First–order Taylor reconstruction at the gridpoint */
        rho = rho0 + gradrho_x * dx + gradrho_y * dy + gradrho_z * dz;
        vx = vx0 + gradvx_x * dx + gradvx_y * dy + gradvx_z * dz;
        vy = vy0 + gradvy_x * dx + gradvy_y * dy + gradvy_z * dz;
        vz = vz0 + gradvz_x * dx + gradvz_y * dy + gradvz_z * dz;
        ie = ie0 + gradie_x * dx + gradie_y * dy + gradie_z * dz;

      } // END else [after if (rc_nn != 0) ]
    } // END else [after if ((fabs(xCart[0]) > MAXHYDROSIZE) || (fabs(xCart[1]) > MAXHYDROSIZE) || (fabs(xCart[2]) > MAXHYDROSIZE))]

    // Compute corresponding pressure
    const REAL press = (commondata->poly_eos_Gamma - 1.0) * rho * ie;

    // Compute rescaled 3-velocity in curvilinear coordinates
    const REAL vCartU[3] = {vx, vy, vz};
    REAL rescaledvU0, rescaledvU1, rescaledvU2;
    compute_rescaledvU_from_vCartU__rfm__SinhCylindrical(commondata, params, vCartU, xx0, xx1, xx2,
                                                         &rescaledvU0, &rescaledvU1, &rescaledvU2);

    // Compute time component of 4-velocity and restrict 3-velocity
    REAL u4Ut;
    const REAL max_lorentz_factor = 10.0;
    compute_u4Ut__rfm__SinhCylindrical(commondata, params, max_lorentz_factor, i0, i1, i2, y_n,
                                       &rescaledvU0, &rescaledvU1, &rescaledvU2, &u4Ut);

    // Compute specific enthalpy
    const REAL h = 1.0 + ie + press / rho;

    // Compute TUU[4][4] in SinhCylindrical coordinates at point (i0, i1, i2)
    compute_T4UU(commondata, params, i0, i1, i2, xx, rho, press, h, u4Ut,
                 rescaledvU0, rescaledvU1, rescaledvU2, y_n, auxevol_gfs);

  } // END LOOP_OMP
  clock_gettime(CLOCK_MONOTONIC, &t1);
  double elapsed_sec = (t1.tv_sec - t0.tv_sec) + 1.0e-9 * (t1.tv_nsec - t0.tv_nsec);
  printf("[Timing] LOOP_OMP: %.6f seconds\n", elapsed_sec);
  free(bin_indices);
  free(bin_offsets);

  // // Load primitives and positions of the hydro particles
  // FILE *fp = fopen("bhah_prims_and_pos_dump.dat", "w");
  // for (int p = 0; p < nParticles; p++) {
  //   int base_idx = p * BHAH_DATA_COMPONENTS;

  //   const double rho = prims_and_pos[base_idx + 0];
  //   const double vx = prims_and_pos[base_idx + 1];
  //   const double vy = prims_and_pos[base_idx + 2];
  //   const double vz = prims_and_pos[base_idx + 3];
  //   const double ie = prims_and_pos[base_idx + 4];
  //   const double x = prims_and_pos[base_idx + 5];
  //   const double y = prims_and_pos[base_idx + 6];
  //   const double z = prims_and_pos[base_idx + 7];

  //   fprintf(fp, "%g %g %g %g %g %g %g %g\n", rho, vx, vy, vz, ie, x, y, z);
  // }
  // fclose(fp);
  // printf("[BHaH_set_PrimsAndPos_particles] wrote %d entries to bhah_prims_and_pos_dump.dat\n", nParticles);

} // END void BHaH_set_PrimsAndPos_particles()
