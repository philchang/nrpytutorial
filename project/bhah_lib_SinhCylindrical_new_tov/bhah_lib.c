#include "bhah_lib.h"
#include "uniform_Lagrange_interp_3D.h"

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

#ifndef BHAH_DATA_COMPONENTS
#define BHAH_DATA_COMPONENTS 8
#endif

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

  // TODO: loop through particles once to identify hydro grid bounds

  // Populate Tmunu
  LOOP_OMP("omp parallel for",
          i0, NGHOSTS, Nxx_plus_2NGHOSTS0 - NGHOSTS,
          i1, NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS,
          i2, NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS) {

    // Compute Cartesian coordinate of the grid point
    REAL xCart[3] = {0.0, 0.0, 0.0};
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    const REAL xx2 = xx[2][i2];

    // Declare variables for hydro particles (poison data for debugging)
    REAL rho = NAN;
    REAL vx = NAN;
    REAL vy = NAN;
    REAL vz = NAN;
    REAL ie = NAN;
    REAL xHydro = NAN;
    REAL yHydro = NAN;
    REAL zHydro = NAN;

    // Initialize distance square to a very high value
    REAL dist2 = 1.0e300;

    // Loop all particles to perform a nearest neighbor search
    for (int p = 0; p < nParticles; p++) {
      int offset = p * BHAH_DATA_COMPONENTS;
      const REAL _rho = prims_and_pos[offset + 0];
      const REAL _vx = prims_and_pos[offset + 1];
      const REAL _vy = prims_and_pos[offset + 2];
      const REAL _vz = prims_and_pos[offset + 3];
      const REAL _ie = prims_and_pos[offset + 4];
      xHydro = prims_and_pos[offset + 5];
      yHydro = prims_and_pos[offset + 6];
      zHydro = prims_and_pos[offset + 7];

      const REAL d2 = (
        (xCart[0] - xHydro) * (xCart[0] - xHydro) +
        (xCart[1] - yHydro) * (xCart[1] - yHydro) +
        (xCart[2] - zHydro) * (xCart[2] - zHydro)
      );

      // Assign primitive variables to the values of the nearest hydro particle
      if (d2 < dist2) {
        dist2 = d2;
        rho = _rho;
        vx = _vx;
        vy = _vy;
        vz = _vz;
        ie = _ie;
      }  // END if (d2 < dist2)
    }  // END for (int p = 0; p < nParticles; p++)

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

  // // Load primitives and positions of the hydro particles
  // FILE *fp = fopen("bhah_prims_and_pos_dump.dat", "w");
  // for (int p = 0; p < nParticles; p++) {
  //   int offset = p * BHAH_DATA_COMPONENTS;

  //   const double rho = prims_and_pos[offset + 0];
  //   const double vx = prims_and_pos[offset + 1];
  //   const double vy = prims_and_pos[offset + 2];
  //   const double vz = prims_and_pos[offset + 3];
  //   const double ie = prims_and_pos[offset + 4];
  //   const double x = prims_and_pos[offset + 5];
  //   const double y = prims_and_pos[offset + 6];
  //   const double z = prims_and_pos[offset + 7];

  //   fprintf(fp, "%g %g %g %g %g %g %g %g\n", rho, vx, vy, vz, ie, x, y, z);
  // }
  // fclose(fp);
  // printf("[BHaH_set_PrimsAndPos_particles] wrote %d entries to bhah_prims_and_pos_dump.dat\n", nParticles);
}  // END void BHaH_set_PrimsAndPos_particles()
