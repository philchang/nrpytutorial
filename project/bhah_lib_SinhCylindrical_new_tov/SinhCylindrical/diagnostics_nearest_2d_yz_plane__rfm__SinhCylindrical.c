#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Output diagnostic quantities at gridpoints closest to yz plane.
 */
void diagnostics_nearest_2d_yz_plane__rfm__SinhCylindrical(commondata_struct *restrict commondata, const params_struct *restrict params,
                                                           REAL *restrict xx[3], MoL_gridfunctions_struct *restrict gridfuncs) {
#include "../set_CodeParameters.h"

  // Unpack grid function pointers from gridfuncs struct
  const REAL *restrict y_n_gfs = gridfuncs->y_n_gfs;
  MAYBE_UNUSED const REAL *restrict auxevol_gfs = gridfuncs->auxevol_gfs;
  const REAL *restrict diagnostic_output_gfs = gridfuncs->diagnostic_output_gfs;

  // Prepare output filename based on 2D plane
  char filename[256];
  sprintf(filename, "out2d-yz-%s-conv_factor%.2f-t%08.2f.txt", CoordSystemName, convergence_factor, time);
  FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
  if (!outfile) {
    fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
    exit(1);
  }

  // Define points for output along the yz-plane in SinhCylindrical coordinates.
  const int numpts_i0 = Nxx0, numpts_i1 = 2, numpts_i2 = Nxx2;
  int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];
#pragma omp parallel for
  for (int i0 = NGHOSTS; i0 < Nxx0 + NGHOSTS; i0++)
    i0_pts[i0 - NGHOSTS] = i0;
  i1_pts[0] = (int)(NGHOSTS + (1.0 / 4.0) * Nxx1 - 1.0 / 2.0);
  i1_pts[1] = (int)(NGHOSTS + (3.0 / 4.0) * Nxx1 - 1.0 / 2.0);
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx2 + NGHOSTS; i2++)
    i2_pts[i2 - NGHOSTS] = i2;
  // Main loop to store data points in the specified plane
  LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
    const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
    const int idx3 = IDX3(i0, i1, i2);
    REAL xCart[3];
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
    {
      // Collect diagnostic data
      const REAL log10HL = log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16));
      const REAL alphaL = y_n_gfs[IDX4pt(ALPHAGF, idx3)];
      const REAL hxx = y_n_gfs[IDX4pt(HDD00GF, idx3)];
      const REAL hxy = y_n_gfs[IDX4pt(HDD01GF, idx3)];
      const REAL hxz = y_n_gfs[IDX4pt(HDD02GF, idx3)];
      const REAL hyy = y_n_gfs[IDX4pt(HDD11GF, idx3)];
      const REAL hyz = y_n_gfs[IDX4pt(HDD12GF, idx3)];
      const REAL hzz = y_n_gfs[IDX4pt(HDD22GF, idx3)];
      const REAL T4UU00 = auxevol_gfs[IDX4pt(T4UU00GF, idx3)];
      fprintf(outfile, "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", xCart[1], xCart[2], log10HL, alphaL, hxx, hxy, hxz, hyy,
              hyz, hzz, T4UU00);
    }
  }

  fclose(outfile);
}
