#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
static inline REAL FD1_arbitrary_upwind_x0_dirn(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                const REAL *restrict gf, const int i0, const int i1, const int i2, const int offset) {
#include "../set_CodeParameters.h"
  switch (offset) {
  case 0:
    return (+1.0 / 12.0 * gf[IDX3(i0 - 2, i1, i2)] - 2.0 / 3.0 * gf[IDX3(i0 - 1, i1, i2)] + 2.0 / 3.0 * gf[IDX3(i0 + 1, i1, i2)] -
            1.0 / 12.0 * gf[IDX3(i0 + 2, i1, i2)]) *
           invdxx0;
  case 1:
    return (-1.0 / 4.0 * gf[IDX3(i0 - 1, i1, i2)] - 5.0 / 6.0 * gf[IDX3(i0, i1, i2)] + 3.0 / 2.0 * gf[IDX3(i0 + 1, i1, i2)] -
            1.0 / 2.0 * gf[IDX3(i0 + 2, i1, i2)] + 1.0 / 12.0 * gf[IDX3(i0 + 3, i1, i2)]) *
           invdxx0;
  case -1:
    return (-1.0 / 12.0 * gf[IDX3(i0 - 3, i1, i2)] + 1.0 / 2.0 * gf[IDX3(i0 - 2, i1, i2)] - 3.0 / 2.0 * gf[IDX3(i0 - 1, i1, i2)] +
            5.0 / 6.0 * gf[IDX3(i0, i1, i2)] + 1.0 / 4.0 * gf[IDX3(i0 + 1, i1, i2)]) *
           invdxx0;
  case 2:
    return (-25.0 / 12.0 * gf[IDX3(i0, i1, i2)] + 4 * gf[IDX3(i0 + 1, i1, i2)] - 3 * gf[IDX3(i0 + 2, i1, i2)] + 4.0 / 3.0 * gf[IDX3(i0 + 3, i1, i2)] -
            1.0 / 4.0 * gf[IDX3(i0 + 4, i1, i2)]) *
           invdxx0;
  case -2:
    return (+1.0 / 4.0 * gf[IDX3(i0 - 4, i1, i2)] - 4.0 / 3.0 * gf[IDX3(i0 - 3, i1, i2)] + 3 * gf[IDX3(i0 - 2, i1, i2)] -
            4 * gf[IDX3(i0 - 1, i1, i2)] + 25.0 / 12.0 * gf[IDX3(i0, i1, i2)]) *
           invdxx0;
  }
  return 0.0 / 0.0; // poison output if offset computed incorrectly
}
/**
 * Compute 1st derivative finite-difference derivative with arbitrary upwind
 */
static inline REAL FD1_arbitrary_upwind_x2_dirn(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                const REAL *restrict gf, const int i0, const int i1, const int i2, const int offset) {
#include "../set_CodeParameters.h"
  switch (offset) {
  case 0:
    return (+1.0 / 12.0 * gf[IDX3(i0, i1, i2 - 2)] - 2.0 / 3.0 * gf[IDX3(i0, i1, i2 - 1)] + 2.0 / 3.0 * gf[IDX3(i0, i1, i2 + 1)] -
            1.0 / 12.0 * gf[IDX3(i0, i1, i2 + 2)]) *
           invdxx2;
  case 1:
    return (-1.0 / 4.0 * gf[IDX3(i0, i1, i2 - 1)] - 5.0 / 6.0 * gf[IDX3(i0, i1, i2)] + 3.0 / 2.0 * gf[IDX3(i0, i1, i2 + 1)] -
            1.0 / 2.0 * gf[IDX3(i0, i1, i2 + 2)] + 1.0 / 12.0 * gf[IDX3(i0, i1, i2 + 3)]) *
           invdxx2;
  case -1:
    return (-1.0 / 12.0 * gf[IDX3(i0, i1, i2 - 3)] + 1.0 / 2.0 * gf[IDX3(i0, i1, i2 - 2)] - 3.0 / 2.0 * gf[IDX3(i0, i1, i2 - 1)] +
            5.0 / 6.0 * gf[IDX3(i0, i1, i2)] + 1.0 / 4.0 * gf[IDX3(i0, i1, i2 + 1)]) *
           invdxx2;
  case 2:
    return (-25.0 / 12.0 * gf[IDX3(i0, i1, i2)] + 4 * gf[IDX3(i0, i1, i2 + 1)] - 3 * gf[IDX3(i0, i1, i2 + 2)] + 4.0 / 3.0 * gf[IDX3(i0, i1, i2 + 3)] -
            1.0 / 4.0 * gf[IDX3(i0, i1, i2 + 4)]) *
           invdxx2;
  case -2:
    return (+1.0 / 4.0 * gf[IDX3(i0, i1, i2 - 4)] - 4.0 / 3.0 * gf[IDX3(i0, i1, i2 - 3)] + 3 * gf[IDX3(i0, i1, i2 - 2)] -
            4 * gf[IDX3(i0, i1, i2 - 1)] + 25.0 / 12.0 * gf[IDX3(i0, i1, i2)]) *
           invdxx2;
  }
  return 0.0 / 0.0; // poison output if offset computed incorrectly
}
/**
 * Compute r(xx0,xx1,xx2) and partial_r x^i.
 */
static inline void r_and_partial_xi_partial_r_derivs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                     const REAL xx0, const REAL xx1, const REAL xx2, REAL *r, REAL *partial_x0_partial_r,
                                                     REAL *partial_x1_partial_r, REAL *partial_x2_partial_r) {
#include "../set_CodeParameters.h"
  const REAL tmp0 = (1.0 / (SINHWRHO));
  const REAL tmp6 = (1.0 / (SINHWZ));
  const REAL tmp5 = ((AMPLRHO) * (AMPLRHO)) / ((exp(tmp0) - exp(-tmp0)) * (exp(tmp0) - exp(-tmp0)));
  const REAL tmp7 = exp(tmp6) - exp(-tmp6);
  const REAL tmp2 = exp(tmp0 * xx0);
  const REAL tmp3 = exp(-tmp0 * xx0);
  const REAL tmp9 = exp(tmp6 * xx2);
  const REAL tmp10 = exp(-tmp6 * xx2);
  const REAL tmp24 = (1.0 / (tmp7));
  const REAL tmp4 = tmp2 - tmp3;
  const REAL tmp11 = -tmp10 + tmp9;
  const REAL tmp21 = tmp4 * tmp5 * (2 * tmp0 * tmp2 + 2 * tmp0 * tmp3);
  const REAL tmp13 = ((AMPLZ) * (AMPLZ)) * ((tmp11) * (tmp11)) / ((tmp7) * (tmp7));
  const REAL tmp18 = ((AMPLZ) * (AMPLZ) * (AMPLZ)) * ((tmp11) * (tmp11)) * (2 * tmp10 * tmp6 + 2 * tmp6 * tmp9) / ((tmp7) * (tmp7) * (tmp7));
  const REAL tmp14 = tmp13 + ((tmp4) * (tmp4)) * tmp5;
  const REAL tmp15 = sqrt(tmp14);
  const REAL tmp19 = (1.0 / sqrt(-tmp13 / tmp14 + 1));
  const REAL tmp25 = (1.0 / 2.0) / pow(tmp14, 3.0 / 2.0);
  const REAL tmp23 = (1.0 / (tmp15));
  const REAL tmp26 = tmp19 * (AMPLZ * tmp23 * tmp24 * (tmp10 * tmp6 + tmp6 * tmp9) - tmp18 * tmp25);
  const REAL tmp27 = (1.0 / ((1.0 / 2.0) * tmp21 * tmp23 * tmp26 + (1.0 / 4.0) * tmp18 * tmp19 * tmp21 / ((tmp14) * (tmp14))));
  *r = tmp15;
  *partial_x0_partial_r = tmp26 * tmp27;
  *partial_x1_partial_r = 0;
  *partial_x2_partial_r = AMPLZ * tmp11 * tmp19 * tmp21 * tmp24 * tmp25 * tmp27;
}
/**
 * Compute \partial_r f
 */
static inline REAL compute_partial_r_f(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                       const REAL *restrict gfs, const int which_gf, const int dest_i0, const int dest_i1, const int dest_i2,
                                       const int FACEi0, const int FACEi1, const int FACEi2, const REAL partial_x0_partial_r,
                                       const REAL partial_x1_partial_r, const REAL partial_x2_partial_r) {
#include "../set_CodeParameters.h"
  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_fd_order/2 = 2
  const int FD1_stencil_radius = 2;

  const int ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  ///////////////////////////////////////////////////////////
  // Next we'll compute partial_xi f, using a maximally-centered stencil.
  //   The {i0,i1,i2}_offset parameters set the offset of the maximally-centered
  //   stencil, such that an offset=0 implies a centered stencil.

  // CHECK: Nxx_plus_2NGHOSTS0=10; FD1_stencil_radius=2. Then Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1 = 7
  //  if dest_i0 = 9, we get i0_offset=7-9=-2, so the (4th order) deriv
  //  stencil is: -4,-3,-2,-1,0

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 1, we get i0_offset = FD1_stencil_radius-dest_i0 = 1,
  //  so the (4th order) deriv stencil is: -1,0,1,2,3

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 0, we get i0_offset = FD1_stencil_radius-1 = 2,
  //  so the (4th order) deriv stencil is: 0,1,2,3,4
  int i0_offset = FACEi0; // Shift stencil away from the face we're updating.
  // Next adjust i0_offset so that FD stencil never goes out of bounds.
  if (dest_i0 < FD1_stencil_radius)
    i0_offset = FD1_stencil_radius - dest_i0;
  else if (dest_i0 > (Nxx_plus_2NGHOSTS0 - FD1_stencil_radius - 1))
    i0_offset = (Nxx_plus_2NGHOSTS0 - FD1_stencil_radius - 1) - dest_i0;
  const REAL partial_x0_f = FD1_arbitrary_upwind_x0_dirn(commondata, params, &gfs[which_gf * ntot], dest_i0, dest_i1, dest_i2, i0_offset);
  const REAL partial_x1_f = 0.0;
  int i2_offset = FACEi2; // Shift stencil away from the face we're updating.
  // Next adjust i2_offset so that FD stencil never goes out of bounds.
  if (dest_i2 < FD1_stencil_radius)
    i2_offset = FD1_stencil_radius - dest_i2;
  else if (dest_i2 > (Nxx_plus_2NGHOSTS2 - FD1_stencil_radius - 1))
    i2_offset = (Nxx_plus_2NGHOSTS2 - FD1_stencil_radius - 1) - dest_i2;
  const REAL partial_x2_f = FD1_arbitrary_upwind_x2_dirn(commondata, params, &gfs[which_gf * ntot], dest_i0, dest_i1, dest_i2, i2_offset);
  return partial_x0_partial_r * partial_x0_f + partial_x1_partial_r * partial_x1_f + partial_x2_partial_r * partial_x2_f;
}

/**
 * *** Apply radiation BCs to all outer boundaries. ***
 *
 */
static inline REAL radiation_bcs(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                 const bc_struct *restrict bcstruct, REAL *restrict xx[3], const REAL *restrict gfs, REAL *restrict gfs_rhss,
                                 const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity, const int dest_i0, const int dest_i1,
                                 const int dest_i2, const short FACEi0, const short FACEi1, const short FACEi2) {
#include "../set_CodeParameters.h"
  // Nearest "interior" neighbor of this gridpoint, based on current face
  const int dest_i0_int = dest_i0 + 1 * FACEi0, dest_i1_int = dest_i1 + 1 * FACEi1, dest_i2_int = dest_i2 + 1 * FACEi2;
  REAL r, partial_x0_partial_r, partial_x1_partial_r, partial_x2_partial_r;
  REAL r_int, partial_x0_partial_r_int, partial_x1_partial_r_int, partial_x2_partial_r_int;
  r_and_partial_xi_partial_r_derivs(commondata, params, xx[0][dest_i0], xx[1][dest_i1], xx[2][dest_i2], &r, &partial_x0_partial_r,
                                    &partial_x1_partial_r, &partial_x2_partial_r);
  r_and_partial_xi_partial_r_derivs(commondata, params, xx[0][dest_i0_int], xx[1][dest_i1_int], xx[2][dest_i2_int], &r_int, &partial_x0_partial_r_int,
                                    &partial_x1_partial_r_int, &partial_x2_partial_r_int);
  const REAL partial_r_f = compute_partial_r_f(commondata, params, xx, gfs, which_gf, dest_i0, dest_i1, dest_i2, FACEi0, FACEi1, FACEi2,
                                               partial_x0_partial_r, partial_x1_partial_r, partial_x2_partial_r);
  const REAL partial_r_f_int = compute_partial_r_f(commondata, params, xx, gfs, which_gf, dest_i0_int, dest_i1_int, dest_i2_int, FACEi0, FACEi1,
                                                   FACEi2, partial_x0_partial_r_int, partial_x1_partial_r_int, partial_x2_partial_r_int);

  const int idx3 = IDX3(dest_i0, dest_i1, dest_i2);
  const int idx3_int = IDX3(dest_i0_int, dest_i1_int, dest_i2_int);

  const REAL partial_t_f_int = gfs_rhss[IDX4pt(which_gf, idx3_int)];

  const REAL c = gf_wavespeed;
  const REAL f_infinity = gf_f_infinity;
  const REAL f = gfs[IDX4pt(which_gf, idx3)];
  const REAL f_int = gfs[IDX4pt(which_gf, idx3_int)];
  const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) / r_int);

  const REAL k = r_int * r_int * r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

  const REAL rinv = 1.0 / r;
  const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

  return partial_t_f_outgoing_wave + k * rinv * rinv * rinv;
}

/**
 * This function is responsible for applying boundary conditions (BCs) to both pure outer and inner
 * boundary points. In the first step, it parallelizes the task using OpenMP and starts by applying BCs to
 * the outer boundary points layer-by-layer, prioritizing the faces in the order x0, x1, x2. The second step
 * applies BCs to the inner boundary points, which may map either to the grid interior or to the outer boundary.
 *
 */
void apply_bcs_outerradiation_and_inner__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                              const bc_struct *restrict bcstruct, REAL *restrict xx[3],
                                                              const REAL custom_wavespeed[NUM_EVOL_GFS], const REAL custom_f_infinity[NUM_EVOL_GFS],
                                                              REAL *restrict gfs, REAL *restrict rhs_gfs) {
#include "../set_CodeParameters.h"

  // Unpack bc_info from bcstruct
  const bc_info_struct *bc_info = &bcstruct->bc_info;

  ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  // Spawn N OpenMP threads, either across all cores, or according to e.g., taskset.
#pragma omp parallel
  {
    for (int which_gz = 0; which_gz < NGHOSTS; which_gz++)
      for (int dirn = 0; dirn < 3; dirn++) {
        // This option results in about 1.6% slower runtime for SW curvilinear at 64x24x24 on 8-core Ryzen 9 4900HS
        // #pragma omp for collapse(2)
        // for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) for(int idx2d=0;idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];idx2d++)
        // {
        //  {
        // Don't spawn a thread if there are no boundary points to fill; results in a nice little speedup.
        if (bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for // threads have been spawned; here we distribute across them
          for (int idx2d = 0; idx2d < bc_info->num_pure_outer_boundary_points[which_gz][dirn]; idx2d++) {
            const short i0 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0;
            const short i1 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1;
            const short i2 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2;
            const short FACEX0 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0;
            const short FACEX1 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1;
            const short FACEX2 = bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2;
            const int idx3 = IDX3(i0, i1, i2);
            for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
              // *** Apply radiation BCs to all outer boundary points. ***
              rhs_gfs[IDX4pt(which_gf, idx3)] = radiation_bcs(commondata, params, bcstruct, xx, gfs, rhs_gfs, which_gf, custom_wavespeed[which_gf],
                                                              custom_f_infinity[which_gf], i0, i1, i2, FACEX0, FACEX1, FACEX2);
            }
          }
        }
      }
  }

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(commondata, params, bcstruct, rhs_gfs); // <- apply inner BCs to RHS gfs only
}
