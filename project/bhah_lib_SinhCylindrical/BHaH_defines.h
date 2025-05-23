#ifndef __BHAH_DEFINES_H__
#define __BHAH_DEFINES_H__
// BHaH core header file, automatically generated from output_BHaH_defines_h within BHaH_defines_h.py,
//    DO NOT EDIT THIS FILE BY HAND.

//********************************************
// Basic definitions for module general:
#include <ctype.h>   // Character type functions, such as isdigit, isalpha, etc.
#include <errno.h>   // Error number definitions
#include <math.h>    // Transcendental functions, etc.
#include <stdbool.h> // bool-typed variables
#include <stdint.h>  // int8_t-typed variables
#include <stdio.h>   // Basic input/output functions, such as *printf, fopen, fwrite, etc.
#include <stdlib.h>  // malloc/free, etc.
#include <string.h>  // String handling functions, such as strlen, strcmp, etc.
#include <time.h>    // Time-related functions and types, such as time(), clock(),
// output_BHaH_defines_h(...,enable_intrinsics=True) was called so we intrinsics headers:
#include "intrinsics/simd_intrinsics.h"
#define REAL double
#define DOUBLE double

// These macros for MIN(), MAX(), and SQR() ensure that if the arguments inside
//   are a function/complex expression, the function/expression is evaluated
//   *only once* per argument. See https://lwn.net/Articles/983965/ for details.
// They are improvements over the original implementations:
// #define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
// #define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
// #define SQR(A) ((A) * (A))
#ifndef MIN
#define MIN(A, B)           \
  ({                        \
    __typeof__(A) _a = (A); \
    __typeof__(B) _b = (B); \
    _a < _b ? _a : _b;      \
  })
#endif
#ifndef MAX
#define MAX(A, B)           \
  ({                        \
    __typeof__(A) _a = (A); \
    __typeof__(B) _b = (B); \
    _a > _b ? _a : _b;      \
  })
#endif
#define SQR(A)              \
  ({                        \
    __typeof__(A) _a = (A); \
    _a *_a;                 \
  })
#ifndef MAYBE_UNUSED
#if __cplusplus >= 201703L
#define MAYBE_UNUSED [[maybe_unused]]
#elif defined(__GNUC__) || defined(__clang__) || defined(__NVCC__)
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif // END check for GCC, Clang, or NVCC
#endif // END MAYBE_UNUSED
// START: CodeParameters declared as #define.
#ifndef MAXNUMGRIDS
#define MAXNUMGRIDS 15 // nrpy.grid
#endif
// END: CodeParameters declared as #define.

//********************************************
// Basic definitions for module nrpy.infrastructures.BHaH.diagnostics.progress_indicator:
#ifdef __linux__
// Timer with nanosecond resolution. Only on Linux.
#define TIMEVAR struct timespec
#define CURRTIME_FUNC(currtime) clock_gettime(CLOCK_REALTIME, currtime)
#define TIME_IN_NS(start, end) (REAL)(1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec);
#else
// Low-resolution timer, 1-second resolution. Widely available.
#define TIMEVAR time_t
#define CURRTIME_FUNC(currtime) time(currtime)
#define TIME_IN_NS(start, end) (REAL)(difftime(end, start) * 1.0e9 + 1e-6) // Round up to avoid divide-by-zero.
#endif

//********************************************
// Basic definitions for module commondata_struct:
typedef struct __commondata_struct__ {
  REAL CFL_FACTOR;                    // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL)
  REAL checkpoint_every;              // (nrpy.infrastructures.BHaH.checkpointing)
  REAL convergence_factor;            // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL diagnostics_output_every;      // (nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library)
  REAL dt;                            // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL)
  REAL eta;                           // (nrpy.equations.general_relativity.BSSN_gauge_RHSs)
  REAL initial_central_density;       // (TOVola)
  REAL initial_ode_step_size;         // (TOVola)
  REAL max_step_size;                 // (TOVola)
  REAL min_step_size;                 // (TOVola)
  REAL ode_error_limit;               // (TOVola)
  REAL poly_eos_Gamma;                // (TOVola)
  REAL poly_eos_K;                    // (TOVola)
  REAL t_0;                           // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL)
  REAL t_final;                       // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL)
  REAL time;                          // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL)
  TIMEVAR start_wallclock_time;       // (nrpy.infrastructures.BHaH.diagnostics.progress_indicator)
  char outer_bc_type[50];             // (nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions)
  int NUMGRIDS;                       // (nrpy.grid)
  int interpolation_stencil_size;     // (TOVola)
  int max_interpolation_stencil_size; // (TOVola)
  int nn;                             // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL)
  int nn_0;                           // (nrpy.infrastructures.BHaH.MoLtimestepping.MoL)
  int ode_max_steps;                  // (TOVola)
  int output_progress_every;          // (nrpy.infrastructures.BHaH.diagnostics.progress_indicator)
} commondata_struct;

//********************************************
// Basic definitions for module params_struct:
typedef struct __params_struct__ {
  REAL AMPLRHO;             // (nrpy.reference_metric_SinhCylindrical)
  REAL AMPLZ;               // (nrpy.reference_metric_SinhCylindrical)
  REAL Cart_originx;        // (nrpy.grid)
  REAL Cart_originy;        // (nrpy.grid)
  REAL Cart_originz;        // (nrpy.grid)
  REAL PI;                  // (nrpy.reference_metric)
  REAL SINHWRHO;            // (nrpy.reference_metric_SinhCylindrical)
  REAL SINHWZ;              // (nrpy.reference_metric_SinhCylindrical)
  REAL dxx0;                // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL dxx1;                // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL dxx2;                // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL f0_of_xx0;           // (nrpy.reference_metric_SinhCylindrical)
  REAL f0_of_xx0__D0;       // (nrpy.reference_metric_SinhCylindrical)
  REAL f0_of_xx0__DD00;     // (nrpy.reference_metric_SinhCylindrical)
  REAL f0_of_xx0__DDD000;   // (nrpy.reference_metric_SinhCylindrical)
  REAL f1_of_xx1;           // (nrpy.reference_metric_SinhCylindrical)
  REAL f1_of_xx1__D1;       // (nrpy.reference_metric_SinhCylindrical)
  REAL f1_of_xx1__DD11;     // (nrpy.reference_metric_SinhCylindrical)
  REAL f1_of_xx1__DDD111;   // (nrpy.reference_metric_SinhCylindrical)
  REAL f2_of_xx0;           // (nrpy.reference_metric_SinhCylindrical)
  REAL f2_of_xx0__D0;       // (nrpy.reference_metric_SinhCylindrical)
  REAL f2_of_xx0__DD00;     // (nrpy.reference_metric_SinhCylindrical)
  REAL f3_of_xx2;           // (nrpy.reference_metric_SinhCylindrical)
  REAL f3_of_xx2__D2;       // (nrpy.reference_metric_SinhCylindrical)
  REAL f3_of_xx2__DD22;     // (nrpy.reference_metric_SinhCylindrical)
  REAL f4_of_xx1;           // (nrpy.reference_metric_SinhCylindrical)
  REAL f4_of_xx1__D1;       // (nrpy.reference_metric_SinhCylindrical)
  REAL f4_of_xx1__DD11;     // (nrpy.reference_metric_SinhCylindrical)
  REAL f4_of_xx1__DDD111;   // (nrpy.reference_metric_SinhCylindrical)
  REAL grid_hole_radius;    // (nrpy.reference_metric)
  REAL grid_physical_size;  // (nrpy.reference_metric)
  REAL invdxx0;             // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL invdxx1;             // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL invdxx2;             // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL xxmax0;              // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL xxmax1;              // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL xxmax2;              // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL xxmin0;              // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL xxmin1;              // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  REAL xxmin2;              // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  bool grid_rotates;        // (nrpy.grid)
  char CoordSystemName[50]; // (nrpy.reference_metric)
  int CoordSystem_hash;     // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  int Nxx0;                 // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  int Nxx1;                 // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  int Nxx2;                 // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  int Nxx_plus_2NGHOSTS0;   // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  int Nxx_plus_2NGHOSTS1;   // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
  int Nxx_plus_2NGHOSTS2;   // (nrpy.infrastructures.BHaH.numerical_grids_and_timestep)
} params_struct;

//********************************************
// Basic definitions for module finite_difference:

// Set the number of ghost zones
// Note that upwinding in e.g., BSSN requires that NGHOSTS = fd_order/2 + 1 <- Notice the +1.
#define NGHOSTS 3

// Declare NO_INLINE macro, used in FD functions. GCC v10+ compilations hang on complex RHS expressions (like BSSN) without this.
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define NO_INLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define NO_INLINE __declspec(noinline)
#else
#define NO_INLINE // Fallback for unknown compilers
#endif

//********************************************
// Basic definitions for module reference_metric:
typedef struct __rfmstruct__ {
  REAL *restrict f0_of_xx0;
  REAL *restrict f0_of_xx0__D0;
  REAL *restrict f0_of_xx0__DD00;
  REAL *restrict f0_of_xx0__DDD000;
  REAL *restrict f3_of_xx2;
  REAL *restrict f3_of_xx2__D2;
  REAL *restrict f3_of_xx2__DD22;
} rfm_struct;

//********************************************
// Basic definitions for module nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions:

// NRPy+ Curvilinear Boundary Conditions: Core data structures
// Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

typedef struct __innerpt_bc_struct__ {
  int dstpt;         // dstpt is the 3D grid index IDX3(i0,i1,i2) of the inner boundary point (i0,i1,i2)
  int srcpt;         // srcpt is the 3D grid index (a la IDX3) to which the inner boundary point maps
  int8_t parity[10]; // parity[10] is a calculation of dot products for the 10 independent parity types
} innerpt_bc_struct;

typedef struct __outerpt_bc_struct__ {
  short i0, i1, i2;              // the outer boundary point grid index (i0,i1,i2), on the 3D grid
  int8_t FACEX0, FACEX1, FACEX2; // 1-byte integers that store
  //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
  //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
  //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
  //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
} outerpt_bc_struct;

typedef struct __bc_info_struct__ {
  int num_inner_boundary_points;                  // stores total number of inner boundary points
  int num_pure_outer_boundary_points[NGHOSTS][3]; // stores number of outer boundary points on each
  //                                                  ghostzone level and direction (update min and
  //                                                  max faces simultaneously on multiple cores)
  int bc_loop_bounds[NGHOSTS][6][6]; // stores outer boundary loop bounds. Unused after bcstruct_set_up()
} bc_info_struct;

typedef struct __bc_struct__ {
  innerpt_bc_struct *restrict inner_bc_array;                   // information needed for updating each inner boundary point
  outerpt_bc_struct *restrict pure_outer_bc_array[NGHOSTS * 3]; // information needed for updating each outer
  //                                                             boundary point
  bc_info_struct bc_info; // stores number of inner and outer boundary points, needed for setting loop
  //                          bounds and parallelizing over as many boundary points as possible.
} bc_struct;

/* PARITY TYPES FOR EVOLVED (plus optional) GRIDFUNCTIONS.
 * SEE "Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb" FOR DEFINITIONS. */
static const int8_t evol_gf_parity[24] = {4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 0, 4, 5, 6, 7, 8, 9, 1, 2, 3, 0, 1, 2, 3};

//********************************************
// Basic definitions for module nrpy.infrastructures.BHaH.MoLtimestepping.MoL:
typedef struct __MoL_gridfunctions_struct__ {
  REAL *restrict y_n_gfs;
  REAL *restrict y_nplus1_running_total_gfs;
  REAL *restrict k_odd_gfs;
  REAL *restrict k_even_gfs;
  REAL *restrict auxevol_gfs;
  REAL *restrict diagnostic_output_gfs;
  REAL *restrict diagnostic_output_gfs2;
} MoL_gridfunctions_struct;

#define LOOP_ALL_GFS_GPS(ii) \
  _Pragma("omp parallel for") for (int(ii) = 0; (ii) < Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_EVOL_GFS; (ii)++)

//********************************************
// Basic definitions for module grid:

// EVOL VARIABLES:
#define NUM_EVOL_GFS 24
#define ADD00GF 0
#define ADD01GF 1
#define ADD02GF 2
#define ADD11GF 3
#define ADD12GF 4
#define ADD22GF 5
#define ALPHAGF 6
#define BETU0GF 7
#define BETU1GF 8
#define BETU2GF 9
#define CFGF 10
#define HDD00GF 11
#define HDD01GF 12
#define HDD02GF 13
#define HDD11GF 14
#define HDD12GF 15
#define HDD22GF 16
#define LAMBDAU0GF 17
#define LAMBDAU1GF 18
#define LAMBDAU2GF 19
#define TRKGF 20
#define VETU0GF 21
#define VETU1GF 22
#define VETU2GF 23

// SET gridfunctions_f_infinity[i] = evolved gridfunction i's value in the limit r->infinity:
static const REAL gridfunctions_f_infinity[NUM_EVOL_GFS] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                                                            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

// SET gridfunctions_wavespeed[i] = evolved gridfunction i's characteristic wave speed:
static const REAL gridfunctions_wavespeed[NUM_EVOL_GFS] = {
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.41421356237310, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

// AUX VARIABLES:
#define NUM_AUX_GFS 2
#define HGF 0
#define MSQUAREDGF 1

// AUXEVOL VARIABLES:
#define NUM_AUXEVOL_GFS 16
#define RBARDD00GF 0
#define RBARDD01GF 1
#define RBARDD02GF 2
#define RBARDD11GF 3
#define RBARDD12GF 4
#define RBARDD22GF 5
#define T4UU00GF 6
#define T4UU01GF 7
#define T4UU02GF 8
#define T4UU03GF 9
#define T4UU11GF 10
#define T4UU12GF 11
#define T4UU13GF 12
#define T4UU22GF 13
#define T4UU23GF 14
#define T4UU33GF 15

// Declare the IDX4(gf,i,j,k) macro, which enables us to store 4-dimensions of
//   data in a 1D array. In this case, consecutive values of "i"
//   (all other indices held to a fixed value) are consecutive in memory, where
//   consecutive values of "j" (fixing all other indices) are separated by
//   Nxx_plus_2NGHOSTS0 elements in memory. Similarly, consecutive values of
//   "k" are separated by Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1 in memory, etc.
#define IDX4(gf, i, j, k) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k) + Nxx_plus_2NGHOSTS2 * (gf))))
#define IDX4pt(gf, idx) ((idx) + (Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2) * (gf))
#define IDX3(i, j, k) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k))))
#define LOOP_REGION(i0min, i0max, i1min, i1max, i2min, i2max) \
  for (int i2 = i2min; i2 < i2max; i2++)                      \
    for (int i1 = i1min; i1 < i1max; i1++)                    \
      for (int i0 = i0min; i0 < i0max; i0++)
#define LOOP_OMP(__OMP_PRAGMA__, i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max)                            \
  _Pragma(__OMP_PRAGMA__) for (int(i2) = (i2min); (i2) < (i2max); (i2)++) for (int(i1) = (i1min); (i1) < (i1max); \
                                                                               (i1)++) for (int(i0) = (i0min); (i0) < (i0max); (i0)++)
#define LOOP_NOOMP(i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max) \
  for (int(i2) = (i2min); (i2) < (i2max); (i2)++)                        \
    for (int(i1) = (i1min); (i1) < (i1max); (i1)++)                      \
      for (int(i0) = (i0min); (i0) < (i0max); (i0)++)
#define LOOP_BREAKOUT(i0, i1, i2, i0max, i1max, i2max) \
  {                                                    \
    i0 = (i0max);                                      \
    i1 = (i1max);                                      \
    i2 = (i2max);                                      \
    break;                                             \
  }
#define IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NG)                                \
  (i0i1i2[0] >= (NG) && i0i1i2[0] < (Nxx_plus_2NGHOSTS0) - (NG) && i0i1i2[1] >= (NG) && i0i1i2[1] < (Nxx_plus_2NGHOSTS1) - (NG) && \
   i0i1i2[2] >= (NG) && i0i1i2[2] < (Nxx_plus_2NGHOSTS2) - (NG))

typedef struct __griddata__ {
  // griddata_struct stores data needed on each grid
  // xx[3] stores the uniform grid coordinates.
  REAL *restrict xx[3];
  // NRPy+ MODULE: nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions
  bc_struct bcstruct; // <- all data needed to apply boundary conditions in curvilinear coordinates
  // NRPy+ MODULE: nrpy.infrastructures.BHaH.MoLtimestepping.MoL
  MoL_gridfunctions_struct gridfuncs; // <- MoL gridfunctions
  // NRPy+ MODULE: params
  params_struct params; // <- BHaH parameters, generated from NRPy+'s CodeParameters
  // NRPy+ MODULE: reference_metric
  char CoordSystemname[100]; // <- the name of the CoordSystem (from reference_metric)
  char gridname[100];        // <- a user-defined alias for describing the grid
  rfm_struct rfmstruct;      // <- includes e.g., 1D arrays of reference metric quantities
} griddata_struct;

//********************************************
// Basic definitions for module nrpy.infrastructures.BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter:
typedef struct __initial_data_struct__ {
  REAL alpha;

  REAL betaSphorCartU0, betaSphorCartU1, betaSphorCartU2;
  REAL BSphorCartU0, BSphorCartU1, BSphorCartU2;

  REAL gammaSphorCartDD00, gammaSphorCartDD01, gammaSphorCartDD02;
  REAL gammaSphorCartDD11, gammaSphorCartDD12, gammaSphorCartDD22;

  REAL KSphorCartDD00, KSphorCartDD01, KSphorCartDD02;
  REAL KSphorCartDD11, KSphorCartDD12, KSphorCartDD22;

  REAL T4SphorCartUU00, T4SphorCartUU01, T4SphorCartUU02, T4SphorCartUU03;
  REAL T4SphorCartUU11, T4SphorCartUU12, T4SphorCartUU13;
  REAL T4SphorCartUU22, T4SphorCartUU23;
  REAL T4SphorCartUU33;

} initial_data_struct;
typedef struct __ID_persist_struct__ {

  // The following arrays store stellar information at all numpoints_arr radii:
  REAL *restrict r_Schw_arr;     // Stellar radial coordinate in units of Schwarzschild radius
  REAL *restrict rho_baryon_arr; // Baryonic mass density
  REAL *restrict rho_energy_arr; // Mass-energy density
  REAL *restrict P_arr;          // Pressure
  REAL *restrict M_arr;          // Integrated rest mass
  REAL *restrict expnu_arr;      // Metric quantity
  REAL *restrict exp4phi_arr;    // Metric quantity
  REAL *restrict r_iso_arr;      // Isotropic radial coordinate, in literature, sometimes called rbar
  int numpoints_arr;             // Number of radii stored in the arrays

} ID_persist_struct;

//********************************************
// Basic definitions for module nrpy.infrastructures.BHaH.rfm_wrapper_functions:
#define SINHCYLINDRICAL 1325205359

//********************************************
// Basic definitions for module BHaH Lib:

typedef struct BHaH_struct {
  commondata_struct *commondata;
  griddata_struct *griddata;
} BHaH_struct;
#endif
