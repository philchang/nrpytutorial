CC ?= gcc  # assigns the value CC to gcc only if environment variable CC is not already set

CFLAGS = -std=gnu99 -O2 -march=native -g -Wall $(shell gsl-config --cflags) -fPIC -shared
VALGRIND_CFLAGS = -std=gnu99 -O2 -g -Wall -Wno-unknown-pragmas $(shell gsl-config --cflags) -fPIC -shared
INCLUDEDIRS =
LDFLAGS = $(shell gsl-config --libs) -lm

# Check for OpenMP support
OPENMP_FLAG = -fopenmp
COMPILER_SUPPORTS_OPENMP := $(shell echo | $(CC) $(OPENMP_FLAG) -E - >/dev/null 2>&1 && echo YES || echo NO)

ifeq ($(COMPILER_SUPPORTS_OPENMP), YES)
    CFLAGS += $(OPENMP_FLAG)
    LDFLAGS += $(OPENMP_FLAG)
endif

OBJ_FILES = apply_bcs_inner_only.o apply_bcs_outerextrap_and_inner.o apply_bcs_outerradiation_and_inner.o bcstruct_set_up.o bhah_diagnostics.o bhah_evolve.o bhah_finalize.o bhah_initialize.o Cart_to_xx_and_nearest_i0i1i2.o cfl_limited_timestep.o commondata_struct_set_to_default.o constraints_eval.o diagnostics.o diagnostics_nearest_1d_y_axis.o diagnostics_nearest_1d_z_axis.o diagnostics_nearest_2d_xy_plane.o diagnostics_nearest_2d_yz_plane.o diagnostics_nearest_grid_center.o enforce_detgammabar_equals_detgammahat.o initial_data.o initial_data_reader__convert_ADM_Spherical_to_BSSN.o MoL_free_memory_non_y_n_gfs.o MoL_free_memory_y_n_gfs.o MoL_malloc_non_y_n_gfs.o MoL_malloc_y_n_gfs.o MoL_step_forward_in_time.o numerical_grid_params_Nxx_dxx_xx.o numerical_grids_and_timestep.o params_struct_set_to_default.o progress_indicator.o read_checkpoint.o rfm_precompute_defines.o rfm_precompute_free.o rfm_precompute_malloc.o rhs_eval.o Ricci_eval.o SinhCylindrical/apply_bcs_outerradiation_and_inner__rfm__SinhCylindrical.o SinhCylindrical/bcstruct_set_up__rfm__SinhCylindrical.o SinhCylindrical/Cart_to_xx_and_nearest_i0i1i2__rfm__SinhCylindrical.o SinhCylindrical/cfl_limited_timestep__rfm__SinhCylindrical.o SinhCylindrical/constraints_eval__rfm__SinhCylindrical.o SinhCylindrical/diagnostics_nearest_1d_y_axis__rfm__SinhCylindrical.o SinhCylindrical/diagnostics_nearest_1d_z_axis__rfm__SinhCylindrical.o SinhCylindrical/diagnostics_nearest_2d_xy_plane__rfm__SinhCylindrical.o SinhCylindrical/diagnostics_nearest_2d_yz_plane__rfm__SinhCylindrical.o SinhCylindrical/diagnostics_nearest_grid_center__rfm__SinhCylindrical.o SinhCylindrical/enforce_detgammabar_equals_detgammahat__rfm__SinhCylindrical.o SinhCylindrical/initial_data_reader__convert_ADM_Spherical_to_BSSN__rfm__SinhCylindrical.o SinhCylindrical/numerical_grid_params_Nxx_dxx_xx__rfm__SinhCylindrical.o SinhCylindrical/rfm_precompute_defines__rfm__SinhCylindrical.o SinhCylindrical/rfm_precompute_free__rfm__SinhCylindrical.o SinhCylindrical/rfm_precompute_malloc__rfm__SinhCylindrical.o SinhCylindrical/rhs_eval__rfm__SinhCylindrical.o SinhCylindrical/Ricci_eval__rfm__SinhCylindrical.o SinhCylindrical/xx_to_Cart__rfm__SinhCylindrical.o TOVola/TOVola_interp.o TOVola/TOVola_solve.o write_checkpoint.o xx_to_Cart.o

all: libbhah_lib.so

%.o: %.c $(COMMON_HEADERS)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

libbhah_lib.so: $(OBJ_FILES)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)


valgrind: clean
	$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all
# Use $(RM) to be cross-platform compatible.
clean:
	$(RM) *.o */*.o *~ */*~ ./#* *.txt *.gp *.dat *.avi *.png libbhah_lib.so
