# CC ?= gcc  # assigns the value CC to gcc only if environment variable CC is not already set

# CFLAGS = -std=gnu99 -O2 -march=native -g -Wall -Wno-unused-variable $(shell gsl-config --cflags) -fPIC -shared
# #CFLAGS = -std=gnu99 -O3 -funroll-loops -march=native -g -Wall -Wno-unused-variable $(shell gsl-config --cflags) -fPIC -shared
# #CFLAGS = -std=gnu99 -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas $(shell gsl-config --cflags) -fPIC -shared

# VALGRIND_CFLAGS = -std=gnu99 -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas $(shell gsl-config --cflags) -fPIC -shared
# INCLUDEDIRS =
# LDFLAGS = $(shell gsl-config --libs) -lm

# # Check for OpenMP support
# OPENMP_FLAG = -fopenmp
# COMPILER_SUPPORTS_OPENMP := $(shell echo | $(CC) $(OPENMP_FLAG) -E - >/dev/null 2>&1 && echo YES || echo NO)

# ifeq ($(COMPILER_SUPPORTS_OPENMP), YES)
#     CFLAGS += $(OPENMP_FLAG)
#     LDFLAGS += $(OPENMP_FLAG)
# endif

# OBJ_FILES = Cart_to_xx_and_nearest_i0i1i2.o MoL_free_memory_non_y_n_gfs.o MoL_free_memory_y_n_gfs.o MoL_malloc_non_y_n_gfs.o MoL_malloc_y_n_gfs.o MoL_step_forward_in_time.o Ricci_eval.o Spherical/Cart_to_xx_and_nearest_i0i1i2__rfm__Spherical.o Spherical/Ricci_eval__rfm__Spherical.o Spherical/apply_bcs_outerradiation_and_inner__rfm__Spherical.o Spherical/bcstruct_set_up__rfm__Spherical.o Spherical/cfl_limited_timestep__rfm__Spherical.o Spherical/constraints_eval__rfm__Spherical.o Spherical/diagnostics_nearest_1d_y_axis__rfm__Spherical.o Spherical/diagnostics_nearest_1d_z_axis__rfm__Spherical.o Spherical/diagnostics_nearest_2d_xy_plane__rfm__Spherical.o Spherical/diagnostics_nearest_2d_yz_plane__rfm__Spherical.o Spherical/diagnostics_nearest_grid_center__rfm__Spherical.o Spherical/enforce_detgammabar_equals_detgammahat__rfm__Spherical.o Spherical/initial_data_reader__convert_ADM_Spherical_to_BSSN__rfm__Spherical.o Spherical/numerical_grid_params_Nxx_dxx_xx__rfm__Spherical.o Spherical/rfm_precompute_defines__rfm__Spherical.o Spherical/rfm_precompute_free__rfm__Spherical.o Spherical/rfm_precompute_malloc__rfm__Spherical.o Spherical/rhs_eval__rfm__Spherical.o Spherical/xx_to_Cart__rfm__Spherical.o TOVola/TOVola_interp.o TOVola/TOVola_solve.o apply_bcs_inner_only.o apply_bcs_outerextrap_and_inner.o apply_bcs_outerradiation_and_inner.o bcstruct_set_up.o bhah_diagnostics.o bhah_evolve.o bhah_finalize.o bhah_initialize.o cfl_limited_timestep.o commondata_struct_set_to_default.o constraints_eval.o diagnostics.o diagnostics_nearest_1d_y_axis.o diagnostics_nearest_1d_z_axis.o diagnostics_nearest_2d_xy_plane.o diagnostics_nearest_2d_yz_plane.o diagnostics_nearest_grid_center.o enforce_detgammabar_equals_detgammahat.o initial_data.o initial_data_reader__convert_ADM_Spherical_to_BSSN.o numerical_grid_params_Nxx_dxx_xx.o numerical_grids_and_timestep.o params_struct_set_to_default.o progress_indicator.o read_checkpoint.o rfm_precompute_defines.o rfm_precompute_free.o rfm_precompute_malloc.o rhs_eval.o write_checkpoint.o xx_to_Cart.o

# all: libbhah_lib.so

# %.o: %.c $(COMMON_HEADERS)
# 	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

# libbhah_lib.so: $(OBJ_FILES)
# 	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# valgrind: clean
# 	$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all

# # Use $(RM) to be cross-platform compatible.
# clean:
# 	$(RM) *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png libbhah_lib.so


CC ?= gcc  # assigns the value CC to gcc only if environment variable CC is not already set

# Remove -shared to create an object suitable for static libraries.
CFLAGS = -std=gnu99 -O2 -march=native -g -Wall -Wno-unused-variable $(shell gsl-config --cflags) -fPIC
#CFLAGS = -std=gnu99 -O3 -funroll-loops -march=native -g -Wall -Wno-unused-variable $(shell gsl-config --cflags) -fPIC
#CFLAGS = -std=gnu99 -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas $(shell gsl-config --cflags) -fPIC

VALGRIND_CFLAGS = -std=gnu99 -O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas $(shell gsl-config --cflags) -fPIC
INCLUDEDIRS =
LDFLAGS = $(shell gsl-config --libs) -lm

# Check for OpenMP support
OPENMP_FLAG = -fopenmp
COMPILER_SUPPORTS_OPENMP := $(shell echo | $(CC) $(OPENMP_FLAG) -E - >/dev/null 2>&1 && echo YES || echo NO)

ifeq ($(COMPILER_SUPPORTS_OPENMP), YES)
    CFLAGS += $(OPENMP_FLAG)
    LDFLAGS += $(OPENMP_FLAG)
endif

OBJ_FILES = bhah_lib.o Cart_to_xx_and_nearest_i0i1i2.o MoL_free_memory_non_y_n_gfs.o MoL_free_memory_y_n_gfs.o MoL_malloc_non_y_n_gfs.o MoL_malloc_y_n_gfs.o MoL_step_forward_in_time.o Ricci_eval.o Spherical/Cart_to_xx_and_nearest_i0i1i2__rfm__Spherical.o Spherical/Ricci_eval__rfm__Spherical.o Spherical/apply_bcs_outerradiation_and_inner__rfm__Spherical.o Spherical/bcstruct_set_up__rfm__Spherical.o Spherical/cfl_limited_timestep__rfm__Spherical.o Spherical/constraints_eval__rfm__Spherical.o Spherical/diagnostics_nearest_1d_y_axis__rfm__Spherical.o Spherical/diagnostics_nearest_1d_z_axis__rfm__Spherical.o Spherical/diagnostics_nearest_2d_xy_plane__rfm__Spherical.o Spherical/diagnostics_nearest_2d_yz_plane__rfm__Spherical.o Spherical/diagnostics_nearest_grid_center__rfm__Spherical.o Spherical/enforce_detgammabar_equals_detgammahat__rfm__Spherical.o Spherical/initial_data_reader__convert_ADM_Spherical_to_BSSN__rfm__Spherical.o Spherical/numerical_grid_params_Nxx_dxx_xx__rfm__Spherical.o Spherical/rfm_precompute_defines__rfm__Spherical.o Spherical/rfm_precompute_free__rfm__Spherical.o Spherical/rfm_precompute_malloc__rfm__Spherical.o Spherical/rhs_eval__rfm__Spherical.o Spherical/xx_to_Cart__rfm__Spherical.o TOVola/TOVola_interp.o TOVola/TOVola_solve.o apply_bcs_inner_only.o apply_bcs_outerextrap_and_inner.o apply_bcs_outerradiation_and_inner.o bcstruct_set_up.o bhah_diagnostics.o bhah_evolve.o bhah_finalize.o bhah_initialize.o cfl_limited_timestep.o commondata_struct_set_to_default.o constraints_eval.o diagnostics.o diagnostics_nearest_1d_y_axis.o diagnostics_nearest_1d_z_axis.o diagnostics_nearest_2d_xy_plane.o diagnostics_nearest_2d_yz_plane.o diagnostics_nearest_grid_center.o enforce_detgammabar_equals_detgammahat.o initial_data.o initial_data_reader__convert_ADM_Spherical_to_BSSN.o numerical_grid_params_Nxx_dxx_xx.o numerical_grids_and_timestep.o params_struct_set_to_default.o progress_indicator.o read_checkpoint.o rfm_precompute_defines.o rfm_precompute_free.o rfm_precompute_malloc.o rhs_eval.o write_checkpoint.o xx_to_Cart.o

# Update target to produce a static library (.a)
all: libbhah_lib.a

%.o: %.c $(COMMON_HEADERS)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

# Use ar to archive the object files into a static library
libbhah_lib.a: $(OBJ_FILES)
	ar rcs $@ $^

valgrind: clean
	$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all

# Use $(RM) to be cross-platform compatible.
clean:
	$(RM) *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png libbhah_lib.a
