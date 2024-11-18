#ifdef __cplusplus
extern "C" {
#endif

#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

void BHaH_setup(const int nxx0, const int nxx1, const int nxx2,
                const REAL cfl, const REAL rmax, BHaH_struct *bhahstruct);
// void BHaH_initialize(BHaH_struct *bhahstruct);
// double BHaH_get_timestep(BHaH_struct *bhahstruct); // this one should be deprecated
// void BHaH_advance_timestep(const REAL t_final, BHaH_struct *bhahstruct);
void BHaH_get_metric_extrinsic_curvature(BHaH_struct *bhahstruct, REAL Cartx, REAL Carty, REAL Cartz,
                                         REAL *alpha_out, REAL *beta0, REAL *beta1, REAL *beta2,
                                         REAL gammaDD[3][3], REAL KDD[3][3]);
// void BHaH_set_Tmunu(const int rgrid, const REAL rmax, REAL *in_Tmunu, BHaH_struct *bhahstruct);  // this function will not be inclulded since it is not used.
void BHaH_output_file(char *filename, BHaH_struct *bhahstruct, int dim);
int BHaH_get_gridpoints(int *indices, REAL *xCartGrid, REAL xCartMax[3], BHaH_struct *bhahstruct);
void BHaH_set_Tmunu_gridpoints(const int nCartGrid, int *indices, REAL *TmunuGrid, BHaH_struct *bhahstruct);
void BHaH_set_TOV(BHaH_struct *bhahstruct);

#ifdef __cplusplus
}
#endif
