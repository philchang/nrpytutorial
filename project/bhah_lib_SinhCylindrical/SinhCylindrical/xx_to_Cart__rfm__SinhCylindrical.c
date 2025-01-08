#include "../BHaH_defines.h"
/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
 * local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center.
 */
void xx_to_Cart__rfm__SinhCylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                      const int i0, const int i1, const int i2, REAL xCart[3]) {
#include "../set_CodeParameters.h"

  const REAL xx0 = xx[0][i0];
  const REAL xx1 = xx[1][i1];
  const REAL xx2 = xx[2][i2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = AMPLRHO*(exp(xx0/SINHWRHO) - exp(-xx0/SINHWRHO))*cos(xx1)/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)) + Cart_originx]"
   *  "[xCart[1] = AMPLRHO*(exp(xx0/SINHWRHO) - exp(-xx0/SINHWRHO))*sin(xx1)/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)) + Cart_originy]"
   *  "[xCart[2] = AMPLZ*(exp(xx2/SINHWZ) - exp(-xx2/SINHWZ))/(exp(1/SINHWZ) - exp(-1/SINHWZ)) + Cart_originz]"
   */
  {
    const REAL tmp0 = (1.0 / (SINHWRHO));
    const REAL tmp3 = (1.0 / (SINHWZ));
    const REAL tmp2 = AMPLRHO * (exp(tmp0 * xx0) - exp(-tmp0 * xx0)) / (exp(tmp0) - exp(-tmp0));
    xCart[0] = Cart_originx + tmp2 * cos(xx1);
    xCart[1] = Cart_originy + tmp2 * sin(xx1);
    xCart[2] = AMPLZ * (exp(tmp3 * xx2) - exp(-tmp3 * xx2)) / (exp(tmp3) - exp(-tmp3)) + Cart_originz;
  }
}
