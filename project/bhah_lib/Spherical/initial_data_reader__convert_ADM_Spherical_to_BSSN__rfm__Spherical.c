#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"

// ADM variables in the Cartesian basis:
typedef struct __ADM_Cart_basis_struct__ {
  REAL alpha, betaU0, betaU1, betaU2, BU0, BU1, BU2;
  REAL gammaDD00, gammaDD01, gammaDD02, gammaDD11, gammaDD12, gammaDD22;
  REAL KDD00, KDD01, KDD02, KDD11, KDD12, KDD22;

  REAL T4UU00, T4UU01, T4UU02, T4UU03;
  REAL T4UU11, T4UU12, T4UU13;
  REAL T4UU22, T4UU23;
  REAL T4UU33;
} ADM_Cart_basis_struct;

// BSSN variables in the Cartesian basis:
typedef struct __BSSN_Cart_basis_struct__ {
  REAL alpha, betaU0, betaU1, betaU2, BU0, BU1, BU2;
  REAL cf, trK;
  REAL gammabarDD00, gammabarDD01, gammabarDD02, gammabarDD11, gammabarDD12, gammabarDD22;
  REAL AbarDD00, AbarDD01, AbarDD02, AbarDD11, AbarDD12, AbarDD22;

  REAL T4UU00, T4UU01, T4UU02, T4UU03;
  REAL T4UU11, T4UU12, T4UU13;
  REAL T4UU22, T4UU23;
  REAL T4UU33;
} BSSN_Cart_basis_struct;

// Rescaled BSSN variables in the rfm basis:
typedef struct __rescaled_BSSN_rfm_basis_struct__ {
  REAL alpha, vetU0, vetU1, vetU2, betU0, betU1, betU2;
  REAL cf, trK;
  REAL hDD00, hDD01, hDD02, hDD11, hDD12, hDD22;
  REAL aDD00, aDD01, aDD02, aDD11, aDD12, aDD22;

  REAL T4UU00, T4UU01, T4UU02, T4UU03;
  REAL T4UU11, T4UU12, T4UU13;
  REAL T4UU22, T4UU23;
  REAL T4UU33;
} rescaled_BSSN_rfm_basis_struct;
/**
 * Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis
 */
static void ADM_SphorCart_to_Cart(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                                  const initial_data_struct *restrict initial_data, ADM_Cart_basis_struct *restrict ADM_Cart_basis) {
#include "../set_CodeParameters.h"
  // Unpack initial_data for ADM vectors/tensors
  const REAL betaSphorCartU0 = initial_data->betaSphorCartU0;
  const REAL betaSphorCartU1 = initial_data->betaSphorCartU1;
  const REAL betaSphorCartU2 = initial_data->betaSphorCartU2;

  const REAL BSphorCartU0 = initial_data->BSphorCartU0;
  const REAL BSphorCartU1 = initial_data->BSphorCartU1;
  const REAL BSphorCartU2 = initial_data->BSphorCartU2;

  const REAL gammaSphorCartDD00 = initial_data->gammaSphorCartDD00;
  const REAL gammaSphorCartDD01 = initial_data->gammaSphorCartDD01;
  const REAL gammaSphorCartDD02 = initial_data->gammaSphorCartDD02;
  const REAL gammaSphorCartDD11 = initial_data->gammaSphorCartDD11;
  const REAL gammaSphorCartDD12 = initial_data->gammaSphorCartDD12;
  const REAL gammaSphorCartDD22 = initial_data->gammaSphorCartDD22;

  const REAL KSphorCartDD00 = initial_data->KSphorCartDD00;
  const REAL KSphorCartDD01 = initial_data->KSphorCartDD01;
  const REAL KSphorCartDD02 = initial_data->KSphorCartDD02;
  const REAL KSphorCartDD11 = initial_data->KSphorCartDD11;
  const REAL KSphorCartDD12 = initial_data->KSphorCartDD12;
  const REAL KSphorCartDD22 = initial_data->KSphorCartDD22;

  const REAL T4SphorCartUU00 = initial_data->T4SphorCartUU00;
  const REAL T4SphorCartUU01 = initial_data->T4SphorCartUU01;
  const REAL T4SphorCartUU02 = initial_data->T4SphorCartUU02;
  const REAL T4SphorCartUU03 = initial_data->T4SphorCartUU03;
  const REAL T4SphorCartUU11 = initial_data->T4SphorCartUU11;
  const REAL T4SphorCartUU12 = initial_data->T4SphorCartUU12;
  const REAL T4SphorCartUU13 = initial_data->T4SphorCartUU13;
  const REAL T4SphorCartUU22 = initial_data->T4SphorCartUU22;
  const REAL T4SphorCartUU23 = initial_data->T4SphorCartUU23;
  const REAL T4SphorCartUU33 = initial_data->T4SphorCartUU33;

  // Perform the basis transform on ADM vectors/tensors from Spherical to Cartesian:

  // Set destination xx[3] based on desired xCart[3]
  REAL xx0, xx1, xx2;
  /*
   *  Original SymPy expressions:
   *  "[xx0 = sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2)]"
   *  "[xx1 = acos(xCart[2]/sqrt(xCart[0]**2 + xCart[1]**2 + xCart[2]**2))]"
   *  "[xx2 = atan2(xCart[1], xCart[0])]"
   */
  {
    const REAL tmp0 = sqrt(((xCart[0]) * (xCart[0])) + ((xCart[1]) * (xCart[1])) + ((xCart[2]) * (xCart[2])));
    xx0 = tmp0;
    xx1 = acos(xCart[2] / tmp0);
    xx2 = atan2(xCart[1], xCart[0]);
  }
  const REAL tmp0 = cos(xx2);
  const REAL tmp1 = sin(xx1);
  const REAL tmp4 = cos(xx1);
  const REAL tmp6 = sin(xx2);
  const REAL tmp12 = ((xx0) * (xx0));
  const REAL tmp3 = tmp0 * xx0;
  const REAL tmp7 = tmp1 * xx0;
  const REAL tmp9 = tmp6 * xx0;
  const REAL tmp10 = ((tmp0) * (tmp0));
  const REAL tmp11 = ((tmp6) * (tmp6));
  const REAL tmp13 = ((tmp1) * (tmp1) * (tmp1));
  const REAL tmp15 = ((tmp4) * (tmp4));
  const REAL tmp16 = tmp1 * tmp12;
  const REAL tmp21 = ((tmp1) * (tmp1) * (tmp1) * (tmp1)) * ((xx0) * (xx0) * (xx0) * (xx0));
  const REAL tmp27 = ((tmp1) * (tmp1));
  const REAL tmp36 = tmp0 * tmp4;
  const REAL tmp42 = tmp0 * tmp6;
  const REAL tmp43 = 2 * tmp4;
  const REAL tmp52 = tmp4 * tmp6;
  const REAL tmp89 = T4SphorCartUU22 * tmp12;
  const REAL tmp18 = tmp10 * tmp15;
  const REAL tmp28 = tmp12 * tmp27;
  const REAL tmp31 = tmp27 * tmp9;
  const REAL tmp37 = KSphorCartDD12 * tmp7;
  const REAL tmp44 = tmp42 * tmp43;
  const REAL tmp56 = tmp36 * tmp7;
  const REAL tmp78 = tmp52 * tmp7;
  const REAL tmp86 = T4SphorCartUU11 * tmp27;
  const REAL tmp19 = (1.0 / ((tmp10 * tmp12 * tmp13 + tmp11 * tmp12 * tmp13 + tmp11 * tmp15 * tmp16 + tmp16 * tmp18) *
                             (tmp10 * tmp12 * tmp13 + tmp11 * tmp12 * tmp13 + tmp11 * tmp15 * tmp16 + tmp16 * tmp18)));
  const REAL tmp24 = 2 * tmp10 * tmp4;
  const REAL tmp33 = -tmp15 * tmp9 - tmp31;
  const REAL tmp40 = KSphorCartDD02 * tmp28;
  const REAL tmp48 = tmp15 * tmp3 + tmp27 * tmp3;
  const REAL tmp58 = tmp11 * tmp27 * xx0;
  const REAL tmp59 = tmp10 * tmp27 * xx0;
  const REAL tmp65 = tmp11 * tmp16 * tmp4;
  const REAL tmp66 = tmp10 * tmp16 * tmp4;
  const REAL tmp70 = tmp0 * tmp28;
  const REAL tmp79 = tmp28 * tmp6;
  const REAL tmp88 = 2 * T4SphorCartUU13 * tmp0 * tmp31;
  const REAL tmp90 = T4SphorCartUU33 * tmp28;
  const REAL tmp91 = T4SphorCartUU23 * tmp16 * tmp44;
  const REAL tmp92 = T4SphorCartUU12 * tmp43 * tmp7;
  const REAL tmp25 = tmp13 * tmp19 * ((xx0) * (xx0) * (xx0));
  const REAL tmp29 = tmp19 * tmp28;
  const REAL tmp35 = KSphorCartDD22 * tmp19;
  const REAL tmp38 = tmp19 * tmp33;
  const REAL tmp49 = tmp19 * tmp48;
  const REAL tmp60 = -tmp58 - tmp59;
  const REAL tmp67 = tmp65 + tmp66;
  const REAL tmp101 = gammaSphorCartDD00 * tmp19 * tmp21;
  const REAL tmp104 = gammaSphorCartDD22 * tmp19;
  const REAL tmp22 = KSphorCartDD00 * tmp19 * tmp21;
  const REAL tmp26 = KSphorCartDD01 * tmp25;
  const REAL tmp30 = KSphorCartDD11 * tmp29;
  const REAL tmp39 = 2 * tmp38;
  const REAL tmp68 = tmp19 * tmp67;
  const REAL tmp71 = KSphorCartDD00 * tmp19 * tmp67;
  const REAL tmp77 = 2 * tmp49;
  const REAL tmp80 = tmp19 * ((tmp60) * (tmp60));
  const REAL tmp81 = tmp19 * ((tmp67) * (tmp67));
  const REAL tmp102 = gammaSphorCartDD01 * tmp25;
  const REAL tmp103 = gammaSphorCartDD11 * tmp29;
  const REAL tmp62 = KSphorCartDD11 * tmp19 * tmp60;
  const REAL tmp64 = KSphorCartDD01 * tmp29 * tmp60;
  const REAL tmp107 = gammaSphorCartDD11 * tmp19 * tmp60;
  const REAL tmp109 = gammaSphorCartDD01 * tmp29 * tmp60;
  ADM_Cart_basis->BU0 = BSphorCartU0 * tmp0 * tmp1 + BSphorCartU1 * tmp3 * tmp4 - BSphorCartU2 * tmp6 * tmp7;
  ADM_Cart_basis->BU1 = BSphorCartU0 * tmp1 * tmp6 + BSphorCartU1 * tmp4 * tmp9 + BSphorCartU2 * tmp0 * tmp7;
  ADM_Cart_basis->BU2 = BSphorCartU0 * tmp4 - BSphorCartU1 * tmp7;
  ADM_Cart_basis->KDD00 =
      2 * tmp0 * tmp38 * tmp40 + tmp10 * tmp22 + tmp18 * tmp30 + tmp24 * tmp26 + ((tmp33) * (tmp33)) * tmp35 + tmp36 * tmp37 * tmp39;
  ADM_Cart_basis->KDD01 = tmp0 * tmp40 * tmp49 + tmp15 * tmp30 * tmp42 + tmp22 * tmp42 + tmp26 * tmp44 + tmp33 * tmp35 * tmp48 +
                          tmp36 * tmp37 * tmp49 + tmp37 * tmp38 * tmp52 + tmp38 * tmp40 * tmp6;
  ADM_Cart_basis->KDD02 =
      KSphorCartDD01 * tmp56 * tmp68 + KSphorCartDD02 * tmp38 * tmp67 + KSphorCartDD12 * tmp38 * tmp60 + tmp0 * tmp64 + tmp56 * tmp62 + tmp70 * tmp71;
  ADM_Cart_basis->KDD11 =
      tmp11 * tmp15 * tmp30 + tmp11 * tmp22 + tmp11 * tmp26 * tmp43 + tmp35 * ((tmp48) * (tmp48)) + tmp37 * tmp52 * tmp77 + tmp40 * tmp6 * tmp77;
  ADM_Cart_basis->KDD12 =
      KSphorCartDD01 * tmp68 * tmp78 + KSphorCartDD02 * tmp49 * tmp67 + KSphorCartDD12 * tmp49 * tmp60 + tmp6 * tmp64 + tmp62 * tmp78 + tmp71 * tmp79;
  ADM_Cart_basis->KDD22 = KSphorCartDD00 * tmp81 + 2 * KSphorCartDD01 * tmp60 * tmp68 + KSphorCartDD11 * tmp80;
  ADM_Cart_basis->T4UU00 = T4SphorCartUU00;
  ADM_Cart_basis->T4UU01 = T4SphorCartUU01 * tmp0 * tmp1 + T4SphorCartUU02 * tmp3 * tmp4 - T4SphorCartUU03 * tmp6 * tmp7;
  ADM_Cart_basis->T4UU02 = T4SphorCartUU01 * tmp1 * tmp6 + T4SphorCartUU02 * tmp4 * tmp9 + T4SphorCartUU03 * tmp0 * tmp7;
  ADM_Cart_basis->T4UU03 = T4SphorCartUU01 * tmp4 - T4SphorCartUU02 * tmp7;
  ADM_Cart_basis->T4UU11 = T4SphorCartUU12 * tmp24 * tmp7 + tmp10 * tmp86 + tmp11 * tmp90 + tmp18 * tmp89 - tmp88 - tmp91;
  ADM_Cart_basis->T4UU12 = -T4SphorCartUU13 * tmp58 + T4SphorCartUU13 * tmp59 - T4SphorCartUU23 * tmp65 + T4SphorCartUU23 * tmp66 +
                           tmp15 * tmp42 * tmp89 + tmp42 * tmp86 - tmp42 * tmp90 + tmp42 * tmp92;
  ADM_Cart_basis->T4UU13 = T4SphorCartUU11 * tmp1 * tmp36 + T4SphorCartUU12 * tmp15 * tmp3 - T4SphorCartUU12 * tmp27 * tmp3 -
                           T4SphorCartUU13 * tmp52 * tmp7 + T4SphorCartUU23 * tmp28 * tmp6 - tmp1 * tmp36 * tmp89;
  ADM_Cart_basis->T4UU22 = tmp10 * tmp90 + tmp11 * tmp15 * tmp89 + tmp11 * tmp86 + tmp11 * tmp92 + tmp88 + tmp91;
  ADM_Cart_basis->T4UU23 = T4SphorCartUU11 * tmp1 * tmp52 + T4SphorCartUU12 * tmp15 * tmp9 - T4SphorCartUU12 * tmp31 +
                           T4SphorCartUU13 * tmp36 * tmp7 - T4SphorCartUU23 * tmp0 * tmp28 - tmp1 * tmp52 * tmp89;
  ADM_Cart_basis->T4UU33 = T4SphorCartUU11 * tmp15 + tmp27 * tmp89 - tmp92;
  ADM_Cart_basis->alpha = initial_data->alpha;
  ADM_Cart_basis->betaU0 = betaSphorCartU0 * tmp0 * tmp1 + betaSphorCartU1 * tmp3 * tmp4 - betaSphorCartU2 * tmp6 * tmp7;
  ADM_Cart_basis->betaU1 = betaSphorCartU0 * tmp1 * tmp6 + betaSphorCartU1 * tmp4 * tmp9 + betaSphorCartU2 * tmp0 * tmp7;
  ADM_Cart_basis->betaU2 = betaSphorCartU0 * tmp4 - betaSphorCartU1 * tmp7;
  ADM_Cart_basis->gammaDD00 = gammaSphorCartDD02 * tmp39 * tmp70 + gammaSphorCartDD12 * tmp39 * tmp56 + tmp10 * tmp101 + tmp102 * tmp24 +
                              tmp103 * tmp18 + tmp104 * ((tmp33) * (tmp33));
  ADM_Cart_basis->gammaDD01 = gammaSphorCartDD02 * tmp0 * tmp28 * tmp49 + gammaSphorCartDD02 * tmp28 * tmp38 * tmp6 +
                              gammaSphorCartDD12 * tmp36 * tmp49 * tmp7 + gammaSphorCartDD12 * tmp38 * tmp52 * tmp7 + tmp101 * tmp42 +
                              tmp102 * tmp44 + tmp103 * tmp15 * tmp42 + tmp104 * tmp33 * tmp48;
  ADM_Cart_basis->gammaDD02 = gammaSphorCartDD00 * tmp68 * tmp70 + gammaSphorCartDD01 * tmp56 * tmp68 + gammaSphorCartDD02 * tmp38 * tmp67 +
                              gammaSphorCartDD12 * tmp38 * tmp60 + tmp0 * tmp109 + tmp107 * tmp56;
  ADM_Cart_basis->gammaDD11 = gammaSphorCartDD02 * tmp77 * tmp79 + gammaSphorCartDD12 * tmp77 * tmp78 + tmp101 * tmp11 + tmp102 * tmp11 * tmp43 +
                              tmp103 * tmp11 * tmp15 + tmp104 * ((tmp48) * (tmp48));
  ADM_Cart_basis->gammaDD12 = gammaSphorCartDD00 * tmp68 * tmp79 + gammaSphorCartDD01 * tmp68 * tmp78 + gammaSphorCartDD02 * tmp49 * tmp67 +
                              gammaSphorCartDD12 * tmp49 * tmp60 + tmp107 * tmp78 + tmp109 * tmp6;
  ADM_Cart_basis->gammaDD22 = gammaSphorCartDD00 * tmp81 + 2 * gammaSphorCartDD01 * tmp60 * tmp68 + gammaSphorCartDD11 * tmp80;
}
/**
 * Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis
 */
static void ADM_Cart_to_BSSN_Cart(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                                  const ADM_Cart_basis_struct *restrict ADM_Cart_basis, BSSN_Cart_basis_struct *restrict BSSN_Cart_basis) {

  // *In the Cartesian basis*, convert ADM quantities gammaDD & KDD
  //   into BSSN gammabarDD, AbarDD, cf, and trK.
  BSSN_Cart_basis->alpha = ADM_Cart_basis->alpha;
  BSSN_Cart_basis->betaU0 = ADM_Cart_basis->betaU0;
  BSSN_Cart_basis->betaU1 = ADM_Cart_basis->betaU1;
  BSSN_Cart_basis->betaU2 = ADM_Cart_basis->betaU2;
  BSSN_Cart_basis->BU0 = ADM_Cart_basis->BU0;
  BSSN_Cart_basis->BU1 = ADM_Cart_basis->BU1;
  BSSN_Cart_basis->BU2 = ADM_Cart_basis->BU2;
  const REAL tmp1 = ADM_Cart_basis->gammaDD00 * ((ADM_Cart_basis->gammaDD12) * (ADM_Cart_basis->gammaDD12));
  const REAL tmp3 = ((ADM_Cart_basis->gammaDD01) * (ADM_Cart_basis->gammaDD01)) * ADM_Cart_basis->gammaDD22;
  const REAL tmp5 = ((ADM_Cart_basis->gammaDD02) * (ADM_Cart_basis->gammaDD02)) * ADM_Cart_basis->gammaDD11;
  const REAL tmp6 = ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD11 * ADM_Cart_basis->gammaDD22 +
                    2 * ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD02 * ADM_Cart_basis->gammaDD12 - tmp1 - tmp3 - tmp5;
  const REAL tmp7 = (1.0 / (tmp6));
  const REAL tmp8 = cbrt(tmp7);
  const REAL tmp9 = 2 * tmp7;
  const REAL tmp10 =
      ADM_Cart_basis->KDD00 * tmp7 *
          (ADM_Cart_basis->gammaDD11 * ADM_Cart_basis->gammaDD22 - ((ADM_Cart_basis->gammaDD12) * (ADM_Cart_basis->gammaDD12))) +
      ADM_Cart_basis->KDD01 * tmp9 *
          (-ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD22 + ADM_Cart_basis->gammaDD02 * ADM_Cart_basis->gammaDD12) +
      ADM_Cart_basis->KDD02 * tmp9 * (ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD12 - ADM_Cart_basis->gammaDD02 * ADM_Cart_basis->gammaDD11) +
      ADM_Cart_basis->KDD11 * tmp7 *
          (ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD22 - ((ADM_Cart_basis->gammaDD02) * (ADM_Cart_basis->gammaDD02))) +
      ADM_Cart_basis->KDD12 * tmp9 *
          (-ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD12 + ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD02) +
      ADM_Cart_basis->KDD22 * tmp7 *
          (ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD11 - ((ADM_Cart_basis->gammaDD01) * (ADM_Cart_basis->gammaDD01)));
  const REAL tmp11 = (1.0 / 3.0) * tmp10;
  BSSN_Cart_basis->AbarDD00 = tmp8 * (ADM_Cart_basis->KDD00 - ADM_Cart_basis->gammaDD00 * tmp11);
  BSSN_Cart_basis->AbarDD01 = tmp8 * (ADM_Cart_basis->KDD01 - ADM_Cart_basis->gammaDD01 * tmp11);
  BSSN_Cart_basis->AbarDD02 = tmp8 * (ADM_Cart_basis->KDD02 - ADM_Cart_basis->gammaDD02 * tmp11);
  BSSN_Cart_basis->AbarDD11 = tmp8 * (ADM_Cart_basis->KDD11 - ADM_Cart_basis->gammaDD11 * tmp11);
  BSSN_Cart_basis->AbarDD12 = tmp8 * (ADM_Cart_basis->KDD12 - ADM_Cart_basis->gammaDD12 * tmp11);
  BSSN_Cart_basis->AbarDD22 = tmp8 * (ADM_Cart_basis->KDD22 - ADM_Cart_basis->gammaDD22 * tmp11);
  BSSN_Cart_basis->T4UU00 = ADM_Cart_basis->T4UU00;
  BSSN_Cart_basis->T4UU01 = ADM_Cart_basis->T4UU01;
  BSSN_Cart_basis->T4UU02 = ADM_Cart_basis->T4UU02;
  BSSN_Cart_basis->T4UU03 = ADM_Cart_basis->T4UU03;
  BSSN_Cart_basis->T4UU11 = ADM_Cart_basis->T4UU11;
  BSSN_Cart_basis->T4UU12 = ADM_Cart_basis->T4UU12;
  BSSN_Cart_basis->T4UU13 = ADM_Cart_basis->T4UU13;
  BSSN_Cart_basis->T4UU22 = ADM_Cart_basis->T4UU22;
  BSSN_Cart_basis->T4UU23 = ADM_Cart_basis->T4UU23;
  BSSN_Cart_basis->T4UU33 = ADM_Cart_basis->T4UU33;
  BSSN_Cart_basis->cf = (1.0 / (sqrt(cbrt(tmp6 / (ADM_Cart_basis->gammaDD00 * ADM_Cart_basis->gammaDD11 * ADM_Cart_basis->gammaDD22 * tmp7 +
                                                  2 * ADM_Cart_basis->gammaDD01 * ADM_Cart_basis->gammaDD02 * ADM_Cart_basis->gammaDD12 * tmp7 -
                                                  tmp1 * tmp7 - tmp3 * tmp7 - tmp5 * tmp7)))));
  BSSN_Cart_basis->gammabarDD00 = ADM_Cart_basis->gammaDD00 * tmp8;
  BSSN_Cart_basis->gammabarDD01 = ADM_Cart_basis->gammaDD01 * tmp8;
  BSSN_Cart_basis->gammabarDD02 = ADM_Cart_basis->gammaDD02 * tmp8;
  BSSN_Cart_basis->gammabarDD11 = ADM_Cart_basis->gammaDD11 * tmp8;
  BSSN_Cart_basis->gammabarDD12 = ADM_Cart_basis->gammaDD12 * tmp8;
  BSSN_Cart_basis->gammabarDD22 = ADM_Cart_basis->gammaDD22 * tmp8;
  BSSN_Cart_basis->trK = tmp10;
}
/**
 * Cartesian -> Spherical basis transformation of BSSN vectors/tensors *except* lambda^i.
 * After the basis transform, all BSSN quantities are rescaled.
 */
static void BSSN_Cart_to_rescaled_BSSN_rfm(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                                           const BSSN_Cart_basis_struct *restrict BSSN_Cart_basis,
                                           rescaled_BSSN_rfm_basis_struct *restrict rescaled_BSSN_rfm_basis) {
#include "../set_CodeParameters.h"

  REAL xx0, xx1, xx2 __attribute__((unused)); // xx2 might be unused in the case of axisymmetric initial data.
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0 = xx[0];
    xx1 = xx[1];
    xx2 = xx[2];
  }
  const REAL tmp0 = sin(xx2);
  const REAL tmp2 = ((xx0) * (xx0));
  const REAL tmp3 = sin(xx1);
  const REAL tmp6 = cos(xx2);
  const REAL tmp8 = cos(xx1);
  const REAL tmp76 = (1.0 / (xx0));
  const REAL tmp1 = ((tmp0) * (tmp0));
  const REAL tmp4 = ((tmp3) * (tmp3) * (tmp3));
  const REAL tmp7 = ((tmp6) * (tmp6));
  const REAL tmp9 = ((tmp8) * (tmp8));
  const REAL tmp10 = tmp2 * tmp3;
  const REAL tmp15 = ((tmp3) * (tmp3));
  const REAL tmp25 = tmp3 * xx0;
  const REAL tmp45 = ((tmp3) * (tmp3) * (tmp3) * (tmp3)) * ((xx0) * (xx0) * (xx0) * (xx0));
  const REAL tmp48 = 2 * tmp6;
  const REAL tmp68 = 2 * tmp0;
  const REAL tmp80 = (1.0 / (tmp3));
  const REAL tmp82 = (1.0 / (tmp2));
  const REAL tmp12 = tmp7 * tmp9;
  const REAL tmp16 = tmp15 * tmp2;
  const REAL tmp22 = tmp7 * tmp8;
  const REAL tmp26 = tmp25 * tmp8;
  const REAL tmp31 = tmp15 * xx0;
  const REAL tmp57 = tmp4 * ((xx0) * (xx0) * (xx0));
  const REAL tmp70 = tmp3 * tmp48 * tmp8;
  const REAL tmp71 = tmp3 * tmp68 * tmp8;
  const REAL tmp74 = BSSN_Cart_basis->AbarDD01 * tmp48;
  const REAL tmp83 = tmp10 * tmp8;
  const REAL tmp87 = tmp2 * tmp9;
  const REAL tmp92 = tmp82 / tmp15;
  const REAL tmp94 = BSSN_Cart_basis->gammabarDD01 * tmp48;
  const REAL tmp13 = tmp1 * tmp10 * tmp9 + tmp1 * tmp2 * tmp4 + tmp10 * tmp12 + tmp2 * tmp4 * tmp7;
  const REAL tmp17 = tmp16 * tmp6;
  const REAL tmp19 = tmp0 * tmp16;
  const REAL tmp21 = tmp1 * tmp10 * tmp8;
  const REAL tmp27 = tmp26 * tmp6;
  const REAL tmp29 = tmp0 * tmp26;
  const REAL tmp32 = tmp1 * tmp31;
  const REAL tmp34 = tmp31 * tmp6;
  const REAL tmp36 = tmp6 * tmp9 * xx0;
  const REAL tmp39 = tmp0 * tmp31;
  const REAL tmp40 = tmp0 * tmp9 * xx0;
  const REAL tmp66 = tmp1 * tmp16;
  const REAL tmp89 = tmp0 * tmp6 * tmp83;
  const REAL tmp14 = (1.0 / (tmp13));
  const REAL tmp23 = tmp10 * tmp22 + tmp21;
  const REAL tmp33 = -tmp31 * tmp7 - tmp32;
  const REAL tmp37 = tmp34 + tmp36;
  const REAL tmp41 = -tmp39 - tmp40;
  const REAL tmp43 = (1.0 / ((tmp13) * (tmp13)));
  const REAL tmp18 = tmp14 * tmp17;
  const REAL tmp20 = tmp14 * tmp19;
  const REAL tmp28 = tmp14 * tmp27;
  const REAL tmp30 = tmp14 * tmp29;
  const REAL tmp38 = tmp14 * tmp37;
  const REAL tmp42 = tmp14 * tmp41;
  const REAL tmp44 = BSSN_Cart_basis->T4UU11 * tmp43;
  const REAL tmp46 = BSSN_Cart_basis->T4UU22 * tmp43;
  const REAL tmp52 = BSSN_Cart_basis->T4UU33 * tmp43;
  const REAL tmp54 = tmp23 * tmp43;
  const REAL tmp59 = tmp33 * tmp43;
  const REAL tmp60 = BSSN_Cart_basis->T4UU12 * tmp37 * tmp43;
  const REAL tmp50 = BSSN_Cart_basis->T4UU12 * tmp43 * tmp48;
  const REAL tmp63 = BSSN_Cart_basis->T4UU12 * tmp41 * tmp43;
  rescaled_BSSN_rfm_basis->T4UU00 = BSSN_Cart_basis->T4UU00;
  rescaled_BSSN_rfm_basis->T4UU01 = BSSN_Cart_basis->T4UU01 * tmp18 + BSSN_Cart_basis->T4UU02 * tmp20 + BSSN_Cart_basis->T4UU03 * tmp14 * tmp23;
  rescaled_BSSN_rfm_basis->T4UU02 = BSSN_Cart_basis->T4UU01 * tmp28 + BSSN_Cart_basis->T4UU02 * tmp30 + BSSN_Cart_basis->T4UU03 * tmp14 * tmp33;
  rescaled_BSSN_rfm_basis->T4UU03 = BSSN_Cart_basis->T4UU01 * tmp42 + BSSN_Cart_basis->T4UU02 * tmp38;
  rescaled_BSSN_rfm_basis->T4UU11 = 2 * BSSN_Cart_basis->T4UU13 * tmp17 * tmp54 + 2 * BSSN_Cart_basis->T4UU23 * tmp19 * tmp54 + tmp0 * tmp45 * tmp50 +
                                    tmp1 * tmp45 * tmp46 + ((tmp23) * (tmp23)) * tmp52 + tmp44 * tmp45 * tmp7;
  rescaled_BSSN_rfm_basis->T4UU12 = BSSN_Cart_basis->T4UU13 * tmp17 * tmp59 + BSSN_Cart_basis->T4UU13 * tmp27 * tmp54 +
                                    BSSN_Cart_basis->T4UU23 * tmp19 * tmp59 + BSSN_Cart_basis->T4UU23 * tmp29 * tmp54 + tmp0 * tmp50 * tmp57 * tmp8 +
                                    tmp1 * tmp46 * tmp57 * tmp8 + tmp22 * tmp44 * tmp57 + tmp23 * tmp33 * tmp52;
  rescaled_BSSN_rfm_basis->T4UU13 = BSSN_Cart_basis->T4UU13 * tmp41 * tmp54 + BSSN_Cart_basis->T4UU23 * tmp37 * tmp54 + tmp17 * tmp41 * tmp44 +
                                    tmp17 * tmp60 + tmp19 * tmp37 * tmp46 + tmp19 * tmp63;
  rescaled_BSSN_rfm_basis->T4UU22 = BSSN_Cart_basis->T4UU13 * tmp26 * tmp48 * tmp59 + BSSN_Cart_basis->T4UU23 * tmp26 * tmp59 * tmp68 +
                                    tmp12 * tmp16 * tmp44 + tmp19 * tmp50 * tmp9 + ((tmp33) * (tmp33)) * tmp52 + tmp46 * tmp66 * tmp9;
  rescaled_BSSN_rfm_basis->T4UU23 = BSSN_Cart_basis->T4UU13 * tmp41 * tmp59 + BSSN_Cart_basis->T4UU23 * tmp37 * tmp59 + tmp27 * tmp41 * tmp44 +
                                    tmp27 * tmp60 + tmp29 * tmp37 * tmp46 + tmp29 * tmp63;
  rescaled_BSSN_rfm_basis->T4UU33 = ((tmp37) * (tmp37)) * tmp46 + ((tmp41) * (tmp41)) * tmp44 + 2 * tmp41 * tmp60;
  rescaled_BSSN_rfm_basis->aDD00 = BSSN_Cart_basis->AbarDD00 * tmp15 * tmp7 + BSSN_Cart_basis->AbarDD02 * tmp70 +
                                   BSSN_Cart_basis->AbarDD11 * tmp1 * tmp15 + BSSN_Cart_basis->AbarDD12 * tmp71 + BSSN_Cart_basis->AbarDD22 * tmp9 +
                                   tmp0 * tmp15 * tmp74;
  rescaled_BSSN_rfm_basis->aDD01 =
      tmp76 * (BSSN_Cart_basis->AbarDD00 * tmp22 * tmp25 - BSSN_Cart_basis->AbarDD02 * tmp34 + BSSN_Cart_basis->AbarDD02 * tmp36 +
               BSSN_Cart_basis->AbarDD11 * tmp1 * tmp26 - BSSN_Cart_basis->AbarDD12 * tmp39 + BSSN_Cart_basis->AbarDD12 * tmp40 -
               BSSN_Cart_basis->AbarDD22 * tmp26 + tmp29 * tmp74);
  rescaled_BSSN_rfm_basis->aDD02 = tmp76 * tmp80 *
                                   (-BSSN_Cart_basis->AbarDD00 * tmp39 * tmp6 + BSSN_Cart_basis->AbarDD01 * tmp15 * tmp7 * xx0 -
                                    BSSN_Cart_basis->AbarDD01 * tmp32 - BSSN_Cart_basis->AbarDD02 * tmp29 +
                                    BSSN_Cart_basis->AbarDD11 * tmp0 * tmp15 * tmp6 * xx0 + BSSN_Cart_basis->AbarDD12 * tmp3 * tmp6 * tmp8 * xx0);
  rescaled_BSSN_rfm_basis->aDD11 =
      tmp82 * (BSSN_Cart_basis->AbarDD00 * tmp12 * tmp2 - BSSN_Cart_basis->AbarDD02 * tmp48 * tmp83 + BSSN_Cart_basis->AbarDD11 * tmp1 * tmp87 -
               BSSN_Cart_basis->AbarDD12 * tmp68 * tmp83 + BSSN_Cart_basis->AbarDD22 * tmp16 + tmp0 * tmp74 * tmp87);
  rescaled_BSSN_rfm_basis->aDD12 = tmp80 * tmp82 *
                                   (-BSSN_Cart_basis->AbarDD00 * tmp89 + BSSN_Cart_basis->AbarDD01 * tmp2 * tmp3 * tmp7 * tmp8 -
                                    BSSN_Cart_basis->AbarDD01 * tmp21 + BSSN_Cart_basis->AbarDD02 * tmp0 * tmp15 * tmp2 +
                                    BSSN_Cart_basis->AbarDD11 * tmp0 * tmp2 * tmp3 * tmp6 * tmp8 - BSSN_Cart_basis->AbarDD12 * tmp17);
  rescaled_BSSN_rfm_basis->aDD22 = tmp92 * (BSSN_Cart_basis->AbarDD00 * tmp66 + BSSN_Cart_basis->AbarDD11 * tmp16 * tmp7 - tmp19 * tmp74);
  rescaled_BSSN_rfm_basis->alpha = BSSN_Cart_basis->alpha;
  rescaled_BSSN_rfm_basis->betU0 = BSSN_Cart_basis->BU0 * tmp18 + BSSN_Cart_basis->BU1 * tmp20 + BSSN_Cart_basis->BU2 * tmp14 * tmp23;
  rescaled_BSSN_rfm_basis->betU1 = xx0 * (BSSN_Cart_basis->BU0 * tmp28 + BSSN_Cart_basis->BU1 * tmp30 + BSSN_Cart_basis->BU2 * tmp14 * tmp33);
  rescaled_BSSN_rfm_basis->betU2 = tmp25 * (BSSN_Cart_basis->BU0 * tmp42 + BSSN_Cart_basis->BU1 * tmp38);
  rescaled_BSSN_rfm_basis->cf = BSSN_Cart_basis->cf;
  rescaled_BSSN_rfm_basis->hDD00 = BSSN_Cart_basis->gammabarDD00 * tmp15 * tmp7 + BSSN_Cart_basis->gammabarDD02 * tmp70 +
                                   BSSN_Cart_basis->gammabarDD11 * tmp1 * tmp15 + BSSN_Cart_basis->gammabarDD12 * tmp71 +
                                   BSSN_Cart_basis->gammabarDD22 * tmp9 + tmp0 * tmp15 * tmp94 - 1;
  rescaled_BSSN_rfm_basis->hDD01 =
      tmp76 * (BSSN_Cart_basis->gammabarDD00 * tmp22 * tmp25 - BSSN_Cart_basis->gammabarDD02 * tmp34 + BSSN_Cart_basis->gammabarDD02 * tmp36 +
               BSSN_Cart_basis->gammabarDD11 * tmp1 * tmp26 - BSSN_Cart_basis->gammabarDD12 * tmp39 + BSSN_Cart_basis->gammabarDD12 * tmp40 -
               BSSN_Cart_basis->gammabarDD22 * tmp26 + tmp29 * tmp94);
  rescaled_BSSN_rfm_basis->hDD02 =
      tmp76 * tmp80 *
      (-BSSN_Cart_basis->gammabarDD00 * tmp39 * tmp6 + BSSN_Cart_basis->gammabarDD01 * tmp15 * tmp7 * xx0 - BSSN_Cart_basis->gammabarDD01 * tmp32 -
       BSSN_Cart_basis->gammabarDD02 * tmp29 + BSSN_Cart_basis->gammabarDD11 * tmp0 * tmp15 * tmp6 * xx0 +
       BSSN_Cart_basis->gammabarDD12 * tmp3 * tmp6 * tmp8 * xx0);
  rescaled_BSSN_rfm_basis->hDD11 = tmp82 * (BSSN_Cart_basis->gammabarDD00 * tmp12 * tmp2 - BSSN_Cart_basis->gammabarDD02 * tmp48 * tmp83 +
                                            BSSN_Cart_basis->gammabarDD11 * tmp1 * tmp87 - BSSN_Cart_basis->gammabarDD12 * tmp68 * tmp83 +
                                            BSSN_Cart_basis->gammabarDD22 * tmp16 + tmp0 * tmp87 * tmp94 - tmp2);
  rescaled_BSSN_rfm_basis->hDD12 = tmp80 * tmp82 *
                                   (-BSSN_Cart_basis->gammabarDD00 * tmp89 + BSSN_Cart_basis->gammabarDD01 * tmp2 * tmp3 * tmp7 * tmp8 -
                                    BSSN_Cart_basis->gammabarDD01 * tmp21 + BSSN_Cart_basis->gammabarDD02 * tmp0 * tmp15 * tmp2 +
                                    BSSN_Cart_basis->gammabarDD11 * tmp0 * tmp2 * tmp3 * tmp6 * tmp8 - BSSN_Cart_basis->gammabarDD12 * tmp17);
  rescaled_BSSN_rfm_basis->hDD22 =
      tmp92 * (BSSN_Cart_basis->gammabarDD00 * tmp66 + BSSN_Cart_basis->gammabarDD11 * tmp16 * tmp7 - tmp16 - tmp19 * tmp94);
  rescaled_BSSN_rfm_basis->trK = BSSN_Cart_basis->trK;
  rescaled_BSSN_rfm_basis->vetU0 = BSSN_Cart_basis->betaU0 * tmp18 + BSSN_Cart_basis->betaU1 * tmp20 + BSSN_Cart_basis->betaU2 * tmp14 * tmp23;
  rescaled_BSSN_rfm_basis->vetU1 =
      xx0 * (BSSN_Cart_basis->betaU0 * tmp28 + BSSN_Cart_basis->betaU1 * tmp30 + BSSN_Cart_basis->betaU2 * tmp14 * tmp33);
  rescaled_BSSN_rfm_basis->vetU2 = tmp25 * (BSSN_Cart_basis->betaU0 * tmp42 + BSSN_Cart_basis->betaU1 * tmp38);
}
/**
 * Compute lambdaU in Spherical coordinates
 */
static void initial_data_lambdaU_grid_interior(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                               REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "../set_CodeParameters.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
        const REAL xx0 = xx[0][i0]; /*
                                     * NRPy+-Generated GF Access/FD Code, Step 1 of 2:
                                     * Read gridfunction(s) from main memory and compute FD stencils as needed.
                                     */
        const REAL hDD00_i2m2 = in_gfs[IDX4(HDD00GF, i0, i1, i2 - 2)];
        const REAL hDD00_i2m1 = in_gfs[IDX4(HDD00GF, i0, i1, i2 - 1)];
        const REAL hDD00_i1m2 = in_gfs[IDX4(HDD00GF, i0, i1 - 2, i2)];
        const REAL hDD00_i1m1 = in_gfs[IDX4(HDD00GF, i0, i1 - 1, i2)];
        const REAL hDD00_i0m2 = in_gfs[IDX4(HDD00GF, i0 - 2, i1, i2)];
        const REAL hDD00_i0m1 = in_gfs[IDX4(HDD00GF, i0 - 1, i1, i2)];
        const REAL hDD00 = in_gfs[IDX4(HDD00GF, i0, i1, i2)];
        const REAL hDD00_i0p1 = in_gfs[IDX4(HDD00GF, i0 + 1, i1, i2)];
        const REAL hDD00_i0p2 = in_gfs[IDX4(HDD00GF, i0 + 2, i1, i2)];
        const REAL hDD00_i1p1 = in_gfs[IDX4(HDD00GF, i0, i1 + 1, i2)];
        const REAL hDD00_i1p2 = in_gfs[IDX4(HDD00GF, i0, i1 + 2, i2)];
        const REAL hDD00_i2p1 = in_gfs[IDX4(HDD00GF, i0, i1, i2 + 1)];
        const REAL hDD00_i2p2 = in_gfs[IDX4(HDD00GF, i0, i1, i2 + 2)];
        const REAL hDD01_i2m2 = in_gfs[IDX4(HDD01GF, i0, i1, i2 - 2)];
        const REAL hDD01_i2m1 = in_gfs[IDX4(HDD01GF, i0, i1, i2 - 1)];
        const REAL hDD01_i1m2 = in_gfs[IDX4(HDD01GF, i0, i1 - 2, i2)];
        const REAL hDD01_i1m1 = in_gfs[IDX4(HDD01GF, i0, i1 - 1, i2)];
        const REAL hDD01_i0m2 = in_gfs[IDX4(HDD01GF, i0 - 2, i1, i2)];
        const REAL hDD01_i0m1 = in_gfs[IDX4(HDD01GF, i0 - 1, i1, i2)];
        const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
        const REAL hDD01_i0p1 = in_gfs[IDX4(HDD01GF, i0 + 1, i1, i2)];
        const REAL hDD01_i0p2 = in_gfs[IDX4(HDD01GF, i0 + 2, i1, i2)];
        const REAL hDD01_i1p1 = in_gfs[IDX4(HDD01GF, i0, i1 + 1, i2)];
        const REAL hDD01_i1p2 = in_gfs[IDX4(HDD01GF, i0, i1 + 2, i2)];
        const REAL hDD01_i2p1 = in_gfs[IDX4(HDD01GF, i0, i1, i2 + 1)];
        const REAL hDD01_i2p2 = in_gfs[IDX4(HDD01GF, i0, i1, i2 + 2)];
        const REAL hDD02_i2m2 = in_gfs[IDX4(HDD02GF, i0, i1, i2 - 2)];
        const REAL hDD02_i2m1 = in_gfs[IDX4(HDD02GF, i0, i1, i2 - 1)];
        const REAL hDD02_i1m2 = in_gfs[IDX4(HDD02GF, i0, i1 - 2, i2)];
        const REAL hDD02_i1m1 = in_gfs[IDX4(HDD02GF, i0, i1 - 1, i2)];
        const REAL hDD02_i0m2 = in_gfs[IDX4(HDD02GF, i0 - 2, i1, i2)];
        const REAL hDD02_i0m1 = in_gfs[IDX4(HDD02GF, i0 - 1, i1, i2)];
        const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
        const REAL hDD02_i0p1 = in_gfs[IDX4(HDD02GF, i0 + 1, i1, i2)];
        const REAL hDD02_i0p2 = in_gfs[IDX4(HDD02GF, i0 + 2, i1, i2)];
        const REAL hDD02_i1p1 = in_gfs[IDX4(HDD02GF, i0, i1 + 1, i2)];
        const REAL hDD02_i1p2 = in_gfs[IDX4(HDD02GF, i0, i1 + 2, i2)];
        const REAL hDD02_i2p1 = in_gfs[IDX4(HDD02GF, i0, i1, i2 + 1)];
        const REAL hDD02_i2p2 = in_gfs[IDX4(HDD02GF, i0, i1, i2 + 2)];
        const REAL hDD11_i2m2 = in_gfs[IDX4(HDD11GF, i0, i1, i2 - 2)];
        const REAL hDD11_i2m1 = in_gfs[IDX4(HDD11GF, i0, i1, i2 - 1)];
        const REAL hDD11_i1m2 = in_gfs[IDX4(HDD11GF, i0, i1 - 2, i2)];
        const REAL hDD11_i1m1 = in_gfs[IDX4(HDD11GF, i0, i1 - 1, i2)];
        const REAL hDD11_i0m2 = in_gfs[IDX4(HDD11GF, i0 - 2, i1, i2)];
        const REAL hDD11_i0m1 = in_gfs[IDX4(HDD11GF, i0 - 1, i1, i2)];
        const REAL hDD11 = in_gfs[IDX4(HDD11GF, i0, i1, i2)];
        const REAL hDD11_i0p1 = in_gfs[IDX4(HDD11GF, i0 + 1, i1, i2)];
        const REAL hDD11_i0p2 = in_gfs[IDX4(HDD11GF, i0 + 2, i1, i2)];
        const REAL hDD11_i1p1 = in_gfs[IDX4(HDD11GF, i0, i1 + 1, i2)];
        const REAL hDD11_i1p2 = in_gfs[IDX4(HDD11GF, i0, i1 + 2, i2)];
        const REAL hDD11_i2p1 = in_gfs[IDX4(HDD11GF, i0, i1, i2 + 1)];
        const REAL hDD11_i2p2 = in_gfs[IDX4(HDD11GF, i0, i1, i2 + 2)];
        const REAL hDD12_i2m2 = in_gfs[IDX4(HDD12GF, i0, i1, i2 - 2)];
        const REAL hDD12_i2m1 = in_gfs[IDX4(HDD12GF, i0, i1, i2 - 1)];
        const REAL hDD12_i1m2 = in_gfs[IDX4(HDD12GF, i0, i1 - 2, i2)];
        const REAL hDD12_i1m1 = in_gfs[IDX4(HDD12GF, i0, i1 - 1, i2)];
        const REAL hDD12_i0m2 = in_gfs[IDX4(HDD12GF, i0 - 2, i1, i2)];
        const REAL hDD12_i0m1 = in_gfs[IDX4(HDD12GF, i0 - 1, i1, i2)];
        const REAL hDD12 = in_gfs[IDX4(HDD12GF, i0, i1, i2)];
        const REAL hDD12_i0p1 = in_gfs[IDX4(HDD12GF, i0 + 1, i1, i2)];
        const REAL hDD12_i0p2 = in_gfs[IDX4(HDD12GF, i0 + 2, i1, i2)];
        const REAL hDD12_i1p1 = in_gfs[IDX4(HDD12GF, i0, i1 + 1, i2)];
        const REAL hDD12_i1p2 = in_gfs[IDX4(HDD12GF, i0, i1 + 2, i2)];
        const REAL hDD12_i2p1 = in_gfs[IDX4(HDD12GF, i0, i1, i2 + 1)];
        const REAL hDD12_i2p2 = in_gfs[IDX4(HDD12GF, i0, i1, i2 + 2)];
        const REAL hDD22_i2m2 = in_gfs[IDX4(HDD22GF, i0, i1, i2 - 2)];
        const REAL hDD22_i2m1 = in_gfs[IDX4(HDD22GF, i0, i1, i2 - 1)];
        const REAL hDD22_i1m2 = in_gfs[IDX4(HDD22GF, i0, i1 - 2, i2)];
        const REAL hDD22_i1m1 = in_gfs[IDX4(HDD22GF, i0, i1 - 1, i2)];
        const REAL hDD22_i0m2 = in_gfs[IDX4(HDD22GF, i0 - 2, i1, i2)];
        const REAL hDD22_i0m1 = in_gfs[IDX4(HDD22GF, i0 - 1, i1, i2)];
        const REAL hDD22 = in_gfs[IDX4(HDD22GF, i0, i1, i2)];
        const REAL hDD22_i0p1 = in_gfs[IDX4(HDD22GF, i0 + 1, i1, i2)];
        const REAL hDD22_i0p2 = in_gfs[IDX4(HDD22GF, i0 + 2, i1, i2)];
        const REAL hDD22_i1p1 = in_gfs[IDX4(HDD22GF, i0, i1 + 1, i2)];
        const REAL hDD22_i1p2 = in_gfs[IDX4(HDD22GF, i0, i1 + 2, i2)];
        const REAL hDD22_i2p1 = in_gfs[IDX4(HDD22GF, i0, i1, i2 + 1)];
        const REAL hDD22_i2p2 = in_gfs[IDX4(HDD22GF, i0, i1, i2 + 2)];
        const REAL FDPart1_Rational_2_3 = 2.0 / 3.0;
        const REAL FDPart1_Rational_1_12 = 1.0 / 12.0;
        const REAL hDD_dD000 = invdxx0 * (FDPart1_Rational_1_12 * (hDD00_i0m2 - hDD00_i0p2) + FDPart1_Rational_2_3 * (-hDD00_i0m1 + hDD00_i0p1));
        const REAL hDD_dD001 = invdxx1 * (FDPart1_Rational_1_12 * (hDD00_i1m2 - hDD00_i1p2) + FDPart1_Rational_2_3 * (-hDD00_i1m1 + hDD00_i1p1));
        const REAL hDD_dD002 = invdxx2 * (FDPart1_Rational_1_12 * (hDD00_i2m2 - hDD00_i2p2) + FDPart1_Rational_2_3 * (-hDD00_i2m1 + hDD00_i2p1));
        const REAL hDD_dD010 = invdxx0 * (FDPart1_Rational_1_12 * (hDD01_i0m2 - hDD01_i0p2) + FDPart1_Rational_2_3 * (-hDD01_i0m1 + hDD01_i0p1));
        const REAL hDD_dD011 = invdxx1 * (FDPart1_Rational_1_12 * (hDD01_i1m2 - hDD01_i1p2) + FDPart1_Rational_2_3 * (-hDD01_i1m1 + hDD01_i1p1));
        const REAL hDD_dD012 = invdxx2 * (FDPart1_Rational_1_12 * (hDD01_i2m2 - hDD01_i2p2) + FDPart1_Rational_2_3 * (-hDD01_i2m1 + hDD01_i2p1));
        const REAL hDD_dD020 = invdxx0 * (FDPart1_Rational_1_12 * (hDD02_i0m2 - hDD02_i0p2) + FDPart1_Rational_2_3 * (-hDD02_i0m1 + hDD02_i0p1));
        const REAL hDD_dD021 = invdxx1 * (FDPart1_Rational_1_12 * (hDD02_i1m2 - hDD02_i1p2) + FDPart1_Rational_2_3 * (-hDD02_i1m1 + hDD02_i1p1));
        const REAL hDD_dD022 = invdxx2 * (FDPart1_Rational_1_12 * (hDD02_i2m2 - hDD02_i2p2) + FDPart1_Rational_2_3 * (-hDD02_i2m1 + hDD02_i2p1));
        const REAL hDD_dD110 = invdxx0 * (FDPart1_Rational_1_12 * (hDD11_i0m2 - hDD11_i0p2) + FDPart1_Rational_2_3 * (-hDD11_i0m1 + hDD11_i0p1));
        const REAL hDD_dD111 = invdxx1 * (FDPart1_Rational_1_12 * (hDD11_i1m2 - hDD11_i1p2) + FDPart1_Rational_2_3 * (-hDD11_i1m1 + hDD11_i1p1));
        const REAL hDD_dD112 = invdxx2 * (FDPart1_Rational_1_12 * (hDD11_i2m2 - hDD11_i2p2) + FDPart1_Rational_2_3 * (-hDD11_i2m1 + hDD11_i2p1));
        const REAL hDD_dD120 = invdxx0 * (FDPart1_Rational_1_12 * (hDD12_i0m2 - hDD12_i0p2) + FDPart1_Rational_2_3 * (-hDD12_i0m1 + hDD12_i0p1));
        const REAL hDD_dD121 = invdxx1 * (FDPart1_Rational_1_12 * (hDD12_i1m2 - hDD12_i1p2) + FDPart1_Rational_2_3 * (-hDD12_i1m1 + hDD12_i1p1));
        const REAL hDD_dD122 = invdxx2 * (FDPart1_Rational_1_12 * (hDD12_i2m2 - hDD12_i2p2) + FDPart1_Rational_2_3 * (-hDD12_i2m1 + hDD12_i2p1));
        const REAL hDD_dD220 = invdxx0 * (FDPart1_Rational_1_12 * (hDD22_i0m2 - hDD22_i0p2) + FDPart1_Rational_2_3 * (-hDD22_i0m1 + hDD22_i0p1));
        const REAL hDD_dD221 = invdxx1 * (FDPart1_Rational_1_12 * (hDD22_i1m2 - hDD22_i1p2) + FDPart1_Rational_2_3 * (-hDD22_i1m1 + hDD22_i1p1));
        const REAL hDD_dD222 = invdxx2 * (FDPart1_Rational_1_12 * (hDD22_i2m2 - hDD22_i2p2) + FDPart1_Rational_2_3 * (-hDD22_i2m1 + hDD22_i2p1));

        /*
         * NRPy+-Generated GF Access/FD Code, Step 2 of 2:
         * Evaluate SymPy expressions and write to main memory.
         */
        const REAL FDPart3tmp0 = sin(xx1);
        const REAL FDPart3tmp2 = 2 * xx0;
        const REAL FDPart3tmp5 = ((xx0) * (xx0) * (xx0));
        const REAL FDPart3tmp6 = ((xx0) * (xx0));
        const REAL FDPart3tmp10 = ((xx0) * (xx0) * (xx0) * (xx0));
        const REAL FDPart3tmp12 = hDD00 + 1;
        const REAL FDPart3tmp27 = hDD_dD012 * xx0;
        const REAL FDPart3tmp28 = cos(xx1);
        const REAL FDPart3tmp63 = -1 / xx0;
        const REAL FDPart3tmp7 = FDPart3tmp6 * hDD11 + FDPart3tmp6;
        const REAL FDPart3tmp11 = ((FDPart3tmp0) * (FDPart3tmp0));
        const REAL FDPart3tmp17 = FDPart3tmp6 * ((hDD01) * (hDD01));
        const REAL FDPart3tmp23 = FDPart3tmp2 * hDD_dD010 + 2 * hDD01 - hDD_dD001;
        const REAL FDPart3tmp26 = FDPart3tmp2 * hDD11 + FDPart3tmp2 + FDPart3tmp6 * hDD_dD110;
        const REAL FDPart3tmp29 = FDPart3tmp28 * hDD02 * xx0;
        const REAL FDPart3tmp33 = FDPart3tmp0 * FDPart3tmp6;
        const REAL FDPart3tmp45 = FDPart3tmp6 * hDD_dD111;
        const REAL FDPart3tmp46 = FDPart3tmp6 * hDD_dD112;
        const REAL FDPart3tmp4 = FDPart3tmp0 * FDPart3tmp2 * hDD_dD020 + 2 * FDPart3tmp0 * hDD02 - hDD_dD002;
        const REAL FDPart3tmp9 = FDPart3tmp0 * FDPart3tmp5 * hDD01 * hDD12 - FDPart3tmp0 * FDPart3tmp7 * hDD02 * xx0;
        const REAL FDPart3tmp13 = FDPart3tmp10 * FDPart3tmp11 * ((hDD12) * (hDD12));
        const REAL FDPart3tmp14 = FDPart3tmp11 * FDPart3tmp6;
        const REAL FDPart3tmp31 = FDPart3tmp0 * hDD_dD021 * xx0;
        const REAL FDPart3tmp35 = FDPart3tmp0 * FDPart3tmp2 * hDD12;
        const REAL FDPart3tmp48 = 2 * FDPart3tmp0 * FDPart3tmp6 * hDD_dD121 + 2 * FDPart3tmp28 * FDPart3tmp6 * hDD12 - FDPart3tmp46;
        const REAL FDPart3tmp49 = -FDPart3tmp26 + 2 * hDD_dD011 * xx0;
        const REAL FDPart3tmp52 = 2 * FDPart3tmp0 * FDPart3tmp28 * FDPart3tmp6;
        const REAL FDPart3tmp55 = -FDPart3tmp12 * FDPart3tmp33 * hDD12 + FDPart3tmp33 * hDD01 * hDD02;
        const REAL FDPart3tmp15 = FDPart3tmp14 * ((hDD02) * (hDD02));
        const REAL FDPart3tmp16 = FDPart3tmp14 * hDD22 + FDPart3tmp14;
        const REAL FDPart3tmp36 = FDPart3tmp33 * hDD_dD120 + FDPart3tmp35;
        const REAL FDPart3tmp41 = 2 * FDPart3tmp11 * xx0;
        const REAL FDPart3tmp53 = FDPart3tmp14 * hDD_dD221 + FDPart3tmp52 * hDD22 + FDPart3tmp52;
        const REAL FDPart3tmp57 = FDPart3tmp14 * hDD_dD222;
        const REAL FDPart3tmp18 = (1.0 / (2 * FDPart3tmp10 * FDPart3tmp11 * hDD01 * hDD02 * hDD12 - FDPart3tmp12 * FDPart3tmp13 +
                                          FDPart3tmp12 * FDPart3tmp16 * FDPart3tmp7 - FDPart3tmp15 * FDPart3tmp7 - FDPart3tmp16 * FDPart3tmp17));
        const REAL FDPart3tmp24 = FDPart3tmp11 * FDPart3tmp5 * hDD02 * hDD12 - FDPart3tmp16 * hDD01 * xx0;
        const REAL FDPart3tmp37 = -FDPart3tmp27 + FDPart3tmp29 + FDPart3tmp31 + FDPart3tmp36;
        const REAL FDPart3tmp42 = FDPart3tmp14 * hDD_dD220 + FDPart3tmp41 * hDD22 + FDPart3tmp41;
        const REAL FDPart3tmp43 = FDPart3tmp27 - FDPart3tmp29 - FDPart3tmp31 + FDPart3tmp36;
        const REAL FDPart3tmp54 = FDPart3tmp27 + FDPart3tmp29 + FDPart3tmp31 - FDPart3tmp33 * hDD_dD120 - FDPart3tmp35;
        const REAL FDPart3tmp58 = 2 * FDPart3tmp0 * FDPart3tmp6 * hDD_dD122 - FDPart3tmp53;
        const REAL FDPart3tmp19 = (1.0 / 2.0) * FDPart3tmp18;
        const REAL FDPart3tmp21 = FDPart3tmp18 * (-FDPart3tmp13 + FDPart3tmp16 * FDPart3tmp7);
        const REAL FDPart3tmp38 = 2 * FDPart3tmp18;
        const REAL FDPart3tmp50 = FDPart3tmp18 * (FDPart3tmp12 * FDPart3tmp16 - FDPart3tmp15);
        const REAL FDPart3tmp59 = 2 * FDPart3tmp0 * hDD_dD022 * xx0 - FDPart3tmp42;
        const REAL FDPart3tmp60 = FDPart3tmp18 * (FDPart3tmp12 * FDPart3tmp7 - FDPart3tmp17);
        const REAL FDPart3tmp20 = FDPart3tmp19 * FDPart3tmp9;
        const REAL FDPart3tmp22 = (1.0 / 2.0) * FDPart3tmp21;
        const REAL FDPart3tmp25 = FDPart3tmp19 * FDPart3tmp24;
        const REAL FDPart3tmp39 = FDPart3tmp24 * FDPart3tmp38;
        const REAL FDPart3tmp44 = FDPart3tmp38 * FDPart3tmp9;
        const REAL FDPart3tmp56 = FDPart3tmp38 * FDPart3tmp55;
        const REAL FDPart3tmp61 = FDPart3tmp19 * FDPart3tmp55;
        const REAL FDPart3tmp62 = (1.0 / 2.0) * FDPart3tmp50;
        const REAL FDPart3tmp64 = (1.0 / 2.0) * FDPart3tmp60;
        in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)] =
            FDPart3tmp21 * (FDPart3tmp20 * FDPart3tmp4 + FDPart3tmp22 * hDD_dD000 + FDPart3tmp23 * FDPart3tmp25) +
            FDPart3tmp39 * (FDPart3tmp20 * FDPart3tmp37 + FDPart3tmp22 * hDD_dD001 + FDPart3tmp25 * FDPart3tmp26) +
            FDPart3tmp44 * (FDPart3tmp20 * FDPart3tmp42 + FDPart3tmp22 * hDD_dD002 + FDPart3tmp25 * FDPart3tmp43) +
            FDPart3tmp50 * (FDPart3tmp20 * FDPart3tmp48 + FDPart3tmp22 * FDPart3tmp49 + FDPart3tmp25 * FDPart3tmp45 + xx0) +
            FDPart3tmp56 * (FDPart3tmp20 * FDPart3tmp53 + FDPart3tmp22 * FDPart3tmp54 + FDPart3tmp25 * FDPart3tmp46) +
            FDPart3tmp60 * (FDPart3tmp11 * xx0 + FDPart3tmp20 * FDPart3tmp57 + FDPart3tmp22 * FDPart3tmp59 + FDPart3tmp25 * FDPart3tmp58);
        in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)] =
            xx0 *
            (FDPart3tmp21 * (FDPart3tmp23 * FDPart3tmp62 + FDPart3tmp25 * hDD_dD000 + FDPart3tmp4 * FDPart3tmp61) +
             FDPart3tmp39 * (FDPart3tmp25 * hDD_dD001 + FDPart3tmp26 * FDPart3tmp62 + FDPart3tmp37 * FDPart3tmp61 + FDPart3tmp63) +
             FDPart3tmp44 * (FDPart3tmp25 * hDD_dD002 + FDPart3tmp42 * FDPart3tmp61 + FDPart3tmp43 * FDPart3tmp62) +
             FDPart3tmp50 * (FDPart3tmp25 * FDPart3tmp49 + FDPart3tmp45 * FDPart3tmp62 + FDPart3tmp48 * FDPart3tmp61) +
             FDPart3tmp56 * (FDPart3tmp25 * FDPart3tmp54 + FDPart3tmp46 * FDPart3tmp62 + FDPart3tmp53 * FDPart3tmp61) +
             FDPart3tmp60 * (FDPart3tmp0 * FDPart3tmp28 + FDPart3tmp25 * FDPart3tmp59 + FDPart3tmp57 * FDPart3tmp61 + FDPart3tmp58 * FDPart3tmp62));
        in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)] =
            FDPart3tmp0 * xx0 *
            (FDPart3tmp21 * (FDPart3tmp20 * hDD_dD000 + FDPart3tmp23 * FDPart3tmp61 + FDPart3tmp4 * FDPart3tmp64) +
             FDPart3tmp39 * (FDPart3tmp20 * hDD_dD001 + FDPart3tmp26 * FDPart3tmp61 + FDPart3tmp37 * FDPart3tmp64) +
             FDPart3tmp44 * (FDPart3tmp20 * hDD_dD002 + FDPart3tmp42 * FDPart3tmp64 + FDPart3tmp43 * FDPart3tmp61 + FDPart3tmp63) +
             FDPart3tmp50 * (FDPart3tmp20 * FDPart3tmp49 + FDPart3tmp45 * FDPart3tmp61 + FDPart3tmp48 * FDPart3tmp64) +
             FDPart3tmp56 * (FDPart3tmp20 * FDPart3tmp54 + FDPart3tmp46 * FDPart3tmp61 + FDPart3tmp53 * FDPart3tmp64 - FDPart3tmp28 / FDPart3tmp0) +
             FDPart3tmp60 * (FDPart3tmp20 * FDPart3tmp59 + FDPart3tmp57 * FDPart3tmp64 + FDPart3tmp58 * FDPart3tmp61));

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  }     // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}

/**
 * Read ADM data in the Spherical basis, and output rescaled BSSN data in the Spherical basis
 */
void initial_data_reader__convert_ADM_Spherical_to_BSSN__rfm__Spherical(
    const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct,
    MoL_gridfunctions_struct *restrict gridfuncs, ID_persist_struct *restrict ID_persist,
    void ID_function(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                     const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data)) {

  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
    // xCart is the global Cartesian coordinate, which accounts for any grid offsets from the origin.
    REAL xCart[3];
    xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

    // Read or compute initial data at destination point xCart
    initial_data_struct initial_data;
    ID_function(commondata, params, xCart, ID_persist, &initial_data);

    ADM_Cart_basis_struct ADM_Cart_basis;
    ADM_SphorCart_to_Cart(commondata, params, xCart, &initial_data, &ADM_Cart_basis);

    BSSN_Cart_basis_struct BSSN_Cart_basis;
    ADM_Cart_to_BSSN_Cart(commondata, params, xCart, &ADM_Cart_basis, &BSSN_Cart_basis);

    rescaled_BSSN_rfm_basis_struct rescaled_BSSN_rfm_basis;
    BSSN_Cart_to_rescaled_BSSN_rfm(commondata, params, xCart, &BSSN_Cart_basis, &rescaled_BSSN_rfm_basis);

    const int idx3 = IDX3(i0, i1, i2);
    gridfuncs->y_n_gfs[IDX4pt(ADD00GF, idx3)] = rescaled_BSSN_rfm_basis.aDD00;
    gridfuncs->y_n_gfs[IDX4pt(ADD01GF, idx3)] = rescaled_BSSN_rfm_basis.aDD01;
    gridfuncs->y_n_gfs[IDX4pt(ADD02GF, idx3)] = rescaled_BSSN_rfm_basis.aDD02;
    gridfuncs->y_n_gfs[IDX4pt(ADD11GF, idx3)] = rescaled_BSSN_rfm_basis.aDD11;
    gridfuncs->y_n_gfs[IDX4pt(ADD12GF, idx3)] = rescaled_BSSN_rfm_basis.aDD12;
    gridfuncs->y_n_gfs[IDX4pt(ADD22GF, idx3)] = rescaled_BSSN_rfm_basis.aDD22;
    gridfuncs->y_n_gfs[IDX4pt(ALPHAGF, idx3)] = rescaled_BSSN_rfm_basis.alpha;
    gridfuncs->y_n_gfs[IDX4pt(BETU0GF, idx3)] = rescaled_BSSN_rfm_basis.betU0;
    gridfuncs->y_n_gfs[IDX4pt(BETU1GF, idx3)] = rescaled_BSSN_rfm_basis.betU1;
    gridfuncs->y_n_gfs[IDX4pt(BETU2GF, idx3)] = rescaled_BSSN_rfm_basis.betU2;
    gridfuncs->y_n_gfs[IDX4pt(CFGF, idx3)] = rescaled_BSSN_rfm_basis.cf;
    gridfuncs->y_n_gfs[IDX4pt(HDD00GF, idx3)] = rescaled_BSSN_rfm_basis.hDD00;
    gridfuncs->y_n_gfs[IDX4pt(HDD01GF, idx3)] = rescaled_BSSN_rfm_basis.hDD01;
    gridfuncs->y_n_gfs[IDX4pt(HDD02GF, idx3)] = rescaled_BSSN_rfm_basis.hDD02;
    gridfuncs->y_n_gfs[IDX4pt(HDD11GF, idx3)] = rescaled_BSSN_rfm_basis.hDD11;
    gridfuncs->y_n_gfs[IDX4pt(HDD12GF, idx3)] = rescaled_BSSN_rfm_basis.hDD12;
    gridfuncs->y_n_gfs[IDX4pt(HDD22GF, idx3)] = rescaled_BSSN_rfm_basis.hDD22;
    gridfuncs->y_n_gfs[IDX4pt(TRKGF, idx3)] = rescaled_BSSN_rfm_basis.trK;
    gridfuncs->y_n_gfs[IDX4pt(VETU0GF, idx3)] = rescaled_BSSN_rfm_basis.vetU0;
    gridfuncs->y_n_gfs[IDX4pt(VETU1GF, idx3)] = rescaled_BSSN_rfm_basis.vetU1;
    gridfuncs->y_n_gfs[IDX4pt(VETU2GF, idx3)] = rescaled_BSSN_rfm_basis.vetU2;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU00GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU00;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU01GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU01;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU02GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU02;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU03GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU03;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU11GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU11;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU12GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU12;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU13GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU13;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU22GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU22;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU23GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU23;
    gridfuncs->auxevol_gfs[IDX4pt(T4UU33GF, idx3)] = rescaled_BSSN_rfm_basis.T4UU33;

    // Initialize lambdaU to zero
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU0GF, idx3)] = 0.0;
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU1GF, idx3)] = 0.0;
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU2GF, idx3)] = 0.0;
  } // END LOOP over all gridpoints on given grid

  // Now we've set all but lambda^i, which will be computed via a finite-difference of hDD.
  //    However, hDD is not correctly set in inner boundary points so we apply inner bcs first.

  // Apply inner bcs to get correct values of all tensor quantities across symmetry boundaries;
  //    BSSN_Cart_to_rescaled_BSSN_rfm() converts each xCart->xx, which guarantees a mapping
  //    to the grid interior. It therefore does not account for parity conditions across
  //    symmetry boundaries being correct.
  apply_bcs_inner_only(commondata, params, bcstruct, gridfuncs->y_n_gfs);

  initial_data_lambdaU_grid_interior(commondata, params, xx, gridfuncs->y_n_gfs);
}
