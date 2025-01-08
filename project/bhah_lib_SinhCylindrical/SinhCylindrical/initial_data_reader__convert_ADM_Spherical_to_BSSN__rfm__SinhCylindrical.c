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
 * Cartesian -> SinhCylindrical basis transformation of BSSN vectors/tensors *except* lambda^i.
 * After the basis transform, all BSSN quantities are rescaled.
 */
static void BSSN_Cart_to_rescaled_BSSN_rfm(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xxL[3],
                                           const BSSN_Cart_basis_struct *restrict BSSN_Cart_basis,
                                           rescaled_BSSN_rfm_basis_struct *restrict rescaled_BSSN_rfm_basis) {
#include "../set_CodeParameters.h"

  const REAL xx0 = xxL[0], xx1 = xxL[1], xx2 = xxL[2];
  const REAL tmp0 = cos(xx1);
  const REAL tmp1 = (1.0 / (SINHWRHO));
  const REAL tmp10 = (1.0 / (SINHWZ));
  const REAL tmp18 = sin(xx1);
  const REAL tmp21 = ((AMPLRHO) * (AMPLRHO));
  const REAL tmp42 = ((AMPLZ) * (AMPLZ));
  const REAL tmp6 = exp(tmp1) - exp(-tmp1);
  const REAL tmp11 = exp(tmp10) - exp(-tmp10);
  const REAL tmp19 = ((tmp18) * (tmp18));
  const REAL tmp27 = ((tmp0) * (tmp0));
  const REAL tmp3 = exp(tmp1 * xx0);
  const REAL tmp4 = exp(-tmp1 * xx0);
  const REAL tmp7 = (1.0 / (tmp6));
  const REAL tmp12 = (1.0 / (tmp11));
  const REAL tmp14 = tmp10 * exp(tmp10 * xx2) + tmp10 * exp(-tmp10 * xx2);
  const REAL tmp52 = 2 * tmp0 * tmp18;
  const REAL tmp5 = tmp3 - tmp4;
  const REAL tmp15 = AMPLZ * tmp12 * tmp14;
  const REAL tmp20 = tmp1 * tmp3 + tmp1 * tmp4;
  const REAL tmp23 = (1.0 / ((tmp6) * (tmp6)));
  const REAL tmp44 = (1.0 / ((tmp11) * (tmp11)));
  const REAL tmp45 = ((tmp14) * (tmp14));
  const REAL tmp63 = ((tmp6) * (tmp6)) / tmp21;
  const REAL tmp69 = tmp11 * tmp6 / (AMPLRHO * AMPLZ * tmp14);
  const REAL tmp9 = AMPLRHO * tmp5 * tmp7;
  const REAL tmp24 = tmp21 * tmp23;
  const REAL tmp33 = AMPLRHO * tmp20 * tmp7;
  const REAL tmp46 = tmp42 * tmp44 * tmp45;
  const REAL tmp65 = (1.0 / (tmp5));
  const REAL tmp66 = (1.0 / (tmp20));
  const REAL tmp25 = tmp20 * tmp24 * tmp5;
  const REAL tmp39 = tmp24 * ((tmp5) * (tmp5));
  const REAL tmp59 = ((tmp20) * (tmp20)) * tmp24;
  const REAL tmp64 = tmp63 / ((tmp20) * (tmp20));
  const REAL tmp67 = tmp63 * tmp65 * tmp66;
  const REAL tmp71 = tmp63 / ((tmp5) * (tmp5));
  const REAL tmp17 = tmp0 * tmp15 * tmp9;
  const REAL tmp26 = tmp19 * tmp25;
  const REAL tmp32 = tmp15 * tmp18 * tmp9;
  const REAL tmp35 = tmp15 * tmp18 * tmp33;
  const REAL tmp40 = tmp27 * tmp39;
  const REAL tmp49 = tmp19 * tmp39;
  const REAL tmp55 = tmp0 * tmp18 * tmp25;
  const REAL tmp60 = tmp19 * tmp59;
  const REAL tmp61 = tmp27 * tmp59;
  const REAL tmp68 = tmp0 * tmp15 * tmp33;
  const REAL tmp29 = tmp15 * tmp25 * tmp27 + tmp15 * tmp26;
  const REAL tmp36 = tmp25 * tmp27 + tmp26;
  const REAL tmp30 = (1.0 / (tmp29));
  const REAL tmp41 = (1.0 / ((tmp29) * (tmp29)));
  const REAL tmp47 = tmp41 * tmp46;
  const REAL tmp48 = BSSN_Cart_basis->T4UU11 * tmp47;
  const REAL tmp57 = BSSN_Cart_basis->T4UU13 * tmp36 * tmp41;
  const REAL tmp74 = tmp15 * tmp30 * tmp36;
  const REAL tmp54 = BSSN_Cart_basis->T4UU12 * tmp47 * tmp52;
  rescaled_BSSN_rfm_basis->T4UU00 = BSSN_Cart_basis->T4UU00;
  rescaled_BSSN_rfm_basis->T4UU01 = BSSN_Cart_basis->T4UU01 * tmp17 * tmp30 + BSSN_Cart_basis->T4UU02 * tmp30 * tmp32;
  rescaled_BSSN_rfm_basis->T4UU02 =
      AMPLRHO * AMPLZ * BSSN_Cart_basis->T4UU02 * tmp0 * tmp12 * tmp14 * tmp20 * tmp30 * tmp7 - BSSN_Cart_basis->T4UU01 * tmp30 * tmp35;
  rescaled_BSSN_rfm_basis->T4UU03 = BSSN_Cart_basis->T4UU03 * tmp30 * tmp36;
  rescaled_BSSN_rfm_basis->T4UU11 = BSSN_Cart_basis->T4UU22 * tmp47 * tmp49 + tmp39 * tmp54 + tmp40 * tmp48;
  rescaled_BSSN_rfm_basis->T4UU12 =
      BSSN_Cart_basis->T4UU12 * tmp20 * tmp21 * tmp23 * tmp27 * tmp41 * tmp42 * tmp44 * tmp45 * tmp5 - BSSN_Cart_basis->T4UU12 * tmp26 * tmp47 +
      BSSN_Cart_basis->T4UU22 * tmp0 * tmp18 * tmp20 * tmp21 * tmp23 * tmp41 * tmp42 * tmp44 * tmp45 * tmp5 - tmp48 * tmp55;
  rescaled_BSSN_rfm_basis->T4UU13 = BSSN_Cart_basis->T4UU23 * tmp32 * tmp36 * tmp41 + tmp17 * tmp57;
  rescaled_BSSN_rfm_basis->T4UU22 = BSSN_Cart_basis->T4UU22 * tmp47 * tmp61 + tmp48 * tmp60 - tmp54 * tmp59;
  rescaled_BSSN_rfm_basis->T4UU23 = AMPLRHO * AMPLZ * BSSN_Cart_basis->T4UU23 * tmp0 * tmp12 * tmp14 * tmp20 * tmp36 * tmp41 * tmp7 - tmp35 * tmp57;
  rescaled_BSSN_rfm_basis->T4UU33 = BSSN_Cart_basis->T4UU33 * ((tmp36) * (tmp36)) * tmp41;
  rescaled_BSSN_rfm_basis->aDD00 =
      tmp64 * (BSSN_Cart_basis->AbarDD00 * tmp61 + BSSN_Cart_basis->AbarDD01 * tmp52 * tmp59 + BSSN_Cart_basis->AbarDD11 * tmp60);
  rescaled_BSSN_rfm_basis->aDD01 =
      tmp67 * (-BSSN_Cart_basis->AbarDD00 * tmp55 + BSSN_Cart_basis->AbarDD01 * tmp20 * tmp21 * tmp23 * tmp27 * tmp5 -
               BSSN_Cart_basis->AbarDD01 * tmp26 + BSSN_Cart_basis->AbarDD11 * tmp0 * tmp18 * tmp20 * tmp21 * tmp23 * tmp5);
  rescaled_BSSN_rfm_basis->aDD02 = tmp66 * tmp69 * (BSSN_Cart_basis->AbarDD02 * tmp68 + BSSN_Cart_basis->AbarDD12 * tmp35);
  rescaled_BSSN_rfm_basis->aDD11 =
      tmp71 * (BSSN_Cart_basis->AbarDD00 * tmp49 - BSSN_Cart_basis->AbarDD01 * tmp39 * tmp52 + BSSN_Cart_basis->AbarDD11 * tmp40);
  rescaled_BSSN_rfm_basis->aDD12 =
      tmp65 * tmp69 * (AMPLRHO * AMPLZ * BSSN_Cart_basis->AbarDD12 * tmp0 * tmp12 * tmp14 * tmp5 * tmp7 - BSSN_Cart_basis->AbarDD02 * tmp32);
  rescaled_BSSN_rfm_basis->aDD22 = BSSN_Cart_basis->AbarDD22;
  rescaled_BSSN_rfm_basis->alpha = BSSN_Cart_basis->alpha;
  rescaled_BSSN_rfm_basis->betU0 = tmp33 * (BSSN_Cart_basis->BU0 * tmp17 * tmp30 + BSSN_Cart_basis->BU1 * tmp30 * tmp32);
  rescaled_BSSN_rfm_basis->betU1 =
      tmp9 * (AMPLRHO * AMPLZ * BSSN_Cart_basis->BU1 * tmp0 * tmp12 * tmp14 * tmp20 * tmp30 * tmp7 - BSSN_Cart_basis->BU0 * tmp30 * tmp35);
  rescaled_BSSN_rfm_basis->betU2 = BSSN_Cart_basis->BU2 * tmp74;
  rescaled_BSSN_rfm_basis->cf = BSSN_Cart_basis->cf;
  rescaled_BSSN_rfm_basis->hDD00 =
      tmp64 * (BSSN_Cart_basis->gammabarDD00 * tmp61 + BSSN_Cart_basis->gammabarDD01 * tmp52 * tmp59 + BSSN_Cart_basis->gammabarDD11 * tmp60 - tmp59);
  rescaled_BSSN_rfm_basis->hDD01 =
      tmp67 * (-BSSN_Cart_basis->gammabarDD00 * tmp55 + BSSN_Cart_basis->gammabarDD01 * tmp20 * tmp21 * tmp23 * tmp27 * tmp5 -
               BSSN_Cart_basis->gammabarDD01 * tmp26 + BSSN_Cart_basis->gammabarDD11 * tmp0 * tmp18 * tmp20 * tmp21 * tmp23 * tmp5);
  rescaled_BSSN_rfm_basis->hDD02 = tmp66 * tmp69 * (BSSN_Cart_basis->gammabarDD02 * tmp68 + BSSN_Cart_basis->gammabarDD12 * tmp35);
  rescaled_BSSN_rfm_basis->hDD11 =
      tmp71 * (BSSN_Cart_basis->gammabarDD00 * tmp49 - BSSN_Cart_basis->gammabarDD01 * tmp39 * tmp52 + BSSN_Cart_basis->gammabarDD11 * tmp40 - tmp39);
  rescaled_BSSN_rfm_basis->hDD12 =
      tmp65 * tmp69 * (AMPLRHO * AMPLZ * BSSN_Cart_basis->gammabarDD12 * tmp0 * tmp12 * tmp14 * tmp5 * tmp7 - BSSN_Cart_basis->gammabarDD02 * tmp32);
  rescaled_BSSN_rfm_basis->hDD22 = ((tmp11) * (tmp11)) * (BSSN_Cart_basis->gammabarDD22 * tmp46 - tmp46) / (tmp42 * tmp45);
  rescaled_BSSN_rfm_basis->trK = BSSN_Cart_basis->trK;
  rescaled_BSSN_rfm_basis->vetU0 = tmp33 * (BSSN_Cart_basis->betaU0 * tmp17 * tmp30 + BSSN_Cart_basis->betaU1 * tmp30 * tmp32);
  rescaled_BSSN_rfm_basis->vetU1 =
      tmp9 * (AMPLRHO * AMPLZ * BSSN_Cart_basis->betaU1 * tmp0 * tmp12 * tmp14 * tmp20 * tmp30 * tmp7 - BSSN_Cart_basis->betaU0 * tmp30 * tmp35);
  rescaled_BSSN_rfm_basis->vetU2 = BSSN_Cart_basis->betaU2 * tmp74;
}
/**
 * Compute lambdaU in SinhCylindrical coordinates
 */
static void initial_data_lambdaU_grid_interior(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                               REAL *restrict xx[3], REAL *restrict in_gfs) {
#include "../set_CodeParameters.h"
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    MAYBE_UNUSED const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
      MAYBE_UNUSED const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
        MAYBE_UNUSED const REAL xx0 = xx[0][i0]; /*
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
        static const REAL FDPart1_Rational_2_3 = 2.0 / 3.0;
        static const REAL FDPart1_Rational_1_12 = 1.0 / 12.0;
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
        const REAL FDPart3tmp0 = ((AMPLRHO) * (AMPLRHO) * (AMPLRHO) * (AMPLRHO));
        const REAL FDPart3tmp2 = (1.0 / (SINHWRHO));
        const REAL FDPart3tmp5 = (1.0 / (SINHWZ));
        const REAL FDPart3tmp82 = (1.0 / ((SINHWZ) * (SINHWZ)));
        const REAL FDPart3tmp94 = (1.0 / ((SINHWRHO) * (SINHWRHO)));
        const REAL FDPart3tmp3 = exp(FDPart3tmp2) - exp(-FDPart3tmp2);
        const REAL FDPart3tmp6 = exp(FDPart3tmp5) - exp(-FDPart3tmp5);
        const REAL FDPart3tmp4 = (1.0 / ((FDPart3tmp3) * (FDPart3tmp3) * (FDPart3tmp3) * (FDPart3tmp3)));
        const REAL FDPart3tmp7 = (1.0 / ((FDPart3tmp6) * (FDPart3tmp6)));
        const REAL FDPart3tmp9 = exp(FDPart3tmp2 * xx0);
        const REAL FDPart3tmp10 = exp(-FDPart3tmp2 * xx0);
        const REAL FDPart3tmp18 = exp(FDPart3tmp5 * xx2);
        const REAL FDPart3tmp19 = exp(-FDPart3tmp5 * xx2);
        const REAL FDPart3tmp27 = (1.0 / ((FDPart3tmp3) * (FDPart3tmp3)));
        const REAL FDPart3tmp40 = (1.0 / (FDPart3tmp6));
        const REAL FDPart3tmp44 = (1.0 / (FDPart3tmp3));
        const REAL FDPart3tmp11 = -FDPart3tmp10 + FDPart3tmp9;
        const REAL FDPart3tmp20 = FDPart3tmp18 * FDPart3tmp5 + FDPart3tmp19 * FDPart3tmp5;
        const REAL FDPart3tmp22 = ((AMPLZ) * (AMPLZ)) * FDPart3tmp7;
        const REAL FDPart3tmp28 = ((AMPLRHO) * (AMPLRHO)) * FDPart3tmp27;
        const REAL FDPart3tmp45 = AMPLRHO * FDPart3tmp44;
        const REAL FDPart3tmp12 = ((FDPart3tmp11) * (FDPart3tmp11));
        const REAL FDPart3tmp15 = FDPart3tmp10 * FDPart3tmp2 + FDPart3tmp2 * FDPart3tmp9;
        const REAL FDPart3tmp42 = AMPLZ * FDPart3tmp20 * FDPart3tmp40;
        const REAL FDPart3tmp56 = FDPart3tmp11 * FDPart3tmp45;
        const REAL FDPart3tmp70 = 2 * FDPart3tmp10 * FDPart3tmp2 + 2 * FDPart3tmp2 * FDPart3tmp9;
        const REAL FDPart3tmp85 = 2 * FDPart3tmp18 * FDPart3tmp82 - 2 * FDPart3tmp19 * FDPart3tmp82;
        const REAL FDPart3tmp88 = AMPLZ * FDPart3tmp40 * (FDPart3tmp18 * FDPart3tmp82 - FDPart3tmp19 * FDPart3tmp82);
        const REAL FDPart3tmp97 = -2 * FDPart3tmp10 * FDPart3tmp94 + 2 * FDPart3tmp9 * FDPart3tmp94;
        const REAL FDPart3tmp100 = -FDPart3tmp10 * FDPart3tmp94 + FDPart3tmp9 * FDPart3tmp94;
        const REAL FDPart3tmp16 = ((FDPart3tmp15) * (FDPart3tmp15));
        const REAL FDPart3tmp23 = ((FDPart3tmp20) * (FDPart3tmp20)) * FDPart3tmp22;
        const REAL FDPart3tmp29 = FDPart3tmp12 * FDPart3tmp28;
        const REAL FDPart3tmp43 = ((AMPLRHO) * (AMPLRHO) * (AMPLRHO)) * FDPart3tmp42 * hDD01 / ((FDPart3tmp3) * (FDPart3tmp3) * (FDPart3tmp3));
        const REAL FDPart3tmp46 = FDPart3tmp15 * FDPart3tmp45;
        const REAL FDPart3tmp57 = FDPart3tmp42 * FDPart3tmp56;
        const REAL FDPart3tmp75 = FDPart3tmp11 * FDPart3tmp28 * FDPart3tmp70;
        const REAL FDPart3tmp86 = FDPart3tmp20 * FDPart3tmp22 * FDPart3tmp85;
        const REAL FDPart3tmp24 = FDPart3tmp23 * hDD22 + FDPart3tmp23;
        const REAL FDPart3tmp25 = FDPart3tmp0 * FDPart3tmp12 * FDPart3tmp16 * FDPart3tmp4 * ((hDD01) * (hDD01));
        const REAL FDPart3tmp30 = FDPart3tmp29 * hDD11 + FDPart3tmp29;
        const REAL FDPart3tmp31 = FDPart3tmp16 * FDPart3tmp28;
        const REAL FDPart3tmp34 = FDPart3tmp23 * FDPart3tmp29 * ((hDD12) * (hDD12));
        const REAL FDPart3tmp39 = FDPart3tmp23 * hDD_dD220;
        const REAL FDPart3tmp47 = FDPart3tmp42 * FDPart3tmp46;
        const REAL FDPart3tmp52 = FDPart3tmp11 * FDPart3tmp15 * FDPart3tmp28;
        const REAL FDPart3tmp65 = FDPart3tmp29 * hDD_dD112;
        const REAL FDPart3tmp66 = FDPart3tmp23 * hDD_dD221;
        const REAL FDPart3tmp72 = FDPart3tmp29 * hDD_dD111;
        const REAL FDPart3tmp76 = FDPart3tmp29 * hDD_dD110 + FDPart3tmp75 * hDD11 + FDPart3tmp75;
        const REAL FDPart3tmp87 = FDPart3tmp23 * hDD_dD222 + FDPart3tmp86 * hDD22 + FDPart3tmp86;
        const REAL FDPart3tmp98 = FDPart3tmp15 * FDPart3tmp28 * FDPart3tmp97;
        const REAL FDPart3tmp32 = FDPart3tmp23 * FDPart3tmp31 * ((hDD02) * (hDD02));
        const REAL FDPart3tmp33 = FDPart3tmp31 * hDD00 + FDPart3tmp31;
        const REAL FDPart3tmp38 = FDPart3tmp31 * hDD_dD002;
        const REAL FDPart3tmp48 = FDPart3tmp12 * FDPart3tmp15 * FDPart3tmp43 * hDD12 - FDPart3tmp30 * FDPart3tmp47 * hDD02;
        const REAL FDPart3tmp53 = FDPart3tmp52 * hDD_dD012;
        const REAL FDPart3tmp54 = FDPart3tmp47 * hDD_dD021;
        const REAL FDPart3tmp61 = FDPart3tmp23 * FDPart3tmp52 * hDD02 * hDD12 - FDPart3tmp24 * FDPart3tmp52 * hDD01;
        const REAL FDPart3tmp73 = 2 * AMPLRHO * AMPLZ * FDPart3tmp11 * FDPart3tmp20 * FDPart3tmp40 * FDPart3tmp44 * hDD_dD121 - FDPart3tmp65;
        const REAL FDPart3tmp77 = 2 * ((AMPLRHO) * (AMPLRHO)) * FDPart3tmp11 * FDPart3tmp15 * FDPart3tmp27 * hDD_dD011 - FDPart3tmp76;
        const REAL FDPart3tmp79 = FDPart3tmp31 * hDD_dD001;
        const REAL FDPart3tmp89 = 2 * FDPart3tmp56 * FDPart3tmp88 * hDD12 + 2 * FDPart3tmp57 * hDD_dD122 - FDPart3tmp66;
        const REAL FDPart3tmp99 = FDPart3tmp31 * hDD_dD000 + FDPart3tmp98 * hDD00 + FDPart3tmp98;
        const REAL FDPart3tmp35 = (1.0 / (2 * ((AMPLZ) * (AMPLZ)) * FDPart3tmp0 * FDPart3tmp12 * FDPart3tmp16 * ((FDPart3tmp20) * (FDPart3tmp20)) *
                                              FDPart3tmp4 * FDPart3tmp7 * hDD01 * hDD02 * hDD12 -
                                          FDPart3tmp24 * FDPart3tmp25 + FDPart3tmp24 * FDPart3tmp30 * FDPart3tmp33 - FDPart3tmp30 * FDPart3tmp32 -
                                          FDPart3tmp33 * FDPart3tmp34));
        const REAL FDPart3tmp59 = FDPart3tmp47 * hDD12 + FDPart3tmp57 * hDD_dD120;
        const REAL FDPart3tmp67 = -FDPart3tmp47 * hDD12 + FDPart3tmp53 + FDPart3tmp54 - FDPart3tmp57 * hDD_dD120;
        const REAL FDPart3tmp68 = FDPart3tmp11 * FDPart3tmp16 * FDPart3tmp43 * hDD02 - FDPart3tmp33 * FDPart3tmp57 * hDD12;
        const REAL FDPart3tmp92 = -FDPart3tmp39 + 2 * FDPart3tmp46 * FDPart3tmp88 * hDD02 + 2 * FDPart3tmp47 * hDD_dD022;
        const REAL FDPart3tmp101 =
            2 * FDPart3tmp52 * hDD_dD010 - FDPart3tmp79 + 2 * hDD01 * (FDPart3tmp100 * FDPart3tmp11 * FDPart3tmp28 + FDPart3tmp31);
        const REAL FDPart3tmp102 = 2 * FDPart3tmp100 * FDPart3tmp42 * FDPart3tmp45 * hDD02 - FDPart3tmp38 + 2 * FDPart3tmp47 * hDD_dD020;
        const REAL FDPart3tmp36 = FDPart3tmp35 * (FDPart3tmp24 * FDPart3tmp30 - FDPart3tmp34);
        const REAL FDPart3tmp49 = (1.0 / 2.0) * FDPart3tmp35;
        const REAL FDPart3tmp60 = FDPart3tmp53 - FDPart3tmp54 + FDPart3tmp59;
        const REAL FDPart3tmp63 = 2 * FDPart3tmp35;
        const REAL FDPart3tmp78 = FDPart3tmp35 * (FDPart3tmp24 * FDPart3tmp33 - FDPart3tmp32);
        const REAL FDPart3tmp80 = -FDPart3tmp53 + FDPart3tmp54 + FDPart3tmp59;
        const REAL FDPart3tmp93 = FDPart3tmp35 * (-FDPart3tmp25 + FDPart3tmp30 * FDPart3tmp33);
        const REAL FDPart3tmp37 = (1.0 / 2.0) * FDPart3tmp36;
        const REAL FDPart3tmp50 = FDPart3tmp48 * FDPart3tmp49;
        const REAL FDPart3tmp62 = FDPart3tmp49 * FDPart3tmp61;
        const REAL FDPart3tmp64 = FDPart3tmp48 * FDPart3tmp63;
        const REAL FDPart3tmp69 = FDPart3tmp63 * FDPart3tmp68;
        const REAL FDPart3tmp81 = FDPart3tmp61 * FDPart3tmp63;
        const REAL FDPart3tmp103 = FDPart3tmp49 * FDPart3tmp68;
        const REAL FDPart3tmp104 = (1.0 / 2.0) * FDPart3tmp78;
        const REAL FDPart3tmp105 = (1.0 / 2.0) * FDPart3tmp93;
        in_gfs[IDX4(LAMBDAU0GF, i0, i1, i2)] =
            FDPart3tmp46 * (FDPart3tmp36 * (FDPart3tmp101 * FDPart3tmp62 + FDPart3tmp102 * FDPart3tmp50 + FDPart3tmp37 * FDPart3tmp99 -
                                            1.0 / 2.0 * FDPart3tmp97 / FDPart3tmp15) +
                            FDPart3tmp64 * (FDPart3tmp37 * FDPart3tmp38 + FDPart3tmp39 * FDPart3tmp50 + FDPart3tmp60 * FDPart3tmp62) +
                            FDPart3tmp69 * (FDPart3tmp37 * FDPart3tmp67 + FDPart3tmp50 * FDPart3tmp66 + FDPart3tmp62 * FDPart3tmp65) +
                            FDPart3tmp78 * ((1.0 / 2.0) * FDPart3tmp11 * FDPart3tmp70 / FDPart3tmp16 + FDPart3tmp37 * FDPart3tmp77 +
                                            FDPart3tmp50 * FDPart3tmp73 + FDPart3tmp62 * FDPart3tmp72) +
                            FDPart3tmp81 * (FDPart3tmp37 * FDPart3tmp79 + FDPart3tmp50 * FDPart3tmp80 + FDPart3tmp62 * FDPart3tmp76) +
                            FDPart3tmp93 * (FDPart3tmp37 * FDPart3tmp92 + FDPart3tmp50 * FDPart3tmp87 + FDPart3tmp62 * FDPart3tmp89));
        in_gfs[IDX4(LAMBDAU1GF, i0, i1, i2)] =
            FDPart3tmp56 * (FDPart3tmp36 * (FDPart3tmp101 * FDPart3tmp104 + FDPart3tmp102 * FDPart3tmp103 + FDPart3tmp62 * FDPart3tmp99) +
                            FDPart3tmp64 * (FDPart3tmp103 * FDPart3tmp39 + FDPart3tmp104 * FDPart3tmp60 + FDPart3tmp38 * FDPart3tmp62) +
                            FDPart3tmp69 * (FDPart3tmp103 * FDPart3tmp66 + FDPart3tmp104 * FDPart3tmp65 + FDPart3tmp62 * FDPart3tmp67) +
                            FDPart3tmp78 * (FDPart3tmp103 * FDPart3tmp73 + FDPart3tmp104 * FDPart3tmp72 + FDPart3tmp62 * FDPart3tmp77) +
                            FDPart3tmp81 * (FDPart3tmp103 * FDPart3tmp80 + FDPart3tmp104 * FDPart3tmp76 + FDPart3tmp62 * FDPart3tmp79 -
                                            1.0 / 2.0 * FDPart3tmp70 / FDPart3tmp11) +
                            FDPart3tmp93 * (FDPart3tmp103 * FDPart3tmp87 + FDPart3tmp104 * FDPart3tmp89 + FDPart3tmp62 * FDPart3tmp92));
        in_gfs[IDX4(LAMBDAU2GF, i0, i1, i2)] =
            FDPart3tmp42 * (FDPart3tmp36 * (FDPart3tmp101 * FDPart3tmp103 + FDPart3tmp102 * FDPart3tmp105 + FDPart3tmp50 * FDPart3tmp99) +
                            FDPart3tmp64 * (FDPart3tmp103 * FDPart3tmp60 + FDPart3tmp105 * FDPart3tmp39 + FDPart3tmp38 * FDPart3tmp50) +
                            FDPart3tmp69 * (FDPart3tmp103 * FDPart3tmp65 + FDPart3tmp105 * FDPart3tmp66 + FDPart3tmp50 * FDPart3tmp67) +
                            FDPart3tmp78 * (FDPart3tmp103 * FDPart3tmp72 + FDPart3tmp105 * FDPart3tmp73 + FDPart3tmp50 * FDPart3tmp77) +
                            FDPart3tmp81 * (FDPart3tmp103 * FDPart3tmp76 + FDPart3tmp105 * FDPart3tmp80 + FDPart3tmp50 * FDPart3tmp79) +
                            FDPart3tmp93 * (FDPart3tmp103 * FDPart3tmp89 + FDPart3tmp105 * FDPart3tmp87 + FDPart3tmp50 * FDPart3tmp92 -
                                            1.0 / 2.0 * FDPart3tmp85 / FDPart3tmp20));

      } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
}

/**
 * Read ADM data in the Spherical basis, and output rescaled BSSN data in the SinhCylindrical basis
 */
void initial_data_reader__convert_ADM_Spherical_to_BSSN__rfm__SinhCylindrical(
    const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct,
    MoL_gridfunctions_struct *restrict gridfuncs, ID_persist_struct *restrict ID_persist,
    void ID_function(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                     const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data)) {

  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
    // xxL are the local coordinates on the destination grid
    const REAL xxL[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};

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
    BSSN_Cart_to_rescaled_BSSN_rfm(commondata, params, xxL, &BSSN_Cart_basis, &rescaled_BSSN_rfm_basis);

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
