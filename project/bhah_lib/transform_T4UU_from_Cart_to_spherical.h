const REAL tmp0 = ((xx0)*(xx0));
const REAL tmp1 = sin(xx1);
const REAL tmp4 = cos(xx2);
const REAL tmp5 = sin(xx2);
const REAL tmp10 = cos(xx1);
const REAL tmp6 = ((tmp5)*(tmp5));
const REAL tmp7 = ((tmp1)*(tmp1)*(tmp1));
const REAL tmp9 = ((tmp4)*(tmp4));
const REAL tmp11 = ((tmp10)*(tmp10));
const REAL tmp23 = tmp1*tmp10*xx0;
const REAL tmp29 = ((tmp1)*(tmp1)*(tmp1)*(tmp1))*((xx0)*(xx0)*(xx0)*(xx0));
const REAL tmp40 = T4CartUU13*tmp4;
const REAL tmp43 = T4CartUU23*tmp5;
const REAL tmp3 = tmp0*((tmp1)*(tmp1));
const REAL tmp13 = tmp0*tmp1*tmp11;
const REAL tmp20 = tmp0*tmp1*tmp10;
const REAL tmp24 = ((tmp1)*(tmp1))*xx0;
const REAL tmp45 = tmp10*tmp7*((xx0)*(xx0)*(xx0));
const REAL tmp14 = tmp0*tmp6*tmp7 + tmp0*tmp7*tmp9 + tmp13*tmp6 + tmp13*tmp9;
const REAL tmp21 = tmp20*tmp6 + tmp20*tmp9;
const REAL tmp25 = -tmp24*tmp6 - tmp24*tmp9;
const REAL tmp27 = tmp11*tmp4*xx0 + tmp24*tmp4;
const REAL tmp28 = -tmp11*tmp5*xx0 - tmp24*tmp5;
const REAL tmp49 = tmp3*tmp4;
const REAL tmp15 = (1.0/(tmp14));
const REAL tmp30 = (1.0/((tmp14)*(tmp14)));
const REAL tmp31 = T4CartUU11*tmp30;
const REAL tmp33 = T4CartUU22*tmp30;
const REAL tmp38 = T4CartUU33*tmp30;
const REAL tmp39 = tmp21*tmp30;
const REAL tmp46 = tmp25*tmp30;
const REAL tmp48 = T4CartUU12*tmp27*tmp30;
const REAL tmp17 = T4CartUU01*tmp15*tmp4;
const REAL tmp19 = T4CartUU02*tmp15*tmp5;
const REAL tmp32 = tmp31*tmp9;
const REAL tmp34 = tmp33*tmp6;
const REAL tmp36 = 2*T4CartUU12*tmp30*tmp5;
const REAL tmp53 = T4CartUU12*tmp28*tmp30;
const REAL tmp57 = 2*tmp23*tmp46;
T4UU00 = T4CartUU00;
T4UU01 = T4CartUU03*tmp15*tmp21 + tmp17*tmp3 + tmp19*tmp3;
T4UU02 = T4CartUU03*tmp15*tmp25 + tmp17*tmp23 + tmp19*tmp23;
T4UU03 = T4CartUU01*tmp15*tmp28 + T4CartUU02*tmp15*tmp27;
T4UU11 = ((tmp21)*(tmp21))*tmp38 + tmp29*tmp32 + tmp29*tmp34 + tmp29*tmp36*tmp4 + 2*tmp3*tmp39*tmp40 + 2*tmp3*tmp39*tmp43;
T4UU12 = tmp21*tmp25*tmp38 + tmp23*tmp39*tmp40 + tmp23*tmp39*tmp43 + tmp3*tmp40*tmp46 + tmp3*tmp43*tmp46 + tmp32*tmp45 + tmp34*tmp45 + tmp36*tmp4*tmp45;
T4UU13 = T4CartUU13*tmp28*tmp39 + T4CartUU23*tmp27*tmp39 + tmp27*tmp3*tmp33*tmp5 + tmp28*tmp31*tmp49 + tmp3*tmp5*tmp53 + tmp48*tmp49;
T4UU22 = tmp11*tmp3*tmp32 + tmp11*tmp3*tmp34 + tmp11*tmp36*tmp49 + ((tmp25)*(tmp25))*tmp38 + tmp40*tmp57 + tmp43*tmp57;
T4UU23 = T4CartUU13*tmp28*tmp46 + T4CartUU23*tmp27*tmp46 + tmp23*tmp27*tmp33*tmp5 + tmp23*tmp28*tmp31*tmp4 + tmp23*tmp4*tmp48 + tmp23*tmp5*tmp53;
T4UU33 = ((tmp27)*(tmp27))*tmp33 + ((tmp28)*(tmp28))*tmp31 + 2*tmp28*tmp48;
