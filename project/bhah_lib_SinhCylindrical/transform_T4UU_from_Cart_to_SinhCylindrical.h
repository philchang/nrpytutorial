const REAL tmp0 = cos(xx1);
const REAL tmp1 = (1.0/(SINHWRHO));
const REAL tmp8 = sin(xx1);
const REAL tmp9 = ((AMPLRHO)*(AMPLRHO));
const REAL tmp15 = (1.0/(SINHWZ));
const REAL tmp33 = ((AMPLZ)*(AMPLZ));
const REAL tmp10 = exp(tmp1) - exp(-tmp1);
const REAL tmp16 = exp(tmp15) - exp(-tmp15);
const REAL tmp3 = exp(tmp1*xx0);
const REAL tmp4 = exp(-tmp1*xx0);
const REAL tmp11 = (1.0/((tmp10)*(tmp10)));
const REAL tmp17 = (1.0/(tmp16));
const REAL tmp19 = tmp15*exp(tmp15*xx2) + tmp15*exp(-tmp15*xx2);
const REAL tmp26 = (1.0/(tmp10));
const REAL tmp34 = (1.0/((tmp16)*(tmp16)));
const REAL tmp5 = tmp3 - tmp4;
const REAL tmp6 = tmp1*tmp3 + tmp1*tmp4;
const REAL tmp12 = tmp11*tmp9;
const REAL tmp20 = AMPLZ*tmp17*tmp19;
const REAL tmp35 = ((tmp19)*(tmp19));
const REAL tmp7 = tmp5*tmp6;
const REAL tmp13 = tmp12*((tmp8)*(tmp8));
const REAL tmp22 = ((tmp0)*(tmp0))*tmp12;
const REAL tmp27 = AMPLRHO*tmp20*tmp26;
const REAL tmp38 = tmp0*tmp12*tmp8;
const REAL tmp14 = tmp13*tmp7;
const REAL tmp39 = 2*T4CartUU12*tmp38;
const REAL tmp24 = tmp14*tmp20 + tmp20*tmp22*tmp7;
const REAL tmp31 = tmp14 + tmp22*tmp7;
const REAL tmp25 = (1.0/(tmp24));
const REAL tmp32 = (1.0/((tmp24)*(tmp24)));
const REAL tmp29 = tmp25*tmp27*tmp5;
const REAL tmp36 = tmp32*tmp33*tmp34*tmp35;
const REAL tmp40 = tmp27*tmp31*tmp32*tmp5;
const REAL tmp37 = tmp36*((tmp5)*(tmp5));
const REAL tmp41 = tmp36*((tmp6)*(tmp6));
T4UU00 = T4CartUU00;
T4UU01 = T4CartUU01*tmp0*tmp29 + T4CartUU02*tmp29*tmp8;
T4UU02 = AMPLRHO*AMPLZ*T4CartUU02*tmp0*tmp17*tmp19*tmp25*tmp26*tmp6 - T4CartUU01*tmp25*tmp27*tmp6*tmp8;
T4UU03 = T4CartUU03*tmp25*tmp31;
T4UU11 = T4CartUU11*tmp22*tmp37 + T4CartUU22*tmp13*tmp37 + tmp37*tmp39;
T4UU12 = -T4CartUU11*tmp36*tmp38*tmp7 + T4CartUU12*((tmp0)*(tmp0))*tmp11*tmp32*tmp33*tmp34*tmp35*tmp5*tmp6*tmp9 - T4CartUU12*tmp14*tmp36 + T4CartUU22*tmp0*tmp11*tmp32*tmp33*tmp34*tmp35*tmp5*tmp6*tmp8*tmp9;
T4UU13 = T4CartUU13*tmp0*tmp40 + T4CartUU23*tmp40*tmp8;
T4UU22 = T4CartUU11*tmp13*tmp41 + T4CartUU22*tmp22*tmp41 - tmp39*tmp41;
T4UU23 = AMPLRHO*AMPLZ*T4CartUU23*tmp0*tmp17*tmp19*tmp26*tmp31*tmp32*tmp6 - T4CartUU13*tmp27*tmp31*tmp32*tmp6*tmp8;
T4UU33 = T4CartUU33*((tmp31)*(tmp31))*tmp32;
