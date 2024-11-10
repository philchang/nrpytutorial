REAL xx0, xx1, xx2;
/*
 *  Original SymPy expressions:
 *  "[xx0 = sqrt(Cartx**2 + Carty**2 + Cartz**2),
 *    xx1 = acos(Cartz/sqrt(Cartx**2 + Carty**2 + Cartz**2)),
 *    xx2 = atan2(Carty, Cartx)]"
 */
{
  const double tmp_0 = sqrt(((Cartx)*(Cartx)) + ((Carty)*(Carty)) + ((Cartz)*(Cartz)));
  xx0 = tmp_0;
  xx1 = acos(Cartz/tmp_0);
  xx2 = atan2(Carty, Cartx);
}
/*
 *  Original SymPy expressions:
 *  "[gammaDD[0][0] = 2*hDD01*xx0**4*sin(xx1)**3*cos(xx1)*cos(xx2)**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + 2*hDD02*xx0**3*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*sin(xx1)**3*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + 2*hDD12*xx0**3*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*sin(xx1)**2*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**4*(hDD00 + 1)*sin(xx1)**4*cos(xx2)**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**2*(hDD11*xx0**2 + xx0**2)*sin(xx1)**2*cos(xx1)**2*cos(xx2)**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + (hDD22*xx0**2*sin(xx1)**2 + xx0**2*sin(xx1)**2)*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    gammaDD[0][1] = 2*hDD01*xx0**4*sin(xx1)**3*sin(xx2)*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD02*xx0**3*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*sin(xx1)**3*sin(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD02*xx0**3*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*sin(xx1)**3*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD12*xx0**3*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*sin(xx1)**2*sin(xx2)*cos(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD12*xx0**3*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*sin(xx1)**2*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**4*(hDD00 + 1)*sin(xx1)**4*sin(xx2)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**2*(hDD11*xx0**2 + xx0**2)*sin(xx1)**2*sin(xx2)*cos(xx1)**2*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + (hDD22*xx0**2*sin(xx1)**2 + xx0**2*sin(xx1)**2)*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    gammaDD[0][2] = hDD01*xx0**3*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)**2*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD01*xx0**2*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD02*xx0*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD12*xx0**2*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**2*(hDD00 + 1)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)**2*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0*(hDD11*xx0**2 + xx0**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    gammaDD[1][0] = 2*hDD01*xx0**4*sin(xx1)**3*sin(xx2)*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD02*xx0**3*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*sin(xx1)**3*sin(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD02*xx0**3*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*sin(xx1)**3*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD12*xx0**3*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*sin(xx1)**2*sin(xx2)*cos(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD12*xx0**3*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*sin(xx1)**2*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**4*(hDD00 + 1)*sin(xx1)**4*sin(xx2)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**2*(hDD11*xx0**2 + xx0**2)*sin(xx1)**2*sin(xx2)*cos(xx1)**2*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + (hDD22*xx0**2*sin(xx1)**2 + xx0**2*sin(xx1)**2)*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    gammaDD[1][1] = 2*hDD01*xx0**4*sin(xx1)**3*sin(xx2)**2*cos(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + 2*hDD02*xx0**3*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*sin(xx1)**3*sin(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + 2*hDD12*xx0**3*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*sin(xx1)**2*sin(xx2)*cos(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**4*(hDD00 + 1)*sin(xx1)**4*sin(xx2)**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**2*(hDD11*xx0**2 + xx0**2)*sin(xx1)**2*sin(xx2)**2*cos(xx1)**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + (hDD22*xx0**2*sin(xx1)**2 + xx0**2*sin(xx1)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    gammaDD[1][2] = hDD01*xx0**3*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)**2*sin(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD01*xx0**2*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)*sin(xx2)*cos(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD02*xx0*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD12*xx0**2*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*sin(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**2*(hDD00 + 1)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)**2*sin(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0*(hDD11*xx0**2 + xx0**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)*sin(xx2)*cos(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    gammaDD[2][0] = hDD01*xx0**3*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)**2*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD01*xx0**2*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD02*xx0*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD12*xx0**2*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**2*(hDD00 + 1)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)**2*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0*(hDD11*xx0**2 + xx0**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)*cos(xx1)*cos(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    gammaDD[2][1] = hDD01*xx0**3*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)**2*sin(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD01*xx0**2*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)*sin(xx2)*cos(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD02*xx0*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + hDD12*xx0**2*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*sin(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0**2*(hDD00 + 1)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)**2*sin(xx2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + xx0*(hDD11*xx0**2 + xx0**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)*sin(xx2)*cos(xx1)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    gammaDD[2][2] = 2*hDD01*xx0*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + (hDD00 + 1)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2) + (hDD11*xx0**2 + xx0**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)**2/(cf**2*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2),
 *    KDD[0][0] = xx0**4*(aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*sin(xx1)**4*cos(xx2)**2/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*xx0**3*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*sin(xx1)**3*cos(xx1)*cos(xx2)**2/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*sin(xx1)**2*cos(xx1)**2*cos(xx2)**2/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*xx0**2*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*sin(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*xx0*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))*sin(xx1)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)**2*(aDD22*xx0**2*sin(xx1)**2/cf**2 + trK*(hDD22*xx0**2*sin(xx1)**2 + xx0**2*sin(xx1)**2)/(3*cf**2))/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    KDD[0][1] = xx0**4*(aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*sin(xx1)**4*sin(xx2)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*xx0**3*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*sin(xx1)**3*sin(xx2)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*sin(xx1)**2*sin(xx2)*cos(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*sin(xx1)**2*sin(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*sin(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))*sin(xx1)*sin(xx2)*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))*sin(xx1)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD22*xx0**2*sin(xx1)**2/cf**2 + trK*(hDD22*xx0**2*sin(xx1)**2 + xx0**2*sin(xx1)**2)/(3*cf**2))/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    KDD[0][2] = xx0**2*(aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    KDD[1][0] = xx0**4*(aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*sin(xx1)**4*sin(xx2)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*xx0**3*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*sin(xx1)**3*sin(xx2)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*sin(xx1)**2*sin(xx2)*cos(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*sin(xx1)**2*sin(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*sin(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))*sin(xx1)*sin(xx2)*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))*sin(xx1)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD22*xx0**2*sin(xx1)**2/cf**2 + trK*(hDD22*xx0**2*sin(xx1)**2 + xx0**2*sin(xx1)**2)/(3*cf**2))/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    KDD[1][1] = xx0**4*(aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*sin(xx1)**4*sin(xx2)**2/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*xx0**3*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*sin(xx1)**3*sin(xx2)**2*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*sin(xx1)**2*sin(xx2)**2*cos(xx1)**2/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*xx0**2*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*sin(xx1)**2*sin(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*xx0*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))*sin(xx1)*sin(xx2)*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))**2*(aDD22*xx0**2*sin(xx1)**2/cf**2 + trK*(hDD22*xx0**2*sin(xx1)**2 + xx0**2*sin(xx1)**2)/(3*cf**2))/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    KDD[1][2] = xx0**2*(aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)**2*sin(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)**2*sin(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)*sin(xx2)*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)*sin(xx2)*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    KDD[2][0] = xx0**2*(aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    KDD[2][1] = xx0**2*(aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)**2*sin(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0**2*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)**2*sin(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)*sin(xx1)*sin(xx2)*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + xx0*(aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*sin(xx1)*sin(xx2)*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD12*xx0**2*sin(xx1)/cf**2 + hDD12*trK*xx0**2*sin(xx1)/(3*cf**2))/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))*(aDD02*xx0*sin(xx1)/cf**2 + hDD02*trK*xx0*sin(xx1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    KDD[2][2] = (aDD00/cf**2 + trK*(hDD00 + 1)/(3*cf**2))*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)**2/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + 2*(aDD01*xx0/cf**2 + hDD01*trK*xx0/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2 + (aDD11*xx0**2/cf**2 + trK*(hDD11*xx0**2 + xx0**2)/(3*cf**2))*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)**2/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)**2,
 *    *beta0 = vetU0*xx0**2*sin(xx1)**2*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2) + vetU1*xx0*sin(xx1)**2*sin(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2) + vetU2*(xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1) + xx0**2*sin(xx1)*cos(xx1)*cos(xx2)**2)/(xx0*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)*sin(xx1)),
 *    *beta1 = vetU0*xx0*sin(xx1)*cos(xx1)*cos(xx2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2) + vetU1*sin(xx1)*sin(xx2)*cos(xx1)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2) + vetU2*(-xx0*sin(xx1)**2*sin(xx2)**2 - xx0*sin(xx1)**2*cos(xx2)**2)/(xx0*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2)*sin(xx1)),
 *    *beta2 = vetU0*(-xx0*sin(xx1)**2*sin(xx2) - xx0*sin(xx2)*cos(xx1)**2)/(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2) + vetU1*(xx0*sin(xx1)**2*cos(xx2) + xx0*cos(xx1)**2*cos(xx2))/(xx0*(xx0**2*sin(xx1)**3*sin(xx2)**2 + xx0**2*sin(xx1)**3*cos(xx2)**2 + xx0**2*sin(xx1)*sin(xx2)**2*cos(xx1)**2 + xx0**2*sin(xx1)*cos(xx1)**2*cos(xx2)**2))]"
 */
{
  const double tmp_0 = cos(xx2);
  const double tmp_1 = ((tmp_0)*(tmp_0));
  const double tmp_2 = ((xx0)*(xx0)*(xx0)*(xx0));
  const double tmp_3 = sin(xx1);
  const double tmp_4 = tmp_2*((tmp_3)*(tmp_3)*(tmp_3)*(tmp_3));
  const double tmp_5 = sin(xx2);
  const double tmp_6 = ((tmp_5)*(tmp_5));
  const double tmp_7 = ((xx0)*(xx0));
  const double tmp_8 = ((tmp_3)*(tmp_3)*(tmp_3));
  const double tmp_10 = cos(xx1);
  const double tmp_11 = ((tmp_10)*(tmp_10));
  const double tmp_12 = tmp_3*tmp_7;
  const double tmp_14 = tmp_1*tmp_11;
  const double tmp_15 = tmp_1*tmp_7*tmp_8 + tmp_11*tmp_12*tmp_6 + tmp_12*tmp_14 + tmp_6*tmp_7*tmp_8;
  const double tmp_16 = (1.0/((tmp_15)*(tmp_15)));
  const double tmp_17 = (1.0/((cf)*(cf)));
  const double tmp_18 = tmp_17*(hDD00 + 1);
  const double tmp_20 = tmp_16*tmp_18*tmp_4;
  const double tmp_23 = 2*tmp_1*tmp_10*tmp_8;
  const double tmp_25 = hDD01*tmp_16*tmp_17*tmp_2;
  const double tmp_26 = tmp_17*(hDD11*tmp_7 + tmp_7);
  const double tmp_27 = ((tmp_3)*(tmp_3));
  const double tmp_28 = tmp_27*tmp_7;
  const double tmp_30 = tmp_16*tmp_26*tmp_28;
  const double tmp_31 = ((xx0)*(xx0)*(xx0));
  const double tmp_33 = tmp_0*tmp_17*tmp_31;
  const double tmp_34 = tmp_5*xx0;
  const double tmp_36 = -tmp_11*tmp_34 - tmp_27*tmp_34;
  const double tmp_37 = tmp_16*tmp_36;
  const double tmp_38 = hDD02*tmp_8;
  const double tmp_41 = 2*tmp_37;
  const double tmp_42 = hDD12*tmp_10;
  const double tmp_44 = tmp_17*(hDD22*tmp_28 + tmp_28);
  const double tmp_45 = tmp_16*tmp_44;
  const double tmp_46 = tmp_0*tmp_5;
  const double tmp_47 = 2*tmp_10*tmp_8;
  const double tmp_49 = tmp_0*tmp_28;
  const double tmp_51 = tmp_11*tmp_16*tmp_49*tmp_5;
  const double tmp_52 = tmp_0*xx0;
  const double tmp_53 = tmp_11*tmp_52 + tmp_27*tmp_52;
  const double tmp_54 = tmp_16*tmp_53;
  const double tmp_56 = tmp_17*tmp_31*tmp_5;
  const double tmp_58 = tmp_27*tmp_42*tmp_56;
  const double tmp_60 = tmp_20*tmp_46 + tmp_25*tmp_46*tmp_47 + tmp_26*tmp_51 + tmp_27*tmp_33*tmp_42*tmp_54 + tmp_33*tmp_38*tmp_54 + tmp_36*tmp_45*tmp_53 + tmp_37*tmp_38*tmp_56 + tmp_37*tmp_58;
  const double tmp_62 = -tmp_1*tmp_27*xx0 - tmp_27*tmp_6*xx0;
  const double tmp_63 = tmp_16*tmp_62;
  const double tmp_65 = tmp_1*tmp_10*tmp_12 + tmp_10*tmp_12*tmp_6;
  const double tmp_67 = tmp_17*tmp_3*tmp_7;
  const double tmp_68 = tmp_10*tmp_65*tmp_67;
  const double tmp_69 = tmp_16*tmp_18*tmp_65;
  const double tmp_70 = tmp_10*tmp_3;
  const double tmp_71 = tmp_52*tmp_70;
  const double tmp_73 = tmp_17*xx0;
  const double tmp_75 = hDD02*tmp_3*tmp_73;
  const double tmp_78 = hDD12*tmp_62*tmp_67;
  const double tmp_79 = hDD01*tmp_0*tmp_16*tmp_68 + hDD01*tmp_27*tmp_33*tmp_63 + tmp_26*tmp_63*tmp_71 + tmp_37*tmp_65*tmp_75 + tmp_37*tmp_78 + tmp_49*tmp_69;
  const double tmp_82 = 2*tmp_54;
  const double tmp_84 = tmp_28*tmp_5;
  const double tmp_85 = tmp_34*tmp_70;
  const double tmp_86 = hDD01*tmp_16*tmp_5*tmp_68 + hDD01*tmp_27*tmp_56*tmp_63 + tmp_26*tmp_63*tmp_85 + tmp_54*tmp_65*tmp_75 + tmp_54*tmp_78 + tmp_69*tmp_84;
  const double tmp_87 = tmp_16*((tmp_65)*(tmp_65));
  const double tmp_88 = tmp_16*((tmp_62)*(tmp_62));
  const double tmp_91 = (1.0/3.0)*trK;
  const double tmp_92 = aDD00*tmp_17 + tmp_18*tmp_91;
  const double tmp_94 = tmp_16*tmp_4*tmp_92;
  const double tmp_95 = aDD01*tmp_73 + hDD01*tmp_73*tmp_91;
  const double tmp_97 = tmp_16*tmp_31*tmp_95;
  const double tmp_98 = aDD11*tmp_17*tmp_7 + tmp_26*tmp_91;
  const double tmp_99 = tmp_16*tmp_28*tmp_98;
  const double tmp_100 = aDD02*tmp_3*tmp_73 + tmp_75*tmp_91;
  const double tmp_102 = aDD12*tmp_67 + hDD12*tmp_67*tmp_91;
  const double tmp_104 = tmp_16*(aDD22*tmp_17*tmp_28 + tmp_44*tmp_91);
  const double tmp_107 = tmp_100*tmp_37*tmp_84 + tmp_100*tmp_49*tmp_54 + tmp_102*tmp_37*tmp_85 + tmp_102*tmp_54*tmp_71 + tmp_104*tmp_36*tmp_53 + tmp_46*tmp_47*tmp_97 + tmp_46*tmp_94 + tmp_51*tmp_98;
  const double tmp_108 = tmp_63*tmp_95;
  const double tmp_109 = tmp_16*tmp_65*tmp_95;
  const double tmp_110 = tmp_16*tmp_65*tmp_92;
  const double tmp_114 = tmp_100*tmp_37*tmp_65 + tmp_102*tmp_37*tmp_62 + tmp_108*tmp_49 + tmp_109*tmp_71 + tmp_110*tmp_49 + tmp_63*tmp_71*tmp_98;
  const double tmp_115 = tmp_100*tmp_54*tmp_65 + tmp_102*tmp_54*tmp_62 + tmp_108*tmp_84 + tmp_109*tmp_85 + tmp_110*tmp_84 + tmp_63*tmp_85*tmp_98;
  const double tmp_116 = (1.0/(tmp_15));
  const double tmp_117 = tmp_116*vetU1;
  const double tmp_118 = tmp_116*vetU0;
  const double tmp_119 = (1.0/(xx0));
  const double tmp_120 = tmp_116*tmp_119*vetU2/tmp_3;
  gammaDD[0][0] = tmp_1*tmp_20 + tmp_14*tmp_30 + tmp_23*tmp_25 + tmp_27*tmp_33*tmp_41*tmp_42 + 2*tmp_33*tmp_37*tmp_38 + ((tmp_36)*(tmp_36))*tmp_45;
  gammaDD[0][1] = tmp_60;
  gammaDD[0][2] = tmp_79;
  gammaDD[1][0] = tmp_60;
  gammaDD[1][1] = tmp_11*tmp_30*tmp_6 + tmp_20*tmp_6 + tmp_25*tmp_47*tmp_6 + tmp_38*tmp_56*tmp_82 + tmp_45*((tmp_53)*(tmp_53)) + tmp_58*tmp_82;
  gammaDD[1][2] = tmp_86;
  gammaDD[2][0] = tmp_79;
  gammaDD[2][1] = tmp_86;
  gammaDD[2][2] = 2*hDD01*tmp_63*tmp_65*tmp_73 + tmp_18*tmp_87 + tmp_26*tmp_88;
  KDD[0][0] = tmp_1*tmp_94 + tmp_100*tmp_41*tmp_49 + tmp_102*tmp_41*tmp_71 + tmp_104*((tmp_36)*(tmp_36)) + tmp_14*tmp_99 + tmp_23*tmp_97;
  KDD[0][1] = tmp_107;
  KDD[0][2] = tmp_114;
  KDD[1][0] = tmp_107;
  KDD[1][1] = tmp_100*tmp_82*tmp_84 + tmp_102*tmp_82*tmp_85 + tmp_104*((tmp_53)*(tmp_53)) + tmp_11*tmp_6*tmp_99 + tmp_47*tmp_6*tmp_97 + tmp_6*tmp_94;
  KDD[1][2] = tmp_115;
  KDD[2][0] = tmp_114;
  KDD[2][1] = tmp_115;
  KDD[2][2] = 2*tmp_108*tmp_65 + tmp_87*tmp_92 + tmp_88*tmp_98;
  *beta0 = tmp_117*tmp_27*tmp_34 + tmp_118*tmp_49 + tmp_120*tmp_65;
  *beta1 = tmp_117*tmp_5*tmp_70 + tmp_118*tmp_71 + tmp_120*tmp_62;
  *beta2 = tmp_117*tmp_119*tmp_53 + tmp_118*tmp_36;
}
