void cons_to_prim(const double *U, double *W) {
 
  W[0] = U[0];
  W[1] = U[1]/U[0];
  W[2] = U[2]/U[0];
  W[3] = U[3]/U[0];
  W[4] = 
  
  
}


double get_slope_x(const double uL, const double u, const double uR) {
  double dqm, dqp;
  dqm =  u-uL;
  dqp = uR-u;
  
  return slope_limiter(dqm,dqp);
  

}

void piecewise_linear_reconstruction(const double *UL, const double *U, const double *UR,
				     double *URL, double *URR) {
  
  double dum, dup, dx
  int i;
  for(i=0;i<NFIELDS;i++) {
    
    dq = get_slope_x(UL[i], U[i], UR[i]);
    
    URL[i] = U[i] - .5*dq;
    URR[i] = U[i] + .5*dq;
  }
  
  return;
}

void reconstruction_evolve(const double *UL, const double *UR, double *newUL, double *newUR, const double dtdx) {
 
  double FL[NFIELDS], FR[NFIELDS];
  
  flux_evaluation(UL, FL);
  flux_evaluation(UR,FR);
  
  for (i=0;i<NFIELDS;i++) {
    newUL[i] = UL[i] + .5*dtdx*(FL[i] - FR[i]);
    newUR[i] = UR[i] + .5*dtdx*(FL[i] - FR[i]);
  }
  
  return;
  
  
}