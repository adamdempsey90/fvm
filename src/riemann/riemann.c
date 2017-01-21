#define QUSER 2.0



void pvrs(const double *WL, const double *WR, double *ps, double *us) {
  /* Primative variable approximate riemman solver.
   * Used in smooth regions to give estimates for pressure and velocity in star region.
   */
  
 
  double rho_L = WL[0];
  double u_L = WL[1];
  double p_L = WL[2];
  double rho_R = WR[0];
  double u_R = WR[1];
  double p_R = WR[2];
  
  double a_L = sound_speed(p_L,rho_L);
  double a_R = sound_speed(p_R,rho_R);
  
  double rho__a_bar = .25*(rho_L + rho_R)*(a_L + a_R);
  
  *ps = .5*(p_L + p_R) - .5*(u_R - u_L) * rho_a__bar;
  *us = .5*(u_L + u_R) - .5*(p_R-p_L)/rho_a_bar;
  
  return;

  
}



void trrs(const double *WL, const double *WR, double *ps, double *us) {
  /* Two rarefaction approximate riemman solver.
   * Used near rarefaction waves to give estimates for pressure and velocity in star region.
   */
  
 
  double rho_L = WL[0];
  double u_L = WL[1];
  double p_L = WL[2];
  double rho_R = WR[0];
  double u_R = WR[1];
  double p_R = WR[2];
  
  double a_L = sound_speed(p_L,rho_L);
  double a_R = sound_speed(p_R,rho_R);
   
  double gamma_1 = (GAMMA - 1)/2.;
  double inv_gamma_1 = 1./gamma_1;
  double gamma_2 = (GAMMA - 1)/(2.*GAMMA);
  double inv_gamma_2 =1/gamma_2;
  double plr = pow(p_L/p_R,gamma_2);
  
  
  *ps = pow( ( a_L + a_R - gamma_1*(u_R - u_L))/(a_l*pow(p_L,-gamma_2) + a_R*pow(p_R,-gamma_2) ), inv_gamma_2);
  *us = (plr*u_L/a_L + u_R/a_R + inv_gamma_1 * (plr - 1) )/(plr/a_L + 1/a_R);
  
  
  return;

  
}



void tsrs(const double *WL, const double *WR, double *ps, double *us) {
  /* Two shock approximate riemman solver.
   * Used near shocks waves to give estimates for pressure and velocity in star region.
   * ps contains the pressure estimate from the pvrs function and is overwritten.
   */
  
 
  double rho_L = WL[0];
  double u_L = WL[1];
  double p_L = WL[2];
  double rho_R = WR[0];
  double u_R = WR[1];
  double p_R = WR[2];
  
  
  double gamma_1 = 2./(GAMMA  + 1);
  double gamma_2 = (GAMMA -1)/(GAMMA + 1);
  
  double p0 = fmax(0.0, *ps);
  
  double Al = gamma_1/rho_L;
  double Ar = gamma_1/rho_R;
  double Bl = gamma_2 * p_L;
  double Br = gamma_2 * p_R;
  double gl = sqrt( Al/(p0 + Bl));
  double gr = sqrt(Ar/(p0 +Br));
  
  *ps = (gl * p_L + gr * pr - (u_R-u_L))/(gl + gr);
  *us = .5*(u_L + u_R) + .5*( (*ps - p_R) *gr - (*ps - p_L)*gl);
  
  
  return;

  
}


void anrs(const double *WL, const double *WR, double *ps, double *us) {
  
  
  double pmax = fmax(WL[2],WR[2]);
  double pmin = fmin(WL[2],WR[2]);
  double Q = pmax/pmin;
  double pstar, ustar;
 
  
  pvrs(WL,WR,&pstar,&ustar);
  
  
  if ( (Q < QUSER) && (pstar < pmax) && (pstar > pmin)) {
    
    *ps = pstar;
    *us = ustar;
    
  }
  else {
    
    if ( pstar < pmin ) {
     
      trrs(WL, WR, &ps, &us);
      
    }
    else {
      tsrs(WL, WR, &pstar, &us);
      *ps = pstar;
    }
  }
  
  return;
  

  
}