void star_region(const double *UL, const double *UR, double *Fhllc) {
  /* Construct the left and right star conservative states with the HLLC solver */
  double WL[5], WR[5];
  double WL1D[3], WR1d[3];
  
  /* Unpack all the conservative variables */
  double rho_L = UL[0];
  double momx_L = UL[1];
  double momy_L = UL[2];
  double momz_L = UL[3];
  double e_L = UL[4];
  
  double rho_R = UR[0];
  double momx_R = UR[1];
  double momy_R = UR[2];
  double momz_R = UR[3];
  double e_R = UR[4];
  
  /* Convert to primative variables */
  cons_to_prim(UL,WL);
  cons_to_prim(UR,WR);
  
  /* Unpack primative variables */
  double u_L = WL[1];
  double v_L = WL[2];
  double w_L = WL[3];
  double p_L = WL[4];
  
  double u_R = WR[1];
  double v_R = WR[2];
  double w_R = WR[3];
  double p_R = WR[4];
  
  double fac_L, fac_R;
  double SL, SR, Sstar;
  double UsL[5], UsR[5];
  
  
  /* First get wave speeds */
  wave_speeds(WL1D, WR1D, &SL, &SR, &Sstar);
  
  if (SL >= 0.0) {
    cons_to_flux(UL,WL, Fhllc);
  }
  else if (Sstar > = 0) {
    fac_L = rho_L * (SL - u_L)/(S_L - Sstar);
    cons_to_flux(UL,WL, Fhllc);
    Fhllc[0] += SL*(fac_L - UL[0]) ;
    Fhllc[1] += SL*(fac_L * Sstar -UL[1]);
    Fhllc[2] += SL*(fac_L * v_L - UL[2]);
    Fhllc[3] += SL*(fac_L * v_L - UL[3]);
    Fhllc[4] += SL*(fac_L * ( e_L/rho_L + (Sstar - u_L) *( Sstar + p_L/(rho_L*(S_L - u_L))) ) - UL[4]);
  }
  else if (SR >= 0) {
    fac_R = rho_R * (SR - u_R)/(S_R - Sstar);
    cons_to_flux(UR,WR, Fhllc);
    Fhllc[0] += SR*(fac_R - UR[0]) ;
    Fhllc[1] += SR*(fac_R * Sstar -UR[1]);
    Fhllc[2] += SR*(fac_R * v_R - UR[2]);
    Fhllc[3] += SR*(fac_R * v_R - UR[3]);
    Fhllc[4] += SR*(fac_R * ( e_R/rho_R + (Sstar - u_R) *( Sstar + p_R/(rho_R*(S_R - u_R))) ) - UR[4]);
  }
  else {
    cons_to_flux(UL,Wl,Fhllc);
  }
  
    
  
  
  

  return;

  
  
  
}