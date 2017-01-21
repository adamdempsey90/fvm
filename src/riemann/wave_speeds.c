void wave_speeds(const double *WL, const double *WR, double *SL, double *SR, double *Sstar) {
  /* Given the left and right state estimate the left, right, and star wave speeds */
  double rho_L = WL[0];
  double u_L = WL[1];
  double p_L = WL[2];
  double rho_R = WR[0];
  double u_R = WR[1];
  double p_R = WR[2];
  
  double a_L = sound_speed(p_L,rho_L);
  double a_R = sound_speed(p_R,rho_R);
  
  double gamma_1 = (GAMMA + 1)/(2*GAMMA);
  
  double ps, us;
  
  
  anrs(WL, WR, &ps, &us);
  
  
  double qL, qR;
  
  if (ps <= p_L) {
    qL = 1.0;
  }
  else {
    qL = sqrt( 1 + gamma_1*(ps/p_L - 1) );
  }
  if (ps <= p_R) {
    qR = 1.0;
  }
  else {
    qR = sqrt( 1 + gamma_1*(ps/p_R - 1) );
  }
  
  *SL = u_L  - a_L*qL;
  *SR = u_R + a_R*qR;
  
  /* For HLL solver we would set Sstar = us here */
  *Sstar = p_R - p_L + rho_L*u_L*( *SL - u_L) - rho_R*u_R*( *SR - u_R);
  *Sstar /= rho_L *( *SL - u_L) - rho_R*( *SR - u_R);
       
    
  
 return; 
}