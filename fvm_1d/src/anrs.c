#include "defs.h"
#define QUSER 2

/* ANRS Riemann solver */
void anrs(const real *WL, const real *WR,real aL, real aR, real g, real *ps, real *us, real gamma_1,real gamma_2,real gamma_3,real gamma_4) {
  real pmax,pmin,Q;
  real pstar, ustar;
  real gl,gr,plr;

  real rhoL = WL[0];
  real uL = WL[1];
  real pL = WL[2];

  real rhoR = WR[0];
  real uR = WR[1];
  real pR = WR[2];

  pmax = fmax(pL,pR);
  pmin = fmin(pL,pR);
  Q = pmax/pmin;
  
 
  
  real rho_a_bar = .25*(rhoL + rhoR)*(aL + aR);
  
  pstar = .5*(pL + pR) - .5*(uR - uL) * rho_a_bar;
  ustar = .5*(uL + uR) - .5*(pR-pL)/rho_a_bar;
  
  pstar = fmax(0.,pstar);
  
  if ( (Q < QUSER) && (pstar < pmax) && (pstar > pmin)) {
    
    *ps = pstar;
    *us = ustar;
    
  }
  else {
    
    if ( pstar < pmin ) {
     
        plr = pow(pL/pR,gamma_2);
  
        *ps = pow( ( aL + aR - .5*gamma_1*(uR - uL))/(aL*pow(pL,-gamma_2) + aR*pow(pR,-gamma_2) ), 1./gamma_2);
        *us = (plr*uL/aL + uR/aR +  (plr - 1)*2/gamma_1 )/(plr/aL + 1./aR);
      
    }
    else {
        /* Two-Shock */
        gl = sqrt(gamma_4/rhoL / (pstar + gamma_3*pL));
        gr = sqrt(gamma_4/rhoR / (pstar + gamma_3*pR));
        *ps = (gl * pL + gr * pR - (uR-uL))/(gl + gr);
        *us = .5*(uL + uR) + .5*( (pstar - pR) *gr - (pstar - pL)*gl);
    }
  }
  
    *ps = fmax(PRESSUREFLOOR,*ps);
    return;
}
