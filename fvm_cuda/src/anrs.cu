#include "defs.h"
#include "cuda_defs.h"
#define QUSER 2

/* ANRS Riemann solver */
__device__ void anrs(const real dL,const real uL, const real pL, const real aL,
        const real dR, const real uR, const real pR, const real aR, 
        real *ps, real *us, real g, real g1,real g2,real g3,real g4) {
  real pmax,pmin,Q;
  real pstar, ustar;
  real gl,gr,plr;


  if (pL < pR) {
      pmax = pR;
      pmin = pL;
  }
  else {
      pmax = pL;
      pmin = pR;
  }

  Q = pmax/pmin;
  
 
  
  real rho_a_bar = .25*(dL + dR)*(aL + aR);
  
  pstar = .5*(pL + pR) - .5*(uR - uL) * rho_a_bar;
  ustar = .5*(uL + uR) - .5*(pR-pL)/rho_a_bar;
  
  pstar = fmax(0.,pstar);
  
  if ( (Q < QUSER) && (pstar < pmax) && (pstar > pmin)) {
    
    *ps = pstar;
    *us = ustar;
    
  }
  else {
    
    if ( pstar < pmin ) {
     
        plr = pow(pL/pR,g2);
  
        *ps = pow( ( aL + aR - .5*g1*(uR - uL))/(aL*pow(pL,-g2) + aR*pow(pR,-g2) ), 1./g2);
        *us = (plr*uL/aL + uR/aR +  (plr - 1)*2/g1 )/(plr/aL + 1./aR);
      
    }
    else {
        /* Two-Shock */
        gl = sqrt(g4/dL / (pstar + g3*pL));
        gr = sqrt(g4/dR / (pstar + g3*pR));
        *ps = (gl * pL + gr * pR - (uR-uL))/(gl + gr);
        *us = .5*(uL + uR) + .5*( (pstar - pR) *gr - (pstar - pL)*gl);
    }
  }
  
  if (*ps < PRESSUREFLOOR) *ps = PRESSUREFLOOR;
  return;
}
