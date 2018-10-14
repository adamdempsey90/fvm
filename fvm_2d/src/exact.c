#include "defs.h"
#define NITERMAX 30


/* Exact solver functions */



real fk_prime(const real p, const real pk, const real dk, const real ak, const real Ak, const real Bk, const real gamma_5) {

    if (p > pk) {
        return ( 1 - .5*(p-pk)/(Bk+p)  ) * sqrt(Ak/(p+Bk));
    }
    else {
        return 1./(ak*dk) * pow(p/pk,-gamma_5);
    }

}

real fk(const real p, const real pk, const real ak, const real Ak, const real Bk, const real gamma_1,const real gamma_2) {

    if (p > pk) {
        return (p-pk)*sqrt(Ak/(p+Bk));
    }
    else {
        return 2*ak/gamma_1 *( pow(p/pk,gamma_2) - 1.);
    }

}

int exact_star_region(const real *WL, const real *WR,  real g, real *ps, real *us, real tol) {
    /* Exact Riemann solver */
    int niter;


    real pstar, ustar;


    real dL = WL[0];
    real uL = WL[1];
    real pL = WL[2];
    real aL = sqrt(g*pL/dL);

    real dR = WR[0];
    real uR = WR[1];
    real pR = WR[2];
    real aR = sqrt(g*pR/dR);

    real du = uR - uL;

    
    real gamma_1 = g-1;
    real gamma_2 = gamma_1/(2*g);
    real gamma_3 = gamma_1/(g+1);
    real gamma_4 = 2./(g+1);
    real gamma_5 = gamma_2/gamma_3;

    /* Use ANRS for initial guess */

    real Al = gamma_4 / dL;
    real Ar = gamma_4 / dR;
    real Bl = gamma_3 * pL;
    real Br = gamma_3 * pR;
    
    /* Newton-Raphson iteration to convergence */
    real pold;
//  pold = .5*(pL + pR);
    anrs(WL,WR,aL,aR,g,&pold,&ustar,gamma_1,gamma_2,gamma_3,gamma_4);
    //printf("Guessing %.3e\n",pold);

    real resid = 1.0;

    niter = 0;

    real fR,fL,dfR,dfL;
    while ((resid > tol)&&(niter < NITERMAX)) {
        fR = fk(pold,pR,aR, Ar,Br,gamma_1,gamma_2);
        fL = fk(pold,pL,aL, Al,Bl,gamma_1,gamma_2);
        dfR = fk_prime(pold,pR,dR,aR,Ar,Br,gamma_5);
        dfL = fk_prime(pold,pL,dL,aL,Al,Bl,gamma_5); 

        pstar = pold - (fL+fR + du)/(dfL+dfR); 
        pstar = fmax(pstar,tol);

        resid = fabs(pstar-pold)*2/(pold+pstar);
        niter += 1;
        pold = pstar;
    }

    *us = .5*(uL + uR+  fR-fL);
    *ps = pstar;

    return niter;

}  

int exact_sample(const real *WL, const real *WR, const real ps, const real us, real *W, real g, real S, real tol) {
    /* Full solution to the Riemann problem with exact solver at x/t 
     * includes support for vacuum.
     */


    /* Sample solution for left and right density in star region */
    real dL = WL[0];
    real uL = WL[1];
    real pL = WL[2];
    real aL = sqrt(g*pL/dL);
    real dR = WR[0];
    real uR = WR[1];
    real pR = WR[2];
    real aR = sqrt(g*pR/dR);
    
    real gp = (g+1)/(2*g);
    real gm = (g-1)/(2*g);
    real g3 = (g-1)/(g+1);


    real d_f, u_f, p_f;

    real SL,SHL,STL,SR, SHR, STR;

    pL = fmax(pL,PRESSUREFLOOR);
    pR = fmax(pR,PRESSUREFLOOR);

    /* Left and right stats are not vacuum */
    aL = sqrt(g*pL/dL);
    aR = sqrt(g*pR/dR);


    if (S < us) {
    // Left of contact
        if (ps > pL) {
            // Shock
            SL = uL - aL*sqrt( gp * ps/pL + gm);
            if (S <= SL) {
                u_f = uL;
                p_f = pL;
                d_f = dL;
            }
            else {
                u_f = us;
                p_f = ps;
                d_f = dL*(ps/pL +g3)/(g3*ps/pL+1); 
            }
        }
        else {
            // Rarefaction
            SHL = uL - aL;
            STL = us - aL* pow(ps/pL,gm);
            if (S <= SHL) {
                u_f = uL;
                p_f = pL;
                d_f = dL;
            }
            else {
                if (S <= STL) {
                    u_f = 2./(g+1) *(aL + (g-1)*uL/2. + S);
                    d_f = dL*pow(2./(g+1) + (g-1)*(uL-S)/(aL*(g+1)),2./(g-1));
                    p_f = pL*pow(2./(g+1) + (uL-S)*(g-1)/(aL*(g+1)),1/gm); 
                }
                else {
                    u_f = us;
                    p_f = ps;
                    d_f = dL *pow(ps/pL,1./g);
                }

            }
        }
    }
    else {
        if (ps > pR) {
            SR = uR + aR*sqrt( gp *ps/pR + gm);
            if (S<=SR) {
                p_f = ps;
                u_f = us;
                d_f = dR*(ps/pR + g3)/(g3 *ps/pR+1);
            }
            else {
                p_f = pR;
                d_f = dR;
                u_f = uR;
            }
        }
        else {
            SHR = uR + aR;
            STR = us + aR* pow(ps/pR,gm);
            if (S <= STR) {
                u_f = us;
                p_f = ps;
                d_f = dR * pow(ps/pR,1./g);
            }
            else {
                if (S <= SHR) {
                    u_f = 2./(g+1) *(-aR + (g-1)*uR/2. + S);
                    p_f = pR*pow(2./(g+1) - g3 *(uR-S)/aR , 1./gm);
                    d_f = dR*pow( 2./(g+1) - g3*(uR-S)/aR, 2./(g-1));
                }
                else {
                    u_f = uR;
                    p_f = pR;
                    d_f = dR;
                }
            }
        }
    }



    /* Done */
    W[0] = d_f;
    W[1] = u_f;
    W[2] = p_f;

    /* For passive scalars and normal velocities
     * we return whether or not to take 
     * the left or right state
     * True = left
     * False = right
     */
    return (us >= S);

}
  
void exact_flux(const real *UL, const real *UR, real *F,real g, real g1, int nf) {
    /* Flux using exact riemann solver.
     * U = (rho, rho*u, rho*v, rho*w, E, rho*si)
     * F = (rho*u, rho*u*u+p, rho*u*v, rho*u*w,u*(E+p), rho*u*si)
     * G = (rho*v, rho*u*v, rho*v*v+p, rho*v*w,u*(E+p), rho*v*si)
     * H = (rho*w, rho*u*w, rho*v*w, rho*w*w+p, rho*w*w,w*(E+p), rho*w*si)
     */
    /* Flux using HLL scheme */
    int n;
    real WL[3], WR[3];
    real WS[3] = {0,0,0};
    real ps, us, ds;
    real SL, SR, rhoL,uL,pL, eL,aL,uL2,uL3;
    real rhoR,uR,pR,eR,aR,uR2,uR3;

    rhoL = UL[0];
    uL = UL[1]/rhoL;
    uL2 = UL[2]/rhoL;
    uL3 = UL[3]/rhoL;
    eL = UL[4];
    pL = (eL - .5*(uL*uL + uL2*uL2 + uL3*uL3)*rhoL)*g1;

    rhoR = UR[0];
    uR = UR[1]/rhoR;
    uR2 = UR[2]/rhoR;
    uR3 = UR[3]/rhoR;
    eR = UR[4];
    pR = (eR - .5*(uR*uR + uR2*uR2 + uR3*uR3)*rhoR)*g1;
    
    pL = fmax(pL,PRESSUREFLOOR);
    pR = fmax(pR,PRESSUREFLOOR);


    

    WL[0] = rhoL; WL[1] = uL; WL[2] = pL;
    WR[0] = rhoR; WR[1] = uR; WR[2] = pR;


    exact_star_region(WL,WR,g,&ps,&us,1e-6);
    int use_left = exact_sample(WL,WR,ps,us,WS,g,0.,1e-6);

    ds = WS[0];
    us = WS[1];
    ps = WS[2];

    F[0] = ds*us;
    F[1] = ds*us*us + ps;
    if (use_left) {
        F[2] = ds*us*uL2 ;
        F[3] = ds*us*uL3 ;
        F[4] = us*( ps/g1 + .5*ds*(us*us + uL2*uL2 + uL3*uL3) + ps);
        for(n=5;n<nf;n++) {
            F[n] = ds*us*UL[n]/rhoL;
        }
    }
    else {
        F[2] = ds*us*uR2 ;
        F[3] = ds*us*uR3 ;
        F[4] = us*( ps/g1 + .5*ds*(us*us + uR2*uR2 + uR3*uR3) + ps);
        for(n=5;n<nf;n++) {
            F[n] = ds*us*UR[n]/rhoR;
        }
    }


    return;
}
