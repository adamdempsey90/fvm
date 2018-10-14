#include "defs.h"
#include "cuda_defs.h"
#define NITERMAX 30


/* Exact solver functions */




__device__ void fk(const real p, const real dk, const real pk, const real ak, const real Ak, const real Bk, const real g1,const real g2, const real g5, real *f, real *df) {

    if (p > pk) {
        *f =  (p-pk)*sqrt(Ak/(p+Bk));
        *df = ( 1 - .5*(p-pk)/(Bk+p)  ) * sqrt(Ak/(p+Bk));
    }
    else {
        *f =  2*ak/g1 *( pow(p/pk,g2) - 1.);
        *df = 1./(ak*dk) * pow(p/pk,-g5);
    }
    return;

}

__device__ int exact_star_region(const real dL, const real uL, 
        const real pL, const real aL,
        const real dR, const real uR, const real pR, const real aR,
        const real g, const real g1, const real g2, 
        const real g3, const real g4, const real g5,
        real *ps, real *us,  real tol) {
    /* Exact Riemann solver */
    int niter;


    real pstar, ustar;


    real du = uR - uL;

    /* Use ANRS for initial guess */

    real Al = g4 / dL;
    real Ar = g4 / dR;
    real Bl = g3 * pL;
    real Br = g3 * pR;
    
    /* Newton-Raphson iteration to convergence */
    real pold;
    anrs(dL,uL,pL,aL,
            dR,uR,pR,aR,
            &pold,&ustar,
            g,g1,g2,g3,g4);

    real resid = 1.0;

    niter = 0;

    real fR,fL,dfR,dfL;
    while ((resid > tol)&&(niter < NITERMAX)) {
        fk(pold,dR, pR, aR, Ar,Br,g1,g2,g5, &fR, &dfR);
        fk(pold,dL, pL, aL, Al,Bl,g1,g2,g5, &fL, &dfL);

        pstar = pold - (fL+fR + du)/(dfL+dfR); 
        if (tol > pstar)
            pstar = tol;

        resid = fabs(pstar-pold)*2/(pold+pstar);
        niter += 1;
        pold = pstar;
    }

    *us = .5*(uL + uR+  fR-fL);
    *ps = pstar;

    return niter;

}  

__device__ int exact_sample(const real dL, const real uL, const real pL, const real aL,
        const real dR, const real uR, const real pR, const real aR,
        real *d_f, real *u_f, real *p_f, 
        real g, real g1, real g2, real g3, real g4, real g5,
        real S, real tol) {
    /* Full solution to the Riemann problem with exact solver at x/t 
     * includes support for vacuum.
     */



    real ps, us;
    exact_star_region(dL,uL,pL,aL,dR,uR,pR,aR,
            g,g1,g2,g3,g4,g5,
            &ps,&us,EXACT_TOL);


    real SL,SHL,STL,SR, SHR, STR;

    /* For passive scalars and normal velocities
     * we return whether or not to take 
     * the left or right state
     * True = left
     * False = right
     */


    if (S < us) {
    // Left of contact
        if (ps > pL) {
            // Shock
            SL = uL - aL*sqrt( g5 * ps/pL + g2);
            if (S <= SL) {
                *u_f = uL;
                *p_f = pL;
                *d_f = dL;
                return TRUE;
            }
            else {
                *u_f = us;
                *p_f = ps;
                *d_f = dL*(ps/pL +g3)/(g3*ps/pL+1); 
                return TRUE;
            }
        }
        else {
            // Rarefaction
            SHL = uL - aL;
            STL = us - aL* pow(ps/pL,g2);
            if (S <= SHL) {
                *u_f = uL;
                *p_f = pL;
                *d_f = dL;
                return TRUE;
            }
            else {
                if (S <= STL) {
                    *u_f = g4 *(aL + g1*uL/2. + S);
                    *d_f = dL*pow(g4 + (uL-S)*g3/aL,2./g1);
                    *p_f = pL*pow(g4 + (uL-S)*g3/aL,1/g2); 
                    return TRUE;
                }
                else {
                    *u_f = us;
                    *p_f = ps;
                    *d_f = dL *pow(ps/pL,1./g);
                    return TRUE;
                }

            }
        }
    }
    else {
        if (ps > pR) {
            SR = uR + aR*sqrt( g5 *ps/pR + g2);
            if (S<=SR) {
                *u_f = us;
                *p_f = ps;
                *d_f = dR*(ps/pR + g3)/(g3 *ps/pR+1);
                return FALSE;
            }
            else {
                *u_f = uR;
                *p_f = pR;
                *d_f = dR;
                return FALSE;
            }
        }
        else {
            SHR = uR + aR;
            STR = us + aR* pow(ps/pR,g2);
            if (S <= STR) {
                *u_f = us;
                *p_f = ps;
                *d_f = dR * pow(ps/pR,1./g);
                return FALSE;
            }
            else {
                if (S <= SHR) {
                    *u_f = g4 *(-aR + g1*uR/2. + S);
                    *p_f = pR*pow(g4 - g3 *(uR-S)/aR , 1./g2);
                    *d_f = dR*pow( g4 - g3*(uR-S)/aR, 2./g1);
                    return FALSE;
                }
                else {
                    *u_f = uR;
                    *p_f = pR;
                    *d_f = dR;
                    return FALSE;
                }
            }
        }
    }

}
  
