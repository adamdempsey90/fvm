#include "defs.h"

/* Piecewise Linear reconstruction */


#define SLOPETHETA 1.

void plm(real *cons, real dt, real *dx, real gamma_1,int dir1, int nf) {
    int n;

    real slopeL, slopeR, slopeC, r;
    real fl, fr, dL,uL,pL,dR,uR,pR,ke;

    real dtdx = .5*dt/dx[1];

    real *WL = &cons[0];
    real *WC = &cons[nf];
    real *WR = &cons[2*nf];


    /* Convert to primitive */

    WL[1] /= WL[0];
    WL[2] /= WL[0];
    WL[3] /= WL[0];
    ke = .5*WL[0]*(WL[1]*WL[1] + WL[2]*WL[2] + WL[3]*WL[3]);
    WL[4] = (WL[4] - ke)*gamma_1;
    for(n=5;n<nf;n++) WL[n] /= WL[0];

    WC[1] /= WC[0];
    WC[2] /= WC[0];
    WC[3] /= WC[0];
    ke = .5*WC[0]*(WC[1]*WC[1] + WC[2]*WC[2] + WC[3]*WC[3]);
    WC[4] = (WC[4] - ke)*gamma_1;
    for(n=5;n<nf;n++) WC[n] /= WC[0];

    WR[1] /= WR[0];
    WR[2] /= WR[0];
    WR[3] /= WR[0];
    ke = .5*WR[0]*(WR[1]*WR[1] + WR[2]*WR[2] + WR[3]*WR[3]);
    WR[4] = (WR[4] - ke)*gamma_1;
    for(n=5;n<nf;n++) WR[n] /= WR[0];


    WL[0] = fmax(PRESSUREFLOOR, WL[0]);
    WC[0] = fmax(PRESSUREFLOOR, WC[0]);
    WR[0] = fmax(PRESSUREFLOOR, WR[0]);
    WL[4] = fmax(PRESSUREFLOOR, WL[4]);
    WC[4] = fmax(PRESSUREFLOOR, WC[4]);
    WR[4] = fmax(PRESSUREFLOOR, WR[4]);

    for(n=0;n<nf;n++) {

        slopeR = (WR[n]-WC[n]) *2/(dx[1] +dx[2]);
        slopeL = (WC[n]-WL[n]) * 2./(dx[0]+dx[1]);
        if (fabs(slopeR) < PRESSUREFLOOR) {
            r = SLOPETHETA;
        }
        else {
            r = slopeL/slopeR;
            r = fmax(0., fmin(fmin( SLOPETHETA*r, .5*(1+r)),SLOPETHETA));
        }

        //r = (r + fabs(r))/(1. + fabs(r));

    
        WL[n] = WC[n] - .5*dx[1]* r *slopeR;
        WR[n] = WC[n] + .5*dx[1]* r *slopeR;

    }
    WL[0] = fmax(WL[0],PRESSUREFLOOR);
    WR[0] = fmax(WR[0],PRESSUREFLOOR);
    WL[4] = fmax(WL[4],PRESSUREFLOOR);
    WR[4] = fmax(WR[4],PRESSUREFLOOR);
    
    dL = WL[0]; uL = WL[1]; pL = WL[4];
    dR = WR[0]; uR = WR[1]; pR = WR[4];


    /* Convert to conservative */

    ke = .5*dL*(WL[1]*WL[1] + WL[2]*WL[2] +WL[3]*WL[3]);
    WL[1] *= dL;
    WL[2] *= dL;
    WL[3] *= dL;
    WL[4] = (pL/gamma_1 + ke);
    for(n=5;n<nf;n++) WL[n] *= dL;

    ke = .5*dR*(WR[1]*WR[1] + WR[2]*WR[2] +WR[3]*WR[3]);
    WR[1] *= dR;
    WR[2] *= dR;
    WR[3] *= dR;
    WR[4] = (pR/gamma_1 + ke);
    for(n=5;n<nf;n++) WR[n] *= dR;

    /* Evolve conservative for dt/2 */

    /* Density */
    fl = dL*uL; fr = dR*uR;
    WL[0] += dtdx*(fl-fr);
    WR[0] += dtdx*(fl-fr);

    /* MX1 */
    fl = uL*WL[1] + pL; fr = uR*WR[1] + pR;
    WL[1] += dtdx*(fl-fr);
    WR[1] += dtdx*(fl-fr);

    /* MX2 */
    fl = uL*WL[2]; fr = uR*WR[2];
    WL[2] += dtdx*(fl-fr);
    WR[2] += dtdx*(fl-fr);

    /* MX3 */
    fl = uL*WL[3]; fr = uR*WR[3];
    WL[3] += dtdx*(fl-fr);
    WR[3] += dtdx*(fl-fr);

    /* Energy */
    fl = uL*(WL[4] + pL); fr = uR*(WR[4]+pR);
    WL[4] += dtdx*(fl-fr);
    WR[4] += dtdx*(fl-fr);

    /* Scalars */
    for(n=5;n<nf;n++) {
        fl = uL*WL[n];
        fr = uR*WR[n];
        WL[n] += dtdx*(fl-fr);
        WR[n] += dtdx*(fl-fr);
    }


    return;
}

