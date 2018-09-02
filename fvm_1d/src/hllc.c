#include "defs.h"

void hll_flux(const real *UL, const real *UR, real *F,real g, real g1, real gp, int nf) {
    /* Flux using HLL scheme */
    int n;
    real WL[3], WR[3];
    real ps, us;
#ifdef HLLC
    real Sstar,Ufac; 
#else
    real sfacL,sfacR,sfacLR;
#endif
    real SL, SR, dL,uL,uL2,uL3,pL, eL,aL;
    real ke;
    real dR,uR,uR2,uR3,pR,eR,aR;

    dL = UL[0];
    uL = UL[1]/dL;
    uL2 = UL[2] / dL;
    uL3 = UL[3] / dL;
    eL = UL[4];
    pL = (eL- .5*dL*(uL*uL + uL2*uL2 + uL3*uL3))*g1;

    dR = UR[0];
    uR = UR[1]/dR;
    uR2 = UR[2] / dR;
    uR3 = UR[3] / dR;
    eR = UR[4];
    pR = (eR- .5*dR*(uR*uR + uR2*uR2 + uR3*uR3))*g1;


    pL = fmax(pL,PRESSUREFLOOR);
    pR = fmax(pR,PRESSUREFLOOR);
    
    aL = sqrt(g*pL/dL);
    aR = sqrt(g*pR/dR);

    

    WL[0] = dL; WL[1] = uL; WL[2] = pL;
    WR[0] = dR; WR[1] = uR; WR[2] = pR;


    real gamma_1 = g-1;
    real gamma_2 = gamma_1/(2*g);
    real gamma_3 = gamma_1/(g+1);
    real gamma_4 = 2./(g+1);
    real gamma_5 = gamma_2/gamma_3;

    anrs(WL,WR,aL,aR,g,&ps,&us,gamma_1,gamma_2,gamma_3,gamma_4);
    int use_left = (us>=0);
    
    SL = (ps<=pL) ? uL-aL : uL-aL*sqrt(1 + gamma_5*(ps/pL-1));
    SR = (ps<=pR) ? uR+aR : uR + aR*sqrt(1 + gamma_5*(ps/pR-1));
#ifdef HLLC
    Sstar = (pR-pL) + dL*uL*(SL-uL) - dR*uR*(SR-uR);
    Sstar /=  dL*(SL-uL) - dR*(SR-uR);
#endif
    if (SL >= 0) {
        /* Flux is due to left state */
        F[0] = dL*uL;
        F[1] = dL*uL*uL + pL;
        F[2] = dL*uL*uL2;
        F[3] = dL*uL*uL3;
        F[4] = (eL+pL)*uL;
        for(n=5;n<nf;n++) {
            F[n] = dL*uL * UL[n]/dL;
        }
    }
    else {
        if (SR <= 0) {
            F[0] = dR*uR;
            F[1] = dR*uR*uR + pR;
            F[2] = dR*uR*uR2;
            F[3] = dR*uR*uR3;
            F[4] = (eR+pR)*uR;
            for(n=5;n<nf;n++) {
                F[n] = dR*uR * UR[n]/dR;
            }

        }
        else {
#ifndef HLLC
            /* HLL Flux 
             * (SR*FL - SL*FR + SL*SR*(UR-UL) ) /(SR-SL)
             */
            sfacL = SR/(SR-SL);
            sfacR = -SL/(SR-SL);
            sfacLR = SL*SR/(SR-SL);

            F[0] = sfacL * dL*uL + sfacR * dR*uR + sfacLR*(dR-dL);
            F[1] = sfacL * (dL*uL*uL + pL) + sfacR*(dR*uR*uR +pR) + sfacLR*(dR*uR-dL*uL);
            F[2] = sfacL * dL*uL*uL2 + sfacR*dR*uR*uR2 + sfacLR*(dR*uR2-dL*uL2);
            F[3] = sfacL * dL*uL*uL3 + sfacR*dR*uR*uR3+ sfacLR*(dR*uR3-dL*uL3);
            F[4] = sfacL * (eL+pL)*uL + sfacR*(eR + pR)*uR + sfacLR*(eR -eL);
            for(n=5;n<nf;n++) {
                F[n] = sfacL *uL * UL[n] + sfacR*uR*UL[n] + sfacLR*(UR[n]-UL[n]);
            }
#else
            if (Sstar >= 0) {
            /* HLLC Flux 
             * FL +SL*(UstarL-UL)
             * UstarL = dL*(SL-uL)/(SL-Sstar) * {
             * 1
             * Sstar
             * vL
             * wL
             * EL/dL + (Sstar-uL)*(Sstar + pL/(dL*(SL-uL)))
             * }
             */
                Ufac = dL*(SL - uL)/(SL-Sstar);
                F[0] = dL*uL;
                F[0] += SL*( Ufac  - UL[0]); 
                
                F[1] = dL*uL*uL + pL;
                F[1] += SL*( Ufac * Sstar  - UL[1]); 

                F[2] = dL*uL*uL2;
                F[2] += SL*( Ufac * uL2  - UL[2]); 

                F[3] = dL*uL*uL3;
                F[3] += SL*( Ufac * uL3 - UL[3]); 

                F[4] = (eL+pL)*uL;
                F[4] += SL*( Ufac * (eL/dL + (Sstar-uL)*(Sstar + pL/(dL*(SL-uL)))) - UL[4]); 

                for(n=5;n<nf;n++) {
                    F[n] = dL*uL * UL[n]/dL;
                    F[n] += SL*( Ufac *UL[n]/dL - UL[n]); 
                }


            }
            else {
            /* HLLC Flux 
             * FR +SR*(UstarR-UR)
             * UstarR = dR*(SR-uR)/(SR-Sstar) * {
             * 1
             * Sstar
             * vR
             * wR
             * ER/dR + (Sstar-uR)*(Sstar + pR/(dR*(SR-uR)))
             * }
             */
                Ufac = dR*(SR - uR)/(SR-Sstar);
                F[0] = dR*uR;
                F[0] += SR*( Ufac  - UR[0]); 
                
                F[1] = dR*uR*uR + pR;
                F[1] += SR*( Ufac * Sstar  - UR[1]); 

                F[2] = dR*uR*uR2;
                F[2] += SR*( Ufac * uR2  - UR[2]); 

                F[3] = dR*uR*uR3;
                F[3] += SR*( Ufac * uR3 - UR[3]); 

                F[4] = (eR+pR)*uR;
                F[4] += SR*( Ufac * (eR/dR + (Sstar-uR)*(Sstar + pR/(dR*(SR-uR)))) - UR[4]); 

                for(n=5;n<nf;n++) {
                    F[n] = dR*uR * UR[n]/dR;
                    F[n] += SR*( Ufac *UR[n]/dR - UR[n]); 
                }
            }
#endif
        }
    }
        



    return;
}

