#include "defs.h"
#include "cuda_defs.h"


__global__ void riemann_fluxes(real *UL, real *UR, real *F, 
        int dir1,int nx1, int nx2, int size_x1, 
        int nf,int ntot, int offset, real g) {
    int i,j,n,indx;
    real dL,uL,pL,eL,aL,uL2,uL3;
    real dR,uR,pR,eR,aR,uR2,uR3;
#ifdef EXACT
    int use_left;
    real ds,us,ps;
#endif

#ifdef HLLC
    real Ufac, SL, SR, Sstar;
#endif
    real g1 = g-1;
    real g2 = g1/(2*g);
    real g3 = g1/(g+1);
    real g4 = 2./(g+1);
    real g5 = g2/g3;

    int il,iu,jl,ju;

    int dir2, dir3; 

    /* 1->2->3 
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;

    if (dir1 == 1) {
        il = -1; iu = nx1+1;
        jl = -NGHX2;ju = nx2+NGHX2;
    }
    else if (dir1 == 2) {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -1; ju = nx2+1;
    }
    else {
        printf("Direction can only be 1 or 2 in 2D!\n");
    }
    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        if ((i>=il)&&(i<iu)&&(j>=jl)&&(j<ju)) {

            dL = UL[indx + 0*ntot];
            uL   = UL[indx + dir1*ntot]/dL;
            uL2  = UL[indx + dir2*ntot]/dL;
            uL3  = UL[indx + dir3*ntot]/dL;
            eL   = UL[indx + 4*ntot];
            pL   = (eL - .5*(uL*uL + uL2*uL2 + uL3*uL3)*dL)*g1;
            aL   = sqrt(g*pL/dL);

            dR = UR[indx + 0*ntot];
            uR   = UR[indx + dir1*ntot]/dR;
            uR2  = UR[indx + dir2*ntot]/dR;
            uR3  = UR[indx + dir3*ntot]/dR;
            eR   = UR[indx + 4*ntot];
            pR   = (eR - .5*(uR*uR + uR2*uR2 + uR3*uR3)*dR)*g1;
            aR   = sqrt(g*pR/dR);

            if (pL < PRESSUREFLOOR) pL = PRESSUREFLOOR;
            if (pR < PRESSUREFLOOR) pR = PRESSUREFLOOR;

#ifdef HLLC

            hllc(dL,uL,pL,aL,
                    dR,uR,pR,aR,
                    &SL,&SR,&Sstar,
                    g,g1,g2,g3,g4,g5);

            if (SL>=0) {
                F[indx + 0*ntot] = dL*uL;
                F[indx + dir1*ntot] = dL*uL*uL+pL;
                F[indx + dir2*ntot] = dL*uL*uL2;
                F[indx + dir3*ntot] = dL*uL*uL3;
                F[indx + 4*ntot] = (eL+pL)*uL;
                for(n=5;n<nf;n++) {
                    F[indx + n*ntot] =dL*uL*UL[indx+n*ntot]/dL;
                } 
            }
            else {
                if (SR <= 0) {
                    F[indx + 0*ntot] = dR*uR;
                    F[indx + dir1*ntot] = dR*uR*uR+pR;
                    F[indx + dir2*ntot] = dR*uR*uR2;
                    F[indx + dir3*ntot] = dR*uR*uR3;
                    F[indx + 4*ntot] = (eR+pR)*uR;
                    for(n=5;n<nf;n++) {
                        F[indx + n*ntot] =dR*uR*UR[indx+n*ntot]/dR;
                    } 
                }
                else {
                    if (Sstar >= 0) {
                        Ufac = dL*(SL-uL)/(SL-Sstar);
                        F[indx + 0*ntot] = dL*uL + 
                            SL*(Ufac - UL[indx+0*ntot]);
                        F[indx + dir1*ntot] = dL*uL*uL + pL +
                            SL*(Ufac * Sstar - UL[indx + dir1*ntot]);
                        F[indx + dir2*ntot] = dL*uL*uL2 + 
                            SL*(Ufac * uL2 - UL[indx + dir2*ntot]);
                        F[indx + dir3*ntot] = dL*uL*uL3 + 
                            SL*(Ufac * uL3 - UL[indx + dir3*ntot]);
                        F[indx + 4*ntot] = (eL+pL)*uL + 
                            SL*(Ufac *(eL/dL + (Sstar-uL)*(Sstar+pL/(dL*(SL-uL))))-UL[indx + 4*ntot]);
                        for(n=5;n<nf;n++) {
                            F[indx + n*ntot] = dL*uL*UL[indx + n*ntot]/dL + 
                                SL*(Ufac * UL[indx + n*ntot]/dL - UL[indx + n*ntot]);
                        }
                    }
                    else {
                        Ufac = dR*(SR-uR)/(SR-Sstar);
                        F[indx + 0*ntot] = dR*uR + 
                            SR*(Ufac - UR[indx+0*ntot]);
                        F[indx + dir1*ntot] = dR*uR*uR + pR +
                            SR*(Ufac * Sstar - UR[indx + dir1*ntot]);
                        F[indx + dir2*ntot] = dR*uR*uR2 + 
                            SR*(Ufac * uR2 - UR[indx + dir2*ntot]);
                        F[indx + dir3*ntot] = dR*uR*uR3 + 
                            SR*(Ufac * uR3 - UR[indx + dir3*ntot]);
                        F[indx + 4*ntot] = (eR+pR)*uR + 
                            SR*(Ufac *(eR/dR + (Sstar-uR)*(Sstar+pR/(dR*(SR-uR))))-UR[indx + 4*ntot]);
                        for(n=5;n<nf;n++) {
                            F[indx + n*ntot] = dR*uR*UR[indx + n*ntot]/dR + 
                                SR*(Ufac * UR[indx + n*ntot]/dR - UR[indx + n*ntot]);
                        }

                    }
                    
                }
            }
#endif
#ifdef EXACT
            use_left = exact_sample(dL,uL,pL,aL,
                    dR,uR,pR,aR,
                    &ds,&us,&ps,
                    g,g1,g2,g3,g4,g5,0.,EXACT_TOL);
            
            F[indx + 0*ntot] = ds*us;
            F[indx + dir1*ntot] = ds*us*us + ps;
            if (use_left) {
                F[indx + dir2*ntot] = ds*us*uL2 ;
                F[indx + dir3*ntot] = ds*us*uL3 ;
                F[indx + 4*ntot] = us*( ps/g1 + .5*ds*(us*us + uL2*uL2 + uL3*uL3) + ps);
                for(n=5;n<nf;n++) {
                    F[indx + n*ntot] = ds*us*UL[indx + n*ntot]/dL;
                }
            }
            else {
                F[indx + dir2*ntot] = ds*us*uR2 ;
                F[indx + dir3*ntot] = ds*us*uR3 ;
                F[indx + 4*ntot] = us*( ps/g1 + .5*ds*(us*us + uR2*uR2 + uR3*uR3) + ps);
                for(n=5;n<nf;n++) {
                    F[indx + n*ntot] = ds*us*UR[indx + n*ntot]/dR;
                }
            }

#endif

            
        }
    }
    

    return;

}


