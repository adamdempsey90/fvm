#include "defs.h"
#include "cuda_defs.h"


#ifdef PPM
#define SLOPETHETA 2.
#else
#define SLOPETHETA 1.
#endif


__device__ __inline__ real slope_limiter(real slopeL, real slopeR) {
    if (fabs(slopeR) < PRESSUREFLOOR) return  SLOPETHETA;
    real r = slopeL/slopeR;
    return MAX2(0., MIN3(SLOPETHETA*r, .5*(1+r),SLOPETHETA));
}




#ifdef PCM
__global__ void reconstruct(real *cons, real *UL, real *UR, real *dx,
        int dir1,int nx1, int nx2, int nx3, int size_x1, int size_x12,
        int nf,int ntot, int offset, real g1, real dt) {
	/*
	 * Piecewise constant reconstruction
	 */
    int i,j,k,n,indx,indxm,indxp;
    int il, iu, jl, ju, kl, ku;
    real dL, uL, uL2,uL3,pL,eL,sL;
    real dC,uC,uC2,uC3,pC,eC,sC;
    real dR,uR,uR2,uR3,pR,eR,sR;
    real dLi,dCi,dRi, dLf, uLf, dRf, uRf;
    real dxm, dxc, dxp;
    real slopeL,slopeR,ke,r, fl,fr;
    real dtdx;
    int dir2, dir3;
    /* 1->2->3
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;

    if (dir1 == 1) {
        il = -2; iu = nx1+2;
        jl = -NGHX2; ju = nx2+NGHX2;
        kl = -NGHX3; ku = nx3 + NGHX3;
    }
    else if (dir1 == 2) {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -2; ju = nx2+2;
        kl = -NGHX3; ku = nx3 + NGHX3;
    }
    else {
    	il = -NGHX1; iu = nx1+NGHX1;
		jl = -NGHX2; ju = nx2+NGHX2;
		kl = -2; ku = nx3 + 2;
    }

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=il)&&(i<iu)&&(j>=jl)&&(j<ju)&&(k>=kl)&&(k<ku)) {
            if (dir1 == 1) {
                indxm = indx - 1 ; // (i-1,j,k)
            }
            else if (dir1 == 2) {
                indxm = indx - size_x1;
            }
            else {
            	indxm = indx - size_x12;
            }
            for(n=0;n<nf;n++) {
                UL[indx + n*ntot] = cons[indx + n*ntot];
                UR[indxm + n*ntot] = cons[indx + n*ntot];
            }

        }
    }
    return;
}
#endif


#ifdef PLM
__global__ void reconstruct(real *cons, real *UL, real *UR, real *dx,
        int dir1,int nx1, int nx2, int nx3, int size_x1, int size_x12,
        int nf,int ntot, int offset, real g1, real dt) {
	/*
	 * Piecewise linear reconstruction
	 */
    int i,j,k,n,indx,indxm,indxp;
    int il, iu, jl, ju, kl, ku;
    real dL, uL, uL2,uL3,pL,eL,sL;
    real dC,uC,uC2,uC3,pC,eC,sC;
    real dR,uR,uR2,uR3,pR,eR,sR;
    real dLi,dCi,dRi, dLf, uLf, dRf, uRf;
    real dxm, dxc, dxp;
    real slopeL,slopeR,ke,r, fl,fr;
    real dtdx;
    int dir2, dir3;
    /* 1->2->3
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;

    if (dir1 == 1) {
        il = -2; iu = nx1+2;
        jl = -NGHX2; ju = nx2+NGHX2;
        kl = -NGHX3; ku = nx3 + NGHX3;
    }
    else if (dir1 == 2) {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -2; ju = nx2+2;
        kl = -NGHX3; ku = nx3 + NGHX3;
    }
    else {
    	il = -NGHX1; iu = nx1+NGHX1;
		jl = -NGHX2; ju = nx2+NGHX2;
		kl = -2; ku = nx3 + 2;
    }

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=il)&&(i<iu)&&(j>=jl)&&(j<ju)&&(k>=kl)&&(k<ku)) {
            if (dir1 == 1) {
                indxm = indx - 1 ; // (i-1,j,k)
                indxp = indx + 1 ; // (i+1,j,k)
                dxm = dx[i-1];
                dxc = dx[i];
                dxp = dx[i+1];
            }
            else if (dir1 == 2) {
                indxm = indx - size_x1;
                indxp = indx + size_x1;
                dxm = dx[j-1];
                dxc = dx[j];
                dxp = dx[j+1];
            }
            else {
            	indxm = indx - size_x12;
            	indxp = indx + size_x12;
            	dxm = dx[k-1];
            	dxc = dx[k];
            	dxp = dx[k+1];
            }
            dtdx = .5*dt/dxc;
            dL  = cons[indxm + 0*ntot];
            dLi = dL;
            uL  = cons[indxm + dir1*ntot]/dL;
            uL2 = cons[indxm + dir2*ntot]/dL;
            uL3 = cons[indxm + dir3*ntot]/dL;
            eL  = cons[indxm + 4*ntot];
            pL = (eL-.5*dL*(uL*uL+uL2*uL2+uL3*uL3))*g1;

            if (dL < PRESSUREFLOOR) dL = PRESSUREFLOOR;
            if (pL < PRESSUREFLOOR) pL = PRESSUREFLOOR;

            dC  = cons[indx + 0*ntot];
            dCi = dC;
            uC  = cons[indx + dir1*ntot]/dC;
            uC2 = cons[indx + dir2*ntot]/dC;
            uC3 = cons[indx + dir3*ntot]/dC;
            eC  = cons[indx + 4*ntot];
            pC = (eC-.5*dC*(uC*uC+uC2*uC2+uC3*uC3))*g1;

            if (dC < PRESSUREFLOOR) dC = PRESSUREFLOOR;
            if (pC < PRESSUREFLOOR) pC = PRESSUREFLOOR;

            dR  = cons[indxp + 0*ntot];
            dRi = dR;
            uR  = cons[indxp + dir1*ntot]/dR;
            uR2 = cons[indxp + dir2*ntot]/dR;
            uR3 = cons[indxp + dir3*ntot]/dR;
            eR  = cons[indxp + 4*ntot];
            pR = (eR-.5*dR*(uR*uR+uR2*uR2+uR3*uR3))*g1;
            if (dR < PRESSUREFLOOR) dR = PRESSUREFLOOR;
            if (pR < PRESSUREFLOOR) pR = PRESSUREFLOOR;


            /* Density */

            slopeR = (dR-dC) * 2./(dxc + dxp);
            slopeL = (dC-dL) * 2./(dxc + dxm);

            r = slope_limiter(slopeL,slopeR);
            dL = dC - .5*dxc * r * slopeR;
            dR = dC + .5*dxc * r * slopeR;
            if (dL < PRESSUREFLOOR) dL = PRESSUREFLOOR;
            if (dR < PRESSUREFLOOR) dR = PRESSUREFLOOR;
            dLf = dL;
            dRf = dR;

            /* ux1 */
            slopeR = (uR-uC) * 2./(dxc + dxp);
            slopeL = (uC-uL) * 2./(dxc + dxm);

            r = slope_limiter(slopeL,slopeR);
            uL = uC - .5*dxc * r * slopeR;
            uR = uC + .5*dxc * r * slopeR;

            uLf = uL;
            uRf = uR;

            /* ux2 */
            slopeR = (uR2-uC2) * 2./(dxc + dxp);
            slopeL = (uC2-uL2) * 2./(dxc + dxm);

            r = slope_limiter(slopeL,slopeR);
            uL2 = uC2 - .5*dxc * r * slopeR;
            uR2 = uC2 + .5*dxc * r * slopeR;

            /* ux3 */
            slopeR = (uR3-uC3) * 2./(dxc + dxp);
            slopeL = (uC3-uL3) * 2./(dxc + dxm);

            r = slope_limiter(slopeL,slopeR);
            uL3 = uC3 - .5*dxc * r * slopeR;
            uR3 = uC3 + .5*dxc * r * slopeR;

            /* pres */
            slopeR = (pR-pC) * 2./(dxc + dxp);
            slopeL = (pC-pL) * 2./(dxc + dxm);

            r = slope_limiter(slopeL,slopeR);
            pL = pC - .5*dxc * r * slopeR;
            pR = pC + .5*dxc * r * slopeR;

            if (pL < PRESSUREFLOOR) pL = PRESSUREFLOOR;
            if (pR < PRESSUREFLOOR) pR = PRESSUREFLOOR;

            /* Evolve cons for dt/2 */

            ke = .5*dL*(uL*uL + uL2*uL2 + uL3*uL3);
            uL *= dL;
            uL2 *= dL;
            uL3 *= dL;
            eL = (pL/g1 + ke);

            ke = .5*dR*(uR*uR + uR2*uR2 + uR3*uR3);
            uR *= dR;
            uR2 *= dR;
            uR3 *= dR;
            eR = (pR/g1 + ke);

            /* Density */
            fl = dLf*uLf; fr = dRf*uRf;
            dL += dtdx*(fl-fr);
            dR += dtdx*(fl-fr);
            UL[indx + 0*ntot] = dR;
            UR[indxm + 0*ntot] = dL;

            /* MX1 */
            fl = uLf*uL + pL; fr = uRf*uR + pR;
            uL += dtdx*(fl-fr);
            uR += dtdx*(fl-fr);
            UL[indx + dir1*ntot] = uR;
            UR[indxm + dir1*ntot] = uL;

            /* MX2 */
            fl = uLf*uL2; fr = uRf*uR2;
            uL2 += dtdx*(fl-fr);
            uR2 += dtdx*(fl-fr);
            UL[indx + dir2*ntot] = uR2;
            UR[indxm + dir2*ntot] = uL2;

            /* MX3 */
            fl = uLf*uL3; fr = uRf*uR3;
            uL3 += dtdx*(fl-fr);
            uR3 += dtdx*(fl-fr);
            UL[indx + dir3*ntot] = uR3;
            UR[indxm + dir3*ntot] = uL3;

            /* Energy */
            fl = uLf*(eL + pL); fr = uRf*(eR+pR);
            eL += dtdx*(fl-fr);
            eR += dtdx*(fl-fr);
            UL[indx + 4*ntot] = eR;
            UR[indxm + 4*ntot] = eL;




            for(n=5;n<nf;n++) {
                sR = cons[indxp + n*ntot]/dRi;
                sL = cons[indxm + n*ntot]/dLi;
                sC = cons[indx + n*ntot]/dCi;
                slopeR = (sR-sC) * 2./(dxc + dxp);
                slopeL = (sC-sL) * 2./(dxc + dxm);
                r = slope_limiter(slopeL,slopeR);
                sL = sC - .5*dxc * r * slopeR;
                sR = sC + .5*dxc * r * slopeR;

                sL *= dLf;
                sR *= dRf;

                fl = sL*uLf;
                fr = sR*uRf;
                sL += dtdx*(fl-fr);
                sR += dtdx*(fl-fr);

                UL[indx + n*ntot] = sR;
                UR[indxm + n*ntot] = sL;

            }
        }
    }

    return;
}
#endif




