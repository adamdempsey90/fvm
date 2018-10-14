#include "defs.h"
#include "cuda_defs.h"

__global__ void update_source(real *cons, real *dhalf, real *F1, real *F2, real *dx1, real *dx2, real *x1, real *x2,
        int nx1, int nx2, int size_x1, 
        int nf,int ntot, int offset, real g1, real dt) {


    int i,j,n,indx;
    int indxp1, indxp2;
    int il, iu, jl, ju;
#ifdef POTENTIAL
    real phil,phir,phic;
#endif


    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)) {
            indxp1 = GINDEX(i+1,j);
            indxp2 = GINDEX(i,j+1);
            indxm1 = GINDEX(i-1,j);
            indxm2 = GINDEX(i,j-1);


#ifdef POTENTIAL
            /* indxm1 ?? */
            /* x1-direction */
            phic = gravpot(x1[i],x2[j]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j]);
            cons[indx + 1*ntot] -= dt/dx1[i] * (phir - phil) * dhalf[indx];

            /* energy */
            cons[indx + 4*ntot] -= dt/dx1[i]*(F1[indx] *(phic-phil) 
                    + F1[indxp1]*(phir-phic));

            /* x2-direction */
            phil = gravpot(x1[i], x2[j] - .5*dx2[j]);
            phir = gravpot(x1[i], x2[j] + .5*dx2[j]);
            cons[indx + 2*ntot] -= dt/dx2[j] * (phir - phil) * dhalf[indx];

            /* energy */
            cons[indx + 4*ntot] -= dt/dx2[j]*(F2[indx] *(phic-phil) 
                    + F2[indxp2]*(phir-phic));



#endif



        }
    }
    return;
}

__global__ void source_terms(real *UL, real *UR, real *dx, real *x1, real *x2,
        int dir1,int nx1, int nx2, int size_x1, 
        int nf,int ntot, int offset, real g1, real dt) {

    int i,j,n,indx;
    int il, iu, jl, ju;
    real dxc,xc,xm,xf,xc2,xc1;
#ifdef POTENTIAL
    real phil,phif,phic,gl,gr;
#endif
    int dir2, dir3;
    /* 1->2->3 
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;

    if (dir1 == 1) {
        il = -2; iu = nx1+2;
        jl = -NGHX2;ju = nx2+NGHX2;
    }
    else if (dir1 == 2) {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -2; ju = nx2+2;
    }
    else {
        printf("Direction can only be 1 or 2 in 2D!\n");
    }
    
    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        if ((i>=il)&&(i<iu)&&(j>=jl)&&(j<ju)) {
            if (dir1 == 1) {
                dxc = dx[i];
                xm = x1[i-1];
                xf = x1[i] - dxc*.5;
                xc = x1[i];
                xc2 = x2[j];
#ifdef POTENTIAL
                phil = gravpot(xm,xc2);
                phif = gravpot(xf,xc2);
                phic = gravpot(xc,xc2);
                gl = (phif - phil)*.5/dxc;
                gr = (phic - phif)*.5/dxc;
#endif
            }
            else {
                dxc = dx[j];
                xm = x2[j-1];
                xf = x2[j] - dxc*.5;
                xc = x2[j];
                xc1 = x1[i];
#ifdef POTENTIAL
                phil = gravpot(xc1,xm);
                phif = gravpot(xc1,xf);
                phic = gravpot(xc1,xc);
                gl = (phif - phil)*.5/dxc;
                gr = (phic - phif)*.5/dxc;
#endif
            }
#ifdef POTENTIAL
            UL[indx + dir1*ntot] -= .5*dt*gl*UL[indx + 0*ntot];
            UR[indx + dir1*ntot] -= .5*dt*gr*UR[indx + 0*ntot];
#endif



        }

    }
    return;
}
