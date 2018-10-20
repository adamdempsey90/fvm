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
            /* x1-direction */
            phic = gravpot(x1[i],x2[j]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j]);
            cons[indx + 1*ntot] -= dt/dx1[i] * (phir - phil) * dhalf[indx];

            /* energy */
            cons[indx + 4*ntot] -= dt/dx1[i]*(F1[indxm1] *(phic-phil) 
                    + F1[indx]*(phir-phic));

            /* x2-direction */
            phil = gravpot(x1[i], x2[j] - .5*dx2[j]);
            phir = gravpot(x1[i], x2[j] + .5*dx2[j]);
            cons[indx + 2*ntot] -= dt/dx2[j] * (phir - phil) * dhalf[indx];

            /* energy */
            cons[indx + 4*ntot] -= dt/dx2[j]*(F2[indxm2] *(phic-phil) 
                    + F2[indx]*(phir-phic));



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
                gl = (phif - phil)*2/dxc;
                gr = (phic - phif)*2/dxc;
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
                gl = (phif - phil)*2/dxc;
                gr = (phic - phif)*2/dxc;
#endif
            }
#ifdef POTENTIAL
            UL[indx + dir1*ntot] -= .5*dt*gr*UL[indx + 0*ntot];
            UR[indxm + dir1*ntot] -= .5*dt*gl*UR[indxm + 0*ntot];
#endif



        }

    }
    return;

__global__ void source_transverse_update(real *UL_1, real *UL_2,
        real *UR_1, real *UR_2, 
        real *F_1, real *F_2, real *dx1, real *dx2, real dt, int nx1, int nx2, int size_x1, int ntot, int offset, int nf) {
    int i,j,n;
    int indx;
    int indxm1,indxp1,indxp2m1, indxp1m2, indxm2, indxp2;
    real dtdx2,dtdx1;

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot;indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        /* X1 - direction */
        if ((i>=-1)&&(i<nx1)&&(j>=-1)&&(j<nx2+1)) {
            indxm1 = GINDEX(i-1,j);
            indxm2 = GINDEX(i,j-1);
            indxp2 = GINDEX(i+1,j); // i,j+1
            indxp1m2 = GINDEX(i+1,j-1);
            dtdx2 = .5*dt/dx2[j];
            phic = gravpot(x1[i],x2[j]);
            phir = gravpot(x1[i],x2[j] + .5*dx2);
            phil = gravpot(x1[i],x2[j] - .5*dx2);
            UL_1[indx + 1*ntot] -= dtdx2*(phir-phic);
            UR_1[indxm1 + 1*ntot] -= dtdx2*(phic-phil);

        }
        /* X2 - direction */
        if ((i>=-1)&&(i<nx1+1)&&(j>=-1)&&(j<nx2)) {
            indxm1 = GINDEX(i-1,j);
            indxm2 = GINDEX(i,j-1);
            indxp1 = GINDEX(i,j+1); // i+1,j
            indxp2m1 = GINDEX(i-1,j+1);
            dtdx1 = .5*dt/dx1[i];
            phic = gravpot(x1[i],x2[j]);
            phir = gravpot(x1[i] + .5*dx1,x2[j]);
            phil = gravpot(x1[i] - .5*dx1,x2[j]);
            UL_2[indx + 2*ntot] -= dtdx1*(phir-phic);
            UR_2[indxm2 + 2*ntot] -= dtdx1*(phic-phil);
        }
    }
    return;

}
