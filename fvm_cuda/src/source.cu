#include "defs.h"
#include "cuda_defs.h"


#ifdef POTENTIAL
    extern __device__ real gravpot(real x, real y);
#endif

__global__ void update_source(real *cons, real *dhalf, real *F1, real *F2, real *dx1, real *dx2, real *x1, real *x2,
        int nx1, int nx2, int size_x1, 
        int nf,int ntot, int offset, real dt) {


    int i,j,n,indx;
    int indxp1, indxp2;
    int il, iu, jl, ju;
    int indxm1,indxm2;
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
    int indxm;
    real dtdx,xc,xp,xm,xc2,xc1;
#ifdef POTENTIAL
    real phir,phil,phic;
#endif
    int dir2, dir3;
    /* 1->2->3 
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;

    if (dir1 == 1) {
        il = -2; iu = nx1+NGHX1;
        jl = -NGHX2;ju = nx2+NGHX2;
    }
    else if (dir1 == 2) {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -2; ju = nx2+NGHX2;
    }
    else {
        printf("Direction can only be 1 or 2 in 2D!\n");
    }
    
    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
#ifdef POTENTIAL
        if ((i>=il)&&(i<iu)&&(j>=jl)&&(j<ju)) {
            if (dir1 == 1) {
                indxm = GINDEX(i-1,j);
                dtdx = dt/dx[i];
                phir = gravpot(x1[i] + .5*dx[i],x2[j]);
                phil = gravpot(x1[i] - .5*dx[i],x2[j]);
                phic = gravpot(x1[i]           ,x2[j]);
            }
            else {
                indxm = GINDEX(i,j-1);
                dtdx = dt/dx[j];
                phir = gravpot(x1[i],x2[j] + .5*dx[j]);
                phil = gravpot(x1[i],x2[j] - .5*dx[j]);
                phic = gravpot(x1[i],x2[j]);
            }
            UL[indx + dir1*ntot] -= dtdx*(phir-phic)*UL[indx + 0*ntot];
            UR[indxm + dir1*ntot] -= dtdx*(phic-phil)*UR[indxm + 0*ntot];

        }
#endif

    }
    return;

}
__global__ void source_transverse_update(real *cons, real *UL_1, real *UL_2, real *UR_1, real *UR_2, real *F_1, real *F_2, real *dx1, real *dx2,real *x1, real *x2, real dt, int nx1, int nx2, int size_x1, int ntot, int offset, int nf) {
    int i,j,n;
    int indx;
    int indxm1,indxp1,indxp2m1, indxp1m2, indxm2, indxp2;
    int indxm1m2;
    real dtdx2,dtdx1;
#ifdef POTENTIAL
    real phil,phir,phic;
#endif

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot;indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        /* X1 - direction */
        if ((i>=-2)&&(i<nx1+1)&&(j>=-2)&&(j<nx2+1)) {
            indxm1 = GINDEX(i-1,j);
            indxm2 = GINDEX(i,j-1);
            indxm1m2 = GINDEX(i-1,j-1);
            indxp2 = GINDEX(i+1,j); // i,j+1
            indxp1m2 = GINDEX(i+1,j-1);
            dtdx2 = .5*dt/dx2[j];
#ifdef POTENTIAL
            phic = gravpot(x1[i],x2[j]);
            phir = gravpot(x1[i],x2[j] + .5*dx2[j]);
            phil = gravpot(x1[i],x2[j] - .5*dx2[j]);
            
            UL_1[indx + 2*ntot] -= dtdx2*(phir-phil)*cons[indx];
            UL_1[indx + 4*ntot] -= dtdx2*((phir-phic)*F_2[indx] + (phic-phil)*F_2[indxm2]);

            UR_1[indxm1 + 2*ntot] -= dtdx2*(phir-phil)*cons[indx];
            UR_1[indxm1 + 4*ntot] -= dtdx2*((phir-phic)*F_2[indx] + (phic-phil)*F_2[indxm2]);

            /*
            phic = gravpot(x1[i-1],x2[j]);
            phir = gravpot(x1[i-1],x2[j] + .5*dx2[j]);
            phil = gravpot(x1[i-1],x2[j] - .5*dx2[j]);
            */

#endif

        }
        /* X2 - direction */
        if ((i>=-1)&&(i<nx1+1)&&(j>=-1)&&(j<nx2+1)) {
            indxm1 = GINDEX(i-1,j);
            indxm2 = GINDEX(i,j-1);
            indxm1m2 = GINDEX(i-1,j-1);
            indxp1 = GINDEX(i,j+1); // i+1,j
            indxp2m1 = GINDEX(i-1,j+1);
            dtdx1 = .5*dt/dx1[i];
#ifdef POTENTIAL
            phic = gravpot(x1[i],x2[j]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j]);
            UL_2[indx + 1*ntot] -= dtdx1*(phir-phil)*cons[indx];
            UL_2[indx + 4*ntot] -= dtdx1*((phir-phic)*F_1[indx] + (phic-phil)*F_1[indxm1]);

            UR_2[indxm2 + 1*ntot] -= dtdx1*(phir-phil)*cons[indx];

            UR_2[indxm2 + 4*ntot] -= dtdx1*((phir-phic)*F_1[indx] + (phic-phil)*F_1[indxm1]);

            /*
            phic = gravpot(x1[i],x2[j-1]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j-1]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j-1]);
            */
            
#endif
        }
    }
    return;

}
