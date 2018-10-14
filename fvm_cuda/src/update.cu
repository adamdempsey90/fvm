#include "defs.h"
#include "cuda_defs.h"


__global__ void update_cons(real *cons, real *intenergy, real *F_1, real *F_2,
        real *dx1, real *dx2, real dt, int nx1, int nx2, int size_x1, int ntot, int offset, int nf) {
    int i,j,n;
    int indx;
    int indxm1, indxm2;
    real dtdx1, dtdx2;


    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)) {
            indxm1 = GINDEX(i-1,j);
            indxm2 = GINDEX(i,j-1);
            dtdx1 = dt/dx1[i];
            dtdx2 = dt/dx2[j];
            for(n=0;n<nf;n++) {
                cons[indx + n*ntot] += dtdx1*(F_1[indxm1 + n*ntot]-F_1[indx + n*ntot]) + 
                              dtdx2*(F_2[indxm2 + n*ntot]-F_2[indx + n*ntot]);

            }
            intenergy[indx] = cons[indx+4*ntot] - .5*(
                    cons[indx + 1*ntot]*cons[indx + 1*ntot] +
                    cons[indx + 2*ntot]*cons[indx + 2*ntot] +
                    cons[indx + 3*ntot]*cons[indx + 3*ntot])/cons[indx];
        }
    }
    return;

}
__global__ void transverse_update(real *UL_1, real *UL_2,
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
            indxm2 = GINDEX(i,j-1);
            indxp2 = GINDEX(i,j+1);
            indxp1m2 = GINDEX(i+1,j-1);
            dtdx2 = .5*dt/dx2[j];
            for(n=0;n<nf;n++) {
                UL_1[indx + n*ntot] += dtdx2*(F_2[indxm2 + n*ntot]-F_2[indx + n*ntot]);
                UR_1[indx + n*ntot] += dtdx2*(F_2[indxp1m2 + n*ntot]-F_2[indxp2 + n*ntot]);

            }
        }
        /* X2 - direction */
        if ((i>=-1)&&(i<nx1+1)&&(j>=-1)&&(j<nx2)) {
            indxm1 = GINDEX(i-1,j);
            indxp1 = GINDEX(i+1,j);
            indxp2m1 = GINDEX(i-1,j+1);
            dtdx1 = .5*dt/dx1[i];
            for(n=0;n<nf;n++) {
                UL_2[indx + n*ntot] += dtdx1*(F_1[indxm1 + n*ntot]-F_1[indx + n*ntot]);
                UR_2[indx + n*ntot] += dtdx1*(F_1[indxp2m1 + n*ntot]-F_1[indxp1 + n*ntot]);

            }
        }
    }
    return;

}
