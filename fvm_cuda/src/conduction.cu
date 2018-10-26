#include "defs.h"
#include "cuda_defs.h"


#ifdef CONDUCTION

extern __device__ real heatcond_func(real dens, real x1, real x2, real delad);

__global__ void conduction_flux(real *cons, real *intenergy, real *F_1, real *F_2,
        real *dx1, real *dx2, real *x1, real *x2, real g, real dt, int nx1, int nx2, int size_x1, int ntot, int offset, int nf) {
    int i,j,n;
    int indx,indxp;
    real dtdx1, dtdx2;

    real cond, tempc, tempr;
    real gamma_1 = g-1;
    real delad = 1. -1./g;
    
    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        
        /* X1 direction */ 
        if ((i>=-NGHX1)&&(i<nx1+2)&&(j>=-NGHX2)&&(j<nx2+NGHX2)) {
            indxp = GINDEX(i+1,j);
            tempc = intenergy[indx] * g / cons[indx]; // Cp = 1 -> rho*e = Cp*T/gamma
            tempr = intenergy[indxp] * g / cons[indxp];
            cond = heatcond_func(.5*(cons[indx]+cons[indxp]),x1[i]+dx1[i]*.5,x2[j],delad);
                
            F_1[indx] = cond * (tempr - tempc) / dx1[i];

        }
        /* X2 direction */
        if (nx2 > 1) {
            if ((i>=-NGHX1)&&(i<nx1+NGHX1)&&(j>=-NGHX2)&&(j<nx2+2)) {
                indxp = GINDEX(i,j+1);
                tempc = intenergy[indx] * g / cons[indx];
                tempr = intenergy[indxp] * g / cons[indxp];
                cond = heatcond_func(.5*(cons[indx]+cons[indxp]),x1[i],x2[j]+.5*dx2[j],delad);
                    
                F_2[indx] = cond * (tempr - tempc) / dx2[j];

            }
        }
        else {
            F_2[indx] = 0.;
        }

            
    }
    return;

}

__global__ void conduction_update(real *cons, real *intenergy, real *F_1, real *F_2,
        real *dx1, real *dx2, real dt, int nx1, int nx2, int size_x1, int ntot, int offset, int nf) {
    int i,j,n;
    int indx;
    int indxm1, indxm2;
    real dtdx1, dtdx2;
    real rhs;


    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)) {
            indxm1 = GINDEX(i-1,j);
            indxm2 = GINDEX(i,j-1);
            dtdx1 = dt/dx1[i];
            dtdx2 = dt/dx2[j];
            rhs = dtdx1*(F_1[indx]-F_1[indxm1]) + dtdx2*(F_2[indx ]-F_2[indxm2]);
            cons[indx + 4*ntot] += rhs;
            intenergy[indx] += rhs;

        }
    }
    return;

}
#endif
