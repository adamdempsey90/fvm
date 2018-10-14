/* Rayleigh-Taylor */
#include "defs.h"
#include <time.h>
#include "cuda_defs.h"


const real g_param = .1;


__device__ real gravpot(real x, real y) {
    return g_param*y;
}

__device__ void hydrostatic_lower(int indxg, int i, int j, real *cons, real *intenergy, real *x1, real *x2, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    int n;
    int indx = GINDEX(i,0);

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + 1*ntot] = cons[indx + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx + 3*ntot];

    intenergy[indxg] = intenergy[indx] - (x2[j]-x2[0])*cons[indx+0*ntot]*g_param/(g-1);
    cons[indxg + 4*ntot] = intenergy[indxg] + .5*(
            cons[indxg + 1*ntot]*cons[indxg + 1*ntot] +
            cons[indxg + 2*ntot]*cons[indxg + 2*ntot] +
            cons[indxg + 3*ntot]*cons[indxg + 3*ntot] 
            )/cons[indxg + 0*ntot];


    for(n=5;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}
__device__ void hydrostatic_upper(int indxg, int i, int j, real *cons, real *intenergy, real *x1, real *x2, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    int n;
    int indx = GINDEX(i,nx2-1);

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + 1*ntot] = cons[indx + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx + 3*ntot];

    intenergy[indxg] = intenergy[indx] - (x2[j]-x2[nx2-1])*cons[indx+0*ntot]*g_param/(g-1);
    cons[indxg + 4*ntot] = intenergy[indxg] + .5*(
            cons[indxg + 1*ntot]*cons[indxg + 1*ntot] +
            cons[indxg + 2*ntot]*cons[indxg + 2*ntot] +
            cons[indxg + 3*ntot]*cons[indxg + 3*ntot] 
            )/cons[indxg + 0*ntot];


    for(n=5;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}
__global__ void boundary_kernel(real *cons, real *intenergy, real *x1, real *x2, int nx1, int nx2, int size_x1, int nf, int ntot, int offset, real g, real time) {

    int i,j,indxg;
    for(indxg = blockIdx.x*blockDim.x + threadIdx.x; indxg<ntot; indxg+=blockDim.x*gridDim.x) {
        j = indxg/size_x1;
        i = indxg -size_x1*j - NGHX1;
        j -= NGHX2;
        if ((i>=-NGHX1)&&(i<0)&&(j>=-NGHX2)&&(j<nx2+NGHX2)) {
        /* Lower x1 */
            periodic_boundary_x1_inner(indxg,i,j,cons,intenergy,nx1,nx2,ntot,nf,size_x1,offset,g,time);

        }
        else if ((j>=-NGHX2)&&(j<0)&&(i>=-NGHX1)&&(i<nx1+NGHX1)) {
        /* Lower x2 */

            hydrostatic_lower(indxg,i,j,cons,intenergy,x1,x2,nx1,nx2,ntot,nf,size_x1,offset,g,time);
        }
        else if ((i>=nx1)&&(i<nx1+NGHX1)&&(j>=-NGHX2)&&(j<nx2+NGHX2))  {
        /* Upper x1 */
            periodic_boundary_x1_outer(indxg,i,j,cons,intenergy,nx1,nx2,ntot,nf,size_x1,offset,g,time);

        }
        else if ((j>=nx2)&&(j<nx2+NGHX2)&&(i>=-NGHX1)&&(i<nx1+NGHX1)) {
        /* Upper x2 */
            hydrostatic_upper(indxg,i,j,cons,intenergy,x1,x2,nx1,nx2,ntot,nf,size_x1,offset,g,time);


        }

    }
    return;
}


//extern "C" {
void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3) {
    *h1 = 1.;
    *h2 = 1.;
    *h3 = 1.;
}
void init_mesh(GridCons *grid, Parameters *params) {
    
    int i,j,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    real xi, xo, dx;
    real *xm ;
    real *xc ;
    real *size;
    

    /* x1 direction */
    n = nx1;
    xm = grid->xm1;
    xc = grid->xc1;
    size= grid->dx1;
    xi = params->x1_min;
    xo = params->x1_max;
    dx = (xo-xi)/(float)n;
    for(i=-NGHX1;i<n+NGHX1+1;i++) {
        xm[i] =  xi + i*dx;
    }
    for(i=-NGHX1;i<n+NGHX1;i++) {
        xc[i] =  (xm[i+1]+xm[i])*.5;
        size[i] = xm[i+1]-xm[i];
    }

    /* x2 direction */
    n = nx2;
    xm = grid->xm2;
    xc = grid->xc2;
    size= grid->dx2;
    xi = params->x2_min;
    xo = params->x2_max;
    dx = (xo-xi)/(float)n;
    for(i=-NGHX2;i<n+NGHX2+1;i++) {
        xm[i] =  xi + i*dx;
    }
    for(i=-NGHX2;i<n+NGHX2;i++) {
        xc[i] =  (xm[i+1]+xm[i])*.5;
        size[i] = xm[i+1]-xm[i];
    }


    real *h1 = &grid->hfac[0*grid->ntot];
    real *h2 = &grid->hfac[1*grid->ntot];
    real *h3 = &grid->hfac[2*grid->ntot];


    for(j=-NGHX2;j<nx2+NGHX2;j++) {
        for(i=-NGHX1;i<nx1+NGHX1;i++) {
         indx = INDEX(i,j);
         scale_factors(grid->xc1[i], grid->xc2[j],0,
                &h1[indx],&h2[indx],&h3[indx]);

        }
    }


    return;

}
void init_gas(GridCons *grid, Parameters *params) {
    int i,j,indx;
    int nx1,nx2,nx3,n,ntot,nf;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;

    real *xm1 = grid->xm1;
    real *x1 = grid->xc1;

    real *xm2 = grid->xm2;
    real *x2 = grid->xc2;

    real *dx1 = grid->dx1;
    real *dx2 = grid->dx2;

    real *rho       = &grid->cons[0*ntot];
    real *mx1       = &grid->cons[1*ntot];
    real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 

    real gamma = params->gamma;
    real gamma_1 = params->gamma_1;
    real pres, ke;

    real u1 = 0;
    real u2 = 0;

    real norm;
    srand(time(NULL));


    for(j=-NGHX2;j<nx2+NGHX2;j++) {
        for(i=-NGHX1;i<nx1+NGHX1;i++) {
            indx = INDEX(i,j); 
            //indx = i + size_x1*j;

            u2 = 0;
            if (fabs(x2[j]) <= .25) {
                u1 = .5;
                rho[indx] = 2.;
            }
            else {
                u1 = -.5;
                rho[indx] = 1.;
            }

            

            pres = 2.5;
            cs = sqrt(gamma * pres/rho[indx]);
            norm =(real)((double)rand() / (double)RAND_MAX );
            u1 += (2*norm-1)*.01;
            u2 += (2*norm-1)*.01;
            mx1[indx] = u1*rho[indx];
            mx2[indx] = u2*rho[indx];
            mx3[indx] = 0.;

            ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
            ke /= 2*rho[indx];
            intenergy[indx] = pres/gamma_1;
            energy[indx] = intenergy[indx] + ke;
            for(n=5;n<nf;n++) {
                grid->cons[n*ntot+indx] = 0;
            }





        }
    }

    return;


}
//}
