/* Sod shock tube */
#include "defs.h"

void init_mesh(GridCons *grid, Parameters *params) {
    
    int i,j,k,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1]; 
    nx3 = grid->nx[2];
    size_x1 = grid->size_x[0];
    size_x12 = grid->size_x12;
    real xi, xo, dx;
    real *xm ;
    real *xc ;
    real *size;
    

    /* x1 direction */
    n = grid->nx[0];
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
    n = grid->nx[1];
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

    /* x3 direction */
    n = grid->nx[2];
    xm = grid->xm3;
    xc = grid->xc3;
    size= grid->dx3;
    xi = params->x3_min;
    xo = params->x3_max;
    dx = (xo-xi)/(float)n;
    for(i=-NGHX3;i<n+NGHX3+1;i++) {
        xm[i] =  xi + i*dx;
    }
    for(i=-NGHX3;i<n+NGHX3;i++) {
        xc[i] =  (xm[i+1]+xm[i])*.5;
        size[i] = xm[i+1]-xm[i];
    }
    


    real *h1 = &grid->hfac[0*grid->ntot + grid->offset];
    real *h2 = &grid->hfac[1*grid->ntot + grid->offset];
    real *h3 = &grid->hfac[2*grid->ntot + grid->offset];


    for(k=-NGHX3;k<nx3+NGHX3;k++) {
        for(j=NGHX2;j<nx2+NGHX2;j++) {
            for(i=-NGHX1;i<nx1+NGHX1;i++) {
                indx = INDEX(i,j,k);
                scale_factors(grid->xc1[i], grid->xc2[j],grid->xc3[k],
                        &h1[indx],&h2[indx],&h3[indx]);

            }
        }
    }


    return;

}

void init_gas(GridCons *grid, Parameters *params) {
    int i,j,k,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1]; 
    nx3 = grid->nx[2];
    size_x1 = grid->size_x[0];
    size_x12 = grid->size_x12;

    real *xm1 = grid->xm1;
    real *xm2 = grid->xm2;
    real *xm3 = grid->xm3;
    real *x1 = grid->xc1;
    real *x2 = grid->xc2;
    real *x3 = grid->xc3;


    real *rho       = &grid->cons[0*grid->ntot + grid->offset];
    real *mx1       = &grid->cons[1*grid->ntot + grid->offset];
    real *mx2       = &grid->cons[2*grid->ntot + grid->offset];
    real *mx3       = &grid->cons[3*grid->ntot + grid->offset];
    real *energy    = &grid->cons[4*grid->ntot + grid->offset];
    real *intenergy = &grid->cons[5*grid->ntot + grid->offset];

    real gamma = params->gamma;
    real gamma_1 = params->gamma_1;
    real pres, ke;

    for(k=-NGHX3;k<nx3+NGHX3;k++) {
        for(j=NGHX2;j<nx2+NGHX2;j++) {
            for(i=-NGHX1;i<nx1+NGHX1;i++) {
                indx = INDEX(i,j,k); 
                pres = (xm1[i] < .5) ? 1. : .1;
                rho[indx] = pres;
                mx1[indx] = 0.;
                mx2[indx] = 0.;
                mx3[indx] = 0.;
                ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
                ke /= 2*rho[indx];
                intenergy[indx] = pres/gamma_1;
                energy[indx] = intenergy[indx] + ke;
            }
        }
    }
    

    return;


}
