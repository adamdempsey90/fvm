/* Sod shock tube */
#include "defs.h"

real sod_dens[2] = {1.,.125};
real sod_pres[2] = {1.,.1};
real sod_vel[2] = {0., 0.};

void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3) {
    *h1 = 1.;
    *h2 = 1.;
    *h3 = 1.;
}

void set_boundary_x3(GridCons *grid, Parameters *params) {
    return;
}
void set_boundary_x2(GridCons *grid, Parameters *params) {
    return;
}
void set_boundary_x1(GridCons *grid, Parameters *params) {
    int i,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;

    real *xm1 = grid->xm1;
    real *x1 = grid->xc1;

    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;

    real pres,ke,gamma_1;
    gamma_1 = params->gamma_1;
    // Inner X1 boundary
    for(i=-NGHX1;i<0;i++) {
        indx = INDEX(i);
        pres = sod_pres[0];
        rho[indx] = sod_dens[0];
        mx1[indx] = sod_dens[0]*sod_vel[0];
        mx2[indx] = 0.;
        mx3[indx] = 0.;
        ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
        ke /= 2*rho[indx];
        intenergy[indx] = pres/gamma_1;
        energy[indx] = intenergy[indx] + ke;
    

        for(n=5;n<nf;n++) {
            grid->cons[n*ntot+indx] = 0;
        }
    }

    // Outer X1 boundary
    for(i=nx1;i<nx1+NGHX1;i++) {
        indx = INDEX(i);
        pres = sod_pres[1];
        rho[indx] = sod_dens[1];
        mx1[indx] = sod_dens[1]*sod_vel[1];
        mx2[indx] = 0.;
        mx3[indx] = 0.;
        ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
        ke /= 2*rho[indx];
        intenergy[indx] = pres/gamma_1;
        energy[indx] = intenergy[indx] + ke;

        for(n=5;n<nf;n++) {
            grid->cons[n*ntot+indx] = 0;
        }
    }

        return;
}

void init_mesh(GridCons *grid, Parameters *params) {
    
    int i,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n;

    nx1 = grid->nx[0];
    size_x1 = grid->size_x[0];
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



    real *h1 = &grid->hfac[0*grid->ntot];
    real *h2 = &grid->hfac[1*grid->ntot];
    real *h3 = &grid->hfac[2*grid->ntot];


    for(i=-NGHX1;i<nx1+NGHX1;i++) {
        indx = INDEX(i);
        scale_factors(grid->xc1[i], 0,0,
                &h1[indx],&h2[indx],&h3[indx]);

    }


    return;

}

void init_gas(GridCons *grid, Parameters *params) {
    int i,indx;
    int nx1,nx2,nx3,n,ntot,nf;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;

    real *xm1 = grid->xm1;
    real *x1 = grid->xc1;


    real *rho       = &grid->cons[0*ntot];
    real *mx1       = &grid->cons[1*ntot];
    real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 

    real gamma = params->gamma;
    real gamma_1 = params->gamma_1;
    real pres, ke;

    for(i=-NGHX1;i<nx1+NGHX1;i++) {
        indx = INDEX(i); 
        if (xm1[i] < .5) {
            pres  = sod_pres[0];
            rho[indx] =sod_dens[0];
            mx1[indx] = sod_dens[0]*sod_vel[0];
        }
        else {
            pres  = sod_pres[1];
            rho[indx] =sod_dens[1];
            mx1[indx] = sod_dens[1]*sod_vel[1];
        }

        mx2[indx] = 0.;
        mx3[indx] = 0.;
        ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
        ke /= 2*rho[indx];
        intenergy[indx] = pres/gamma_1;
        energy[indx] = intenergy[indx] + ke;
        for(n=5;n<nf;n++) {
            grid->cons[n*ntot+indx] = 0;
        }
    }
    

    return;


}
