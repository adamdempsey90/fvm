/* Kelvin-Helmholtz  */
#include "defs.h"

void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3) {
    *h1 = 1.;
    *h2 = 1.;
    *h3 = 1.;
}

void set_boundary_x3(GridCons *grid, Parameters *params) {
    /* Set boundary conditions in x3direction */
    return;
}
void set_boundary_x2(GridCons *grid, Parameters *params) {
    /* Set boundary conditions in x2 direction */
    reflecting_boundary_x2_inner(grid,params);
    reflecting_boundary_x2_outer(grid,params);
    return;
}
void set_boundary_x1(GridCons *grid, Parameters *params) {
    /* Set boundary conditions in x1 direction */
    reflecting_boundary_x1_inner(grid,params);
    reflecting_boundary_x1_outer(grid,params);
    return;
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

    printf("NTOT %d\n",ntot);
    real *rho       = &grid->cons[0*ntot];
    real *mx1       = &grid->cons[1*ntot];
    real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 

    real gamma = params->gamma;
    real gamma_1 = params->gamma_1;
    real pres, ke;

    for(j=0;j<nx2;j++) {
        for(i=0;i<nx1;i++) {
            //indx = INDEX(i,j); 
            indx = i + size_x1*j;

            if (xm2[j+1]  <=  .15 - xm1[i+1]) {
            //if (xm2[j+1]  <= .1) {
     //           printf("%d\t%d\t%lg\t%lg%lg\n",i,j,x1[i],x2[j],.15-x1[i]);
                rho[indx] = .125;
                pres = .14;
            }
            else {
                rho[indx] = 1.;
                pres = 1.;
            }

            mx1[indx] = 0.;
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
    }

   
    set_boundary_x1(grid,params);
    set_boundary_x2(grid,params);


    return;


}
