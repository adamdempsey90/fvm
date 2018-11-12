/* Implosion */
#include "defs.h"
#include <time.h>
#include "cuda_defs.h"




__device__ void x2_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	if (x1[i] > 1./6) {
		reflecting_boundary_inner(2,indxg,i,j,k,cons,intenergy,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	}
	return;
}
__device__ void x2_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	real u1,u2;
	real *rho       = &cons[0*ntot];
	real *mx1       = &cons[1*ntot];
	real *mx2       = &cons[2*ntot];
	real *mx3       = &cons[3*ntot];
	real *energy    = &cons[4*ntot];

	if (x1[i] < 1./6 + (1. + 20*time)/sqrt(3.)) {
		rho[indxg] = 8.;
		intenergy[indxg] = 291.25;
		u1 = 8.25*sqrt(3.)/2.;
		u2 = -8.25*.5;
	}
	else {
		rho[indxg] = 1.4;
		u1 = 0.;
		u2 = 0.;
		intenergy[indxg] = 2.5;
	}
	mx1[indxg] = u1*rho[indxg];
	mx2[indxg] = u2*rho[indxg];
	mx3[indxg] = 0;

	energy[indxg] = intenergy[indxg] + (mx1[indxg]*mx1[indxg] + mx2[indxg]*mx2[indxg] + mx3[indxg]*mx3[indxg])/(2*rho[indxg]);
	return;
}
__device__ void x1_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	/* ICs */
	return;
}

__device__ void x1_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	outflow_boundary_outer(1,indxg,i,j,k,cons,intenergy,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}


void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3) {
	*h1 = 1.;
	*h2 = 1.;
	*h3 = 1.;
}
void init_mesh(GridCons *grid, Parameters *params) {

	init_uniform_mesh(grid,params);

	return;

}
void init_gas(GridCons *grid, Parameters *params) {
    int i,j,k,indx;
    int nx1,nx2,nx3,n,ntot,nf;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    nx3 = grid->nx[2];
    size_x1 = grid->size_x1;
    size_x12 = grid->size_x12;
    ntot = grid->ntot;
    nf = grid->nf;

    real *x1 = grid->xc1;

	real *x2 = grid->xc2;


	real *rho       = &grid->cons[0*ntot];
	real *mx1       = &grid->cons[1*ntot];
	real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 

    real gamma_1 = params->gamma - 1;
    real ke;

    real u1 = 0;
    real u2 = 0;
    real u3 = 0;
    real energ;



    for(k=-NGHX3;k<nx3+NGHX3;k++) {
		for(j=-NGHX2;j<nx2+NGHX2;j++) {
			for(i=-NGHX1;i<nx1+NGHX1;i++) {
				indx = INDEX(i,j,k);

			    if (x1[i] < 1./6 + x2[j]/sqrt(3.0) + 1e-8) {
			        rho[indx] = 8.0;
			        energ = 291.25;
			        u1 = 8.25*sqrt(3.)/2.;
			        u2 = -8.25*.5;
			    }
			    else {
			        rho[indx] = 1.4;
			        energ = 2.5;
			        u1 = 0.;
			        u2 = 0.;
			    }


				mx1[indx] = u1*rho[indx];
				mx2[indx] = u2*rho[indx];
				mx3[indx] = u3*rho[indx];

				ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
				ke /= 2*rho[indx];
				intenergy[indx] = energ;
				energy[indx] = intenergy[indx] + ke;
				for(n=5;n<nf;n++) {
					grid->cons[n*ntot+indx] = 0;
				}

			}
		}
    }
    return;


}
