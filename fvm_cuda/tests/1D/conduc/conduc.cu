/* 1D Sod Shock Tube */
#include "defs.h"
#include <time.h>
#include "cuda_defs.h"

__device__ __managed__ real cond = .01;


__device__ real heatcond_func(real dens, real x1, real x2, real x3,real delad) {
    return cond;
}
__device__ real thermal_diff(real dens, real x1, real x2, real x3,real delad) {
    return heatcond_func(dens,x1,x2,x3,delad)/dens/delad;
}

/* Select boundary conditions */
__device__ void x1_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	//outflow_boundary_inner(1,indxg,i,j,k,cons,intenergy,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}
__device__ void x1_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	//outflow_boundary_outer(1,indxg,i,j,k,cons,intenergy,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}

//extern "C" {
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

    cond = params->kappa;
    real w = params->width;
    real tinit = params->tinit;



    real *x1 = grid->xc1;

    real *rho       = &grid->cons[0*ntot];
    real *mx1       = &grid->cons[1*ntot];
    real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 

    real gamma = params->gamma ;
    real delad = 1. - 1./gamma;

    real ke;

    real u1 = 0;
    real u2 = 0;
    real u3 = 0;
    real temp;
    real chi = cond / (1 - 1./gamma);
    real t = .1;
    real pres = 1.;



    for(k=-NGHX3;k<nx3+NGHX3;k++) {
		for(j=-NGHX2;j<nx2+NGHX2;j++) {
			for(i=-NGHX1;i<nx1+NGHX1;i++) {
				indx = INDEX(i,j,k);


				// Set gaussian temp ic.
				temp = (tinit-1)*exp(-x1[i]*x1[i]/(4*chi*t))/sqrt(4*M_PI*chi*t) + 1.;

				rho[indx] = pres/(temp*delad);

				mx1[indx] = u1*rho[indx];
				mx2[indx] = u2*rho[indx];
				mx3[indx] = u3*rho[indx];

				ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
				ke /= 2*rho[indx];
				intenergy[indx] = temp *rho[indx]/ gamma;
				energy[indx] = intenergy[indx] + ke;
				for(n=5;n<nf;n++) {
					grid->cons[n*ntot+indx] = 0;
				}





			}
		}
    }
    return;


}
//}
