/* Rayleigh-Taylor */
#include "defs.h"
#include <time.h>
#include "cuda_defs.h"


__device__ __managed__ real g_param = .1;
__device__ __managed__ real nu = .01;


__device__ real kinematic_viscosity(real x1, real x2, real x3) {
	return nu;
}

__device__ real gravpot(real x, real y, real z) {
    return g_param*y;
}

__device__ void hydrostatic_lower(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    int n;
    int indx = GINDEX(i,-j-1,k);


//    inten = inten0 - (x2[j]-x2[0])*cons[indx0+0*ntot]*.1/(g-1);




    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + 1*ntot] = cons[indx + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx + 3*ntot];
    cons[indxg + 4*ntot] = cons[indx + 4*ntot];

    cons[indxg + 4*ntot] -= g_param*cons[indx + 0*ntot]/(g) * (x2[j] - x2[-j-1]);


    for(n=5;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}
__device__ void hydrostatic_upper(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    int n;
    int indx = GINDEX(i,nx2+ -(j-nx2)-1,k);

//    inten = inten0 - (x2[j]-x2[nx2-1])*cons[indx0+0*ntot]*.1/(g-1);

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + 1*ntot] = cons[indx + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx + 3*ntot];
    cons[indxg + 4*ntot] = cons[indx + 4*ntot];

    cons[indxg + 4*ntot] -= g_param*cons[indx + 0*ntot]/(g) * (x2[j] - x2[nx2+-(j-nx2)-1]);



    for(n=5;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}
__device__ void x2_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	hydrostatic_lower(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}
__device__ void x2_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	hydrostatic_upper(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}
__device__ void x1_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	periodic_boundary_inner(1,indxg,i,j,k,cons,intenergy,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}
__device__ void x1_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	periodic_boundary_outer(1,indxg,i,j,k,cons,intenergy,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
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

    g_param = params->g;
    nu = params->nu;

    real *x1 = grid->xc1;

	real *xm2 = grid->xm2;
	real *x2 = grid->xc2;


	real *rho       = &grid->cons[0*ntot];
	real *mx1       = &grid->cons[1*ntot];
	real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 

    real gamma = params->gamma;
    real gamma_1 = gamma-1;
    real ke;

    real u1 = 0;
    real u2 = 0;
    real u3 = 0;
    real pres = 1.;
    real amp = params->amp;



    for(k=-NGHX3;k<nx3+NGHX3;k++) {
		for(j=-NGHX2;j<nx2+NGHX2;j++) {
			for(i=-NGHX1;i<nx1+NGHX1;i++) {
				indx = INDEX(i,j,k);

				u2 = 0.;

				/*
				if (xm2[j] >= .01*cos(2*M_PI*x1[i])) {
					pres = 1./gamma - .1*2*(x2[j]- 0.0);
					rho[indx] = 2.;
				}
				else {
					pres = 1./gamma - .1*1*(x2[j]- 0.0);
					rho[indx] = 1.;
				}
				*/

				if (fabs(xm2[j]) <=1e-8) {
					u2 = amp*cos(2*M_PI*x1[i]);
				}



				if (xm2[j] >= 0) {

					pres = 1./gamma - g_param*2*(x2[j]- 0.0);
					rho[indx] = 2.;
				}
				else {
					pres = 1./gamma - g_param*1*(x2[j]- 0.0);
					rho[indx] = 1.;
				}


				mx1[indx] = u1*rho[indx];
				mx2[indx] = u2*rho[indx];
				mx3[indx] = u3*rho[indx];

				ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
				ke /= 2*rho[indx];
				intenergy[indx] = pres/gamma_1;
				energy[indx] = intenergy[indx] + ke;
				for(n=5;n<nf;n++) {
					grid->cons[n*ntot+indx] = 0;
            	}

			}
		}
	}
	return;


	}

