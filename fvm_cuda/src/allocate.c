#include "defs.h"
extern void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3);
void allocate(GridCons *grid, FluxCons *fluxes,Parameters *params) {
    int nx1,nx2,nx3, ntot, nf, nscalars;
    int size_x1, size_x2, size_x3, size_x12;


    size_x1 = params->nx1;
    size_x2 = params->nx2;
    size_x3 = params->nx3;
    nx1 = size_x1- 2*NGHX1;
    nx2 = size_x2- 2*NGHX2;
    nx3 = size_x3- 2*NGHX3;


    size_x12 = size_x1*size_x2;
    ntot = size_x1*size_x2*size_x3;

    nscalars = params->nscalars;
    nf = 5 + nscalars; // rho/mx1/mx2/mx3/energy  + num scalars

    grid->nscalars = nscalars;
    grid->nf = nf;

    grid->nx[0] = nx1;
    grid->nx[1] = nx2;
    grid->nx[2] = nx3;
    grid->size_x1 = size_x1;
    grid->size_x2 = size_x2;
    grid->size_x3 = size_x3;
    grid->size_x12 = size_x1*size_x2;
    grid->ntot = ntot;
    grid->offset = NGHX1 + size_x1*NGHX2  + size_x12*NGHX3;

    /* 1D arrays */
    grid->xm1 = (real *)malloc(sizeof(real)*(size_x1 + 1));
    grid->xc1 = (real *)malloc(sizeof(real)*size_x1);
    grid->dx1 = (real *)malloc(sizeof(real)*size_x1);

    grid->xm2 = (real *)malloc(sizeof(real)*(size_x2 + 1));
    grid->xc2 = (real *)malloc(sizeof(real)*size_x2);
    grid->dx2 = (real *)malloc(sizeof(real)*size_x2);

    grid->xm3 = (real *)malloc(sizeof(real)*(size_x3 + 1));
    grid->xc3 = (real *)malloc(sizeof(real)*size_x3);
    grid->dx3 = (real *)malloc(sizeof(real)*size_x3);

    grid->xm1 = &grid->xm1[NGHX1];
    grid->xc1 = &grid->xc1[NGHX1];
    grid->dx1 = &grid->dx1[NGHX1];

    grid->xm2 = &grid->xm2[NGHX2];
    grid->xc2 = &grid->xc2[NGHX2];
    grid->dx2 = &grid->dx2[NGHX2];

    grid->xm3 = &grid->xm3[NGHX3];
    grid->xc3 = &grid->xc3[NGHX3];
    grid->dx3 = &grid->dx3[NGHX3];

    /* 3D arrays */
    grid->hfac  = (real *)malloc(sizeof(real)*ntot*3);
    grid->cons = (real *)malloc(sizeof(real)*ntot*nf);
    grid->intenergy  = (real *)malloc(sizeof(real)*ntot);

    grid->hfac = &grid->hfac[grid->offset];
    grid->cons = &grid->cons[grid->offset];
    grid->intenergy = &grid->intenergy[grid->offset];



    fluxes->UL_1 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->UR_1 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->Fstar_1 = (real *)malloc(sizeof(real)*ntot*nf);

    fluxes->UL_1 = &fluxes->UL_1[grid->offset];
    fluxes->UR_1 = &fluxes->UR_1[grid->offset];
    fluxes->Fstar_1 = &fluxes->Fstar_1[grid->offset];

    fluxes->UL_2 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->UR_2 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->Fstar_2 = (real *)malloc(sizeof(real)*ntot*nf);

    fluxes->UL_2 = &fluxes->UL_2[grid->offset];
    fluxes->UR_2 = &fluxes->UR_2[grid->offset];
    fluxes->Fstar_2 = &fluxes->Fstar_2[grid->offset];


    fluxes->UL_3 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->UR_3 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->Fstar_3 = (real *)malloc(sizeof(real)*ntot*nf);

    fluxes->UL_3 = &fluxes->UL_3[grid->offset];
    fluxes->UR_3 = &fluxes->UR_3[grid->offset];
    fluxes->Fstar_3 = &fluxes->Fstar_3[grid->offset];
    return;

}

void init_uniform_mesh(GridCons *grid, Parameters *params) {

    int i,j,k,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12;
    int n;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    nx3 = grid->nx[2];
    size_x1 = grid->size_x1;
    size_x12 = grid->size_x12;
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

    /* x3 direction */
	n = nx3;
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


    real *h1 = &grid->hfac[0*grid->ntot];
    real *h2 = &grid->hfac[1*grid->ntot];
    real *h3 = &grid->hfac[2*grid->ntot];

    for(k=-NGHX3;k<nx3+NGHX3;k++) {
		for(j=-NGHX2;j<nx2+NGHX2;j++) {
			for(i=-NGHX1;i<nx1+NGHX1;i++) {
			 indx = INDEX(i,j,k);
			 scale_factors(grid->xc1[i], grid->xc2[j],grid->xc3[k],
					&h1[indx],&h2[indx],&h3[indx]);

			}
		}
    }


    return;

}
