#include "defs.h"
void allocate(GridCons *grid, FluxCons *fluxes,Parameters *params) {
    int nx1,nx2,ntot, nf, nscalars;
    int size_x1, size_x2;


    nx1 = params->nx1;
    nx2 = params->nx2;
    size_x1 = nx1 + 2*NGHX1; 
    size_x2 = nx2 + 2*NGHX2;
    ntot = size_x1*size_x2;

    nscalars = params->nscalars;
    nf = 5 + nscalars; // rho/mx1/mx2/mx3/energy  + num scalars

    grid->nscalars = nscalars;
    grid->nf = nf;

    grid->nx[0] = nx1;
    grid->nx[1] = nx2;
    grid->size_x[0] = size_x1;
    grid->size_x[1] = size_x2;
    grid->ntot = ntot;
    grid->offset = NGHX1 + size_x1*NGHX2; 

    /* 1D arrays */
    grid->xm1 = (real *)malloc(sizeof(real)*(size_x1 + 1));
    grid->xc1 = (real *)malloc(sizeof(real)*size_x1);
    grid->dx1 = (real *)malloc(sizeof(real)*size_x1);

    grid->xm2 = (real *)malloc(sizeof(real)*(size_x2 + 1));
    grid->xc2 = (real *)malloc(sizeof(real)*size_x2);
    grid->dx2 = (real *)malloc(sizeof(real)*size_x2);

    grid->xm1 = &grid->xm1[NGHX1];
    grid->xc1 = &grid->xc1[NGHX1];
    grid->dx1 = &grid->dx1[NGHX1];

    grid->xm2 = &grid->xm2[NGHX2];
    grid->xc2 = &grid->xc2[NGHX2];
    grid->dx2 = &grid->dx2[NGHX2];

    /* 3D arrays */
    grid->hfac  = (real *)malloc(sizeof(real)*ntot*3);
    grid->cons = (real *)malloc(sizeof(real)*ntot*nf);
    grid->intenergy  = (real *)malloc(sizeof(real)*ntot);

    grid->hfac = &grid->hfac[grid->offset];
    grid->cons = &grid->cons[grid->offset];
    grid->intenergy = &grid->intenergy[grid->offset];

    grid->time = 0.;

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
    return;

}
