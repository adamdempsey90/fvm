#include "defs.h"


void update_cons(GridCons *grid, FluxCons *fluxes,Parameters *params, real dt) {
    int i,j,indx;
    int nx1,nx2,nx3;
    int n,nf,ntot;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;

    real dtdx1,dtdx2,dtdx3;
    real *dx1 = grid->dx1;
    real *dx2 = grid->dx2;

    real *h1 = &grid->hfac[0*grid->ntot];
    real *h2 = &grid->hfac[1*grid->ntot];
    real *h3 = &grid->hfac[2*grid->ntot];

    real *cons = grid->cons;
    real *intenergy = grid->intenergy;
    real *F_1 = fluxes->Fstar_1;
    real *F_2 = fluxes->Fstar_2;
    real dU1;
    int indxm1, indxm2;
#pragma omp parallel for private(indx,indxm1,indxm2,dtdx1,dtdx2);
    for(n=0;n<nf;n++) {
        for(j=0;j<nx2;j++) {
            for(i=0;i<nx1;i++) {
                indx = INDEX(i,j) + n*ntot;
                indxm1 = INDEX(i-1,j) + n*ntot;
                indxm2 = INDEX(i,j-1) + n*ntot;
                dtdx1 = dt/dx1[i];
                dtdx2 = dt/dx2[j];
                cons[indx] += dtdx1*(F_1[indxm1]-F_1[indx]) + 
                              dtdx2*(F_2[indxm2]-F_2[indx]);
            }
        }
    }
    /* Sync internal energy */
    real ke;
#pragma omp parallel for private(indx);
    for(j=0;j<nx2;j++) {
        for(i=0;i<nx1;i++) {
            indx = INDEX(i,j);
            intenergy[indx] = cons[indx+4*ntot] - .5*(
                    cons[indx + ntot]*cons[indx + ntot] +
                    cons[indx + 2*ntot]*cons[indx + 2*ntot] +
                    cons[indx + 3*ntot]*cons[indx + 3*ntot])/cons[indx];
        }
    }

    return;



}
void transverse_update(GridCons *grid, FluxCons *fluxes,Parameters *params, real dt) {
    int i,j,indx;
    int nx1,nx2,nx3;
    int n,nf,ntot;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;

    real dtdx1,dtdx2,dtdx3;
    real *dx1 = grid->dx1;
    real *dx2 = grid->dx2;

    real *h1 = &grid->hfac[0*grid->ntot];
    real *h2 = &grid->hfac[1*grid->ntot];
    real *h3 = &grid->hfac[2*grid->ntot];


    real *UL_1  = fluxes->UL_1;
    real *UR_1  = fluxes->UR_1;
    real *F_1 = fluxes->Fstar_1;
    real *UL_2  = fluxes->UL_2;
    real *UR_2  = fluxes->UR_2;
    real *F_2 = fluxes->Fstar_2;

    int indxm1, indxm2;
    int indxp1,indxp2, indxp1m2,indxp2m1;
    /* X1 direction */
#pragma omp parallel for private(indx,indxm2,indxp2,indxp1m2,dtdx2);           
    for(n=0;n<nf;n++) {
        for(j=-1;j<nx2+1;j++) {
            for(i=-1;i<nx1;i++) {
                indx = INDEX(i,j) + n*ntot;
                indxm2 = INDEX(i,j-1) + n*ntot;
                indxp2 = INDEX(i,j+1) + n*ntot;
                indxp1m2 = INDEX(i+1,j-1) + n*ntot;

                dtdx2 = .5*dt/dx2[j];
                UL_1[indx] += dtdx2*(F_2[indxm2]-F_2[indx]);
                UR_1[indx] += dtdx2*(F_2[indxp1m2]-F_2[indxp2]);
            }
        }
    }
    /* X2 direction */
#pragma omp parallel for private(indx,indxm1,indxp1,indxp2m1,dtdx1);           
    for(n=0;n<nf;n++) {
        for(j=-1;j<nx2;j++) {
            for(i=-1;i<nx1+1;i++) {
                indx = INDEX(i,j) + n*ntot;
                indxm1 = INDEX(i-1,j) + n*ntot;
                indxp1 = INDEX(i+1,j) + n*ntot;
                indxp2m1 = INDEX(i-1,j+1) + n*ntot;

                dtdx1 = .5*dt/dx1[i];
                UL_2[indx] += dtdx1*(F_1[indxm1]-F_1[indx]);
                UR_2[indx] += dtdx1*(F_1[indxp2m1]-F_1[indxp1]);
            }
        }
    }
  
    return;



}
