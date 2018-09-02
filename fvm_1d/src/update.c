#include "defs.h"


void update_cons(GridCons *grid, FluxCons *fluxes,Parameters *params, real dt) {
    int i,indx,indxp1,indxp2,indxp3;
    int nx1,nx2,nx3;
    int n,nf,ntot;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;

    real dtdx1,dtdx2,dtdx3;
    real *dx1 = grid->dx1;

    real *h1 = &grid->hfac[0*grid->ntot];
    real *h2 = &grid->hfac[1*grid->ntot];
    real *h3 = &grid->hfac[2*grid->ntot];

    real *cons = grid->cons;
    real *F_1 = fluxes->Fstar_1;
    real dU1;
    int indxm1;
    for(n=0;n<nf;n++) {
        for(i=0;i<nx1;i++) {
            indx = INDEX(i) + n*ntot;
            indxm1 = INDEX(i-1) + n*ntot;
            dtdx1 = dt/dx1[i];
            cons[indx] += dtdx1*(F_1[indxm1]-F_1[indx]);
        }
    }

    return;



}
