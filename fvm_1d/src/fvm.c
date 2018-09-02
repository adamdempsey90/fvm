#include "defs.h"


void allocate(GridCons *grid,FluxCons *fluxes, Parameters *params);
int main(int argc, char *argv[]) {
    int i,indx,ntot;
    int nx1;
    int size_x1; 



    GridCons *grid = (GridCons *)malloc(sizeof(GridCons));
    FluxCons *fluxes = (FluxCons *)malloc(sizeof(FluxCons));
    Parameters *params = (Parameters *)malloc(sizeof(Parameters));

    read_pars(params);
    allocate(grid,fluxes,params);
    init_mesh(grid,params);
    init_gas(grid,params);

    ntot = grid->ntot;

    printf("IC\n");

    real *x1        = grid->xc1;
    real *rho       = &grid->cons[0*ntot];
    real *mx1       = &grid->cons[1*ntot];
    real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy;

    nx1 = grid->nx[0];
    size_x1 = grid->size_x[0];


    

    int Nout = params->Nout;
    int step = 0;
    real dtout = params->dtout;

    output(step,grid,fluxes,params);

    if (params->one_step) {
        printf("Executing one step then exiting.\n");
        algogas_single(params->tend,grid,fluxes,params);
        step += 1;
        printf("Output %d at t=%.2e\n",step,grid->time);
        output(step,grid,fluxes,params); // Output 
    }

    else {
        printf("Evolving from %.2e to %.2e\n",grid->time,params->tend);

        for(step=1;step<=Nout;step++) {
            algogas_dt(dtout,grid,fluxes,params); // Evolve for a time of dtout
            printf("Output %d at t=%.2e\n",step,grid->time);
            output(step,grid,fluxes,params); // Output 
        }
    }




    printf("Exiting.\n");

    return 1;

}


void allocate(GridCons *grid, FluxCons *fluxes,Parameters *params) {
    int nx1,ntot, nf, nscalars;
    int size_x1;


    nx1 = params->nx1;
    size_x1 = nx1 + 2*NGHX1; 
    ntot = size_x1;

    nscalars = params->nscalars;
    nf = 5 + nscalars; // rho/mx1/mx2/mx3/energy  + num scalars

    grid->nscalars = nscalars;
    grid->nf = nf;

    grid->nx[0] = nx1;
    grid->size_x[0] = size_x1;
    grid->ntot = ntot;
    grid->offset = NGHX1; 

    /* 1D arrays */
    grid->xm1 = (real *)malloc(sizeof(real)*(size_x1 + 1));
    grid->xc1 = (real *)malloc(sizeof(real)*size_x1);
    grid->dx1 = (real *)malloc(sizeof(real)*size_x1);


    grid->xm1 = &grid->xm1[NGHX1];
    grid->xc1 = &grid->xc1[NGHX1];
    grid->dx1 = &grid->dx1[NGHX1];

    /* 3D arrays */
    grid->hfac  = (real *)malloc(sizeof(real)*ntot*3);
    grid->cons = (real *)malloc(sizeof(real)*ntot*nf);
    grid->intenergy  = (real *)malloc(sizeof(real)*ntot*3);

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
    return;

}
