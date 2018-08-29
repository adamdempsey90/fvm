#include "defs.h"


void allocate(GridCons *grid,FluxCons *fluxes, Parameters *params);
void read_pars(Parameters *params);
int main(int argc, char *argv[]) {

    int i,j,k,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 

    GridCons *grid = (GridCons *)malloc(sizeof(GridCons));
    FluxCons *fluxes = (FluxCons *)malloc(sizeof(FluxCons));
    Parameters *params = (Parameters *)malloc(sizeof(Parameters));

    read_pars(params);
    allocate(grid,fluxes,params);
    init_mesh(grid,params);
    init_gas(grid,params);

    printf("IC\n");

    real *x1        = grid->xc1;
    real *x2        = grid->xc2;
    real *x3        = grid->xc3;
    real *rho       = &grid->cons[0*grid->ntot + grid->offset];
    real *mx1       = &grid->cons[1*grid->ntot + grid->offset];
    real *mx2       = &grid->cons[2*grid->ntot + grid->offset];
    real *mx3       = &grid->cons[3*grid->ntot + grid->offset];
    real *energy    = &grid->cons[4*grid->ntot + grid->offset];
    real *intenergy = &grid->cons[5*grid->ntot + grid->offset];

    nx1 = grid->nx[0];
    nx2 = grid->nx[1]; 
    nx3 = grid->nx[2];
    size_x1 = grid->size_x[0];
    size_x12 = grid->size_x12;


    for(k=0;k<nx3;k++) {
        for(j=0;j<nx2;j++) {
            for(i=0;i<nx1;i++) {
                indx = INDEX(i,j,k); 
                printf("%.5f\t%.3f\t%.3f\t%.3f\n",x1[i],
                        rho[indx], energy[indx],intenergy[indx]);
            }
        }
    }
    

    int Nout = params->Nout;
    int step = 0;
    real dtout = params->dtout;

    output(step,grid,fluxes,params);

    printf("Evolving from %.2e to %.2e\n",grid->time,params->tend);

    for(step=0;step<Nout;step++) {
        algogas(dtout,grid,fluxes,params); // Evolve for a time of dtout
        printf("Output %d at t=%.2e\n",step,grid->time);
        output(step,grid,fluxes,params); // Output 
        
    }




    printf("Exiting.\n");

    return 1;

}
void read_pars(Parameters *params) {
    params->nx1 = 100;
    params->nx2 = 1;
    params->nx3 = 1;
    params->gamma = 1.4;
    params->cfl = .2;

    params->x1_min = 0.;
    params->x1_max = 1.;
    params->x2_min = 0.;
    params->x2_max = 1.;
    params->x3_min = 0.;
    params->x3_max = 1.;


    params->gamma_1 = params->gamma-1.;
    params->gamma_c = params->gamma * params->gamma_1;

    params->tend = .2;
    params->Nout = 2;
    params->dtout = (params->tend)/(float)params->Nout;

    strcpy(params->outputname ,"out/test_out");
    return;

}


void allocate(GridCons *grid, FluxCons *fluxes,Parameters *params) {
    int nx1,nx2,nx3, ntot, nf, nscalars;
    int size_x1, size_x2, size_x3;


    nx1 = params->nx1;
    nx2 = params->nx2;
    nx3 = params->nx3;
    size_x1 = nx1 + 2*NGHX1; 
    size_x2 = nx2 + 2*NGHX2; 
    size_x3 = nx3 + 2*NGHX3; 
    ntot = size_x1*size_x2*size_x3;

    nscalars = params->nscalars;
    nf = 5 + 1 + nscalars; // rho/mx1/mx2/mx3/energy + internal energy + num scalars

    grid->nscalars = nscalars;
    grid->nf = nf;

    grid->nx[0] = nx1;
    grid->nx[1] = nx2;
    grid->nx[2] = nx3;
    grid->size_x[0] = size_x1;
    grid->size_x[1] = size_x2;
    grid->size_x[2] = size_x3;
    grid->size_x12 = size_x1*size_x2;
    grid->ntot = ntot;
    grid->offset = NGHX1 + size_x1 * NGHX2 + size_x1*size_x2 *NGHX3;

    /* 1D arrays */
    grid->xm1 = (real *)malloc(sizeof(real)*(size_x1 + 1));
    grid->xm2 = (real *)malloc(sizeof(real)*(size_x2 + 1));
    grid->xm3 = (real *)malloc(sizeof(real)*(size_x3 + 1));
    grid->xc1 = (real *)malloc(sizeof(real)*size_x1);
    grid->xc2 = (real *)malloc(sizeof(real)*size_x2);
    grid->xc3 = (real *)malloc(sizeof(real)*size_x3);
    grid->dx1 = (real *)malloc(sizeof(real)*size_x1);
    grid->dx2 = (real *)malloc(sizeof(real)*size_x2);
    grid->dx3 = (real *)malloc(sizeof(real)*size_x3);


    grid->xm1 = &grid->xm1[NGHX1];
    grid->xm2 = &grid->xm2[NGHX2];
    grid->xm3 = &grid->xm3[NGHX3];
    grid->xc1 = &grid->xc1[NGHX1];
    grid->xc2 = &grid->xc2[NGHX2];
    grid->xc3 = &grid->xc3[NGHX3];
    grid->dx1 = &grid->dx1[NGHX1];
    grid->dx2 = &grid->dx2[NGHX2];
    grid->dx3 = &grid->dx3[NGHX3];

    /* 3D arrays */
    grid->hfac  = (real *)malloc(sizeof(real)*ntot*3);
    grid->cons = (real *)malloc(sizeof(real)*ntot*nf);


    grid->time = 0.;

    fluxes->Ustar_1 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->Ustar_2 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->Ustar_3 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->Fstar_1 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->Fstar_2 = (real *)malloc(sizeof(real)*ntot*nf);
    fluxes->Fstar_3 = (real *)malloc(sizeof(real)*ntot*nf);
    return;

}
