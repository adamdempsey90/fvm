#include "defs.h"
#include "cuda_defs.h"


int main(int argc, char *argv[]) {
    int Nout, step;



    GridCons *grid = (GridCons *)malloc(sizeof(GridCons));
    FluxCons *fluxes = (FluxCons *)malloc(sizeof(FluxCons));
    Parameters *params = (Parameters *)malloc(sizeof(Parameters));

    printf("Reading pars\n");
    read_pars(params,argc,argv);
    printf("Allocating\n");
    allocate(grid,fluxes,params);
    printf("Init mesh\n");
    init_mesh(grid,params);
    printf("Init gas\n");
    init_gas(grid,params);

    printf("Outputting results to %s\n",params->outputname);


    size_t totbytes = sizeof(real)*(7*(grid->ntot)*(grid->nf) // cons, flux, UL/R
            + grid->ntot // intenergy 
            + grid->ntot // d_half
            + grid->size_x[0] 
            + grid->size_x[1]
            + 1024);
    printf("%.2f GB will be used on device\n",totbytes/(real)1e9);



    Nout = params->Nout;
    step = 0;
    real dtout = params->dtout;

#ifndef PROF
    output(step,grid,fluxes,params);
#endif

    if (params->one_step) {
        printf("Executing one step then exiting.\n");
        algogas_onestep(params->tend,grid,fluxes,params);
        step += 1;
        printf("Output %d at t=%.2e\n",step,grid->time);
#ifndef PROF
        output(step,grid,fluxes,params); // Output 
#endif
    }

    else {
        printf("Evolving from %.2e to %.2e\n",grid->time,params->tend);

        for(step=1;step<=Nout;step++) {
            algogas_dt(dtout,grid,fluxes,params); // Evolve for a time of dtout
            printf("Output %d at t=%.2e\n",step,grid->time);
#ifndef PROF
            output(step,grid,fluxes,params); // Output 
#endif
        }
    }




    printf("Exiting.\n");

    return 0;

}


