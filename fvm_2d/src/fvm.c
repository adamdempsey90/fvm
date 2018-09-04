#include "defs.h"


int main(int argc, char *argv[]) {
    int Nout, step;



    GridCons *grid = (GridCons *)malloc(sizeof(GridCons));
    FluxCons *fluxes = (FluxCons *)malloc(sizeof(FluxCons));
    Parameters *params = (Parameters *)malloc(sizeof(Parameters));

    read_pars(params,argc,argv);
    allocate(grid,fluxes,params);
    init_mesh(grid,params);
    init_gas(grid,params);




    Nout = params->Nout;
    step = 0;
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


