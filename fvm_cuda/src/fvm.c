#include "defs.h"
#include "cuda_defs.h"


int main(int argc, char *argv[]) {
    int Nout, step;
    int restart = FALSE;



    GridCons *grid = (GridCons *)malloc(sizeof(GridCons));
    FluxCons *fluxes = (FluxCons *)malloc(sizeof(FluxCons));
    Parameters *params = (Parameters *)malloc(sizeof(Parameters));

    printf("Reading pars\n");
    read_pars(params,argc,argv);
    printf("Allocating\n");
    allocate(grid,fluxes,params);
    grid->time = 0.;
    printf("Init mesh\n");
    init_mesh(grid,params);

    if (restart) {
    	read_restart("out/restart.h5",grid,fluxes,params);
    }
    else {
    	printf("Init gas\n");
    	init_gas(grid,params);
    }

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
    real dtout = (params->tend - grid->time)/(float)params->Nout;
    real dt_curr;
    
    int threads = 256;
    int blocks = min((grid->ntot+threads-1)/threads,1024);
    printf("Threads %d, Blocks %d\n",threads,blocks);
    
#ifndef PROF
    output(step,grid,fluxes,params);
#endif
    
   dt_curr = algogas_firststep(params->tend,threads, blocks,grid, fluxes,params);

    if (params->one_step) {
        printf("Executing one step then exiting.\n");
        step += 1;
        printf("Output %d at t=%.2e\n",step,grid->time);
#ifndef PROF
        output(step,grid,fluxes,params); // Output 
#endif
    }
    else {
        printf("Evolving from %.2e to %.2e\n",grid->time,params->tend);

        for(step=1;step<=Nout;step++) {
            // Evolve for a time of dtout
            dt_curr = algogas_dt(dt_curr, dtout,threads, blocks,grid,fluxes,params); 
            printf("Output %d at t=%.2e\n",step,grid->time);
#ifndef PROF
            output(step,grid,fluxes,params); // Output 
#endif
        }
    }




    printf("Exiting.\n");

    return 0;

}


