#include "defs.h"
#include "cuda_defs.h"
#include <ctype.h>

#define MAXBLOCKS 1024
#define MINBLOCKS 32
#define MAXTHREADS 256
void default_pars(Parameters *params);

void set_threads_blocks(int n, int *threads, int *blocks) {
	*threads = MAXTHREADS;
	*blocks = (n+*threads-1)/(*threads);
    if (*blocks > MAXBLOCKS) *blocks = MAXBLOCKS;
	return;
    /*
	if (n <= MAXTHREADS) {
		*threads = MAXTHREADS;
		*blocks = 1;
		return;
	}
	for(*threads=MAXTHREADS; *threads >= 64; *threads /= 2) {
		*blocks = min((n+*threads-1)/(*threads),MAXBLOCKS);
		if (*blocks >= MINBLOCKS) return;

	}
	printf("Can't determine threads/blocks from problem size %d\n",n);
	return;
    */

}

int main(int argc, char *argv[]) {
    int Nout, step;
    int threads, blocks;
    int restart = FALSE;
    int onestep = FALSE;
    char parfile[512];
    char restartfile[512];



    GridCons *grid = (GridCons *)malloc(sizeof(GridCons));
    FluxCons *fluxes = (FluxCons *)malloc(sizeof(FluxCons));
    Parameters *params = (Parameters *)malloc(sizeof(Parameters));
    
    /* Set default parameters first */
    default_pars(params);

    /* Check for command line arguments */
    if (argc > 1) {
        strcpy(parfile,argv[argc-1]);
        argc--;
        printf("Parameter file is %s\n",parfile);

        /* Point to first argument */
        argc--;
        argv = &argv[1];


        while(argc > 0) {
            /* Test for restart and one-step flags */
            if((*argv[0] == '-') && (*(argv[0]+1) != '\0') && (*(argv[0]+2) == '\0')) {
                if ( (char)tolower(*(argv[0]+1)) == 'r' ){
                    /* Restart in form
                     * -R restartfile
                     */
                    restart = TRUE;
                    strcpy(restartfile,argv[1]);
                    argc -= 2;
                    argv = &argv[2];
                }
                else if ( (char)tolower(*(argv[0]+1)) == 's' ){
                    /* When -S flag is present we just 
                     * take one step and exit
                     */
                    onestep = TRUE;
                    argc -= 1;
                    argv = &argv[1];
                }
                else {
                    printf("Unknown command line argument %s\n",argv[0]);
                    argc=0;
                }
            }
            else {
                /* No -? commands but argc > 0 
                 * Assume this is redefined arguments 
                 */
                break;
            }

        }

        /* Read parameter file */
        read_param_file(parfile,argc,argv,params);

    }
    else {
    	printf("No parameter file given, assuming in.par\n");
    	strcpy(parfile,"in.par");
    	read_param_file(parfile,argc,argv,params);
    }

#ifdef DIMS3
    printf("Compiled with 3D\n");
#endif
#ifdef HLLC
    printf("Compiled with HLLC\n");
#endif
#ifdef EXACT
    printf("Compiled with EXACT\n");
#endif
#ifdef POTENTIAL
    printf("Compiled with POTENTIAL\n");
#endif
#ifdef CONDUCTION
    printf("Compiled with CONDUCTION\n");
#endif
#ifdef VISCOSITY
    printf("Compiled with VISCOSITY\n");
#endif
    allocate(grid,fluxes,params);
    grid->time = 0.;

    init_mesh(grid,params);

    if (restart) {
    	read_restart(restartfile,grid,fluxes,params);
    }
    else {
    	init_gas(grid,params);
    }

    Nout = params->nout;
    step = 0;
    real dtout = (params->tend - grid->time)/(float)Nout;
    real dt_curr;
    
    set_threads_blocks(grid->ntot,&threads,&blocks);
#ifndef SILENT
    printf("Outputting results to %s\n",params->outputname);

    size_t totbytes = sizeof(real)*(7*(grid->ntot)*(grid->nf) // cons, flux, UL/R
            + grid->ntot // intenergy 
            + grid->ntot // d_half
            + grid->ntot*3 // scale factors
            + grid->size_x1*3 +1 // xm, xc, dx
            + grid->size_x2*3 +1
            + grid->size_x3*3 +1
            + 1024); // Shared mem block
    printf("%.2f GB will be used on device\n",totbytes/(real)1e9);
    printf("Threads %d, Blocks %d\n",threads,blocks);
#endif

#ifndef PROF
    output(step,grid,fluxes,params);
#endif
    
   dt_curr = algogas_firststep(params->tend,threads, blocks,grid, fluxes,params);

    if (onestep) {
#ifndef SILENT
        printf("Executing one step then exiting.\n");
#endif
        step += 1;
#ifndef SILENT
        printf("Output %d at t=%.2e\n",step,grid->time);
#endif
#ifndef PROF
        output(step,grid,fluxes,params); // Output 
#endif
    }
    else {
#ifndef SILENT
        printf("Evolving from %.2e to %.2e\n",grid->time,params->tend);
#endif
        for(step=1;step<=Nout;step++) {
            // Evolve for a time of dtout
            dt_curr = algogas_dt(dt_curr, dtout,threads, blocks,grid,fluxes,params); 
#ifndef SILENT
            printf("Output %d at t=%.2e\n",step,grid->time);
#endif
#ifndef PROF
            output(step,grid,fluxes,params); // Output 
#endif
        }
    }



#ifndef SILENT
    printf("Exiting.\n");
#endif

    return 0;

}

void default_pars(Parameters *params) {
    params->nx1 = 128;
    params->nx2 = 1;
    params->nx3 = 1;
    params->nscalars = 0;
    params->gamma = 1.4;
    params->cfl = .2;

    params->x1_min = 0.;
    params->x1_max = 1;

    params->x2_min = 0.;
    params->x2_max = 1;

    params->x3_min = 0.;
    params->x3_max = 1;


    params->tend = .3;
    params->nout = 10;


    strcpy(params->outputname ,"out/default");

}

