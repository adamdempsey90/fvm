#include "defs.h"
#include "cuda_defs.h"
#include <ctype.h>
#include <time.h>


void cons_to_prim_grid(GridCons *grid, Parameters *params);
void prim_to_cons_grid(GridCons *grid, Parameters *params);


int main(int argc, char *argv[]) {
    struct timespec tic, toc; 
    double elapsed_sec;
    int Nout, step;
    int threads, blocks;
    int restart = FALSE;
    int onestep = FALSE;
    int nostep = FALSE;
    char parfile[512];
    char restartfile[512];

    

    clock_gettime(CLOCK_MONOTONIC, &tic);


    GridCons *grid = (GridCons *)malloc(sizeof(GridCons));
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
                else if ( (char)tolower(*(argv[0]+1)) == 'b' ){
                    /* When -B flag is present we just
                     * set boundary conditions and exit
                     */
                    nostep = TRUE;
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
    fflush(stdout);
    allocate(grid,params);
    grid->time = 0.;

    init_mesh(grid,params);

    if (restart) {
    	read_restart(restartfile,grid,params);
    	prim_to_cons_grid(grid,params);
    }
    else {
    	init_gas(grid,params);
    	cons_to_prim_grid(grid,params);


    }



    config_kernels(&threads,&blocks,grid,params);


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
    fflush(stdout);

#endif
    
    if (onestep) {
        params->maxsteps = 0;
#ifndef SILENT
        printf("Executing one step then exiting.\n");
        fflush(stdout);
#endif
    }

    driver(grid,params);



    clock_gettime(CLOCK_MONOTONIC, &toc);
    elapsed_sec = (double)(toc.tv_sec - tic.tv_sec) + 1e-9*(toc.tv_nsec-tic.tv_nsec);

#ifndef SILENT
    int hrs,mins;
    double secs;
    printf("Total execution time of ");
    if (elapsed_sec < 60.) {
        printf("%.2fs\n",elapsed_sec);
    }
    else if (elapsed_sec< 3600.) {
        mins = (int)(elapsed_sec/60.);
        secs = elapsed_sec - 60.*mins;
        printf("%dm%.2fs\n",mins,secs);
    }
    else {
        hrs = (int)(elapsed_sec/3600.);
        mins = (int)((elapsed_sec - 3600*hrs)/60.);
        secs = elapsed_sec - 60.*mins - 3600*hrs;
        printf("%dh%dm%.2fs\n",hrs,mins,secs);
    }
    printf("Exiting.\n");
    fflush(stdout);
#endif


    return 0;

}


void cons_to_prim_grid(GridCons *grid, Parameters *params) {
	int indx,n;
	int ntot = grid->ntot;
	int nf = grid->nf;
	int offset = grid->offset;
	real *cons = grid->cons;
	real *prim = grid->prim;
	real *intenergy = grid->intenergy;
	real g1 = params->gamma-1;
	for(indx=-offset;indx<ntot-offset;indx++) {
		prim[indx] =cons[indx];
		prim[indx + 1*ntot] = cons[indx + 1*ntot]/cons[indx];
		prim[indx + 2*ntot] = cons[indx + 2*ntot]/cons[indx];
		prim[indx + 3*ntot] = cons[indx + 3*ntot]/cons[indx];
		prim[indx + 4*ntot] = intenergy[indx] * g1;
		for(n=5;n<nf;n++) prim[indx + n*ntot] = cons[indx + n*ntot]/cons[indx];
	}
	return;
}
void prim_to_cons_grid(GridCons *grid, Parameters *params) {
	int indx,n;
	int ntot = grid->ntot;
	int nf = grid->nf;
	int offset = grid->offset;
	real ke;
	real *cons = grid->cons;
	real *prim = grid->prim;
	real *intenergy = grid->intenergy;
	real g1 = params->gamma-1;
	for(indx=-offset;indx<ntot-offset;indx++) {
		cons[indx] = prim[indx];
		cons[indx + 1*ntot] = prim[indx + 1*ntot]*prim[indx];
		cons[indx + 2*ntot] = prim[indx + 2*ntot]*prim[indx];
		cons[indx + 3*ntot] = prim[indx + 3*ntot]*prim[indx];
		ke = .5*prim[indx]*( prim[indx + 1*ntot]*prim[indx + 1*ntot] +
				             prim[indx + 2*ntot]*prim[indx + 2*ntot] +
				             prim[indx + 3*ntot]*prim[indx + 3*ntot]);
		intenergy[indx] = prim[indx + 4*ntot]/g1;
		cons[indx + 4*ntot] = intenergy[indx] + ke;
		for(n=5;n<nf;n++) cons[indx + n*ntot] = prim[indx + n*ntot]*prim[indx];
	}
	return;
}
