
void init_gas(GridCons *grid, Parameters *params); 
void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3);
void init_mesh(GridCons *grid, Parameters *params);
extern "C" {
void algogas_dt(real dtout, GridCons *grid, FluxCons *fluxes, Parameters *params);
void algogas_onestep(real dtout, GridCons *grid, FluxCons *fluxes, Parameters *params);
void output(int step, GridCons *grid, FluxCons *fluxes, Parameters *params);
void allocate(GridCons *grid,FluxCons *fluxes, Parameters *params);
void read_pars(Parameters *params, int argc, char *argv[]);
}





