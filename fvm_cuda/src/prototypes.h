
void init_gas(GridCons *grid, Parameters *params); 
void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3);
void init_mesh(GridCons *grid, Parameters *params);
extern "C" {
real algogas_dt(real dt,real dtout, int threads, int blocks, GridCons *grid, FluxCons *fluxes, Parameters *params);
real algogas_firststep(real dtout, int threads, int blocks, GridCons *grid, FluxCons *fluxes, Parameters *params);
void output(int step, GridCons *grid, FluxCons *fluxes, Parameters *params);
void allocate(GridCons *grid,FluxCons *fluxes, Parameters *params);
void read_pars(Parameters *params, int argc, char *argv[]);
}





