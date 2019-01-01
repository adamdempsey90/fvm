
void init_gas(GridCons *grid, Parameters *params); 
void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3);
void init_mesh(GridCons *grid, Parameters *params);
extern "C" {
real algogas_dt(real dt,real dtout, int threads, int blocks, GridCons *grid, Parameters *params);
real algogas_firststep(real dtout, int threads, int blocks, int restart, int nostep, GridCons *grid, Parameters *params);
void snapshot(char *name, GridCons *grid, Parameters *params);
void snapshot_1d(char *name, GridCons *grid, Parameters *params);
void snapshot_2d(char *name, GridCons *grid, Parameters *params);
void snapshot_3d(char *name, GridCons *grid, Parameters *params);
void read_restart(const char *fname, GridCons *grid, Parameters *params);
void allocate(GridCons *grid, Parameters *params);
void init_uniform_mesh(GridCons *grid, Parameters *params);
}





