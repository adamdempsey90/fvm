
void init_gas(GridCons *grid, Parameters *params); 
void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3);
void init_mesh(GridCons *grid, Parameters *params);
void output(int step, GridCons *grid, FluxCons *fluxes, Parameters *params);
void algogas_single(real dt_max, GridCons *grid, FluxCons *fluxes, Parameters *params);
void algogas_dt(real dtout, GridCons *grid, FluxCons *fluxes, Parameters *params);
void reconstruct(GridCons *grid,FluxCons *fluxes,Parameters *params, real dt);
void riemann_fluxes(const real *UL, const real *UR, real *flux, int dir1,int nx[3],int size_x1, int size_x12,int nf,int ntot ,real gamma);
void update_cons(GridCons *grid, FluxCons *fluxes,Parameters *params, real dt);
void set_boundary_x3(GridCons *grid, Parameters *params);
void set_boundary_x2(GridCons *grid, Parameters *params);
void set_boundary_x1(GridCons *grid, Parameters *params);
void read_pars(Parameters *params);

void hll_flux(const real *UL, const real *UR, real *F,real g, real g1, real gp, int nf);

void anrs(const real *WL, const real *WR,real aL, real aR, real g, real *ps, real *us, real gamma_1,real gamma_2,real gamma_3,real gamma_4);
int exact_sample(const real *WL, const real *WR, const real ps, const real us, real *W, real g, real S, real tol);
int exact_star_region(const real *WL, const real *WR,  real g, real *ps, real *us, real tol);

void exact_flux(const real *UL, const real *UR, real *F,real g, real g1, int nf);
//void plm(real *cons_m, real *cons_c, real *cons_p, real *UL, real *UR, real dt, real *dx, real gamma_1,int dir1, int nf);
void plm(real *cons, real dt, real *dx, real gamma_1,int dir1, int nf);
