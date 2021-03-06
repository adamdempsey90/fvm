typedef struct GridCons {
    int nx[1]; // Does not include ghost zones
    int size_x[1]; // Includes ghost zones
    //int size_x12;
    int ntot; // Total size with ghost_zones;
    int nf; // Number of fields in cons array, 5 + internal_energy + N-scalars
    int nscalars;
    int offset; // Index offset for cons
    real time;
    real *xm1;
    real *xc1;
    real *dx1;
    real *hfac;
    real *cons; 
    real *intenergy;

} GridCons;

typedef struct FluxCons {
    real *UL_1;
    real *UR_1;
    real *Fstar_1;
} FluxCons;
typedef struct Parameters {
    int nx1, nscalars;
    real x1_min, x1_max;
    real gamma;
    real gamma_1;
    real gamma_c;
    real tend;
    real Nout;
    real dtout;
    real cfl;
    char outputname[512];
    int one_step;
} Parameters;
