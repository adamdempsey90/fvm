typedef struct GridCons {
    int nx[3]; // Does not include ghost zones
    int size_x[3]; // Includes ghost zones
    int size_x12;
    int ntot; // Total size with ghost_zones;
    int nf; // Number of fields in cons array, 5 + internal_energy + N-scalars
    int nscalars;
    int offset; // Index offset for cons
    real time;
    real *xm1;
    real *xm2;
    real *xm3;
    real *xc1;
    real *xc2;
    real *xc3;
    real *dx1;
    real *dx2;
    real *dx3;
    real *hfac;
    real *cons; 
    //real *rho;
    //real *mx1;
    //real *mx2;
    //real *mx3;
    //real *energy;
    //real *intenergy;

} GridCons;

typedef struct FluxCons {
    real *Ustar_1;
    real *Ustar_2;
    real *Ustar_3;
    real *Fstar_1;
    real *Fstar_2;
    real *Fstar_3;
} FluxCons;
typedef struct Parameters {
    int nx1,nx2,nx3, nscalars;
    real x1_min, x1_max;
    real x2_min, x2_max;
    real x3_min, x3_max;
    real gamma;
    real gamma_1;
    real gamma_c;
    real tend;
    real Nout;
    real dtout;
    real cfl;
    char outputname[512];
} Parameters;
