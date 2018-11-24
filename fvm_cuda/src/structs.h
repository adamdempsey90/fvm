typedef struct GridCons {
    int nx[3]; // Does not include ghost zones
    int size_x1;
    int size_x2;
    int size_x3;// Includes ghost zones
    int size_x12;
    int ntot; // Total size with ghost_zones;
    int nf; // Number of fields in cons array, 5 + internal_energy + N-scalars
    int nscalars;
    int offset; // Index offset for cons
    real time;
    real *xm1;
    real *xc1;
    real *dx1;
    real *xm2;
    real *xc2;
    real *dx2;
    real *xm3;
    real *xc3;
    real *dx3;
    real *hfac;
    real *cons; 
    real *intenergy;
    real *prim;

} GridCons;
