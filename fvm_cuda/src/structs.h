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

	int gridSize_update_cons, blockSize_update_cons;
	int gridSize_transverse, blockSize_transverse;

	int gridSize_riemann , blockSize_riemann ;
	int gridSize_reconstruct, blockSize_reconstruct;
	int gridSize_reduc, blockSize_reduc;
	int gridSize_boundary, blockSize_boundary;

#ifdef POTENTIAL
	int gridSize_source_transverse, blockSize_source_transverse;
	int gridSize_source, blockSize_source;
	int gridSize_update_source, blockSize_update_source;
#endif
#ifdef VISCOSITY
	int gridSize_viscous_flux, blockSize_viscous_flux;
	int gridSize_divergence , blockSize_divergence ;
#endif
#ifdef CONDUCTION
	int gridSize_conduction_flux, blockSize_conduction_flux;
#endif

} GridCons;
