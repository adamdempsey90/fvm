typedef real (*pfunc)(real rho, real vx1, real vx2, real vx3, real pres,
		real x1, real x2, real x3,
		real dx1, real dx2, real dx3,
		real gamma, real time);


typedef struct Avg1D {
    int avgx1,avgx2,avgx3;
    int nx1,nx2;
    real *res;
    pfunc func;
    char *name;
} Avg1D;

typedef struct Avg2D {
    int avgx1,avgx2,avgx3;
    int nx1;
    real *res;
    pfunc avg_func;
    char *name;
} Avg2D;

typedef struct Snap1D {
    int nx1;
    real x1_lims[2];
    real x2_lims[2];
    real x3_lims[2];
    real *res;
    pfunc func;
    char *name;
} Snap1D;

typedef struct Snap2D {
    int nx1,nx2;
    real x1_lims[2];
    real x2_lims[2];
    real x3_lims[2];
    real *res;
    pfunc snap_func;
    char *name;
} Snap2D;

typedef struct Snap3D {
    int nx1,nx2,nx3;
    real x1_lims[2];
    real x2_lims[2];
    real x3_lims[2];
    real *res;
    pfunc snap_func;
    char *name;
} Snap3D;




void enroll_1D_average(char *name, pfunc func, int avg_dims[3], Avg1D *avg_arr, int *count);
void enroll_1D_snap(char *name, pfunc func, Snap1D *avg_arr, int *count);
void free_1D_averages(Avg1D *avg_arr, int count);
void free_1D_snaps(Snap1D *avg_arr, int count);
void compute_1D_averages(Avg1D *avg_arr, int count,
		real *cons, real *intenergy,
		real *x1, real *x2, real *x3,
		real *dx1, real *dx3, real *dx3,
		GridCons *grid, Parameters *params);
void output_1D_averages(char *fname, Avg1D *avg_arr, int count, GridCons *grid, Parameters *params);
void compute_1D_snapshots(Snap1D *snap_arr, int count,
		real *cons, real *intenergy,
		real *x1, real *x2, real *x3,
		real *dx1, real *dx3, real *dx3,
		GridCons *grid, Parameters *params);
void output_1D_snapshots(char *fname,Snap1D *snap_arr, int count, GridCons *grid, Parameters *params);
