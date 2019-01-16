typedef real (*pfunc)(real rho, real vx1, real vx2, real vx3, real pres,
		real x1, real x2, real x3,
		real dx1, real dx2, real dx3,
		real gamma, real time);


typedef struct Avg1D {
    int avgx1,avgx2,avgx3;
    int nx1,nx2;
    real *res;
    pfunc avg_func;
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
    pfunc snap_func;
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
