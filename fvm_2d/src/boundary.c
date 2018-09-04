#include "defs.h"

void reflecting_boundary_x2_inner(GridCons *grid, Parameters *params) {
    /* Reflecting inner boundary.
     * This sets  dU/dx = 0 and u2 = 0 
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    il = -NGHX1; iu = nx1 + NGHX1;
    jl = -NGHX2; ju = 0;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /* -1 -> 0
             * -2 -> 1 
             * -3 -> 2
             *------------
             *    x | o
             *   x  |  o
             *  x   |   o
             *---------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(i,-j-1);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = -mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}
void reflecting_boundary_x1_inner(GridCons *grid, Parameters *params) {
    /* Reflecting inner boundary.
     * This sets  dU/dx = 0 and u2 = 0 
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    jl = -NGHX2; ju = nx2 + NGHX2;
    il = -NGHX1; iu = 0;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            indxg = INDEX(i,j);
            indx = INDEX((-i-1),j);

            rho[indxg] = rho[indx];
            mx1[indxg] = -mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}

void reflecting_boundary_x2_outer(GridCons *grid, Parameters *params) {
    /* Reflecting outer boundary.
     * This sets  dU/dx = 0 and u2 = 0 
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    il = -NGHX1; iu = nx1 + NGHX1;
    jl = nx2; ju = nx2+NGHX2;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /* nx -> nx-1     (0 -> -1)
             * nx+1 -> nx-2   (1 -> -2)
             * nx+2 -> nx-3   (2 -> -3)
             *------------
             *    o | x
             *   o  |  x
             *  o   |   x
             *---------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(i,nx2 + -(j-nx2)-1);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = -mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}
void reflecting_boundary_x1_outer(GridCons *grid, Parameters *params) {
    /* Reflecting inner boundary.
     * This sets  dU/dx = 0 and u2 = 0 
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    jl = -NGHX2; ju = nx2 + NGHX2;
    il = nx1 ; iu = nx1 + NGHX1;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            indxg = INDEX(i,j);
            indx = INDEX(nx1 + -(i-nx1)-1,j);

            rho[indxg] = rho[indx];
            mx1[indxg] = -mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}

void outflow_boundary_x2_inner(GridCons *grid, Parameters *params) {
    /* Simple outflow inner boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    il = -NGHX1; iu = nx1 + NGHX1;
    jl = -NGHX2; ju = 0;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /* -1 -> 0
             * -2 -> 1 
             * -3 -> 2
             *------------
             *  xxx | o
             *      |  o
             *      |   o
             *---------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(i,0);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}
void outflow_boundary_x1_inner(GridCons *grid, Parameters *params) {
    /* Simple outflow inner boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    il = -NGHX1; iu = 0;
    jl = -NGHX2; ju = nx2 + NGHX2;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /* 
             *------------
             *  xxx | o
             *      |  o
             *      |   o
             *---------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(0,j);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}


void outflow_boundary_x2_outer(GridCons *grid, Parameters *params) {
    /* Simple outflow outer boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    il = -NGHX1; iu = nx1 + NGHX1;
    jl = nx2; ju = nx2+NGHX2;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /* 
             *------------
             *    o | xxx
             *   o  |  
             *  o   |   
             *---------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(i,nx2-1);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}
void outflow_boundary_x1_outer(GridCons *grid, Parameters *params) {
    /* Simple outflow inner boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    il = nx1; iu = nx1+NGHX1;
    jl = -NGHX2; ju = nx2 + NGHX2;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /* -1 -> 0
             * -2 -> 1 
             * -3 -> 2
             *------------
             *  xxx | o
             *      |  o
             *      |   o
             *---------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(nx1-1,j);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}

void periodic_boundary_x2_inner(GridCons *grid, Parameters *params) {
    /* Periodic inner boundary.
     * 
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    il = -NGHX1; iu = nx1 + NGHX1;
    jl = -NGHX2; ju = 0;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /* -1 -> nx-1
             * -2 -> nx-2 
             * -3 -> nx-3
             *------------------
             *    x |         o|
             *   x  |        o |
             *  x   | o o o o  |
             *------------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(i,nx2+j);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}
void periodic_boundary_x1_inner(GridCons *grid, Parameters *params) {
    /* Periodic inner boundary.
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    jl = -NGHX2; ju = nx2 + NGHX2;
    il = -NGHX1; iu = 0;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /*
             *------------------
             *  | o       |x
             *  |  o      | x
             *  |   o o o |  x
             *------------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(nx1+i,j);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}
void periodic_boundary_x2_outer(GridCons *grid, Parameters *params) {
    /* Periodic outer boundary.
     * 
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    il = -NGHX1; iu = nx1 + NGHX1;
    jl = nx2; ju = nx2+NGHX2;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            /*
             * nx -> 0
             * nx+1 -> 1
             * nx+2 -> 2
             *------------------
             *  | o       |x
             *  |  o      | x
             *  |   o o o |  x
             *------------------
             */
            indxg = INDEX(i,j);
            indx = INDEX(i,j-nx2);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}
void periodic_boundary_x1_outer(GridCons *grid, Parameters *params) {
    /* Periodic outer boundary.
     */
    int i,j,indx,indxg;
    int jl,ju,il,iu;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    int n,ntot,nf;

    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;


    real *rho    = &grid->cons[0*ntot]; 
    real *mx1    = &grid->cons[1*ntot]; 
    real *mx2    = &grid->cons[2*ntot]; 
    real *mx3    = &grid->cons[3*ntot]; 
    real *energy = &grid->cons[4*ntot]; 
    real *intenergy = grid->intenergy;


    jl = -NGHX2; ju = nx2 + NGHX2;
    il = nx1; iu = nx1+NGHX1;

    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            indxg = INDEX(i,j);
            indx = INDEX(i-nx1,j);

            rho[indxg] = rho[indx];
            mx1[indxg] = mx1[indx];
            mx2[indxg] = mx2[indx];
            mx3[indxg] = mx3[indx];
            energy[indxg] = energy[indxg];
            intenergy[indxg] = intenergy[indx];

        }
    }



    return;
}
void ic_boundary_x2_inner(GridCons *grid, Parameters *params) {
    /* Initial condition inner boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
void ic_boundary_x1_inner(GridCons *grid, Parameters *params) {
    /* Initial condition inner boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
void ic_boundary_x2_outer(GridCons *grid, Parameters *params) {
    /* Initial condition outer boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
void ic_boundary_x1_outer(GridCons *grid, Parameters *params) {
    /* Initial condition outer boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
