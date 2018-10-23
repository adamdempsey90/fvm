#include "defs.h"
#include "cuda_defs.h"


__device__ void reflecting_boundary_x2_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Reflecting inner boundary.
     * This sets  dU/dx = 0 and u2 = 0 
     */
    /* -1 -> 0
     * -2 -> 1 
     * -3 -> 2
     *------------
     *    x | o
     *   x  |  o
     *  x   |   o
     *---------------
     */
    int n;
    int indx = GINDEX(i,-j-1);

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + 1*ntot] = cons[indx + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx + 3*ntot];
    intenergy[indxg] = intenergy[indx];

    for(n=4;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}
__device__ void reflecting_boundary_x1_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Reflecting inner boundary.
     * This sets  dU/dx = 0 and u2 = 0 
     */
    int n;
    int indx = GINDEX((-i-1),j);

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + 1*ntot] = -cons[indx + 1*ntot];
    cons[indxg + 2*ntot] = cons[indx + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx + 3*ntot];
    intenergy[indxg] = intenergy[indx];

    for(n=4;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}

__device__ void reflecting_boundary_x2_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Reflecting outer boundary.
     * This sets  dU/dx = 0 and u2 = 0 
     */
    /* n -> n-1
     * n+1 -> n-2 
     * n+2-> n-3
     *------------
     *    o | x
     *   o  |  x
     *  o   |   x
     *---------------
     */
    int n;
    int indx = GINDEX(i,nx2 + -(j-nx2)-1);

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + 1*ntot] = cons[indx + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx + 3*ntot];
    intenergy[indxg] = intenergy[indx];

    for(n=4;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }



    return;
}
__device__ void reflecting_boundary_x1_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Reflecting inner boundary.
     * This sets  dU/dx = 0 and u2 = 0 
     */
    int n;
    int indx = GINDEX(nx1 + -(i-nx1)-1,j); 

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + 1*ntot] = -cons[indx + 1*ntot];
    cons[indxg + 2*ntot] = cons[indx + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx + 3*ntot];
    intenergy[indxg] = intenergy[indx];

    for(n=4;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}

__device__ void outflow_boundary_x2_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Simple outflow inner boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */

    int n;
    int indx = GINDEX(i,0); 

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}
__device__ void outflow_boundary_x1_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Simple outflow inner boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
    int n;
    int indx = GINDEX(0,j); 

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }



    return;
}


__device__ void outflow_boundary_x2_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Simple outflow outer boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
    int n;
    int indx = GINDEX(i,nx2-1); 

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }



    return;
}
__device__ void outflow_boundary_x1_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Simple outflow inner boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
    int n;
    int indx = GINDEX(nx1-1,j); 

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}

__device__ void periodic_boundary_x2_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Periodic inner boundary.
     * 
     * -1 -> nx-1
     * -2 -> nx-2 
     * -3 -> nx-3
     *------------------
     *    x |         o|
     *   x  |        o |
     *  x   | o o o o  |
     *------------------
     */
    int n;
    int indx = GINDEX(i, nx2+j); 

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }




    return;
}
__device__ void periodic_boundary_x1_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Periodic inner boundary.
     */
    int n;
    int indx = GINDEX(nx1+i,j); 

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}
__device__ void periodic_boundary_x2_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Periodic outer boundary.
     * 
     * nx -> 0
     * nx+1 -> 1
     * nx+2 -> 2
     *------------------
     *  | o       |x
     *  |  o      | x
     *  |   o o o |  x
     *------------------
     */
    int n;
    int indx = GINDEX(i,j-nx2); 

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }




    return;
}
__device__ void periodic_boundary_x1_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Periodic outer boundary.
     */
    int n;
    int indx = GINDEX(i-nx1,j); 

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }




    return;
}
__device__ void ic_boundary_x2_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Initial condition inner boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
__device__ void ic_boundary_x1_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Initial condition inner boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
__device__ void ic_boundary_x2_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Initial condition outer boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
__device__ void ic_boundary_x1_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) {
    /* Initial condition outer boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
