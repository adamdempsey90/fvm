#include "defs.h"
#include "cuda_defs.h"

__device__ void reflecting_boundary_inner(int dir1,int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
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
    int indx;
    int dir2, dir3;
    /* 1->2->3
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;


    switch (dir1) {
    	case 1:
    		indx = GINDEX(-i-1,j,k);
    		break;
    	case 2:
    		indx = GINDEX(i,-j-1,k);
    		break;
    	case 3:
    		indx = GINDEX(i,j,-k-1);
    		break;
    	default:
    		indx = GINDEX(-i-1,j,k);
    		break;
    }

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + dir1*ntot] = -cons[indx + dir1*ntot];
    cons[indxg + dir2*ntot] =  cons[indx + dir2*ntot];
    cons[indxg + dir3*ntot] =  cons[indx + dir3*ntot];
    intenergy[indxg] = intenergy[indx];

    for(n=4;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}

__device__ void reflecting_boundary_outer(int dir1, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
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
    int indx;
    int dir2, dir3;
    /* 1->2->3
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;
    switch (dir1) {
    	case 1:
    		indx = GINDEX(nx1 + -(i-nx1)-1,j,k);
    		break;
    	case 2:
    		indx = GINDEX(i,nx2 + -(j-nx2)-1,k);
    		break;
    	case 3:
    		indx = GINDEX(i,j,nx3 + -(k-nx3)-1);
    		break;
    	default:
    		indx = GINDEX(nx1 + -(i-nx1)-1,j,k);
    		break;
    }

    cons[indxg + 0*ntot] = cons[indx + 0*ntot];
    cons[indxg + dir1*ntot] = -cons[indx + dir1*ntot];
    cons[indxg + dir2*ntot] =  cons[indx + dir2*ntot];
    cons[indxg + dir3*ntot] =  cons[indx + dir3*ntot];
    intenergy[indxg] = intenergy[indx];

    for(n=4;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }



    return;
}

__device__ void outflow_boundary_inner(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    /* Simple outflow inner boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
    /*
     *------------
     *    x | o o o
     *   x  |
     *  x   |
     *---------------
     */
    int n;
    int indx;
    switch (dir) {
    	case 1:
    		indx = GINDEX(0,j,k);
    		break;
    	case 2:
    		indx = GINDEX(i,0,k);
    		break;
    	case 3:
    		indx = GINDEX(i,j,0);
    		break;
    	default:
    		indx = GINDEX(0,j,k);
    		break;
    }

    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}

__device__ void outflow_boundary_outer(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    /* Simple outflow outer boundary.
     * This copies data from last cell into all ghost cells.
     * Note that this is NOT a non-reflecting b.c (especially for non-cartesian domains
     * or strong shearing flows).
     */
	/*
     *---------------
     * o o o | x
     *       |  x
     *       |   x
     *---------------
	 */
    int n;
    int indx;
    switch (dir) {
    	case 1:
    		indx = GINDEX(nx1-1,j,k);
    		break;
    	case 2:
    		indx = GINDEX(i,nx2-1,k);
    		break;
    	case 3:
    		indx = GINDEX(i,j,nx3-1);
    		break;
    	default:
    		indx = GINDEX(nx1-1,j,k);
    		break;
    }
    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }



    return;
}
__device__ void periodic_boundary_inner(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
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
    int indx;
    switch (dir) {
    	case 1:
    		indx = GINDEX(nx1+i,j,k);
    		break;
    	case 2:
    		indx = GINDEX(i,nx2+j,k);
    		break;
    	case 3:
    		indx = GINDEX(i,j,nx3+k);
    		break;
    	default:
    		indx = GINDEX(nx1+i,j,k);
    		break;
    }
    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }




    return;
}
__device__ void periodic_boundary_outer(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
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
    int indx;
    switch (dir) {
    	case 1:
    		indx = GINDEX(i-nx1,j,k);
    		break;
    	case 2:
    		indx = GINDEX(i,j-nx2,k);
    		break;
    	case 3:
    		indx = GINDEX(i,j,k-nx3);
    		break;
    	default:
    		indx = GINDEX(i-nx1,j,k);
    		break;
    }
    intenergy[indxg] = intenergy[indx];

    for(n=0;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }




    return;
}

__device__ void ic_boundary(int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    /* Initial condition inner boundary.
     * The update steps do not alter the ghost zone values, 
     * so they should remain at their initial values.
     */
    return;
}
