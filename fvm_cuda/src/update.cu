#include "defs.h"
#include "cuda_defs.h"

__global__ void compute_dhalf(real *cons, real *dhalf, real *F_1, real *F_2,real *F_3,
        real *dx1, real *dx2, real *dx3, real dt, int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {
    int i,j,k;
    int indx;


    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=-2)&&(i<nx1+2)&&(j>=-2)&&(j<nx2+2)&&(k>=-2)&&(k<nx3+2)) {
            dhalf[indx] = cons[indx] +  .5*dt/dx1[i]*(F_1[indx - 1]        - F_1[indx]);
#ifdef DIMS2
            dhalf[indx] += .5*dt/dx2[j]*(F_2[indx - size_x1]  - F_2[indx]);
#endif
#ifdef DIMS3
            dhalf[indx] += .5*dt/dx3[k]*(F_3[indx - size_x12] - F_3[indx]);
#endif
        }
    }
    return;

}

__global__ void update_cons(real *cons, real *intenergy, real *F_1, real *F_2, real *F_3,
        real *dx1, real *dx2, real *dx3, real dt, int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {
    int i,j,k,n;
    int indx;
    int nan_check;
    real dtdx1;
#ifdef DIMS2
    real dtdx2;
#endif
#ifdef DIMS3
    real dtdx3;
#endif

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)&&(k>=0)&&(k<nx3)) {
            dtdx1 = dt/dx1[i];
#ifdef DIMS2
            dtdx2 = dt/dx2[j];
#endif
#ifdef DIMS3
            dtdx3 = dt/dx3[k];
#endif
            for(n=0;n<nf;n++) {
            	//printf("%d (%d,%d), %lg \n",n,i,j,F_1[indx + n*ntot]);
                cons[indx + n*ntot] += dtdx1*(F_1[indx - 1        + n*ntot]- F_1[indx + n*ntot]);
#ifdef DIMS2
                cons[indx + n*ntot] += dtdx2*(F_2[indx - size_x1  + n*ntot]- F_2[indx + n*ntot]);
#endif
#ifdef DIMS3
                cons[indx + n*ntot] += dtdx3*(F_3[indx - size_x12 + n*ntot]- F_3[indx + n*ntot]);
#endif

            }
            intenergy[indx] = cons[indx+4*ntot] - .5*(
                    cons[indx + 1*ntot]*cons[indx + 1*ntot] +
                    cons[indx + 2*ntot]*cons[indx + 2*ntot] +
                    cons[indx + 3*ntot]*cons[indx + 3*ntot])/cons[indx];

        }
    }
    return;

}
__global__ void transverse_update(real *UL_1, real *UL_2, real *UL_3,
        real *UR_1, real *UR_2, real *UR_3,
        real *F_1, real *F_2, real *F_3, real *dx1, real *dx2, real *dx3, real dt,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {

	/*
	 * 			G(i,j)			G(i+1,j)
	 *		+	+	+	+	+	+	+	+	+
	 *		+				+				+
	 *		+		UL(i,j)	+ UR(i,j)		+
	 *		+	  			+	 			+
	 *		+	  (i,j)		+ 	  (i+1,j)	+
	 *		+				+ 				+
	 * 		+	+	+	+	+	+	+	+	+
	 *			G(i,j-1)		G(i+1,j-1)
	 *
	 *
	 */
    int i,j,k,n;
    int indx;
    real dtdx;

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot;indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        /* X1 - direction */
        if ((i>=-2)&&(i<nx1+2)&&(j>=-2)&&(j<nx2+3)&&(k>=-2)&&(k<nx3+3)) {
            dtdx = .5*dt/dx2[j];
            for(n=0;n<nf;n++) {
                UL_1[indx + n*ntot] += dtdx*(F_2[indx - size_x1     + n*ntot] -F_2[indx     + n*ntot]);
                UR_1[indx + n*ntot] += dtdx*(F_2[indx - size_x1 + 1 + n*ntot] -F_2[indx + 1 + n*ntot]);

            }
            /* Add X3 flux */
#ifdef DIMS3
			dtdx = .5*dt/dx3[k];
			for(n=0;n<nf;n++) {
				UL_1[indx + n*ntot] += dtdx*(F_3[indx - size_x12     + n*ntot] -F_3[indx     + n*ntot]);
				UR_1[indx + n*ntot] += dtdx*(F_3[indx - size_x12 + 1 + n*ntot] -F_3[indx + 1 + n*ntot]);

			}

#endif
        }
        /* X2 - direction */
        if ((i>=-2)&&(i<nx1+3)&&(j>=-2)&&(j<nx2+2)&&(k>=-2)&&(k<nx3+3)) {
            /* Add X1 flux */
            dtdx = .5*dt/dx1[i];
            for(n=0;n<nf;n++) {
                UL_2[indx + n*ntot] += dtdx*(F_1[indx - 1           + n*ntot] -F_1[indx           + n*ntot]);
                UR_2[indx + n*ntot] += dtdx*(F_1[indx - 1 + size_x1 + n*ntot] -F_1[indx + size_x1 + n*ntot]);

            }
            /* Add X3 flux */
#ifdef DIMS3
                /* Add X1 flux */
			dtdx = .5*dt/dx3[k];
			for(n=0;n<nf;n++) {
				UL_2[indx + n*ntot] += dtdx*(F_3[indx - size_x12           + n*ntot] -F_3[indx           + n*ntot]);
				UR_2[indx + n*ntot] += dtdx*(F_3[indx - size_x12 + size_x1 + n*ntot] -F_3[indx + size_x1 + n*ntot]);

			}

#endif
        }
        /* X3 - direction */
#ifdef DIMS3
        if ((i>=-2)&&(i<nx1+3)&&(j>=-2)&&(j<nx2+3)&&(k>=-2)&&(k<nx3+2)) {
            /* Add X1 flux */
            dtdx = .5*dt/dx1[i];
            for(n=0;n<nf;n++) {
                UL_3[indx + n*ntot] += dtdx*(F_1[indx - 1            + n*ntot] -F_1[indx            + n*ntot]);
                UR_3[indx + n*ntot] += dtdx*(F_1[indx - 1 + size_x12 + n*ntot] -F_1[indx + size_x12 + n*ntot]);

            }
            /* Add X2 flux */
			dtdx = .5*dt/dx2[j];
			for(n=0;n<nf;n++) {
				UL_3[indx + n*ntot] += dtdx*(F_2[indx - size_x1            + n*ntot] -F_2[indx            + n*ntot]);
				UR_3[indx + n*ntot] += dtdx*(F_2[indx - size_x1 + size_x12 + n*ntot] -F_2[indx + size_x12 + n*ntot]);

			}


        }
#endif
    }
    return;

}
