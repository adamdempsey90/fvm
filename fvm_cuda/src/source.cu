#include "defs.h"
#include "cuda_defs.h"
#ifdef POTENTIAL



extern __device__ real gravpot(real x, real y, real z);


__global__ void update_source(real *cons, real *dhalf, real *F1, real *F2, real *F3,
		real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int nf,int ntot, int offset, real dt) {


    int i,j,k,indx;
    real phil,phir,phic;


    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)&&(k>=0)&&(k<nx3)) {


            /* x1-direction */
            phic = gravpot(x1[i]            ,x2[j],x3[k]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j],x3[k]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j],x3[k]);
            cons[indx + 1*ntot] -= dt/dx1[i] * (phir - phil) * dhalf[indx];

            /* energy */
            cons[indx + 4*ntot] -= dt/dx1[i]*(F1[indx - 1] *(phic-phil)
                    + F1[indx]*(phir-phic));

#ifdef DIMS2
            /* x2-direction */
            phil = gravpot(x1[i], x2[j] - .5*dx2[j], x3[k]);
            phir = gravpot(x1[i], x2[j] + .5*dx2[j], x3[k]);
            cons[indx + 2*ntot] -= dt/dx2[j] * (phir - phil) * dhalf[indx];

            /* energy */
            cons[indx + 4*ntot] -= dt/dx2[j]*(F2[indx - size_x1] *(phic-phil)
                    + F2[indx]*(phir-phic));
#endif
#ifdef DIMS3
            /* x3-direction */
            phil = gravpot(x1[i], x2[j], x3[k] - .5*dx3[k]);
            phir = gravpot(x1[i], x2[j], x3[k] + .5*dx3[k]);
            cons[indx + 3*ntot] -= dt/dx3[k] * (phir - phil) * dhalf[indx];

            /* energy */
            cons[indx + 4*ntot] -= dt/dx3[k]*(F3[indx - size_x12] *(phic-phil)
                    + F3[indx]*(phir-phic));
#endif

        }
    }
    return;
}

__global__ void source_terms(real *UL, real *UR, real *dx, real *x1, real *x2, real *x3,
        int dir1,int nx1, int nx2, int nx3, int size_x1, int size_x12,
        int nf,int ntot, int offset, real g1, real dt) {
	/*
	 *		+	+	+	+	+	+	+	+	+
	 *		+		UL(i,j)	+ UR(i,j)		+
	 *		+				+ 				+
	 *		+	  	   Phi(i+1/2,j)	 		+
	 *		+				+				+
	 *		+	  Phi(i,j)	+ 	Phi(i+1,j)  +
	 *		+				+ 				+
	 * 		+	+	+	+	+	+	+	+	+
	 *
	 */
    int i,j,k,indx;
    int il, iu, jl, ju, kl,ku;
    real dtdx;
    real phir,phil,phic;

    if (dir1 == 1) {
        il = -NGHX1; iu = nx1+2;
        jl = -NGHX2;ju = nx2+NGHX2;
        kl = -NGHX3; ku = nx3+NGHX3;
    }
    else if (dir1 == 2) {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -NGHX2; ju = nx2+2;
        kl = -NGHX3; ku = nx3+NGHX3;
    }
    else {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -NGHX2; ju = nx2+NGHX2;
        kl = -NGHX3; ku = nx3+2;
    }
    
    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=il)&&(i<iu)&&(j>=jl)&&(j<ju)&&(k>=kl)&&(k<ku)) {
            if (dir1 == 1) {
                dtdx = dt/dx[i];
                phil = gravpot(x1[i]           ,x2[j], x3[k]);
                phic = gravpot(x1[i] + .5*dx[i],x2[j], x3[k]);
                phir = gravpot(x1[i+1]         ,x2[j], x3[k]);
            }
            else if (dir1 == 2) {
                dtdx = dt/dx[j];
                phil = gravpot(x1[i], x2[j]           , x3[k]);
                phic = gravpot(x1[i], x2[j] + .5*dx[j], x3[k]);
                phir = gravpot(x1[i], x2[j+1]         , x3[k]);
            }
            else {
                dtdx = dt/dx[k];
                phil = gravpot(x1[i], x2[j], x3[k]           );
                phic = gravpot(x1[i], x2[j], x3[k] + .5*dx[k]);
                phir = gravpot(x1[i], x2[j], x3[k+1]         );

            }
            UL[indx + dir1*ntot] -= dtdx*(phic-phil)*UL[indx];
            UR[indx + dir1*ntot] -= dtdx*(phir-phic)*UR[indx];

        }

    }
    return;

}
__global__ void source_transverse_update(real *cons, real *UL_1, real *UL_2, real *UL_3, real *UR_1, real *UR_2, real *UR_3,
		real *F_1, real *F_2, real *F_3, real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3, real dt,
		int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {
    int i,j,k;
    int indx;
    real phil,phir,phic;
    real dtdx2,dtdx1;
#ifdef DIMS3
    real dtdx3;
#endif


    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot;indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        /* X1 - direction */
        if ((i>=-2)&&(i<nx1+2)&&(j>=-2)&&(j<nx2+1)&&(k>=-2)&&(k<nx3+1)) {

            dtdx2 = .5*dt/dx2[j];

            /* Add X2 contribution */
            phic = gravpot(x1[i],x2[j]            , x3[k]);
            phir = gravpot(x1[i],x2[j] + .5*dx2[j], x3[k]);
            phil = gravpot(x1[i],x2[j] - .5*dx2[j], x3[k]);
            
            UL_1[indx + 2*ntot] -= dtdx2*(phir-phil)*cons[indx];
            UL_1[indx + 4*ntot] -= dtdx2*((phir-phic)*F_2[indx] + (phic-phil)*F_2[indx - size_x1]);

            phic = gravpot(x1[i+1],x2[j]            , x3[k]);
            phir = gravpot(x1[i+1],x2[j] + .5*dx2[j], x3[k]);
            phil = gravpot(x1[i+1],x2[j] - .5*dx2[j], x3[k]);

            UR_1[indx + 2*ntot] -= dtdx2*(phir-phil)*cons[indx+1];
            UR_1[indx + 4*ntot] -= dtdx2*((phir-phic)*F_2[indx+1] + (phic-phil)*F_2[indx - size_x1 + 1]);

#ifdef DIMS3
            /* Add X3 contribution */
            	dtdx3 = .5*dt/dx3[k];
				phic = gravpot(x1[i],x2[j], x3[k]);
				phir = gravpot(x1[i],x2[j], x3[k] + .5*dx3[k]);
				phil = gravpot(x1[i],x2[j], x3[k] - .5*dx3[k]);

				UL_1[indx + 3*ntot] -= dtdx3*(phir-phil)*cons[indx];
				UL_1[indx + 4*ntot] -= dtdx3*((phir-phic)*F_3[indx] + (phic-phil)*F_3[indx - size_x12]);

				phic = gravpot(x1[i+1],x2[j]            , x3[k]);
				phir = gravpot(x1[i+1],x2[j], x3[k] + .5*dx3[k]);
				phil = gravpot(x1[i+1],x2[j], x3[k] - .5*dx3[k]);

				UR_1[indx + 3*ntot] -= dtdx3*(phir-phil)*cons[indx+1];
				UR_1[indx + 4*ntot] -= dtdx3*((phir-phic)*F_3[indx+1] + (phic-phil)*F_3[indx - size_x12 + 1]);

#endif

        }
        /* X2 - direction */
        if ((i>=-2)&&(i<nx1+1)&&(j>=-2)&&(j<nx2+1)&&(k>=-2)&&(k<nx3+1)) {
            dtdx1 = .5*dt/dx1[i];
            /* Add X1 contribution */
            phic = gravpot(x1[i]            ,x2[j], x3[k]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j], x3[k]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j], x3[k]);
            UL_2[indx + 1*ntot] -= dtdx1*(phir-phil)*cons[indx];
            UL_2[indx + 4*ntot] -= dtdx1*((phir-phic)*F_1[indx] + (phic-phil)*F_1[indx - 1]);

            phic = gravpot(x1[i]            ,x2[j+1], x3[k]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j+1], x3[k]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j+1], x3[k]);

            UR_2[indx + 1*ntot] -= dtdx1*(phir-phil)*cons[indx + size_x1];
            UR_2[indx + 4*ntot] -= dtdx1*((phir-phic)*F_1[indx + size_x1] + (phic-phil)*F_1[indx - 1 + size_x1]);

#ifdef DIMS3
            	/* Add X3 contribution */
            	dtdx3 = .5*dt/dx3[k];
				phic = gravpot(x1[i],x2[j], x3[k]            );
				phir = gravpot(x1[i],x2[j], x3[k] + .5*dx3[k]);
				phil = gravpot(x1[i],x2[j], x3[k] - .5*dx3[k]);

				UL_2[indx + 3*ntot] -= dtdx3*(phir-phil)*cons[indx];
				UL_2[indx + 4*ntot] -= dtdx3*((phir-phic)*F_3[indx] + (phic-phil)*F_3[indx - size_x12]);

				phic = gravpot(x1[i],x2[j+1], x3[k]            );
				phir = gravpot(x1[i],x2[j+1], x3[k] + .5*dx3[k]);
				phil = gravpot(x1[i],x2[j+1], x3[k] - .5*dx3[k]);

				UR_2[indx + 3*ntot] -= dtdx3*(phir-phil)*cons[indx+size_x1];
				UR_2[indx + 4*ntot] -= dtdx3*((phir-phic)*F_3[indx+size_x1] + (phic-phil)*F_3[indx - size_x12 + size_x1]);

#endif


        }
        /* X3 - direction */
#ifdef DIMS3
        if ((i>=-2)&&(i<nx1+1)&&(j>=-2)&&(j<nx2+1)&&(k>=-2)&&(k<nx3+1)) {
            dtdx1 = .5*dt/dx1[i];
            /* Add X1 contribution */
            phic = gravpot(x1[i]            ,x2[j], x3[k]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j], x3[k]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j], x3[k]);
            UL_3[indx + 1*ntot] -= dtdx1*(phir-phil)*cons[indx];
            UL_3[indx + 4*ntot] -= dtdx1*((phir-phic)*F_1[indx] + (phic-phil)*F_1[indx - 1]);

            phic = gravpot(x1[i]            ,x2[j], x3[k+1]);
            phir = gravpot(x1[i] + .5*dx1[i],x2[j], x3[k+1]);
            phil = gravpot(x1[i] - .5*dx1[i],x2[j], x3[k+1]);

            UR_3[indx + 1*ntot] -= dtdx1*(phir-phil)*cons[indx + size_x12];
            UR_3[indx + 4*ntot] -= dtdx1*((phir-phic)*F_1[indx + size_x12] + (phic-phil)*F_1[indx - 1 + size_x12]);


			/* Add X2 contribution */
            dtdx2 = .5*dt/dx2[j];
            phic = gravpot(x1[i],x2[j]            , x3[k]);
            phir = gravpot(x1[i],x2[j] + .5*dx2[j], x3[k]);
            phil = gravpot(x1[i],x2[j] - .5*dx2[j], x3[k]);

            UL_3[indx + 2*ntot] -= dtdx2*(phir-phil)*cons[indx];
            UL_3[indx + 4*ntot] -= dtdx2*((phir-phic)*F_2[indx] + (phic-phil)*F_2[indx - size_x1]);

            phic = gravpot(x1[i],x2[j]            , x3[k+1]);
            phir = gravpot(x1[i],x2[j] + .5*dx2[j], x3[k+1]);
            phil = gravpot(x1[i],x2[j] - .5*dx2[j], x3[k+1]);

            UR_3[indx + 2*ntot] -= dtdx2*(phir-phil)*cons[indx+size_x12];
            UR_3[indx + 4*ntot] -= dtdx2*((phir-phic)*F_2[indx+size_x12] + (phic-phil)*F_2[indx - size_x1 + size_x12]);


            


        }
#endif
    }
    return;

}
#endif
