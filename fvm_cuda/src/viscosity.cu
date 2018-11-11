#include "defs.h"
#include "cuda_defs.h"

#ifdef VISCOSITY

extern __device__ real kinematic_viscosity(real x1, real x2, real x3);

__global__ void compute_divergence(real *cons, real *vel, 
        real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {
    int i,j,k,indx;

    real vm,vc,vp;

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=-NGHX1+1)&&(i<nx1+2)&&(j>=-NGHX2+1)&&(j<nx2+2)&&(k>=-NGHX3+1)&&(k<nx3+2)) {

            /* X1 direction */
            vm = cons[indx - 1 + 1*ntot]/cons[indx - 1];
            vc = cons[indx     + 1*ntot]/cons[indx];
            vp = cons[indx + 1 + 1*ntot]/cons[indx + 1];
            vel[indx + 1*ntot] = vc; 
            vel[indx]  = (vp-vm)/(2*dx1[i]);

            /* X2 direction */
#ifdef DIMS2
            vm = cons[indx - size_x1 + 2*ntot]/cons[indx - size_x1];
            vc = cons[indx           + 2*ntot]/cons[indx];
            vp = cons[indx + size_x1 + 2*ntot]/cons[indx + size_x1];
            vel[indx + 2*ntot] = vc; 
            vel[indx]  += (vp-vm)/(2*dx2[j]);
#endif


            /* X3 direction */
#ifdef DIMS3
            vm = cons[indx - size_x12 + 3*ntot]/cons[indx - size_x12];
            vc = cons[indx            + 3*ntot]/cons[indx];
            vp = cons[indx + size_x12 + 3*ntot]/cons[indx + size_x12];
            vel[indx + 3*ntot] = vc; 
            vel[indx]  += (vp-vm)/(2*dx3[k]);
#endif
        }



    }
    return;
}

__global__ void viscous_flux(real *vel, real *rho, real *F_1, real *F_2, real *F_3,
        real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {
	/*
	 *  *  *  *  *  *  *  *  *  *  *  *  *
	 *                                   *
	 *                 *                 *
	 *     (i,j+1)           (i+1,j+1)   *
	 *                 *                 *
 	 *                                   *
	 *  *  F2[i,j]  *  *  *  *  *  *  *  *
	 *                                   *
	 *                 *                 *
 F1[i-1,j] (i,j)    F1[i,j]  (i+1,j)     *
	 *                 *                 *
 	 *                                   *
	 *  *  F2[i,j-1]*  *  *  *  *  *  *  *
 	 *                                   *
	 *                 *                 *
	 *     (i,j-1)           (i+1,j-1)   *
	 *                 *                 *
 	 *                                   *
	 *  *  *  *  *  *  *  *  *  *  *  *  *
	 */

    int i,j,k;
    int indx;
    real s1,s2,s3,visc;
    real *divv = &vel[0];
    real *vx1 = &vel[ntot];
    real *vx2 = &vel[2*ntot];
    real *vx3 = &vel[3*ntot];

    real *F11 = &F_1[1*ntot];
    real *F12 = &F_1[2*ntot];
    real *F13 = &F_1[3*ntot];
    real *F1e = &F_1[4*ntot];
#ifdef DIMS2
    real *F21 = &F_2[1*ntot];
	real *F22 = &F_2[2*ntot];
	real *F23 = &F_2[3*ntot];
    real *F2e = &F_2[4*ntot];

#endif
#ifdef DIMS3
	real *F31 = &F_3[1*ntot];
	real *F32 = &F_3[2*ntot];
	real *F33 = &F_3[3*ntot];
    real *F3e = &F_3[4*ntot];

#endif
    s1=0; s2=0; s3=0;
    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);

        /* X1 direction */
        if ((i>=-NGHX1)&&(i<nx1+2)&&(j>=-NGHX2+1)&&(j<nx2+2)&&(k>=-NGHX3+1)&&(k<nx3+2)) {
        	visc=  .5*(rho[indx] + rho[indx+1]) * kinematic_viscosity(x1[i] +.5*dx1[i], x2[j], x3[k]);

            s1 = 2*(vx1[indx+1] - vx1[indx])/dx1[i] + (divv[indx] + divv[indx+1])/3.;
            s2 = (vx2[indx+1]-vx2[indx])/dx1[i];
#ifdef DIMS2
            s2 += .25*( (vx1[indx+1 + size_x1] - vx1[indx+1 - size_x1])/dx2[j+1] + (vx1[indx+size_x1]-vx1[indx-size_x1])/dx2[j]);
#endif
            s3 = (vx3[indx+1]-vx3[indx])/dx1[i];
#ifdef DIMS3
            s3 += .25*( (vx1[indx+1 + size_x12] - vx1[indx+1 - size_x12])/dx3[k+1] + (vx1[indx+size_x12]-vx1[indx-size_x12])/dx3[k]);
#endif
            F11[indx] -= visc*s1;
            F12[indx] -= visc*s2;
            F13[indx] -= visc*s3;
            F1e[indx] -= .5*visc*( (vx1[indx + 1] + vx1[indx])*s1
            		 	 	 	+(vx2[indx + 1] + vx2[indx])*s2
            		 	 	 	+(vx3[indx + 1] + vx3[indx])*s3);
        }
        /* X2 direction */
#ifdef DIMS2
		if ((i>=-NGHX1+1)&&(i<nx1+2)&&(j>=-NGHX2)&&(j<nx2+2)&&(k>=-NGHX3+1)&&(k<nx3+2)) {
			visc=  .5*(rho[indx] + rho[indx+size_x1])* kinematic_viscosity(x1[i], x2[j] + .5*dx2[j], x3[k]);
			s1 = .25*( (vx2[indx+size_x1 + 1] - vx2[indx+size_x1 - 1])/dx1[i+1] + (vx2[indx+1]-vx2[indx-1])/dx1[i])
					  + (vx1[indx+size_x1]-vx1[indx])/dx2[j];
			s2 = 2*(vx2[indx+size_x1] - vx2[indx])/dx2[j] + (divv[indx] +divv[indx+size_x1])/3;
			s3 = (vx3[indx+size_x1]-vx3[indx])/dx2[j];
#ifdef DIMS3
			s3 += .25*( (vx2[indx+size_x1 + size_x12] - vx1[indx+size_x1 - size_x12])/dx3[k+1]
			                   + (vx2[indx+size_x12]-vx1[indx-size_x12])/dx3[k]);
#endif
            F21[indx] -= visc*s1;
            F22[indx] -= visc*s2;
            F23[indx] -= visc*s3;
            F2e[indx] -= .5*visc*( (vx1[indx + size_x1] + vx1[indx])*s1
            		 	 	 	+(vx2[indx + size_x1] + vx2[indx])*s2
            		 	 	 	+(vx3[indx + size_x1] + vx3[indx])*s3);
		}
#endif
        /* X3 direction */
#ifdef DIMS3
		if ((i>=-NGHX1+1)&&(i<nx1+2)&&(j>=-NGHX2+1)&&(j<nx2+2)&&(k>=-NGHX3)&&(k<nx3+2)) {
			visc=  .5*(rho[indx] + rho[indx+size_x12]) * kinematic_viscosity(x1[i], x2[j], x3[k] + .5*dx3[k]);
			s1 = .25*( (vx3[indx+size_x12 + 1] - vx3[indx+size_x12 - 1])/dx1[i+1] + (vx3[indx+1]-vx3[indx-1])/dx1[i])
					  + (vx1[indx+size_x12]-vx1[indx])/dx3[k];
			s2 = .25*( (vx3[indx+size_x12 + size_x1] - vx3[indx+size_x12 - size_x1])/dx2[j+1]
			                   + (vx3[indx+size_x1]-vx3[indx-size_x1])/dx2[j])
					  + (vx2[indx+size_x12]-vx2[indx])/dx3[k];
			s3 = 2*(vx3[indx+size_x12] - vx3[indx])/dx3[k] + (divv[indx] +divv[indx+size_x12])/3;
            F31[indx] -= visc*s1;
            F32[indx] -= visc*s2;
            F33[indx] -= visc*s3;
            F3e[indx] -= .5*visc*( (vx1[indx + size_x12] + vx1[indx])*s1
            		 	 	 	+(vx2[indx + size_x12] + vx2[indx])*s2
            		 	 	 	+(vx3[indx + size_x12] + vx3[indx])*s3);
		}
#endif

    }
    return;

}

#endif
