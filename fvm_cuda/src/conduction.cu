#include "defs.h"
#include "cuda_defs.h"


#ifdef CONDUCTION

extern __host__ __device__ real heatcond_func(real dens, real x1, real x2, real x3, real delad);

__global__ void conduction_flux(real *cons, real *intenergy, real *F_1, real *F_2, real *F_3,
        real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3, real g,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {
    int i,j,k;
    int indx,indxp;

    real cond, tempc, tempr;
    real delad = 1. -1./g;
    
    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        
        /* X1 direction */ 
        if ((i>=-NGHX1)&&(i<nx1+2)&&(j>=-NGHX2)&&(j<nx2+NGHX2)&&(k>=-NGHX3)&&(k<nx3+NGHX3)) {
            indxp = indx + 1;
            tempc = intenergy[indx] * g / cons[indx]; // Cp = 1 -> rho*e = Cp*T/gamma
            tempr = intenergy[indxp] * g / cons[indxp];
            cond = heatcond_func(.5*(cons[indx]+cons[indxp]),x1[i]+.5*dx1[i],x2[j],x3[k],delad);
                
            F_1[indx + 4*ntot] -= cond * (tempr - tempc) / dx1[i];

        }
        /* X2 direction */
#ifdef DIMS2
		if ((i>=-NGHX1)&&(i<nx1+NGHX1)&&(j>=-NGHX2)&&(j<nx2+2)&&(k>=-NGHX3)&&(k<nx3+NGHX3)) {
			indxp = indx + size_x1;
			tempc = intenergy[indx] * g / cons[indx];
			tempr = intenergy[indxp] * g / cons[indxp];
			cond = heatcond_func(.5*(cons[indx]+cons[indxp]),x1[i],x2[j]+.5*dx2[j],x3[k],delad);

			F_2[indx + 4*ntot] -= cond * (tempr - tempc) / dx2[j];

		}
#endif
        /* X3 direction */
#ifdef DIMS3
		if ((i>=-NGHX1)&&(i<nx1+NGHX1)&&(j>=-NGHX2)&&(j<nx2+NGHX2)&&(k>=-NGHX3)&&(k<nx3+2)) {
			indxp = indx + size_x12;
			tempc = intenergy[indx] * g / cons[indx];
			tempr = intenergy[indxp] * g / cons[indxp];
			cond = heatcond_func(.5*(cons[indx]+cons[indxp]),x1[i],x2[j],x3[k]+.5*dx3[k],delad);

			F_3[indx + 4*ntot] -= cond * (tempr - tempc) / dx3[k];

		}
#endif
            
    }
    return;

}

//__global__ void conduction_update(real *cons, real *intenergy, real *F_1, real *F_2, real *F_3,
//        real *dx1, real *dx2, real *dx3, real dt,
//        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {
//    int i,j,k;
//    int indx;
//    real rhs;
//
//
//    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
//    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
//        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)&&(k>=0)&&(k<nx3)) {
//
//            rhs =  dt/dx1[i]*(F_1[indx +4*ntot] -F_1[indx - 1       +4*ntot]);
//#ifdef DIMS2
//            rhs += dt/dx2[j]*(F_2[indx +4*ntot] -F_2[indx - size_x1 +4*ntot]);
//#endif
//#ifdef DIMS3
//            rhs += dt/dx3[k]*(F_3[indx +4*ntot] -F_3[indx - size_x12 +4*ntot]);
//#endif
//
//            cons[indx + 4*ntot] += rhs;
//            intenergy[indx] += rhs;
//            printf("(%d,%d)\t%e\t%e\n",i,j,intenergy[indx],cons[indx+4*ntot]-intenergy[indx]);
//
//        }
//    }
//    return;
//
//}
//__global__ void diffusion_update(real *cons, real *intenergy, real *F_1, real *F_2, real *F_3,
//        real *dx1, real *dx2, real *dx3, real dt, int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf) {
//    int i,j,k,n;
//    int indx;
//    int nan_check;
//    real dtdx1;
//#ifdef DIMS2
//    real dtdx2;
//#endif
//#ifdef DIMS3
//    real dtdx3;
//#endif
//
//    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
//    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
//        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)&&(k>=0)&&(k<nx3)) {
//            dtdx1 = dt/dx1[i];
//#ifdef DIMS2
//            dtdx2 = dt/dx2[j];
//#endif
//#ifdef DIMS3
//            dtdx3 = dt/dx3[k];
//#endif
//            for(n=0;n<nf;n++) {
////            	printf("%d (%d,%d), %e %e %e\n",n,i,j,
////            			F_2[indx + n*ntot],F_2[indx - size_x1        + n*ntot],
////            			F_2[indx - size_x1        + n*ntot]-F_2[indx + n*ntot]);
//
//                cons[indx + n*ntot] += dtdx1*(F_1[indx - 1        + n*ntot]- F_1[indx + n*ntot]);
//#ifdef DIMS2
//                cons[indx + n*ntot] += dtdx2*(F_2[indx - size_x1  + n*ntot]- F_2[indx + n*ntot]);
//#endif
//#ifdef DIMS3
//                cons[indx + n*ntot] += dtdx3*(F_3[indx - size_x12 + n*ntot]- F_3[indx + n*ntot]);
//#endif
//
//            }
//
//            intenergy[indx] = cons[indx+4*ntot] - .5*(
//                    cons[indx + 1*ntot]*cons[indx + 1*ntot] +
//                    cons[indx + 2*ntot]*cons[indx + 2*ntot] +
//                    cons[indx + 3*ntot]*cons[indx + 3*ntot])/cons[indx];
////            printf("%e\t%e\t%e\t%e\n",cons[indx+4*ntot] ,  .5*(
////                    cons[indx + 1*ntot]*cons[indx + 1*ntot] +
////                    cons[indx + 2*ntot]*cons[indx + 2*ntot] +
////                    cons[indx + 3*ntot]*cons[indx + 3*ntot])/cons[indx],
////                    intenergy[indx], cons[indx + 4*ntot]-intenergy[indx]);
//
//        }
//    }
//    return;
//
//}
#endif
