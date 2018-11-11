#include "defs.h"
#include "cuda_defs.h"


#ifdef CONDUCTION
extern __device__ real thermal_diff(real rho, real x1, real x2, real x3, real delad);
#endif
#ifdef VISCOSITY
extern __device__ real kinematic_viscosity(real x1, real x2, real x3);
#endif

__inline__ __device__ real warpReduceMin(real val) {
    for(int offset = 16; offset >0; offset /= 2) {
#if __CUDACC_VER_MAJOR__ < 9
        real tmp_val = __shfl_down(val,offset);
        if (tmp_val < val)
            val = tmp_val;
#else
        real tmp_val = __shfl_down_sync(FULL_MASK,val,offset);
        if (tmp_val < val)
            val = tmp_val;
#endif


    }
    return val;
}
__inline__ __device__ real blockReduceMin(real val) {
    static __shared__ real shared[32];
    int lane = threadIdx.x % 32;
    int wid = threadIdx.x / 32;

    val = warpReduceMin(val);

    if (lane == 0) shared[wid] = val;
    __syncthreads();

    val = (threadIdx.x < blockDim.x / 32) ? shared[lane] : FLOATMAX;
    if (wid ==0) {

    	val = warpReduceMin(val);
    }
    return val;
}

__global__ void timestep_kernel(real *cons, real *dx1, real *dx2,real *dx3, real *x1, real *x2,real *x3, real *out ,
		int nx1, int nx2, int nx3, int size_x1, int size_x12,int ntot,int offset, real g) {
    int i,j,k,indx;
    real curr_min = FLOATMAX ;
    real pres,cs,dt1,dt2,dt3,dt;
    dt1 = FLOATMAX;
    dt2 = FLOATMAX;
    dt3 = FLOATMAX;
#ifdef CONDUCTION
    real chi;
    real delad = 1. - 1./g;
#endif
#ifdef VISCOSITY
    real nu;
#endif

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot;indx +=blockDim.x*gridDim.x) {
    	unpack_indices(indx,&i,&j,&k,size_x1,size_x12);
        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)&&(k>=0)&&(k<nx3)) {

            pres = (g-1)*(cons[indx + 4*ntot]-
                    .5*(  cons[indx + 1*ntot] * cons[indx + 1*ntot]
                        + cons[indx + 2*ntot] * cons[indx + 2*ntot]
                        + cons[indx + 3*ntot] * cons[indx + 3*ntot] )/cons[indx]);

            cs = sqrt( g* pres/cons[indx]);
            dt1 = dx1[i]/(fabs(cons[indx + 1*ntot]/cons[indx]) + cs);
#ifdef DIMS2
            dt2 = dx2[j]/(fabs(cons[indx + 2*ntot]/cons[indx]) + cs);
#endif
#ifdef DIMS3
            dt3 = dx3[j]/(fabs(cons[indx + 3*ntot]/cons[indx]) + cs);
#endif

            dt = MIN3(dt1,dt2,dt3);
            

#ifdef CONDUCTION
            chi = thermal_diff(cons[indx],x1[i],x2[j],x3[k], delad);
            dt1 = dx1[i]*dx1[i]/chi;
            dt = MIN2(dt,dt1);
#ifdef DIMS2
            dt2 = dx2[j]*dx2[j]/chi;
            dt = MIN2(dt,dt2);
#endif
#ifdef DIMS3
            dt3 = dx3[k]*dx3[k]/chi;
            dt = MIN2(dt,dt3);
#endif
#endif


#ifdef VISCOSITY
            nu = kinematic_viscosity(x1[i],x2[j],x3[k]);
            dt1 = dx1[i]*dx1[i]/nu;
            dt = MIN2(dt,dt1);
#ifdef DIMS2
            dt2 = dx2[j]*dx2[j]/nu;
            dt = MIN2(dt,dt2);
#endif
#ifdef DIMS3
            dt3 = dx3[k]*dx3[k]/nu;
            dt = MIN2(dt,dt3);
#endif
#endif

            if (dt < curr_min) curr_min = dt;


        }


    }
    curr_min = blockReduceMin(curr_min);
    if (threadIdx.x ==0) out[blockIdx.x]=curr_min;
    return;
}
__global__ void timestep_kernel_final(real *in, real *out ,int n, real cfl) {
    real val = FLOATMAX;
    for(int i = blockIdx.x*blockDim.x + threadIdx.x; i<n;i+=blockDim.x*gridDim.x) {
        if (in[i] < val)
            val = in[i];
    }
    val = blockReduceMin(val);
    if (threadIdx.x ==0) out[blockIdx.x]=val;
    return;
}
