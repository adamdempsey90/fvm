#include "defs.h"
#include "cuda_defs.h"


#ifdef CONDUCTION
extern __device__ real thermal_diff(real rho, real x1, real x2, real delad);
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
    if (wid ==0) val = warpReduceMin(val);
    return val;
}

__global__ void timestep_kernel(real *cons, real *dx1, real *dx2, real *x1, real *x2, real *out ,int nx1, int nx2, int size_x1, int ntot,int offset, real g, real g1) {
    int i,j,indx;
    real curr_min = FLOATMAX;
    real pres,cs,u1,dt1,u2,dt2,dt;
    
#ifdef CONDUCTION
    real kappa;
    real delad = 1. - 1./g;
#endif
    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot;indx +=blockDim.x*gridDim.x) {
        j = indx/size_x1;
        i = indx -size_x1*j - NGHX1;
        j -= NGHX2;
        if ((i>=0)&&(i<nx1)&&(j>=0)&&(j<nx2)) {

            pres = g1*(cons[indx + 4*ntot]-
                    .5*(  cons[indx + 1*ntot] * cons[indx + 1*ntot]
                        + cons[indx + 2*ntot] * cons[indx + 2*ntot]
                        + cons[indx + 3*ntot] * cons[indx + 3*ntot] )/cons[indx]);

            cs = sqrt( g* pres/cons[indx]);
            u1 = fabs(cons[indx + 1*ntot]/cons[indx]);
            dt1 = dx1[i]/(u1 + cs);
            u2 = fabs(cons[indx + 2*ntot]/cons[indx]);
            dt2 = dx2[j]/(u2 + cs);

            dt = (dt1 < dt2) ? dt1 : dt2;
            
            printf("%d %d dt %e %e %e\n",i,j,dt1,dt2,dt);

#ifdef CONDUCTION
            kappa = thermal_diff(cons[indx],x1[i],x2[j],delad);
            dt1 = dx1[i]*dx1[i]/kappa;
            dt2 = dx2[j]*dx2[j]/kappa;
            dt1 = (dt1 < dt2) ? dt1 : dt2;
            dt = (dt < dt1) ? dt : dt1;
           // printf("%d %d new dt %e %e %e\n",i,j,dt1,dt2,dt);
#endif

            if (dt < curr_min) curr_min = dt;


        }


    }
    printf("%d %e\n",blockIdx.x,curr_min);
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
    if (threadIdx.x ==0) out[blockIdx.x]=cfl*val;
    return;
}
