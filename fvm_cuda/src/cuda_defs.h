#include <cuda.h>
#include <cuda_runtime.h>

#define cudaCheckError() {                                          \
 cudaError_t cer=cudaGetLastError();                                 \
 if(cer!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(cer));           \
   exit(0); \
 }                                                                 \
}




#define FULL_MASK 0xffffffff


#define GINDEX(i,j,k) ( (offset) + (i) + (j)*size_x1 + (k)*size_x12)


#ifdef DIMS3
#define VOL(i,j,k) dx1[i]*dx2[j]*dx3[k]
#define AREA_12(i,j,k) dx1[i]*dx2[j]
#define AREA_13(i,j,k) dx1[i]*dx3[k]
#define AREA_23(i,j,k) dx2[j]*dx3[k]
#else
#ifdef DIMS2
#define VOL(i,j,k) dx1[i]*dx2[j]
#define AREA_12(i,j,k) dx1[i]
#define AREA_13(i,j,k) dx1[i]
#define AREA_23(i,j,k) dx2[j]

#else
#define VOL(i,j,k) dx1[i]
#define AREA_12(i,j,k) 1 
#define AREA_13(i,j,k) 1 
#define AREA_23(i,j,k) 1 
#endif
#endif




__inline__ __device__ void unpack_indices(int indx, int *i, int *j, int *k, int size_x1, int size_x12) {
	*k = indx/size_x12;
    *j = (indx-(*k)*size_x12)/size_x1;
    *i = indx -size_x1*(*j) -size_x12*(*k) - NGHX1;
    *j -= NGHX2;
    *k -= NGHX3;
    return;
}
__inline __device__ real MIN2(real a, real b) {return (a<b) ? a : b;}
__inline __device__ real MAX2(real a, real b) {return (a>b) ? a : b;}
__inline __device__ real MIN3(real a, real b, real c) {return (a<b) ? ( (a<c) ? a : c ) : ( (b<c) ? b : c );}
__inline __device__ real MAX3(real a, real b, real c) {return (a>b) ? ( (a>c) ? a : c ) : ( (b>c) ? b : c );}


/* driver.cu */
 void driver(GridCons *grid, Parameters *params);

/* algogas.cu */
real set_bc_timestep(real dt_max,
        real *d_cons,
        real *d_intenergy,
        real *d_dx1,
        real *d_dx2,
        real *d_dx3,
        real *d_x1,
        real *d_x2,
        real *d_x3,
        real *dt_arr,
        int *nan_arr,
        int *nan_res,
        GridCons *grid, Parameters *params);
void algogas_single(real dt,
        real *d_cons,
        real *d_intenergy,
        real *d_UL_1,
        real *d_UR_1,
        real *d_F_1,
        real *d_UL_2,
        real *d_UR_2,
        real *d_F_2,
        real *d_UL_3,
        real *d_UR_3,
        real *d_F_3,
        real *d_dhalf,
        real *d_dx1,
        real *d_dx2,
        real *d_dx3,
        real *d_x1,
        real *d_x2,
        real *d_x3,
        real *dt_arr,
        GridCons *grid, Parameters *params);

__global__ void zero_flux_array(real *F1, real *F2, real *F3, int ntot, int nf);

/* update.cu */
__global__ void compute_dhalf(real *cons, real *dhalf, real *F_1, real *F_2,real *F_3,
        real *dx1, real *dx2, real *dx3, real dt, int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);

__global__ void update_cons(real *cons, real *intenergy, real *F_1, real *F_2, real *F_3,
        real *dx1, real *dx2, real *dx3, real dt, int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);

__global__ void transverse_update(real *UL_1, real *UL_2, real *UL_3,
        real *UR_1, real *UR_2, real *UR_3,
        real *F_1, real *F_2, real *F_3, real *dx1, real *dx2, real *dx3, real dt,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);
__global__ void cons_to_prim(real *cons, real *intenergy, real *prim, real g1,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);
__global__ void prim_to_cons(real *cons, real *intenergy, real *prim, real g1,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);

/* riemann.cu */
__global__ void riemann_fluxes(real *UL, real *UR, real *F, 
        int dir1,int nx1, int nx2, int nx3, int size_x1, int size_x12,
        int nf,int ntot, int offset, real g);

__device__ int exact_sample(const real dL, const real uL, const real pL, const real aL,
        const real dR, const real uR, const real pR, const real aR,
        real *d_f, real *u_f, real *p_f, 
        real g, real g1, real g2, real g3, real g4, real g5,
        real S, real tol);
__device__ void anrs(const real dL,const real uL, const real pL, const real aL,
        const real dR, const real uR, const real pR, const real aR, 
        real *ps, real *us, real g, real g1,real g2,real g3,real g4);

/* reconstruct.cu */
__global__ void reconstruct(real *cons, real *UL, real *UR, real *dx,
        int dir1,int nx1, int nx2, int nx3, int size_x1, int size_x12,
        int nf,int ntot, int offset, real g1, real dt);

/* calc_dt.cu */
__global__ void timestep_kernel_final(real *in, real *out ,int n, real cfl);
__global__ void timestep_kernel(real *cons, real *dx1, real *dx2,real *dx3, real *x1, real *x2,real *x3, real *out ,
		int nx1, int nx2, int nx3, int size_x1, int size_x12,int ntot,int offset, real g);


/* hllc.cu */
__device__ void hllc(real dL, real uL, real pL, real aL,
        real dR, real uR, real pR, real aR,
        real *SL, real *SR, real *Sstar,
        real g, real g1, real g2, real g3, real g4, real g5);

/* source.cu */
#ifdef POTENTIAL
__global__ void update_source(real *cons, real *dhalf, real *F1, real *F2, real *F3,
		real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int nf,int ntot, int offset, real dt);
__global__ void source_terms(real *UL, real *UR, real *dx, real *x1, real *x2, real *x3,
        int dir1,int nx1, int nx2, int nx3, int size_x1, int size_x12,
        int nf,int ntot, int offset, real g1, real dt);
__global__ void source_transverse_update(real *cons, real *UL_1, real *UL_2, real *UL_3, real *UR_1, real *UR_2, real *UR_3,
		real *F_1, real *F_2, real *F_3, real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3, real dt,
		int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);
#endif

/* conduction.cu */
#ifdef CONDUCTION
__global__ void conduction_flux(real *cons, real *intenergy, real *F_1, real *F_2, real *F_3,
        real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3, real g,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);

#endif

/* viscosity.cu */
#ifdef VISCOSITY
__global__ void compute_divergence(real *cons, real *vel,
        real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);
__global__ void viscous_flux(real *vel, real *rho, real *F_1, real *F_2, real *F_3,
        real *dx1, real *dx2, real *dx3, real *x1, real *x2, real *x3,
        int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);
#endif

/* prob.cu */
__device__ void x1_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void x2_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void x3_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void x1_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void x2_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void x3_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);

/* reduction.cu */


/* Boundary conditions */
__global__ void boundary_kernel(real *cons, real *intenergy, real *x1, real *x2, real *x3,
		int nx1, int nx2, int nx3, int size_x1, int size_x12, int nf, int ntot, int offset, real g, real time);

/* boundary.cu */
__device__ void reflecting_boundary_inner(int dir,int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void reflecting_boundary_outer(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void outflow_boundary_inner(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void outflow_boundary_outer(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void periodic_boundary_inner(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void periodic_boundary_outer(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void ic_boundary(int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);


/* config.cu */
void config_kernels(int *threads, int *blocks, GridCons *grid, Parameters *params);

/* Reduction operations */
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

    val = (threadIdx.x < blockDim.x / 32) ? shared[lane] : 1e99;
    if (wid ==0) val = warpReduceMin(val);
    return val;
}

__inline__ __device__ real warpReduceSum(real val) {
    for(int offset = 16; offset >0; offset /= 2) {
#if __CUDACC_VER_MAJOR__ < 9
        val +=__shfl_down(val,offset);
#else
        val += __shfl_down_sync(FULL_MASK,val,offset);
#endif


    }
    return val;
}
__inline__ __device__ real blockReduceSum(real val) {
    static __shared__ real shared[32];
    int lane = threadIdx.x % 32;
    int wid = threadIdx.x / 32;

    val = warpReduceSum(val);

    if (lane == 0) shared[wid] = val;
    __syncthreads();

    val = (threadIdx.x < blockDim.x / 32) ? shared[lane] : 0.;
    if (wid ==0) val = warpReduceSum(val);

    return val;
}



__inline__ __device__ int warpReduceBoolOR(int val) {
    for(int offset = 16; offset >0; offset /= 2) {
#if __CUDACC_VER_MAJOR__ < 9
        val |=__shfl_down(val,offset);
#else
        val |= __shfl_down_sync(FULL_MASK,val,offset);
#endif


    }
    return val;
}
__inline__ __device__ int blockReduceBoolOR(int val) {
    static __shared__ int shared[32];
    int lane = threadIdx.x % 32;
    int wid = threadIdx.x / 32;

    val = warpReduceBoolOR(val);

    if (lane == 0) shared[wid] = val;
    __syncthreads();

    val = (threadIdx.x < blockDim.x / 32) ? shared[lane] : 0;
    if (wid ==0) val = warpReduceBoolOR(val);

    return val;
}





