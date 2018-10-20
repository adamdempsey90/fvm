#include <cuda.h>
#include <cuda_runtime.h>

#define POTENTIAL 

#define cudaCheckError() {                                          \
 cudaError_t cer=cudaGetLastError();                                 \
 if(cer!=cudaSuccess) {                                              \
   printf("Cuda failure %s:%d: '%s'\n",__FILE__,__LINE__,cudaGetErrorString(cer));           \
   exit(0); \
 }                                                                 \
}

#define FULL_MASK 0xffffffff

#define GINDEX(i,j) ( (offset) + (i) + (j)*size_x1)

__global__ void compute_dhalf(real *cons, real *dhalf, real *F_1, real *F_2,
        real *dx1, real *dx2, real dt, int nx1, int nx2, int size_x1, int ntot, int offset, int nf);
__global__ void update_cons(real *cons, real *intenergy, real *F_1, real *F_2,real *dx1, real *dx2, real dt, int nx, int ny, int size_x1, int ntot,int offset, int nf);
__global__ void source_transverse_update(real *UL_1, real *UL_2,
        real *UR_1, real *UR_2, 
        real *F_1, real *F_2, real *dx1, real *dx2, real dt, int nx1, int nx2, int size_x1, int ntot, int offset, int nf);
__global__ void riemann_fluxes(real *UL, real *UR, real *F, 
        int dir1,int nx1, int nx2, int size_x1, 
        int nf,int ntot, int offset, real g);
__global__ void transverse_update(real *UL_1, real *UL_2,
        real *UR_1, real *UR_2, 
        real *F_1, real *F_2, real *dx1, real *dx2, real dt, int nx, int ny, int size_x1, int ntot, int offset, int nf);

__device__ int exact_sample(const real dL, const real uL, const real pL, const real aL,
        const real dR, const real uR, const real pR, const real aR,
        real *d_f, real *u_f, real *p_f, 
        real g, real g1, real g2, real g3, real g4, real g5,
        real S, real tol);
__device__ void anrs(const real dL,const real uL, const real pL, const real aL,
        const real dR, const real uR, const real pR, const real aR, 
        real *ps, real *us, real g, real g1,real g2,real g3,real g4);
__global__ void plm(real *cons, real *UL, real *UR, real *dx,
        int dir1,int nx1, int nx2, int size_x1, 
        int nf,int ntot, int offset, real g1, real dt);
__global__ void timestep_kernel_final(real *in, real *out ,int n, real cfl);
__global__ void timestep_kernel(real *cons, real *dx1, real *dx2, real *out ,int nx1, int nx2, int size_x1, int ntot,int offset, real g, real g1);
__global__ void boundary_kernel(real *cons, real *intenergy, real *x1, real *x2, int nx1, int nx2, int size_x1, int nf, int ntot, int offset, real g, real time);

__device__ void hllc(real dL, real uL, real pL, real aL,
        real dR, real uR, real pR, real aR,
        real *SL, real *SR, real *Sstar,
        real g, real g1, real g2, real g3, real g4, real g5);
__global__ void source_terms(real *UL, real *UR, real *dx, real *x1, real *x2,
        int dir1,int nx1, int nx2, int size_x1, 
        int nf,int ntot, int offset, real g1, real dt);


/* Boundary conditions */
__device__ void reflecting_boundary_x2_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time);
__device__ void reflecting_boundary_x1_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void reflecting_boundary_x2_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void reflecting_boundary_x1_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void outflow_boundary_x2_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void outflow_boundary_x1_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void outflow_boundary_x2_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void outflow_boundary_x1_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void periodic_boundary_x2_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void periodic_boundary_x1_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void periodic_boundary_x2_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void periodic_boundary_x1_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void ic_boundary_x2_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void ic_boundary_x1_inner(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void ic_boundary_x2_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
__device__ void ic_boundary_x1_outer(int indxg, int i, int j, real *cons, real *intenergy, int nx1, int nx2, int ntot, int nf, int size_x1, int offset, real g, real time) ;
