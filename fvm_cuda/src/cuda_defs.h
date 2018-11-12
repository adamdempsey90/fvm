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


/* update.cu */
__global__ void compute_dhalf(real *cons, real *dhalf, real *F_1, real *F_2,real *F_3,
        real *dx1, real *dx2, real *dx3, real dt, int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);

__global__ void update_cons(real *cons, real *intenergy, real *F_1, real *F_2, real *F_3,
        real *dx1, real *dx2, real *dx3, real dt, int nx1, int nx2, int nx3, int size_x1, int size_x12, int ntot, int offset, int nf);

__global__ void transverse_update(real *UL_1, real *UL_2, real *UL_3,
        real *UR_1, real *UR_2, real *UR_3,
        real *F_1, real *F_2, real *F_3, real *dx1, real *dx2, real *dx3, real dt,
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
__global__ void plm(real *cons, real *UL, real *UR, real *dx,
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

/* Boundary conditions */
/* boundary.cu */
__device__ void reflecting_boundary_inner(int dir,int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void reflecting_boundary_outer(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void outflow_boundary_inner(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void outflow_boundary_outer(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void periodic_boundary_inner(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void periodic_boundary_outer(int dir, int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);
__device__ void ic_boundary(int indxg, int i, int j,int k, real *cons, real *intenergy, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time);