#include "cuda_defs.h"


__global__ void sum_2d_over_x1(const real *in, real *out, const real dx, const int nx,const int ny) {
 	 /* Sum along first index in in[j][i] = in[i + nx*j]
 	  * Each block is responsible for a fixed j, and we use a parallel
 	  * reduction of the threads from i=0,nx-1
 	  * output is 1D array out[j]
 	  * NB, TPB = ny, nx
 	  */
 	  __shared__ real sdata[1024];


 	  sdata[threadIdx.x] = 0;
 	  for (int i = threadIdx.x; i < nx; i += blockDim.x) {
 			  sdata[threadIdx.x] += in[ i + nx*blockIdx.x]*dx;
 	  }
 	  __syncthreads();
 	  for (int i = blockDim.x>>1; i > 0; i>>=1){
 		if (threadIdx.x < i) sdata[threadIdx.x] += sdata[threadIdx.x+i];
 		__syncthreads();}
 	  if (threadIdx.x == 0) out[blockIdx.x] = sdata[0];
 }

__global__ void sum_2d_over_x2(const real *in, real *out, const real dx, const int nx, const int ny) {
	 /* Sum along first index for in[j][i] = in[i + nx*j]
	  * Each thread is responsible for a fixed  i, and
	  * sums through all j values.
	  * out is a 1D array out[i]
	  * NB, TPB = nx/TPB, ny
	  */
	  for(int indx = threadIdx.x+blockDim.x*blockIdx.x; indx < nx; indx += blockDim.x*gridDim.x) {
		  real temp = 0;
		  for (int j = 0; j < ny; j++) temp += in[indx + j*nx]*dx;
		  out[indx] = temp;
	  }
}

 /*
  *  Reductions over 1 dimension, 3D array -> 2D array
  */

 __global__ void sum_3d_over_x1(const real *in, real *out, const real dx, const int nx,const int ny,const int nz) {
 	 /* Sum along first index in in[k][j][i] = in[i + nx*j + nx*ny*k]
 	  * Each block is responsible for a fixed (k,j), and we use a parallel
 	  * reduction of the threads from i=0,nx-1
 	  * output is 2D array out[k][j]
 	  */
 	  __shared__ real sdata[1024];


 	  sdata[threadIdx.x] = 0;
 	  for (int i = threadIdx.x; i < nx; i += blockDim.x) {
 			  sdata[threadIdx.x] += in[ i + nx*blockIdx.x]*dx;
 	  }
 	  __syncthreads();
 	  for (int i = blockDim.x>>1; i > 0; i>>=1){
 		if (threadIdx.x < i) sdata[threadIdx.x] += sdata[threadIdx.x+i];
 		__syncthreads();}
 	  if (threadIdx.x == 0) out[blockIdx.x] = sdata[0];
 }
 __global__ void sum_3d_over_x3(const real *in, real *out, const real dx, const int nx, const int ny, const int nz) {

	 /* Sum along third index for in[k][j][i] = in[i + nx*j + nx*ny*k]
	  * Each thread is responsible for a fixed  (i,j), and
	  * sums through all k values.
	  * out is a 2D array out[j][i] = out[i + j*nx]
	  */
      int size = nx*ny;
	  for(int indx = threadIdx.x+blockDim.x*blockIdx.x; indx < size; indx += blockDim.x*gridDim.x) {
		  real temp = 0;
		  for (int k = 0; k < nz; k++) temp += in[indx + k*size]*dx;
		  out[indx] = temp;
	  }
}
__global__ void sum_3d_over_x2(const real *in, real *out,const real dx, const int nx, const int ny, const int nz) {
	 /* Sum along second index for in[k][j][i] = in[i + nx*j + nx*ny*k]
	  * Each thread is responsible for a fixed  (i,k), and
	  * sums through all k values.
	  * out is a 2D array out[k][i] = out[i + nx*k]
	  */
    int i,k;
	  for(int indx = threadIdx.x+blockDim.x*blockIdx.x; indx < nx*nz; indx += blockDim.x*gridDim.x) {
          k = indx/nx;
          i = indx - nx*k;
		  real temp = 0;
		  for (int j = 0; j < ny; j++) temp += in[i + k*nx*ny + j*nx]*dx;
		  out[indx] = temp;
	  }
}

/*
 * Reductions over 2 dimensions, 3D array -> 1D array
 *
 */
__global__ void sum_3d_over_x2_x3(const real *in, real *out, const real dx, const int nx, const int ny, const int nz) {
	 /* Sum along second and third indices for in[k][j][i] = in[i + nx*j + nx*ny*k]
	  * Each thread is responsible for a fixed  i, and
	  * sums through all k values.
	  * out is a 1D array out[i]
	  */
	  for(int i = threadIdx.x+blockDim.x*blockIdx.x; i < nx; i += blockDim.x*gridDim.x) {
		  real temp = 0;
		  for (int k = 0; k < nz; k++) {
			  for(int j=0; j< ny; j++) {
				  temp += in[i + j*nx + k*nx*ny]*dx;
			  }
		  }
		  out[i] = temp;
	  }
}
__global__ void sum_3d_over_x1_x2(const real *in, real *out, const real dx, const int nx, const int ny, const int nz) {
	 /* Sum along first and second index in in[k][j][i] = in[i + nx*j + nx*ny*k]
	  * Each block is responsible for a fixed k, and we use a parallel
	  * reduction of the threads from i=0,nx-1 and j=0,ny-1
	  * output is 1D array out[k]
	  */
	  __shared__ real sdata[1024];
	  sdata[threadIdx.x] = 0;
	  for (int i = threadIdx.x; i < nx*ny; i += blockDim.x) {
			  sdata[threadIdx.x] += in[ i + blockIdx.x * nx*ny]*dx;
	  }
	  __syncthreads();
	  for (int i = blockDim.x>>1; i > 0; i>>=1){
		if (threadIdx.x < i) sdata[threadIdx.x] += sdata[threadIdx.x+i];
		__syncthreads();}
	  if (threadIdx.x == 0) out[blockIdx.x] = sdata[0];
}
__global__ void sum_3d_over_x1_x3(const real *in, real *out, const real dx, const int nx, const int ny, const int nz) {
	 /* Sum along first and third index in in[k][j][i] = in[i + nx*j + nx*ny*k]
	  * Each block is responsible for a fixed j, and we use a parallel
	  * reduction of the threads from i=0,nx-1 and k=0,nz-1
	  * output is 1D array out[j]
	  */
	  __shared__ real sdata[1024];
	  sdata[threadIdx.x] = 0;
	  for (int i = threadIdx.x; i < nx; i += blockDim.x) {
		  for(int k = 0 ; k < nz ; k++) {
			  sdata[threadIdx.x] += in[ i + k*nx*ny + blockIdx.x *nx]*dx;
		  }
	  }
	  __syncthreads();
	  for (int i = blockDim.x>>1; i > 0; i>>=1){
		if (threadIdx.x < i) sdata[threadIdx.x] += sdata[threadIdx.x+i];
		__syncthreads();}
	  if (threadIdx.x == 0) out[blockIdx.x] = sdata[0];
}
