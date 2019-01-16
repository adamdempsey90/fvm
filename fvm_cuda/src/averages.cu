#include "cuda_defs.h"
#include "averages.h"

__device__ real s_func(real rho, real vx1, real vx2, real vx3, real pres,
		real x1, real x2, real x3,
		real dx1, real dx2, real dx3,
		real gamma, real time) {

    return log(pow(pres,1./gamma)/rho);
}

void free_1D_averages(Avg1D *avg_arr, int count) {
	if (count > 0) {
		for(int i=0;i<count;i++) {
			free(avg_arr[i]);
		}
		free(avg_arr);
	}
}
void free_1D_snaps(Snap1D *avg_arr, int count) {
	if (count > 0) {
		for(int i=0;i<count;i++) {
			free(avg_arr[i]);
		}
		free(avg_arr);
	}
}
//void free_2D_averages(Avg2D *avg_arr, int count) {
//	if (count > 0) {
//		for(int i=0;i<count;i++) {
//			free(avg_arr[i]);
//		}
//		free(avg_arr);
//	}
//}
//
//void free_2D_snaps(Snap2D *avg_arr, int count) {
//	if (count > 0) {
//		for(int i=0;i<count;i++) {
//			free(avg_arr[i]);
//		}
//		free(avg_arr);
//	}
//}
//void free_3D_snaps(Snap3D *avg_arr, int count) {
//	if (count > 0) {
//		for(int i=0;i<count;i++) {
//			free(avg_arr[i]);
//		}
//		free(avg_arr);
//	}
//}

void enroll_1D_average(char *name, pfunc func, int avg_dims[3], Avg1D *avg_arr, int *count) {

	Avg1D *avg = (Avg1D *)malloc(sizeof(Avg1D));
	avg->nx1 = grid->nx1;
	avg->nx2 = grid->nx2;
	avg->avgx1 = avg_dims[0];
	avg->avgx2 = avg_dims[1];
	avg->avgx3 = avg_dims[2];
	strcpy(avg->name,name);
	dev_func_ptr = func;
    cudaMemcpyFromSymbol(&(avg->func), dev_func_ptr,sizeof(pfunc));
	cudaCheckError();
	avg_arr[*count] = avg;
	*count += 1;

}

void enroll_1D_snap(char *name, pfunc func, Snap1D *avg_arr, int *count) {
	Snap1D *avg = (Snap1D *)malloc(sizeof(Snap1D));
	avg->nx1 = grid->nx1;
	avg->x1_lims[0] = params->x1_min;
	avg->x1_lims[1] = params->x1_max;

	strcpy(avg->name,name);
	dev_func_ptr = func;
    cudaMemcpyFromSymbol(&(avg->func), dev_func_ptr,sizeof(pfunc));
	cudaCheckError();
	avg_arr[*count] = avg;
	*count += 1;

}


//void enroll_2D_average(char *name, pfunc func, int avg_dims[3], Avg2D *avg_arr, int *count) {
//
//	Avg2D *avg = (Avg2D *)malloc(sizeof(Avg2D));
//	avg->nx1 = grid->nx1;
//	avg->avgx1 = avg_dims[0];
//	avg->avgx2 = avg_dims[1];
//	avg->avgx3 = avg_dims[2];
//	strcpy(avg->name,name);
//	dev_func_ptr = func;
//    cudaMemcpyFromSymbol(&(avg->func), dev_func_ptr,sizeof(pfunc));
//	cudaCheckError();
//	*count += 1;
//	avg_arr[*count-1] = avg;
//}
//void enroll_2D_snap(char *name, pfunc func, Snap2D *avg_arr, int *count) {
//
//	Snap2D *avg = (Snap2D *)malloc(sizeof(Snap2D));
//	avg->nx1 = grid->nx1;
//	avg->nx2 = grid->nx2;
//	strcpy(avg->name,name);
//	dev_func_ptr = func;
//    cudaMemcpyFromSymbol(&(avg->func), dev_func_ptr,sizeof(pfunc));
//	cudaCheckError();
//	*count += 1;
//	avg_arr[*count-1] = avg;
//}
//void enroll_3D_snap(char *name, pfunc func, Snap3D *avg_arr, int *count) {
//
//	Snap3D *avg = (Snap3D *)malloc(sizeof(Snap3D));
//	avg->nx1 = grid->nx1;
//	avg->nx2 = grid->nx2;
//	avg->nx3 = grid->nx3;
//	strcpy(avg->name,name);
//	dev_func_ptr = func;
//    cudaMemcpyFromSymbol(&(avg->func), dev_func_ptr,sizeof(pfunc));
//	cudaCheckError();
//	*count += 1;
//	avg_arr[*count-1] = avg;
//}



void compute_1D_averages(Avg1D *avg_arr, int count,
		real *cons, real *intenergy,
		real *x1, real *x2, real *x3,
		real *dx1, real *dx3, real *dx3,
		GridCons *grid, Parameters *params) {

	/* Compute user defined 1D profiles */

	int i, size;
	Avg1D *avg;
	real *work;
	pfunc *func_ptr;
	int avgx1, avgx2, avgx3;


	if (count > 0) {
		for(i=0 ; i<count ; i++) {
			/* Allocate arrays on host and device */
			avgx1 = FALSE;
			avgx2 = FALSE;
			avgx3 = FALSE;

			avg = &avg_arr[i];
			res = avg->res;
			size = avg->nx1 * avg->nx2;
			func = avg->avg_func;

			res = (real *)malloc(sizeof(real)*size);
			for(i=0;i<size;i++) res[i] = 0.;
			cudaMalloc((void**)&work,sizeof(real)*size);
			cudaCheckError();
            cudaMemcpy(work,res,sizeof(real)*size,cudaMemcpyDeviceToHost);
            cudaCheckError();

			/* Decide which axes to average over */
			if (avg->avgx1) avgx1 = TRUE;
			if (avg->avgx2) avgx2 = TRUE;
			if (avg->avgx3) avgx3 = TRUE;

			/* Launch Kernel */
			average_3D<<<NB[avgx1][avgx2][avgx3], TPB[avgx1][avgx2][avgx3] >>>(cons,intenergy,work,func,avgx1,avgx2,avgx3);

			/* Copy result to host */
            cudaMemcpy(res,work,sizeof(real)*size,cudaMemcpyDeviceToHost);
            cudaCheckError();
			/* Free device pointer */
			cudaFree(work); cudaCheckError();
			/* Done */
		}
	}
	return;
}
void output_1D_averages(char *fname, Avg1D *avg_arr, int count, GridCons *grid, Parameters *params) {
	/* Output user defined 1D profiles */
	int i, size;
	Avg1D *avg;
	real *res;
	if (count > 0) {

		for(i=0 ; i<count ; i++) {
			res = avg->res;
			size = avg->nx1 * avg->nx2;
			/* Write to file */

			/* Free result array for next time */
			free(res); res=NULL;

		}
	}
	return;
}

//
//void compute_2D_averages(Avg2D *avg_arr, int count,
//		real *cons, real *intenergy,
//		real *x1, real *x2, real *x3,
//		real *dx1, real *dx3, real *dx3,
//		GridCons *grid, Parameters *params) {
//
//	/* Compute user defined 2D profiles */
//
//	int i, size;
//	Avg2D *avg;
//	real *work;
//	pfunc *func;
//	int avgx1, avgx2, avgx3;
//
//	if (count > 0) {
//		for(i=0 ; i<count ; i++) {
//			/* Allocate arrays on host and device */
//			avgx1 = FALSE;
//			avgx2 = FALSE;
//			avgx3 = FALSE;
//
//			avg = &avg_arr[i];
//			res = avg->res;
//			size = avg->nx1;
//			func = avg->avg_func;
//
//			res = (real *)malloc(sizeof(real)*size);
//			for(i=0;i<size;i++) res[i] = 0.;
//			cudaMalloc((void**)&work,sizeof(real)*size);
//			cudaCheckError();
//            cudaMemcpy(work,res,sizeof(real)*size,cudaMemcpyDeviceToHost);
//            cudaCheckError();
//
//			/* Decide which axes to average over */
//			if (avg->avgx1) avgx1 = TRUE;
//			if (avg->avgx2) avgx2 = TRUE;
//			if (avg->avgx3) avgx3 = TRUE;
//
//			/* Launch Kernel */
//			average_3D<<<NB[avgx1][avgx2][avgx3], TPB[avgx1][avgx2][avgx3] >>>(cons,intenergy,work,func,avgx1,avgx2,avgx3);
//
//			/* Copy result to host */
//            cudaMemcpy(res,work,sizeof(real)*size,cudaMemcpyDeviceToHost);
//            cudaCheckError();
//			/* Free device pointer */
//			cudaFree(work); cudaCheckError();
//			/* Done */
//		}
//	}
//	return;
//}
//void output_2D_averages(char *fname, Avg2D *avg_arr, int count, GridCons *grid, Parameters *params) {
//	/* Output user defined 2D profiles */
//	int i, size;
//	Avg2D *avg;
//	real *res;
//	if (count > 0) {
//
//		for(i=0 ; i<count ; i++) {
//			res = avg->res;
//			size = avg->nx1;
//			/* Write to file */
//
//			/* Free result array for next time */
//			free(res); res=NULL;
//
//		}
//	}
//	return;
//}


void compute_1D_snapshots(Snap1D *snap_arr, int count,
		real *cons, real *intenergy,
		real *x1, real *x2, real *x3,
		real *dx1, real *dx3, real *dx3,
		GridCons *grid, Parameters *params) {

	/* Compute user defined 1D profiles */

	int i, size;
	Snap1D *snap;
	real *work;
	pfunc *snap_ptr;

	if (count > 0) {
		for(i=0 ; i<count ; i++) {
			/* Allocate arrays on host and device */

			snap = &snap_arr[i];
			res = snap->res;
			size = snap->nx1 ;
			func = snap->snap_func;

			res = (real *)malloc(sizeof(real)*size);
			for(i=0;i<size;i++) res[i] = 0.;
			cudaMalloc((void**)&work,sizeof(real)*size);
			cudaCheckError();
            cudaMemcpy(work,res,sizeof(real)*size,cudaMemcpyDeviceToHost);
            cudaCheckError();


			/* Launch Kernel */
			snap_3D<<< snap_nb, snap_tpb >>>(cons,intenergy,work,func);

			/* Copy result to host */
            cudaMemcpy(res,work,sizeof(real)*size,cudaMemcpyDeviceToHost);
            cudaCheckError();
			/* Free device pointer */
			cudaFree(work); cudaCheckError();
			/* Done */
		}
	}
	return;
}
void output_1D_snapshots(char *fname,Snap1D *snap_arr, int count, GridCons *grid, Parameters *params) {
	/* Output user defined 1D profiles */
	int i, size;
	Snap1D *snap;
	real *res;
	if (count > 0) {

		for(i=0 ; i<count ; i++) {
			res = snap->res;
			size = snap->nx1;
			/* Write to file */

			/* Free result array for next time */
			free(res); res=NULL;

		}
	}
	return;
}

//void compute_2D_snapshots(Snap2D *snap_arr, int count,
//		real *cons, real *intenergy,
//		real *x1, real *x2, real *x3,
//		real *dx1, real *dx3, real *dx3,
//		GridCons *grid, Parameters *params) {
//
//	/* Compute user defined 2D profiles */
//
//	int i, size;
//	Snap2D *snap;
//	real *work;
//	pfunc *snap_ptr;
//
//	if (count > 0) {
//		for(i=0 ; i<count ; i++) {
//			/* Allocate arrays on host and device */
//			snap = &snap_arr[i];
//			res = snap->res;
//			size = snap->nx1 * snap->nx2;
//			func = snap->snap_func;
//
//			res = (real *)malloc(sizeof(real)*size);
//			for(i=0;i<size;i++) res[i] = 0.;
//			cudaMalloc((void**)&work,sizeof(real)*size);
//			cudaCheckError();
//            cudaMemcpy(work,res,sizeof(real)*size,cudaMemcpyDeviceToHost);
//            cudaCheckError();
//
//
//			/* Launch Kernel */
//			snap_3D<<< snap_nb, snap_tpb >>>(cons,intenergy,work,func);
//
//			/* Copy result to host */
//            cudaMemcpy(res,work,sizeof(real)*size,cudaMemcpyDeviceToHost);
//            cudaCheckError();
//			/* Free device pointer */
//			cudaFree(work); cudaCheckError();
//			/* Done */
//		}
//	}
//	return;
//}
//void output_2D_snapshots(char *fname,Snap2D *snap_arr, int count, GridCons *grid, Parameters *params) {
//	/* Output user defined 2D profiles */
//	int i, size;
//	Snap2D *snap;
//	real *res;
//	if (count > 0) {
//
//		for(i=0 ; i<count ; i++) {
//			res = snap->res;
//			size = snap->nx1;
//			/* Write to file */
//
//			/* Free result array for next time */
//			free(res); res=NULL;
//
//		}
//	}
//	return;
//}
//void compute_3D_snapshots(Snap3D *snap_arr, int count,
//		real *cons, real *intenergy,
//		real *x1, real *x2, real *x3,
//		real *dx1, real *dx3, real *dx3,
//		GridCons *grid, Parameters *params) {
//
//	/* Compute user defined 3D profiles */
//
//	int i, size;
//	Snap3D *snap;
//	real *work;
//	pfunc *snap_ptr;
//
//	if (count > 0) {
//		for(i=0 ; i<count ; i++) {
//			/* Allocate arrays on host and device */
//			snap = &snap_arr[i];
//			res = snap->res;
//			size = snap->nx1 * snap->nx2 * snap->nx3;
//			func = snap->snap_func;
//
//			res = (real *)malloc(sizeof(real)*size);
//			for(i=0;i<size;i++) res[i] = 0.;
//			cudaMalloc((void**)&work,sizeof(real)*size);
//			cudaCheckError();
//            cudaMemcpy(work,res,sizeof(real)*size,cudaMemcpyDeviceToHost);
//            cudaCheckError();
//
//
//			/* Launch Kernel */
//			snap_3D<<< snap_nb, snap_tpb >>>(cons,intenergy,work,func);
//
//			/* Copy result to host */
//            cudaMemcpy(res,work,sizeof(real)*size,cudaMemcpyDeviceToHost);
//            cudaCheckError();
//			/* Free device pointer */
//			cudaFree(work); cudaCheckError();
//			/* Done */
//		}
//	}
//	return;
//}
//void output_3D_snapshots(char *fname,Snap3D *snap_arr, int count, GridCons *grid, Parameters *params) {
//	/* Output user defined 3D profiles */
//	int i, size;
//	Snap3D *snap;
//	real *res;
//	if (count > 0) {
//
//		for(i=0 ; i<count ; i++) {
//			res = snap->res;
//			size = snap->nx1;
//			/* Write to file */
//
//			/* Free result array for next time */
//			free(res); res=NULL;
//
//		}
//	}
//	return;
//}


void set_averages_config(int (*NB)[2][2], int (*TPB)[2][2], const int nx, const int ny, const int nz, const int tpb0) {
    int nb, tpb;

    NB[FALSE][FALSE][FALSE] = 1;
    TPB[FALSE][FALSE][FALSE] = 1;

    nb = ny*nz;
    tpb = get_tpb(nx);
    NB[TRUE][FALSE][FALSE] = nb;
    TPB[TRUE][FALSE][FALSE] = tpb;



    tpb  = tpb0;
    nb = (nx*nz + tpb -1)/tpb;
    NB[FALSE][TRUE][FALSE] = nb;
    TPB[FALSE][TRUE][FALSE] = tpb;

    tpb = tpb0;
    nb = (nx*ny + tpb-1)/tpb;
    NB[FALSE][FALSE][TRUE] = nb;
    TPB[FALSE][FALSE][TRUE] = tpb;

    nb = nz;
    tpb = tpb0;
    NB[TRUE][TRUE][FALSE] = nb;
    TPB[TRUE][TRUE][FALSE] = tpb;

    nb = ny;
    tpb = get_tpb(nx);
    NB[TRUE][FALSE][TRUE] = nb;
    TPB[TRUE][FALSE][TRUE] = tpb;

    tpb = get_tpb(nx);
    nb = (nx + tpb-1)/tpb;
    NB[FALSE][TRUE][TRUE] = nb;
    TPB[FALSE][TRUE][TRUE] = tpb;

    tpb = tpb0;
    nb = (nx*ny*nz + tpb-1)/tpb;
    NB[TRUE][TRUE][TRUE] = nb;
    TPB[TRUE][TRUE][TRUE] = tpb;
}



__global__ void sum_3d(const real *in, real *out,
		const real dx1, const real dx2, const real dx3,
		const int nx1, const int nx2, const int nx3,
		const int avgx1, const int avgx2 ,const int avgx3) {
    /* Average over directions indicated by avgx1, avgx2, avgx3 */
    int i,j,k,indx,size, stride;
    real res = 0.;
    /* x1 */
    if (avgx1) {
        /* These require a shared mem reduction */
        if (avgx2) {
            if (avgx3) {
                /* x1 x2 x3 */
                size = nx1*nx2*nx3;
                for(i = blockIdx.x*blockDim.x + threadIdx.x; i<size;i +=blockDim.x*gridDim.x) {
                   res += in[i]*dx1*dx2*dx3;
                }
            }
            else {
                /* x1 x2 */
                size = nx1*nx2;
                stride = size*blockIdx.x;
                for (i = threadIdx.x; i < size; i += blockDim.x) {
                   res += in[ i + stride]*dx1*dx2;
                }
            }
        }
        else {
            if (avgx3) {
                /* x1 x3 */
                size = nx1*nx2;
                for(k=0;k<nx3;k++) {
                    stride = blockIdx.x*nx1 + k*size;
                    for (i = threadIdx.x; i < nx1; i += blockDim.x) {
                        res += in[ i + stride]*dx1*dx3;
                    }
                }
            }
            else {
                /* x1 */
                stride = nx1*blockIdx.x;
                for (i = threadIdx.x; i < nx1; i += blockDim.x) {
                   res += in[ i + stride]*dx1;
                }
            }

        }
    /* Do reductions for each block
     * Note that if all axes are selected then out needs to be summed as well
     */
        res = blockReduceSum(res);
        if (threadIdx.x == 0) out[blockIdx.x] = res;
        return;
    }

    if (avgx2) {
        if (avgx3) {
            /* x2 x3 */
            size = nx2*nx3;
            for(indx = threadIdx.x+blockDim.x*blockIdx.x; indx < nx1; indx += blockDim.x*gridDim.x) {
                res = 0;
                for(j=0;j<size;j++) res += in[indx + j*nx1]*dx2*dx3;
                out[indx] = res;
            }
            return;
        }
        else {
            /* x2 */
            for(indx = threadIdx.x+blockDim.x*blockIdx.x; indx < nx1*nx3; indx += blockDim.x*gridDim.x) {
                k = indx/nx1;
                i = indx - nx1*k;
                stride = i + k*nx1*nx2;
                res = 0;
                for (j = 0; j < nx2; j++) res += in[stride + j*nx1]*dx2;
                out[indx] = res;
            }
            return;
        }
    }

    if (avgx3) {
        /* x3 */
        size = nx1*nx2;
        for(indx = threadIdx.x+blockDim.x*blockIdx.x; indx < size; indx += blockDim.x*gridDim.x) {
            res = 0;
            for (k = 0; k < nx3; k++) res += in[indx + k*size]*dx3;
            out[indx] = res;
        }
        return;
    }

    /* No averaging selected */
    return;
}
