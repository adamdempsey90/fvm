int NB[2][2][2], TPB[2][2][2];

typedef real (*pfunc)(int indx, real *prim, real *intenergy,
		real x1, real x2, real x3,
		real dx1, real dx2, real dx3,
		real gamma, real time);


//typedef struct Avg1D {
//    int avgx1,avgx2,avgx3;
//    int nx1,nx2;
//    real *res;
//    pfunc func;
//    char *name;
//} Avg1D;

typedef struct Avg0D {
    real res;
    pfunc func;
    char *name;
} Avg0D;


Avg0D* allocate_Avg0D(char *name, pfunc *func) {
	Avg0D *avg = (Avg0D *)malloc(sizeof(Avg0D));
	strcpy(avg->name,name);
	avg->func = func;
	avg->res = 0.;
	return avg;
}


void enroll_Avg0D(Avg0D *avgs, Avg0D *newavg, int *count) {
	avgs[*count] = newavg;
	*count++;
}

void compute_avg123(const real *prim, pfunc *func, real *res) {
	real ans =0.;
	avg_123<<<TPB,NPB>>>(prim, out, avg->func);
	for(int i=0;i<NPB;i++) ans += out[i];
	*res = ans;
	return;
}

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
__global__ void avg_3d(const real *prim, real *out, pfunc *func,
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
                   res += func(i,prim)*dx1*dx2*dx3;
                }
            }
            else {
                /* x1 x2 */
                size = nx1*nx2;
                stride = size*blockIdx.x;
                for (i = threadIdx.x; i < size; i += blockDim.x) {
                   res += func( i + stride,prim)*dx1*dx2;
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
                        res += func( i + stride,prim)*dx1*dx3;
                    }
                }
            }
            else {
                /* x1 */
                stride = nx1*blockIdx.x;
                for (i = threadIdx.x; i < nx1; i += blockDim.x) {
                   res += func( i + stride,prim)*dx1;
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
                for(j=0;j<size;j++) res += func(indx + j*nx1,prim)*dx2*dx3;
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
                for (j = 0; j < nx2; j++) res += func(stride + j*nx1,prim)*dx2;
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
            for (k = 0; k < nx3; k++) res += func(indx + k*size,prim)*dx3;
            out[indx] = res;
        }
        return;
    }

    /* No averaging selected */
    return;
}
