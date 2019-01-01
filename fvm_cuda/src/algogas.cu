#include "defs.h"
#include "cuda_defs.h"

#define DTMIN 1e-8

__global__ void nancheck_kernel(real *cons, int *out, int ntot,int nf);





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
        GridCons *grid, Parameters *params) {


    int nx1 = grid->nx[0];
    int nx2 = grid->nx[1];
    int nx3 = grid->nx[2];
    int ntot = grid->ntot;
    int size_x1 = grid->size_x1;
    int size_x12 = grid->size_x12;
    int nf = grid->nf;
    int offset = grid->offset;

   /* Add in operator split effects here. */


#ifdef CONDUCTION
    /* Add conduction */
     conduction_flux<<<grid->gridSize_conduction_flux, grid->blockSize_conduction_flux>>>(d_cons,
            d_intenergy,
            d_F_1,
            d_F_2,
            d_F_3,
            d_dx1,
            d_dx2,
            d_dx3,
            d_x1,
            d_x2,
            d_x3,
            params->gamma,
            nx1,
            nx2,
            nx3,
            size_x1,
            size_x12,
            ntot,
            offset,
            nf);
    cudaCheckError();

#endif

#ifdef VISCOSITY
    /* Add viscosity */
    /* Store velocities and divergence in one of
     * the reconstruction arrays
     */
     compute_divergence<<<grid->gridSize_divergence, grid->blockSize_divergence>>>(d_cons,
            d_UL_1,
            d_dx1,
            d_dx2,
            d_dx3,
            d_x1,
            d_x2,
            d_x3,
            nx1,
            nx2,
            nx3,
            size_x1,
            size_x12,
            ntot,
            offset,
            nf);
    cudaCheckError();
    viscous_flux<<<grid->gridSize_viscous_flux, grid->blockSize_viscous_flux>>>(d_UL_1,
    		d_cons,
           d_F_1,
           d_F_2,
           d_F_3,
           d_dx1,
           d_dx2,
           d_dx3,
           d_x1,
           d_x2,
           d_x3,
           nx1,
           nx2,
           nx3,
           size_x1,
           size_x12,
           ntot,
           offset,
           nf);
   cudaCheckError();
#endif
#if defined(CONDUCTION) || defined(VISCOSITY)
   /* Update conservative variables with diffusive fluxes */
	update_cons<<<grid->gridSize_update_cons, grid->blockSize_update_cons>>>(d_cons,
		   d_intenergy,
		   d_F_1,
		   d_F_2,
		   d_F_3,
		   d_dx1,
		   d_dx2,
		   d_dx3,
		   dt,
		   nx1,
		   nx2,
		   nx3,
		   size_x1,
		   size_x12,
		   ntot,
		   offset,
		   nf);
	cudaCheckError();
#endif




    /* X1 reconstruction */
	reconstruct<<<grid->gridSize_reconstruct, grid->blockSize_reconstruct>>>(d_cons ,
		d_UL_1,
		d_UR_1,
		d_dx1,
		1,
		nx1,
		nx2,
		nx3,
		size_x1,
		size_x12,
		nf,
		ntot,
		offset,
		params->gamma-1,
		dt);
	cudaCheckError();
#ifdef POTENTIAL
	source_terms<<<grid->gridSize_source, grid->blockSize_source>>>(d_UL_1 ,
		d_UR_1,
		d_dx1,
		d_x1,
		d_x2,
		d_x3,
		1,
		nx1,
		nx2,
		nx3,
		size_x1,
		size_x12,
		nf,
		ntot,
		offset,
		params->gamma-1,
		dt);
	cudaCheckError();

#endif
	riemann_fluxes<<<grid->gridSize_riemann, grid->blockSize_riemann>>>(d_UL_1 ,
			d_UR_1 ,
			d_F_1 ,
			1,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			nf,
			ntot,
			offset,
			params->gamma);
	cudaCheckError();

	/* x2 reconstruction */
#ifdef DIMS2
	reconstruct<<<grid->gridSize_reconstruct, grid->blockSize_reconstruct>>>(d_cons ,
            d_UL_2,
            d_UR_2,
            d_dx2,
            2,
            nx1,
            nx2,
            nx3,
            size_x1,
            size_x12,
            nf,
            ntot,
            offset,
            params->gamma-1,
            dt);
    cudaCheckError();
#ifdef POTENTIAL
	source_terms<<<grid->gridSize_source, grid->blockSize_source>>>(d_UL_2 ,
		d_UR_2,
		d_dx2,
		d_x1,
		d_x2,
		d_x3,
		2,
		nx1,
		nx2,
		nx3,
		size_x1,
		size_x12,
		nf,
		ntot,
		offset,
		params->gamma-1,
		dt);
	cudaCheckError();
    
#endif
	riemann_fluxes<<<grid->gridSize_riemann, grid->blockSize_riemann>>>(d_UL_2 ,
			d_UR_2 ,
			d_F_2 ,
			2,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			nf,
			ntot,
			offset,
			params->gamma);
	cudaCheckError();
#endif
#ifdef DIMS3
	reconstruct<<<grid->gridSize_reconstruct, grid->blockSize_reconstruct>>>(d_cons ,
            d_UL_3,
            d_UR_3,
            d_dx3,
            3,
            nx1,
            nx2,
            nx3,
            size_x1,
            size_x12,
            nf,
            ntot,
            offset,
            params->gamma-1,
            dt);
    cudaCheckError();
#ifdef POTENTIAL
	source_terms<<<grid->gridSize_source, grid->blockSize_source>>>(d_UL_3 ,
		d_UR_3,
		d_dx3,
		d_x1,
		d_x2,
		d_x3,
		3,
		nx1,
		nx2,
		nx3,
		size_x1,
		size_x12,
		nf,
		ntot,
		offset,
		params->gamma-1,
		dt);
	cudaCheckError();

#endif
	riemann_fluxes<<<grid->gridSize_riemann, grid->blockSize_riemann>>>(d_UL_3 ,
			d_UR_3 ,
			d_F_3 ,
			3,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			nf,
			ntot,
			offset,
			params->gamma);
	cudaCheckError();
#endif


    /* Evolve interface states with transverse fluxes */

#ifdef CTU
#ifdef DIMS2
	transverse_update<<<grid->gridSize_transverse, grid->blockSize_transverse>>>(d_UL_1,
			d_UL_2,
			d_UL_3,
			d_UR_1,
			d_UR_2,
			d_UR_3,
			d_F_1 ,
			d_F_2 ,
			d_F_3 ,
			d_dx1,
			d_dx2,
			d_dx3,
			dt,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			ntot,
			offset,
			nf);
	cudaCheckError();
#ifdef POTENTIAL
	source_transverse_update<<<grid->gridSize_source_transverse, grid->blockSize_source_transverse>>>(d_cons,
			d_UL_1,
			d_UL_2,
			d_UL_3,
			d_UR_1,
			d_UR_2,
			d_UR_3,
			d_F_1 ,
			d_F_2 ,
			d_F_3 ,
			d_dx1,
			d_dx2,
			d_dx3,
			d_x1,
			d_x2,
			d_x3,
			dt,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			ntot,
			offset,
			nf);
	cudaCheckError();
#endif
    /* Compute new fluxes */
	riemann_fluxes<<<grid->gridSize_riemann, grid->blockSize_riemann>>>(d_UL_1 ,
			d_UR_1 ,
			d_F_1 ,
			1,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			nf,
			ntot,
			offset,
			params->gamma);
	cudaCheckError();
	riemann_fluxes<<<grid->gridSize_riemann, grid->blockSize_riemann>>>(d_UL_2 ,
			d_UR_2 ,
			d_F_2 ,
			2,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			nf,
			ntot,
			offset,
			params->gamma);
	cudaCheckError();
#ifdef DIMS3
	riemann_fluxes<<<grid->gridSize_riemann, grid->blockSize_riemann>>>(d_UL_3 ,
				d_UR_3 ,
				d_F_3 ,
				3,
				nx1,
				nx2,
				nx3,
				size_x1,
				size_x12,
				nf,
				ntot,
				offset,
				params->gamma);
		cudaCheckError();
#endif // DIMS3
#endif // DIMS2
#endif // CTU
#ifdef POTENTIAL
	compute_dhalf<<<grid->gridSize_update_cons, grid->blockSize_update_cons>>>(d_cons,
			d_dhalf,
			d_F_1,
			d_F_2,
			d_F_3,
			d_dx1,
			d_dx2,
			d_dx3,
			dt,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			ntot,
			offset,
			nf);
	cudaCheckError();

	update_source<<<grid->gridSize_update_source, grid->blockSize_update_source>>>(d_cons,
			d_dhalf,
			d_F_1,
			d_F_2,
			d_F_3,
			d_dx1,
			d_dx2,
			d_dx3,
			d_x1,
			d_x2,
			d_x3,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			nf,
			ntot,
			offset,dt);
	cudaCheckError();
    
#endif
    /* Final update */
    update_cons<<<grid->gridSize_update_cons, grid->blockSize_update_cons>>>(d_cons,
            d_intenergy,
            d_F_1,
            d_F_2,
            d_F_3,
            d_dx1,
            d_dx2,
            d_dx3,
            dt,
            nx1,
            nx2,
            nx3,
            size_x1,
            size_x12,
            ntot,
            offset,
            nf);
    cudaCheckError();


    return;
}



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
        GridCons *grid, Parameters *params) {


    int nx1 = grid->nx[0];
    int nx2 = grid->nx[1];
    int nx3 = grid->nx[2];
    int ntot = grid->ntot;
    int size_x1 = grid->size_x1;
    int size_x12 = grid->size_x12;
    int nf = grid->nf;
    int offset = grid->offset;
    real dt;
    real h_dt_arr[1024];
    int h_nan_arr[1024];
    int blocks = grid->gridSize_reduc;
//    cudaMemcpy(&grid->cons[-offset],d_cons,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
//	cudaCheckError();
//
//	cudaMemcpy(&grid->intenergy[-offset],d_intenergy,sizeof(real)*ntot,cudaMemcpyDeviceToHost);
//	cudaCheckError();

//	printf("%lg\n",dt);
    /* Calculate new timestep */

    timestep_kernel<<<grid->gridSize_reduc , grid->blockSize_reduc >>>(d_cons,
    		d_dx1,
    		d_dx2,
    		d_dx3,
    		d_x1,
    		d_x2,
    		d_x3,
    		dt_arr,
    		nx1,nx2,nx3,
    		size_x1,size_x12,
            ntot,offset,params->gamma);
    cudaCheckError();

    /* Do final reduction on host */
    cudaMemcpy(h_dt_arr,dt_arr,sizeof(real)*blocks,cudaMemcpyDeviceToHost);
    cudaCheckError();

    dt = FLOATMAX;
    for(int i=0;i<blocks;i++) {
    	if (h_dt_arr[i] < dt) dt = h_dt_arr[i];
    }
    dt *= params->cfl;
//    timestep_kernel_final<<<1,blocks>>>(dt_arr,dt_arr,blocks,params->cfl);
//    cudaCheckError();
//
//    cudaMemcpy(&dt,dt_arr,sizeof(real),cudaMemcpyDeviceToHost);
//    cudaCheckError();

//	curr_min = FLOATMAX;
//	real pres,dt1,cs;
//	for(int i=0;i<nx1;i++) {
//        pres = grid->intenergy[i] * params->gamma_1;
//
//        cs = sqrt( params->gamma* pres/grid->cons[i]);
//        dt1 = grid->dx1[i]/(fabs(grid->cons[i + 1*ntot]/grid->cons[i]) + cs);
//        //printf("%lg\t%lg\t%lg\n",cs,dt1,curr_min);
//        if (dt1 < curr_min) curr_min = dt1;
//	}
 //   printf("dt_max %lg \n",dt_max);

//    if (dt > dt_max) dt = dt_max;
    if (dt < DTMIN){
        printf("Timestep %.4e fell below minimum value of %.1e\n",dt,DTMIN);
        exit(0);
    }




    /* Set boundaries */

    boundary_kernel<<<grid->gridSize_boundary , grid->blockSize_boundary>>>(d_cons,d_intenergy,d_x1,d_x2,d_x3,nx1,nx2,nx3,size_x1,size_x12,nf,ntot,offset,params->gamma,grid->time);
    cudaCheckError();
    cudaDeviceSynchronize();
    cudaCheckError();
    
    /* Check for NaN */
    nancheck_kernel<<<grid->gridSize_reduc , grid->blockSize_reduc >>>(d_cons, nan_arr, ntot,nf);
    cudaCheckError();

    /* Do final reduction on host */
    cudaMemcpy(h_nan_arr,nan_arr,sizeof(int)*blocks,cudaMemcpyDeviceToHost);
    cudaCheckError();

    *nan_res = FALSE;
    for(int i=0;i<blocks;i++) *nan_res |= h_nan_arr[i];

    return dt;
}
__global__ void boundary_kernel(real *cons, real *intenergy, real *x1, real *x2, real *x3,
		int nx1, int nx2, int nx3, int size_x1, int size_x12, int nf, int ntot, int offset, real g, real time) {

    int i,j,k,indxg;
    for(indxg = blockIdx.x*blockDim.x + threadIdx.x; indxg<ntot; indxg+=blockDim.x*gridDim.x) {
    	unpack_indices(indxg,&i,&j,&k,size_x1,size_x12);

        if ((i>=-NGHX1)&&(i<0)&&(j>=-NGHX2)&&(j<nx2+NGHX2)&&(k>=-NGHX3)&&(k<nx3+NGHX3)) {
        /* Lower x1 */
        	x1_boundary_inner(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
        }
        else if ((i>=nx1)&&(i<nx1+NGHX1)&&(j>=-NGHX2)&&(j<nx2+NGHX2)&&(k>=-NGHX3)&&(k<nx3+NGHX3))  {
         /* Upper x1 */
         	x1_boundary_outer(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
        }
#ifdef DIMS2
        else if ((j>=-NGHX2)&&(j<0)&&(i>=-NGHX1)&&(i<nx1+NGHX1)&&(k>=-NGHX3)&&(k<nx3+NGHX3)) {
        /* Lower x2 */
        	x2_boundary_inner(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
        }
        else if ((j>=nx2)&&(j<nx2+NGHX2)&&(i>=-NGHX1)&&(i<nx1+NGHX1)&&(k>=-NGHX3)&&(k<nx3+NGHX3)) {
        /* Upper x2 */
            x2_boundary_outer(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
        }
#endif
#ifdef DIMS3
        else if ((k>=-NGHX2)&&(k<0)&&(i>=-NGHX1)&&(i<nx1+NGHX1)&&(j>=-NGHX2)&&(j<nx2+NGHX2)) {
        /* Lower x3 */
        	x3_boundary_inner(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
        }
        else if ((k>=nx3)&&(k<nx3+NGHX3)&&(i>=-NGHX1)&&(i<nx1+NGHX1)&&(j>=-NGHX2)&&(j<nx2+NGHX2)) {
        /* Upper x3 */
            x3_boundary_outer(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
        }
#endif
    }
    return;
}

__global__ void nancheck_kernel(real *cons, int *out, int ntot,int nf) {
    int indx,n;

    int curr_res = FALSE;

    for(indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot;indx +=blockDim.x*gridDim.x) {
    	for(n=0;n<nf;n++) curr_res |= (cons[indx +n*ntot] != cons[indx + n*ntot]);

    }
    curr_res = blockReduceBoolOR(curr_res);
    if (threadIdx.x ==0) out[blockIdx.x]=curr_res;
    return;
}


__global__ void zero_flux_array(real *F1, real *F2, real *F3, int ntot, int nf) {
	for(int indx = blockIdx.x*blockDim.x + threadIdx.x; indx<ntot; indx+=blockDim.x*gridDim.x) {
		for(int n=0;n<nf;n++) {
			F1[indx + n*ntot] = 0.;
			F2[indx + n*ntot] = 0.;
			F3[indx + n*ntot] = 0.;
		}
	}
	return;
}
