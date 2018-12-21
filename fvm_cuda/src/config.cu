#include "defs.h"
#include "cuda_defs.h"


#define OccuPre(kern) { \
	cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,kern, 0, 0);              \
	gridSize = (ntot + blockSize - 1) / blockSize;                  \
}
#define OccuPost(kern,name) { \
	cudaDeviceSynchronize();           \
	cudaOccupancyMaxActiveBlocksPerMultiprocessor( &maxActiveBlocks,kern, blockSize,0);        \
	cudaGetDevice(&device);           \
	cudaGetDeviceProperties(&props, device);           \
	occupancy = (maxActiveBlocks * blockSize / props.warpSize) / (float)(props.maxThreadsPerMultiProcessor / props.warpSize);           \
	printf("*****************************************************************\n");           \
	printf("* Occupancy results for %s kernel:                        *\n",name);           \
	printf("*    MinGridsize: %d, gridsize: %d, maxactiveblocks: %d          *\n",minGridSize,gridSize,maxActiveBlocks);            \
	printf("*    <<<%d,%d>>>                                                *\n",gridSize,blockSize);           \
	printf("*    Launched blocks of size %d. Theoretical occupancy: %.3f  *\n",blockSize, occupancy);           \
	printf("*****************************************************************\n");           \
}


#define MAXBLOCKS 1024 // From shared mem requirements in reductions

/*
void set_threads_blocks(int n, int *threads, int *blocks) {
	*threads = MAXTHREADS;
	*blocks = (n+*threads-1)/(*threads);
    if (*blocks > MAXBLOCKS) *blocks = MAXBLOCKS;
	return;

}
*/
void config_kernels(int *threads, int *blocks, GridCons *grid, Parameters *params)
{ 
	int blockSize;   // The launch configurator returned block size
	int minGridSize; // The minimum grid size needed to achieve the
				   // maximum occupancy for a full device launch
	int gridSize;    // The actual grid size needed, based on input size
	int nx1 = grid->nx[0];
	int nx2 = grid->nx[1];
	int nx3 = grid->nx[2];
    int size_x2 = grid->size_x2;
    int size_x3 = grid->size_x3;


	int size_x12 = grid->size_x12;

	int ntot = grid->ntot;
	int size_x1 = grid->size_x1;
	int nf = grid->nf;
	int offset = grid->offset;
	real dt = 1e-6;

	real *d_cons, *d_intenergy;
	real *d_F_1, *d_UL_1, *d_UR_1;
	real *d_F_2, *d_UL_2, *d_UR_2;
	real *d_F_3, *d_UL_3, *d_UR_3;
	real *d_dx1, *d_dx2, *d_dx3;
	real *d_x1, *d_x2, *d_x3;
	real *d_dhalf;



	cudaMalloc((void**)&d_dx1,sizeof(real)*size_x1);
	cudaCheckError();
	cudaMemcpy(d_dx1,&grid->dx1[-NGHX1],sizeof(real)*size_x1,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_dx2,sizeof(real)*size_x2);
	cudaCheckError();
	cudaMemcpy(d_dx2,&grid->dx2[-NGHX2],sizeof(real)*size_x2,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_dx3,sizeof(real)*size_x3);
	cudaCheckError();
	cudaMemcpy(d_dx3,&grid->dx3[-NGHX3],sizeof(real)*size_x3,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_x1,sizeof(real)*size_x1);
	cudaCheckError();
	cudaMemcpy(d_x1,&grid->xc1[-NGHX1],sizeof(real)*size_x1,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_x2,sizeof(real)*size_x2);
	cudaCheckError();
	cudaMemcpy(d_x2,&grid->xc2[-NGHX2],sizeof(real)*size_x2,cudaMemcpyHostToDevice);
	cudaCheckError();
	cudaMalloc((void**)&d_x3,sizeof(real)*size_x3);
	cudaCheckError();
	cudaMemcpy(d_x3,&grid->xc3[-NGHX3],sizeof(real)*size_x3,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_cons,sizeof(real)*ntot*nf);
	cudaCheckError();
	cudaMemcpy(d_cons,&grid->cons[-offset],sizeof(real)*ntot*nf,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_intenergy,sizeof(real)*ntot);
	cudaCheckError();
	cudaMemcpy(d_intenergy,&grid->intenergy[-offset],sizeof(real)*ntot,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_UL_1,sizeof(real)*ntot*nf);
	cudaCheckError();

	cudaMalloc((void**)&d_UR_1,sizeof(real)*ntot*nf);
	cudaCheckError();

	cudaMalloc((void**)&d_F_1,sizeof(real)*ntot*nf);
	cudaCheckError();

	cudaMalloc((void**)&d_UL_2,sizeof(real)*ntot*nf);
	cudaCheckError();

	cudaMalloc((void**)&d_UR_2,sizeof(real)*ntot*nf);
	cudaCheckError();

	cudaMalloc((void**)&d_F_2,sizeof(real)*ntot*nf);
	cudaCheckError();

	cudaMalloc((void**)&d_UL_3,sizeof(real)*ntot*nf);
	cudaCheckError();

	cudaMalloc((void**)&d_UR_3,sizeof(real)*ntot*nf);
	cudaCheckError();

	cudaMalloc((void**)&d_F_3,sizeof(real)*ntot*nf);
	cudaCheckError();
   	cudaMalloc((void**)&d_dhalf,sizeof(real)*ntot);
   	cudaCheckError();
	cudaMemcpy(d_dhalf,d_cons,sizeof(real)*ntot,cudaMemcpyDeviceToDevice);
	cudaCheckError();

	int device, maxActiveBlocks;
	float occupancy;
	cudaDeviceProp props;

	OccuPre(boundary_kernel)
    boundary_kernel<<<gridSize, blockSize>>>(d_cons,
    		d_intenergy,
            d_x1  + NGHX1,
            d_x2  + NGHX2,
            d_x3  + NGHX3,
    		nx1,nx2,nx3,
    		size_x1,size_x12,
    		nf,ntot,offset,
    		params->gamma,0.0);
    cudaCheckError();
    OccuPost(boundary_kernel,"boundary_kernel")
	grid->gridSize_boundary = gridSize;
	grid->blockSize_boundary  = blockSize;

#ifdef CONDUCTION
    /* Add conduction */
	OccuPre(conduction_flux)
     conduction_flux<<< gridSize, blockSize>>>(d_cons,
            d_intenergy,
            d_F_1,
            d_F_2,
            d_F_3,
            d_dx1 + NGHX1,
            d_dx2 + NGHX2,
            d_dx3 + NGHX3,
            d_x1  + NGHX1,
            d_x2  + NGHX2,
            d_x3  + NGHX3,
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
    OccuPost(conduction_flux,"conduction_flux")
	grid->gridSize_conduction_flux = gridSize;
	grid->blockSize_conduction_flux  = blockSize;

#endif

#ifdef VISCOSITY
    /* Add viscosity */
    /* Store velocities and divergence in one of
     * the reconstruction arrays
     */
	OccuPre(compute_divergence)

     compute_divergence<<< gridSize, blockSize>>>(d_cons,
            d_UL_1,
            d_dx1 + NGHX1,
            d_dx2 + NGHX2,
            d_dx3 + NGHX3,
            d_x1  + NGHX1,
            d_x2  + NGHX2,
            d_x3  + NGHX3,
            nx1,
            nx2,
            nx3,
            size_x1,
            size_x12,
            ntot,
            offset,
            nf);
    cudaCheckError();
    OccuPost(compute_divergence,"divergence")
	grid->gridSize_divergence = gridSize;
	grid->blockSize_divergence  = blockSize;

	OccuPre(viscous_flux)

    viscous_flux<<< gridSize, blockSize>>>(d_UL_1,
    		d_cons,
           d_F_1,
           d_F_2,
           d_F_3,
           d_dx1 + NGHX1,
           d_dx2 + NGHX2,
           d_dx3 + NGHX3,
           d_x1  + NGHX1,
           d_x2  + NGHX2,
           d_x3  + NGHX3,
           nx1,
           nx2,
           nx3,
           size_x1,
           size_x12,
           ntot,
           offset,
           nf);
   cudaCheckError();
   OccuPost(viscous_flux,"viscous_flux")
	grid->gridSize_viscous_flux = gridSize;
	grid->blockSize_viscous_flux = blockSize;
#endif


    /* X1 reconstruction */
	OccuPre(plm)
	plm<<< gridSize, blockSize>>>(d_cons ,
		d_UL_1,
		d_UR_1,
		d_dx1 +NGHX1,
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
	OccuPost(plm,"plm x1")
	grid->gridSize_plm = gridSize;
	grid->blockSize_plm  = blockSize;

#ifdef POTENTIAL
	OccuPre(source_terms)

	source_terms<<< gridSize, blockSize>>>(d_UL_1 ,
		d_UR_1,
		d_dx1 + NGHX1,
        d_x1  + NGHX1,
        d_x2  + NGHX2,
        d_x3  + NGHX3,
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
	OccuPost(source_terms,"source_terms x1")
	grid->gridSize_source = gridSize;
	grid->blockSize_source  = blockSize;


#endif
	OccuPre(riemann_fluxes)

	riemann_fluxes<<< gridSize, blockSize>>>(d_UL_1 ,
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
	OccuPost(riemann_fluxes,"riemann_fluxes x1")
	grid->gridSize_riemann = gridSize;
	grid->blockSize_riemann  = blockSize;




    /* Evolve interface states with transverse fluxes */

#ifdef CTU
#ifdef DIMS2
	OccuPre(transverse_update)

	transverse_update<<< gridSize, blockSize>>>(d_UL_1,
			d_UL_2,
			d_UL_3,
			d_UR_1,
			d_UR_2,
			d_UR_3,
			d_F_1 ,
			d_F_2 ,
			d_F_3 ,
			d_dx1 + NGHX1,
			d_dx2 + NGHX2,
			d_dx3 + NGHX3,
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
	OccuPost(transverse_update,"transverse_update")
	grid->gridSize_transverse = gridSize;
	grid->blockSize_transverse = blockSize;

#ifdef POTENTIAL
	OccuPre(source_transverse_update)

	source_transverse_update<<< gridSize, blockSize>>>(d_cons,
			d_UL_1,
			d_UL_2,
			d_UL_3,
			d_UR_1,
			d_UR_2,
			d_UR_3,
			d_F_1 ,
			d_F_2 ,
			d_F_3 ,
			d_dx1 + NGHX1,
			d_dx2 + NGHX2,
			d_dx3 + NGHX3,
			d_x1 + NGHX1,
			d_x2 + NGHX2,
			d_x3 + NGHX3,
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
	OccuPost(source_transverse_update,"source_transverse")
	grid->gridSize_source_transverse = gridSize;
	grid->blockSize_source_transverse = blockSize;
#endif
#endif // DIMS2
#endif // CTU
#ifdef POTENTIAL
	OccuPre(update_source)

	update_source<<< gridSize, blockSize>>>(d_cons,
			d_dhalf,
			d_F_1,
			d_F_2,
			d_F_3,
			d_dx1 + NGHX1,
			d_dx2 + NGHX2,
			d_dx3 + NGHX3,
			d_x1 + NGHX1,
			d_x2 + NGHX2,
			d_x3 + NGHX3,
			nx1,
			nx2,
			nx3,
			size_x1,
			size_x12,
			nf,
			ntot,
			offset,dt);
	cudaCheckError();
	OccuPost(update_source,"update_source ")
	grid->gridSize_update_source = gridSize;
	grid->blockSize_update_source = blockSize;


#endif
    /* Final update */
	OccuPre(update_cons)

    update_cons<<< gridSize, blockSize>>>(d_cons,
            d_intenergy,
            d_F_1,
            d_F_2,
            d_F_3,
			d_dx1 + NGHX1,
			d_dx2 + NGHX2,
			d_dx3 + NGHX3,
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
	OccuPost(update_cons,"update_cons ")
	grid->gridSize_update_cons = gridSize;
	grid->blockSize_update_cons = blockSize;

	/* Take update_cons as default */
	*blocks = gridSize;
	*threads = blockSize;




	OccuPre(timestep_kernel)
	real *dt_arr;
	if (gridSize > MAXBLOCKS) gridSize = MAXBLOCKS;
   	cudaMalloc((void**)&dt_arr,sizeof(real)*blockSize);
   	cudaCheckError();
    timestep_kernel<<<gridSize, blockSize>>>(d_cons,
			d_dx1 + NGHX1,
			d_dx2 + NGHX2,
			d_dx3 + NGHX3,
			d_x1 + NGHX1,
			d_x2 + NGHX2,
			d_x3 + NGHX3,
    		dt_arr,
    		nx1,nx2,nx3,size_x1,size_x12,
            ntot,offset,params->gamma);
    cudaCheckError();
	OccuPost(timestep_kernel,"timestep")
	grid->gridSize_reduc = gridSize;
	grid->blockSize_reduc = blockSize;



    cudaFree(d_cons); cudaCheckError();
    cudaFree(d_intenergy); cudaCheckError();
    cudaFree(d_F_1); cudaCheckError();
    cudaFree(d_UL_1); cudaCheckError();
    cudaFree(d_UR_1); cudaCheckError();
    cudaFree(d_F_2); cudaCheckError();
    cudaFree(d_UL_2); cudaCheckError();
    cudaFree(d_UR_2); cudaCheckError();
    cudaFree(d_F_3); cudaCheckError();
    cudaFree(d_UL_3); cudaCheckError();
    cudaFree(d_UR_3); cudaCheckError();
    cudaFree(d_dx1); cudaCheckError();
    cudaFree(d_dx2); cudaCheckError();
    cudaFree(d_dx3); cudaCheckError();
    cudaFree(d_x1); cudaCheckError();
    cudaFree(d_x2); cudaCheckError();
    cudaFree(d_x3); cudaCheckError();
    cudaFree(dt_arr); cudaCheckError();
    cudaFree(d_dhalf); cudaCheckError();


	return;
}

