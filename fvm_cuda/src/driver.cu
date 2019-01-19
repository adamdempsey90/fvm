#include "defs.h"
#include "cuda_defs.h"
#include <ctype.h>
#include <time.h>
#include <unistd.h>



int check_termination(void) {
    if ( access( "STOP", F_OK ) != -1 ) return TRUE;
    return FALSE;
}
int check_dump(void) {
    if ( access( "DUMP", F_OK ) != -1 ) return TRUE;
    return FALSE;
}

void driver(GridCons *grid, Parameters *params) {
    struct timespec tic, toc; 
    real size_in_gb;
    real elapsed, elapsed_sec, elapsed_nsec, rate;
    real dt, dt_min, end_time, dt_max; 
    real tstart = grid->time;
    real tstop = params->tend;
    char fname[512];
    int dump_count = 0;



    /* Setup step counters */
    long long  unsigned int nsteps, nsteps_dt;
    nsteps = 0;



    int n_0d = 1; // volume averaged outputs
    int n_1d = 1; // 1D outputs
    real t_0d, t_1d;
    real dt_0d = (tstop - tstart)/(real)(params->nout0d);
    if (dt_0d < 0) dt_0d = 1e99; // No 0d outputs
    else printf("Time between %d, 0D outputs is %.3e\n",params->nout0d,dt_0d);

    real dt_1d = (tstop - tstart)/(real)(params->nout1d);
    if (dt_1d < 0) dt_1d = 1e99; // No 1d outputs
    else printf("Time between %d, 1D outputs is %.3e\n",params->nout1d,dt_1d);


#ifdef DIMS2
    int n_2d = 1; // 2D outputs
    real t_2d;
    real dt_2d = (tstop - tstart)/(real)(params->nout2d);
    if (dt_2d < 0) dt_2d = 1e99; // No 2d outputs
    else printf("Time between %d, 2D outputs is %.3e\n",params->nout2d,dt_2d);

#endif
#ifdef DIMS3
    int n_3d = 1; // 3D outputs
    real t_3d;
    real dt_3d = (params->tend - tstart)/(real)(params->nout3d);
    if (dt_3d < 0) dt_3d = 1e99; // No 3d outputs
    else printf("Time between %d, 3D outputs is %.3e\n",params->nout3d,dt_3d);

#endif






    /* Allocate GPU arrays */
    int ntot = grid->ntot;
    int size_x1 = grid->size_x1;
    int size_x2 = grid->size_x2;
    int size_x3 = grid->size_x3;
    int nf = grid->nf;
    int offset = grid->offset;

    int nan_res;
    real *d_cons, *d_intenergy;
    real *d_F_1, *d_UL_1, *d_UR_1;
    real *d_F_2, *d_UL_2, *d_UR_2;
    real *d_F_3, *d_UL_3, *d_UR_3;
    real *d_dx1, *d_dx2, *d_dx3;
    real *d_x1, *d_x2, *d_x3;
    real *dt_arr;
    real *d_dhalf;
    int *nan_arr;
    
    clock_gettime(CLOCK_MONOTONIC, &tic);
    
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


	cudaMalloc((void**)&dt_arr,sizeof(real)*(grid->gridSize_reduc));
	cudaCheckError();


	cudaMalloc((void**)&nan_arr,sizeof(int)*(grid->gridSize_reduc));
	cudaCheckError();

    clock_gettime(CLOCK_MONOTONIC, &toc);
    elapsed_sec = (double)(toc.tv_sec - tic.tv_sec) + 1e-9*(toc.tv_nsec-tic.tv_nsec);
#ifndef SILENT
    printf("Device allocation took %.2e ms\n", elapsed_sec*1e3);
#endif

    printf("Evolving from %.3e to %.3e\n",tstart,tstop);
    dt_max = tstop - tstart;
    /* Set initial bc and get first timestep */
    dt = set_bc_timestep(dt_max, 
            d_cons,
            d_intenergy,
            d_dx1 + NGHX1,
            d_dx2 + NGHX2,
            d_dx3 + NGHX3,
            d_x1  + NGHX1,
            d_x2  + NGHX2,
            d_x3  + NGHX3,
            dt_arr,
            nan_arr,
            &nan_res,
            grid,params);

    sprintf(fname, "%s_%d",params->outputname, 0);
#ifdef DIMS3
    snapshot_3d(fname, grid, params);
#else
#ifdef DIMS2
    snapshot_2d(fname, grid, params);
#else
    snapshot_1d(fname, grid, params);
#endif
#endif
    /* Start main loop 
       Use do-while to ensure at least one step is exectuted
     */
    do {
        /* compute when next outputs are */
    	t_0d = tstart + n_0d * dt_0d;
    	t_1d = tstart + n_1d * dt_1d;

#ifdef DIMS2
    	t_2d = tstart + n_2d * dt_2d;

#endif
#ifdef DIMS3
    	t_3d = tstart + n_3d * dt_3d;
#endif

        /* Compute smallest time until nearest outptu */
    	dt_min = fmin(t_0d - grid->time, t_1d - grid->time);
#ifdef DIMS2
        dt_min = fmin(dt_min, t_2d - grid->time);
#endif
#ifdef DIMS3
        dt_min = fmin(dt_min, t_3d - grid->time);
#endif


        /* Measure time to evolve dt_min */
    	clock_gettime(CLOCK_MONOTONIC, &tic);
    	nsteps_dt = 0;

        /* Evolve for dt_min */
        end_time = grid->time + dt_min;
        while ((grid->time < end_time)&&(nsteps < params->maxsteps)) {
    	    /* Zero flux arrays */
    	    zero_flux_array<<<grid->gridSize_update_cons, grid->blockSize_update_cons>>>(d_F_1,d_F_2,d_F_3,ntot,nf);
    	cudaCheckError();

            /* Ensure we end exactly on output time */
            dt_max = end_time - grid->time;
            if (dt > dt_max) dt = dt_max;

            /* Advance by dt */
            algogas_single(dt, 
                d_cons,
                d_intenergy,
                d_UL_1,
                d_UR_1,
                d_F_1,
                d_UL_2,
                d_UR_2,
                d_F_2,
                d_UL_3,
                d_UR_3,
                d_F_3,
                d_dhalf,
                d_dx1 + NGHX1,
                d_dx2 + NGHX2,
                d_dx3 + NGHX3,
                d_x1  + NGHX1,
                d_x2  + NGHX2,
                d_x3  + NGHX3,
                dt_arr,
                grid, params);

            grid->time += dt;
            nsteps += 1;
            nsteps_dt += 1;

            
            /* Set new bc and timestep */
            dt = set_bc_timestep(dt_max, 
                    d_cons,
                    d_intenergy,
                    d_dx1 + NGHX1,
                    d_dx2 + NGHX2,
                    d_dx3 + NGHX3,
                    d_x1  + NGHX1,
                    d_x2  + NGHX2,
                    d_x3  + NGHX3,
                    dt_arr,
                    nan_arr,
                    &nan_res,
                    grid,params);
            if (dt < 0) {
                /* NaN */
            	printf("NaN detected in domain!\n");
            	fflush(stdout);
            	exit(-1);
            }
        } 

        /* Measure time that took */
        //dt_curr = algogas_dt(dt_curr, dt_min, &nsteps_elapsed,grid,params);
        clock_gettime(CLOCK_MONOTONIC, &toc);
        elapsed_sec = (double)(toc.tv_sec - tic.tv_sec) + 1e-9*(toc.tv_nsec-tic.tv_nsec);
        elapsed = elapsed_sec/nsteps_dt;
        rate = (real)(ntot) / elapsed;
#ifndef SILENT
        printf("Step %.2e \t Time %.2e \t avg DT %.2e\n\tTook %.3e s for %llu steps, avg. %.3e zones per sec.\n", (double)nsteps, grid->time, dt_min/(real)nsteps_dt,elapsed_sec,nsteps_dt, rate);
        fflush(stdout);
#endif

        /* Check for STOP, DUMP */

        if (check_termination()) {
            printf("Detected STOP file ...\n");
            cudaMemcpy(&grid->cons[-offset],d_cons,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();

            cudaMemcpy(&grid->prim[-offset],d_UL_1,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();


            cudaMemcpy(&grid->intenergy[-offset],d_intenergy,sizeof(real)*ntot,cudaMemcpyDeviceToHost);
            cudaCheckError();
            sprintf(fname, "dump_%d",dump_count);
            dump_count++;
        	snapshot(fname, grid, params);
            remove("STOP");
            grid->time = 1e99;
            break;

        }
        if (check_dump()) {
            printf("Detected DUMP file ...\n");
            cudaMemcpy(&grid->cons[-offset],d_cons,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();

            cudaMemcpy(&grid->prim[-offset],d_UL_1,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();


            cudaMemcpy(&grid->intenergy[-offset],d_intenergy,sizeof(real)*ntot,cudaMemcpyDeviceToHost);
            cudaCheckError();
            sprintf(fname, "dump_%d",dump_count);
            dump_count++;
        	snapshot(fname, grid, params);
            remove("DUMP");
        }

        /* Convert to primitives */

        cons_to_prim<<<grid->gridSize_update_cons, grid->blockSize_update_cons>>>(d_cons,d_intenergy,d_UL_1,params->gamma-1,
    	 grid->nx[0],grid->nx[1],grid->nx[2], size_x1, grid->size_x12, ntot, offset, nf);
        cudaCheckError();

        /* Check which output to do */
        if (grid->time >= t_0d) {
        	//output_0d(n_0d, dt_curr, grid, params);
        	n_0d += 1;
        }
        if (grid->time >= t_1d) {
#ifndef DIMS2
            /* 1d outputs are snapshots*/
            clock_gettime(CLOCK_MONOTONIC, &tic);
            cudaMemcpy(&grid->cons[-offset],d_cons,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();

            cudaMemcpy(&grid->prim[-offset],d_UL_1,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();


            cudaMemcpy(&grid->intenergy[-offset],d_intenergy,sizeof(real)*ntot,cudaMemcpyDeviceToHost);
            cudaCheckError();
            clock_gettime(CLOCK_MONOTONIC, &toc);
            elapsed_sec = (double)(toc.tv_sec - tic.tv_sec) + 1e-9*(toc.tv_nsec-tic.tv_nsec);
            size_in_gb = sizeof(real)*ntot*(2*nf+1)/1e9;
            rate = size_in_gb/elapsed_sec;
#ifndef SILENT
            printf("Device to host copy of %.2e GB took %.2e ms, %.2e GB/s\n", elapsed_sec*1e3, size_in_gb,rate);
        
#endif
            sprintf(fname, "%s_%d",params->outputname, n_1d);
            snapshot(fname,grid,params);
#else
            /* DIMS > 1 so 1d outputs are averages */
        	//output_1d(n_1d, dt_curr, grid, params);
#endif
        	n_1d += 1;
        }
#ifdef DIMS2
        if (grid->time  >= t_2d) {
#ifndef DIMS3
            /* 2d outputs ares snapshots*/
            clock_gettime(CLOCK_MONOTONIC, &tic);
            cudaMemcpy(&grid->cons[-offset],d_cons,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();

            cudaMemcpy(&grid->prim[-offset],d_UL_1,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();


            cudaMemcpy(&grid->intenergy[-offset],d_intenergy,sizeof(real)*ntot,cudaMemcpyDeviceToHost);
            cudaCheckError();
            clock_gettime(CLOCK_MONOTONIC, &toc);
            elapsed_sec = (double)(toc.tv_sec - tic.tv_sec) + 1e-9*(toc.tv_nsec-tic.tv_nsec);
            size_in_gb = sizeof(real)*ntot*(2*nf+1)/1e9;;;;
            rate = size_in_gb/elapsed_sec;
#ifndef SILENT
            printf("Device to host copy of %.2e GB took %.2e ms, %.2e GB/s\n", elapsed_sec*1e3, size_in_gb,rate);
        
#endif
            sprintf(fname, "%s_%d",params->outputname, n_2d);
            snapshot(fname,grid,params);
#else
            /* DIMS > 2 so 2d outputs are averages */

            //output_2d(n_2d, dt_curr,grid,params);
#endif
        	n_2d += 1;
        }
#endif
#ifdef DIMS3
        if (grid->time  >= t_3d) {
            /* 3d outputs are snapshots*/
            clock_gettime(CLOCK_MONOTONIC, &tic);
            cudaMemcpy(&grid->cons[-offset],d_cons,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();

            cudaMemcpy(&grid->prim[-offset],d_UL_1,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
            cudaCheckError();


            cudaMemcpy(&grid->intenergy[-offset],d_intenergy,sizeof(real)*ntot,cudaMemcpyDeviceToHost);
            cudaCheckError();
            clock_gettime(CLOCK_MONOTONIC, &toc);
            elapsed_sec = (double)(toc.tv_sec - tic.tv_sec) + 1e-9*(toc.tv_nsec-tic.tv_nsec);
            size_in_gb = sizeof(real)*ntot*(2*nf+1)/1e9;;;;
            rate = size_in_gb/elapsed_sec;
#ifndef SILENT
            printf("Device to host copy of %.2e GB took %.2e ms, %.2e GB/s\n", elapsed_sec*1e3, size_in_gb,rate);
        
#endif
            sprintf(fname, "%s_%d",params->outputname, n_3d);
            snapshot(fname,grid,params);

            //output_3d(n_3d, dt_curr,grid,params);
        	n_3d += 1;
        }
#endif

    /*End of main loop */
    } while ((grid->time <= params->tend) && (nsteps < params->maxsteps ));

    if (nsteps >= params->maxsteps) {
        printf("Hit maximum number of steps (%.2e)\n",(double)params->maxsteps);
        sprintf(fname, "dump_%d",dump_count);
        dump_count++;
        snapshot(fname,grid,params);
    }

    /* Free GPU arrays */
    cudaFree(d_cons); cudaCheckError();
    cudaFree(d_intenergy);cudaCheckError();
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
    cudaFree(nan_arr); cudaCheckError();

    return;
}


