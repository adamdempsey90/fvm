#include "defs.h"
#include "cuda_defs.h"

#define DTMIN 1e-8

//extern "C" {
void algogas_single(real dt,
        real *d_cons,
        real *d_intenergy,
        real *d_UL_1,
        real *d_UR_1,
        real *d_F_1,
        real *d_UL_2,
        real *d_UR_2,
        real *d_F_2,
        real *d_dhalf,
        real *d_dx1,
        real *d_dx2,
        real *d_x1,
        real *d_x2,
        real *dt_arr,
        int blocks,
        int threads,
        GridCons *grid, Parameters *params) {


    int nx1 = grid->nx[0];
    int nx2 = grid->nx[1];
    int ntot = grid->ntot;
    int size_x1 = grid->size_x[0];
    int nf = grid->nf;
    int offset = grid->offset;

    
#ifdef CONDUCTION
     conduction_flux<<<blocks,threads>>>(d_cons,
            d_intenergy,
            d_F_1,
            d_F_2,
            d_dx1,
            d_dx2,
            d_x1,
            d_x2,
            params->gamma,
            dt,
            nx1,
            nx2,
            size_x1,
            ntot,
            offset,
            nf);
    cudaCheckError();
    conduction_update<<<blocks,threads>>>(d_cons,
            d_intenergy,
            d_F_1,
            d_F_2,
            d_dx1,
            d_dx2,
            dt,
            nx1,
            nx2,
            size_x1,
            ntot,
            offset,
            nf);
    cudaCheckError();
#endif




    if (nx1 > 1) {
        plm<<<blocks,threads>>>(d_cons ,
            d_UL_1, 
            d_UR_1,
            d_dx1,
            1,
            nx1,
            nx2,
            size_x1,
            nf,
            ntot,
            offset,
            params->gamma_1,
            dt);
        cudaCheckError();
#ifdef POTENTIAL
        source_terms<<<blocks,threads>>>(d_UL_1 ,
            d_UR_1,
            d_dx1,
            d_x1,
            d_x2,
            1,
            nx1,
            nx2,
            size_x1,
            nf,
            ntot,
            offset,
            params->gamma_1,
            dt);
        cudaCheckError();
    
#endif
        riemann_fluxes<<<blocks,threads>>>(d_UL_1 ,
                d_UR_1 ,
                d_F_1 ,
                1,
                nx1,
                nx2,
                size_x1,
                nf,
                ntot,
                offset,
                params->gamma);
        cudaCheckError();
    }
    if (nx2 > 1) {

    plm<<<blocks,threads>>>(d_cons ,
            d_UL_2,
            d_UR_2,
            d_dx2,
            2,
            nx1,
            nx2,
            size_x1,
            nf,
            ntot,
            offset,
            params->gamma_1,
            dt);
    cudaCheckError();
#ifdef POTENTIAL
        source_terms<<<blocks,threads>>>(d_UL_2 ,
            d_UR_2,
            d_dx2,
            d_x1,
            d_x2,
            2,
            nx1,
            nx2,
            size_x1,
            nf,
            ntot,
            offset,
            params->gamma_1,
            dt);
        cudaCheckError();
    
#endif
        riemann_fluxes<<<blocks,threads>>>(d_UL_2 ,
                d_UR_2 ,
                d_F_2 ,
                2,
                nx1,
                nx2,
                size_x1,
                nf,
                ntot,
                offset,
                params->gamma);
        cudaCheckError();
    }


    /* Evolve interface states with transverse fluxes */

#ifdef CTU
    if ((nx1 > 1)&&(nx2 > 1)) {
        transverse_update<<<blocks,threads>>>(d_UL_1,
                d_UL_2,
                d_UR_1,
                d_UR_2,
                d_F_1 ,
                d_F_2 ,
                d_dx1,
                d_dx2,
                dt,
                nx1,
                nx2,
                size_x1,
                ntot,
                offset,
                nf);
        cudaCheckError();
        source_transverse_update<<<blocks,threads>>>(d_cons,
                d_UL_1,
                d_UL_2,
                d_UR_1,
                d_UR_2,
                d_F_1 ,
                d_F_2 ,
                d_dx1,
                d_dx2,
                d_x1,
                d_x2,
                dt,
                nx1,
                nx2,
                size_x1,
                ntot,
                offset,
                nf);
        cudaCheckError();
    /* Compute new fluxes */
        riemann_fluxes<<<blocks,threads>>>(d_UL_1 , 
                d_UR_1 ,
                d_F_1 ,
                1,
                nx1,
                nx2,
                size_x1,
                nf,
                ntot,
                offset,
                params->gamma);
        cudaCheckError();
        riemann_fluxes<<<blocks,threads>>>(d_UL_2 ,
                d_UR_2 ,
                d_F_2 ,
                2,
                nx1,
                nx2,
                size_x1,
                nf,
                ntot,
                offset,
                params->gamma);
        cudaCheckError();
    }
#endif 
#ifdef POTENTIAL
    compute_dhalf<<<blocks,threads>>>(d_cons,
            d_dhalf,
            d_F_1,
            d_F_2,
            d_dx1,
            d_dx2,
            dt,
            nx1,
            nx2,
            size_x1,
            ntot,
            offset,
            nf);
    update_source<<<blocks,threads>>>(d_cons,
            d_dhalf,
            d_F_1,
            d_F_2,
            d_dx1,
            d_dx2,
            d_x1,
            d_x2,
            nx1,
            nx2,
            size_x1,
            nf,
            ntot,
            offset,dt);
        cudaCheckError();
    
#endif
    /* Final update */
    update_cons<<<blocks,threads>>>(d_cons,
            d_intenergy,
            d_F_1,
            d_F_2,
            d_dx1,
            d_dx2,
            dt,
            nx1,
            nx2,
            size_x1,
            ntot,
            offset,
            nf);
    cudaCheckError();

   
    




    return;
}

real algogas_dt(real dt, real dtout, int threads, int blocks, GridCons *grid, FluxCons *fluxes, Parameters *params) {
    real end_time = grid->time + dtout;
    real dt_max;
    int ntot = grid->ntot;
    int nf = grid->nf;
    int size_x1 = grid->size_x[0];
    int size_x2 = grid->size_x[1];
    int offset = grid->offset;


    real *d_cons, *d_intenergy;
    real *d_F_1, *d_UL_1, *d_UR_1;
    real *d_F_2, *d_UL_2, *d_UR_2;
    real *d_dx1, *d_dx2;
    real *d_x1, *d_x2;
    real *dt_arr;
    real *d_dhalf;
    
    
    cudaMalloc((void**)&d_dx1,sizeof(real)*size_x1);
	cudaCheckError();
	cudaMemcpy(d_dx1,&grid->dx1[-NGHX1],sizeof(real)*size_x1,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_dx2,sizeof(real)*size_x2);
	cudaCheckError();
	cudaMemcpy(d_dx2,&grid->dx2[-NGHX2],sizeof(real)*size_x2,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_x1,sizeof(real)*size_x1);
	cudaCheckError();
	cudaMemcpy(d_x1,&grid->xc1[-NGHX1],sizeof(real)*size_x1,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_x2,sizeof(real)*size_x2);
	cudaCheckError();
	cudaMemcpy(d_x2,&grid->xc2[-NGHX2],sizeof(real)*size_x2,cudaMemcpyHostToDevice);
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

	cudaMalloc((void**)&d_dhalf,sizeof(real)*ntot);
	cudaCheckError();


	cudaMalloc((void**)&dt_arr,sizeof(real)*blocks);
	cudaCheckError();

   
    
    while (grid->time < end_time) { 
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
            d_dhalf,
            d_dx1 + NGHX1,
            d_dx2 + NGHX2,
            d_x1  + NGHX1,
            d_x2  + NGHX2,
            dt_arr,
            blocks,
            threads,
            grid, params);
    /* Set new timestep and bcs */
        grid->time += dt;
        dt_max = end_time - grid->time;
        dt = set_bc_timestep(dt_max, 
                d_cons,
                d_intenergy,
                d_dx1 + NGHX1,
                d_dx2 + NGHX2,
                d_x1  + NGHX1,
                d_x2  + NGHX2,
                dt_arr,
                blocks,
                threads,
                grid,params);

    }

    /* Copy to host */
	cudaMemcpy(&grid->cons[-offset],d_cons,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
    cudaCheckError();

	cudaMemcpy(&grid->intenergy[-offset],d_intenergy,sizeof(real)*ntot,cudaMemcpyDeviceToHost);
    cudaCheckError();

    /* Free device arrays */
    cudaFree(d_cons);
    cudaFree(d_intenergy);
    cudaFree(d_F_1);
    cudaFree(d_UL_1);
    cudaFree(d_UR_1);
    cudaFree(d_F_2);
    cudaFree(d_UL_2);
    cudaFree(d_UR_2);
    cudaFree(d_dx1);
    cudaFree(d_dx2);
    cudaFree(d_x1);
    cudaFree(d_x2);
    cudaFree(dt_arr);
    cudaFree(d_dhalf);

    return dt;
}
real algogas_firststep(real dtout, int threads, int blocks, GridCons *grid, FluxCons *fluxes, Parameters *params) {
    real end_time = grid->time + dtout;
    real dt_max,dt;
    int ntot = grid->ntot;
    int nf = grid->nf;
    int size_x1 = grid->size_x[0];
    int size_x2 = grid->size_x[1];
    int offset = grid->offset;


    real *d_cons, *d_intenergy;
    real *d_F_1, *d_UL_1, *d_UR_1;
    real *d_F_2, *d_UL_2, *d_UR_2;
    real *d_dx1, *d_dx2;
    real *d_x1, *d_x2;
    real *dt_arr;
    real *d_dhalf;


    cudaMalloc((void**)&d_dx1,sizeof(real)*size_x1);
	cudaCheckError();
	cudaMemcpy(d_dx1,&grid->dx1[-NGHX1],sizeof(real)*size_x1,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_dx2,sizeof(real)*size_x2);
	cudaCheckError();
	cudaMemcpy(d_dx2,&grid->dx2[-NGHX2],sizeof(real)*size_x2,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_x1,sizeof(real)*size_x1);
	cudaCheckError();
	cudaMemcpy(d_x1,&grid->xc1[-NGHX1],sizeof(real)*size_x1,cudaMemcpyHostToDevice);
	cudaCheckError();

	cudaMalloc((void**)&d_x2,sizeof(real)*size_x2);
	cudaCheckError();
	cudaMemcpy(d_x2,&grid->xc2[-NGHX2],sizeof(real)*size_x2,cudaMemcpyHostToDevice);
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

	cudaMalloc((void**)&d_dhalf,sizeof(real)*ntot);
	cudaCheckError();


	cudaMalloc((void**)&dt_arr,sizeof(real)*blocks);
	cudaCheckError();

   

    dt_max = end_time - grid->time;

    /* Set new timestep and bcs */

    dt = set_bc_timestep(dt_max, 
        d_cons,
        d_intenergy,
        d_dx1 + NGHX1,
        d_dx2 + NGHX2,
        d_x1  + NGHX1,
        d_x2  + NGHX2,
        dt_arr,
        blocks,
        threads,
        grid,params);

    /* Take one step */

    algogas_single(dt,
        d_cons,
        d_intenergy,
        d_UL_1,
        d_UR_1,
        d_F_1,
        d_UL_2,
        d_UR_2,
        d_F_2,
        d_dhalf,
        d_dx1 + NGHX1,
        d_dx2 + NGHX2,
        d_x1 + NGHX1,
        d_x2 + NGHX2,
        dt_arr,
        blocks,
        threads,
        grid, params);

    grid->time += dt;
    /* Get new timestep */
    dt_max = end_time - grid->time;

    dt = set_bc_timestep(dt_max,
        d_cons,
        d_intenergy,
        d_dx1 + NGHX1,
        d_dx2 + NGHX2,
        d_x1  + NGHX1,
        d_x2  + NGHX2,
        dt_arr,
        blocks,
        threads,
        grid,params);
    
    /* Copy results to host */
	cudaMemcpy(&grid->cons[-offset],d_cons,sizeof(real)*ntot*nf,cudaMemcpyDeviceToHost);
    cudaCheckError();

	cudaMemcpy(&grid->intenergy[-offset],d_intenergy,sizeof(real)*ntot,cudaMemcpyDeviceToHost);
    cudaCheckError();

    /* Free device arrays */
    cudaFree(d_cons); cudaCheckError();
    cudaFree(d_intenergy);cudaCheckError();
    cudaFree(d_F_1); cudaCheckError();
    cudaFree(d_UL_1); cudaCheckError();
    cudaFree(d_UR_1); cudaCheckError();
    cudaFree(d_F_2); cudaCheckError();
    cudaFree(d_UL_2); cudaCheckError();
    cudaFree(d_UR_2); cudaCheckError();
    cudaFree(d_dx1); cudaCheckError();
    cudaFree(d_dx2); cudaCheckError();
    cudaFree(d_x1); cudaCheckError();
    cudaFree(d_x2); cudaCheckError();
    cudaFree(dt_arr); cudaCheckError();
    cudaFree(d_dhalf); cudaCheckError();

    return dt;
}
//}


real set_bc_timestep(real dt_max, 
        real *d_cons,
        real *d_intenergy,
        real *d_dx1,
        real *d_dx2,
        real *d_x1,
        real *d_x2,
        real *dt_arr,
        int blocks,
        int threads,
        GridCons *grid, Parameters *params) {


    int nx1 = grid->nx[0];
    int nx2 = grid->nx[1];
    int ntot = grid->ntot;
    int size_x1 = grid->size_x[0];
    int nf = grid->nf;
    int offset = grid->offset;
    real dt;

    /* Calculate new timestep */
    timestep_kernel<<<blocks,threads>>>(d_cons,d_dx1,d_dx2,d_x1,d_x2,dt_arr,nx1,nx2,size_x1,
            ntot,offset,params->gamma,params->gamma_1);
    cudaCheckError();
    
    timestep_kernel_final<<<1,blocks>>>(dt_arr,dt_arr,blocks,params->cfl);
    cudaCheckError();

    cudaMemcpy(&dt,dt_arr,sizeof(real),cudaMemcpyDeviceToHost);
    cudaCheckError();
   
    if (dt < DTMIN){
        printf("Timestep %.4e fell below minimum value of %.1e\n",dt,DTMIN);
        //exit(1);
    }
    if (dt > dt_max) dt = dt_max;



    /* Set boundaries */

    boundary_kernel<<<blocks,threads>>>(d_cons,d_intenergy,d_x1,d_x2,nx1,nx2,size_x1,nf,ntot,offset,params->gamma,grid->time);
    cudaCheckError();
    
    return dt;
}
