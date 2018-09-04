#include "defs.h"

#define DTMIN 1e-8

real calc_dt(GridCons *grid, Parameters *params);
void algogas_single(real dt_max, GridCons *grid, FluxCons *fluxes, Parameters *params) {
    int i,ntot;
    real dt;
    ntot = grid->ntot;
    dt = calc_dt(grid,params);
    if (dt < DTMIN){
        printf("Timestep %.4e fell below minimum value of %.1e\n",dt,DTMIN);
        exit(1);
    }
    dt = fmin(dt,dt_max);


    /* Set boundaries */

    if (grid->nx[0]>1) {
        set_boundary_x1(grid,params);
    }
    if (grid->nx[1]>1) { 
        set_boundary_x2(grid,params);
    }

    /* Reconstruction pt. 1 */
    /* Reconstruct and then evolve boundary values for dt/2 */
    reconstruct(grid,fluxes,params,dt);


    
    /* Riemann solve with boundary values */   
    if (grid->nx[0] > 1) {
        riemann_fluxes(fluxes->UL_1, fluxes->UR_1, fluxes->Fstar_1, 1,
                grid->nx,
                grid->size_x[0],
                0,
                grid->nf,
                ntot,
                params->gamma);
    }
    if (grid->nx[1] > 1) {
        riemann_fluxes(fluxes->UL_2, fluxes->UR_2, fluxes->Fstar_2, 2,
                grid->nx,
                grid->size_x[0],
                0,
                grid->nf,
                ntot,
                params->gamma);
    }


    /* Evolve interface states with transverse fluxes */

#ifdef CTU
    if ((grid->nx[0] > 1)&&(grid->nx[1]>1)) {
        transverse_update(grid,fluxes,params,dt);
    }

    /* Compute new fluxes */
    if (grid->nx[0] > 1) {
        riemann_fluxes(fluxes->UL_1, fluxes->UR_1, fluxes->Fstar_1, 1,
                grid->nx,
                grid->size_x[0],
                0,
                grid->nf,
                ntot,
                params->gamma);
    }
    if (grid->nx[1] > 1) {
        riemann_fluxes(fluxes->UL_2, fluxes->UR_2, fluxes->Fstar_2, 2,
                grid->nx,
                grid->size_x[0],
                0,
                grid->nf,
                ntot,
                params->gamma);
    }

#endif 
    /* Final update */
    update_cons(grid,fluxes,params,dt);


    grid->time += dt;

    return;
}
void algogas_dt(real dtout, GridCons *grid, FluxCons *fluxes, Parameters *params) {
    real end_time = grid->time + dtout;
    real dt_max;

    while (grid->time < end_time) { 
        /* Time-step */
        dt_max = end_time - grid->time;
        algogas_single(dt_max,grid,fluxes,params);

    }
    return;
}
real calc_dt(GridCons *grid, Parameters *params) {
    int i,j,indx;
    int nx1,nx2,nx3,ntot;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;

    real dt = FLOATMAX;
    real dt_cell; 

    real *x1 = grid->xc1;
    real *x2 = grid->xc2;

    real *h1 = &grid->hfac[0*ntot ];
    real *h2 = &grid->hfac[1*ntot ];
    real *h3 = &grid->hfac[2*ntot ];

    real *rho       = &grid->cons[0*ntot];
    real *mx1       = &grid->cons[1*ntot];
    real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 

    real cs;
    real gamma = params->gamma;
    real gamma_1 = params->gamma_1; // (gamma-1)*gamma
    real dx1,dx2,dx3,dt1,dt2,dt3,u1,u2,u3;
    real pres;
    for(j=0;j<nx2;j++) {
        for(i=0;i<nx1;i++) {
            indx = INDEX(i,j); 
            pres = gamma_1*(energy[indx]-.5*(mx1[indx]*mx1[indx]
                    +mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx])/rho[indx]);
            cs = sqrt( gamma * pres/rho[indx]);
            dx1 = h1[indx]*grid->dx1[i];
            u1 = fabs(mx1[indx]/rho[indx]);
            dt1 = dx1/(u1 + cs);

            dx2 = h2[indx]*grid->dx2[j];
            u2 = fabs(mx2[indx]/rho[indx]);
            dt2 = dx2/(u2 + cs);

            dt_cell = fmin(dt1,dt2);
            dt = fmin(dt,dt_cell);
        }
    }
    return dt*(params->cfl);
}
