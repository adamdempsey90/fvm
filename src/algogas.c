#include "defs.h"

#define DTMIN 1e-8

real calc_dt(GridCons *grid, Parameters *params);
void algogas(real dtout, GridCons *grid, FluxCons *fluxes, Parameters *params) {
    real end_time = grid->time + dtout;
    real dt_max;
    real dt;

    while (grid->time < end_time) { 
        /* Time-step */
        dt_max = end_time - grid->time;
        dt = calc_dt(grid,params);
        if (dt < DTMIN){
            printf("Timestep %.4e fell below minimum value of %.1e\n",dt,DTMIN);
            exit(1);
        }
        dt = fmin(dt,dt_max);

        /* Reconstruction pt. 1 */
        /* Reconstruct and then evolve boundary values for dt/2 */
 //       reconstruct_and_evolve(grid,fluxes,dt*.5,params);
        
        /* Riemann solve with boundary values */
//        riemann_boundary(fluxes);

//#ifndef 1D
        /* Transverse flux updates */
  //      transverse_update(grid,fluxes,dt,params);
        /* Riemann solve with boundary values */
 //       riemann_boundary(fluxes);
//#endif
        /* Final update */
 //       update_cons(grid,fluxes,dt,params);


        grid->time += dt;
    }
    return;
}
real calc_dt(GridCons *grid, Parameters *params) {
    int i,j,k,indx;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1]; 
    nx3 = grid->nx[2];
    size_x1 = grid->size_x[0];
    size_x12 = grid->size_x12;

    real dt = FMAX;
    real dt_cell; 

    real *x1 = grid->xc1;
    real *x2 = grid->xc2;
    real *x3 = grid->xc3;

    real *h1 = &grid->hfac[0*grid->ntot + grid->offset];
    real *h2 = &grid->hfac[1*grid->ntot + grid->offset];
    real *h3 = &grid->hfac[2*grid->ntot + grid->offset];

    real *rho       = &grid->cons[0*grid->ntot + grid->offset];
    real *mx1       = &grid->cons[1*grid->ntot + grid->offset];
    real *mx2       = &grid->cons[2*grid->ntot + grid->offset];
    real *mx3       = &grid->cons[3*grid->ntot + grid->offset];
    real *energy    = &grid->cons[4*grid->ntot + grid->offset];
    real *intenergy = &grid->cons[5*grid->ntot + grid->offset];

    real cs;
    real gamma = params->gamma;
    real gamma_c = params->gamma_c; // (gamma-1)*gamma
    real dx1,dx2,dx3,dt1,dt2,dt3,u1,u2,u3;

    for(k=0;k<nx3;k++) {
        for(j=0;j<nx2;j++) {
            for(i=0;i<nx1;i++) {
                indx = INDEX(i,j,k); 
                cs = sqrt( gamma_c * intenergy[indx]/rho[indx]);
                dx1 = h1[indx]*grid->dx1[i];
                //dx2 = h2[indx]*(x2[j+1]-x2[j]);
                //dx3 = h3[indx]*(x3[k+1]-x3[k]);
                //printf("%.2e\t%.2e\n",cs,dx1);
                u1 = fabs(mx1[indx]/rho[indx]);
                //u2 = fabs(mx2[indx]/rho[indx]);
                //u3 = fabs(mx3[indx]/rho[indx]);
                dt1 = dx1/(u1 + cs);
                //printf("%.2e\t%.2e\t%.2e\t%.2e\n",cs,dx1,u1,dt1);
                dt_cell = dt1;
                //dt2 = dx2/(u2 + cs);
                //dt3 = dx3/(u3 + cs);
                //dt_cell = fmin(dt3,fmin(dt1,dt2));
                dt = fmin(dt,dt_cell);
            }
        }
    }
    return dt*(params->cfl);
}
