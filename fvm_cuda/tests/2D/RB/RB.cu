/* Plane-Parallel */
#include "defs.h"
#include <time.h>
#include "cuda_defs.h"




__device__ __managed__ real cond = .01;
__device__ __managed__ real nu = .01;
__device__ __managed__ real Tlower = 1000.;
__device__ __managed__ real Tupper = 100.;
__device__ __managed__ real g_param = .1;



__device__ real gravpot(real x, real y, real z) {
    return g_param*y;
}

__device__ real heatcond_func(real dens, real x1, real x2, real x3,real delad) {
	real y = 1./(1 + exp(-(x2-1.4)/.01)) + 1./(1 + exp((x2-.6)/.01));

    return cond * exp(3*y);
}
__device__ real thermal_diff(real dens, real x1, real x2, real x3,real delad) {
    return heatcond_func(dens,x1,x2,x3,delad)/dens/delad;
}
__device__ real kinematic_viscosity(real x1, real x2, real x3) {
	return nu;
}

__device__ void fixed_temp_lower(int indxg, int i, int j, int k, real *cons, real *intenergy, real *x1, real *x2, real *x3,
		int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    /* Fixed temp with hydrostatic balance
     * 
     * dP/dy = - d g
     * rho*e = Cp*T*rho/gamma
     * delT = -dT/dy *delad/g
     * T = T0 +dT/dy (y-y0)
     * P/P0 = (T/T0)^(1./delT)
     * d/d0 = (T/T0)^(1./delT-1)
     */
    int n;
    int indx_r = GINDEX(i,-j-1,k);
    int indx = GINDEX(i,0,k);
    
    real T0,P0,d0,delad, delT,temp;
    delad = 1 - 1./g;

    T0 = intenergy[indx] * g/cons[indx];
    d0 = cons[indx];
    P0 = d0*T0*delad;
    real x0 = .5*(x2[0] + x2[-1]);
    
    
    delT = (T0 - Tlower)/(x2[0]-x0);
    temp = Tlower + delT * (x2[j] - x0);
  //  printf("%lg %lg %lg %lg %lg %lg\n",T0,Tlower,delT,x2[j],x0,temp);
    delT *= -delad/g_param;
    
    

    /* Velocities are wall & no slip */
    cons[indxg + 1*ntot] = -cons[indx_r + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx_r + 2*ntot];
    cons[indxg + 3*ntot] = -cons[indx_r + 3*ntot];
    
    

    cons[indxg] = d0* pow(temp/T0,1./delT-1);
    intenergy[indxg] = temp*cons[indxg]/g;
    
    cons[indxg + 4*ntot] = intenergy[indxg]  + .5*(cons[indxg+1*ntot]*cons[indxg+1*ntot]
    +cons[indxg + 2*ntot]*cons[indxg + 2*ntot]
    +cons[indxg + 3*ntot]*cons[indxg + 3*ntot])/cons[indxg];
    

    /* outflow for scalars */
    for(n=5;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}

__device__ void fixed_temp_upper(int indxg, int i, int j, int k, real *cons, real *intenergy, real *x1, real *x2, real *x3,
		int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    /* Fixed temp with hydrostatic balance
     *
     * dP/dy = - d g
     * rho*e = Cp*T*rho/gamma
     * delT = -dT/dy *delad/g
     * T = T0 +dT/dy (y-y0)
     * P/P0 = (T/T0)^(1./delT)
     * d/d0 = (T/T0)^(1./delT-1)
     */
    int n;
    int indx_r = GINDEX(i,nx2 + -(j-nx2)-1,k);
    int indx = GINDEX(i,nx2-1,k);

    real T0,P0,d0,delad, delT,temp;
    delad = 1 - 1./g;

    T0 = intenergy[indx] * g/cons[indx];
    d0 = cons[indx];
    P0 = d0*T0*delad;
    real x0 = .5*(x2[nx2] + x2[nx2-1]);


    delT = (Tupper - T0)/(x0-x2[nx2-1]);
    temp = Tupper + delT * (x2[j] - x0);
  //  printf("%lg %lg %lg %lg %lg %lg\n",T0,Tlower,delT,x2[j],x0,temp);
    delT *= -delad/g_param;
    


    /* Velocities are wall & no slip */
    cons[indxg + 1*ntot] = -cons[indx_r + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx_r + 2*ntot];
    cons[indxg + 3*ntot] = -cons[indx_r + 3*ntot];
    
    

    cons[indxg] = d0* pow(temp/T0,1./delT-1);
    intenergy[indxg] = temp*cons[indxg]/g;
    
    cons[indxg + 4*ntot] = intenergy[indxg]  + .5*(cons[indxg+1*ntot]*cons[indxg+1*ntot]
    +cons[indxg + 2*ntot]*cons[indxg + 2*ntot]
    +cons[indxg + 3*ntot]*cons[indxg + 3*ntot])/cons[indxg];
    

    /* outflow for scalars */
    for(n=5;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}

__device__ void x2_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	fixed_temp_lower(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}
__device__ void x2_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    fixed_temp_upper(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}
__device__ void x1_boundary_inner(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	periodic_boundary_inner(1,indxg,i,j,k,cons,intenergy,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}
__device__ void x1_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
	periodic_boundary_outer(1,indxg,i,j,k,cons,intenergy,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}


void scale_factors(real x1, real x2, real x3, real *h1, real *h2, real *h3) {
	*h1 = 1.;
	*h2 = 1.;
	*h3 = 1.;
}
void init_mesh(GridCons *grid, Parameters *params) {

	init_uniform_mesh(grid,params);

	return;

}
void init_gas(GridCons *grid, Parameters *params) {
    int i,j,k,indx;
    int nx1,nx2,nx3,n,ntot,nf;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    nx3 = grid->nx[2];
    size_x1 = grid->size_x1;
    size_x12 = grid->size_x12;
    ntot = grid->ntot;
    nf = grid->nf;


    Tlower = params->tlower;
    Tupper = params->tupper;
    cond = params->kappa;
    nu = params->nu;
    g_param = params->g;

	real *x2 = grid->xc2;
	real *xm2 = grid->xm2;


	real *rho       = &grid->cons[0*ntot];
	real *mx1       = &grid->cons[1*ntot];
	real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 
    real gamma = params->gamma;
    real ke;

    real u1 = 0;
    real u2 = 0;
    real u3 = 0;
    //real pres,P0;

    real delad = 1. - 1./gamma;
    //real delrad = F*delad/(cond*g_param);
    real D0,temp;
    D0 = 1.;
    real delT = (Tupper - Tlower)/(xm2[nx2]-xm2[0]);

    real delrho = delT*-delad/g_param;
    for(k=-NGHX3;k<nx3+NGHX3;k++) {
		for(j=-NGHX2;j<nx2+NGHX2;j++) {
			for(i=-NGHX1;i<nx1+NGHX1;i++) {
				indx = INDEX(i,j,k);

            

				temp = Tlower + delT *(x2[j]-xm2[0]);


				//pres = P0 *pow(temp/Tlower,1./delT*-g_param/delad);
				rho[indx] = D0*pow(temp/Tlower,1./delrho - 1.);


				mx1[indx] = u1*rho[indx];
				mx2[indx] = u2*rho[indx];
				mx3[indx] = u3*rho[indx];

				ke = mx1[indx]*mx1[indx] + mx2[indx]*mx2[indx] + mx3[indx]*mx3[indx];
				ke /= 2.*rho[indx];
				intenergy[indx] = temp *rho[indx]/ gamma;
				energy[indx] = intenergy[indx] + ke;
				for(n=5;n<nf;n++) {
					grid->cons[n*ntot+indx] = 0;
				}




		}
    }
}


return;


}
//}
