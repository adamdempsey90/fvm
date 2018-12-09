/* Plane-Parallel */
#include "defs.h"
#include <time.h>
#include "cuda_defs.h"



__device__ __managed__ real delad = .4;
__device__ __managed__ real Ftot = 1.;
__device__ __managed__ real Tbot = 1.;
__device__ __managed__ real T0 = 1.;
__device__ __managed__ real T1 = 1.;
__device__ __managed__ real Ttop = .01;
__device__ __managed__ real Pbot = 1.;
__device__ __managed__ real P0 = 1.;
__device__ __managed__ real P1 = 1.;
__device__ __managed__ real Ptop = .01;
__device__ __managed__ real xi = .01;
__device__ __managed__ real delta = .01;
__device__ __managed__ real loz = .01;
__device__ __managed__ real slope = .01;
__device__ __managed__ real nu = .01;
__device__ __managed__ real g_param = .1;
__device__ __managed__ real minF = .1;
__device__ __managed__ real ksmooth = .1;

real Tfunc(real z);
real Pfunc(real z);
__host__ __device__ static real eix(const real x);
__host__ __device__ static real e1xb(const real x);


__host__ __device__ real gravpot(real x, real y, real z) {
    return g_param*y;
}
__host__ __device__ real kinematic_viscosity(real x1, real x2, real x3) {
	return nu;
}
__host__ __device__ real heatcond_func(real dens, real x1, real x2, real x3,real delad2) {
	  real zcval = (1.-minF)/slope;


	  real result;
	  real logfac;

	  real xval = (x2 - zcval)/ksmooth;
	  real xmval = (-1. - zcval)/ksmooth;
	  real x1val = (1-x2 - zcval)/ksmooth;
	  real xm1val = (1.-(-1.) - zcval)/ksmooth;

	  logfac = (1 + exp(xval))*(1+exp(-x1val));
	  logfac /= (1+exp(xmval))*(1+exp(-xm1val));
	  result =  1. - slope*x2 + slope*ksmooth*log(logfac);

	return result*Ftot;
}
__host__ __device__ real thermal_diff(real dens, real x1, real x2, real x3,real delad2) {
    return heatcond_func(dens,x1,x2,x3,delad2)/dens/delad2;
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
    
    real Ts,Ds, Ps,Tval,Dval,Pval,delT;

    Pval = intenergy[indx_r] *(g-1);
    Dval = cons[indx_r];

    Ts = intenergy[indx] * g/cons[indx];
	Ds = cons[indx];
	Ps = Ts*delad*Ds;

    real x0 = .5*(x2[0] + x2[-1]);
    
    
    delT = (Ts - Tbot)/(x2[0]-x0);

    Tval = 2*Tbot - Pval/(delad*Dval);
//temp = Tbot + delT * (x2[j] - x0);
    printf("(%d,%d) %lg %lg %lg %lg %lg %lg\n",i,j,Pval,Dval,Pval,Ts,Ds,Ps);
    delT *= -delad/g_param;
    
    
    


    /* Velocities are reflecting */
    cons[indxg + 1*ntot] = cons[indx_r + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx_r + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx_r + 3*ntot];
    
    

    cons[indxg] = Ds* pow(Tval/Ts,1./delT-1);
    intenergy[indxg] = Tval*cons[indxg]/g;
    
    cons[indxg + 4*ntot] = intenergy[indxg]  + .5*(cons[indxg+1*ntot]*cons[indxg+1*ntot]
    +cons[indxg + 2*ntot]*cons[indxg + 2*ntot]
    +cons[indxg + 3*ntot]*cons[indxg + 3*ntot])/cons[indxg];
    

    /* outflow for scalars */
    for(n=5;n<nf;n++) {
        cons[indxg + n*ntot] = cons[indx + n*ntot];
    }


    return;
}
__device__ void fixed_flux_upper(int indxg, int i, int j, int k, real *cons, real *intenergy, real *x1, real *x2, real *x3,
		int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    /* Fixed heat flux with hydrostaic balance 
     * Constant conductivity K.
     * 
     * F = -K dT/dy  
     * dP/dy = - d g
     * rho*e = Cp*T*rho/gamma
     * delrad = F*delad/(K*g)
     * T = T0 - F/K (y-y0)
     * P/P0 = (T/T0)^(1./delrad)
     * d/d0 = (T/T0)^(1./delrad-1)
     */
    int n;
    int indx_r = GINDEX(i,nx2+ -(j-nx2)-1,k);
    int indx = GINDEX(i,nx2-1,k);
    real Tval,Pval;
    real Ts,Ds,Ps,x20;
    
    Ts = intenergy[indx] * g/cons[indx];
	Ds = cons[indx];
	Ps = Ds*Ts*delad;
	x20 = 1-x2[nx2-1];


	Tval = Ts - log((1-slope*(1-x2[j]))/(1-slope*x20))/slope; // Exact conductive T(z)
	Pval = Ps*exp( (1-slope*x20)/delad*exp(slope*Ts)*(eix(-slope*Tval)-eix(-slope*Ts))); // Exact hydrostatic pressure for T(z)
	cons[indxg] = Pval/(delad*Tval);

    /* Velocities are reflecting */
    cons[indxg + 1*ntot] = cons[indx_r + 1*ntot];
    cons[indxg + 2*ntot] = -cons[indx_r + 2*ntot];
    cons[indxg + 3*ntot] = cons[indx_r + 3*ntot];
    
    

    intenergy[indxg] = Pval/(g-1);
    

    
    
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
	//fixed_temp_lower(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
	return;
}
__device__ void x2_boundary_outer(int indxg, int i, int j,int k, real *cons, real *intenergy, real *x1, real *x2, real *x3, int nx1, int nx2, int nx3, int ntot, int nf, int size_x1, int size_x12, int offset, real g, real time) {
    fixed_flux_upper(indxg,i,j,k,cons,intenergy,x1,x2,x3,nx1,nx2,nx3,ntot,nf,size_x1,size_x12,offset,g,time);
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


    loz = params->loz;
    slope = params->slope;
    delta = params->delta;
    Ftot = params->ftot;
    xi = params->xi;
    nu = params->nu;
    g_param = params->g;
    ksmooth = params->ksmooth;
    minF = params->minf;

	real *x2 = grid->xc2;
	real *x2m = grid->xm2;


	real *rho       = &grid->cons[0*ntot];
	real *mx1       = &grid->cons[1*ntot];
	real *mx2       = &grid->cons[2*ntot];
    real *mx3       = &grid->cons[3*ntot];
    real *energy    = &grid->cons[4*ntot];
    real *intenergy = grid->intenergy; 
    real gamma = params->gamma;
    real ke;
    real temp,pres;

    real u1 = 0;
    real u2 = 0;
    real u3 = 0;
    //real pres,P0;

    delad = 1.-1./gamma;
	Ttop = xi/delad;
	T1 = Ttop + log((1+slope)/(1+slope*loz))/slope;
	T0 = T1 + (1. + delta/delad)*(1+2*loz);

	Ptop = xi;
	P1 = Ptop *exp( (1+slope)/delad *exp(slope*Ttop)*(eix(-slope*T1)-eix(-slope*Ttop)));
	P0 = P1 *pow(T0/T1,1./(delad+delta));


	Tbot = Tfunc(-1.0);
	Pbot = Pfunc(-1.0);

	printf("Delad:%lg\nDelta:%lg\nTtop:%lg\nPtop:%lg\nT1:%lg\nT0:%lg\nP1:%lg\nP0:%lg\n",delad,delta,Ttop,Ptop,T1,T0,P1,P0);

	real norm;
	srand(time(NULL));

    for(k=-NGHX3;k<nx3+NGHX3;k++) {
		for(j=-NGHX2;j<nx2+NGHX2;j++) {
			for(i=-NGHX1;i<nx1+NGHX1;i++) {
				indx = INDEX(i,j,k);

				temp = Tfunc(x2[j]);
				pres = Pfunc(x2[j]);
				rho[indx] =  pres/(delad*temp);
				//if (j < 0) printf("IC\t(%d,%d)\t%lg\t%lg\t%lg\n",i,j,temp,pres,rho[indx]);

				if ((x2[j]>0)&&(x2[j]<1)) {
					norm =(real)((double)rand() / (double)RAND_MAX );
					u1 += (norm-.5)*.001;
					norm =(real)((double)rand() / (double)RAND_MAX );
					u2 += (norm-.5)*.001;
				}
				else {
					u1 = 0;
					u2 = 0;
				}
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
//    FILE *f = fopen("out/kappa.dat","w");
//    for(j=-NGHX2;j<nx2+NGHX2;j++) fprintf(f,"%.8e\t%.8e\t%.8e\t%.8e\n",x2m[j],heatcond_func(1.,0.,x2m[j],0.,delad),Tfunc(x2m[j]),Pfunc(x2m[j]));
//    fclose(f);


    return;


}

__host__ __device__ static real e1xb(const real x) {
  int maxit;
  real eps,euler,fpmin,fpmax;
  eps=6e-8;
  euler=.57721566;
  maxit=20;
  fpmin=1e-30;
  fpmax=1e30;
  int k;
  real sum,term;

  if (x == 0.0) return fpmax;
  if (x <= 1.) {
    sum = 1.;
    term = 1.;
    k=1;
    while ((fabs(term)>fabs(sum)*eps)&&(k<=maxit)){
      term = -term*k*x/( (k+1.)*(k+1.));
      sum = sum + term;
      k++;
    }
    return -euler - log(x) + x*sum;
  }

  term = 0.;
  for (k=maxit+(int)floor(80/x);k>=1;k--) {
    term = k/(1. + k/(x+term));
  }
  return exp(-x) / (x+term);
}

__host__ __device__ static real eix(const real x) {
  int maxit;
  real eps,euler,fpmin,fpmax;
  eps=6e-8;
  euler=.57721566;
  maxit=100;
  fpmin=1e-30;
  fpmax=1e30;
  int k;
  real sum,term;
  if (x == 0.0) return -fpmax;

  if (x < 0.0) return -e1xb(-x);
  if (fabs(x) <= 40.) {
    // Power series around x=0

    sum = 1.0;
    term =1.0;
    k = 1;
    while ((fabs(term/sum) > eps)&&(k<=maxit)) {
      term = term*k*x/( (k+1.)*(k+1.));
      sum += term;
      k++;
    }
    return euler + log(x) + x*sum;
  }

  // Asymptotic expansion (the series is not convergent)
  sum = 1.0;
  term = 1.0;

  for (k=1;k<=20;k++) {
    term = term*k/x;
    sum += term;
  }
  return exp(x)/x * sum;
}


real Tfunc(real z) {

    if (z < -loz) {
        return T0 + log((1-slope*z)/(1+slope*loz))/slope;
    }

    if (z < 1+loz) {
        return T1 - (1. + delta/delad)*(z -(1+loz));
    }

    return Ttop - log((1-slope*(1-z))/(1+slope))/slope;
}
real Pfunc(real z) {
	real Tval = Tfunc(z);
    if (z < -loz) {
        return P0*exp((1+slope*loz)/delad *exp(-slope*T0)*(eix(slope*Tval)-eix(slope*T0)));
    }

    if (z < 1+loz) {
        return P1*pow(Tval/T1,1./(delad+delta));
    }

    return Ptop *exp((1+slope)/delad *exp(slope*Ttop)*(eix(-slope*Tval)-eix(-slope*Ttop)));
}


