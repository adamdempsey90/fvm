#include "defs.h"



void riemann_fluxes(const real *UL, const real *UR, real *flux, int dir1,int nx[3],int size_x1, int size_x12,int nf,int ntot, real gamma) {
    int i,j,k,indx,indxp;
    int nx1,nx2,nx3,n;
    nx1 = nx[0];
    nx2 = nx[1];
    real rho,mx1,mx2,mx3,energy;

    int il,iu,jl,ju;

    int dir2, dir3; 

    /* 1->2->3 
     * 2->3->1
     * 3->1->2
     */
    dir2 = (dir1)%3 + 1;
    dir3 = (dir2)%3 + 1;

    if (dir1 == 1) {
        il = -1; iu = nx1+1;
        jl = -NGHX2;ju = nx2+NGHX2;
    }
    else if (dir1 == 2) {
        il = -NGHX1; iu = nx1+NGHX1;
        jl = -1; ju = nx2+1;
    }
    else {
        printf("Direction can only be 1 or 2 in 2D!\n");
        exit(1);
    }

	real *UL_cell, *UR_cell, *F;
#pragma omp threadprivate(UL_cell,UR_cell,F)
#pragma omp parallel num_threads(nth)
{
    UL_cell = (real *)malloc(sizeof(real)*nf);
    UR_cell = (real *)malloc(sizeof(real)*nf);
    F = (real *)malloc(sizeof(real)*nf);
}

#pragma omp parallel for num_threads(nth) private(indx,indxm2,indxp2,indxp1m2,dtdx2)           
    for(j=jl;j<ju;j++) {
        for(i=il;i<iu;i++) {
            indx = INDEX(i,j); 
            UL_cell[0] = UL[0*ntot + indx];
            UL_cell[1] = UL[dir1*ntot + indx];
            UL_cell[2] = UL[dir2*ntot + indx];
            UL_cell[3] = UL[dir3*ntot + indx];
            for(n=4;n<nf;n++) {
                UL_cell[n] = UL[n*ntot + indx];
            }
            UR_cell[0] = UR[0*ntot + indx];
            UR_cell[1] = UR[dir1*ntot + indx];
            UR_cell[2] = UR[dir2*ntot + indx];
            UR_cell[3] = UR[dir3*ntot + indx];
            for(n=4;n<nf;n++) {
                UR_cell[n] = UR[n*ntot + indx];
            }
#ifdef EXACT_RIEMANN
            exact_flux(UL_cell,UR_cell,F,gamma,gamma-1,nf);
#endif
#if defined(HLLC)||defined(HLL)
            hll_flux(UL_cell,UR_cell,F,gamma,gamma-1, (gamma+1)/(2*gamma),nf);
#endif
            flux[0*ntot + indx] = F[0];
            flux[dir1*ntot + indx] = F[1];
            flux[dir2*ntot + indx] = F[2];
            flux[dir3*ntot + indx] = F[3];
            flux[4*ntot + indx] = F[4];
            for(n=5;n<nf;n++) {
                flux[n*ntot + indx] = F[n];
            }
        }
    }
    

#pragma omp parallel num_threads(nth)
{
    free(UL_cell);free(UR_cell);free(F);
}
    return;

}


