#include "defs.h"



void riemann_fluxes(const real *UL, const real *UR, real *flux, int dir1,int nx[3],int size_x1, int size_x12,int nf,int ntot, real gamma) {
    int i,j,k,indx,indxp;
    int nx1,nx2,nx3,n;
    nx1 = nx[0];
    real rho,mx1,mx2,mx3,energy;

    int dir2, dir3; 

    /* 1->2->3 
     * 2->3->1
     * 3->1->2
     */
    //dir2 = (dir1)%3 + 1;
    //dir3 = (dir1+1)%3 + 1;
      dir2 = 2;
      dir3 = 3;
  //  printf("%d\n",nf);
  //  printf("%d\t%d\t%d\n",nx1,size_x1,size_x12);


  //  for(i=-NGHX1;i<nx1+NGHX1;i++) printf("%d\t%lg\t%lg\n",i,UL[4*ntot+i],UR[4*ntot+i]);
    real *UL_cell = (real *)malloc(sizeof(real)*nf);
    real *UR_cell = (real *)malloc(sizeof(real)*nf);
    real *F = (real *)malloc(sizeof(real)*nf);

    for(i=-1;i<nx1+NGHX1;i++) {
        //printf("kj %d %d \n",k,j);
        indx = INDEX(i) ;
        //printf("UL_Cell %d\n",indx);
        //for(n=0;n<nf;n++) printf("(%d,%lg)\t",n,UL[n*ntot+indx]);
        //printf("\n");
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
    

    free(UL_cell);free(UR_cell);free(F);

    return;

}


