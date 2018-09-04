#include "defs.h"




void reconstruct(GridCons *grid,FluxCons *fluxes,Parameters *params, real dt) {
    int i,j,k,indx,indxp,indxm,n;
    int ntot,nf;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
    nx2 = grid->nx[1];
    size_x1 = grid->size_x[0];
    ntot = grid->ntot;
    nf = grid->nf;
#ifdef PLM
    real dx_arr[3];
   // real temp_L[5];
   // real temp_R[5];
   // real plm_cons_m[5], plm_cons_c[5], plm_cons_p[5];
    real *plm_cons = (real *)malloc(sizeof(real)*nf*3);
#endif


    real *UL_1 = fluxes->UL_1;
    real *UR_1 = fluxes->UR_1;
    real *cons = grid->cons;
    real *dx1 = grid->dx1;

    real *UL_2 = fluxes->UL_2;
    real *UR_2 = fluxes->UR_2;
    real *dx2 = grid->dx2;

/* X1 - Direction */
    if (nx1 > 1) {
        for(j=-NGHX2; j < nx2 + NGHX2;j++) {
            for(i=-2;i<nx1+2;i++) {
    
                indxm = INDEX(i-1,j);
                indx  = INDEX(i,j); 
                indxp = INDEX(i+1,j);

#ifdef PCM
                for(n=0;n<nf;n++) {
                    UL_1[indx + n*ntot] = cons[indx + n*ntot];
                    UR_1[indx + n*ntot] = cons[indxp + n*ntot];
                }
#endif
#ifdef PLM

                dx_arr[0] = dx1[i-1];
                dx_arr[1] = dx1[i];
                dx_arr[2] = dx1[i+1];
                for(n=0;n<nf;n++) {
                    plm_cons[n]        = cons[indxm + n*ntot];
                    plm_cons[n + nf]   = cons[indx + n*ntot];
                    plm_cons[n + 2*nf] = cons[indxp + n*ntot];
                }
                //plm(&plm_cons_m[0], &plm_cons_c[0], &plm_cons_p[0],
                //        &temp_L[0],&temp_R[0],dt,dx_arr, params->gamma_1,1,nf);
                plm(plm_cons,dt,dx_arr,params->gamma_1, 1, nf);
                for(n=0;n<nf;n++) {
                    UL_1[indx + n*ntot] = plm_cons[n + 2*nf]; //temp_R[n];
                    UR_1[indxm + n*ntot] = plm_cons[n]; //temp_L[n];
                }
#endif
            }
        }
    }


/* X2 - Direction */
    if (nx2 > 1) {
        for(j=-2; j < nx2 + 2;j++) {
            for(i=-NGHX1;i<nx1+NGHX1;i++) {
    
                indxm = INDEX(i,j-1);
                indx  = INDEX(i,j); 
                indxp = INDEX(i,j+1);

#ifdef PCM
                for(n=0;n<nf;n++) {
                    UL_2[indx + n*ntot] = cons[indx + n*ntot];
                    UR_2[indx + n*ntot] = cons[indxp + n*ntot];
                }
#endif
#ifdef PLM

                dx_arr[0] = dx2[j-1];
                dx_arr[1] = dx2[j];
                dx_arr[2] = dx2[j+1];
                for(n=0;n<nf;n++) {
                    plm_cons[n]        = cons[indxm + n*ntot];
                    plm_cons[n + nf]   = cons[indx + n*ntot];
                    plm_cons[n + 2*nf] = cons[indxp + n*ntot];
                }
                plm(plm_cons,dt,dx_arr,params->gamma_1, 2, nf);
                for(n=0;n<nf;n++) {
                    UL_2[indx + n*ntot] = plm_cons[n + 2*nf]; //temp_R[n];
                    UR_2[indxm + n*ntot] = plm_cons[n]; //temp_L[n];
                }
#endif
            }
        }
    }
#ifdef PLM
    free(plm_cons);
#endif

    return;
}

