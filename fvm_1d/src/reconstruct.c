#include "defs.h"




void reconstruct(GridCons *grid,FluxCons *fluxes,Parameters *params, real dt) {
    int i,j,k,indx,indxp,indxm,n;
    int ntot,nf;
    int nx1,nx2,nx3;
    int size_x1,size_x12; 
    nx1 = grid->nx[0];
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


    real *UL = fluxes->UL_1;
    real *UR = fluxes->UR_1;
    real *cons = grid->cons;
    real *dx1 = grid->dx1;


/* X1 - Direction */
    if (nx1 > 1) {
        for(i=-1;i<nx1+1;i++) {
    
            indxm = INDEX(i-1);
            indx  = INDEX(i); 
            indxp = INDEX(i+1);

#ifdef PCM
            for(n=0;n<nf;n++) {
                UL[indx + n*ntot] = cons[indx + n*ntot];
                UR[indx + n*ntot] = cons[indxp + n*ntot];
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
                UL[indx + n*ntot] = plm_cons[n + 2*nf]; //temp_R[n];
                UR[indxm + n*ntot] = plm_cons[n]; //temp_L[n];
            }
#endif
        }
    }

#ifdef PLM
    free(plm_cons);
#endif

    return;
}

