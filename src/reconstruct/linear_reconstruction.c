#include <math.h>

void linear_reconstruction(double dt) {

#ifdef X
    linear_reconstruction_x(dt);
#endif
#ifdef Y
    linear_reconstruction_y(dt);
#endif
#ifdef Z
    linear_reconstruction_z(dt);
#endif

}

void linear_reconstruction_x(double dt) {
    int i,j,k;
    for(k=0;k<Nx;k++) {
        for(j=0;j<Ny;j++) {
            for(i=0;i<Nx-1;i++) {
                indx = i + Nx*j + Nx*Ny*k;
                dx = xmin[i+1] - xmin[i];
                rho_L[indx] = rho[indx] + .5*(1 - dt/dx* 

                 



            }
        }
    }



}
