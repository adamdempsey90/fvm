#include "defs.h"
#include "cuda_defs.h"

__device__ void hllc(real dL, real uL, real pL, real aL,
        real dR, real uR, real pR, real aR,
        real *SL, real *SR, real *Sstar,
        real g, real g1, real g2, real g3, real g4, real g5) {

    real ps, us;
    anrs(dL,uL,pL,aL,
            dR,uR,pR,aR,
            &ps,&us,
            g,g1,g2,g3,g4);

    *SL = (ps<=pL) ? uL-aL : uL-aL*sqrt(1 + g5*(ps/pL-1));
    *SR = (ps<=pR) ? uR+aR : uR+aR*sqrt(1 + g5*(ps/pR-1));
    *Sstar = (pR-pL) + dL*uL*(*SL-uL) - dR*uR*(*SR-uR);
    *Sstar /= dL*(*SL-uL) - dR*(*SR-uR);
    
    return;
}
