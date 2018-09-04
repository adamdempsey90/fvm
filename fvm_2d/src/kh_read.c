#include "defs.h"
void read_pars(Parameters *params, int arc, char *argv[]) {
    params->nx1 = 320;
    params->nx2 = 180;
    params->gamma = 1.4;
    params->cfl = .2;

    params->x1_min = 0.;
    params->x1_max = 1.;

    params->x2_min = 0.;
    params->x2_max = .5625;

    params->gamma_1 = params->gamma-1.;
    params->gamma_c = params->gamma * params->gamma_1;

    params->tend = 1.;
    params->Nout = 20;
    params->dtout = (params->tend)/(float)params->Nout;

    params->one_step = FALSE;

    strcpy(params->outputname ,"out/kh");
    return;

}
