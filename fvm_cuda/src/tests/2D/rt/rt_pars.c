#include "defs.h"
void read_pars(Parameters *params, int arc, char *argv[]) {
    params->nx1 = 100;
    params->nx2 = 400;
    params->gamma = 1.4;
    params->cfl = .2;

    params->x1_min = -.5;
    params->x1_max = .5;

    params->x2_min = -2;
    params->x2_max = 2;

    params->gamma_1 = params->gamma-1.;
    params->gamma_c = params->gamma * params->gamma_1;

    params->tend = 100;
    params->Nout = 100;
    params->dtout = (params->tend)/(float)params->Nout;

    params->one_step = FALSE;

    strcpy(params->outputname ,"out/test");
    return;

}
