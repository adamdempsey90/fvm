#include "defs.h"
void read_pars(Parameters *params, int arc, char *argv[]) {
    params->nx1 = 400;
    params->nx2 = 400;
    params->gamma = 1.4;
    params->cfl = .2;

    params->x1_min = -.5;
    params->x1_max = .5;

    params->x2_min = -.5;
    params->x2_max = .5;

    params->gamma_1 = params->gamma-1.;
    params->gamma_c = params->gamma * params->gamma_1;

    params->tend = 5.;
    params->Nout = 50;
    params->dtout = (params->tend)/(float)params->Nout;

    params->one_step = FALSE;

    strcpy(params->outputname ,"out/kh");
    return;

}
