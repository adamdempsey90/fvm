#include "defs.h"
void read_pars(Parameters *params) {
    params->nx1 = 100;
    params->gamma = 1.4;
    params->cfl = .2;

    params->x1_min = 0.;
    params->x1_max = 1.;


    params->gamma_1 = params->gamma-1.;
    params->gamma_c = params->gamma * params->gamma_1;

    params->tend = .25;
    params->Nout = 20;
    params->dtout = (params->tend)/(float)params->Nout;

    params->one_step = TRUE;

    strcpy(params->outputname ,"out/sod");
    return;

}
