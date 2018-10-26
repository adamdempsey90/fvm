#include "defs.h"
void read_pars(Parameters *params, int arc, char *argv[]) {
    params->nx1 = 64;
    params->nx2 = 64;
    params->gamma = 1.6666666667;
    params->cfl = .2;

    params->x1_min = 0;
    params->x1_max = 1;

    params->x2_min = 0;
    params->x2_max = 1;

    params->gamma_1 = params->gamma-1.;
    params->gamma_c = params->gamma * params->gamma_1;

    params->tend = 4;
    params->Nout = 50; //10;
    params->dtout = (params->tend)/(float)params->Nout;

    params->one_step = TRUE;

    strcpy(params->outputname ,"out/test");
    return;

}
