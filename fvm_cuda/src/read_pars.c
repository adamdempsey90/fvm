#include "defs.h"
void read_pars(Parameters *params, int arc, char *argv[]) {
    params->nx1 = 128;
    params->nx2 = 128;
    params->gamma = 1.4;
    params->cfl = .2;

    params->x1_min = 0.;
    params->x1_max = 1;

    params->x2_min = 0.;
    params->x2_max = 1;

    params->gamma_1 = params->gamma-1.;
    params->gamma_c = params->gamma * params->gamma_1;

    params->tend = .3;
    params->Nout = 10; //10;

    params->one_step = FALSE;

    strcpy(params->outputname ,"out/test");
    return;

}
