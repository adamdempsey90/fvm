#include <math.h>

extern void boundary(void);
extern void source_step(double);
extern double cfl(void);
extern double linear_reconstruction(double);
extern void riemann(double);
extern void viscosity(double);
extern void update(double);

void algogas(double dt_max) {

    boundary();
    double dt = fmin(cfl(),dt_max);
    ctu(dt);

    return;
}
