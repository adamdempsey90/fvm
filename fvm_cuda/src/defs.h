#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define TRUE 1
#define FALSE 0
#define real double 
//#define ISFLOAT 
#define FLOATMAX 1e8
#define PRESSUREFLOOR  1e-8
#include "structs.h"
#include "prototypes.h"

#define NGHX1 3
#define NGHX2 3

//#define INDEX(i,j,k) ( (i) + (j)*size_x1 + (k)*size_x12)
#define INDEX(i,j) ( (i) + (j)*size_x1)


#define CTU
//#define PCM
#define PLM
//#define MINMOD

//#define POTENTIAL

#define EXACT
#define EXACT_TOL 1e-6

//#define HLL
//#define HLLC

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif



