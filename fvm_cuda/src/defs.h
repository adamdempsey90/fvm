#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "prob_defs.h"



#define TRUE 1
#define FALSE 0

//#define ISFLOAT 
#ifdef ISFLOAT
  #define real float
#else
  #define real double
#endif

#define FLOATMAX 1e8
#define PRESSUREFLOOR  1e-8
#include "structs.h"
#include "prototypes.h"


//#define INDEX(i,j,k) ( (i) + (j)*size_x1 + (k)*size_x12)
#define INDEX(i,j) ( (i) + (j)*size_x1)

//#define PROF

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif



