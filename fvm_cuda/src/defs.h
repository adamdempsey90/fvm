#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "prob.h"





#define TRUE 1
#define FALSE 0

#ifndef real
#define real double
#endif


#ifndef EXACT_TOL
#define EXACT_TOL 1e-6
#endif

#ifndef FLOATMAX
#define FLOATMAX 1e8
#endif

#ifndef PRESSUREFLOOR
#define PRESSUREFLOOR  1e-8
#endif


#ifdef PPM
#define NGHX1 4
#else
#define NGHX1 3
#endif



#ifdef DIMS3
	#define NGHX3 NGHX1
#ifndef DIMS2
	#define DIMS2
#endif
#else
	#define NGHX3 0
#endif


#ifdef DIMS2
	#define NGHX2 NGHX1
#ifndef DIMS1
	#define DIMS1
#endif
#else
	#define NGHX2 NGHX1
#endif



#include "structs.h"
#include "prototypes.h"


#define INDEX(i,j,k) ( (i) + (j)*size_x1 + (k)*size_x12)


#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif



