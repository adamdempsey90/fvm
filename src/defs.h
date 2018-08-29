#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define TRUE 1
#define FALSE 0
#define real float
#define ISFLOAT 
#define FMAX 1e8

#include "structs.h"
#include "prototypes.h"

#define NGHX1 3
#define NGHX2 0
#define NGHX3 0

#define INDEX(i,j,k) (i + j*size_x1 + k*size_x12)


