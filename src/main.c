#include "fvm.h"

int main(int argc , char *argv[]) {
    int n;
    input(argc,argv);
    allocate_fields();
    
    initialize();
    time = 0;
    for(n=0;n<DT;n++) {
        printf("t = %.5f (%d DT)\n", time, n);
        algogas(time,dt);
        output(time);
    }

    free_fields();
    return 1; 
}
