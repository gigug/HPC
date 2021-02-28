#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <math.h>
#include <time.h>
#include "OMPlib.h"

#define MSGLEN 2048

int main(int argc, char *argv[]){
int m = atoi(argv[1]), n = atoi(argv[2]), kern_type = 1, temp = atoi(argv[3]), reps, rep = atoi(argv[4]);
int procs = atoi(argv[4]);
int symm = 0;
double param;
double *kernel;

kernel = (double *)calloc(m * n, sizeof(double));

if (m%2 == 0){
    printf("Please choose an odd value for the kernel height\n");
    exit(1);
}

if (n%2 == 0){
    printf("Please choose an odd value for the kernel width\n");
    exit(1);
}

if (kern_type == 2 || kern_type == 3){
    if(argc < 5){
        printf("Please choose a value for the parameter\n");
        exit(1);
    }
    param = atof(argv[4]);
}

if (kern_type == 3){
    if(argc < 6){
        printf("Please choose a symmetry option\n");
        exit(1);
    }
    symm = atoi(argv[5]);
}

int xsize      = XWIDTH;
int ysize      = YWIDTH;
int maxval     = MAXVAL;

switch (kern_type)
{
    case 1:
    meankernel(m, n, kernel);
    break;
    case 2:
    weightkernel(m, n, param, kernel);
    break;
    case 3:
    gaussiankernel(m, n, param, symm, kernel);
    break;
}
for (reps = 0; reps < rep; reps ++){
    void * ptr = generate_gradient( maxval, xsize, ysize );
    read_pgm_image(&ptr, &maxval, &xsize, &ysize, "earth-large.pgm");

    // swap the endianism
    if (I_M_LITTLE_ENDIAN)
        swap_image(ptr, xsize, ysize, maxval);

    void * res = generate_gradient( maxval, xsize, ysize );

    elaborate(ptr, xsize, ysize, maxval, kernel, m, n, res, procs);

    if ( I_M_LITTLE_ENDIAN )
            swap_image(res, xsize, ysize, maxval);

    char *name = (char*)malloc(30);
    snprintf(name, 30, "earth-large.b_%d_%dx%d.pgm", kern_type,m,n);
    write_pgm_image(res, maxval, xsize, ysize, name);

    free(ptr);
    free(kernel);
    free(res);
    printf("The image has been written back\n");
}
}
