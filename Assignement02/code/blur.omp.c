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
int kern_type = atoi(argv[1]), m = atoi(argv[2]), n = atoi(argv[2]);
double param;

int symm;
if (argc < 7){
    symm = 0;
}
else{
    symm = atoi(argv[6]);
}
double *kernel;
char *name_i_com = (char*)malloc(40), *name_i, *name_o_com = (char*)malloc(40), *name_o;
if (argc < 6){
    name_i = argv[3];
    name_o = argv[4];
    param = 0;
}
else{
    param = atof(argv[3]);
    name_i = argv[4];
    name_o = argv[5];
}

if (argc > 7){
    symm = atoi(argv[6]);
    if (argc > 8){
        m = atoi(argv[7]);
    }
}
snprintf(name_i_com, 40, "%s.pgm", name_i);
if(kern_type == 1){
    double param_dec = (param - (int)param) * pow(10, strlen(argv[3])-2);
    int param_dec_int = (int)param_dec;
    snprintf(name_o_com, 40, "%s.b_%d_%dx%d_0%d.omp.pgm", name_o, kern_type, m, n, param_dec_int);
}
else{
    snprintf(name_o_com, 40, "%s.b_%d_%dx%d.omp.pgm", name_o, kern_type,m,n);
}


kernel = (double *)calloc(m * n, sizeof(double));
int xsize      = XWIDTH;
int ysize      = YWIDTH;
int maxval     = MAXVAL;

switch (kern_type)
{
    case 0:
    meankernel(m, n, kernel);
    break;
    case 1:
    weightkernel(m, n, param, symm, kernel);
    break;
    case 2:
    gaussiankernel(m, n, param, symm, kernel);
    break;
}

int layer = 1;
if (maxval > 255){
    layer = 2;
}

void * ptr = (unsigned short int*)malloc(xsize * ysize * layer);
read_pgm_image(&ptr, &maxval, &xsize, &ysize, name_i_com);

// swap the endianism
if (I_M_LITTLE_ENDIAN)
    swap_image(ptr, xsize, ysize, maxval);

void * res = (unsigned short int*)malloc(xsize * ysize * layer);

elaborate(ptr, xsize, ysize, maxval, kernel, m, n, res);

if ( I_M_LITTLE_ENDIAN )
        swap_image(res, xsize, ysize, maxval);

write_pgm_image(res, maxval, xsize, ysize, name_o_com);
free(name_o_com);
free(name_i_com);
free(ptr);
free(kernel);
free(res);
printf("The image has been written back\n");
return 0;
}
