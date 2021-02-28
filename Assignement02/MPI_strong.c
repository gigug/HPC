#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <math.h>
#include <time.h>
#include "MPIlib.h"
#include <mpi.h>

#define MSGLEN 2048


int main(int argc, char *argv[]){

MPI_Init(&argc, &argv);

int m = atoi(argv[1]), n = atoi(argv[2]), kern_type = atoi(argv[3]);
double kernel[m*n];
int i_rank, ranks;
int param, symm;
MPI_Status status;
MPI_Comm_rank( MPI_COMM_WORLD, &i_rank);
MPI_Comm_size( MPI_COMM_WORLD, &ranks);
//g
int xsize, ysize, maxval;
xsize = 0;
ysize = 0;
maxval = 0;

void * ptr;

switch (kern_type){
    case 1:
    meankernel(m, n, kernel);
    break;
    case 2:
    weightkernel(m, n, param, symm, kernel);
    break;
    case 3:
    gaussiankernel(m, n, param, symm, kernel);
    break;
}

if (i_rank == 0){
    read_pgm_image(&ptr, &maxval, &xsize, &ysize, "earth-large.pgm");
}


MPI_Bcast(&xsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&ysize, 1, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(&maxval, 1, MPI_INT, 0, MPI_COMM_WORLD);

int flo, start, end, i;
flo = floor(ysize/ranks);

int first, last;

start = i_rank * flo;
end = (i_rank + 1) * flo; 

first = start - (m - 1)/2;
last = end + (m - 1)/2;

if (start == 0){
    first = 0;
}
if (end == ysize){
    last = ysize;
}

int sendcounts[ranks];
int displs[ranks];

int first2[ranks];
int last2[ranks];
int c_start2[ranks];
int c_end2[ranks];

int num;
num = (ranks - 1) * (m-1);

unsigned short int megapic[xsize*(ysize + num)];


if (i_rank == 0){
    for(i = 0; i < ranks; i++){
        c_start2[i] = i * flo;
        c_end2[i] = (i + 1) * flo; 
        if ( i == ranks - 1){
            c_end2[i] = ysize;
        }
        first2[i] = c_start2[i] - (m - 1)/2;
        last2[i] = c_end2[i] + (m - 1)/2;
        if (c_start2[i] == 0){
            first2[i] = 0;
        }
        if (c_end2[i] == ysize){
            last2[i] = ysize;
        }
        sendcounts[i] = (last2[i] - first2[i]) * xsize; 
    }

    int i, j, k, index, index_disp = 0;
    index = 0;
    displs[0] = 0;

    for (k = 0; k < ranks; k++){
        for (i = first2[k]*xsize; i < last2[k]*xsize; i++){
            megapic[index] = ((unsigned short int *)ptr)[i];
            index++;
        }
        index_disp++;
        displs[index_disp] = index;
    }

}

MPI_Bcast(displs, ranks, MPI_INT, 0, MPI_COMM_WORLD);
MPI_Bcast(sendcounts, ranks, MPI_INT, 0, MPI_COMM_WORLD);

//unsigned short int minipic[xsize*(last-first)];
void *minipic = (unsigned short int*)malloc(xsize*(last-first) * sizeof(short int));

MPI_Barrier(MPI_COMM_WORLD);
MPI_Scatterv(megapic, sendcounts, displs, MPI_UNSIGNED_SHORT, minipic, (last-first)*xsize, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

if (I_M_LITTLE_ENDIAN)
    swap_image(minipic, xsize, last-first, maxval);

void * res = elaborate(minipic, xsize, ysize, maxval, kernel, m, n, start, end, first, last);
void * final = (unsigned short int*)malloc(xsize*ysize* sizeof(short int)); // POTREBBE ESSERE NECESSARIO PRIVATIZZARE DISP E PERIOD

int sendcounts2[ranks];
int displs2[ranks];

int c_start3[ranks];
int c_end3[ranks];

if (i_rank == 0){
    for(i = 0; i < ranks; i++){
        c_start3[i] = i * flo;
        c_end3[i] = (i + 1) * flo; 
        if ( i == ranks - 1){
            c_end3[i] = ysize;
        }
        sendcounts2[i] = (c_end3[i] - c_start3[i]) * xsize;
        displs2[i] = c_start3[i] * xsize;
    }
}



MPI_Barrier(MPI_COMM_WORLD);
MPI_Gatherv(res, (end-start)*xsize, MPI_UNSIGNED_SHORT, final, sendcounts2, displs2, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

if ( I_M_LITTLE_ENDIAN )
    swap_image(final, xsize, ysize, maxval);

if (i_rank == 0){
    char *name = (char*)malloc(30);
    snprintf(name, 25, "earth-large.b_%d_%dx%d.pgm", i_rank,m,n);
    write_pgm_image(final, maxval, xsize, ysize, "qui.pgm");
}
MPI_Finalize();
