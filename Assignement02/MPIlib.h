#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define XWIDTH 256
#define YWIDTH 256
#define MAXVAL 65535


#if ((0x100 & 0xf) == 0x0)
#define I_M_LITTLE_ENDIAN 1
#define swap(mem) (( (mem) & (short int)0xff00) >> 8) +	\
  ( ((mem) & (short int)0x00ff) << 8)
#else
#define I_M_LITTLE_ENDIAN 0
#define swap(mem) (mem)
#endif


// =============================================================
//  utilities for managinf pgm files
//
//  * write_pgm_image
//  * read_pgm_image
//  * swap_image
//
// =============================================================
void center(int m, int n, int * x_center, int * y_center);
void meankernel(int m, int n, double * kernel){
  double det = 1.0/(m*n);
  int i, j;
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      kernel[i + m * j] = det;
    }
  }
}

void weightkernel(int m, int n, double param, int symm, double * kernel){
  // per asimmetria chiedi da quale parte mettere i valori negativi
  double w;
  int i, j, x_center, y_center;
  if (symm == 1){
    center(m, n, &x_center, &y_center);
  }
  else{
    x_center = (m-1)/2;
    y_center = (n-1)/2;
  }
  w = (1-param)/(m * n - 1);
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      kernel[i + m * j] = w;
    }
  }
  //kernel[(int)floor((m-1)/2.0) + m * (int)floor((n-1)/2.0)] = param;
  kernel[x_center + m * y_center] = param;
}

void gaussiankernel(int m, int n, double param, int symm, double * kernel){
  double sum = 0;
  int i, j, x_center, y_center;
  double den;
  if (symm == 1){
    center(m, n, &x_center, &y_center);
  }
  else{
    x_center = (m-1)/2;
    y_center = (n-1)/2;
  }
  if (param == 0){
    kernel[x_center + m * y_center] = 1;
    return;
  }
  den = 1/(sqrt(2 * M_PI)*param);
  param = param * param;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      double x = i - x_center;
      double y = j - y_center;
      kernel[i + m * j] = den * exp(-(x * x + y * y) / (2 * param));
      sum += kernel[i + m * j];
    }
  }
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            kernel[i + m * j] /= sum;
        }
    }
}

void center(int m, int n, int * x_center, int * y_center){
  printf("Enter an x value from 1 to %d \n", m);
  while(*x_center < 1 || *x_center > m){
    scanf("%d", x_center);
  }
  *x_center -= 1;
  printf("Enter a y value from 1 to %d \n", n);
  while(*y_center < 1 || *y_center > n){
    scanf("%d", y_center);
  }
  *y_center -= 1;
  return;
}


void * elaborate(void * ptr, int xsize, int ysize, int maxval, double * kernel, int m, int n, int start, int end, int first, int last)
{

short int *sImage;   // the image when a two bytes are used for each pixel
void * res;
sImage = (unsigned short int*)calloc(xsize*(end-start), sizeof(short int));
unsigned short int _maxval = swap((unsigned short int)maxval);
int idx = 0;



unsigned short int ind_i = (m - 1)/2;
unsigned short int ind_j = (n - 1)/2;

double count = 0;
double val = 0;

unsigned short int val_u, val_s;


int i, j, ii, jj, iii, jjj;

// scorro le righe
iii = 0;
jjj = 0;

for (i = 0 + (start - first); i < end - first; i++){

  //scorro le colonne
  
  for (j = 0; j < xsize; j++){

    // scorro le righe (kernel)
    
    for (ii = -ind_i; ii < ind_i + 1; ii++){

      // scorro le colonne (kernel)
      
      for (jj = -ind_j; jj < ind_j + 1; jj++){

        if( (i + ii) < 0 || (i + ii) >= last - first || (j + jj) < 0 || (j + jj) >= xsize ){
          count += kernel[(jj + ind_j) + n * (ii + ind_i)];
        }
        else{
          val += ((unsigned short int *)ptr)[(j + jj) + xsize * (i + ii)] * kernel[(jj + ind_j) + n * (ii + ind_i)];
          //printf("%d\n", ((unsigned short int *)ptr)[(j + jj) + xsize * (i + ii)]);
        }

      }

    }
    if (count != 0 && count != 1){
      val /= (1-count);
    }
    val_u = (short int)val; // CONTROLLA SE FUNZIONA
    val_s = swap(val_u);
    sImage[jjj + xsize * iii] = val_u;
    //((unsigned short int *)res)[jjj + xsize * iii] = val_s;
    val = 0;
    count = 0;
    jjj++;
    //printf("iii is %d, jjj is %d\n", iii, jjj);

  }
iii++;
jjj = 0;
//printf("Index %d, start is %d, first is %d, end is %d, last is %d, value is %d\n", i, start, first, end, last, sImage[i]);
} 
res = (void*)sImage;
return res;
}



void write_pgm_image( void *image, int maxval, int xsize, int ysize, const char *image_name)
/*
 * image        : a pointer to the memory region that contains the image
 * maxval       : either 255 or 65536
 * xsize, ysize : x and y dimensions of the image
 * image_name   : the name of the file to be written
 *
 */
{
  FILE* image_file; 
  image_file = fopen(image_name, "w"); 
  
  // Writing header
  // The header's format is as follows, all in ASCII.
  // "whitespace" is either a blank or a TAB or a CF or a LF
  // - The Magic Number (see below the magic numbers)
  // - the image's width
  // - the height
  // - a white space
  // - the image's height
  // - a whitespace
  // - the maximum color value, which must be between 0 and 65535
  //
  // if he maximum color value is in the range [0-255], then
  // a pixel will be expressed by a single byte; if the maximum is
  // larger than 255, then 2 bytes will be needed for each pixel
  //

  int color_depth = 1 + ( maxval > 255 );

  fprintf(image_file, "P5\n# generated by\n# put here your name\n%d %d\n%d\n", xsize, ysize, maxval);
  
  // Writing file
  fwrite( image, 1, xsize*ysize*color_depth, image_file);  

  fclose(image_file); 
  return ;

  /* ---------------------------------------------------------------

     TYPE    MAGIC NUM     EXTENSION   COLOR RANGE
           ASCII  BINARY

     PBM   P1     P4       .pbm        [0-1]
     PGM   P2     P5       .pgm        [0-255]
     PPM   P3     P6       .ppm        [0-2^16[
  
  ------------------------------------------------------------------ */
}


void read_pgm_image( void **image, int *maxval, int *xsize, int *ysize, const char *image_name)
/*
 * image        : a pointer to the pointer that will contain the image
 * maxval       : a pointer to the int that will store the maximum intensity in the image
 * xsize, ysize : pointers to the x and y sizes
 * image_name   : the name of the file to be read
 *
 */
{
  FILE* image_file; 
  image_file = fopen(image_name, "r"); 

  *image = NULL;
  *xsize = *ysize = *maxval = 0;
  
  char    MagicN[2];
  char   *line = NULL;
  size_t  k, n = 0;
  
  // get the Magic Number
  k = fscanf(image_file, "%2s%*c", MagicN );

  // skip all the comments
  k = getline( &line, &n, image_file);
  while ( (k > 0) && (line[0]=='#') ){
    k = getline( &line, &n, image_file);
  }

  if (k > 0)
    {
      k = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);
      if ( k < 3 )
    fscanf(image_file, "%d%*c", maxval);
    }
  else
    {
      *maxval = -1;         // this is the signal that there was an I/O error
                // while reading the image header
      free( line );
      return;
    }
  free( line );
  
  int color_depth = 1 + ( *maxval > 255 );
  unsigned int size = *xsize * *ysize * color_depth;
  
  if ( (*image = (char*)malloc( size )) == NULL )
    {
      fclose(image_file);
      *maxval = -2;         // this is the signal that memory was insufficient
      *xsize  = 0;
      *ysize  = 0;
      return;
    }
  
  if ( fread( *image, 1, size, image_file) != size )
    {
      free( image );
      image   = NULL;
      *maxval = -3;         // this is the signal that there was an i/o error
      *xsize  = 0;
      *ysize  = 0;
    }  

  fclose(image_file);
  return;
}


void swap_image( void *image, int xsize, int ysize, int maxval )
/*
 * This routine swaps the endianism of the memory area pointed
 * to by ptr, by blocks of 2 bytes
 *
 */
{
  if ( maxval > 255 )
    {
      // pgm files has the short int written in
      // big endian;
      // here we swap the content of the image from
      // one to another
      //
      int i;
      unsigned int size = xsize * ysize;
      for (i = 0; i < size; i++){
        ((unsigned short int*)image)[i] = swap(((unsigned short int*)image)[i]);
      }
    }
  return;
}



// =============================================================
//

void * generate_gradient( int maxval, int xsize, int ysize )
/*
 * just and example about how to generate a vertical gradient
 * maxval is either 255 or 65536, xsize and ysize are the
 * x and y dimensions of the image to be generated.
 * The memory region that will contain the image is returned
 * by the function as a void *

 */
{
  char      *cImage;   // the image when a single byte is used for each pixel
  short int *sImage;   // the image when a two bytes are used for each pixel
  void      *ptr;
  
  int minval      = 0; 
  int delta       = (maxval - minval) / ysize;
  int yy, xx;
  if(delta < 1)
    delta = 1;
  
  if( maxval < 256 )
    // generate a gradient with 1 byte of color depth
    {
      cImage = (char*)calloc( xsize*ysize, sizeof(char) );
      unsigned char _maxval = (char)maxval;
      int idx = 0;
      for (yy = 0; yy < ysize; yy++ )
    {
      unsigned char value = minval + yy*delta;
      for(xx = 0; xx < xsize; xx++ )
        cImage[idx++] = (value > _maxval)?_maxval:value;
    }
      ptr = (void*)cImage;
    }
  else
    // generate a gradient with 2 bytes of color depth
    {
      sImage = (unsigned short int*)calloc(xsize*ysize, sizeof(short int));
      unsigned short int _maxval = swap((unsigned short int)maxval);
      int idx = 0;
      for (yy = 0; yy < ysize; yy++ ){
        unsigned short int value  = (short int) (minval+ yy*delta);
        unsigned short int _value = swap( value );    // swap high and low bytes, the format expect big-endianism
          for(xx = 0; xx < xsize; xx++ ){
            sImage[idx++] = (value > maxval)?_maxval:_value;
          }
      }
      ptr = (void*)sImage;  
    }

  return ptr;
}