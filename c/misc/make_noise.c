#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define DIMX 10
#define DIMY 10

int main() {
  double arr[DIMX*DIMY];
  FILE *outFile;
  outFile = fopen("randarray.txt","w");

  srand(time(NULL));
  for (int i = 0; i < DIMX*DIMY; i++) 
    fprintf(outFile, "%f ", (double) rand() / (double) RAND_MAX);

  fclose(outFile);

   return 0;
}


