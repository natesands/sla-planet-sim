#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define DIMX 256
#define DIMY 256

int main() {
  double arr[DIMX*DIMY];
  FILE *outFile;
  outFile = fopen("real_noise.txt","w");

  srand(time(NULL));
  for (int i = 0; i < DIMX*DIMY; i++) 
    fprintf(outFile, "%f ", (double) rand() / (double) RAND_MAX);

  fclose(outFile);

   return 0;
}


