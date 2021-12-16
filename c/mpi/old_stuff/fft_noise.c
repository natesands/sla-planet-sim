/*-------------------------------------------------------------------------------- 
Toy program to test parallelize FFTW3's 2d DFT, real -> complex.  
Must compile with linker flags gcc test1d.c -lfftw3 -lm 
Will read in file "real_noise.txt" containing 256 * 256 doubles 
representing a 2D grid and execute DFT.  Output is stored in array 
cmplx_out.  The file "complex_noise.txt" contains the expected
output.

Change DIMX, DIMY arrays to test smaller grids.  Uncomment lines
to print input/output.

Notes:
See
https://www.fftw.org/fftw3_doc/Multi_002dthreaded-FFTW.html

-------------------------------------------------------------------------------*/

#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIMX 256
#define DIMY 256

char *cfs(fftw_complex c);
char buff[100];
int main(int argc, char *argv[]) 
{
  FILE *real_in, *fout;
  fftw_plan p;
  fftw_complex *cmplx_out, *noise;
  double fp;
  int i;
  noise = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)* DIMX * DIMY);
  cmplx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DIMX * DIMY);
  real_in = fopen("real_noise.txt", "r"); 
  for (i=0; i < DIMX * DIMY; i++) {
    fscanf(real_in, "%lf", &fp);
    noise[i] = fp;
  }
  fclose(real_in);

//  for (i=0; i < DIMX * DIMY; i++)
//    printf("%f ", creal(noise[i]))

  p = fftw_plan_dft_2d(DIMX, DIMY, noise, cmplx_out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);


  fout = fopen("complex_noise.txt", "w");

  for (i=0; i < DIMX * DIMY ; i++)
    fprintf(fout, "%s ", cfs(cmplx_out[i]));

  fclose(fout);
  return 0;

}

/* complex format string for printing complex numbers*/
char *cfs(fftw_complex c) {
  int n;
  n = snprintf(buff, 100, "%f + %fi",creal(c), cimag(c));
  return buff;
}
