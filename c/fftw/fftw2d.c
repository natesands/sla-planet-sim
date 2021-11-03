/*-------------------------------------------------------------------------------- 
Test FFTW3's 2d DFT.  Must compile with linker flags gcc test1d.c -lfftw3 -lm 
Notes:
-in and out arrays can be the same for an in-place transform
-if program performs many transforms of the same size and initialization time is not important, use FFTW_MEASURE
-must create the plan before initializing the input
-if <complex.h> is included before <fftw3.h> then fftw_complex is the native double-precision complex type 
and you can manipulate it with ordinary arithmetic
-------------------------------------------------------------------------------*/

#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIMX 3
#define DIMY 3

char buff[100];

/* complex format string for printing complex numbers*/
char *cfs(fftw_complex c) {
  int n;
  n = snprintf(buff, 100, "%f + %fi",creal(c), cimag(c));
  return buff;
}

int main() {
  fftw_plan p;
  int i, j;
  double complex f[]  ={ 1.0,  2.0,  1.0,
                         0.0,  0.0,  0.0,
                        -1.0, -2.0, -1.0};
  double complex *out = (double complex*) fftw_malloc(sizeof(double complex)*DIMX*DIMY);


  p = fftw_plan_dft_2d(DIMX, DIMY,  f, out, FFTW_FORWARD, FFTW_ESTIMATE);

 fftw_execute(p);
 printf("plan executed\n");
 printf("f:\n");
 for(i=0; i<DIMX; i++)
  for(j=0; j<DIMY; j++)
    printf("%s%c", cfs(f[DIMX*i+j]), j == DIMY-1 ? '\n' : ' ');
 printf("\n");
 printf("out:\n");
 for(i=0; i<DIMX; i++)
  for(j=0; j<DIMY; j++)
    printf("%s%c", cfs(out[DIMX*i+j]), j == DIMY-1 ? '\n' : ' ');
 printf("\n");
 return 0;
}
/*--------------------------------------------------------------------------------
 OUTPUT:
0.000000 + 0.000000i 0.000000 + 0.000000i 0.000000 + 0.000000i
6.000000 + -3.464102i -1.500000 + -0.866025i -0.000000 + 1.732051i
6.000000 + 3.464102i -0.000000 + -1.732051i -1.500000 + 0.866025i
--------------------------------------------------------------------------------*/  

  

