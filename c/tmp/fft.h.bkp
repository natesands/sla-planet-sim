/*--------------------------------------------------------------------------------
 * 2D fast fourier transforms using FFTW2.
 *
 * Requires gcc flags -lfftw3 -lm
 *------------------------------------------------------------------------------*/


#include <complex.h>
#include <fftw3.h>
#include <time.h>    // random seed>
#include <stdlib.h>  // rand(); use rand() / RAND_MAX for uniform rv [0,1]
#include <math.h>    // M_PI, M_E

char buff[100];      // for cfs()

double complex* fft2d_r2c(double *f, int dimx, int dimy);
double * fft2d_c2r(double complex *f, int dimx, int dimy);
char *cfs(complex double c);


