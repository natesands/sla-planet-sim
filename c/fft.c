/*--------------------------------------------------------------------------------
 * Implementation of 2D fast fourier transforms using FFTW2.
 *
 * Note: all 2D arrays are in row-major format.
 *
 * FFTW3 docs:
 * "The multi-dimensional transforms of FFTW, in general, compute simply the separable product of the given 1d transform along each dimension of the array. Since each of these transforms is unnormalized, computing the forward followed by the backward/inverse multi-dimensional transform will result in the original array scaled by the product of the normalization factors for each dimension (e.g. the product of the dimension sizes, for a multi-dimensional DFT)." 
 *
 * Requires gcc flags -lfftw3 -lm
 *------------------------------------------------------------------------------*/

#include "fft.h"

/*******************************************************************************/
double complex* fft2d_r2c(double *f, int dimx, int dimy) {
/********************************************************************************
Performs FFT from 2D array of reals to array of complex values. Returns 
pointer to new array. 
********************************************************************************/
  fftw_plan p;
  int i, j;
  double complex *ff, *g;    // return array
  ff = (double complex*) fftw_malloc(sizeof(double complex)*dimx*dimy); 
  g = (double complex*) fftw_malloc(sizeof(double complex)*dimx*dimy); 

  for (i=0; i<dimx*dimy; i++)    // copy f to ff
    ff[i] = (double complex) f[i];

  p = fftw_plan_dft_2d(dimx, dimy, ff, g, 
                       FFTW_FORWARD,
                       FFTW_ESTIMATE);
  fftw_execute(p);
  for(i=0; i<dimx; i++) {
    for(j=0; j<dimy; j++) {
      if (i==dimx/2 || j==dimy/2)   // TODO: why are these values set to zero?
        g[dimy*i+j] = 0.0;
    } // end for j
  } // end for i
  fftw_free(ff);
  return g;
} //fft2d_r2c


/*******************************************************************************/
double* fft2d_c2r(double complex *f, int dimx, int dimy) {
/********************************************************************************
Performs FFT from 2D array of complex values to array of reals. Returns 
pointer to new array. 

-Input: complex array f
  Copy f to temp complex array ff;
  Copy real part of ff to g;
  Free ff;
  return g;
********************************************************************************/
  fftw_plan p;
  int i, j;
  double complex *ff, *fff;  /* f is copied to ff; ff is modified, then transformed to fff */
  double *g;    /* return array */
  ff = (double complex*) fftw_malloc(sizeof(double complex)*dimx*dimy);
  fff = (double complex*) fftw_malloc(sizeof(double complex)*dimx*dimy);
  g = (double *) fftw_malloc(sizeof(double)*dimx*dimy); 

  for (i=0; i<dimx; i++) {
    for (j=0; j<dimy; j++) {
      if (i==dimx/2 || j==dimy/2)   // TODO: why are these values set to zero?
        ff[dimy*i+j] = 0.0 + I*0.0;
      else ff[dimy*i+j] = f[dimy*i+j];
    }
  }
  p = fftw_plan_dft_2d(dimx, dimy, ff, fff, FFTW_BACKWARD,FFTW_ESTIMATE);
  fftw_execute(p);             
  for (i=0; i<dimx*dimy; i++)
    g[i] = creal(fff[i]) / (double) (dimx*dimy);   // ff needs to be normalized
  return g;
} //fft2d_c2r

