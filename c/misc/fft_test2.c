/* Requires gcc flags -lfftw3 -lm 
 *  NOTE - ff2d_r2c and fft2_c2r zeros out columns and rows*
 * */
#include "fft.h"
#include "sla.h"
#include <stdio.h>

int main() {
  int dimx, dimy,i,j;
  dimx=6;
  dimy=4;
  fftw_plan p;
  double real_in[] = {0.7749, 0.3998, 0.9106, 0.1361, 
               0.8173, 0.2599, 0.1818, 0.8693, 
               0.8687, 0.8001, 0.2638, 0.5797,
              0.0844, 0.4314, 0.1455, 0.5499,
              0.3, .55, .234, .6765,
              .4567, .673, .788, .988};
  //double complex *cmplx_out = (double complex*) fftw_malloc(sizeof(double complex)*dimx*dimy);
  double complex *cmplx_out;
  double *real_out;
  double *real_in2 = malloc(sizeof(double)*dimx*dimy);
  for (i=0; i<dimx*dimy; i++)
    real_in2[i] = real_in[i%24];
  cmplx_out = fft2d_r2c(real_in, dimx, dimy);


  /* 

  printf("real in:\n");
  for (i=0; i<dimx; i++) 
    for (j=0; j<dimy; j++) 
      printf(j < dimy-1 ? "%lf " : "%lf\n", real_in[dimy*i+j]);
  cmplx_out = fft2d_r2c(real_in, dimx, dimy);
  printf("\ncomplex out:\n");
  for (i=0; i<dimx; i++) {
    for (j=0; j<dimy; j++) 
      printf(creal(cmplx_out[dimy*i+j]) < 0.0  ? "%s " : " %s  ", cfs(cmplx_out[dimy*i + j]));
    printf("\n");
  }
  real_out = fft2d_c2r(cmplx_out, dimx, dimy);
  printf("real out:\n");
  for (i=0; i<dimx; i++) 
    for (j=0; j<dimy; j++) 
      printf(j < dimy-1 ? "%lf " : "%lf\n", real_out[dimy*i+j]);
*/
}




