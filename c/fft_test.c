#include "fft.h"
#include <stdio.h>
/* Requires gcc flags -lfftw3 -lm */

int main() {
  int dimx, dimy,i,j;
  dimx=4; dimy=4;
  fftw_plan p;
  double complex *cmplx_a = (double complex *) fftw_malloc(sizeof(double complex)*dimx*dimy);
  double complex *real_a = (double complex *) fftw_malloc(sizeof(double complex)*dimx*dimy);
  double k[] = {0.7749, 0.3998, 0.9106, 0.1361, 
               0.8173, 0.2599, 0.1818, 0.8693, 
               0.8687, 0.8001, 0.2638, 0.5797, 
               0.0844, 0.4314, 0.1455, 0.5499};

  double complex a[]  = {0.7749+I*0.0, 0.3998+I*0.0, 0.9106+I*0.0, 0.1361+I*0.0, 
               0.8173+I*0.0, 0.2599+I*0.0, 0.1818+I*0.0, 0.8693+I*0.0, 
               0.8687+I*0.0, 0.8001+I*0.0, 0.2638+I*0.0, 0.5797+I*0.0, 
               0.0844+I*0.0, 0.4314+I*0.0, 0.1455+I*0.0, 0.5499+I*0.0};
  printf("INVERSE TEST:\n");
  p = fftw_plan_dft_2d(dimx, dimy, a, cmplx_a, 
                       FFTW_FORWARD,
                       FFTW_ESTIMATE);

  fftw_execute(p);

  for (i=0; i<dimx; i++) {
    for (j=0; j<dimy; j++) 
      printf(creal(cmplx_a[dimx*i+j]) < 0.0  ? "%s " : " %s  ", cfs(cmplx_a[dimx*i + j]));
    printf("\n");
  }

  
  p = fftw_plan_dft_2d(dimx, dimy, cmplx_a, real_a,
                       FFTW_BACKWARD,
                       FFTW_ESTIMATE);

  fftw_execute(p);
  printf("real_a:\n");
  for (i=0; i<dimx; i++) {
    for (j=0; j<dimy; j++) 
      printf(creal(real_a[dimx*i+j]) < 0.0  ? "%s " : " %s  ", cfs(real_a[dimx*i + j]));
    printf("\n");
  }
  double complex kc[] = {8.0732 + I*0.0000, 1.0436 + I*0.2438, 0.0208 + I*0.0000, 1.0436 - I*0.2438,
                 -0.2909 - I*0.9171,-0.2497 -I* 0.7399, 1.3969 - I*0.6213,-1.2315 - I*0.6533,
                 1.3942 + I*0.0000,-0.1052 - I*1.2120, 1.7838 + I*0.0000,-0.1052 + I*1.2120,
                 -0.2909 + I*0.9171,-1.2315 + I*0.6533, 1.3969 + I*0.6213,-0.2497 + I*0.7399};
  double *kk;
  double complex *kfft;
  kfft = fft2d_r2c(k, dimx, dimy);
  printf("kft:\n");
  for (i=0; i<dimx; i++) {
    for (j=0; j<dimy; j++) 
      printf(creal(kfft[dimx*i+j]) < 0.0  ? "%s " : " %s  ", cfs(kfft[dimx*i + j]));
    printf("\n");
  }



  
  kk = fft2d_c2r(kc, dimx, dimy);
  printf("kk:\n");
  for (i=0; i<dimx; i++) {
    for (j=0; j<dimy; j++) 
      printf(kk[dimx*i+j] < 0.0  ? "%lf " : " %lf  ", kk[dimx*i + j]);
    printf("\n");
  }

  

  return 0;
}
/* 
kfft:
 8.073200 + 0.000000i, 1.043600 + 0.243800i, 0.000000 + 0.000000i, 1.043600 + -0.243800i
-0.290900 + -0.917100i -0.249700 + -0.739900i,0.000000 + 0.000000i,-1.231500 + -0.653300i
 0.000000 + 0.000000i, 0.000000 + 0.000000i, 0.000000 + 0.000000i, 0.000000 + 0.000000i
-0.290900 + 0.917100i -1.231500 + 0.653300i,0.000000 + 0.000000i,-0.249700 + 0.739900i
*/



