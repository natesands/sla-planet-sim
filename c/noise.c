#include <complex.h>
#include <stdio.h>
#include <time.h>   // for random seed
#include <stdlib.h>  // rand(); use rand() / RAND_MAX for uniform [0,1]
#include <math.h>  // M_PI, M_E

void noise(double k[10][10], int nx, int ny, double kmin, double kmax, int kpow, double rms_noise) {
  double complex **cnoise = (double complex **) malloc(sizeof(double complex*)*nx);
  for (int i=0;i<nx;i++)
    cnoise[i] = (double complex *) malloc(sizeof(double complex) *ny);
  int ii, jj;
  for (jj=0; jj<ny; jj++)
    for (ii=0; ii<nx; ii++) 
      if ((k[ii][jj] >=kmin) && (k[ii][jj] <=kmax))
        cnoise[ii][jj] = cpow(M_E, 2*I*M_PI*((double) rand()/ (double) RAND_MAX)) /
            k[ii][jj]*k[ii][jj];
  for (int i=0; i<ny; i++) {
    for (int j=0; j<nx; j++)
      printf("%f ", creal(cnoise[i][j]));
    printf("\n");
  }
  printf("\n\n");
  for (int i=0; i<ny; i++) {
    for (int j=0; j<nx; j++)
      printf("%f ", cimag(cnoise[i][j]));
    printf("\n");
  }
  printf("\n");
  return;
}

void noise1d(double *k, int n) {
  for (int i=0; i<n; i++)
    k[i] = (double) rand() / (double) RAND_MAX;
}
    
int main() {
  srand(time(NULL));
  int n=10;
  double *k = (double *) malloc(sizeof(double)*n);
  double kk[n][n];
  noise1d(k,n);
  for (int i=0; i<n; i++)
    printf("%f\n", k[i]);
 
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      kk[i][j] = (double) rand() / (double) RAND_MAX;
  
  noise(kk, n, n, .001, .8, 1, 1);


  return 0;
}



