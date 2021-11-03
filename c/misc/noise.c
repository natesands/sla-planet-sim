#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <time.h>   // for random seed
#include <stdlib.h>  // rand(); use rand() / RAND_MAX for uniform [0,1]
#include <math.h>  // M_PI, M_E

char buff[100];

/* complex format string for printing complex numbers*/
char *cfs(fftw_complex c) {
  int n;
  n = snprintf(buff, 100, "%f + %fi",creal(c), cimag(c));
  return buff;
}
/*
void fft2d(double complex *f, double complex *g, int nx, int ny, int dir) {
  fftw_plan p;
  int i, j;
  double complex *ff= (double complex*) fftw_malloc(sizeof(double complex)*nx*ny);
  for (i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      ff[i*nx+j] = f[i*nx+j];

  if (dir < 0) {
    p = fftw_plan_dft_2d(nx, ny, ff, g, 
                       FFTW_FORWARD,
                       FFTW_ESTIMATE);
    fftw_execute(p);
    for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
        if (i==nx/2 || j==ny/2)   // TODO: why are these values set to zero?
          g[nx*i+j] = 0.0;
        else
          g[nx*i+j] = creal(g[nx*i+j]);
      } // end for j
    } // end for i
  } // end if
  else {
    for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
        if (i==nx/2 || j==ny/2)
          ff[nx*i+j] = 0.0;
      } // end for j
    } // end for i
    p = fftw_plan_dft_2d(nx, ny, ff, g, 
                       FFTW_BACKWARD,
                       FFTW_ESTIMATE);
    fftw_execute(p);
  } // end else
} // fft2d
*/

double complex cmean1D(double complex *k, int size) {
  double complex sum = 0;
  for(int i=0; i<size; sum+=k[i++]);
  return sum / size;
}

double complex cmean2d(double complex **k, int nx,int ny) {
  double complex sum = 0;
  for (int i=0; i<nx; i++)
    sum += cmean1D(k[i],ny);
  return sum / nx;
}

double complex** fft2d(double complex **f, int nx, int ny, int dir) {
  /* Performs fft2d in either direction depending on sign of dir.  
     If dir < 0: real physical space --> complex frequency space.
     If dir > 0: complex frequency space --> real physical space.
     Takes a 2D array ptr f and returns a 2D array array ptr g. */

  fftw_plan p;
  int i, j;
  double complex *gg, *ff;  // FFTW only works on single dimensional arrays in row-major order. 
  double complex **g;    // This will be the 2D array pointer returned.

  ff= (double complex*) fftw_malloc(sizeof(double complex)*nx*ny);
  gg= (double complex*) fftw_malloc(sizeof(double complex)*nx*ny);
  g = (double complex**) fftw_malloc(sizeof(double complex*)*nx);
  for (i=0; i<nx; i++)
    g[i] = (double complex*) fftw_malloc(sizeof(double complex)*ny);
 
  // copy f to ff 
  for (i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      ff[i*nx+j] = f[i][j];

  if (dir < 0) {
    p = fftw_plan_dft_2d(nx, ny, ff, gg, 
                       FFTW_FORWARD,
                       FFTW_ESTIMATE);
    fftw_execute(p);
    for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
        if (i==nx/2 || j==ny/2)   // TODO: why are these values set to zero?
          gg[nx*i+j] = 0.0;
        else
          gg[nx*i+j] = creal(gg[nx*i+j]);
      } // end for j
    } // end for i
  } // end if
  else {
    for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
        if (i==nx/2 || j==ny/2)
          ff[nx*i+j] = 0.0;
      } // end for j
    } // end for i
    p = fftw_plan_dft_2d(nx, ny, ff, gg, 
                       FFTW_BACKWARD,
                       FFTW_ESTIMATE);
    fftw_execute(p);
  } // end else
  // copy gg to g and return
  for (i=0; i<nx; i++)
    for(j=0; j<ny; j++)
      g[i][j] = f[i][j];

  return g;
} // fft2d

/* element-wise square of 2D array of complex doubles */ 
// TODO: modifies original array - don't use?
void csq2d(double complex **k, int dimx, int dimy) {
  for (int i=0; i<dimx; i++)
    for (int j=0; j<dimy; j++)
      k[i][j] *= k[i][j];
}

/* returns root mean squared of 2D array of complex doubles */
double complex crms2d(double complex **k, int dimx, int dimy) {
  double complex rms = 0;
  // sum square of each element
  for (int i=0; i<dimx; i++)
   for (int j=0; j<dimy; j++)
    rms += k[i][j]*k[i][j]; 
  // average
  rms /= dimx*dimy;
  return csqrt(rms);
}


/* Initialize grid of random noise in Fourier space */
// TODO: what is the type of the array that is passed to this function?
double complex** noise2d(double **k, int nx, int ny, double kmin, double kmax, int kpow, double rms_noise) {
 // TODO:  What function does rms_noise have? 
  double complex rms;
  double complex **noise = (double complex **) malloc(sizeof(double complex*)*nx);
  for (int i=0;i<nx;i++)
    noise[i] = (double complex *) malloc(sizeof(double complex) *ny);
  int ii, jj;
  for (jj=1; jj<ny; jj++)                         // TODO: why is indexing starting at 1?
    for (ii=1; ii<nx; ii++) 
      if ((k[ii][jj] >=kmin) && (k[ii][jj] <=kmax))
        noise[ii][jj] = cpow(M_E, 2*I*M_PI*((double) rand()/ (double) RAND_MAX)) /
            pow(k[ii][jj],kpow);
  noise = fft2d(noise, nx, ny, FFTW_FORWARD);
/*  for (int i=0; i<ny; i++) {
    for (int j=0; j<nx; j++)
      printf("%f ", creal(noise[i][j]));
    printf("\n");
  }
  printf("\n\n");
  for (int i=0; i<ny; i++) {
    for (int j=0; j<nx; j++)
      printf("%f ", cimag(noise[i][j]));
    printf("\n");
  }
  printf("\n");
  */
  rms = crms2d(noise, nx, ny);
  // noise = noise/rms * rms_noise
  for (int i=0; i<nx; i++)
    for (int j=0; j<ny; j++)
      noise[i][j] = noise[i][j]/rms * rms_noise; 
  return noise;
}

void noise1d(double *k, int n) {
  for (int i=0; i<n; i++)
    k[i] = (double) rand() / (double) RAND_MAX;
}
/*
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
*/  

int main() {
  srand(time(NULL));
  int n=10;
  double *k = (double *) malloc(sizeof(double)*n);
  double **kk = (double **) malloc(sizeof(double*)*n);
  for (int i=0; i<n; i++) 
    kk[i] = (double*) malloc(sizeof(double)*n);
  noise1d(k,n);
  for (int i=0; i<n; i++)
    printf("%f\n", k[i]);
   
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
      kk[i][j] = (double) rand() / (double) RAND_MAX;
  
  noise2d(kk, n, n, .001, .8, 1, 1);
   
  double complex a[] =  {32.0000 + 0.0000*I, -2.5000 - 0.8660*I,  -2.5000 + 0.8660*I,
    -7.0000 - 3.4641*I,  0.5000 + 0.8660*I,  -2.5000 - 2.5981*I,
    -7.0000 + 3.4641*I, -2.5000 + 2.5981*I,   0.5000 - 0.8660*I };

  double complex **aa = (double complex**) malloc(sizeof(double complex*)*3);
  for (int i=0; i<3; i++) {
    aa[i] = (double complex *) malloc(sizeof(double complex));
    for (int j=0; j<3; j++)
      aa[i][j] = a[3*i+j];
  }
  /* test crms2d */
  printf("rms = %s\n", cfs(crms2d(aa,3,3)));

  return 0;
}



