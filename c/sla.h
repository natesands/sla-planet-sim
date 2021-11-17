/*-------------------------------------------------------------------------------
sla.h is an include file for planetary formation simulation, sla.c,
based on the research and Matlab source code contained in 
SIMULATING THE BIRTH OF PLANETS: A SPECTRAL SEMI-LAGRANGIAN HYDRODYNAMIC APPROACH,
a Master's thesis by Wendy Crumrine.
--------------------------------------------------------------------------------*/
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>   // fast fourier transforms 
#include <stdio.h>
#include <math.h>
#include <time.h>   // for random seed
#include "fft.h"
/* Constants----------------------------------------------------------------------
NX, NY = number of grid points in the X,Y directions
LX, LY = dimensions of grid
NT = number of time steps
BUFX, BUFY = size of buffer space for 2d arrays
--------------------------------------------------------------------------------*/

#define NX 16  // TODO: set to 16 for testing
#define NY 16
#define LX 4.0
#define LY 4.0
#define NT 4096
#define BUFX (NX/16)
#define BUFY (NY/16)

/* Variables---------------------------------------------------------------------
dx, dy = grid cell dimensions
x, y, kx, ky = intermediate variables for X,Y,KX,KY
X,Y = physical coordinate mesh
KX, KY = wave number mesh
K2 = kx^2 + ky^2 for (kx,ky) in [KX,KY]
hypvisc
wzq : vorticity over time
-------------------------------------------------------------------------------*/
double dx = LX/NX;
double dy = LY/NY;
double *x,*y,*kx,*ky;
double **X,**Y, **KX, **KY, **K2;
double hypvisc[NX][NY];
double wzq[NT+1][NX*NY];
double complex *ipsiky, *negipsikx;
double *vx;
double *vy;
double *V2;
double vxb[NX][NY];
double complex *psi; 
double delx[NX][NY];
double dely[NX][NY];
double xi[NX][NY];  
double yi[NX][NY];  
double vx_buf[NX+2*BUFX][NY+2*BUFY];
double vy_buf[NX+2*BUFX][NY+2*BUFY];
double wz_buf[NX+2*BUFX][NY+2*BUFY];
//double x_buf[NX+2*BUFX][NY+2*BUFY]; 
double **x_buf;
//double y_buf[nx+2*bufx][ny+2*bufy]; 
double **y_buf;
//double rho[NT+1][NX][NY];     
double ***rho;
double complex *vxw_x;
double complex *vxw_y;;
double *temp_real;
double complex *temp_complex;
double nlxf[NX][NY]; 
double nlyf[NX][NY]; 
double hf[NX][NY];   
double dPdx[NX][NY]; 
double dPdy[NX][NY]; 
double qx[NX][NY];   
double qy[NX][NY];   
double divq[NX][NY]; 
double crlq[NX][NY]; 
double rho_buf[NX+2*BUFX][NY+2*BUFY]; 
double divq_buf[NX+2*BUFX][NY+2*BUFY];
double crlq_buf[NX+2*BUFX][NY+2*BUFY];
char buff[100];

/* Parameters-------------------------------------------------------------------
------------------------------------------------------------------------------*/

//double rho0, omega, shear, dPdR, tau, dt;
// int num_sl_disp_iter, num_pressure_iter, hypviscpow;

double rho0 = 1.0;
double omega = 1.0;
double shear = -1.5*0;
double dPdR = -0.10;
double tau = 0.005;
double dt = 1.0/8.0;
int num_sl_disp_iter = 4;
int num_pressure_iter = 4;
int hypviscpow = 8;

/* Functions-------------------------------------------------------------------
grid2d:     initialize mesh grids -- equiv. to Matlab's ndgrid for d=2
sq:         square value
add_buffer: pads Nx x Ny array by bufx and bufy, and tiles values. e.g.
            
            A B C
            D E F
            G H I

            becomes

         H I G H I G H
         B C A B C A B
         E F D E F D E
         H I G H I G H
         B C A B C A B

         for bufx=1, bufy=2.  
         returns pointer to (Nx+2*bufx) x (Ny+2*bufy) array.
cfs:  Returns c string for printing complex numbers. 
cmean1d / cmean2d:  mean of 1D/2D array of complex doubles
fft2d:  Performs fast fourier transform real->freq and freq->real
------------------------------------------------------------------------------*/

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

/* returns root mean squared of 2D array of doubles */
double complex rms2d(double **k, int dimx, int dimy) {
  double rms = 0;
  // sum square of each element
  for (int i=0; i<dimx; i++)
   for (int j=0; j<dimy; j++)
    rms += k[i][j]*k[i][j]; 
  // average
  rms /= dimx*dimy;
  return sqrt(rms);
}

/* take sqrt of elements of real matrix */
double** msqrt(double **M, int dimx, int dimy) {
  double **sqrtM = (double **) malloc(sizeof(double*)*dimx);
  for (int i=0; i<dimx; i++)
    sqrtM[i] = (double *) malloc(sizeof(double)*dimy);

  for (int i=0; i<dimx; i++)
    for (int j=0; j<dimy; j++)
      sqrtM[i][j] = sqrt(M[i][j]);

  return sqrtM;
}

void grid2d (double **X, double **Y, int Nx, int Ny, double *x, double *y) {
  // Initialize mesh grids X,Y each Nx by Ny
  for (int j=0; j < Ny; j++)
    for (int i=0; i < Nx; i++)
      X[i][j] = x[i];

  for (int i=0; i < Nx; i++)
    for (int j=0; j < Ny; j++)
      Y[i][j] = y[j];
}

double sq(double x) { return x*x;}

double** add_buffer(double **X, int Nx, int Ny, int bufx, int bufy) {
  double **Xb = (double **) malloc(sizeof(double *)*(Nx+2*bufx));
  for (int i=0; i<Nx+2*bufx; i++)
    Xb[i] = (double *) malloc(sizeof(double)*(Ny+2*bufy));
  for (int i=0; i<Nx; i++)  // pad X with zero-filled buffer
    for (int j=0; j<Ny; j++) 
      Xb[i+bufx][j+bufy] = X[i][j];
  for (int i=0; i<bufx; i++)
    for (int j=0; j<Ny; j++)
      Xb[i][bufy+j] = X[Nx-bufx+i][j];
  for (int i=0; i<bufx; i++)
    for (int j=0; j<Ny; j++)
      Xb[Nx+bufx+i][j+bufy] = X[i][j];
  for (int j=0; j<bufy; j++)
    for (int i=0; i<Nx+2*bufx; i++)
      Xb[i][j] = Xb[i][Ny+j];
  for (int j=0; j<bufy; j++)
    for (int i=0; i<Nx+2*bufx; i++)
      Xb[i][Ny+bufy+j] = Xb[i][bufy+j];
  return Xb;
}

/* complex format string for printing complex numbers*/
char *cfs(fftw_complex c) {
  int n;
  n = snprintf(buff, 100, "%f + %fi",creal(c), cimag(c));
  return buff;
}

/* print real matrix*/
void printrm(double **M, int dimx, int dimy)
{
  for (int i=0; i<dimx; i++) {
    for (int j=0; j<dimy; j++)
      printf(j != dimy-1 ? "%-6lf\t" : "%-6lf\n", M[i][j]);
  }
}

/* print real matrix in row major format */
void printrm_rowmaj(double *M, int dimx, int dimy)
{
  for (int i=0; i<dimx; i++) {
    for (int j=0; j<dimy; j++)
      printf(j != dimy-1 ? "%-6lf\t" : "%-6lf\n", M[i*dimx +j]);
  }
}
          
/* print complex matrix*/
void printcm(double complex **M, int dimx, int dimy)
{
  for (int i=0; i<dimx; i++) {
    for (int j=0; j<dimy; j++)
      printf(j != dimy-1 ? "%-6s\t" : "%-6s\n", cfs(M[i][j]));
  }
}

/* print complex matrix in row major format*/
void printcm_rowmaj(double complex *M, int dimx, int dimy)
{
  for (int i=0; i<dimx; i++) {
    for (int j=0; j<dimy; j++)
      printf(j != dimy-1 ? "%-6s\t" : "%-6s\n", cfs(M[i*dimx+j]));
  }
}

/* mean of array of complex values */
double complex cmean1D(double complex *k, int size) {
  double complex sum = 0;
  for(int i=0; i<size; sum+=k[i++]);
  return sum / size;
}
/* mean of matrix of complex values */
double complex cmean2d(double complex **k, int nx,int ny) {
  double complex sum = 0;
  for (int i=0; i<nx; i++)
    sum += cmean1D(k[i],ny);
  return sum / nx;
}

/* element-wise square of 2D array of complex doubles */ 
// TODO: modifies original array - don't use?
void csq2d(double complex **k, int dimx, int dimy) {
  for (int i=0; i<dimx; i++)
    for (int j=0; j<dimy; j++)
      k[i][j] *= k[i][j];
}

/* Initialize grid of random noise in Fourier space */
// TODO:  row major problem with fftw
// TODO: what is the type of the array that is passed to this function?
double ** noise2d(double **k, int nx, int ny, double kmin, double kmax, int kpow, double rms_noise) {
  int i, j;
 // TODO:  What function does rms_noise have? 
  double complex rms;
  double *noise_return_flat = (double *) malloc(sizeof(double)*nx*ny);
  double **noise_return = (double **) malloc(sizeof(double *)*nx);
  for (i=0; i<nx; i++)
    noise_return[i] = (double *) malloc(sizeof(double)*ny);
  double complex *noise = (double complex *) fftw_malloc(sizeof(double complex) * nx*ny);
  int ii, jj;
  for (jj=1; jj<ny; jj++) {                        // TODO: why is indexing starting at 1? 
    for (ii=1; ii<nx; ii++) { 
      if ((k[ii][jj] >=kmin) && (k[ii][jj] <=kmax)) {
        noise[ii*nx + jj] = cpow(M_E, 2*I*M_PI*((double) rand()/ (double) RAND_MAX)) /
            pow(k[ii][jj],kpow);
      }
    }
  }
  
  noise_return_flat = fft2d_c2r(noise, nx, ny);
  for (i=0; i<nx; i++) {
    for(j=0; j<ny; j++) {
      noise_return[i][j] = noise_return_flat[nx*i + j];
    }
  }

  rms = rms2d(noise_return, nx, ny); // TODO: rms2d takes rms of squared entries
  // noise = noise/rms * rms_noise
  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      noise_return[i][j] = noise_return[i][j] / rms * rms_noise;
  return noise_return;
}
