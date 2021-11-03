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

/* Constants----------------------------------------------------------------------
NX, NY = number of grid points in the X,Y directions
LX, LY = dimensions of grid
NT = number of time steps
BUFX, BUFY = size of buffer space for 2d arrays
--------------------------------------------------------------------------------*/

#define NX 6
#define NY 6
#define LX 4.0
#define LY 4.0
#define NT 4096
#define BUFX (NX/3)
#define BUFY (NY/3)

/* Variables---------------------------------------------------------------------
dx, dy = grid cell dimensions
x, y, kx, ky = intermediate variables for X,Y,KX,KY
X,Y = physical coordinate mesh
KX, KY = wave number mesh
K2 = kx^2 + ky^2 for (kx,ky) in [KX,KY]
-------------------------------------------------------------------------------*/
double dx = LX/NX;
double dy = LY/NY;
double *x,*y,*kx,*ky;
double **X,**Y, **KX, **KY, **K2;
double hypvisc[NX][NY];
double wzq[NX][NY][NT+1];
double vx[NX][NY];
double vxb[NX][NY];
double vy[NX][NY];  
double psi[NX][NY]; 
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
double rho[NX][NY][NT+1];     
double vxw_x[NX][NY];
double vxw_y[NX][NY];
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

double complex cmean1D(double complex *k, int size) {
  double complex sum = 0;
  for(int i=0; i<size; sum+=k[i++]);
  return sum / size;
}

double complex cmean2d(double complex k[3][3], int nx,int ny) {
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
