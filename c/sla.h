#ifndef SLA_H
#define SLA_H
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
#include <string.h> // memset

/* Constants----------------------------------------------------------------------
NX, NY = number of grid points in the X,Y directions
LX, LY = dimensions of grid
NT = number of time steps
BUFX, BUFY = size of buffer space for 2d arrays
--------------------------------------------------------------------------------*/

#define NX 32  // TODO: NX, NY reset to 256
#define NY 16 
#define LX 4.0
#define LY 4.0
#define NT 4096
#define BUFX (NX/16)
#define BUFY (NY/16)

/* Parameters-------------------------------------------------------------------
------------------------------------------------------------------------------*/

//double rho0, omega, shear, dPdR, tau, dt;
// int num_sl_disp_iter, num_pressure_iter, hypviscpow;

double rho0 = 1.0;
double omega = 1.0;
double shear = -1.5*1;  // set to zero in orig implementation
double dPdR = -0.10;
double tau = 0.005;
double dt = 1.0/8.0;
int num_sl_disp_iter = 4;
int num_pressure_iter = 4;
int hypviscpow = 8;
int bufx = BUFX;   // TODO: get rid of this
int bufy = BUFY;
int ufx = NX/16;


/* Variables---------------------------------------------------------------------
dx, dy = grid cell dimensions
x, y, kx, ky = intermediate variables for X,Y,KX,KY
X,Y = physical coordinate mesh
KX, KY = wave number mesh
K2 = kx^2 + ky^2 for (kx,ky) in [KX,KY]
hypvisc
wzq : vorticity 
rho : dust density
-------------------------------------------------------------------------------*/
double dx = LX/NX;
double dy = LY/NY;
double *x,*y,*kx,*ky;
//double **X,**Y, **KX, **KY, **K2;
double *X, *Y, *KX, *KY, *K2;
double *hypvisc;
double *vxb;
double* wzq[NT+1];  //  wzq, rho both (NT+1) x (NX*NY) 2D grids
double* rho[NT+1];
double fac1[NX*NY];
double *vx;
double *vy;
double *V2;
double **vx_plus_vxb;  // TODO: get rid of this
fftw_complex *psi; 
fftw_complex *tmp_cplx_arr;
double *tmp_real_arr;
double delx[NX][NY];
double dely[NX][NY];
double xi[NX][NY];
double yi[NX][NY];  
double **vx_buf;
double vy_buf[NX+2*BUFX][NY+2*BUFY];
double wz_buf[NX+2*BUFX][NY+2*BUFY];
//double x_buf[NX+2*BUFX][NY+2*BUFY]; 
double *x_buf;
//double y_buf[nx+2*bufx][ny+2*bufy]; 
double *y_buf;
fftw_complex vxw_x[NX*NY];
fftw_complex vxw_y[NX*NY];
fftw_complex qx[NX*NY], qy[NX*NY];
fftw_complex nlxf[NX*NY], nlyf[NX*NY];
fftw_complex hf[NX*NY];
double divq[NX*NY];
double crlq[NX*NY];

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
fftw_complex crms2d(fftw_complex *a, int dimx, int dimy) {
  fftw_complex rms = 0;
  // sum square of each element
  for (int i=0; i<dimx; i++)
   for (int j=0; j<dimy; j++)
    rms += a[i*NY+j]*a[i*NY+j]; 
  // average
  rms /= dimx*dimy;
  return csqrt(rms);
}

/* returns root mean squared of 2D array of doubles */
fftw_complex rms2d(double *a, int dimx, int dimy) {
  double rms = 0;
  // sum square of each element
  for (int i=0; i<dimx; i++)
   for (int j=0; j<dimy; j++)
    rms += a[i*dimy+j]*a[i*dimy+j]; 
  // average
  rms /= dimx*dimy;
  return sqrt(rms);
}

/*  square root of entrys of an array */
double* msqrt(double *M, int size) {
  double *sqrtM = (double*) malloc(sizeof(double*)*size);
  for (int i=0; i<size; i++)
    sqrtM[i] = sqrt(M[i]);
  return sqrtM;
}

/* take sqrt of elements of real matrix */
/*
double** msqrt(double **M, int dimx, int dimy) {
  double **sqrtM = (double **) malloc(sizeof(double*)*dimx);
  for (int i=0; i<dimx; i++)
    sqrtM[i] = (double *) malloc(sizeof(double)*dimy);

  for (int i=0; i<dimx; i++)
    for (int j=0; j<dimy; j++)
      sqrtM[i][j] = sqrt(M[i][j]);

  return sqrtM;
}
*/

void grid2d(double *X, double *Y, int dimx, int dimy, double *x, double *y) {
  int i,j;
  for (i=0; i<dimx; i++)
    for (j=0; j<dimy; j++) {
      X[i*dimy+j] = x[i];
      Y[i*dimy+j] = y[j];
    }
 // printf("\n");
}
/*

void grid2d (double **X, double **Y, int Nx, int Ny, double *x, double *y) {
  // Initialize mesh grids X,Y each Nx by Ny
  for (int j=0; j < Ny; j++)
    for (int i=0; i < Nx; i++)
      X[i][j] = x[i];

  for (int i=0; i < Nx; i++)
    for (int j=0; j < Ny; j++)
      Y[i][j] = y[j];
}
*/

double sq(double x) { return x*x;}

double* add_buffer(double *X, int dimx, int dimy, int bufx, int bufy) {
  int i,j;
  int new_dimx = 2*bufx + dimx;
  int new_dimy = 2*bufy + dimy;
  double *Xb = (double *) malloc(sizeof(double) * new_dimx * new_dimy);
  memset(Xb, 0.0, sizeof(double) * new_dimx * new_dimy);
  for (i=0; i<dimx; i++)  // pad X with zero-filled buffer
    for (j=0; j<dimy; j++) 
      Xb[new_dimy * (i + bufx)+bufy+j] = X[i*dimy + j];
  for (i=0; i<bufx; i++)
    for (j=0; j<dimy; j++) {
      Xb[i*new_dimy+bufy+j] = X[(dimx-bufx+i)*dimy+j];
      Xb[new_dimy*(i+bufx+dimx) + bufy+j] = X[i*dimy + j]; 
    }
  for (i=0; i<new_dimx; i++)
    for (j=0; j<bufy; j++) {
      Xb[i*new_dimy+j] = Xb[i*new_dimy+dimy+j];
      Xb[i*new_dimy+dimy+bufy+j] = Xb[i*new_dimy+bufy+j]; 
    }
  return Xb;
}
/*
double** add_buffer(double **X, int Nx, int Ny, int bufx, int bufy) {
  double **Xb = (double **) malloc(sizeof(double *)*(Nx+2*bufx));
  for (int i=0; i<Nx+2*bufx; i++)
    Xb[i] = (double *) malloc(sizeof(double)*(Ny+2*bufy));
  for (int i=0; i<Nx; i++)
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
*/

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
void printrmat(double *M, int dimx, int dimy)
{
  for (int i=0; i<dimx; i++) {
    for (int j=0; j<dimy; j++)
      printf(j != dimy-1 ? "%-6lf\t" : "%-6lf\n", M[i*dimy +j]);
  }
}
          
/* print complex matrix*/
void printcm(fftw_complex **M, int dimx, int dimy)
{
  for (int i=0; i<dimx; i++) {
    for (int j=0; j<dimy; j++)
      printf(j != dimy-1 ? "%-6s\t" : "%-6s\n", cfs(M[i][j]));
  }
}

/* print complex matrix in row major format*/
void printcmat(fftw_complex *M, int dimx, int dimy)
{
  for (int i=0; i<dimx; i++) {
    for (int j=0; j<dimy; j++)
      printf(j != dimy-1 ? "%-6s\t" : "%-6s\n", cfs(M[i*dimy+j]));
  }
}

/* mean of array of complex values */
fftw_complex cmean1D(fftw_complex *k, int size) {
  fftw_complex sum = 0;
  for(int i=0; i<size; sum+=k[i++]);
  return sum / size;
}
/* mean of matrix of complex values */
fftw_complex cmean2d(fftw_complex **k, int dimx,int dimy) {
  fftw_complex sum = 0;
  for (int i=0; i<dimx; i++)
    sum += cmean1D(k[i],dimy);
  return sum / dimx;
}

/* element-wise square of 2D array of complex doubles */ 
// TODO: modifies original array - don't use?
void csq2d(fftw_complex **k, int dimx, int dimy) {
  for (int i=0; i<dimx; i++)
    for (int j=0; j<dimy; j++)
      k[i][j] *= k[i][j];
}


/* creates dimx-by-dimy row-major array of noise */
double* noise2d(double *k, int dimx, int dimy, double kmin, double kmax, int kpow, double rms_noise) {
  int i, j;
  fftw_complex rms;  // TODO:  What function does rms_noise have? 
  double *noise_return;
  fftw_complex *noise = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * dimx*dimy);
  int ii, jj;
  for (jj=1; jj<dimy; jj++) {                        // TODO: why is indexing starting at 1? 
    for (ii=1; ii<dimx; ii++) { 
      if ( (k[ii*dimy + jj] >= kmin) && (k[ii*dimy + jj] <=kmax)) {
        noise[ii*dimy + jj] = cpow(M_E, 2*I*M_PI*((double) rand()/ (double) RAND_MAX)) /
            pow(k[ii*dimy + jj],kpow);
      }
    }
  }
  noise_return = fft2d_c2r(noise, dimx, dimy);
  rms = rms2d(noise_return, dimx, dimy); // TODO: rms2d takes rms of squared entries
  for (i=0; i<dimx; i++)
    for (j=0; j<dimy; j++)
      noise_return[i*dimy+j] = noise_return[i*dimy+j] / rms * rms_noise;
  return noise_return;
}

/* find velocity from vorticity via streamfunction */
void update_velocity_via_streamfunc(int timestep) {
  int i;
  //fftw_free(vxw_x); 
  //fftw_free(vxw_y);
  fftw_complex *vxw_x_tmp, *vxw_y_tmp;
  psi = fft2d_r2c(wzq[timestep], NX, NY);
  for (i=0; i<NX*NY; i++) 
    psi[i] /= K2[i];

  
  tmp_cplx_arr = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
  for (i=0; i<NX*NY; i++)
    tmp_cplx_arr[i] = I * KY[i] * psi[i];

  vx = fft2d_c2r(tmp_cplx_arr, NX, NY);

  for (i=0; i<NX*NY; i++)
    tmp_cplx_arr[i] = -I * KX[i] * psi[i];
  vy = fft2d_c2r(tmp_cplx_arr, NX, NY);

  tmp_real_arr = (double*) malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX*NY; i++)
    tmp_real_arr[i] = vx[i]*vx[i] + vy[i]*vy[i];
   
  //vxw_x = fft2d_r2c(tmp_real_arr, NX, NY); 
  vxw_x_tmp = fft2d_r2c(tmp_real_arr, NX, NY); 
  for (i=0; i<NX*NY; i++)
    vxw_x[i] = vxw_x_tmp[i];
 // vxw_y = (fftw_complex*) fftw_malloc(sizeof(fftw_complex*) * NX * NY);  // TODO:  shouldn't need to reallocate each time
  for (i=0; i<NX*NY; i++) {
  vxw_x[i] *= -0.5;
  vxw_y[i] = I*KY[i]*vxw_x[i];
  vxw_x[i] = I*KX[i]*vxw_x[i];
  tmp_real_arr[i] = vy[i] * (wzq[timestep][i] + 2.0 * (omega + shear));
  }

  printf("tmp_real_arr:\n");
  //printrmat(tmp_real_arr, NX, NY);

  //free(tmp_cplx_arr);
  tmp_cplx_arr = fft2d_r2c(tmp_real_arr, NX, NY);
  
  for (i=0; i<NX*NY; i++) {
    vxw_x[i] += tmp_cplx_arr[i];
    tmp_real_arr[i] = vx[i] * (wzq[timestep][i] + 2.0 * omega);
  } 
  printf("tmp_real_arr:\n");
  //printrmat(tmp_real_arr, NX, NY);

  printf("vxw_x:\n");
  //printcmat(vxw_x, NX, NY);
  // OK tmp_real_arr, vxw_x
  //free(tmp_cplx_arr);
  //tmp_cplx_arr = fft2d_r2c(tmp_real_arr, NX, NY);
  //printf("tmp_cplx_arr:\n");
  //printcmat(tmp_cplx_arr, NX, NY);
  //for (i=0; i<NX*NY; i++)
  //  vxw_y[i] -= tmp_cplx_arr[i];
  //fftw_free(tmp_cplx_arr);
  //free(tmp_real_arr);

} 
/* updates qx, qy, divq, crlq */ 
void update_drift_vel_gas_P(int timestep) {
  double* rho_frame = rho[timestep];
  fftw_complex *i_kx_hf, *i_ky_hf, *tmp_complex, *qx_tmp, *qy_tmp;
  double *dPdx, *dPdy, *tmp_real, *tmp_divq, *tmp_crlq;   
  int i, j, k;
  i_kx_hf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
  i_ky_hf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
  tmp_real = (double*) fftw_malloc(sizeof(double)*NX*NY);
  tmp_complex= (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++) 
      fac1[i*NY+j] = (rho_frame[i*NY+j] / rho0) / ((1.0 + rho_frame[i*NY+j] / rho0)*(1.0 + rho_frame[i*NY+j] / rho0)
          + (2.0*omega*tau)*(2.0*omega*tau));    // fac1 off from matlab ~0.001
  printf("fac1:\n");
  printrmat(fac1, NX, NY);
  for (k=0; k < num_pressure_iter; k++) {
    for (i=0; i<NX*NY; i++) {
      nlxf[i] = vxw_x[i] + qx[i];
      nlyf[i] = vxw_y[i] + qy[i];
      hf[i] = -I*(KX[i] * nlxf[i] + KY[i] * nlyf[i]) / K2[i];
      i_kx_hf[i] = I*KX[i]*hf[i];
      i_ky_hf[i] = I*KY[i]*hf[i];
    }
    dPdx = fft2d_c2r(i_kx_hf, NX, NY);
    dPdy = fft2d_c2r(i_ky_hf, NX, NY);
    for (i=0; i<NX*NY; i++)
      dPdy[i] += dPdR;
    for (i=0; i<NX*NY; i++)
      tmp_real[i] = fac1[i]*( (1.0 + rho_frame[i] / rho0) * dPdx[i] + 2.0*omega*tau*dPdy[i]);
    qx_tmp = fft2d_r2c(tmp_real, NX, NY);
    for (i=0; i<NX*NY; i++)
      tmp_real[i] = fac1[i]*( (1.0 + rho_frame[i] / rho0) * dPdy[i] - 2.0 * omega * tau * dPdx[i]);
    qy_tmp = fft2d_r2c(tmp_real, NX, NY);
    for (i=0; i<NX*NY; i++) {
      qx[i] = qx_tmp[i];
      qy[i] = qy_tmp[i];
    }
  }
  for (i=0; i<NX*NY; i++)
    tmp_complex[i] = I*(KX[i] * qx[i] + KY[i] * qy[i]);
  tmp_divq = fft2d_c2r(tmp_complex, NX, NY);
  for (i=0; i<NX*NY; i++)
    divq[i] = dt * tmp_divq[i] * rho0 * tau;
  for (i=0; i<NX*NY; i++)
    tmp_complex[i] = I*(KX[i] * qy[i] - KY[i] * qx[i]);
  tmp_crlq = fft2d_c2r(tmp_complex, NX, NY);
  for (i=0; i<NX*NY; i++)
    crlq[i] = dt * tmp_crlq[i];

  fftw_free(i_kx_hf);
  fftw_free(i_ky_hf);
  fftw_free(tmp_real);
  fftw_free(tmp_complex);
  fftw_free(dPdx);
  fftw_free(dPdy);
  fftw_free(qx_tmp);
  fftw_free(qy_tmp);
  fftw_free(tmp_divq);
  fftw_free(tmp_crlq);
}
/*
% compute gas pressure and dust drift velocity
fac1 = (rho(:,:,tt)/rho0)./((1.0+rho(:,:,tt)/rho0).^2+(2.0*omega*tau)^2);
for jj=1:num_pressure_iter
    nlxf = vxw_x+qx;
    nlyf = vxw_y+qy;
    hf   = -1i*(kx.*nlxf+ky.*nlyf)./k2;
    dPdx = fft2d(1i*kx.*hf,+2);
    dPdy = fft2d(1i*ky.*hf,+2)+dPdR;
    qx   = fft2d(fac1.*((1.0+rho(:,:,tt)/rho0).*dPdx+2.0*omega*tau*dPdy),-2);
    qy   = fft2d(fac1.*((1.0+rho(:,:,tt)/rho0).*dPdy-2.0*omega*tau*dPdx),-2);
end
divq = dt*fft2d(1i*(kx.*qx+ky.*qy),+2)*(rho0*tau);
crlq = dt*fft2d(1i*(kx.*qy-ky.*qx),+2);
*/
#endif /* SLA_H */
