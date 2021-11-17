/*--------------------------------------------------------------------------------
sla.c
Planetary formation simulation based on the research and Matlab source code contained in 
SIMULATING THE BIRTH OF PLANETS: A SPECTRAL SEMI-LAGRANGIAN HYDRODYNAMIC APPROACH,
a Master's thesis by Wendy Crumrine.

TODOS:  
-- put everything in row major
-- get rid of i=0's
-- 
--------------------------------------------------------------------------------*/
#include "sla.h"
#include <stdio.h>
#include <math.h>
int main() {

  int i, j; 
  srand(time(NULL));

  // Initialize spatial and wave number grids
  // X,Y = physical coordinate mesh
  X = (double **) malloc(sizeof(double *) * NX);
  for (i=0; i<NX; i++)
    X[i] = (double *) malloc(sizeof(double) * NY);
  Y = (double **) malloc(sizeof(double *) * NX);
  for (i=0; i<NX; i++)
    Y[i] = (double *) malloc(sizeof(double) * NY); x = (double *) malloc(sizeof(double)*NX); y = (double *) malloc(sizeof(double)*NY);
  for (i=0; i<NX; i++) 
    x[i] = (double) i*dx-.5*LX;
  for (i=0; i<NY; i++)
    y[i] = (double) i*dy-.5*LY;

  grid2d(X, Y, NX, NY, x, y);

  // KX, KY = wave number mesh
  KX = (double **) malloc(sizeof(double *) * NX);
  for (i=0; i<NX; i++)
    KX[i] = (double *) malloc(sizeof(double) * NY);
  KY = (double **) malloc(sizeof(double *) * NX);
  for (i=0; i<NX; i++)
    KY[i] = (double *) malloc(sizeof(double) * NY);
  kx = (double *) malloc(sizeof(double)*NX);
  ky = (double *) malloc(sizeof(double)*NY);
  for (i=0; i<NX; i++) {                         // TODO: double check indexing
    kx[i] = (double) (i<NX/2+1) ? i : i-NX;
    kx[i] = kx[i]*(2*M_PI/LX);
  }
  for (i=0; i<NY; i++) {
    ky[i] = (double) (i<NY/2+1) ? i : i-NY;
    ky[i] = ky[i]*(2*M_PI/LY);
  }

  grid2d(KX, KY, NX, NY, kx, ky);

  // K2 = kx^2 + ky^2 for (kx,ky) in [KX,KY]
  K2 = (double **) malloc(sizeof(double *) * NX);
  for (i=0; i<NX; i++)
    K2[i] = (double *) malloc(sizeof(double) * NY);
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++)
      K2[i][j] = KX[i][j]*KX[i][j] + KY[i][j]*KY[i][j];
  K2[0][0] = 1.0;

  printf("X:\n");
  for (i=0; i<NX; i++) {
    for (j=0; j<NY; j++) 
      printf("%f\t", X[i][j]);
    printf("\n");
  }
  
  printf("Y:\n");
  for (i=0; i<NX; i++) {
    for (j=0; j<NY; j++) 
      printf("%f\t", Y[i][j]);
    printf("\n");
  }

  printf("KX:\n");
  for (i=0; i<NX; i++) {
    for (j=0; j<NY; j++) 
      printf("%f\t", KX[i][j]);
    printf("\n");
  }
  
  printf("KY:\n");
  for (i=0; i<NX; i++) {
    for (j=0; j<NY; j++) 
      printf("%f\t", KY[i][j]);
    printf("\n");
  }

  printf("K2:\n");
  for (i=0; i<NX; i++) {
    for (j=0; j<NY; j++) 
      printf("%f\t", K2[i][j]);
    printf("\n");
  }

//  rho0 = 1.0;
//  omega = 1.0;
    shear = -1.5*1;
//  dPdR = -0.10;
//  tau = 0.005;
//  dt = 1.0/8.0;
//  nt = 4096;
 // ufx = NX/16;
 // bufy = NY/16;
//  num_sl_disp_iter = 4;
//  num_pressure_iter = 4;
//  hypviscpow = 8;
  
//  hypvisc = (double **) malloc(sizeof(double *) * NX);
//  for (i=0; i<NX; i++)
//    hypvisc[i] = (double *) malloc(sizeof(double) * NY);
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++) 
      hypvisc[i][j] = exp(-8.0*(pow(K2[i][j]/sq(M_PI/dx),(hypviscpow/2))));
  hypvisc[0][0]=1.0;


  printf("hypvisc:\n");
  for (i=0; i<NX; i++) {
    for (j=0; j<NY; j++) 
      printf("%f\t", hypvisc[i][j]);
    printf("\n");
  }

  
//  printf("qx:\n");
//  for (i=0; i<NX; i++) {
//    for (j=0; j<NY; j++) 
//      printf("%f\t", qx[i][j]);
//    printf("\n");
//  }
//  printf("BUFY: %d, BUFX: %d\n", BUFY, BUFX);
//
//  double **XX = (double **) malloc(sizeof(double *) * 3);
//  for (i=0; i<3; i++) 
//    XX[i] = (double *) malloc(sizeof(double)*3);
//  for (i=0; i<9; i++) 
//    //printf("%d, %d\n", i/3, i%3);
//   XX[i/3][i%3] = i+1;
//
//  int tmpbufx = 2;
//  int tmpbufy = 2; 
//  XX = add_buffer(XX, 3, 3, tmpbufx,tmpbufy);
//
//  printf("XX:\n");
//  for (i=0; i<3+2*tmpbufx; i++) {
//    for (j=0; j<3+2*tmpbufy; j++) 
//      printf("%f\t", XX[i][j]);
//    printf("\n");
//  }
//
  x_buf = add_buffer(X, NX, NY, BUFX, BUFY); // x_buf is (NX+2*BUFX) x (NY+2*BUFY)
  for (i=0; i<BUFX; i++)
    for (j=0; j<NY+2*BUFY; j++) {
      x_buf[i][j] -= LX;
      x_buf[NX+BUFX+i][j] += LX; 
    }

  printf("x_buf:\n");
  for (i=0; i < NX+2*BUFX;i++) {
    for(j=0; j<NY+2*BUFY;j++)
      printf("%f\t",x_buf[i][j]);
    printf("\n");
  }
  y_buf = add_buffer(Y, NX, NY, BUFX, BUFY); // y_buf is (NX+2*BUFX) x (NY+2*BUFY)
  for (j=0; j<BUFY; j++)
    for (i=0; i<NX+2*BUFX; i++) {
      y_buf[i][j] -= LY;
      y_buf[i][NY+BUFY+j] += LY;
    }

  printf("y_buf:\n");
  for (i=0; i < NX+2*BUFX;i++) {
    for(j=0; j<NY+2*BUFY;j++)
      printf("%f\t",y_buf[i][j]);
    printf("\n");
  }
  /* initialize background shear */
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++) 
      vxb[i][j] = -shear*Y[i][j]*(0.5*(tanh(8.0*(Y[i][j]+0.4*LY))-tanh(8.0*(Y[i][j]-0.4*LY))));

  printf("vxb:\n");
  for (i=0; i < NX;i++) {
    for(j=0; j<NY;j++)
      printf("%f\t",vxb[i][j]);
    printf("\n");
  }
  
  /* initialize vorticity */
  double scalar_mult = 0.1;    /* TODO: why is this set to 0 in the matlab code?? */
  for (i=0; i < NX; i++)
    for (j=0; j < NY; j++)
      wzq[0][i*NX + j] = scalar_mult * X[i][j];
  

  /* initialize random dust density */
  rho = (double ***) malloc(sizeof(double **) * NT);
  rho[0]=noise2d(msqrt(K2,NX,NY), NX, NY, M_PI / dx / 64.0, M_PI / dx / 2.0, 1.0, 1.0);
  for (i=0; i< NX; i++) 
    for (j=0; j<NY; j++)
      rho[0][i][j] = 0.1 * (1.0 + .01*rho[0][i][j]);

  printf("rho t=0:\n");
  for (i=0; i < NX;i++) {
    for(j=0; j<NY;j++)
      printf("%f\t",rho[0][i][j]);
    printf("\n");
  }

  /* find velocity from vorticity via streamfunction */
  psi = fft2d_r2c(wzq[0], NX, NY);
  ipsiky = (double complex*) fftw_malloc(sizeof(double complex)*NX*NY);
  negipsikx = (double complex*) fftw_malloc(sizeof(double complex)*NX*NY);
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++) {
      psi[i*NX+j] = psi[i*NX+j] / K2[i][j];
      ipsiky[i*NX+j] = I*psi[i*NX+j]*KY[i][j];
      negipsikx[i*NX+j] = -I*psi[i*NX+j]*KX[i][j];
    }
  vx = fft2d_c2r(ipsiky, NX, NY);
  vy = fft2d_c2r(negipsikx, NX, NY);
  free(ipsiky);
  free(negipsikx);
  printf("psi\n"); 
  printcm_rowmaj(psi, NX, NY);
  printf("vx:\n");
  printrm_rowmaj(vx, NX, NY);
  printf("vy:\n");
  printrm_rowmaj(vy, NX, NY);
  V2 = (double *) malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      V2[i*NX+j] = vx[i*NX+j]*vx[i*NX+j] + vy[i*NX+j]*vy[i*NX+j];
  vxw_x = fft2d_r2c(V2, NX, NY);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      vxw_x[i*NX+j] = -0.5 * vxw_x[i*NX+j];
  vxw_y = (double complex*) fftw_malloc(sizeof(double complex)*NX*NY);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++) {
      vxw_y[i*NX+j] = I*KY[i][j]*vxw_x[i*NX+j];
      vxw_x[i*NX+j] = I*KX[i][j]*vxw_x[i*NX+j];
    }
  temp_real = (double *) malloc(sizeof(double)*NX*NY);
  temp_complex = (double complex *) fftw_malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      temp_real[i*NX+j] = vy[i*NX+j]*(wzq[0][i*NX+j]+2.0*(omega+shear));
  temp_complex = fft2d_r2c(temp_real, NX, NY);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      vxw_x[i*NX+j] = vxw_x[i*NX+j] + temp_complex[i*NX+j];
  free(temp_real);
  free(temp_complex);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      temp_real[i*NX+j] = vx[i*NX+j]*(wzq[0][i*NX+j]+2.0*(omega));
  temp_complex = fft2d_r2c(temp_real, NX, NY);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      vxw_y[i*NX+j] = vxw_y[i*NX+j] - temp_complex[i*NX+j];
  free(temp_real);
  free(temp_complex);

  return 0;

}


