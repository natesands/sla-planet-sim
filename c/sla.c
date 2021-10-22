/*--------------------------------------------------------------------------------
sla.c
Planetary formation simulation based on the research and Matlab source code contained in 
SIMULATING THE BIRTH OF PLANETS: A SPECTRAL SEMI-LAGRANGIAN HYDRODYNAMIC APPROACH,
a Master's thesis by Wendy Crumrine.
--------------------------------------------------------------------------------*/
#include "sla.h"
#include <stdio.h>
//#include <stdlib.h> 
#include <math.h>

int main() {
  // Initialize spatial and wave number grids

  
  // X,Y = physical coordinate mesh
  X = (double **) malloc(sizeof(double *) * NX);
  for (int i=0; i<NX; i++)
    X[i] = (double *) malloc(sizeof(double) * NY);
  Y = (double **) malloc(sizeof(double *) * NX);
  for (int i=0; i<NX; i++)
    Y[i] = (double *) malloc(sizeof(double) * NY);
  x = (double *) malloc(sizeof(double)*NX);
  y = (double *) malloc(sizeof(double)*NY);
  for (int i=0; i<NX; i++) 
    x[i] = (double) i*dx-.5*LX;
  for (int i=0; i<NY; i++)
    y[i] = (double) i*dy-.5*LY;

  grid2d(X, Y, NX, NY, x, y);

  // KX, KY = wave number mesh
  KX = (double **) malloc(sizeof(double *) * NX);
  for (int i=0; i<NX; i++)
    KX[i] = (double *) malloc(sizeof(double) * NY);
  KY = (double **) malloc(sizeof(double *) * NX);
  for (int i=0; i<NX; i++)
    KY[i] = (double *) malloc(sizeof(double) * NY);
  kx = (double *) malloc(sizeof(double)*NX);
  ky = (double *) malloc(sizeof(double)*NY);
  for (int i=0; i<NX; i++) {                         // TODO: double check indexing
    kx[i] = (double) (i<NX/2+1) ? i : i-NX;
    kx[i] = kx[i]*(2*M_PI/LX);
  }
  for (int i=0; i<NY; i++) {
    ky[i] = (double) (i<NY/2+1) ? i : i-NY;
    ky[i] = ky[i]*(2*M_PI/LY);
  }

  grid2d(KX, KY, NX, NY, kx, ky);


  // K2 = kx^2 + ky^2 for (kx,ky) in [KX,KY]
  K2 = (double **) malloc(sizeof(double *) * NX);
  for (int i=0; i<NX; i++)
    K2[i] = (double *) malloc(sizeof(double) * NY);
  for (int i=0; i<NX; i++)
    for (int j=0; j<NY; j++)
      K2[i][j] = KX[i][j]*KX[i][j] + KY[i][j]*KY[i][j];
  K2[0][0] = 1.0;

  printf("X:\n");
  for (int i=0; i<NX; i++) {
    for (int j=0; j<NY; j++) 
      printf("%f\t", X[i][j]);
    printf("\n");
  }
  
  printf("Y:\n");
  for (int i=0; i<NX; i++) {
    for (int j=0; j<NY; j++) 
      printf("%f\t", Y[i][j]);
    printf("\n");
  }

  printf("KX:\n");
  for (int i=0; i<NX; i++) {
    for (int j=0; j<NY; j++) 
      printf("%f\t", KX[i][j]);
    printf("\n");
  }
  
  printf("KY:\n");
  for (int i=0; i<NX; i++) {
    for (int j=0; j<NY; j++) 
      printf("%f\t", KY[i][j]);
    printf("\n");
  }

  printf("K2:\n");
  for (int i=0; i<NX; i++) {
    for (int j=0; j<NY; j++) 
      printf("%f\t", K2[i][j]);
    printf("\n");
  }

//  rho0 = 1.0;
//  omega = 1.0;
//  shear = -1.5*0;
//  dPdR = -0.10;
//  tau = 0.005;
//  dt = 1.0/8.0;
//  nt = 4096;
 // bufx = NX/16;
 // bufy = NY/16;
//  num_sl_disp_iter = 4;
//  num_pressure_iter = 4;
//  hypviscpow = 8;
  
//  hypvisc = (double **) malloc(sizeof(double *) * NX);
//  for (int i=0; i<NX; i++)
//    hypvisc[i] = (double *) malloc(sizeof(double) * NY);
  for (int i=0; i<NX; i++)
    for (int j=0; j<NY; j++) 
      hypvisc[i][j] = exp(-8.0*(pow(K2[i][j]/sq(M_PI/dx),(hypviscpow/2))));
  hypvisc[0][0]=1.0;


  printf("hypvisc:\n");
  for (int i=0; i<NX; i++) {
    for (int j=0; j<NY; j++) 
      printf("%f\t", hypvisc[i][j]);
    printf("\n");
  }

  
//  printf("qx:\n");
//  for (int i=0; i<NX; i++) {
//    for (int j=0; j<NY; j++) 
//      printf("%f\t", qx[i][j]);
//    printf("\n");
//  }
//  printf("BUFY: %d, BUFX: %d\n", BUFY, BUFX);
//
//  double **XX = (double **) malloc(sizeof(double *) * 3);
//  for (int i=0; i<3; i++) 
//    XX[i] = (double *) malloc(sizeof(double)*3);
//  for (int i=0; i<9; i++) 
//    //printf("%d, %d\n", i/3, i%3);
//   XX[i/3][i%3] = i+1;
//
//  int tmpbufx = 2;
//  int tmpbufy = 2; 
//  XX = add_buffer(XX, 3, 3, tmpbufx,tmpbufy);
//
//  printf("XX:\n");
//  for (int i=0; i<3+2*tmpbufx; i++) {
//    for (int j=0; j<3+2*tmpbufy; j++) 
//      printf("%f\t", XX[i][j]);
//    printf("\n");
//  }
//
  x_buf = add_buffer(X, NX, NY, BUFX, BUFY); // x_buf is (NX+2*BUFX) x (NY+2*BUFY)
  for (int i=0; i<BUFX; i++)
    for (int j=0; j<NY+2*BUFY; j++) {
      x_buf[i][j] -= LX;
      x_buf[NX+BUFX+i][j] += LX; 
    }

  printf("x_buf:\n");
  for (int i=0; i < NX+2*BUFX;i++) {
    for(int j=0; j<NY+2*BUFY;j++)
      printf("%f\t",x_buf[i][j]);
    printf("\n");
  }
  y_buf = add_buffer(Y, NX, NY, BUFX, BUFY); // y_buf is (NX+2*BUFX) x (NY+2*BUFY)
  for (int j=0; j<BUFY; j++)
    for (int i=0; i<NX+2*BUFX; i++) {
      y_buf[i][j] -= LY;
      y_buf[i][NY+BUFY+j] += LY;
    }

  printf("y_buf:\n");
  for (int i=0; i < NX+2*BUFX;i++) {
    for(int j=0; j<NY+2*BUFY;j++)
      printf("%f\t",y_buf[i][j]);
    printf("\n");
  }
 return 0;


}


