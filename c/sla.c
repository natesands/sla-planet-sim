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


  /* Initialize space grids X, Y using points x, y */
  /*
  X = (double **) malloc(sizeof(double *) * NX);
  for (i=0; i<NX; i++)
    X[i] = (double *) malloc(sizeof(double) * NY);
  Y = (double **) malloc(sizeof(double *) * NX);
  for (i=0; i<NX; i++)
    Y[i] = (double *) malloc(sizeof(double) * NY); 
  x = (double *) malloc(sizeof(double)*NX); 
  y = (double *) malloc(sizeof(double)*NY);
  for (i=0; i<NX; i++) 
    x[i] = (double) i*dx-.5*LX;
  for (i=0; i<NY; i++)
    y[i] = (double) i*dy-.5*LY;
  grid2d(X, Y, NX, NY, x, y);
  */
  
  /* Initialize space grids X, Y using points x, y */
  X = (double*) malloc(sizeof(double)*NX*NY);
  Y = (double*) malloc(sizeof(double)*NX*NY);
  x = (double *) malloc(sizeof(double)*NX); 
  y = (double *) malloc(sizeof(double)*NY);
  for (i=0; i<NX; i++) {
    x[i] = (double) i*dx-.5*LX;
    printf("x:%f\n", x[i]); 
  }
  for (i=0; i<NY; i++) {
    y[i] = (double) i*dy-.5*LY;
    printf("y:%f\n", y[i]);
  }
  grid2d(X, Y, NX, NY, x, y);

  /*  Initialize wave number grids KX, KY using points kx, ky */
  /*
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
  */

  
  /*  Initialize wave number grids KX, KY using points kx, ky */
  KX = (double*) malloc(sizeof(double)*NX*NY);
  KY = (double*) malloc(sizeof(double)*NX*NY);
  kx = (double*) malloc(sizeof(double)*NX);
  ky = (double*) malloc(sizeof(double)*NY);
  for (i=0; i<NX; i++) {                         // TODO: double check indexing
    kx[i] = (double) (i<NX/2+1) ? i : i-NX;
    kx[i] = kx[i]*(2*M_PI/LX);
  }
  for (i=0; i<NY; i++) {
    ky[i] = (double) (i<NY/2+1) ? i : i-NY;
    ky[i] = ky[i]*(2*M_PI/LY);
  }
  grid2d(KX, KY, NX, NY, kx, ky);

  /* Initialize grid of squared wave numbers K2 = kx^2 + ky^2 for (kx,ky) in [KX,KY] */
  /*
  K2 = (double **) malloc(sizeof(double *) * NX);
  for (i=0; i<NX; i++)
    K2[i] = (double *) malloc(sizeof(double) * NY);
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++)
      K2[i][j] = KX[i][j]*KX[i][j] + KY[i][j]*KY[i][j];
  K2[0][0] = 1.0;
  */
  /* Initialize grid of squared wave numbers K2 = kx^2 + ky^2 for (kx,ky) in [KX,KY] */
  K2 = (double*) malloc(sizeof(double*)*NX*NY);
  for (i=0; i<NX*NY; i++)
    K2[i] = KX[i]*KX[i] + KY[i]*KY[i];
  K2[0] = 1.0E+64;;   // TODO:  this was done in MATLAB code to prevent division by zero?
  
  printf("X:\n");
  printrm_rowmaj(X, NX, NY);
 /* for (i=0; i<NX; i++) {
    for (j=0; j<NY; j++) 
      printf("%f\t", X[i][j]);
    printf("\n");
  } */
  
  printf("Y:\n");
  printrm_rowmaj(Y, NX, NY);
  /*for (i=0; i<NX; i++) {
    for (j=0; j<NY; j++) 
      printf("%f\t", Y[i][j]);
    printf("\n");
  }*/
  printf("KX:\n");
  printrm_rowmaj(KX,NX,NY);  
 
  printf("KY:\n");
  printrm_rowmaj(KY,NX,NY);  

  printf("K2:\n");
  printrm_rowmaj(K2,NX,NY);  
/*
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
*/

//    hypvisc[i] = (double *) malloc(sizeof(double) * NY);
/*
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
*/
  hypvisc = (double*) malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX*NY; i++)
    hypvisc[i] = exp(-8.0*(pow(K2[i]/sq(M_PI/dx), hypviscpow/2)));
  printf("hypvisc:\n");
  printrm_rowmaj(hypvisc, NX, NY);
  hypvisc[0] = 1.0;

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
  x_buf = add_buffer(X, NX, NY, bufx, bufy);
  for (i=0; i < bufx; i++)
    for (j=0; j < NY+2*bufy; j++) {
      x_buf[i*(NY+2*bufy) + j] -= LX;
      x_buf[(NX + bufx)*(NY+2*bufy) + i*(NY+2*bufy) + j] += LX;
    }
  y_buf = add_buffer(Y, NX, NY, bufx, bufy);

  printf("x_buf:\n");
  printrm_rowmaj(x_buf, NX+2*bufx, NY+2*bufy);
 
  y_buf = add_buffer(Y, NX, NY, bufx, bufy);
  for (j=0; j < bufy; j++)
    for (i=0; i < NX+2*bufx; i++) {
      y_buf[i*(NY+2*bufy) + j] -= LY;
      y_buf[i*(NY+2*bufy) + NY+bufy+j] += LY;;
    } 

  printf("y_buf:\n");
  printrm_rowmaj(y_buf, NX+2*bufx, NY+2*bufy);

  /*
  x_buf = add_buffer(X, NX, NY, BUFX, BUFY); // x_buf is (NX+2*BUFX) x (NY+2*BUFY)
  for (i=0; i<BUFX; i++)
    for (j=0; j<NY+2*BUFY; j++) {
      x_buf[i][j] -= LX;
      x_buf[NX+BUFX+i][j] += LX; 
    }
  
  y_buf = add_buffer(Y, NX, NY, BUFX, BUFY); 
  for (j=0; j<BUFY; j++)
    for (i=0; i<NX+2*BUFX; i++) {
      y_buf[i][j] -= LY;
      y_buf[i][NY+BUFY+j] += LY;
    }

  
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++) 
      vxb[i][j] = -shear*Y[i][j]*(0.5*(tanh(8.0*(Y[i][j]+0.4*LY))-tanh(8.0*(Y[i][j]-0.4*LY))));
  */

  /* Initialize background shear */ 
  vxb = (double*) malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX*NY; i++)
    vxb[i] = -shear * Y[i] * (0.5 * (tanh(8.0 * (Y[i] + 0.4*LY))-tanh(8.0* (Y[i] - 0.4*LY))));    // NB: shear=0 in MATLAB code

  printf("vxb:\n");
  printrm_rowmaj(vxb, NX, NY);

  
  /* Initialize vorticity */
  wzq[0] = (double*) malloc(sizeof(double) * NX * NY);   // TODO: why is wzq[0] set to 0*x in the MATLAB code?
  //memset(wzq[0], 0.0, sizeof(double) * NX * NY);

  /* TODO: why is wzq set to 0*x in the matlab code? The values are already 0*/

 
  double scalar_mult = 1;      
  for (i=0; i < NX; i++)
    for (j=0; j < NY; j++)
      wzq[0][i*NY + j] = scalar_mult * X[i*NY+j];
  printf("wzq:\n");
  printrm_rowmaj(wzq[0], NX, NY);
  
 
  /* initialize random dust density */
  rho[0]=noise2d(msqrt(K2,NX*NY), NX, NY, M_PI / dx / 64.0, M_PI / dx / 2.0, 1.0, 1.0);
  for (i=0; i< NX*NY; i++) 
    rho[0][i] = 0.1 * (1.0 + .01*rho[0][i]);

  printf("rho t=0:\n");
  printrm_rowmaj(rho[0], NX, NY);
  /* find velocity from vorticity via streamfunction at t=0 */
  update_velocity_via_streamfunc(0);
  
  printf("psi\n");
  printcm_rowmaj(psi, NX, NY);
/*
  printf("vxw_x:\n");
  printcm_rowmaj(vxw_x, NX, NY);

  printf("vxw_y:\n");
  printcm_rowmaj(vxw_y, NX, NY);

  printf("vx:\n");
  printrm_rowmaj(vx, NX, NY);

  printf("vy:\n");
  printrm_rowmaj(vy, NX, NY);
*/



  // TODO:  the values of the matrices are odd... (in MATLAB as well)
  /*
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
  */
  // fftw_free(ipsiky);
  // fftw_free(negipsikx);
  /*
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
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      temp_real[i*NX+j] = vy[i*NX+j]*(wzq[0][i*NX+j]+2.0*(omega+shear));
  temp_complex = fft2d_r2c(temp_real, NX, NY);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      vxw_x[i*NX+j] = vxw_x[i*NX+j] + temp_complex[i*NX+j];
      */
  //fftw_free(temp_complex);
  /*
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      temp_real[i*NX+j] = vx[i*NX+j]*(wzq[0][i*NX+j]+2.0*(omega));
  temp_complex = fft2d_r2c(temp_real, NX, NY);
  for (i=0; i<NX; i++)
    for(j=0; j<NY; j++)
      vxw_y[i*NX+j] = vxw_y[i*NX+j] - temp_complex[i*NX+j];
      */
  //free(temp_real);
  //fftw_free(temp_complex);
  /*
  printf("vxw_x\n"); 
  printcm_rowmaj(vxw_x, NX, NY);
  printf("vxw_y\n");
  printcm_rowmaj(vxw_y, NX, NY);
*/
  /* compute gas pressure and dust drift velocity */
  /*
  fac1 = (double *) malloc(sizeof(double)*NX*NY);
  double* rho_frame = rho[0];
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++) 
      fac1[i*NX+j] = (rho_frame[i*NX+j] / rho0) / ((1.0 + rho_frame[i*NX+j] / rho0)*(1.0 + rho_frame[i*NX+j] / rho0)
          + (2.0*omega*tau)*(2.0*omega*tau));
 
  printf("fac1:\n"); 
  printrm_rowmaj(fac1, NX, NY);
  double complex *tmp_dPdx = (double complex*) fftw_malloc(sizeof(double complex)*NX*NY);
  double complex *tmp_dPdy = (double complex*) fftw_malloc(sizeof(double complex)*NX*NY);
  double *tmp_dPdx_r = (double *) malloc(sizeof(double)*NX*NY);
  double *tmp_dPdy_r = (double *) malloc(sizeof(double)*NX*NY);
  double complex *tmp_crlq = (double complex*) fftw_malloc(sizeof(double complex)*NX*NY);
  double complex *tmp_divq = (double complex*) fftw_malloc(sizeof(double complex)*NX*NY);
  for (int k=0; k < num_pressure_iter; k++) {
    for (i=0; i<NX; i++)
      for(j=0; j<NY; j++) {
        nlxf[i*NX+j] = vxw_x[i*NX+j] + (k==0 ? 0 : qx[i*NX+j]);
        nlyf[i*NX+j] = vxw_y[i*NX+j] + (k==0 ? 0 : qy[i*NX+j]);
        hf[i*NX+j] = -I*(KX[i][j]*nlxf[i*NX+j]+KY[i][j]*nlyf[i*NX+j])/K2[i][j];
        tmp_dPdx[i*NX+j] = I*KX[i][j]*hf[i*NX+j];
        tmp_dPdy[i*NX+j] = I*KY[i][j]*hf[i*NX+j];
      }
    dPdx = fft2d_c2r(tmp_dPdx, NX, NY);
    dPdy = fft2d_c2r(tmp_dPdy, NX, NY);
    for (i=0; i<NX; i++)
      for(j=0; j<NY; j++)
        dPdy[i*NX+j] += dPdR;
    for (i=0; i<NX; i++)
      for(j=0; j<NY; j++) {
        tmp_dPdx_r[i*NX+j] = fac1[i*NX+j]*((1.0 + rho_frame[i*NX+j] / rho0)*dPdx[i*NX+j] + 2.0*omega*tau*dPdy[i*NX+j]);
        tmp_dPdy_r[i*NX+j] = fac1[i*NX+j]*((1.0 + rho_frame[i*NX+j] / rho0)*dPdy[i*NX+j] - 2.0*omega*tau*dPdx[i*NX+j]);
      }
    fftw_free(qx);
    fftw_free(qy);
    qx = fft2d_r2c(tmp_dPdx_r, NX, NY);
    qy = fft2d_r2c(tmp_dPdy_r, NX, NY);
  }
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++) {
      tmp_divq[i*NX+j] = I*(KX[i][j]*qx[i*NX+j] + KY[i][j]*qy[i*NX+j]);
      tmp_crlq[i*NX+j] = I*(KX[i][j]*qy[i*NX+j] - KY[i][j]*qx[i*NX+j]);
    }
  divq = fft2d_c2r(tmp_divq, NX, NY);
  crlq = fft2d_c2r(tmp_crlq, NY, NY);
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++) {
      divq[i*NX+j] = dt * divq[i*NX+j]*rho0*tau;
      crlq[i*NX+j] = dt * crlq[i*NX+j];
    }
  vx_plus_vxb = (double **) malloc(sizeof(double *)*NX);
  for (i=0; i<NX; i++)
    vx_plus_vxb[i] = (double *) malloc(sizeof(double)*NY);
  for (i=0; i<NX; i++)
    for (j=0; j<NY; j++)
      vx_plus_vxb[i][j] = vx[i*NX+j] + vxb[i][j];
  vx_buf = add_buffer(vx_plus_vxb, NX, NY, bufx, bufy);
  printf("vx_buf:\n");
  printrm(vx_plus_vxb,NX+2*bufx,NY+2*bufy);

    

  fftw_free(tmp_dPdx);
  fftw_free(tmp_dPdy);
  free(tmp_dPdx_r);
  free(tmp_dPdy_r);
  fftw_free(tmp_divq);
  fftw_free(tmp_crlq);
  printf("qx:\n"); 
  printcm_rowmaj(qx, NX, NY);
  printf("qy:\n"); 
  printcm_rowmaj(qy, NX, NY);
  printf("divq:\n");
  printrm_rowmaj(divq, NX, NY);
  printf("crlq:\n");
  printrm_rowmaj(crlq, NX, NY);
  */
  return 0;

}


