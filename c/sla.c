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
  X = (double*) malloc(sizeof(double)*NX*NY);
  Y = (double*) malloc(sizeof(double)*NX*NY);
  x = (double *) malloc(sizeof(double)*NX); 
  y = (double *) malloc(sizeof(double)*NY);
  for (i=0; i<NX; i++) {
    x[i] = (double) i*dx-.5*LX;
  }
  for (i=0; i<NY; i++) {
    y[i] = (double) i*dy-.5*LY;
  }
  grid2d(X, Y, NX, NY, x, y);

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
  K2 = (double*) malloc(sizeof(double*)*NX*NY);
  for (i=0; i<NX*NY; i++)
    K2[i] = KX[i]*KX[i] + KY[i]*KY[i];
  K2[0] = 1.0E+64;;   // TODO:  this was done in MATLAB code to prevent division by zero?
  
  printf("X:\n");
 // printrm_rowmaj(X, NX, NY);
  
  printf("Y:\n");
  //printrm_rowmaj(Y, NX, NY);
  printf("KX:\n");
  //printrm_rowmaj(KX,NX,NY);  
  printf("KY:\n");
  //printrm_rowmaj(KY,NX,NY);  
  printf("K2:\n");
  //printrm_rowmaj(K2,NX,NY);  


  hypvisc = (double*) malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX*NY; i++)
    hypvisc[i] = exp(-8.0*(pow(K2[i]/sq(M_PI/dx), hypviscpow/2)));
  printf("hypvisc:\n");
  //printrm_rowmaj(hypvisc, NX, NY);
  hypvisc[0] = 1.0;

  x_buf = add_buffer(X, NX, NY, bufx, bufy);
  for (i=0; i < bufx; i++)
    for (j=0; j < NY+2*bufy; j++) {
      x_buf[i*(NY+2*bufy) + j] -= LX;
      x_buf[(NX + bufx)*(NY+2*bufy) + i*(NY+2*bufy) + j] += LX;
    }
  y_buf = add_buffer(Y, NX, NY, bufx, bufy);

  printf("x_buf:\n");
  //printrm_rowmaj(x_buf, NX+2*bufx, NY+2*bufy);
 
  y_buf = add_buffer(Y, NX, NY, bufx, bufy);
  for (j=0; j < bufy; j++)
    for (i=0; i < NX+2*bufx; i++) {
      y_buf[i*(NY+2*bufy) + j] -= LY;
      y_buf[i*(NY+2*bufy) + NY+bufy+j] += LY;;
    } 

  printf("y_buf:\n");
  //printrm_rowmaj(y_buf, NX+2*bufx, NY+2*bufy);

  /* Initialize background shear */ 
  vxb = (double*) malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX*NY; i++)
    vxb[i] = -shear * Y[i] * (0.5 * (tanh(8.0 * (Y[i] + 0.4*LY))-tanh(8.0* (Y[i] - 0.4*LY))));    // NB: shear=0 in MATLAB code

  printf("vxb:\n");
  printrm_rowmaj(vxb, NX, NY);
  
  /* Initialize vorticity */
  wzq[0] = (double*) malloc(sizeof(double) * NX * NY);   // TODO: why is wzq[0] set to 0*x in the MATLAB code?
  memset(wzq[0], 0.0, sizeof(double) * NX * NY);

  /* TODO: why is wzq set to 0*x in the matlab code? The values are already 0*/

  double scalar_mult = 1;      // <-- set to 0 in MATLAB code
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

  printf("vx:\n");
  printrm_rowmaj(vx, NX, NY);

  printf("vy:\n");
  printrm_rowmaj(vy, NX, NY);
  printf("vxw_x:\n");
  printcm_rowmaj(vxw_x, NX, NY);
  printf("vxw_y:\n");
  printcm_rowmaj(vxw_y, NX, NY);


  /* compute gas pressure and dust drift velocity */
 // qx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
 // qy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
 for (i=0; i<NX*NY; i++) { 
   qx[i]=0.0; 
   qy[i]=0.0; 
 }

  update_drift_vel_gas_P(0);

  // TODO:  the values of the matrices are odd... (in MATLAB as well)

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
  fftw_complex *tmp_dPdx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
  fftw_complex *tmp_dPdy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
  double *tmp_dPdx_r = (double *) malloc(sizeof(double)*NX*NY);
  double *tmp_dPdy_r = (double *) malloc(sizeof(double)*NX*NY);
  fftw_complex *tmp_crlq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
  fftw_complex *tmp_divq = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
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


