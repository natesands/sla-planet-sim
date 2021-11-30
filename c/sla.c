/*--------------------------------------------------------------------------------
sla.c
Planetary formation simulation based on the research and Matlab source code contained in 
SIMULATING THE BIRTH OF PLANETS: A SPECTRAL SEMI-LAGRANGIAN HYDRODYNAMIC APPROACH,
a Master's thesis by Wendy Crumrine.

TODOS:  
-- malloc -> fftw_malloc
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
  
  //printf("X:\n");
 // printrmat(X, NX, NY);
  //printf("Y:\n");
  //printrmat(Y, NX, NY);
  //printf("KX:\n");
  //printrmat(KX,NX,NY);  
  //printf("KY:\n");
  //printrmat(KY,NX,NY);  
  //printf("K2:\n");
  //printrmat(K2,NX,NY);  

  hypvisc = (double*) malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX*NY; i++)
    hypvisc[i] = exp(-8.0*(pow(K2[i]/sq(M_PI/dx), hypviscpow/2)));
  //printf("hypvisc:\n");
  //printrmat(hypvisc, NX, NY);
  hypvisc[0] = 1.0;

  x_buf = add_buffer(X, NX, NY, bufx, bufy);
  for (i=0; i < bufx; i++)
    for (j=0; j < NY+2*bufy; j++) {
      x_buf[i*(NY+2*bufy) + j] -= LX;
      x_buf[(NX + bufx)*(NY+2*bufy) + i*(NY+2*bufy) + j] += LX;
    }
  y_buf = add_buffer(Y, NX, NY, bufx, bufy);

  //printf("x_buf:\n");
  //printrmat(x_buf, NX+2*bufx, NY+2*bufy);
 
  y_buf = add_buffer(Y, NX, NY, bufx, bufy);
  for (j=0; j < bufy; j++)
    for (i=0; i < NX+2*bufx; i++) {
      y_buf[i*(NY+2*bufy) + j] -= LY;
      y_buf[i*(NY+2*bufy) + NY+bufy+j] += LY;;
    } 

  //printf("y_buf:\n");
  //printrmat(y_buf, NX+2*bufx, NY+2*bufy);

  /* Initialize background shear */ 
  vxb = (double*) malloc(sizeof(double)*NX*NY);
  for (i=0; i<NX*NY; i++)
    vxb[i] = -shear * Y[i] * (0.5 * (tanh(8.0 * (Y[i] + 0.4*LY))-tanh(8.0* (Y[i] - 0.4*LY))));    // NB: shear=0 in MATLAB code

  printf("vxb:\n");
  printrmat(vxb, NX, NY);
  
  /* Initialize vorticity */
  wzq[0] = (double*) malloc(sizeof(double) * NX * NY);   // TODO: why is wzq[0] set to 0*x in the MATLAB code?
  memset(wzq[0], 0.0, sizeof(double) * NX * NY);

  /* TODO: why is wzq set to 0*x in the matlab code? The values are already 0*/

  double scalar_mult = 1;      // <-- set to 0 in MATLAB code
  for (i=0; i < NX; i++)
    for (j=0; j < NY; j++)
      wzq[0][i*NY + j] = scalar_mult * X[i*NY+j];
  printf("wzq:\n");
  printrmat(wzq[0], NX, NY);
 
  /* initialize random dust density */
  rho[0]=noise2d(msqrt(K2,NX*NY), NX, NY, M_PI / dx / 64.0, M_PI / dx / 2.0, 1.0, 1.0);
  for (i=0; i< NX*NY; i++) 
    rho[0][i] = 0.1 * (1.0 + .01*rho[0][i]);

  printf("rho t=0:\n");
  printrmat(rho[0], NX, NY);
  /* find velocity from vorticity via streamfunction at t=0 */
  update_velocity_via_streamfunc(0);
 
  printf("psi\n");
  printcmat(psi, NX, NY);

  printf("vx:\n");
  printrmat(vx, NX, NY);

  printf("vy:\n");
  printrmat(vy, NX, NY);
  printf("vxw_x:\n");
  printcmat(vxw_x, NX, NY);
  printf("vxw_y:\n");
  printcmat(vxw_y, NX, NY);

  /* compute gas pressure and dust drift velocity */
 // qx = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
 // qy = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*NX*NY);
 for (i=0; i<NX*NY; i++) {   // probably not necessary
   qx[i]=0.0; 
   qy[i]=0.0; 
 }

  update_drift_vel_gas_P(0);

  printf("qx:\n");
  printcmat(qx, NX, NY); 
  printf("qy:\n");
  printcmat(qy, NX, NY); 
  printf("crlq:\n");
  printrmat(crlq, NX, NY);
  printf("divq:\n");
  printrmat(divq, NX, NY);
  // TODO:  the values of the matrices are odd... (in MATLAB as well)

  /* initialize delx, dely, xi, yi, update displacements */
  real_mat_scalar_mult(delx, vx, dt, NX, NY);
  real_mat_scalar_mult(dely, vy, dt, NX, NY);
  update_xi_yi();
  iterate_displacements();
  printf("delx:\n");
  printrmat(delx, NX, NY);
  printf("dely:\n");
  printrmat(dely, NX, NY);
  /* initial advect forward step for wzq... */

  fftw_complex *tmp_wzq_cmplx;
  double *wz_buf, *tmp_wzq_real;
  wz_buf = add_buffer(wzq[0], NX, NY, bufx, bufy);
  printf("wz_buf 0:\n");
  printrmat(wz_buf, NX + 2*bufx, NY + 2*bufy);
  wzq[1] = (double*) fftw_malloc(sizeof(double) * NX * NY);  
  interpolate_grid(wzq[1], x_buf, y_buf, wz_buf, xi, yi);
  printf("interp wzq:\n");
  printrmat(wzq[1], NX, NY);
  add_real_mats(wzq[1], wzq[1], crlq, NX, NY);
  tmp_wzq_cmplx = fft2d_r2c(wzq[1], NX, NY);
  for (i=0; i < NX * NY; i++)
    tmp_wzq_cmplx[i] *= hypvisc[i];
  tmp_wzq_real = fft2d_c2r(tmp_wzq_cmplx, NX, NY);
  for (i=0; i < NX * NY; i++)
    wzq[1][i] = tmp_wzq_real[i];
  fftw_free(tmp_wzq_cmplx);
  fftw_free(tmp_wzq_real);
  fftw_free(wz_buf);
  printf("wzq %d\n", 1);
  printrmat(wzq[1], NX, NY);

  /* ...and rho. */
  fftw_complex *tmp_rho_cmplx;
  double *rho_buf, *tmp_rho_real;
  rho_buf = add_buffer(rho[0], NX, NY, bufx, bufy);
  printf("rho_buf:\n");
  printrmat(rho_buf, NX + 2*bufx, NY + 2*bufy);
  rho[1] = (double*) fftw_malloc(sizeof(double) * NX * NY);
  interpolate_grid(rho[1], x_buf, y_buf, rho_buf, xi, yi);
  printf("interp rho:\n");
  printrmat(rho[1], NX, NY);
  subtract_real_mats(rho[1], rho[1], divq, NX, NY);
  tmp_rho_cmplx = fft2d_r2c(rho[1], NX, NY);
  for (i=0; i < NX * NY; i++)
    tmp_rho_cmplx[i] *= hypvisc[i];
  tmp_rho_real = fft2d_c2r(tmp_rho_cmplx, NX, NY);
  for (i=0; i < NX * NY; i++)
    rho[1][i] = tmp_rho_real[i];
  fftw_free(tmp_rho_cmplx);
  fftw_free(tmp_rho_real);
  fftw_free(rho_buf);
  printf("rho %d\n", 1);
  printrmat(rho[1], NX, NY);
 

  /* Main loop */ /*starting t=1*/
  for (int timestep = 1; timestep < NT; timestep++) {
    printf("T: %d\n", timestep);
    update_velocity_via_streamfunc(timestep);
    printf("vx %d:\n", timestep);
    printrmat(vx, NX, NY);
    printf("vy %d:\n", timestep);
    printrmat(vy, NX, NY);
    update_drift_vel_gas_P(timestep);
    iterate_displacements();
    update_xi2_yi2();
    wzq_advect_step(timestep);
    rho_advect_step(timestep);

    printf("wzq %d\n", timestep);
    printrmat(wzq[timestep], NX, NY);
    printf("rho %d\n", timestep);
    printrmat(rho[timestep], NX, NY);

  }
  
  return 0;

}


