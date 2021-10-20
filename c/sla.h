/*--------------------------------------------------------------------------------
sla.h is an include file for planetary formation simulation, sla.c
--------------------------------------------------------------------------------*/
#include <stdlib.h>

/* Constants----------------------------------------------------------------------
NX, NY = number of grid points in the X,Y directions
LX, LY = dimensions of grid
NT = number of time steps
--------------------------------------------------------------------------------*/

#define NX 6
#define NY 6
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
double x_buf[NX+2*BUFX][NY+2*BUFY]; 
double y_buf[NX+2*BUFX][NY+2*BUFY]; 
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


/* Parameters-------------------------------------------------------------------
------------------------------------------------------------------------------*/

double rho0, omega, shear, dPdR, tau, dt;
int num_sl_disp_iter, num_pressure_iter, hypviscpow;


/* Functions-------------------------------------------------------------------
grid2d: initialize mesh grids -- equiv. to Matlab's ndgrid for d=2
sq: square 
add_buffer: 
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
  for (int i=0; i<Nx; i++)
    for (int j=0; j<Ny; j++) 
      Xb[i+bufx][j+bufy] = X[i][j];
  return Xb;
}





