/*-------------------------------------------------------------------------------
sla.h is an include file for planetary formation simulation, sla.c,
based on the research and Matlab source code contained in 
SIMULATING THE BIRTH OF PLANETS: A SPECTRAL SEMI-LAGRANGIAN HYDRODYNAMIC APPROACH,
a Master's thesis by Wendy Crumrine.
--------------------------------------------------------------------------------*/
#include <stdlib.h>

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
grid2d: initialize mesh grids -- equiv. to Matlab's ndgrid for d=2
sq: square value
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





