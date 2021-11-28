
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

int main()
{

  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  const size_t N = 21;  /* num interp points */
  //const double xa[] = { -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 }; /* define unit square */
 // const double ya[] = { -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 }; /* define unit square */
  double xa[21];
  double ya[21];
  const size_t nx = sizeof(xa) / sizeof(double); /* x grid points */
  const size_t ny = sizeof(ya) / sizeof(double); /* y grid points */
  double *za = malloc(nx * ny * sizeof(double));
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  size_t i, j;
  
  for (i=0; i < N; i++) {
    xa[i] = -5.0 + (double) i * (10.0 / (N - 1.0));
    ya[i] = xa[i];
  }

  double R, z, eps;
  double X[ny][nx];
  double Y[ny][nx];
  double Z[ny][nx];
  eps = 2.2204e-16;
  for (j=0; j < N; j++)
    for (i=0; i < N; i++) {
      X[j][i] = xa[i];
      Y[j][i] = ya[j];
    }
 printf("X:\n"); 
  for (i=0; i < ny; i++) {
    for (j=0; j < nx; j++) 
      printf("%lf ", X[i][j]);
    printf("\n");
  }

 printf("Y:\n"); 
  for (i=0; i < ny; i++) {
    for (j=0; j < nx; j++) 
      printf("%lf \t", Y[i][j]);
    printf("\n");
  }

  gsl_interp2d *bicubic = gsl_interp2d_alloc(T, nx, ny);
  for(i=0; i<nx; i++)
    for(j=0; j<ny; j++) {
      R = sqrt(xa[i] * xa[i] + ya[j] * ya[j]);
      z = sin(R) / (R+eps);
      gsl_interp2d_set(bicubic, za,i, j, z);
    }
  printf("za:\n");
  for (i=0; i < ny; i++) {
    for (j=0; j < nx; j++) 
      printf("%lf ", za[j*nx+i]);
    printf("\n");
  }

  gsl_interp2d_init(bicubic, xa, ya, za, nx, ny);
  for (i=0; i < ny; i++)
    for (j=0; j < nx; j++)
      Z[i][j] = gsl_interp2d_eval(bicubic, xa, ya, za, X[i][j], Y[i][j], xacc, yacc);

  gsl_interp2d_free(bicubic);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  free(za);
 
 printf("Z:\n"); 
  for (i=0; i < ny; i++) {
    for (j=0; j < nx; j++) 
      printf("%lf ", Z[i][j]);
    printf("\n");
  }

  return 0;
}

/* given m x_i values, and n y_j values, 
 * X is an n-by-m grid whose rows are x1, x2, ..., xm
 * Y is an n-by-m grid whose columns are y1, y2, ... yn
 * Z is an n-by-m grid whose entries are zij = f(xi, yj)
 */
