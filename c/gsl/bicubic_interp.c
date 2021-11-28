
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

int main()
{

  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  const size_t N = 21;  /* num interp points */
  const double xa[] = { -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 }; 
  const double ya[] = { -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 }; 

  const size_t nx = sizeof(xa) / sizeof(double); /* x grid points */
  const size_t ny = sizeof(ya) / sizeof(double); /* y grid points */
  double *za = malloc(nx * ny * sizeof(double));
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();
  
  double R, z, eps;
  double X[N][N];
  double Y[N][N];
  double Z[N][N];
  eps = 2.2204e-16;
  size_t i, j;

  /* create mesh grid of query points,  including midpoints between sample values */  

  double query_pts[N];
  for (i = 0; i < N; i++)
    query_pts[i] = -5.0 + (double) i * (10.0 / (N -1.0));

 
  for (i=0; i < N; i++) 
    for (j=0; j < N; j++) {
      X[i][j] = query_pts[j];
      Y[i][j] = query_pts[i];
    }

   printf("X:\n"); 
    for (i=0; i < N; i++) {
      for (j=0; j < N; j++) 
        printf("%lf ", X[i][j]);
      printf("\n");
    }

   printf("Y:\n"); 
    for (i=0; i < N; i++) {
      for (j=0; j < N; j++) 
        printf("%lf \t", Y[i][j]);
      printf("\n");
    }

  /* initialize interpolator */
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

  for (i=0; i < N; i++)
    for (j=0; j < N; j++)
      Z[i][j] = gsl_interp2d_eval(bicubic, xa, ya, za, X[i][j], Y[i][j], xacc, yacc);

 printf("Z:\n"); 
  for (i=0; i < N; i++) {
    for (j=0; j < N; j++) 
      printf("%lf ", Z[i][j]);
    printf("\n");
  }

  gsl_interp2d_free(bicubic);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  free(za);
  return 0;

}




  
  
