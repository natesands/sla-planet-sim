#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

int
main()
{
  const gsl_interp2d_type *T = gsl_interp2d_bilinear;
  const size_t N = 11;             /* number of points to interpolate */
  size_t i, j;
  double *xa = malloc(sizeof(double)*N);
  double *ya = malloc(sizeof(double)*N);
  for (i=0; i < N; i++) {
    xa[i] = -5.0 + i;
    ya[i] = -5.0 + i;
  } 

  double *za = malloc(N * N * sizeof(double));
  gsl_spline2d *spline = gsl_spline2d_alloc(T, N, N);
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();


  /* set z grid values */
  for (i=0 ; i<N; i++)
    for (j=0; j<N; j++) {
      double R = sqrt(xa[i]*xa[i] + ya[i]*ya[i]);
      gsl_spline2d_set(spline, za, i, j, sin(R) / (R + 0.000000001)) ;
    }
  /* initialize interpolation */
  gsl_spline2d_init(spline, xa, ya, za, N, N);

  /* interpolate N values in x and y and print out grid for plotting */

  for (i = 0; i < 10; i++)
    {
      double xi = -5 + (double) i*0.5;

      for (j = 0; j < 10; j++) 
        {
          double yj = -5 + (double) j*0.5;
          double zij = gsl_spline2d_eval(spline, xi, yj, xacc, yacc);

          printf("%f %f %f\n", xi, yj, zij);
        }
      printf("\n");
    }
  gsl_spline2d_free(spline);
  gsl_interp_accel_free(xacc);
  gsl_interp_accel_free(yacc);
  free(za);

  return 0;
}
