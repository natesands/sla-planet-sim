#include <complex.h>
#include <stdio.h>

char buff[100];

/* complex format string for printing complex numbers*/
char *cfs(double complex c) {
  int n;
  n = snprintf(buff, 100, "%f + %fi",creal(c), cimag(c));
  return buff;
}

double complex cmean1D(double complex *k, int size) {
  double complex sum = 0;
  for(int i=0; i<size; sum+=k[i++]);
  return sum / size;
}

double complex cmean2d(double complex k[3][3], int nx,int ny) {
  double complex sum = 0;
  for (int i=0; i<nx; i++)
    sum += cmean1D(k[i],ny);
  return sum / nx;
}


int main(){

  double complex a[3]= {1.0+3.0*I, 2.0, 1.0-5.0*I};
  double complex f[3][3]  ={ {1.0 +5*I,  2.0+6*I,  1.0+7*I},
                        {5.0,  6.0,  7.0},
                        { -1.0, -2.0, -1.0}};
  /* add */
  printf("1+3i + 1-3i = %s\n", cfs(1+3*I + 1-3*I));
  /* divide */
  printf("1+3i / 3 = %s\n", cfs((1+3*I)/3));
  /* mult */
  printf("1+3i * 1-5i = %s\n", cfs((1+3*I)*(1-5*I)));
  /* mean 1D */
  printf("mean{1+3*I, 2, 1-5*I}=%s\n", cfs(cmean1D(a,3)));
  /* mean 2d */
  printf("%s\n", cfs(cmean2d(f,3,3)));
  return 0;
}



