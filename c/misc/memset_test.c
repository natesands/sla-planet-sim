#include <stdio.h>
#include <string.h>
#include <stdlib.h>
int main() {
  size_t N = 10;
  int i;
  double *p = (double *) malloc(sizeof(double)*N);
  for (i=0; i < N; i++)
    p[i] = (double) i;
  memset(p, 0, 3*sizeof(double));
  for (i = 0; i<N; i++)
    printf("%f\n", p[i]);
  return 0;
}

