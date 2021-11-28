#include <stdio.h>
#include <stdlib.h>


#define NX 32
#define NY 16

int main()
{
  double *a = malloc(sizeof(double)*NX*NY);

  for (int i=0; i<NX; i++)
    for (int j=0; j<NY; j++)
      a[i*NY+j] = (double) i;

  for (int i=0; i<NX*NY; i++)
    printf("%f\n", a[i]);
  return 0;
}
