#include <stdio.h>
#include <gsl/gsl_vector.h>

int main() {
  int i;
  gsl_vector *v = gsl_vector_alloc(3);
  for (i=0; i<3; i++)
    gsl_vector_set(v,i,1.23+i);

  for (i=0; i<v->size; i++)
    printf("v_%d = %g", i, gsl_vector_get(v,i));
}

