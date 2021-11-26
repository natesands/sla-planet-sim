#include <stdio.h>
#include "/Users/ironchefnate/Documents/USC/CSCI_596/project/sla-planet-sim/c/sla.h"
#include <stdlib.h>

#define DIMX 4
#define DIMY 3
#define XBUF 2
#define YBUF 3
int main() {
  double *X = (double*) malloc(sizeof(double)*DIMX*DIMY);
  double *Xb;
  int i, j;
  for (i=0; i<DIMX; i++)
    for (j=0; j<DIMY; j++)
      X[i*DIMY+j] = (double) i*DIMY+j;
  printf("X:\n"); 
  printrm_rowmaj(X,DIMX,DIMY);
 
  Xb = add_buffer(X,DIMX,DIMY,XBUF,YBUF);
  printf("Xb:\n");
  printrm_rowmaj(Xb,DIMX+2*XBUF, DIMY+2*YBUF);
  return 0;
}


