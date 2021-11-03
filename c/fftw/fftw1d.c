
/*-------------------------------------------------------------------------------- 
Test FFTW3's 1d DFT.  Must compile with linker flags gcc test1d.c -lfftw3 -lm 
Notes:
-in and out arrays can be the same for an in-place transform
-if program performs many transforms of the same size and initialization time is not important, use FFTW_MEASURE
-must create the plan before initializing the input
-The DFT results are stored in-order in the array out, with the zero-frequency (DC) component in out[0]
-if <complex.h> is included before <fftw3.h> then fftw_complex is the native double-precision complex type 
and you can manipulate it with ordinary arithmetic
-------------------------------------------------------------------------------*/
#include <fftw3.h>
#include <stdio.h>

#define c_re(c)   (c[0])
#define c_im(c)   (c[1])

char buff[100];
/* print complex */
void pc(fftw_complex c) {
  printf("%f + %fi", c_re(c), c_im(c));
}
/* complex format string for printing complex numbers*/
char *cfs(fftw_complex c) {
  int n;
  n = snprintf(buff, 100, "%f + %fi",c_re(c), c_im(c));
  return buff;
}

int main() {
  fftw_complex in[] = {{3.0,0.0}, {2.0,0}};
  fftw_complex out[2];
  fftw_plan p;
  int i;

  p = fftw_plan_dft_1d(2, in, out,FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_execute(p);
  printf("IN:\n");
  for (i=0; i<2; i++)
    printf("%s ", cfs(in[i]));
  printf("\nOUT:\n");
  for (i=0; i<2; i++)
    printf("%s ", cfs(out[i]));
  printf("\n");

  fftw_destroy_plan(p);
  
  return 0;
}

