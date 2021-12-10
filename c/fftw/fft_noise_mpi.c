/*-------------------------------------------------------------------------------- 
Toy program to test parallelize FFTW3's 2d DFT, real -> complex.  
Must compile with linker flags gcc test1d.c -lfftw3 -lm 
Will read in file "real_noise.txt" containing 256 * 256 doubles 
representing a 2D grid and execute DFT.  Output is stored in array 
cmplx_out.  The file "complex_noise.txt" contains the expected
output.

Change DIMX, DIMY arrays to test smaller grids.  Uncomment lines
to print input/output.

Notes:
See
https://www.fftw.org/fftw3_doc/Multi_002dthreaded-FFTW.html

-------------------------------------------------------------------------------*/

#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fftw3-mpi.h>

#define DIMX 256
#define DIMY 256

char *cfs(fftw_complex c);
char buff[100];
int main(int argc, char *argv[]) 
{
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fftw_mpi_init();
  FILE *real_in, *fout;
  fftw_plan p;
  fftw_complex *cmplx_out, *noise, *data;
  double fp;
  int i,j;
  ptrdiff_t alloc_local, local_n0, local_0_start
 


  if (rank == 0){
    real_in = fopen("real_noise.txt", "r"); 
    for (i=0; i < DIMX * DIMY; i++) {
      fscanf(real_in, "%lf", &fp);
      data[i] = fp;
      fclose(real_in);
    }
  }

  MPI_Bcast(data, DIMX * DIMY, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

  alloc_local = fftw_mpi_local_size_2d(DIMX, DIMY, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
  noise = fftw_alloc_complex(alloc_local);
  MPI_Barrier(MPI_COMM_WORLD);
  for (i = 0; i < local_n0; ++i) for (j = 0; j < DIMY; ++j)
        noise[i*DIMY + j] = data[i*DIMY + j];
  

//  for (i=0; i < DIMX * DIMY; i++)
//    printf("%f ", creal(noise[i]))

  p = fftw_mpi_plan_dft_2d(DIMX, DIMY, noise, cmplx_out, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  MPI_Finalize();

  
  fout = fopen("complex_noise.txt", "w");

//  for (i=0; i < DIMX * DIMY ; i++)
//    printf("%s ", cfs(cmplx_out[i]))

  fclose(fout);
  return 0;

}

/* complex format string for printing complex numbers*/
char *cfs(fftw_complex c) {
  int n;
  n = snprintf(buff, 100, "%f + %fi",creal(c), cimag(c));
  return buff;
}
