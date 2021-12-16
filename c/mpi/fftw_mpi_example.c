/*-------------------------------------------------------------------------------- 
Program to test parallelization ofFFTW3's 2d DFT, real -> complex.  

Will read required doubles from file "noise.txt" containing 256 * 256 doubles.

Change N0, N1 arrays to test smaller grids.

Notes:
See
https://www.fftw.org/fftw3_doc/Multi_002dthreaded-FFTW.html

% mpicc -o fft_noise_mpi fft_noise_mpi.c -lmpi -lfftw3 -lfftw3_mpi -lm
% mpirun -np 2 fft_noise_mpi
-------------------------------------------------------------------------------*/

#include <complex.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <stdlib.h>
#include <fftw3-mpi.h>
#include <unistd.h>
/* grid dimensions */

#define N0 256
#define N1 256

fftw_complex dft_in[N0 * N1];
fftw_complex dft_out[N0 * N1];
char buff[100];

/* complex format string for printing complex numbers*/
char *cfs(fftw_complex c) {
  int n;
  n = snprintf(buff, 100, "%.3f + %.3fi",creal(c), cimag(c));
  return buff;
}

/* loads data from file into array dft_in */
void load_data(char *file_name, int dimx, int dimy) {
  FILE* fin = fopen(file_name, "r");
  double fp;
  int i;
  for (i = 0; i < dimx * dimy; i++) {
    fscanf(fin, "%lf", &fp);
    dft_in[i] = (fftw_complex) fp;
  } 
  fclose(fin);
}

/* copies each worker's alloted data into their array */
void get_data(fftw_complex *arr, ptrdiff_t dimx, ptrdiff_t dimy,  
              ptrdiff_t local_0_start, ptrdiff_t local_n0) {
  int i;
  printf("start : %td, rows: %td\n", local_0_start, local_n0);
  for (i = 0; i < local_n0 * dimy; i++) {
    arr[i] = dft_in[local_0_start * dimy + i];
  }
}

/* write transformed data to file */
void write_output(char *file_name, int dimx, int dimy) {    
  FILE* fout = fopen(file_name, "w");
  for (int i=0; i < N0 * N1; i++)
    fprintf(fout, "%s ", cfs(dft_out[i]));
  fclose(fout);
}

int main(int argc, char **argv)
{
    fftw_plan plan;
    fftw_complex *data;
    int myid, numprocs, manager = 0, worker_done, begin = 0, rows_completed, chunk_size;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Status status; 
    fftw_mpi_init();

    /* fix this... requires N0*N1 be divisible by numprocs */
    chunk_size = (N0 * N1) / numprocs;

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    data = fftw_alloc_complex(alloc_local);

    /* create plan for in-place forward DFT */
    plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
                                FFTW_FORWARD, FFTW_ESTIMATE);    

    /* Manager */
    if (myid == manager)  {

      /* load data from file into array dft_in*/ 
      fftw_complex *cbuf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * chunk_size);
      load_data("noise.txt", N0, N1);
  //    printf("data in: \n");
  //    for (i=0; i < N0 * N1; i++)
  //      printf("%.3f+%.3fi\n", creal(dft_in[i]), cimag(dft_in[i]));
      
      /* copy first chunk to data array */
      for (i = 0; i < local_n0 * N1; i++) 
        data[i] = dft_in[local_0_start * N1 + i];
      
      /* send data to workers */  
      for (i = 1; i < numprocs; i++) {
        for(int j = 0; j < chunk_size; j++) 
          cbuf[j] = dft_in[i*chunk_size + j];
        MPI_Send(cbuf, chunk_size, MPI_C_DOUBLE_COMPLEX, i, i, MPI_COMM_WORLD);
      }
      /* compute in-place transform */
      fftw_execute(plan);
      printf("proc %d completed transform of %d rows\n", myid, chunk_size / N1);
      /* copy to output array dft_out */
      for (int i=0; i < chunk_size; i++)
        dft_out[local_0_start*N1 + i] = data[i];

      /* receive transformed data from workers and copy to output array */
      for (i = 1; i < numprocs; i++) {
        MPI_Recv(cbuf, chunk_size, MPI_C_DOUBLE_COMPLEX, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        for (j = 0; j < chunk_size; j++)
          dft_out[i*chunk_size + j] = cbuf[j];
      } 

      /* write output array to file */
      write_output("noise_out.txt", N0, N1);    

      fftw_free(cbuf);
    }

    /* workers receive data from manager */
    else {
      MPI_Recv(data, chunk_size, MPI_C_DOUBLE_COMPLEX, manager, MPI_ANY_TAG,
               MPI_COMM_WORLD, &status); 

      /* compute transform in-place */
      fftw_execute(plan);
      printf("proc %d completed transform of %d rows\n", myid, chunk_size / N1);

      /* send manager transformed data */
      MPI_Send(data, chunk_size, MPI_C_DOUBLE_COMPLEX, manager, myid, MPI_COMM_WORLD);
    }

    fftw_destroy_plan(plan);
    MPI_Finalize();
}
