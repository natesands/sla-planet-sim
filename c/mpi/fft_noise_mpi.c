/*-------------------------------------------------------------------------------- 
Toy program to test parallelize FFTW3's 2d DFT, real -> complex.  
Must compile with linker flags gcc test1d.c -lfftw3 -lm 
Will read in file "real_noise.txt" containing 256 * 256 doubles 
representing a 2D grid and execute DFT.  Output is stored in array 
cmplx_out.  The file "complex_noise.txt" contains the expected
output.

Change N0, N1 arrays to test smaller grids.  Uncomment lines
to print input/output.

Notes:
See
https://www.fftw.org/fftw3_doc/Multi_002dthreaded-FFTW.html

% mpicc -o fft_noise_mpi fft_noise_mpi.c -lmpi -lfftw3 -lfftw3_mpi -lm
% mpirun -np 2 fft_noise_mpi
-------------------------------------------------------------------------------*/

#include <complex.h>
#include <fftw3.h>
#include <stdio.h>
#include <math.h>
#include <fftw3-mpi.h>
#include <stdlib.h>

#define N0 4
#define N1 4

char *cfs(fftw_complex c);
char buff[100];
int main(int argc, char **argv) 
{

  fftw_plan plan;
  fftw_complex *mydata, *noise;
  double fp;
  int myid, manager, numprocs;
  ptrdiff_t i,j,  local_n0, local_0_start, alloc_local;
  manager = 0;
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  fftw_mpi_init();
  /* manager */
  if (myid == manager ) {
    int worker_rows, worker_0_start, worker; 
    fftw_complex *buf_ptr, buf[N1];
    fftw_complex *cmplx_in, *cmplx_out;    /* array to hold final dft */
    FILE *noise;                /* input file ptr */
    FILE *fout;                 /* file to write output of dft */
    cmplx_in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)* N0 * N1);
    cmplx_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N0 * N1);
    noise = fopen("noise_small.txt", "r");   /* read noise from file */
    for (i=0; i < N0 * N1; i++) {
      fscanf(noise, "%lf", &fp);
      cmplx_in[i] = (fftw_complex) fp;
    }

    fclose(noise);
    
    for (i=0; i < N0*N1; i++)
      printf("%.3f+%.3fi\n", creal(cmplx_in[i]), cimag(cmplx_in[i]));

    alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD, &local_n0, 
                                         &local_0_start);

    /* reserve space for manager's share of work and copy data */
    mydata = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * local_n0 * N1);
    for (i = 0; i < local_n0 * N1; i++)
      mydata[i] = cmplx_in[local_0_start*N0 + i];

    printf("manager's data:\n");
    for (i=0; i < local_n0 * N1; i++)
      printf("%.3f+%.3fi\n", creal(mydata[i]), cimag(mydata[i]));

    /* receive requests for work from workers and send them rows */
    for (i = 1; i < numprocs; i++) {
      MPI_Recv(&worker_rows, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      worker = status.MPI_SOURCE;
      worker_0_start = status.MPI_TAG;
      printf("worker %d requested %d rows starting at index %d\n", worker, 
          worker_rows, worker_0_start);
      //buf_ptr = (fftw_complex *) malloc(sizeof(fftw_complex)*worker_rows*N1);
      for (j=0; j < worker_rows * N1; j++) {
        buf[j] = cmplx_in[worker_0_start*N0+j];
        //buf_ptr[i] = cmplx_in[worker_0_start*N0 + i];
      }
      printf("sending worker %d:\n", worker);
      for (j=0; j < worker_rows * N1; j++)
        printf("%.3f+%.3fi\n", creal(buf[j]), cimag(buf[j]));
      MPI_Send(buf, worker_rows*N1, MPI_C_DOUBLE_COMPLEX, worker, worker,
               MPI_COMM_WORLD);
      
    //  for (i=0; i < 5; i++)
    //    printf("%.1f+%.1fi\n", creal(buf_ptr[i]), cimag(buf_ptr[i]));
      // fftw_free(buf_ptr);
    }
  }
  /* worker */
  else {
    mydata = (fftw_complex *) malloc(sizeof(fftw_complex)*local_n0*N1);
    alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD, &local_n0, 
                                         &local_0_start);
    printf("myid = %d, local_n0 = %li, local_0_start = %li\n", myid, local_n0, local_0_start);
    /* request rows from manager */
    MPI_Send(&local_n0, 1, MPI_INT, manager, local_0_start, MPI_COMM_WORLD);
    /* receive rows */
    MPI_Recv(mydata, local_n0 * N1, MPI_C_DOUBLE_COMPLEX, manager, MPI_ANY_TAG, 
                MPI_COMM_WORLD, &status);
    printf("worker %d received %zu elements from manager\n", myid, status._ucount);
    for (j=0; j < local_n0*N1; j++)
      printf("%.3f+%.3fi\n", creal(mydata[j]), cimag(mydata[j]));
  }
  fftw_destroy_plan(plan);
  MPI_Finalize();
  return 0;
}



  /* broadcast 

  alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD, &local_n0,
                                       &local_0_start);
  plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
                              FFTW_FORWARD, FFTW_ESTIMATE);

  printf("myid = %d, local_n0 = %li, local_0_start = %li\n", myid, local_n0, local_0_start);

//  p = fftw_plan_dft_2d(N0, N1, noise, cmplx_out, FFTW_FORWARD, FFTW_ESTIMATE);
//  fftw_execute(p);
//
//
//  fout = fopen("complex_noise.txt", "w");
//
//  for (i=0; i < N0 * N1 ; i++)
//    fprintf(fout, "%s ", cfs(cmplx_out[i]));
//
//  fclose(fout);
*/  


/* complex format string for printing complex numbers*/
char *cfs(fftw_complex c) {
  int n;
  n = snprintf(buff, 100, "%f + %fi",creal(c), cimag(c));
  return buff;
}
