/*-------------------------------------------------------------------------------- 
Toy program to test parallelize FFTW3's 2d DFT, real -> complex.  
Must compile with linker flags gcc test1d.c -lfftw3 -lm 
Will read required doubles from file "noise.txt" containing 256 * 256 doubles.

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
#include <fftw3-mpi.h>
#include <stdlib.h>
#include <fftw3-mpi.h>

/* grid dimensions */

#define N0 256
#define N1 256

fftw_complex dft_in[N0 * N1];
fftw_complex dft_out[N0 * N1];

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
  for (i = 0; i < local_n0 * dimy; i++)
    arr[i] = dft_in[local_0_start * dimy + i];
}


int main(int argc, char **argv)
{
    fftw_plan plan;
    fftw_complex *data;
    int myid, numprocs, manager = 0, worker_done, begin, rows_completed;
    ptrdiff_t alloc_local, local_n0, local_0_start, i, j;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Status status; 
    fftw_mpi_init();

    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
                                         &local_n0, &local_0_start);
    data = fftw_alloc_complex(alloc_local);

    /* create plan for in-place forward DFT */
    plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
                                FFTW_FORWARD, FFTW_ESTIMATE);    

    /* manager loads data from file */ 
    if (myid == manager)  {
      load_data("noise.txt", N0, N1);
      begin = 1;
      MPI_Bcast(&begin, 1, MPI_INT, manager, MPI_COMM_WORLD);
    }
    /* workers receive broadcast to begin */
    else {
      MPI_Bcast(&begin, 1, MPI_INT, manager, MPI_COMM_WORLD);
    }

    /* grab share of data */
    get_data(data, N0, N1, local_0_start, local_n0);

    /* compute transforms in-place */
    fftw_execute(plan);

    /* write data to output array */
    for (i = 0; i < local_n0 * N1; i++)
      dft_out[local_0_start * N1 + i] = data[i];
  
    /* Manager polls workers */
    if (myid == manager ) {
      rows_completed = local_n0;  /* rows computed by manager */;
      printf("Manager computed rows 0 through %d\n", rows_completed-1);
      for (i = 1; i < numprocs; i++) {
        MPI_Recv(&worker_done, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        printf("Worker %d computed rows %d through %d\n", status.MPI_SOURCE,
                rows_completed, rows_completed + status.MPI_TAG -1);
        rows_completed += status.MPI_TAG;
      }
    }
    else {
      worker_done = 1;
      MPI_Send(&worker_done, 1, MPI_INT, manager, local_n0, MPI_COMM_WORLD);
    }
    
//    for (i=0; i < local_n0 * N1; i++)
//      printf("%.3f+%.3fi\n", creal(data[i]), cimag(data[i]));
//
    fftw_destroy_plan(plan);

    MPI_Finalize();
}
