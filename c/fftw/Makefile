fftw: 
	gcc -I${FFTW_ROOT}/include fft_noise.c -o fft_noise -L${FFTW_ROOT}/lib -lfftw3 -lm

fftw-mpi:
	mpicc -I${FFTW_ROOT}/include fft_noise_mpi.c -o fft_noise_mpi -L${FFTW_ROOT}/lib -lfftw3_mpi -lfftw3 -lm

test-mpi:
	mpicc -I${FFTW_ROOT}/include test.c -o test -L${FFTW_ROOT}/lib -lfftw3_mpi -lfftw3 -lm