fftw_complex *bgy3d_fft_malloc (DA da);
void bgy3d_fft_free (fftw_complex *ptr);

fftw_complex *bgy3d_fft_axpby (DA da, fftw_complex *restrict y,
                               double alpha, double beta,
                               const fftw_complex *x);

fftw_complex *bgy3d_fft_set (DA da, fftw_complex *y, double alpha);


fftw_complex *ComputeFFTfromVec_fftw(DA da, fftwnd_mpi_plan fft_plan, Vec g,
				fftw_complex *g_fft, fftw_complex *work);
void ComputeVecfromFFT_fftw(DA da, fftwnd_mpi_plan fft_plan, Vec g,
			    fftw_complex *g_fft, fftw_complex *work);
