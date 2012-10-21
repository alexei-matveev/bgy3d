/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 et sta ai: */
void bgy3d_fft_init_da (const int N[3],
                        fftwnd_mpi_plan *fw, fftwnd_mpi_plan *bw,
                        DA *da, DA *da_mg);

fftw_complex *bgy3d_fft_malloc (DA da);
void bgy3d_fft_free (fftw_complex *ptr);

fftw_complex *bgy3d_fft_axpby (DA da, fftw_complex *restrict y,
                               double alpha, double beta,
                               const fftw_complex *x);

fftw_complex *bgy3d_fft_set (DA da, fftw_complex *y, double alpha);


void ComputeFFTfromVec_fftw (Mat fft_mat, Vec g, fftw_complex *g_fft);
void ComputeVecfromFFT_fftw (Mat fft_mat, Vec g, fftw_complex *g_fft);

FFT_DATA *ComputeFFTfromVec(DA da, struct fft_plan_3d *fft_plan, Vec g,
			    FFT_DATA *g_fft);
void ComputeVecfromFFT(DA da, struct fft_plan_3d *fft_plan, Vec g,
		       FFT_DATA *g_fft);

double bgy3d_fft_test (int m, int n, int k);
