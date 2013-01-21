void bgy3d_fft_mat_create (const int N[3], Mat *A, DA *da, DA *dc);

void bgy3d_fft_interp (const Mat A,
                       const Vec Y, /* complex, intent(in) */
                       int np, double x[np][3], /* intent(in) */
                       double y[np]);           /* intent(out) */

