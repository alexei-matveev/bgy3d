/* Used in bgy3dtest.c: */
typedef struct BGY3dFourierStruct
{
  DA da;
  Vec f[3];
  Vec v[3];

  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;

  Vec g_ini;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;

  ProblemData *PD;
} *BGY3dFourierData;

BGY3dFourierData BGY3dFourierData_malloc(ProblemData *PD);
BGY3dFourierData BGY3dFourierData_kirk_malloc(ProblemData *PD);
void BGY3dFourierData_free(BGY3dFourierData BDD);
void ExtractAxis(BGY3dFourierData BDD, Vec g, int axis);
Vec BGY3dDiv_solve_Fourier(ProblemData *PD, Vec g_ini);
