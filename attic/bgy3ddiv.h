/* -*- mode: c; c-basic-offset: 2; -*- vim: set sw=2 tw=70 et sta ai: */
/* Used in bgy3d-test.c: */
typedef struct BGY3dDivStruct
{
  DA da;
  Vec f[3];
  Vec ddU;
  Vec boundary;
  Vec v[3], i[3], v2[3];
  Mat M;
  Mat FD[3];

  real LJ_params[2];            /* sigma and epsilon  */
  real beta, rho;

  Vec g_ini, g_SA;

  /* Parallel FFT */
  struct fft_plan_3d *fft_plan;
  FFT_DATA *(fg2_fft[3]), *g_fft, *gfg2_fft;

  const ProblemData *PD;
} *BGY3dDivData;

BGY3dDivData BGY3dDivData_malloc(const ProblemData *PD, PetscTruth flg);
BGY3dDivData BGY3dDivData_kirk_malloc(const ProblemData *PD, PetscTruth flg);
void BGY3dDivData_free(BGY3dDivData BDD);
void AssembleMatrix(BGY3dDivData BDD, DA da, Mat M);
Vec BGY3dDiv_solve(const ProblemData *PD, Vec g_ini);
Vec BGY3dDiv_solve2(const ProblemData *PD, Vec g_ini);
void AssembleFDMatrix(BGY3dDivData BDD, DA da, Mat M, int vdim);
void ComputeIntegralPart(BGY3dDivData BDD, Vec g, Vec f);
void ComputeIntegralPart_kirk(BGY3dDivData BDD, Vec g, Vec f);
void AssembleSystemMatrix(BGY3dDivData BDD, Mat SM, Vec f);
void AssembleSystemMatrix_part2(BGY3dDivData BDD, Mat SM);
void ComputeRHS(BGY3dDivData BDD, Vec b, Vec g0, Vec f);
void ComputeBGY3dDiv_F(BGY3dDivData BDD, Mat SM, Vec g0, Vec dg, Vec g,
		       Vec b, Vec f);
void ComputeRHS2(BGY3dDivData BDD, Vec b);
void ShiftVec(DA da, Vec g, Vec scratch, const int N[3]);
void AssembleSystemMatrix_part2b(BGY3dDivData BDD, Mat SM);
