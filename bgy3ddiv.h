BGY3dDivData BGY3dDivData_malloc(ProblemData *PD, PetscTruth flg);
BGY3dDivData BGY3dDivData_kirk_malloc(ProblemData *PD, PetscTruth flg);
void BGY3dDivData_free(BGY3dDivData BDD);
void AssembleMatrix(BGY3dDivData BDD, DA da, Mat M);
Vec BGY3dDiv_solve(ProblemData *PD, Vec g_ini, int vdim);
Vec BGY3dDiv_solve2(ProblemData *PD, Vec g_ini, int vdim);
void AssembleFDMatrix(BGY3dDivData BDD, DA da, Mat M, int vdim);
void ComputeIntegralPart(BGY3dDivData BDD, Vec g, Vec f);
void ComputeIntegralPart_kirk(BGY3dDivData BDD, Vec g, Vec f);
void AssembleSystemMatrix(BGY3dDivData BDD, Mat SM, Vec f);
void AssembleSystemMatrix_part2(BGY3dDivData BDD, Mat SM);
void ComputeRHS(BGY3dDivData BDD, Vec b, Vec g0, Vec f);
void ComputeBGY3dDiv_F(BGY3dDivData BDD, Mat SM, Vec g0, Vec dg, Vec g,
		       Vec b, Vec f);
void ComputeRHS2(BGY3dDivData BDD, Vec b);
void ShiftVec(DA da, Vec g, Vec scratch, int N[3]);
void AssembleSystemMatrix_part2b(BGY3dDivData BDD, Mat SM);
