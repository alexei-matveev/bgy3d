Vec BGY3dDiv_test(ProblemData *PD, Vec g_ini, int vdim);
void InitializeTestData(BGY3dDivData BDD, Vec g, real sigma_g, real sigma_K);
void ComputeRHStest(BGY3dDivData BDD, Vec g, Vec rhs, real sigma_g,
		    real sigma_K);
Vec BGY3dDiv_solve_FourierTest(ProblemData *PD, Vec g_ini, int vdim);
Vec BGY3d_Convolution_Test(ProblemData *PD, Vec g_ini, int vdim);
