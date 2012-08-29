BGY3dFourierData BGY3dFourierData_malloc(ProblemData *PD);
BGY3dFourierData BGY3dFourierData_kirk_malloc(ProblemData *PD);
void BGY3dFourierData_free(BGY3dFourierData BDD);
void ExtractAxis(BGY3dFourierData BDD, Vec g, int axis);
Vec BGY3dDiv_solve_Fourier(ProblemData *PD, Vec g_ini, int vdim);
