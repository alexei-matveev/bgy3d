BGY3dParameterVec BGY3dParameterVec_malloc(ProblemData *PD);
void BGY3dParameterVec_free(BGY3dParameterVec par_vec);

/* Vec */
Vec BGY3d_vec_solve(ProblemData *PD, Vec g_ini, int vdim);
void CreateInitialGuess_vec(BGY3dParameterVec par_vec, Vec g);
PetscErrorCode ComputeVec_F(SNES snes, Vec g, Vec f, void *pa);
