#ifndef STELLA_MAT_H
#define STELLA_MAT_H

#include "petscmat.h"

#ifdef WITH_CEDAR
#include <cedar/capi.h>

typedef struct
{
	int nd;
	cedar_mat so;
	cedar_vec solvec;
	cedar_vec rhsvec;
} stella_bmg_mat;
#endif

PetscErrorCode stella_bmg_SetValuesStencil(Mat mat, PetscInt m, const MatStencil idxm[], PetscInt n,
                                           const MatStencil idxn[], const PetscScalar v[],
                                           InsertMode addv);


PetscErrorCode stella_bmg_mult(Mat mat, Vec, Vec);

#endif
