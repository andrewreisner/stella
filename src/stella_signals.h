#ifndef STELLA_SIGNALS_H
#define STELLA_SIGNALS_H

#include "stella_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Signal boundary conditions have changed
 */
PetscErrorCode stella_changed_bc(stella * slv);


/**
 * Signal rhs has changed
 */
PetscErrorCode stella_changed_rhs(stella *slv);

#ifdef __cplusplus
}
#endif

#endif
