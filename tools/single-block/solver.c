#include <stdlib.h>

#include "solver.h"
#include "problem.h"
#include "grid.h"


solver *solver_create(grid *grd, problem *pb)
{
	solver *slv = (solver*) malloc(sizeof(solver));
	slv->bnd = boundary_create(grd, pb);
	slv->state = state_create(grd, pb);
	slv->axisymmetric = 0;

	return slv;
}


PetscErrorCode solver_init(solver *slv, grid *grd)
{
	PetscErrorCode ierr;
	int i;
	int stride[3], offset[3];

	blist *bnd = slv->bnd;

	for (i = 0; i < 3; i++) offset[i] = 0;
	stride[0] = grd->nx; stride[1] = grd->ny; stride[2] = grd->nz;

	ierr = PetscInitialize(NULL, NULL, "PETScOptions.txt", NULL);CHKERRQ(ierr);
	int overlap_periodic = 1;
	ierr = stella_init(&slv->ptr, grd->comm, grd->num_global, grd->num_procs,
	                   grd->num_local,
	                   offset, stride,
	                   grd->cart_coord, grd->periodic, overlap_periodic, grd->nd,
	                   slv->axisymmetric);CHKERRQ(ierr);
	ierr = stella_set_grid(slv->ptr, grd->is, grd->ie, grd->num_pts, grd->xyz);CHKERRQ(ierr);
	ierr = stella_set_external(slv->ptr, slv->state->phi, slv->state->eps,
	                           slv->state->debye, slv->state->jc);CHKERRQ(ierr);
	ierr = stella_set_sol(slv->ptr, slv->state->sol);CHKERRQ(ierr);

	for(i = 0; i < bnd->len; i++) {
		ierr = stella_boundary_add(slv->ptr->boundary, (stella_bctype) bnd->v[i].type, bnd->v[i].norm_dir,
		                           bnd->v[i].is, bnd->v[i].ie, bnd->v[i].dirichlet);CHKERRQ(ierr);
	}

	ierr = stella_setup_op(slv->ptr);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode solver_run(solver *slv)
{
	PetscErrorCode ierr;

	ierr = stella_set_rhs(slv->ptr, slv->state->rhs);CHKERRQ(ierr);
	ierr = stella_solve(slv->ptr);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode solver_destroy(solver *slv)
{
	PetscErrorCode ierr;

	ierr = stella_cleanup(slv->ptr);CHKERRQ(ierr);
	boundary_destroy(slv->bnd);
	free(slv->bnd);
	state_destroy(slv->state);
	free(slv->state);

	ierr = PetscFinalize();CHKERRQ(ierr);

	return 0;
}
