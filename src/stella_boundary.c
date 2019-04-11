#include <stdlib.h>

#include "stella_boundary.h"


PetscErrorCode stella_boundary_create(stella_boundary **stella_bnd, stella_level *level,
                                  stella_dmap *dmap, stella_state *state, stella_fd *fd)
{
	stella_boundary *bnd = (stella_boundary*) malloc(sizeof(stella_boundary));

	bnd->bc = NULL;
	bnd->level = level;
	bnd->slv_dmap = dmap;
	bnd->state = state;
	bnd->fd = fd;
	bnd->num_bc = 0;
	bnd->axisymmetric = 0;

	*stella_bnd = bnd;
	return 0;
}


PetscErrorCode stella_boundary_apply(stella_boundary *bnd, Mat A, DM da)
{
	PetscErrorCode ierr;
	int i;

	ierr = DMGlobalToLocalBegin(da, bnd->level->add_cont, INSERT_VALUES, bnd->level->ladd_cont);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(da, bnd->level->add_cont, INSERT_VALUES, bnd->level->ladd_cont);CHKERRQ(ierr);

	for (i = 0; i < bnd->num_bc; i++) {
		if (bnd->bc[i]->btype == STELLA_NEUMANN) {
			ierr = stella_bc_apply(bnd->bc[i], A, da);CHKERRQ(ierr);
		}
	}
	for (i = 0; i < bnd->num_bc; i++) {
		if (bnd->bc[i]->btype == STELLA_DIRICHLET ||
		    bnd->bc[i]->btype == STELLA_SCHWARZ) {
			ierr = stella_bc_apply(bnd->bc[i], A, da);CHKERRQ(ierr);
		}
	}

	for (i = 0; i < bnd->num_bc; i++) {
		ierr = stella_bc_symmetric(bnd->bc[i], A, da);CHKERRQ(ierr);
	}

	ierr = DMLocalToGlobalBegin(da, bnd->level->ladd_cont, ADD_VALUES, bnd->level->add_cont);CHKERRQ(ierr);
	ierr = DMLocalToGlobalEnd(da, bnd->level->ladd_cont, ADD_VALUES, bnd->level->add_cont);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_boundary_apply_rhs(stella_boundary *bnd, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	int i;

	for (i = 0; i < bnd->num_bc; i++) {
		ierr = stella_bc_apply_rhs(bnd->bc[i], da, rhs);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_boundary_set(stella_boundary *bnd, char *norm_dir, double *values)
{
	PetscErrorCode ierr;

	stella_bctype bcs[2];
	bcs[0] = STELLA_DIRICHLET;
	bcs[1] = STELLA_NEUMANN;

	bnd->num_bc = 2;
	bnd->bc = malloc(bnd->num_bc*sizeof(stella_bc*));

	int i;
	for (i = 0; i < bnd->num_bc; i++) {
		ierr = stella_bc_create(&bnd->bc[i], bcs[i], norm_dir, values,
		                        bnd->level, bnd->slv_dmap, bnd->state, bnd->fd);CHKERRQ(ierr);
		bnd->bc[i]->axisymmetric = bnd->axisymmetric;
	}

	return 0;
}


PetscErrorCode stella_boundary_destroy(stella_boundary *bnd)
{
	PetscErrorCode ierr;
	int i;

	for (i = 0; i < bnd->num_bc; i++) {
		ierr = stella_bc_destroy(bnd->bc[i]);CHKERRQ(ierr);
	}

	free(bnd->bc);

	return 0;
}
