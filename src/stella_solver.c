#include <stdlib.h>

#include "stella_solver.h"
#include "stella_io.h"
#include "stella_signals.h"
#include "stella_util.h"

#include "stella_mat.h"
#include "stella_pc.h"
#include "stella_classify.h"

#ifdef WITH_CEDAR
#include <cedar/capi.h>
#include "stella_cedar.h"
#endif

#include "petscmat.h"

static inline PetscErrorCode petsc_options_has_name(const char name[], PetscBool *set) {

	PetscErrorCode ierr;

	#ifdef PETSC_OPTIONS_TAKES_DATABASE
	ierr = PetscOptionsHasName(PETSC_NULL, PETSC_NULL, name, set);CHKERRQ(ierr);
	#else
	ierr = PetscOptionsHasName(PETSC_NULL, name, set);CHKERRQ(ierr);
	#endif

	return 0;
}

static inline PetscErrorCode petsc_options_set_value(const char name[], const char value[]) {

	PetscErrorCode ierr;

	#ifdef PETSC_OPTIONS_TAKES_DATABASE
	ierr = PetscOptionsSetValue(PETSC_NULL, name, value);CHKERRQ(ierr);
	#else
	ierr = PetscOptionsSetValue(name, value);CHKERRQ(ierr);
	#endif

	return 0;
}


/**
 * Fill external array with updated solution
 */
static PetscErrorCode update_solution(stella *slv)
{
	PetscErrorCode ierr;
	PetscInt i, j, k, xs, ys, zs, xm, ym, zm;
	PetscInt ngx,ngy,ngz;
	DMBoundaryType bx, by, bz;

	ierr = DMDAGetInfo(slv->dm, 0, &ngx, &ngy, &ngz,0,0,0,0,0, &bx, &by, &bz,0);CHKERRQ(ierr);
	ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	if (slv->grid.nd == 2) {
		Vec loc_x;
		PetscScalar **xvec;
		double **phi;

		ierr = stella_dmap_get(slv->dmap, slv->state.phi, &phi);CHKERRQ(ierr);
		ierr = DMGetLocalVector(slv->dm, &loc_x);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(slv->dm, slv->x, INSERT_VALUES, loc_x);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(slv->dm, slv->x, INSERT_VALUES, loc_x);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, loc_x, &xvec);CHKERRQ(ierr);

		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				phi[j][i] = xvec[j][i];
			}
		}

		if ((by == DM_BOUNDARY_PERIODIC) && (ys + ym == ngy) && slv->grid.overlap_periodic) {
			for (i = xs; i < xs + xm; i++) {
				phi[ngy][i] = xvec[ngy][i];
			}
		}

		if ((bx == DM_BOUNDARY_PERIODIC) && (xs + xm == ngx) && slv->grid.overlap_periodic) {
			for (j = ys; j < ys + ym; j++) {
				phi[j][ngx] = xvec[j][ngx];
			}
		}

		ierr = stella_dmap_restore(slv->dmap, &phi);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, loc_x, &xvec);CHKERRQ(ierr);
		if (stella_log(slv, STELLA_LOG_VTK)) {
			ierr = DMLocalToGlobalBegin(slv->dm, loc_x, INSERT_VALUES, slv->x);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(slv->dm, loc_x, INSERT_VALUES, slv->x);CHKERRQ(ierr);
			ierr = stella_io_vtkwrite(slv->dm, slv->x, "phi", slv->grid.id, slv->ts);CHKERRQ(ierr);
		}
		ierr = DMRestoreLocalVector(slv->dm, &loc_x);CHKERRQ(ierr);
	} else if (slv->grid.nd == 3) {
		PetscScalar ***xvec;
		double ***phi;
		Vec loc_x;

		ierr = stella_dmap_get(slv->dmap, slv->state.phi, &phi);CHKERRQ(ierr);
		ierr = DMGetLocalVector(slv->dm, &loc_x);CHKERRQ(ierr);
		ierr = DMGlobalToLocalBegin(slv->dm, slv->x, INSERT_VALUES, loc_x);CHKERRQ(ierr);
		ierr = DMGlobalToLocalEnd(slv->dm, slv->x, INSERT_VALUES, loc_x);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, loc_x, &xvec);CHKERRQ(ierr);

		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					phi[k][j][i] = xvec[k][j][i];
				}
			}
		}

		if ((bz == DM_BOUNDARY_PERIODIC) && (zs + zm == ngz) && slv->grid.overlap_periodic) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					phi[ngz][j][i] = xvec[ngz][j][i];
				}
			}
		}
		if ((by == DM_BOUNDARY_PERIODIC) && (ys + ym == ngy) && slv->grid.overlap_periodic) {
			for (k = zs; k < zs + zm; k++) {
				for (i = xs; i < xs + xm; i++) {
					phi[k][ngy][i] = xvec[k][ngy][i];
				}
			}
		}
		if ((bx == DM_BOUNDARY_PERIODIC) && (xs + xm == ngx) && slv->grid.overlap_periodic) {
			for (k = zs; k < zs + zm; k++) {
				for (j = ys; j < ys + ym; j++) {
					phi[k][j][ngx] = xvec[k][j][ngx];
				}
			}
		}

		ierr = stella_dmap_restore(slv->dmap, &phi);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, loc_x, &xvec);CHKERRQ(ierr);
		if (stella_log(slv, STELLA_LOG_VTK)) {
			ierr = DMLocalToGlobalBegin(slv->dm, loc_x, INSERT_VALUES, slv->x);CHKERRQ(ierr);
			ierr = DMLocalToGlobalEnd(slv->dm, loc_x, INSERT_VALUES, slv->x);CHKERRQ(ierr);
			ierr = stella_io_vtkwrite(slv->dm, slv->x, "phi", slv->grid.id, slv->ts);CHKERRQ(ierr);
		}
		ierr = DMRestoreLocalVector(slv->dm, &loc_x);CHKERRQ(ierr);
	}

	return 0;
}

#ifdef WITH_CEDAR
/**
 * Sets a custom cedar config file name. If not set, stella defaults to config.json.
 */
PetscErrorCode stella_set_cedar_config(stella *slv, char fname[]) {
	slv->cedar_config_filename = fname;
}
#endif


/**
 * Creates and initializes data needed for setting up elliptic discretization.
 */
static PetscErrorCode stella_setup(stella *slv, int offset[], int stride[])
{
	PetscErrorCode ierr;
	PetscInt xs, ys, zs, xm, ym, zm, ngx, ngy, ngz;
	stella_pc *shell;
	PCType pc_type;

	if (stella_log(slv, STELLA_LOG_STATUS)) {
		ierr = stella_io_print(slv->comm, "Setting up electric field solver");CHKERRQ(ierr);
	}

	slv->ts = 0;

	if (stella_log(slv, STELLA_LOG_RESIDUAL)) {
		ierr = petsc_options_set_value("-ksp_monitor_short", NULL);CHKERRQ(ierr);
	}

	ierr = DMDASetFieldName(slv->dm, 0,"potential");CHKERRQ(ierr);

	if (slv->grid.nd == 2) {
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);
		slv->dmap = stella_dmap_create_2d(xs - offset[0], ys - offset[1], xm, ym, stride);
	} else if (slv->grid.nd == 3) {
		ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
		slv->dmap = stella_dmap_create_3d(xs - offset[0], ys - offset[1], zs - offset[2], xm, ym, zm, stride);
	} else {
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Unsupported dimmension: %d", slv->grid.nd);
	}
	ierr = DMDAGetInfo(slv->dm, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

	ierr = KSPCreate(slv->comm, &slv->ksp);CHKERRQ(ierr);
	if (stella_log(slv, STELLA_LOG_EIGS)) {
		ierr = KSPSetComputeEigenvalues(slv->ksp, PETSC_TRUE);CHKERRQ(ierr);
	}

	ierr = KSPSetFromOptions(slv->ksp);CHKERRQ(ierr);

	if (!slv->options.algebraic) {
		ierr = KSPGetPC(slv->ksp, &slv->pc);
		ierr = PCSetType(slv->pc, PCSHELL);CHKERRQ(ierr);
		ierr = stella_pc_create(&shell);CHKERRQ(ierr);
		ierr = PCShellSetApply(slv->pc, stella_pc_apply);CHKERRQ(ierr);
		ierr = PCShellSetDestroy(slv->pc, stella_pc_destroy);CHKERRQ(ierr);
		ierr = PCShellSetContext(slv->pc, shell);CHKERRQ(ierr);
		ierr = PCShellSetName(slv->pc, "cedar");CHKERRQ(ierr);
		ierr = PCShellSetSetUp(slv->pc, stella_pc_setup);CHKERRQ(ierr);
	}

	ierr = stella_fd_create(&slv->fd, STELLA_FD_STANDARD_O2);CHKERRQ(ierr);
	ierr = stella_operator_create(&slv->op, &slv->level, slv->fd, slv->grid.nd);CHKERRQ(ierr);
	ierr = stella_boundary_create(&slv->boundary, &slv->level,
	                              slv->dmap, &slv->state, slv->fd);CHKERRQ(ierr);
	slv->op->axisymmetric = slv->options.axisymmetric;
	slv->boundary->axisymmetric = slv->options.axisymmetric;

	ierr = DMCreateGlobalVector(slv->dm, &slv->level.dcoef);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(slv->dm, &slv->level.ldcoef);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->level.bcoef);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(slv->dm, &slv->level.lbcoef);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->level.add_cont);CHKERRQ(ierr);
	ierr = DMCreateLocalVector(slv->dm, &slv->level.ladd_cont);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->level.nscale);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->level.pm);CHKERRQ(ierr);
	ierr = VecSet(slv->level.pm, 1.0);CHKERRQ(ierr);
	ierr = VecSet(slv->level.add_cont, 0.0);CHKERRQ(ierr);
	ierr = VecSet(slv->level.nscale, 1.0);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->rhs);CHKERRQ(ierr);
	ierr = DMCreateGlobalVector(slv->dm, &slv->x);CHKERRQ(ierr);

	if (!slv->options.algebraic) {
		#ifdef WITH_CEDAR
		stella_bmg_mat *mat_ctx = (stella_bmg_mat*) malloc(sizeof(stella_bmg_mat));
		int cedar_err;
		if (slv->grid.nd == 2) {
			cedar_err = cedar_mat_create2d(slv->grid.topo, CEDAR_STENCIL_NINE_PT, &mat_ctx->so);chkerr(cedar_err);
			cedar_err = cedar_vec_create2d(slv->grid.topo, &mat_ctx->solvec);chkerr(cedar_err);
			cedar_err = cedar_vec_create2d(slv->grid.topo, &mat_ctx->rhsvec);chkerr(cedar_err);
			mat_ctx->nd = 2;
			ierr = MatCreateShell(slv->comm, xm*ym, xm*ym, ngx*ngy, ngx*ngy, mat_ctx, &slv->A);CHKERRQ(ierr);
		} else {
			cedar_err = cedar_mat_create3d(slv->grid.topo, CEDAR_STENCIL_XXVII_PT, &mat_ctx->so);chkerr(cedar_err);
			cedar_err = cedar_vec_create3d(slv->grid.topo, &mat_ctx->solvec);chkerr(cedar_err);
			cedar_err = cedar_vec_create3d(slv->grid.topo, &mat_ctx->rhsvec);chkerr(cedar_err);
			mat_ctx->nd = 3;
			ierr = MatCreateShell(slv->comm, xm*ym*zm, xm*ym*zm, ngx*ngy*ngz, ngx*ngy*ngz, mat_ctx, &slv->A);CHKERRQ(ierr);
		}
		ierr = MatShellSetOperation(slv->A, MATOP_MULT, (void(*)(void)) stella_bmg_mult);CHKERRQ(ierr);

		// If not custom cedar config filename was set fall back to config.json
		if(slv->cedar_config_filename == NULL || strlen(slv->cedar_config_filename) == 0)
			cedar_err = cedar_config_create("config.json", &mat_ctx->conf);
		else
			cedar_err = cedar_config_create(slv->cedar_config_filename, &mat_ctx->conf);

		// If the specified config does not exist, inform user and fall back to config.json
		if (cedar_err == CEDAR_ERR_FNAME) {
			char err_str[80];
			sprintf(err_str, "Cedar config not found: %s - Falling back to config.json.", slv->cedar_config_filename);
			stella_io_print(PETSC_COMM_WORLD, err_str);
			cedar_config_create("config.json", &mat_ctx->conf);
		}

		#endif
	} else {
		ierr = DMCreateMatrix(slv->dm, &slv->A);CHKERRQ(ierr);
	}

	if (stella_log(slv, STELLA_LOG_STATUS)) {
		ierr = stella_io_print(slv->comm, "set-up is finished");CHKERRQ(ierr);
	}
	return 0;
}


PetscErrorCode stella_init(stella **solver_ctx, MPI_Comm comm,
                           int nGlobal[], int nProcs[], int nLocal[],
                           int offset[], int stride[],
                           int cartCoord[], int periodic[], int periodic_storage,
                           int nd, int axisymmetric)
{
	PetscErrorCode ierr;

	stella *slv = (stella*) malloc(sizeof(stella));
	slv->comm = comm;
	slv->options.axisymmetric = axisymmetric;
	{
		PetscBool flg;
		ierr = petsc_options_has_name("-cedar", &flg);CHKERRQ(ierr);
		if (flg) {
		#ifdef WITH_CEDAR
			slv->options.algebraic = 0;
		#endif
		} else {
			slv->options.algebraic = 1;
		}
	}

	stella_set_log(slv, STELLA_LOG_STATUS | STELLA_LOG_RESIDUAL);
	ierr = stella_grid_setup(&slv->grid, &slv->dm, &slv->comm, nGlobal, nProcs,
	                         nLocal, cartCoord, periodic, periodic_storage, nd);CHKERRQ(ierr);
	slv->level.dm = slv->dm;
	ierr = stella_setup(slv, offset, stride);CHKERRQ(ierr);

	*solver_ctx = slv;

	return 0;
}


PetscErrorCode stella_solve(stella *slv)
{
	PetscErrorCode ierr;
	PetscBool cedar_direct = 0;
	slv->ts++;

	ierr = petsc_options_has_name("-cedar_direct", &cedar_direct);CHKERRQ(ierr);

	// Temporarily signal bcs have changed to be safe
	ierr = stella_changed_bc(slv);CHKERRQ(ierr);
	// Reuse previous solution at subsequent iteration
	if(slv->ts == 1) {
		ierr = VecSet(slv->x, 0);CHKERRQ(ierr);
	}
	if (!slv->options.algebraic && cedar_direct) {
		PC bmg_pc;
		ierr = KSPGetPC(slv->ksp, &bmg_pc);CHKERRQ(ierr);
		ierr = PCApply(bmg_pc, slv->rhs, slv->x);CHKERRQ(ierr);
	} else {
		ierr = KSPSolve(slv->ksp, slv->rhs, slv->x);CHKERRQ(ierr);
	}

	if (slv->state.sol) {
	  ierr = stella_io_vtkwrite(slv->dm, slv->x, "phi", 0, slv->ts);CHKERRQ(ierr);
	  ierr = stella_io_vtkwrite(slv->dm, slv->sol, "sol", 0, 0);CHKERRQ(ierr);
	  Vec err;
	  ierr = VecDuplicate(slv->sol, &err);CHKERRQ(ierr);
	  ierr = VecWAXPY(err, -1, slv->x, slv->sol);CHKERRQ(ierr);
	  ierr = VecAbs(err);CHKERRQ(ierr);
	  ierr = stella_io_vtkwrite(slv->dm, err, "err", 0, 0);CHKERRQ(ierr);
	  ierr = VecDestroy(&err);CHKERRQ(ierr);
	}

	// Update external array with solution
	ierr = update_solution(slv);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_set_external(stella *slv, double phi[], double dcoef[],
                                   double bcoef[], double jump_x[], double jump_y[], double jump_z[])
{
	PetscErrorCode ierr;

	slv->state.phi = phi;
	slv->state.dcoef = dcoef;
	slv->state.bcoef = bcoef;
	slv->state.jump[0] = jump_x;
	slv->state.jump[1] = jump_y;
	slv->state.jump[2] = jump_z;
	slv->state.sol = NULL;

	ierr = stella_store_external_array(slv, slv->state.dcoef, slv->level.dcoef);CHKERRQ(ierr);
	// populate halo region
	ierr = DMGlobalToLocalBegin(slv->dm, slv->level.dcoef, INSERT_VALUES, slv->level.ldcoef);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(slv->dm, slv->level.dcoef, INSERT_VALUES, slv->level.ldcoef);CHKERRQ(ierr);

	ierr = stella_store_external_array(slv, slv->state.bcoef, slv->level.bcoef);CHKERRQ(ierr);
	// populate halo region
	ierr = DMGlobalToLocalBegin(slv->dm, slv->level.bcoef, INSERT_VALUES, slv->level.lbcoef);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(slv->dm, slv->level.bcoef, INSERT_VALUES, slv->level.lbcoef);CHKERRQ(ierr);

	return 0;
}



PetscErrorCode stella_set_rhs(stella *slv, double rhs[])
{
	PetscErrorCode ierr;

	slv->state.rhs = rhs;
	ierr = stella_changed_rhs(slv);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_set_sol(stella *slv, double sol[])
{
	PetscErrorCode ierr;

	ierr = DMCreateGlobalVector(slv->dm, &slv->sol);
	slv->state.sol = sol;

	ierr = stella_store_external_array(slv, slv->state.sol, slv->sol);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_set_grid(stella *slv, int is[], int ie[], int num_cells, double xyz[])
{
	PetscErrorCode ierr;
	PetscInt i, j, k;
	PetscInt xs, ys, zs, xm, ym, zm;
	DM cda;
	Vec gc;

	for (i = 0; i < slv->grid.nd; i++) {
		slv->grid.is[i] = is[i];
		slv->grid.ie[i] = ie[i];
	}
	slv->grid.num_cells = num_cells;
	slv->grid.xyz = xyz;

	ierr = DMGetCoordinateDM(slv->dm, &cda);CHKERRQ(ierr);
	ierr = DMGetCoordinates(slv->dm, &gc);CHKERRQ(ierr);
	ierr = DMDAGetCorners(cda, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	if (slv->grid.nd == 2) {
		DMDACoor2d **coors;
		double **x, **y;
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz, &x);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz + num_cells, &y);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda, gc, &coors);CHKERRQ(ierr);

		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				coors[j][i].x = x[j][i];
				coors[j][i].y = y[j][i];
			}
		}

		ierr = stella_dmap_restore(slv->dmap, &x);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &y);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(cda, gc, &coors);CHKERRQ(ierr);
	} else if (slv->grid.nd == 3) {
		DMDACoor3d ***coors;
		double ***x, ***y, ***z;
		ierr = DMDAVecGetArray(cda, gc, &coors);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz, &x);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz + num_cells, &y);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->grid.xyz + 2*num_cells, &z);CHKERRQ(ierr);

		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {
					coors[k][j][i].x = x[k][j][i];
					coors[k][j][i].y = y[k][j][i];
					coors[k][j][i].z = z[k][j][i];
				}
			}
		}

		ierr = stella_dmap_restore(slv->dmap, &x);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &y);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &z);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(cda, gc, &coors);CHKERRQ(ierr);

	}

	ierr = stella_metric_create(&slv->level.metric, slv->dm, slv->fd);CHKERRQ(ierr);

	if (stella_log(slv, STELLA_LOG_PROBLEM)) {
		PetscViewer vout;
		ierr = PetscViewerBinaryOpen(slv->comm, "grid.dat", FILE_MODE_WRITE, &vout);CHKERRQ(ierr);
		ierr = VecView(gc, vout);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&vout);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_set_boundary(stella *slv, stella_ptypes ptypes,
                                   char classify[], char norm_dir[], double values[])
{
	PetscErrorCode ierr;

	ierr = stella_classify_create(&slv->level.classify, slv->dm, slv->dmap, classify, ptypes);CHKERRQ(ierr);
	ierr = stella_boundary_set(slv->boundary, norm_dir, values);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_setup_op(stella *slv)
{
	PetscErrorCode ierr;
    Vec Vec_tmp;
    PetscScalar    **array;
    PetscInt         m,n;

	ierr = stella_operator_assemble(slv->op, slv->A, slv->dm);CHKERRQ(ierr);
	ierr = stella_boundary_apply(slv->boundary, slv->A, slv->dm);CHKERRQ(ierr);

	if (stella_log(slv, STELLA_LOG_STATUS)) {
		ierr = stella_io_print(slv->comm, "Matrix assembled");CHKERRQ(ierr);
	}

	if (slv->options.algebraic) {
		ierr = MatAssemblyBegin(slv->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		ierr = MatAssemblyEnd(slv->A, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
		if (stella_log(slv, STELLA_LOG_PROBLEM)) {
			PetscViewer vout;
			ierr = PetscViewerBinaryOpen(slv->comm, "a.dat", FILE_MODE_WRITE, &vout);CHKERRQ(ierr);
			ierr = MatView(slv->A, vout);CHKERRQ(ierr);
			ierr = PetscViewerDestroy(&vout);CHKERRQ(ierr);
		}
	} else {
		#ifdef WITH_CEDAR
		if (stella_log(slv, STELLA_LOG_PROBLEM)) {
			stella_bmg_mat *ctx;
			ierr = MatShellGetContext(slv->A, (void**) &ctx);CHKERRQ(ierr);
			int cedar_err;
			cedar_err = cedar_mat_dump(ctx->so);chkerr(cedar_err);
		}
		#endif
	}

	ierr = KSPSetOperators(slv->ksp, slv->A, slv->A);CHKERRQ(ierr);

	return 0;
}


int stella_log(stella *slv, stella_log_level level)
{
	return slv->options.log_level & level;
}


void stella_set_log(stella *slv, stella_log_level level)
{
	slv->options.log_level = level;
}


PetscErrorCode stella_cleanup(stella *slv)
{
	PetscErrorCode ierr;
	int i;

	ierr = stella_boundary_destroy(slv->boundary);CHKERRQ(ierr);
	ierr = stella_operator_destroy(slv->op);CHKERRQ(ierr); free(slv->op);
	ierr = stella_fd_destroy(slv->fd);CHKERRQ(ierr); free(slv->fd);
	stella_dmap_destroy(slv->dmap); free(slv->dmap);
	ierr = stella_metric_destroy(slv->level.metric);CHKERRQ(ierr);
	free(slv->level.metric);

	ierr = DMDestroy(&slv->dm);CHKERRQ(ierr);
	ierr = KSPDestroy(&slv->ksp);CHKERRQ(ierr);

	return 0;
}
