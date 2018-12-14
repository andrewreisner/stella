#include "stella_io.h"
#include "stella_util.h"
#include "stella_signals.h"
#include "stella_mat.h"

// private signal helpers

static inline int smallerValue(int i1, int i2, double v1, double v2)
{
	return   (v1==v2) ? i1 : ( v1<v2 ? i1 : i2 );
}


/**
 * Contribution to rhs needed for axisymmetric solve
 */
static PetscErrorCode contribute_axisymmetric(stella *slv)
{
	PetscErrorCode ierr;
	PetscInt i, j, k, xs, ys, zs, xm, ym, zm, ngx, ngy, ngz;
	Vec b;

	b = slv->rhs;

	ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(slv->dm, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

	if (slv->grid.nd == 2) {
		double **rhs;
		PetscScalar **bvec;
		DM cda;
		DMDACoor2d **coors;
		Vec gc;

		ierr = DMGetCoordinateDM(slv->dm, &cda);CHKERRQ(ierr);
		ierr = DMGetCoordinates(slv->dm, &gc);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda, gc, &coors);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, b, &bvec);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.rhs, &rhs);CHKERRQ(ierr);

		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				bvec[j][i] = rhs[j][i] * coors[j][i].x;
			}
		}

		ierr = DMDAVecRestoreArray(cda, gc, &coors);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, b, &bvec);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &rhs);CHKERRQ(ierr);
	} else {
		SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP,"Axisymmetric contributions not supported for 3D");
	}

	return 0;
}


/**
 * Contribution to rhs needed for interface (jump) condition
 */
static PetscErrorCode contribute_interface(stella *slv)
{
	PetscErrorCode ierr;
	PetscInt i, j, k, xs, ys, zs, xm, ym, zm, ngx, ngy, ngz;
	Vec b;

	b = slv->rhs;

	ierr = DMDAGetCorners(slv->dm, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAGetInfo(slv->dm, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);

	if (slv->grid.nd == 2) {
		double **rhs, **jump_x, **jump_y;
		PetscScalar **bvec, **dcoef;

		ierr = stella_dmap_get(slv->dmap, slv->state.rhs, &rhs);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.jump[0], &jump_x);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.jump[1], &jump_y);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, slv->level.ldcoef, &dcoef);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, b, &bvec);CHKERRQ(ierr);

		double fxp, fxm, fyp, fym;
		PetscScalar **jac[4];
		PetscScalar **x_s, **x_t, **y_s, **y_t;
		stella_metric *met;
		met = slv->level.metric;

		for (i = 0; i < 4; i++) {
			ierr = DMDAVecGetArray(slv->dm, met->ljac_v[i], &jac[i]);CHKERRQ(ierr);
		}

		x_s = jac[met->t2map[0][0]];
		x_t = jac[met->t2map[0][1]];
		y_s = jac[met->t2map[1][0]];
		y_t = jac[met->t2map[1][1]];

		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {

				if ((i != ngx-1) && (dcoef[j][i] != dcoef[j][i+1]))
					fxp = 2.0*dcoef[j][i]*jump_x[j][i] / (dcoef[j][i] + dcoef[j][i+1]) / (x_s[j][i]+x_s[j][i+1]);
				else fxp = 0;
				if ((i != 0) && (dcoef[j][i] != dcoef[j][i-1]))
					fxm = 2.0*dcoef[j][i]*jump_x[j][i] / (dcoef[j][i] + dcoef[j][i-1]) / (x_s[j][i]+x_s[j][i-1]);
				else fxm = 0;
				if ((j != ngy-1) && (dcoef[j][i] != dcoef[j+1][i]))
					fyp = 2.0*dcoef[j][i]*jump_y[j][i] / (dcoef[j][i] + dcoef[j+1][i]) / (y_t[j][i]+y_t[j+1][i]);
				else fyp = 0;
				if ((j != 0) && (dcoef[j][i] != dcoef[j-1][i]))
					fym = 2.0*dcoef[j][i]*jump_y[j][i] / (dcoef[j][i] + dcoef[j-1][i]) / (y_t[j][i]+y_t[j-1][i]);
				else fym = 0;

				bvec[j][i] = rhs[j][i] - fxp - fxm - fyp - fym;
			}
		}

		for (i = 0; i < 4; i++) {
			ierr = DMDAVecRestoreArray(slv->dm, met->ljac_v[i], &jac[i]);CHKERRQ(ierr);
		}

		ierr = stella_dmap_restore(slv->dmap, &jump_x);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &jump_y);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &rhs);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, b, &bvec);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, slv->level.ldcoef, &dcoef);CHKERRQ(ierr);
	} else { // 3D
		double ***rhs, ***jump_x, ***jump_y, ***jump_z;
		PetscScalar ***bvec, ***dcoef;

		double fxp, fxm, fyp, fym, fzp, fzm;
		PetscScalar ***jac[9];
		PetscScalar ***x_r, ***y_s, ***z_t;
		stella_metric *met;
		met = slv->level.metric;

		for (i = 0; i < 9; i++) {
			ierr = DMDAVecGetArray(slv->dm, met->ljac_v[i], &jac[i]);CHKERRQ(ierr);
		}

		x_r = jac[met->t3map[0][0]];
		y_s = jac[met->t3map[1][1]];
		z_t = jac[met->t3map[2][2]];

		ierr = DMDAVecGetArray(slv->dm, slv->level.ldcoef, &dcoef);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.rhs, &rhs);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.jump[0], &jump_x);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.jump[1], &jump_y);CHKERRQ(ierr);
		ierr = stella_dmap_get(slv->dmap, slv->state.jump[2], &jump_z);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(slv->dm, b, &bvec);CHKERRQ(ierr);
		for (k = zs; k < zs + zm; k++) {
			for (j = ys; j < ys + ym; j++) {
				for (i = xs; i < xs + xm; i++) {

					if ((i != ngx-1) && (dcoef[k][j][i] != dcoef[k][j][i+1]))
						fxp = 2.0*dcoef[k][j][i]*jump_x[k][j][i] / (dcoef[k][j][i] + dcoef[k][j][i+1]) / (x_r[k][j][i]+x_r[k][j][i+1]);
					else fxp = 0;
					if ((i != 0) && (dcoef[k][j][i] != dcoef[k][j][i-1]))
						fxm = 2.0*dcoef[k][j][i]*jump_x[k][j][i] / (dcoef[k][j][i] + dcoef[k][j][i-1]) / (x_r[k][j][i]+x_r[k][j][i-1]);
					else fxm = 0;
					if ((j != ngy-1) && (dcoef[k][j][i] != dcoef[k][j+1][i]))
						fyp = 2.0*dcoef[k][j][i]*jump_y[k][j][i] / (dcoef[k][j][i] + dcoef[k][j+1][i]) / (y_s[k][j][i]+y_s[k][j+1][i]);
					else fyp = 0;
					if ((j != 0) && (dcoef[k][j][i] != dcoef[k][j-1][i]))
						fym = 2.0*dcoef[k][j][i]*jump_y[k][j][i] / (dcoef[k][j][i] + dcoef[k][j-1][i]) / (y_s[k][j][i]+y_s[k][j-1][i]);
					else fym = 0;
					if ((k != ngz-1) && (dcoef[k][j][i] != dcoef[k+1][j][i]))
						fzp = 2.0*dcoef[k][j][i]*jump_z[k][j][i] / (dcoef[k][j][i] + dcoef[k+1][j][i]) / (z_t[k][j][i]+z_t[k+1][j][i]);
					else fzp = 0;
					if ((k != 0) && (dcoef[k][j][i] != dcoef[k-1][j][i]))
						fzm = 2.0*dcoef[k][j][i]*jump_z[k][j][i] / (dcoef[k][j][i] + dcoef[k-1][j][i]) / (z_t[k][j][i]+z_t[k-1][j][i]);
					else fzm = 0;


					bvec[k][j][i] = rhs[k][j][i] - fxp - fxm - fyp - fym - fzp - fzm;
				}
			}
		}

		for (i = 0; i < 9; i++) {
			ierr = DMDAVecRestoreArray(slv->dm, met->ljac_v[i], &jac[i]);CHKERRQ(ierr);
		}
		ierr = stella_dmap_restore(slv->dmap, &rhs);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &jump_x);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &jump_y);CHKERRQ(ierr);
		ierr = stella_dmap_restore(slv->dmap, &jump_z);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, b, &bvec);CHKERRQ(ierr);
		ierr = DMDAVecRestoreArray(slv->dm, slv->level.ldcoef, &dcoef);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_changed_bc(stella * slv)
{
	PetscErrorCode ierr;
	ierr = VecSet(slv->level.pm, 1.0);CHKERRQ(ierr);
	ierr = VecSet(slv->level.add_cont, 0.0);CHKERRQ(ierr);
	ierr = VecSet(slv->level.nscale, 1.0);CHKERRQ(ierr);

	ierr = stella_boundary_apply(slv->boundary, slv->A, slv->dm);CHKERRQ(ierr);

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
			cedar_mat_dump(ctx->so);
		}
		#endif
	}

	ierr = KSPSetOperators(slv->ksp, slv->A, slv->A);CHKERRQ(ierr);

	ierr = stella_changed_rhs(slv);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode stella_changed_rhs(stella *slv)
{
	PetscErrorCode ierr;
	Vec b;

	b = slv->rhs;

	if (slv->op->axisymmetric)
		contribute_axisymmetric(slv);
	else
		contribute_interface(slv);

	// 1/jac term moved to rhs
	ierr = VecPointwiseMult(b, slv->level.metric->jac, b);CHKERRQ(ierr);

	ierr = stella_boundary_apply_rhs(slv->boundary, slv->dm, slv->rhs);CHKERRQ(ierr);

	ierr = VecPointwiseMult(slv->level.add_cont, slv->level.pm, slv->level.add_cont);CHKERRQ(ierr);
	ierr = VecAXPY(b, 1, slv->level.add_cont);CHKERRQ(ierr);
	ierr = VecPointwiseMult(b, slv->level.nscale, b);CHKERRQ(ierr);

	if (stella_log(slv, STELLA_LOG_PROBLEM)) {
		PetscViewer bout;
		PetscViewerASCIIOpen(slv->comm, "b.txt", &bout);
		ierr = VecView(b, bout);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_changed_dcoef(stella *slv)
{
	PetscErrorCode ierr;

	ierr = stella_store_external_array(slv, slv->state.dcoef, slv->level.dcoef);CHKERRQ(ierr);
	// populate halo region
	ierr = DMGlobalToLocalBegin(slv->dm, slv->level.dcoef, INSERT_VALUES, slv->level.ldcoef);CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(slv->dm, slv->level.dcoef, INSERT_VALUES, slv->level.ldcoef);CHKERRQ(ierr);

	ierr = stella_operator_assemble(slv->op, slv->A, slv->dm);CHKERRQ(ierr);

	if (stella_log(slv, STELLA_LOG_STATUS)) {
		ierr = stella_io_print(slv->comm, "Matrix assembled");CHKERRQ(ierr);
	}

	ierr = stella_changed_bc(slv);CHKERRQ(ierr);

	return 0;
}
