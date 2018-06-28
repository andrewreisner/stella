#include "stella_dirichlet.h"
#include "stella_mat.h"
#include "stella_operator.h"


/**
 * Deletes connections from interior points that neighbor a Dirichlet boundary.
 * This step is needed to get a symmetric operator.
 */
static PetscErrorCode update_dirichlet_sym(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt i,j,cnt;
	PetscInt ngx, ngy;
	PetscScalar v[9];
	MatStencil row[9], col;
	stella_classify *cls = bc->level->classify;
	char **classify;

	PetscInt xs, ys, xm, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	char dirichlet = cls->ptypes.dirichlet;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);

	for (i = 0; i < 8; i++) v[i] = 0;

	for (j = ys; j < ys + ym; j++) {
		for (i = xs; i < xs + xm; i++) {
			if (classify[j][i] == dirichlet) {
				col.i = i; col.j = j;
				cnt = 0;
				if (i-1 >= 0 && classify[j][i-1] != dirichlet) {
					row[cnt].i = i-1; row[cnt].j = j;
					cnt++;
					if (j-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j-1;
						cnt++;
					}
					if (j+1 <= ngy-1) {
						row[cnt].i = i-1; row[cnt].j = j+1;
						cnt++;
					}
				}
				if (i+1 <= ngx-1 && classify[j][i+1] != dirichlet) {
					row[cnt].i = i+1; row[cnt].j = j;
					cnt++;
					if (j-1 >= 0) {
						row[cnt].i = i+1; row[cnt].j = j-1;
						cnt++;
					}
					if (j+1 <= ngy-1) {
						row[cnt].i = i+1; row[cnt].j = j+1;
						cnt++;
					}
				}
				if (j-1 >= 0 && classify[j-1][i] != dirichlet) {
					row[cnt].i = i; row[cnt].j = j-1;
					cnt++;
					if (i-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j-1;
						cnt++;
					}
					if (i+1 <= ngx - 1) {
						row[cnt].i = i+1; row[cnt].j = j-1;
						cnt++;
					}
				}
				if (j+1 <= ngy-1 && classify[j+1][i] != dirichlet) {
					row[cnt].i = i; row[cnt].j = j+1;
					cnt++;
					if (i-1 >= 0) {
						row[cnt].i = i-1; row[cnt].j = j+1;
						cnt++;
					}
					if (i+1 <= ngx - 1) {
						row[cnt].i = i+1; row[cnt].j = j+1;
						cnt++;
					}
				}

				ierr = MatSetValuesStencil(A, cnt, row, 1, &col, v, INSERT_VALUES);CHKERRQ(ierr);
			}
		}
	}

	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode update_dirichlet_sym_3d(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt i, j, k, cnt;
	PetscInt ngx, ngy, ngz;
	PetscScalar v[64];
	MatStencil row[64], col;
	stella_classify *cls = bc->level->classify;
	char ***classify;

	PetscInt xs, ys, zs, xm, ym, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	char dirichlet = cls->ptypes.dirichlet;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);

	for (i = 0; i < 19; i++) v[i] = 0;
	for (k = zs; k < zs + zm; k++) {
		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				if (classify[k][j][i] == dirichlet) {
					col.i = i; col.j = j; col.k = k;
					cnt = 0;

					if (i-1 >= 0 && classify[k][j][i-1] != dirichlet) {
						row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k;
						cnt++;
						if (j-1 >= 0) {
							row[cnt].i = i-1; row[cnt].j = j-1; row[cnt].k = k;
							cnt++;
						}
						if (j+1 <= ngy-1) {
							row[cnt].i = i-1; row[cnt].j = j+1; row[cnt].k = k;
							cnt++;
						}
						if (k-1 >= 0) {
							row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k-1;
							cnt++;
						}
						if (k+1 <= ngz - 1) {
							row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k+1;
							cnt++;
						}
					}
					if (i+1 <= ngx - 1 && classify[k][j][i+1] != dirichlet) {
						row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k;
						cnt++;
						if (j-1 >= 0) {
							row[cnt].i = i+1; row[cnt].j = j-1; row[cnt].k = k;
							cnt++;
						}
						if (j+1 <= ngy-1) {
							row[cnt].i = i+1; row[cnt].j = j+1; row[cnt].k = k;
							cnt++;
						}
						if (k-1 >= 0) {
							row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k-1;
							cnt++;
						}
						if (k+1 <= ngz - 1) {
							row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k+1;
							cnt++;
						}
					}
					if (j-1 >= 0 && classify[k][j-1][i] != dirichlet) {
						row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k;
						cnt++;
						if (i-1 >= 0) {
							row[cnt].i = i-1; row[cnt].j = j-1; row[cnt].k = k;
							cnt++;
						}
						if (i+1 <= ngx - 1) {
							row[cnt].i = i+1; row[cnt].j = j-1; row[cnt].k = k;
							cnt++;
						}
						if (k-1 >= 0) {
							row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k-1;
							cnt++;
						}
						if (k+1 <= ngz - 1) {
							row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k+1;
							cnt++;
						}
					}
					if (j+1 <= ngy - 1 && classify[k][j+1][i] != dirichlet) {
						row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k;
						cnt++;
						if (i-1 >= 0) {
							row[cnt].i = i-1; row[cnt].j = j+1; row[cnt].k = k;
							cnt++;
						}
						if (i+1 <= ngx - 1) {
							row[cnt].i = i+1; row[cnt].j = j+1; row[cnt].k = k;
							cnt++;
						}
						if (k-1 >= 0) {
							row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k-1;
							cnt++;
						}
						if (k+1 <= ngz - 1) {
							row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k+1;
							cnt++;
						}
					}
					if (k-1 >= 0 && classify[k-1][j][i] != dirichlet) {
						row[cnt].i = i; row[cnt].j = j; row[cnt].k = k-1;
						cnt++;
						if (i-1 >= 0) {
							row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k-1;
							cnt++;
						}
						if (i+1 <= ngx - 1) {
							row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k-1;
							cnt++;
						}
						if (j-1 >= 0) {
							row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k-1;
							cnt++;
						}
						if (j+1 <= ngy - 1) {
							row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k-1;
							cnt++;
						}
					}
					if (k+1 <= ngz - 1 && classify[k+1][j][i] != dirichlet) {
						row[cnt].i = i; row[cnt].j = j; row[cnt].k = k+1;
						cnt++;
						if (i-1 >= 0) {
							row[cnt].i = i-1; row[cnt].j = j; row[cnt].k = k+1;
							cnt++;
						}
						if (i+1 <= ngx - 1) {
							row[cnt].i = i+1; row[cnt].j = j; row[cnt].k = k+1;
							cnt++;
						}
						if (j-1 >= 0) {
							row[cnt].i = i; row[cnt].j = j-1; row[cnt].k = k+1;
							cnt++;
						}
						if (j+1 <= ngy - 1) {
							row[cnt].i = i; row[cnt].j = j+1; row[cnt].k = k+1;
							cnt++;
						}
					}

					ierr = MatSetValuesStencil(A, cnt, row, 1, &col, v, INSERT_VALUES);CHKERRQ(ierr);
				}
			}
		}
	}

	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);

	return 0;
}



/**
 * Eliminates connections from boundary to neighboring points.
 * Saves the weight of this connection as additive contribution to RHS.
 */
static PetscErrorCode apply_dirichlet(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	DM cda;
	DMDACoor2d **coors;
	Vec lc;
	PetscInt i,j,cnt;
	PetscInt ngx, ngy;
	PetscScalar v[9];
	MatStencil row, col[9];
	MatType mtype;
	PetscBool is_boxmg;
	PetscScalar **acont;
	stella_metric *met;
	PetscScalar **dcoef;
	PetscScalar **coef[5];

	stella_classify *cls = bc->level->classify;
	char **classify;
	double **bcvals;
	char dirichlet = cls->ptypes.dirichlet;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);
	ierr = stella_dmap_get(bc->slv_dmap, bc->values, &bcvals);CHKERRQ(ierr);

	PetscInt xs, ys, xm, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);

	met = bc->level->metric;

	ierr = MatGetType(A, &mtype);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->ladd_cont, &acont);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);

	if (bc->axisymmetric) {
		ierr = DMGetCoordinateDM(da, &cda);CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da, &lc);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda, lc, &coors);CHKERRQ(ierr);
	}

	for (i = 0; i < 5; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	v[0] = 1.0;
	for (j = ys; j < ys + ym; j++) {
		for (i = xs; i < xs + xm; i++) {
			if (classify[j][i] == dirichlet) {
				row.i = i; row.j = j;
				col[0].i = i; col[0].j = j;
				cnt = 1;

				if (i != 0) {
					col[cnt].i = i-1; col[cnt].j = j;
					v[cnt] = 0; cnt++;
					if (classify[j][i-1] != dirichlet) {
						double weight = coef[MET_W][j][i]*have(dcoef[j][i],dcoef[j][i-1])*bcvals[j][i];
						if (bc->axisymmetric) weight *= (coors[j][i].x + coors[j][i-1].x) * .5;
						acont[j][i-1] += weight;
					}
				}
				if (j != 0) {
					col[cnt].i = i; col[cnt].j = j-1;
					v[cnt] = 0; cnt++;
					if (classify[j-1][i] != dirichlet) {
						double weight = coef[MET_S][j][i]*have(dcoef[j][i],dcoef[j-1][i])*bcvals[j][i];
						if (bc->axisymmetric) weight *= coors[j][i].x;
						acont[j-1][i] += weight;
					}
				}
				// SW
				if (i != 0 && j != 0) {
					col[cnt].i = i-1; col[cnt].j = j-1;
					v[cnt] = 0; cnt++;
					if ((classify[j][i-1] != dirichlet) || (classify[j-1][i] != dirichlet)) {
						double weight = coef[MET_SW][j][i]*have(have(dcoef[j][i], dcoef[j][i-1]),
						                                        have(dcoef[j-1][i], dcoef[j-1][i-1])) * bcvals[j][i];
						acont[j-1][i-1] += weight;
					}
				}
				if (i != ngx - 1) {
					col[cnt].i = i+1; col[cnt].j = j;
					v[cnt] = 0; cnt++;
					if (classify[j][i+1] != dirichlet) {
						double weight = coef[MET_W][j][i+1]*have(dcoef[j][i],dcoef[j][i+1])*bcvals[j][i];
						if (bc->axisymmetric) weight *= (coors[j][i+1].x + coors[j][i].x) * .5;
						acont[j][i+1] += weight;
					}
				}
				// SE
				if (i != ngx - 1 && j != 0) {
					col[cnt].i = i+1; col[cnt].j = j-1;
					v[cnt] = 0; cnt++;
					if ((classify[j][i+1] != dirichlet) || (classify[j-1][i] != dirichlet)) {
						double weight = coef[MET_SE][j][i]*have(have(dcoef[j][i], dcoef[j][i+1]),
						                                        have(dcoef[j-1][i], dcoef[j-1][i+1]))*bcvals[j][i];
						acont[j-1][i+1] += weight;
					}
				}
				// NW
				if (i !=0 && j != ngy - 1) {
					col[cnt].i = i-1; col[cnt].j = j+1;
					v[cnt] = 0; cnt++;
					if  ((classify[j][i-1] != dirichlet) || (classify[j+1][i] != dirichlet)) {
						double weight = coef[MET_SE][j+1][i-1]*have(have(dcoef[j][i], dcoef[j][i-1]),
						                                            have(dcoef[j+1][i], dcoef[j+1][i-1]))*bcvals[j][i];
						acont[j+1][i-1] += weight;
					}
				}
				if (j != ngy - 1) {
					col[cnt].i = i; col[cnt].j = j+1;
					v[cnt] = 0; cnt++;
					if (classify[j+1][i] != dirichlet) {
						double weight = coef[MET_S][j+1][i]*have(dcoef[j][i],dcoef[j+1][i])*bcvals[j][i];
						if (bc->axisymmetric) weight *= coors[j][i].x;
						acont[j+1][i] += weight;
					}
				}
				// NE
				if ((i != ngx - 1) && (j != ngy - 1)) {
					col[cnt].i = i+1; col[cnt].j = j+1;
					v[cnt] = 0; cnt++;
					if ((classify[j][i+1] != dirichlet) || (classify[j+1][i] != dirichlet)) {
						double weight = coef[MET_SW][j+1][i+1]*have(have(dcoef[j][i], dcoef[j][i+1]),
						                                            have(dcoef[j+1][i], dcoef[j+1][i+1]))*bcvals[j][i];
						acont[j+1][i+1] += weight;
					}
				}

				ierr = PetscStrcmp(mtype, MATSHELL, &is_boxmg);CHKERRQ(ierr);
				if (is_boxmg) {
					ierr = stella_bmg_SetValuesStencil(A, 1, &row, cnt, col, v, INSERT_VALUES);CHKERRQ(ierr);
				} else {
					ierr = MatSetValuesStencil(A, 1, &row, cnt, col, v, INSERT_VALUES);CHKERRQ(ierr);
				}
			}
		}
	}

	ierr = stella_dmap_restore(bc->slv_dmap, &bcvals);CHKERRQ(ierr);
	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->ladd_cont, &acont);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);
	for (i = 0; i < 5; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	if (bc->axisymmetric) {
		ierr = DMDAVecRestoreArray(cda, lc, &coors);CHKERRQ(ierr);
	}

	return 0;
}


static PetscErrorCode apply_dirichlet_3d(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,cnt;
	PetscInt ngx, ngy, ngz;
	PetscScalar v[19];
	MatStencil row, col[19];
	MatType mtype;
	PetscBool is_boxmg;
	double dcoefh[19];
	PetscScalar ***acont;
	stella_metric *met;
	PetscScalar ***dcoef;
	PetscScalar ***coef[10];

	ierr = MatGetType(A, &mtype);
	PetscErrorCode (*set_stencil)(Mat mat, PetscInt m, const MatStencil idxm[], PetscInt n,
	                              const MatStencil idxn[], const PetscScalar v[],
	                              InsertMode addv);
	ierr = PetscStrcmp(mtype, MATSHELL, &is_boxmg);CHKERRQ(ierr);
	if (is_boxmg) {
		set_stencil = &stella_bmg_SetValuesStencil;
	} else {
		set_stencil = &MatSetValuesStencil;
	}

	stella_classify *cls = bc->level->classify;
	char ***classify;
	double ***bcvals;
	char dirichlet = cls->ptypes.dirichlet;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);
	ierr = stella_dmap_get(bc->slv_dmap, bc->values, &bcvals);CHKERRQ(ierr);


	met = bc->level->metric;
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	PetscInt xs, ys, zs, xm, ym, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->ladd_cont, &acont);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);

	enum dir {
		O = 0, N, S, E, W, F, B, NW, NE, SW, SE, NF, NB, SF, SB, FE, FW, BE, BW
	};

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	int il, jl, ir, jr, kl, kr;

	v[0] = 1.0;
	for (k = zs; k < zs + zm; k++) {
		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				if (classify[k][j][i] == dirichlet) {

					il = i>0     ? i-1 : 0;
					jl = j>0     ? j-1 : 0;
					kl = k>0     ? k-1 : 0;
					ir = i<ngx-1 ? i+1 : ngx-1;
					jr = j<ngy-1 ? j+1 : ngy-1;
					kr = k<ngz-1 ? k+1 : ngz-1;

					dcoefh[N] = have(dcoef[k][j][i], dcoef[k][jr][i]);
					dcoefh[S] = have(dcoef[k][j][i], dcoef[k][jl][i]);
					dcoefh[W] = have(dcoef[k][j][i], dcoef[k][j][il]);
					dcoefh[E] = have(dcoef[k][j][i], dcoef[k][j][ir]);
					dcoefh[F] = have(dcoef[k][j][i], dcoef[kr][j][i]);
					dcoefh[B] = have(dcoef[k][j][i], dcoef[kl][j][i]);
					dcoefh[SW] = have(dcoefh[W],
					                  have(dcoef[k][jl][i], dcoef[k][jl][il]));
					dcoefh[NW] = have(dcoefh[W],
					                  have(dcoef[k][jr][i], dcoef[k][jr][il]));
					dcoefh[SE] = have(dcoefh[E],
					                  have(dcoef[k][jl][i], dcoef[k][jl][ir]));
					dcoefh[NE] = have(dcoefh[E],
					                  have(dcoef[k][jr][i], dcoef[k][jr][ir]));
					dcoefh[FE] = have(dcoefh[E],
					                  have(dcoef[kr][j][i], dcoef[kr][j][ir]));
					dcoefh[BE] = have(dcoefh[E],
					                  have(dcoef[kl][j][i], dcoef[kl][j][ir]));
					dcoefh[FW] = have(dcoefh[W],
					                  have(dcoef[kr][j][i], dcoef[kr][j][il]));
					dcoefh[BW] = have(dcoefh[W],
					                  have(dcoef[kl][j][i], dcoef[kl][j][il]));
					dcoefh[NF] = have(dcoefh[F],
					                  have(dcoef[k][jr][i], dcoef[kr][jr][i]));
					dcoefh[SF] = have(dcoefh[F],
					                  have(dcoef[k][jl][i], dcoef[kr][jl][i]));
					dcoefh[NB] = have(dcoefh[B],
					                  have(dcoef[k][jr][i], dcoef[kl][jr][i]));
					dcoefh[SB] = have(dcoefh[B],
					                  have(dcoef[k][jl][i], dcoef[kl][jl][i]));
					row.i = i; row.j = j; row.k = k;
					col[0].i = i; col[0].j = j; col[0].k = k;
					cnt = 1;

					if (i != 0) {
						col[cnt].i = i-1; col[cnt].j = j; col[cnt].k = k;
						v[cnt] = 0; cnt++;
						if (classify[k][j][i-1] != dirichlet) {
							double weight = coef[MET_W][k][j][i] * dcoefh[W] * bcvals[k][j][i];
							acont[k][j][i-1] += weight;
						}
					}
					if (j != 0) {
						col[cnt].i = i; col[cnt].j = j-1; col[cnt].k = k;
						v[cnt] = 0; cnt++;
						if (classify[k][j-1][i] != dirichlet) {
							double weight = coef[MET_S][k][j][i] * dcoefh[S] * bcvals[k][j][i];
							acont[k][j-1][i] += weight;
						}
					}
					if (k != 0) {
						col[cnt].i = i; col[cnt].j = j; col[cnt].k = k-1;
						v[cnt] = 0; cnt++;
						if (classify[k-1][j][i] != dirichlet) {
							double weight = coef[MET_B][k][j][i] * dcoefh[B] * bcvals[k][j][i];
							acont[k-1][j][i] += weight;
						}
					}
					if (i != ngx - 1) {
						col[cnt].i = i+1; col[cnt].j = j; col[cnt].k = k;
						v[cnt] = 0; cnt++;
						if (classify[k][j][i+1] != dirichlet) {
							double weight = coef[MET_W][k][j][i+1] * dcoefh[E] * bcvals[k][j][i];
							acont[k][j][i+1] += weight;
						}
					}
					if (j != ngy - 1) {
						col[cnt].i = i; col[cnt].j = j+1; col[cnt].k = k;
						v[cnt] = 0; cnt++;
						if (classify[k][j+1][i] != dirichlet) {
							double weight = coef[MET_S][k][j+1][i] * dcoefh[N] * bcvals[k][j][i];
							acont[k][j+1][i] += weight;
						}
					}
					if (k != ngz - 1) {
						col[cnt].i = i; col[cnt].j = j; col[cnt].k = k+1;
						v[cnt] = 0; cnt++;
						if (classify[k+1][j][i] != dirichlet) {
							double weight = coef[MET_B][k+1][j][i] * dcoefh[F] * bcvals[k][j][i];
							acont[k+1][j][i] += weight;
						}
					}
					if (i != 0 && j != 0) {
						col[cnt].i = i-1; col[cnt].j = j-1; col[cnt].k = k;
						v[cnt] = 0; cnt++;
						if ((classify[k][j][i-1] != dirichlet) || (classify[k][j-1][i] != dirichlet)) {
							double weight = coef[MET_SW][k][j][i] * dcoefh[SW] * bcvals[k][j][i];
							acont[k][j-1][i-1] += weight;
						}
					}
					if (i != ngx - 1 && j != 0) {
						col[cnt].i = i+1; col[cnt].j = j-1; col[cnt].k = k;
						v[cnt] = 0; cnt++;
						if ((classify[k][j][i+1] != dirichlet) || (classify[k][j-1][i] != dirichlet)) {
							double weight = coef[MET_SE][k][j][i] * dcoefh[SE] * bcvals[k][j][i];
							acont[k][j-1][i+1] += weight;
						}
					}
					if (i != 0 && j != ngy-1) {
						col[cnt].i = i-1; col[cnt].j = j+1; col[cnt].k = k;
						v[cnt] = 0; cnt++;
						if ((classify[k][j][i-1] != dirichlet) || (classify[k][j+1][i] != dirichlet)) {
							double weight = coef[MET_SE][k][j+1][i-1] * dcoefh[NW] * bcvals[k][j][i];
							acont[k][j+1][i-1] += weight;
						}
					}
					if ((i != ngx -1) && (j != ngy - 1)) {
						col[cnt].i = i+1; col[cnt].j = j+1; col[cnt].k = k;
						v[cnt] = 0; cnt++;
						if ((classify[k][j][i+1] != dirichlet) || (classify[k][j+1][i] != dirichlet)) {
							double weight = coef[MET_SW][k][j+1][i+1] * dcoefh[NW] * bcvals[k][j][i];
							acont[k][j+1][i+1] += weight;
						}
					}
					if (i != 0 && k != 0) {
						col[cnt].i = i-1; col[cnt].j = j; col[cnt].k = k-1;
						v[cnt] = 0; cnt++;
						if ((classify[k][j][i-1] != dirichlet) || (classify[k-1][j][i] != dirichlet)) {
							double weight = coef[MET_WB][k][j][i] * dcoefh[BW] * bcvals[k][j][i];
							acont[k-1][j][i-1] += weight;
						}
					}
					if (i != ngx - 1 && k != 0) {
						col[cnt].i = i+1; col[cnt].j = j; col[cnt].k = k-1;
						v[cnt] = 0; cnt++;
						if ((classify[k][j][i+1] != dirichlet) || (classify[k-1][j][i] != dirichlet)) {
							double weight = coef[MET_EB][k][j][i] * dcoefh[BE] * bcvals[k][j][i];
							acont[k-1][j][i+1] += weight;
						}
					}
					if (i != 0 && k != ngz-1) {
						col[cnt].i = i-1; col[cnt].j = j; col[cnt].k = k+1;
						v[cnt] = 0; cnt++;
						if ((classify[k][j][i-1] != dirichlet) || (classify[k+1][j][i] != dirichlet)) {
							double weight = coef[MET_EB][k+1][j][i-1]  * dcoefh[BE] * bcvals[k][j][i];
							acont[k+1][j][i-1] += weight;
						}
					}
					if ((i != ngx-1) && (k != ngz -1)) {
						col[cnt].i = i+1; col[cnt].j = j; col[cnt].k = k+1;
						v[cnt] = 0; cnt++;
						if ((classify[k][j][i+1] != dirichlet) || (classify[k+1][j][i] != dirichlet)) {
							double weight = coef[MET_WB][k+1][j][i+1] * dcoefh[FE] * bcvals[k][j][i];
							acont[k+1][j][i+1] += weight;
						}
					}
					if (j != 0 && k != 0) {
						col[cnt].i = i; col[cnt].j = j-1; col[cnt].k = k-1;
						v[cnt] = 0; cnt++;
						if ((classify[k][j-1][i] != dirichlet) || (classify[k-1][j][i] != dirichlet)) {
							double weight = coef[MET_SB][k][j][i] * dcoefh[SB] * bcvals[k][j][i];
							acont[k-1][j-1][i] += weight;
						}
					}
					if (j != ngy - 1 && k != 0) {
						col[cnt].i = i; col[cnt].j = j+1; col[cnt].k = k-1;
						v[cnt] = 0; cnt++;
						if ((classify[k][j+1][i] != dirichlet) || (classify[k-1][j][i] != dirichlet)) {
							double weight = coef[MET_NB][k][j][i] * dcoefh[NB] * bcvals[k][j][i];
							acont[k-1][j+1][i] += weight;
						}
					}
					if (j != 0 && k != ngz-1) {
						col[cnt].i = i; col[cnt].j = j-1; col[cnt].k = k+1;
						v[cnt] = 0; cnt++;
						if ((classify[k][j-1][i] != dirichlet) || (classify[k+1][j][i] != dirichlet)) {
							double weight = coef[MET_NB][k+1][j-1][i]  * dcoefh[SF] * bcvals[k][j][i];
							acont[k+1][j-1][i] += weight;
						}
					}
					if ((j != ngy-1) && (k != ngz -1)) {
						col[cnt].i = i; col[cnt].j = j+1; col[cnt].k = k+1;
						v[cnt] = 0; cnt++;
						if ((classify[k][j+1][i] != dirichlet) || (classify[k+1][j][i] != dirichlet)) {
							double weight = coef[MET_SB][k+1][j+1][i] * dcoefh[NF] * bcvals[k][j][i];
							acont[k+1][j+1][i] += weight;
						}
					}
					ierr = set_stencil(A, 1, &row, cnt, col, v, INSERT_VALUES);CHKERRQ(ierr);
				}
			}
		}
	}

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	ierr = stella_dmap_restore(bc->slv_dmap, &bcvals);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->ladd_cont, &acont);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);

	return 0;
}




static PetscErrorCode apply_dirichlet_rhs(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	PetscInt i,j;
	PetscScalar **rhs_vec, **pm;
	double **bcvals;
	stella_classify *cls = bc->level->classify;
	char **classify;
	char dirichlet = cls->ptypes.dirichlet;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);

	PetscInt xs, ys, xm, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);

	ierr = stella_dmap_get(bc->slv_dmap, bc->values, &bcvals);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, rhs, &rhs_vec);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->pm, &pm);CHKERRQ(ierr);
	for (j = ys; j < ys + ym; j++) {
		for (i = xs; i < xs + xm; i++) {
			if (classify[j][i] == dirichlet) {
				rhs_vec[j][i] = bcvals[j][i];
				pm[j][i] = 0.0;
			}
		}
	}
	ierr = stella_dmap_restore(bc->slv_dmap, &bcvals);CHKERRQ(ierr);
	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, rhs, &rhs_vec);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->pm, &pm);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode apply_dirichlet_rhs_3d(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	PetscInt i,j,k;
	PetscScalar ***rhs_vec, ***pm;
	double ***bcvals;
	stella_classify *cls = bc->level->classify;
	char ***classify;
	char dirichlet = cls->ptypes.dirichlet;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);

	PetscInt xs, ys, zs, xm, ym, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	ierr = stella_dmap_get(bc->slv_dmap, bc->values, &bcvals);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, rhs, &rhs_vec);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->pm, &pm);CHKERRQ(ierr);
	for (k = zs; k < zs + zm; k++) {
		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				if (classify[k][j][i] == dirichlet) {
					rhs_vec[k][j][i] = bcvals[k][j][i];
					pm[k][j][i] = 0.0;
				}
			}
		}
	}
	ierr = stella_dmap_restore(bc->slv_dmap, &bcvals);CHKERRQ(ierr);
	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, rhs, &rhs_vec);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->pm, &pm);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode apply_op(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt dim;

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	if (dim == 2) {
		ierr = apply_dirichlet(bc, A, da);CHKERRQ(ierr);
	} else if (dim == 3) {
		ierr = apply_dirichlet_3d(bc, A, da);CHKERRQ(ierr);
	}

	return 0;
}


static PetscErrorCode symmetric(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt dim;

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	if (dim == 2) {
		ierr = update_dirichlet_sym(bc, A, da);CHKERRQ(ierr);
	} else if (dim == 3) {
		ierr = update_dirichlet_sym_3d(bc, A, da);CHKERRQ(ierr);
	}

	return 0;
}


static PetscErrorCode apply_rhs(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	PetscInt dim;

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	if (dim == 2) {
		ierr = apply_dirichlet_rhs(bc, da, rhs);CHKERRQ(ierr);
	} else if (dim == 3) {
		ierr = apply_dirichlet_rhs_3d(bc, da, rhs);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_dirichlet_create(stella_bc *bc)
{
	bc->apply = &apply_op;
	bc->apply_rhs = &apply_rhs;
	bc->symmetric = &symmetric;

	return 0;
}
