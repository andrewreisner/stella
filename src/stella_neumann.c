#include "stella_neumann.h"
#include "stella_mat.h"
#include "stella_operator.h"

#include "assert.h"

/**
 * Second order discretization of Neumann boundary conditions.
 */
static PetscErrorCode apply_neumann(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,ngx,ngy;
	PetscScalar v[9];
	MatStencil row, col[9];
	PetscScalar **dcoef;
	PetscScalar **bcoef;
	PetscScalar **coef[5];
	DM cda;
	PetscScalar **rootg;
	DMDACoor2d **coors;
	Vec lc;

	MatType mtype;
	PetscBool is_cedar;
	double acc;
	stella_metric *met;

	stella_classify *cls = bc->level->classify;
	char **classify, **norm_dirs;
	char neumann = cls->ptypes.neumann;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);
	ierr = stella_dmap_get_char(bc->slv_dmap, bc->norm_dirs, &norm_dirs);CHKERRQ(ierr);

	PetscInt xs, ys, xm, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0,0,0,0,0);CHKERRQ(ierr);
	ierr = MatGetType(A, &mtype);CHKERRQ(ierr);

	met = bc->level->metric;

	ierr = DMDAVecGetArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->lbcoef, &bcoef);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, met->jac, &rootg);CHKERRQ(ierr);

	if (bc->axisymmetric) {
		ierr = DMGetCoordinateDM(da, &cda);CHKERRQ(ierr);
		ierr = DMGetCoordinatesLocal(da, &lc);CHKERRQ(ierr);
		ierr = DMDAVecGetArray(cda, lc, &coors);CHKERRQ(ierr);
	}

	for (i = 0; i < 5; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	// Stencil indexes
	PetscInt si, sj;
	PetscInt sis, sie, sjs, sje;

	// Left and right indexes. Used to avoid indexing out of the
	// array
	PetscInt il, ir;
	PetscInt jl, jr;

	// Harmonic averages of dcoef
	PetscScalar DCOEF[3][3];
	PetscInt ind;

	for (j = ys; j < ys + ym; j++) {
		for (i = xs; i < xs + xm; i++) {
			if (classify[j][i] == neumann) {
				int normdOut = -1 * norm_dirs[j][i]; // normal pointing out
				row.i = i; row.j = j;
				col[0].i = i; col[0].j = j;

				il = i>0     ? i-1 : 0;
				jl = j>0     ? j-1 : 0;
				ir = i<ngx-1 ? i+1 : ngx-1;
				jr = j<ngy-1 ? j+1 : ngy-1;

				DCOEF[0][0] = coef[MET_SW][j][i]*have(have(dcoef[j][i], dcoef[j][il]),
				                                      have(dcoef[jl][i], dcoef[jl][il]));
				DCOEF[0][1] = coef[MET_S][j][i]*have(dcoef[j][i],dcoef[jl][i]);
				DCOEF[0][2] = coef[MET_SE][j][i]*have(have(dcoef[j][i], dcoef[j][ir]),
				                                      have(dcoef[jl][i], dcoef[jl][ir]));
				DCOEF[1][0] = coef[MET_W][j][i]*have(dcoef[j][i],dcoef[j][il]);
				DCOEF[1][1] = 0.0;
				DCOEF[1][2] = coef[MET_W][j][ir]*have(dcoef[j][i],dcoef[j][ir]);
				DCOEF[2][0] = coef[MET_SE][jr][il]*have(have(dcoef[j][i], dcoef[j][il]),
				                                        have(dcoef[jr][i], dcoef[jr][il]));
				DCOEF[2][1] = coef[MET_S][jr][i]*have(dcoef[j][i],dcoef[jr][i]);
				DCOEF[2][2] = coef[MET_SW][jr][ir]*have(have(dcoef[j][i], dcoef[j][ir]),
				                                        have(dcoef[jr][i], dcoef[jr][ir]));

				if (bc->axisymmetric) {
					DCOEF[0][1] *= coors[j][i].x;
					DCOEF[2][1] *= coors[j][i].x;
					DCOEF[1][0] *= (coors[j][i].x + coors[j][i-1].x) * .5;
					DCOEF[1][2] *= (coors[j][i+1].x + coors[j][i].x) * .5;
				}

				sjs=-1; sje=+1;
				sis=-1; sie=+1;
				if(i==0)     sis = 0;
				if(i==ngx-1) sie = 0;
				if(j==0)     sjs = 0;
				if(j==ngy-1) sje = 0;

				// Handle corners
				if( ((i==0) || (i==ngx-1)) &&
				    ((j==0) || (j==ngy-1)) ) {

					DCOEF[0][0]=0.0;DCOEF[0][2]=0.0;
					DCOEF[2][0]=0.0;DCOEF[2][2]=0.0;
					// Scale corners by 1.0 for symmetry
					/* DCOEF[0][1] *= 2.0/2.0; */
					/* DCOEF[1][0] *= 2.0/2.0; */
					/* DCOEF[1][2] *= 2.0/2.0; */
					/* DCOEF[2][1] *= 2.0/2.0; */

				}
				// Handle edges
				else {

					if(abs(normdOut)==2)
						for(si=-1; si<=1; si++) {
							DCOEF[2][si+1] *= 2.0;
							DCOEF[0][si+1] *= 2.0;
						}
					else // abs(normdOut)==1
						for(sj=-1; sj<=1; sj++) {
							DCOEF[sj+1][0] *= 2.0;
							DCOEF[sj+1][2] *= 2.0;
						}

				}

				ind = 1;
				for(sj=sjs; sj<=sje; sj++) {
					for(si=sis; si<=sie; si++) {

						if((sj!=0) || (si!=0)) {
							v[ind] = -1*DCOEF[sj+1][si+1];
							col[ind].i = i+si;
							col[ind].j = j+sj;
							ind++;
						}

					}
				}

				acc = 0;
				for (k = 1; k < ind; k++){
					v[k] = .5 * v[k];
					acc += v[k];
				}
				v[0] = -1*acc;

				PetscScalar diag_cont;
				diag_cont = rootg[j][i]*bcoef[j][i];
				if (bc->axisymmetric) {
					diag_cont *= coors[j][i].x;
				}
				PetscScalar bscale = .5;
				if (ind == 4)
					bscale = .25;
				v[0] += PetscSign(v[0]) * bscale*PetscAbsScalar(diag_cont);

				ierr = PetscStrcmp(mtype, MATSHELL, &is_cedar);CHKERRQ(ierr);
				if (is_cedar) {
					ierr = stella_bmg_SetValuesStencil(A, 1, &row, ind, col, v, INSERT_VALUES);CHKERRQ(ierr);
				} else {
					ierr = MatSetValuesStencil(A, 1, &row, ind, col, v, INSERT_VALUES);CHKERRQ(ierr);
				}
			}

		}
	}

	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	ierr = stella_dmap_restore_char(bc->slv_dmap, &norm_dirs);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->lbcoef, &bcoef);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, met->jac, &rootg);CHKERRQ(ierr);
	for (i = 0; i < 5; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}

	if (bc->axisymmetric) {
		ierr = DMDAVecRestoreArray(cda, lc, &coors);CHKERRQ(ierr);
	}

	return 0;
}


/**
 * Second order discretization of Neumann boundary conditions.
 */
static PetscErrorCode apply_neumann_3d(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,ngx,ngy,ngz;
	PetscScalar v[27];
	PetscInt cnt;
	MatStencil row, col[27];
	MatType mtype;
	PetscBool is_cedar;
	double dcoefh[19];
	PetscScalar ***acont;
	stella_metric *met;
	PetscScalar ***dcoef;
	PetscScalar ***bcoef;
	PetscScalar ***rootg;
	PetscScalar ***coef[10];
	double acc, scale;
	DMBoundaryType bx, by, bz;
	int per_x, per_y, per_z;

	PetscInt xs, ys, zs, xm, ym, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	stella_classify *cls = bc->level->classify;
	char ***classify, ***norm_dirs;
	char neumann = cls->ptypes.neumann;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);
	ierr = stella_dmap_get_char(bc->slv_dmap, bc->norm_dirs, &norm_dirs);CHKERRQ(ierr);

	ierr = MatGetType(A, &mtype);
	PetscErrorCode (*set_stencil)(Mat mat, PetscInt m, const MatStencil idxm[], PetscInt n,
	                              const MatStencil idxn[], const PetscScalar v[],
	                              InsertMode addv);
	ierr = PetscStrcmp(mtype, MATSHELL, &is_cedar);CHKERRQ(ierr);
	if (is_cedar) {
		set_stencil = &stella_bmg_SetValuesStencil;
	} else {
		set_stencil = &MatSetValuesStencil;
	}

	met = bc->level->metric;
	enum dir {
		O = 0, N, S, E, W, F, B, NW, NE, SW, SE, NF, NB, SF, SB, FE, FW, BE, BW
	};

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecGetArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}
	ierr = DMDAVecGetArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, bc->level->lbcoef, &bcoef);CHKERRQ(ierr);
	ierr = DMDAVecGetArray(da, met->jac, &rootg);CHKERRQ(ierr);
	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0,&bx,&by,&bz,0);CHKERRQ(ierr);
	per_x = bx==DM_BOUNDARY_PERIODIC;
	per_y = by==DM_BOUNDARY_PERIODIC;
	per_z = bz==DM_BOUNDARY_PERIODIC;

	// Stencil indexes
	PetscInt si, sj, sk;
	PetscInt sis, sie;
	PetscInt sjs, sje;
	PetscInt sks, ske;

	// Left and right indexes. Used to avoid indexing out of the
	// array
	PetscInt il, ir;
	PetscInt jl, jr;
	PetscInt kl, kr;

	// Harmonic averages of dcoef
	PetscInt ind;
	PetscScalar DCOEF[3][3][3];
	for(sk=-1; sk<=1; sk++)
		for(sj=-1; sj<=1; sj++)
			for(si=-1; si<=1; si++)
				DCOEF[sk+1][sj+1][si+1] = 0.0;

	for (k = zs; k < zs + zm; k++) {
		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				if (classify[k][j][i] == neumann) {
					int normdOut = -1 * norm_dirs[k][j][i];
					il = (i>0 || per_x)    ? i-1 : 0;
					jl = (j>0 || per_y)    ? j-1 : 0;
					kl = (k>0 || per_z)    ? k-1 : 0;
					ir = (i<ngx-1 || per_x) ? i+1 : ngx-1;
					jr = (j<ngy-1 || per_y) ? j+1 : ngy-1;
					kr = (k<ngz-1 || per_z) ? k+1 : ngz-1;

					// N
					DCOEF[1][2][1]  = coef[MET_S][k][jr][i];
					dcoefh[N]       = have(dcoef[k][j][i], dcoef[k][jr][i]);
					DCOEF[1][2][1] *= dcoefh[N];
					// S
					DCOEF[1][0][1]  = coef[MET_S][k][j][i];
					dcoefh[S]       = have(dcoef[k][j][i], dcoef[k][jl][i]);
					DCOEF[1][0][1] *= dcoefh[S];
					// W
					DCOEF[1][1][0]  = coef[MET_W][k][j][i];
					dcoefh[W]       = have(dcoef[k][j][i], dcoef[k][j][il]);
					DCOEF[1][1][0] *= dcoefh[W];
					// E
					DCOEF[1][1][2]  = coef[MET_W][k][j][ir];
					dcoefh[E]       = have(dcoef[k][j][i], dcoef[k][j][ir]);
					DCOEF[1][1][2] *= dcoefh[E];
					// F
					DCOEF[2][1][1]  = coef[MET_B][kr][j][i];
					dcoefh[F]       = have(dcoef[k][j][i], dcoef[kr][j][i]);
					DCOEF[2][1][1] *= dcoefh[F];
					// B
					DCOEF[0][1][1]  = coef[MET_B][k][j][i];
					dcoefh[B]       = have(dcoef[k][j][i], dcoef[kl][j][i]);
					DCOEF[0][1][1] *= dcoefh[B];
					// SW
					DCOEF[1][0][0]  = coef[MET_SW][k][j][i];
					DCOEF[1][0][0] *= have(dcoefh[W], have(dcoef[k][jl][i], dcoef[k][jl][il]));
					// NW
					DCOEF[1][2][0]  = coef[MET_SE][k][jr][il];
					DCOEF[1][2][0] *= have(dcoefh[W], have(dcoef[k][jr][i], dcoef[k][jr][il]));
					// SE
					DCOEF[1][0][2]  = coef[MET_SE][k][j][i];
					DCOEF[1][0][2] *= have(dcoefh[E], have(dcoef[k][jl][i], dcoef[k][jl][ir]));
					// NE
					DCOEF[1][2][2]  = coef[MET_SW][k][jr][ir];
					DCOEF[1][2][2] *= have(dcoefh[E], have(dcoef[k][jr][i], dcoef[k][jr][ir]));
					// FE
					DCOEF[2][1][2]  = coef[MET_WB][kr][j][ir];
					DCOEF[2][1][2] *= have(dcoefh[E], have(dcoef[kr][j][i], dcoef[kr][j][ir]));
					// BE
					DCOEF[0][1][2]  = coef[MET_EB][k][j][i];
					DCOEF[0][1][2] *= have(dcoefh[E], have(dcoef[kl][j][i], dcoef[kl][j][ir]));
					// FW
					DCOEF[2][1][0]  = coef[MET_EB][kr][j][il];
					DCOEF[2][1][0] *= have(dcoefh[W], have(dcoef[kr][j][i], dcoef[kr][j][il]));
					// BW
					DCOEF[0][1][0]  = coef[MET_WB][k][j][i];
					DCOEF[0][1][0] *= have(dcoefh[W], have(dcoef[kl][j][i], dcoef[kl][j][il]));
					// NF
					DCOEF[2][2][1]  = coef[MET_SB][kr][jr][i];
					DCOEF[2][2][1] *= have(dcoefh[F], have(dcoef[k][jr][i], dcoef[kr][jr][i]));
					// SF
					DCOEF[2][0][1]  = coef[MET_NB][kr][jl][i];
					DCOEF[2][0][1] *= have(dcoefh[F], have(dcoef[k][jl][i], dcoef[kr][jl][i]));
					// NB
					DCOEF[0][2][1]  = coef[MET_NB][k][j][i];
					DCOEF[0][2][1] *= have(dcoefh[B], have(dcoef[k][jr][i], dcoef[kl][jr][i]));
					// SB
					DCOEF[0][0][1]  = coef[MET_SB][k][j][i];
					DCOEF[0][0][1] *= have(dcoefh[B], have(dcoef[k][jl][i], dcoef[kl][jl][i]));

					row.i = i; row.j = j; row.k = k;
					col[0].i = i; col[0].j = j; col[0].k = k;


					// Stencil start and end points
					sks=-1; ske=+1;
					sjs=-1; sje=+1;
					sis=-1; sie=+1;
					if(i==0 && !per_x)     sis = 0;
					if(i==ngx-1 && !per_x) sie = 0;
					if(j==0 && !per_y)     sjs = 0;
					if(j==ngy-1 && !per_y) sje = 0;
					if(k==0 && !per_z)     sks = 0;
					if(k==ngz-1 && !per_z) ske = 0;

					// Handle corners
					if( ((i==0 && !per_x) || (i==ngx-1 && !per_x)) &&
					    ((j==0 && !per_y) || (j==ngy-1 && !per_y)) &&
					    ((k==0 && !per_z) || (k==ngz-1 && !per_z))) {

						scale = 0.125;

						DCOEF[1][0][0]=0.0;DCOEF[1][0][2]=0.0;
						DCOEF[1][2][0]=0.0;DCOEF[1][2][2]=0.0;
						DCOEF[0][1][0]=0.0;DCOEF[0][1][2]=0.0;
						DCOEF[2][1][0]=0.0;DCOEF[2][1][2]=0.0;
						DCOEF[0][0][1]=0.0;DCOEF[0][2][1]=0.0;
						DCOEF[2][0][1]=0.0;DCOEF[2][2][1]=0.0;

						for(sk=0; sk<=2; sk++)
							for(sj=0; sj<=2; sj++)
								for(si=0; si<=2; si++)
									DCOEF[sk][sj][si] *= 2;

					}
					// Handle edges
					else if (( ((i==0 && !per_x) || (i==ngx-1 && !per_x)) &&
					           ((j==0 && !per_y) || (j==ngy-1 && !per_y)) )
					         ||
					         ( ((i==0 && !per_x) || (i==ngx-1 && !per_x)) &&
					           ((k==0 && !per_z) || (k==ngz-1 && !per_z)) )
					         ||
					         ( ((j==0 && !per_y) || (j==ngy-1 && !per_y)) &&
					           ((k==0 && !per_z) || (k==ngz-1 && !per_z)) )
						) {

						scale = 0.25;

						// Intersection of y-z faces
						if( ((j==0 && !per_y) || (j==ngy-1 && !per_y)) &&
						    ((k==0 && !per_z) || (k==ngz-1 && !per_z)) )
							for(si=-1; si<=1; si++) {
								DCOEF[0][2][si+1] *= 4.0;
								DCOEF[0][0][si+1] *= 4.0;
								DCOEF[2][2][si+1] *= 4.0;
								DCOEF[2][0][si+1] *= 4.0;

								DCOEF[1][2][si+1] *= 2.0;
								DCOEF[1][0][si+1] *= 2.0;
								DCOEF[2][1][si+1] *= 2.0;
								DCOEF[0][1][si+1] *= 2.0;
							}
						// Intersection of x-z faces
						else if( ((i==0 && !per_x) || (i==ngx-1 && !per_x)) &&
						         ((k==0 && !per_z) || (k==ngz-1 && !per_z)) )
							for(sj=-1; sj<=1; sj++) {
								DCOEF[0][sj+1][0] *= 4.0;
								DCOEF[0][sj+1][2] *= 4.0;
								DCOEF[2][sj+1][0] *= 4.0;
								DCOEF[2][sj+1][2] *= 4.0;

								DCOEF[2][sj+1][1] *= 2.0;
								DCOEF[0][sj+1][1] *= 2.0;
								DCOEF[1][sj+1][0] *= 2.0;
								DCOEF[1][sj+1][2] *= 2.0;
							}

						// Intersection of x-y faces
						else
							for(sk=-1; sk<=1; sk++) {
								DCOEF[sk+1][0][0] *= 4.0;
								DCOEF[sk+1][0][2] *= 4.0;
								DCOEF[sk+1][2][0] *= 4.0;
								DCOEF[sk+1][2][2] *= 4.0;

								DCOEF[sk+1][1][0] *= 2.0;
								DCOEF[sk+1][1][2] *= 2.0;
								DCOEF[sk+1][0][1] *= 2.0;
								DCOEF[sk+1][2][1] *= 2.0;
							}

					}
					// Handle faces
					else {

						scale = .5;

						if(abs(normdOut)==3) {
							for(sj=-1; sj<=1; sj++)
								for(si=-1; si<=1; si++) {
									DCOEF[0][sj+1][si+1] *= 2.0;
									DCOEF[2][sj+1][si+1] *= 2.0;
								}
						}
						else if(abs(normdOut)==2)
							for(sk=-1; sk<=1; sk++)
								for(si=-1; si<=1; si++) {
									DCOEF[sk+1][0][si+1] *= 2.0;
									DCOEF[sk+1][2][si+1] *= 2.0;
								}
						else // abs(normdOut)==1
							for(sk=-1; sk<=1; sk++)
								for(sj=-1; sj<=1; sj++) {
									DCOEF[sk+1][sj+1][0] *= 2.0;
									DCOEF[sk+1][sj+1][2] *= 2.0;
								}

					}

					ind = 1;
					for(sk=sks; sk<=ske; sk++)
						for(sj=sjs; sj<=sje; sj++)
							for(si=sis; si<=sie; si++)
								if((sk!=0) || (sj!=0) || (si!=0)) {
									v[ind] = -1*DCOEF[sk+1][sj+1][si+1];
									col[ind].i = i+si;
									col[ind].j = j+sj;
									col[ind].k = k+sk;
									ind++;
								}

					acc = 0;
					for (sk=1; sk<ind; sk++){
						v[sk] = scale*v[sk];
						acc += v[sk];
					}
					v[0] = -1*acc;

					PetscScalar diag_cont;
					diag_cont = rootg[k][j][i]*bcoef[k][j][i];
					PetscScalar bscale = .5;
					if (ind == 12) {
						bscale = .25;
					} else if (ind == 8) {
						bscale = 1.0/6.0;
					}
					v[0] += PetscSign(v[0]) * bscale*PetscAbsScalar(diag_cont);

					ierr = set_stencil(A, 1, &row, ind, col, v, INSERT_VALUES);CHKERRQ(ierr);
				}
			}
		}
	}

	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	ierr = stella_dmap_restore_char(bc->slv_dmap, &norm_dirs);CHKERRQ(ierr);

	for (i = 0; i < 10; i++) {
		ierr = DMDAVecRestoreArray(da, met->lcoef[i], &coef[i]);CHKERRQ(ierr);
	}
	ierr = DMDAVecRestoreArray(da, bc->level->ldcoef, &dcoef);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->lbcoef, &bcoef);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, met->jac, &rootg);CHKERRQ(ierr);

	return 0;
}


PetscErrorCode apply_neumann_rhs(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	PetscInt i,j,ngx,ngy;
	PetscScalar **scale;
	DMBoundaryType bx, by;
	int per_x, per_y;

	stella_classify *cls = bc->level->classify;
	char **classify;
	char neumann = cls->ptypes.neumann;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);

	PetscInt xs, ys, xm, ym;
	ierr = DMDAGetCorners(da, &xs, &ys, 0, &xm, &ym, 0);CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, 0,0,0,0,0,0,&bx,&by,0,0);CHKERRQ(ierr);
	per_x = bx==DM_BOUNDARY_PERIODIC;
	per_y = by==DM_BOUNDARY_PERIODIC;

	ierr = DMDAVecGetArray(da, bc->level->nscale, &scale);CHKERRQ(ierr);
	for (j = ys; j < ys + ym; j++) {
		for (i = xs; i < xs + xm; i++) {
			if (classify[j][i] == neumann) {
				if( ((i==0 && !per_x) || (i==ngx-1 && !per_x)) &&
				    ((j==0 && !per_y) || (j==ngy-1 && !per_y)) )
					scale[j][i] = .25;
				else
					scale[j][i] = .5;
			}
		}
	}

	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);
	ierr = DMDAVecRestoreArray(da, bc->level->nscale, &scale);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode apply_neumann_rhs_3d(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;
	PetscInt i,j,k,ngx,ngy,ngz;
	PetscScalar ***scale;
	int nb;
	DMBoundaryType bx, by, bz;
	int per_x, per_y, per_z;

	PetscInt xs, ys, zs, xm, ym, zm;
	ierr = DMDAGetCorners(da, &xs, &ys, &zs, &xm, &ym, &zm);CHKERRQ(ierr);

	stella_classify *cls = bc->level->classify;
	char ***classify;
	char neumann = cls->ptypes.neumann;
	ierr = stella_classify_get(cls, &classify);CHKERRQ(ierr);

	ierr = DMDAGetInfo(da, 0, &ngx, &ngy, &ngz,0,0,0,0,0,&bx,&by,&bz,0);CHKERRQ(ierr);
	per_x = bx==DM_BOUNDARY_PERIODIC;
	per_y = by==DM_BOUNDARY_PERIODIC;
	per_z = bz==DM_BOUNDARY_PERIODIC;

	ierr = DMDAVecGetArray(da, bc->level->nscale, &scale);CHKERRQ(ierr);

	for (k = zs; k < zs + zm; k++) {
		for (j = ys; j < ys + ym; j++) {
			for (i = xs; i < xs + xm; i++) {
				if (classify[k][j][i] == neumann) {
					// Corners
					if( ((i==0 && !per_x) || (i==ngx-1 && !per_x)) &&
					    ((j==0 && !per_y) || (j==ngy-1 && !per_y)) &&
					    ((k==0 && !per_z) || (k==ngz-1 && !per_z)) )
						scale[k][j][i] = 0.125;
					// Edges
					else if (( ((i==0 && !per_x) || (i==ngx-1 && !per_x)) &&
					           ((j==0 && !per_y) || (j==ngy-1 && !per_y)) )
					         ||
					         ( ((i==0 && !per_x) || (i==ngx-1 && !per_x)) &&
					           ((k==0 && !per_z) || (k==ngz-1 && !per_z)) )
					         ||
					         ( ((j==0 && !per_y) || (j==ngy-1 && !per_y)) &&
					           ((k==0 && !per_z) || (k==ngz-1 && !per_z)) )
						)
						scale[k][j][i] = 0.25;
					// Faces
					else
						scale[k][j][i] = 0.5;

				}
			}
		}
	}

	ierr = DMDAVecRestoreArray(da, bc->level->nscale, &scale);CHKERRQ(ierr);
	ierr = stella_classify_restore(cls, &classify);CHKERRQ(ierr);

	return 0;
}


static PetscErrorCode apply_op(stella_bc *bc, Mat A, DM da)
{
	PetscErrorCode ierr;
	PetscInt dim;

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	if (dim == 2) {
		ierr = apply_neumann(bc, A, da);CHKERRQ(ierr);
	} else if (dim == 3) {
		ierr = apply_neumann_3d(bc, A, da);CHKERRQ(ierr);
	}

	return 0;
}


static PetscErrorCode apply_rhs(stella_bc *bc, DM da, Vec rhs)
{
	PetscErrorCode ierr;

	PetscInt dim;

	ierr = DMDAGetInfo(da, &dim, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);CHKERRQ(ierr);

	if (dim == 2) {
		ierr = apply_neumann_rhs(bc, da, rhs);CHKERRQ(ierr);
	} else {
		ierr = apply_neumann_rhs_3d(bc, da, rhs);CHKERRQ(ierr);
	}

	return 0;
}


PetscErrorCode stella_neumann_create(stella_bc *bc)
{
	bc->apply = &apply_op;
	bc->apply_rhs = &apply_rhs;

	return 0;
}
