#include <stdlib.h>
#include <stdio.h>

#include "boundary.h"

#ifndef max
#define max(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

static int contains(int *proc_is, int *proc_ie, int *hole_is, int *hole_ie, int nd)
{
	int flag = 1;
	int i;
	for (i = 0; i < nd; i++) {
		flag = flag && (hole_is[i] < proc_ie[i]);
		flag = flag && (hole_ie[i] >= proc_is[i]);
	}

	return flag;
}


static void add_2d(boundary *bnd, grid *grd, problem *pb)
{
	char *classify = bnd->classify;
	char *norm_dirs = bnd->norm_dir;
	double *values = bnd->values;

	int i;
	for (i = 0; i < 4; i++) {
		int stride;
		int offset;
		int len;
		char norm_dir;

		len = 0;
		if (i == NORTH
		    && grd->cart_coord[1] == grd->num_procs[1]-1
		    && (!grd->periodic[1])) {
			stride = 1;
			offset = (grd->ibeg[1] + (grd->ny-1))*grd->len[0];
			len = grd->len[0];
			norm_dir = -2;
		} else if (i == SOUTH
		           && grd->cart_coord[1] == 0
		           && (!grd->periodic[1])) {
			stride = 1;
			offset = grd->ibeg[1]*grd->len[0];
			len = grd->len[0];
			norm_dir = 2;
		} else if (i == EAST
		           && grd->cart_coord[0] == grd->num_procs[0]-1
		           && (!grd->periodic[0])) {
			stride = grd->len[0];
			offset = grd->ibeg[0] + grd->nx - 1;
			len = grd->len[1];
			norm_dir = -1;
		} else if (i == WEST
		           && grd->cart_coord[0] == 0
		           && (!grd->periodic[0])) {
			stride = grd->len[0];
			offset = grd->ibeg[0];
			len = grd->len[1];
			norm_dir = 1;
		}

		int j;
		for (j = 0; j < len; j++) {
			int ind = j*stride + offset;
			// Give Dirichlet bc precedence
			if (!((classify[ind] == DIRICHLET) && (pb->boundary[i] != DIRICHLET))) {
				classify[ind] = pb->boundary[i];
				if (pb->boundary[i] == DIRICHLET)
					values[ind] = pb->sol(grd->x[ind], grd->y[ind], 0);
				else
					values[ind] = 0.0;
				norm_dirs[ind] = norm_dir;
			}
		}
	}
}


static void add_3d(boundary *bnd, grid *grd, problem *pb)
{
	char *classify = bnd->classify;
	char *norm_dirs = bnd->norm_dir;
	double *values = bnd->values;

	int i;
	for (i = 0; i < 6; i++) {
		int ii, jj;
		int stride[2];
		int len[2];
		int offset;
		char norm_dir;

		len[0] = 0;
		len[1] = 0;
		if (i == NORTH
		    && grd->cart_coord[1] == grd->num_procs[1]-1
		    && (!grd->periodic[1])) {
			stride[0] = grd->len[0]*grd->len[1];
			stride[1] = 1;
			offset = (grd->ibeg[1] + (grd->ny-1))*grd->len[0];
			norm_dir = -2;
			len[0] = grd->len[2];
			len[1] = grd->len[0];
		} else if (i == SOUTH
		           && grd->cart_coord[1] == 0
		           && (!grd->periodic[1])) {
			stride[0] = grd->len[0]*grd->len[1];
			stride[1] = 1;
			offset = grd->ibeg[1]*grd->len[0];
			norm_dir = 2;
			len[0] = grd->len[2];
			len[1] = grd->len[0];
		} else if (i == EAST
		           && grd->cart_coord[0] == grd->num_procs[0]-1
		           && (!grd->periodic[0])) {
			stride[0] = grd->len[0]*grd->len[1];
			stride[1] = grd->len[0];
			offset = grd->ibeg[0] + grd->nx - 1;
			norm_dir = -1;
			len[0] = grd->len[2];
			len[1] = grd->len[1];
		} else if (i == WEST
		           && grd->cart_coord[0] == 0
		           && (!grd->periodic[0])) {
			stride[0] = grd->len[0]*grd->len[1];
			stride[1] = grd->len[0];
			offset = grd->ibeg[0];
			norm_dir = 1;
			len[0] = grd->len[2];
			len[1] = grd->len[1];
		} else if (i == FRONT
		           && grd->cart_coord[2] == grd->num_procs[2] - 1
		           && (!grd->periodic[2])) {
			stride[0] = grd->len[0];
			stride[1] = 1;
			offset = (grd->ibeg[2] + (grd->nz-1))*grd->len[0]*grd->len[1];
			norm_dir = -3;
			len[0] = grd->len[1];
			len[1] = grd->len[2];
		} else if (i == BACK
		           && grd->cart_coord[2] == 0
		           && (!grd->periodic[2])) {
			stride[0] = grd->len[0];
			stride[1] = 1;
			offset = grd->ibeg[2] * grd->len[0] * grd->len[1];
			norm_dir = 3;
			len[0] = grd->len[1];
			len[1] = grd->len[0];
		}

		for (ii = 0; ii < len[0]; ii++) {
			for (jj = 0; jj < len[1]; jj++) {
				int ind = ii*stride[0] + jj*stride[1] + offset;
				// Give Dirichlet bc precedence
				if (!((classify[ind] == DIRICHLET) && (pb->boundary[i] != DIRICHLET))) {
					classify[ind] = pb->boundary[i];
					if (pb->boundary[i] == DIRICHLET)
						values[ind] = pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
					else
						values[ind] = 0.0;
					norm_dirs[ind] = norm_dir;
				}
			}
		}
	}
}


static void add_electrodes(boundary *bnd, grid* grd, problem *pb)
{
	char *classify = bnd->classify;
	char *norm_dir = bnd->norm_dir;
	double *values = bnd->values;
	int i;
	for (i = 0; i < pb->nholes; i++) {
		int bnd_is[3], el_is[3];
		int bnd_ie[3], el_ie[3];

		int k;
		for (k = 0; k < grd->nd; k++) {
			double esize = pb->holes[i].rel_size[k] * grd->num_global[k];
			double offset = pb->holes[i].rel_offset[k] * grd->num_global[k];
			el_is[k] = offset - esize * 0.5;
			el_ie[k] = offset + esize * 0.5;
		}

		if (contains(grd->is, grd->ie,
		             el_is, el_ie, grd->nd)) {
			int ne[3];
			ne[2] = 1;
			for (k = 0; k < grd->nd; k++) {
				bnd_is[k] = max(el_is[k], grd->is[k]);
				bnd_ie[k] = min(el_ie[k], grd->ie[k]);
				ne[k] = bnd_ie[k] - bnd_is[k] + 1;
			}

			printf("Electrode: (%d %d %d) -> (%d %d %d)\n",
			       bnd_is[0], bnd_is[1], bnd_is[2],
			       bnd_ie[0], bnd_ie[1], bnd_ie[2]);

			if (grd->nd == 2) {
				int ii, jj;
				for (jj = 0; jj < ne[1]; jj++) {
					for (ii = 0; ii < ne[0]; ii++) {
						int ind = (bnd_is[1] - grd->is[1] + jj)*grd->nx + (bnd_is[0] - grd->is[0] + ii);
						values[ind] = pb->sol(grd->x[ind], grd->y[ind], 0);
						classify[ind] = DIRICHLET;
						norm_dir[ind] = 0.0;
					}
				}
			} else {
				int ii,jj,kk;
				for (kk = 0; kk < ne[2]; kk++) {
					for (jj = 0; jj < ne[1]; jj++) {
						for (ii = 0; ii < ne[0]; ii++) {
							int ind = (((bnd_is[2] - grd->is[2] + kk)*grd->nx*grd->ny) +
							           ((bnd_is[1] - grd->is[1] + jj)*grd->nx) +
							           ((bnd_is[0] - grd->is[0] + ii)));
							values[ind] = pb->sol(grd->x[ind], grd->y[ind], grd->z[ind]);
							classify[ind] = DIRICHLET;
							norm_dir[ind] = 0.0;
						}
					}
				}
			}
		}
	}
}


boundary *boundary_create(grid *grd, problem *pb)
{
	int i, j;

	boundary *bnd = (boundary*) malloc(sizeof(boundary));

	bnd->classify = (char*) malloc(grd->num_pts*sizeof(char));
	bnd->norm_dir = (char*) malloc(grd->num_pts*sizeof(char));
	bnd->values = (double*) malloc(grd->num_pts*sizeof(double));

	for (i = 0; i < grd->num_pts; i++)
		bnd->classify[i] = 0;

	if (grd->nd == 3) {
		add_3d(bnd, grd, pb);
	} else {
		add_2d(bnd, grd, pb);
	}

	if (pb->id == ELECTRODE) {
		add_electrodes(bnd, grd, pb);
	}

	return bnd;
}


void boundary_destroy(boundary *bnd)
{
	free(bnd->classify);
	free(bnd->norm_dir);
	free(bnd->values);
}
