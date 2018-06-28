#include <stdlib.h>
#include "state.h"

state *state_create(grid *grd, problem *pb)
{
	int i, j, k;
	int ind;

	state *st = (state*) malloc(sizeof(state));
	st->rhs = (double*) malloc(grd->num_pts*sizeof(double));
	st->phi = (double*) malloc(grd->num_pts*sizeof(double));
	st->eps = (double*) malloc(grd->num_pts*sizeof(double));
	st->debye = (double*) malloc(grd->num_pts*sizeof(double));
	st->sol = (double*) malloc(grd->num_pts*sizeof(double));
	st->jc = (double*) malloc(grd->num_pts*sizeof(double));

	for (i = 0; i < grd->num_pts; i++) {
		st->debye[i] = 0.0;
	}

	if (pb->id == JUMP) {
		if (grd->nd == 2) {
			for (j = 0; j < grd->len[1]; j++) {
				for (i = 0; i < grd->len[0]; i++) {
					ind = j*grd->len[0] + i;
					if (grd->x[ind] <= 0)
						st->eps[ind] = 4;
					else
						st->eps[ind] = 2;
				}
			}
		} else {
			for (k = 0; k < grd->len[2]; k++) {
				for (j = 0; j < grd->len[1]; j++) {
					for (i = 0; i < grd->len[0]; i++) {
						ind = k*grd->len[0]*grd->len[1] + j*grd->len[0] + i;
						if (grd->x[ind] <= 0)
							st->eps[ind] = 4;
						else
							st->eps[ind] = 2;
					}
				}
			}
		}
		for (i = 0; i < grd->num_pts; i++) {
			st->jc[i] = -5;
		}
	} else {
		for (i = 0; i < grd->num_pts; i++) {
			st->eps[i] = 1.0;
			st->jc[i] = 0;
		}
	}

	grid_eval(grd, pb->rhs, st->rhs);
	grid_eval(grd, pb->sol, st->sol);

	return st;
}


void state_destroy(state* st)
{
	free(st->rhs);
	free(st->phi);
	free(st->eps);
	free(st->debye);
}
