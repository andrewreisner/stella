#include <stdlib.h>
#include <math.h>
#include "state.h"

state *state_create(grid *grd, problem *pb)
{
	int i, j, k;
	int ind, indpx, indpy, indmx, indmy;

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
	} else if (pb->id == CBOARD) {
		if (grd->nd == 2) {
			for (j = 0; j < grd->len[1]; j++) {
				for (i = 0; i < grd->len[0]; i++) {
					ind = j*grd->len[0] + i;
					if ((grd->x[ind] <= 0 && grd->y[ind] < 0) || (grd->x[ind] > 0 && grd->y[ind] >= 0))
						st->eps[ind] = 2;
					else
						st->eps[ind] = 4;
				}
			}
		} else {
			for (k = 0; k < grd->len[2]; k++) {
				for (j = 0; j < grd->len[1]; j++) {
					for (i = 0; i < grd->len[0]; i++) {
						ind = k*grd->len[0]*grd->len[1] + j*grd->len[0] + i;
						if ((grd->x[ind] <= 0 && grd->y[ind] < 0) || (grd->x[ind] > 0 && grd->y[ind] >= 0))
							st->eps[ind] = 4;
						else
							st->eps[ind] = 2;
					}
				}
			}
		}
		if (grd->nd == 2) {
			for (j = 0; j < grd->len[1]; j++) {
				for (i = 0; i < grd->len[0]; i++) {
					ind = j*grd->len[0] + i;
					indpy = (j+1)*grd->len[0] + i;
					indpx = j*grd->len[0] + (i+1);
					indmy = (j-1)*grd->len[0] + i;
					indmx = j*grd->len[0] + (i-1);
					if (grd->x[ind] <= 0 && grd->y[ind] >= 0)
						st->jc[ind] = -5;
					else if (grd->x[ind] > 0 && grd->y[ind] >= 0)
						st->jc[ind] = -5;
					else if (grd->x[ind] < 0 && grd->y[ind] < 0)
						st->jc[ind] = 9.5;
					else if (grd->x[ind] >= 0 && grd->y[ind] < 0)
						st->jc[ind] = 9.5;

					//borders
					if ((grd->x[ind] < 0) && (j != grd->len[0]-1) && (st->eps[ind] != st->eps[indpy]))
						st->jc[ind] = 5;
					if ((grd->x[ind] < 0) && (j != 0) && (st->eps[ind] != st->eps[indmy]))
						st->jc[ind] = 5;
					if ((grd->x[ind] >= 0) && (j != grd->len[0]-1) && (st->eps[ind] != st->eps[indpy]))
						st->jc[ind] = -9.5;
					if ((grd->x[ind] >= 0) && (j != 0) && (st->eps[ind] != st->eps[indmy]))
						st->jc[ind] = -9.5;
				}
			}

		} else {
			for (k = 0; k < grd->len[2]; k++) {
				for (j = 0; j < grd->len[1]; j++) {
					for (i = 0; i < grd->len[0]; i++) {
						ind = k*grd->len[0]*grd->len[1] + j*grd->len[0] + i;
						if ((grd->x[ind] <= 0 && grd->y[ind] < 0) || (grd->x[ind] > 0 && grd->y[ind] >= 0))
							st->jc[ind] = -5;
						else if (grd->x[ind] <= 0 && grd->y[ind] >= 0)
							st->jc[ind] = 5;
						else
							st->jc[ind] = -10;
					}
				}
			}
		}
	} else if (pb->id == JSINE) {
		if (grd->nd == 2) {
			for (j = 0; j < grd->len[1]; j++) {
				for (i = 0; i < grd->len[0]; i++) {
					ind = j*grd->len[0] + i;
					if (grd->x[ind] <= 0)
						st->eps[ind] = 5;
					else
						st->eps[ind] = 2;
				}
			}
		}
		for (i = 0; i < grd->num_pts; i++) {
			st->jc[i] = -sin(2*M_PI*grd->y[i])*M_PI;
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
