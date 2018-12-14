#include <stdlib.h>
#include <math.h>
#include "state.h"

state *state_create(grid *grd, problem *pb)
{
	int i;

	state *st = (state*) malloc(sizeof(state));
	st->rhs = (double*) malloc(grd->num_pts*sizeof(double));
	st->phi = (double*) malloc(grd->num_pts*sizeof(double));
	st->debye = (double*) malloc(grd->num_pts*sizeof(double));
	st->sol = (double*) malloc(grd->num_pts*sizeof(double));
	st->eps = (double*) malloc(grd->num_pts*sizeof(double));
	st->jump[0] = (double*) malloc(grd->num_pts*sizeof(double));
	st->jump[1] = (double*) malloc(grd->num_pts*sizeof(double));
	st->jump[2] = (double*) malloc(grd->num_pts*sizeof(double));

	for (i = 0; i < grd->num_pts; i++) {
		st->debye[i] = 0.0;
	}

	grid_eval(grd, pb->rhs, st->rhs);
	grid_eval(grd, pb->sol, st->sol);
	grid_eval(grd, pb->eps, st->eps);
	grid_eval(grd, pb->jc[0], st->jump[0]);
	grid_eval(grd, pb->jc[1], st->jump[1]);
	grid_eval(grd, pb->jc[2], st->jump[2]);

	return st;
}


void state_destroy(state* st)
{
	free(st->rhs);
	free(st->phi);
	free(st->eps);
	free(st->debye);
	free(st->jump[0]);
	free(st->jump[1]);
	free(st->jump[2]);
}
