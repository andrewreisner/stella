#include <iostream>
#include <vector>
#include <array>

#include <mpi.h>
#include "gtest/gtest.h"

#include "stella_signals.h"
#include "problem.h"
#include "map.h"
#include "grid.h"
#include "solver.h"


namespace {

static double norm(grid *grd, double *e_sol, double *c_sol)
{
	double *diff = new double[grd->num_pts];

	int i;
	for (i = 0; i < grd->num_pts; i++)
		diff[i] = 0;

	if (grd->nd == 2) {
		int i, j;
		for (j = grd->ibeg[1]; j <= grd->iend[1]; j++) {
			for (i = grd->ibeg[0]; i <= grd->iend[0]; i++) {
				int ind = j*grd->len[0] + i;
				diff[ind] = e_sol[ind] - c_sol[ind];
			}
		}
	} else {
		int i, j, k;
		for (k = grd->ibeg[2]; k <= grd->iend[2]; k++) {
			for (j = grd->ibeg[1]; j <= grd->iend[1]; j++) {
				for (i = grd->ibeg[0]; i <= grd->iend[0]; i++) {
					int ind = k*grd->len[1]*grd->len[0] + j*grd->len[0] + i;
					diff[ind] = e_sol[ind] - c_sol[ind];
				}
			}
		}
	}

	double curr = 0;
	for (i = 0; i < grd->num_pts; i++) {
		if (fabs(diff[i]) > curr) curr = fabs(diff[i]);
	}

	delete[] diff;

	return curr;
}


static double mpi_norm(grid *grd, double *e_sol, double *c_sol)
{
	int rank;
	double n;

	double mynorm = norm(grd, e_sol, c_sol);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Reduce(&mynorm, &n, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	if (rank == 0) return n;
	else return 0;
}


TEST(Signals, dcoef2D) {
	// expected norm for correct dcoef
	PetscErrorCode ierr;
	int nx = 200;
	int i, j, ind;
	grid *grd = grid_create(-1, 1, nx,
	                        -1, 1, nx,
	                         0, 0, 0);
	problem *pb = problem_create(JSINE, 2, 0);
	solver *sol = solver_create(grd, pb);

	double *dcoef = sol->state->eps;

	ierr = solver_init(sol, grd);
	ASSERT_EQ(ierr, 0);
	ierr = solver_run(sol);
	ASSERT_EQ(ierr, 0);

	double *u = new double[grd->num_pts];
	grid_eval(grd, pb->sol, u);

	double nrm_exp = mpi_norm(grd, sol->state->phi, u);

	//new grid for dcoef test
	grd = grid_create(-1, 1, nx,
	                  -1, 1, nx,
	                   0, 0, 0);
	pb = problem_create(JSINE, 2, 0);
	sol = solver_create(grd, pb);

	// inject bad dcoef
	dcoef = sol->state->eps;
	for (int i = 0; i < grd->num_pts; i++)
		dcoef[i] = 0;

	ierr = solver_init(sol, grd);
	ASSERT_EQ(ierr, 0);
	ierr = solver_run(sol);
	ASSERT_EQ(ierr, 0);

	u = new double[grd->num_pts];
	grid_eval(grd, pb->sol, u);
	{
		double nrm_bad = mpi_norm(grd, sol->state->phi, u);
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank > 0)
			nrm_bad = 1;
		ASSERT_GT(nrm_bad, nrm_exp);
	}

	// inject correct dcoef
	for (j = 0; j < grd->len[1]; j++) {
		for (i = 0; i < grd->len[0]; i++) {
			ind = j*grd->len[0] + i;
			if (grd->x[ind] <= 0)
				dcoef[ind] = 5;
			else
				dcoef[ind] = 2;
		}
	}

	ierr = stella_changed_dcoef(sol->ptr);
	ASSERT_EQ(ierr, 0);
	ierr = stella_solve(sol->ptr);
	ASSERT_EQ(ierr, 0);

	double nrm = mpi_norm(grd, sol->state->phi, u);
	ASSERT_EQ(nrm, nrm_exp);
}


TEST(Signals, dcoef3D) {
	// expected norm for correct dcoef
	PetscErrorCode ierr;
	int nx = 40;
	int i, j, k, ind;
	grid *grd = grid_create(-1, 1, nx,
	                        -1, 1, nx,
	                        -1, 1, nx);
	problem *pb = problem_create(JSINE, 3, 0);
	solver *sol = solver_create(grd, pb);

	double *dcoef = sol->state->eps;

	ierr = solver_init(sol, grd);
	ASSERT_EQ(ierr, 0);
	ierr = solver_run(sol);
	ASSERT_EQ(ierr, 0);

	double *u = new double[grd->num_pts];
	grid_eval(grd, pb->sol, u);

	double nrm_exp = mpi_norm(grd, sol->state->phi, u);

	//new grid for dcoef test
	grd = grid_create(-1, 1, nx,
	                  -1, 1, nx,
	                  -1, 1, nx);
	pb = problem_create(JSINE, 3, 0);
	sol = solver_create(grd, pb);

	// inject bad dcoef
	dcoef = sol->state->eps;
	for (int i = 0; i < grd->num_pts; i++)
		dcoef[i] = 0;

	ierr = solver_init(sol, grd);
	ASSERT_EQ(ierr, 0);
	ierr = solver_run(sol);
	ASSERT_EQ(ierr, 0);

	u = new double[grd->num_pts];
	grid_eval(grd, pb->sol, u);
	{
		double nrm_bad = mpi_norm(grd, sol->state->phi, u);
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank > 0)
			nrm_bad = 1;
		ASSERT_GT(nrm_bad, nrm_exp);
	}

	// inject correct dcoef
	for (k = 0; k < grd->len[2]; k++) {
		for (j = 0; j < grd->len[1]; j++) {
			for (i = 0; i < grd->len[0]; i++) {
				ind = k*grd->len[0]*grd->len[1] + j*grd->len[0] + i;
				if ((grd->x[ind] > 0 && grd->y[ind] >= 0 && grd->z[ind] >= 0) || (grd->x[ind] <= 0 && grd->y[ind] < 0 && grd->z[ind] >= 0))
			        dcoef[ind] = 3;
			    else if ((grd->x[ind] <= 0 && grd->y[ind] >= 0 && grd->z[ind] >= 0) || (grd->x[ind] > 0 && grd->y[ind] < 0 && grd->z[ind] >= 0))
			        dcoef[ind] = 1;
			    else if ((grd->x[ind] > 0 && grd->y[ind] >= 0 && grd->z[ind] < 0) || (grd->x[ind] <= 0 && grd->y[ind] < 0 && grd->z[ind] < 0))
			        dcoef[ind] = 7;
			    else
			        dcoef[ind] = 5;
			}
		}
	}

	ierr = stella_changed_dcoef(sol->ptr);
	ASSERT_EQ(ierr, 0);
	ierr = stella_solve(sol->ptr);
	ASSERT_EQ(ierr, 0);

	double nrm = mpi_norm(grd, sol->state->phi, u);
	ASSERT_EQ(nrm, nrm_exp);
}

}
