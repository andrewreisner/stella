#include <limits.h>
#include <math.h>
#include <cmath>
#include <mpi.h>
#include <vector>
#include <array>
#include "gtest/gtest.h"

#include "problem.h"
#include "map.h"
#include "grid.h"
#include "solver.h"


namespace {


static double estimate_order(std::vector<double> & hs, std::vector<double> & norms)
{
	return std::log(norms[1] / norms[0]) / std::log(hs[1] / hs[0]);
}

static double norm(grid *grd, double *e_sol, double *c_sol)
{
	double *diff = new double[grd->num_pts];

	int i;
	for (i = 0; i < grd->num_pts; i++)
		diff[i] = 0;

	{
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


TEST(Models3, UniformSine)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 20;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);
		problem *pb = problem_create(SIN, 3, 0);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		double nrm = mpi_norm(grd, sol->state->phi, u);
		norms.push_back(nrm);

		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}


TEST(Models3, UniformTSine)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 53;
	double tol = 1e-3;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);
		problem *pb = problem_create(TSINE, 3, 0);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		double nrm = mpi_norm(grd, sol->state->phi, u);
		norms.push_back(nrm);

		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}


TEST(Models3, UniformMixed)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 40;
	double tol = 1e-2;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);
		problem *pb = problem_create(MIXED, 3, 0);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		double nrm = mpi_norm(grd, sol->state->phi, u);
		norms.push_back(nrm);

		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}


TEST(Models3, Periodic)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 40;
	double tol = 1e-2;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);

		grd->periodic[1] = 1;
		mapping *mp = map_create(MAP_POLAR, grd->nd);
		grid_apply_map(grd, mp);

		problem *pb = problem_create(PERIODIC, 3, mp->id);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		double nrm = mpi_norm(grd, sol->state->phi, u);
		norms.push_back(nrm);

		map_destroy(mp); free(mp);
		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}

TEST(Models3, Jump)
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
		std::vector<double> norms;
		std::vector<double> hs;
		std::array<int, 2> nvals = {40, 80};
		double tol = 1e-1;

		for (int k = 0; k < 2; k++) {
			int i, ierr;

			auto nx = nvals[k];

			hs.push_back(2.0/(nx-1));

			grid *grd = grid_create(-1, 1, nx,
			                        -1, 1, nx,
			                        -1, 1, nx);
			problem *pb = problem_create(JUMP, 3, 0);
			solver *sol = solver_create(grd, pb);

			ierr = solver_init(sol, grd);
			ASSERT_EQ(ierr, 0);
			ierr = solver_run(sol);
			ASSERT_EQ(ierr, 0);

			double *u = new double[grd->num_pts];
			grid_eval(grd, pb->sol, u);

			double nrm = mpi_norm(grd, sol->state->phi, u);
			norms.push_back(nrm);

			problem_destroy(pb); free(pb);
			ierr = solver_destroy(sol); free(sol);
			ASSERT_EQ(ierr, 0);
			grid_destroy(grd); free(grd);
			delete[] u;
		}

		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank == 0) {
			double order = estimate_order(hs, norms);
			ASSERT_GT(order, 2 - tol);
		}
}

TEST(Models, Warping) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 21;
	double tol = 5e-1;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 1, nx);

		mapping *mp = map_create(MAP_DISTORT, grd->nd);
		grid_apply_map(grd, mp);

		problem *pb = problem_create(MTUT, 3, mp->id);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		double nrm = mpi_norm(grd, sol->state->phi, u);
		norms.push_back(nrm);

		map_destroy(mp); free(mp);
		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;

		nx = (nx-1)*2 + 1;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}

TEST(Models3, JumpSine)
{
	std::vector<double> norms;
	std::vector<double> hs;
	std::array<int, 2> nvals = {40, 80};
	double tol = 1e-1;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		auto nx = nvals[k];

		hs.push_back(2.0/(nx-1));

		grid *grd = grid_create(-1, 1, nx,
		                        -1, 1, nx,
		                        -1, 1, nx);
		problem *pb = problem_create(JSINE, 3, 0);
		solver *sol = solver_create(grd, pb);

		ierr = solver_init(sol, grd);
		ASSERT_EQ(ierr, 0);
		ierr = solver_run(sol);
		ASSERT_EQ(ierr, 0);

		double *u = new double[grd->num_pts];
		grid_eval(grd, pb->sol, u);

		double nrm = mpi_norm(grd, sol->state->phi, u);
		norms.push_back(nrm);

		problem_destroy(pb); free(pb);
		ierr = solver_destroy(sol); free(sol);
		ASSERT_EQ(ierr, 0);
		grid_destroy(grd); free(grd);
		delete[] u;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}

}  // namespace
