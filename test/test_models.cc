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
		int i, j;
		for (j = grd->ibeg[1]; j <= grd->iend[1]; j++) {
			for (i = grd->ibeg[0]; i <= grd->iend[0]; i++) {
				int ind = j*grd->len[0] + i;
				diff[ind] = e_sol[ind] - c_sol[ind];
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


TEST(Models, UniformSine)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);
		problem *pb = problem_create(SIN, 2, 0);
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


TEST(Models, UniformTSine)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 301;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);
		problem *pb = problem_create(TSINE, 2, 0);
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


TEST(Models, UniformMixed) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);
		problem *pb = problem_create(MIXED, 2, 0);
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


TEST(Models, UniformMixed1) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);
		problem *pb = problem_create(MIXED_1, 2, 0);
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


TEST(Models, UniformMixed2) {
	PetscErrorCode ierr;
	int i;
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 101;
	double tol = 1e-9;

	grid *grd = grid_create(0, 1, nx,
	                        0, 1, nx,
	                        0, 0, 0);
	problem *pb = problem_create(MIXED_2, 2, 0);
	solver *sol = solver_create(grd, pb);

	ierr = solver_init(sol, grd);
	ASSERT_EQ(ierr, 0);
	ierr = solver_run(sol);
	ASSERT_EQ(ierr, 0);

	double *u = new double[grd->num_pts];
	grid_eval(grd, pb->sol, u);

	double nrm = mpi_norm(grd, sol->state->phi, u);

	ASSERT_LT(nrm, tol);

	problem_destroy(pb); free(pb);
	ierr = solver_destroy(sol); free(sol);
	ASSERT_EQ(ierr, 0);
	grid_destroy(grd); free(grd);
	delete[] u;
}


TEST(Models, StretchSine) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	int ny = 101;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(ny-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, ny,
		                        0, 0, 0);
		problem *pb = problem_create(SIN, 2, 0);
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
		ny = (ny-1)*2 + 1;
	}
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}


TEST(Models, StretchMixed) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	int ny = 101;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(ny-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, ny,
		                        0, 0, 0);

		problem *pb = problem_create(MIXED, 2, 0);
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
		ny = (ny-1)*2 + 1;
	}
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		double order = estimate_order(hs, norms);
		ASSERT_GT(order, 2 - tol);
	}
}


TEST(Models, CurlSine) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);

		mapping *mp = map_create(MAP_CURL, grd->nd);
		grid_apply_map(grd, mp);

		problem *pb = problem_create(SIN, 2, mp->id);
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


TEST(Models, SecondTutorial) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-3;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);

		mapping *mp = map_create(MAP_SECOND, grd->nd);
		grid_apply_map(grd, mp);

		problem *pb = problem_create(MTUT, 2, mp->id);
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


TEST(Models, Rot) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);
		mapping *mp = map_create(MAP_ROT, grd->nd);
		grid_apply_map(grd, mp);

		problem *pb = problem_create(ROT, 2, mp->id);
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


TEST(Models, Axisymmetric) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(.1, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);
		problem *pb = problem_create(AXISYMMETRIC, 2, 0);
		solver *sol = solver_create(grd, pb);
		sol->axisymmetric = 1;

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


TEST(Models, Electrode)
{
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-5;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);
		problem *pb = problem_create(ELECTRODE, 2, 0);
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


TEST(Models, Jump)
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	std::vector<double> norms;
	std::vector<double> hs;
	std::array<int, 2> nvals = {200, 400};
	double tol = 1e-2;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		auto nx = nvals[k];

		hs.push_back(2.0/(nx-1));

		grid *grd = grid_create(-1, 1, nx,
		                        -1, 1, nx,
		                        -1, 1, 0);
		problem *pb = problem_create(JUMP, 2, 0);
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


TEST(Models, HalfAnnulus) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-3;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);

		mapping *mp = map_create(MAP_HANNULUS, grd->nd);
		grid_apply_map(grd, mp);

		problem *pb = problem_create(MTUT, 2, mp->id);
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


TEST(Models, Warping) {
	std::vector<double> norms;
	std::vector<double> hs;
	int nx = 201;
	double tol = 1e-3;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		hs.push_back(1.0/(nx-1));

		grid *grd = grid_create(0, 1, nx,
		                        0, 1, nx,
		                        0, 0, 0);

		mapping *mp = map_create(MAP_DISTORT, grd->nd);
		grid_apply_map(grd, mp);

		problem *pb = problem_create(MTUT, 2, mp->id);
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
TEST(Models, Cboard)
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
		std::vector<double> norms;
		std::vector<double> hs;
		std::array<int, 2> nvals = {200, 400};
		double tol = 1.3e-2;

		for (int k = 0; k < 2; k++) {
			int i, ierr;

			auto nx = nvals[k];

			hs.push_back(2.0/(nx-1));

			grid *grd = grid_create(-1, 1, nx,
			                        -1, 1, nx,
			                        -1, 1, 0);
			problem *pb = problem_create(CBOARD, 2, 0);
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

TEST(Models, JumpSine)
{
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	std::vector<double> norms;
	std::vector<double> hs;
	std::array<int, 2> nvals = {200, 400};
	double tol = 1e-1;

	for (int k = 0; k < 2; k++) {
		int i, ierr;

		auto nx = nvals[k];

		hs.push_back(2.0/(nx-1));

		grid *grd = grid_create(-1, 1, nx,
		                        -1, 1, nx,
		                        -1, 1, 0);
		problem *pb = problem_create(JSINE, 2, 0);
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
