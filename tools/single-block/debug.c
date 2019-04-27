#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <mpi.h>

#include "problem.h"
#include "map.h"
#include "grid.h"
#include "solver.h"
#include "option.h"



double debug_map_x(double r, double s, double t)
{
	double r1 = .5;
	double r2 = 1;
	return (r1 + (r2 - r1)*s)*sin(M_PI*r);
}



double debug_map_y(double r, double s, double t)
{
	double r1 = .5;
	double r2 = 1;
	return (r1 + (r2 - r1)*s)*cos(M_PI*r);
}


static double dcoef(double x, double y, double z)
{
	//return 10;
	return 10 - x;
	/* return 10 + (double)rand() / (double)RAND_MAX; */
}


static double nojump(double x, double y, double z)
{
	return 0;
}


static double test_rhs(double x, double y, double z)
{
	return 1.0;
}


static double test_sol(double x, double y, double z)
{
	return 0.0;
}


int main(int argc, char *argv[])
{
	PetscErrorCode ierr;

	MPI_Init(&argc,&argv);

	grid *grd = grid_create(0, 1, 50,
	                        0, 1, 20,
	                        0, 1, 0);

	/* mapping *mp = map_create(MAP_HANNULUS, grd->nd); */
	mapping *mp = (mapping*) malloc(sizeof(mapping));
	mp->id = MAP_HANNULUS;
	mp->x = &debug_map_x;
	mp->y = &debug_map_y;
	grid_apply_map(grd, mp);

	problem *pb = (problem*) malloc(sizeof(problem));
	{ // setup problem
		pb->id = SIN;
		pb->nd = grd->nd;
		pb->nholes = 0;

		pb->boundary[NORTH] = DIRICHLET;
		pb->boundary[SOUTH] = DIRICHLET;
		pb->boundary[EAST] =  NEUMANN;
		pb->boundary[WEST] =  NEUMANN;
		pb->eps = &dcoef;
		for (int i = 0; i < 3; i++) pb->jc[i] = &nojump;
		pb->rhs = &test_rhs;
		pb->sol = &test_sol;
	}

	solver *sol = solver_create(grd, pb);

	ierr = solver_init(sol, grd);CHKERRQ(ierr);
	ierr = solver_run(sol);CHKERRQ(ierr);

	grid_destroy(grd); free(grd);
	map_destroy(mp); free(mp);
	free(pb);
	ierr = solver_destroy(sol);CHKERRQ(ierr); free(sol);

	MPI_Finalize();

	return 0;
}
