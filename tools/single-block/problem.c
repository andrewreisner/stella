#include <stdlib.h>
#include <math.h>

#include "problem.h"
#include "map.h"
#include "boundary.h"


const char *problem_name[NUM_PROBLEMS] = {"Multigrid Tutorial", "Mixed Boundaries",
                                          "Mixed Boundaries 1", "Mixed Boundaries 2",
                                          "Sine", "Translated Sine",
                                          "Rotation", "Electrode", "Jump",
                                          "Axisymmetric", "Periodic", "Checkerboard", "Jump Sine"};
const char *problem_key[NUM_PROBLEMS] = {"tutorial", "mixed", "mixed1", "mixed2", "sin",
                                         "tsine", "rotation",
                                         "electrode", "jump","axisymmetric", "periodic", "cboard", "jsine"};

double electrode_rhs(double x, double y, double z)
{
	return 2*(M_PI*M_PI)*cos(M_PI*x)*cos(M_PI*y);
}


double electrode_sol(double x, double y, double z)
{
	return cos(M_PI*x)*cos(M_PI*y);
}


double electrode_rhs_3d(double x, double y, double z)
{
	return 3*(M_PI*M_PI)*cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*z);
}


double electrode_sol_3d(double x, double y, double z)
{
	return cos(M_PI*x)*cos(M_PI*y)*cos(M_PI*z);
}


double axisymmetric_rhs(double x, double y, double z)
{
	return 1./x*(8*pow(M_PI, 2)*x*sin(2*M_PI*x)*sin(2*M_PI*y) - 2*M_PI*sin(2*M_PI*y)*cos(2*M_PI*x));
}


double axisymmetric_sol(double x, double y, double z)
{
	return sin(2*M_PI*x)*sin(2*M_PI*y);
}


double mtut_rhs(double x, double y, double z)
{
	return 2*((1-(6*(x*x)))*(y*y)*(1-(y*y)) + (1 - (6*(y*y)))*(x*x)*(1-(x*x)));
}


double mtut_sol(double x, double y, double z)
{
	return (x*x - x*x*x*x)*(y*y*y*y - y*y);
}

double mtut_rhs_3d(double x, double y, double z)
{
    return 12*(M_PI*M_PI)*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
}

double mtut_sol_3d(double x, double y, double z)
{
    return sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
}

double sin_rhs(double x, double y, double z)
{
	return 8*(M_PI*M_PI)*sin(2*M_PI*x)*sin(2*M_PI*y);
}


double sin_sol(double x, double y, double z)
{
	return sin(2*M_PI*x)*sin(2*M_PI*y);
}


double tsin_rhs(double x, double y, double z)
{
	return 8*(M_PI*M_PI)*sin(2*M_PI*x + 1.4)*sin(2*M_PI*y+1);
}


double tsin_sol(double x, double y, double z)
{
	return sin(2*M_PI*x+1.4)*sin(2*M_PI*y+1);
}


double polar_tsin_rhs(double x, double y, double z)
{
	return M_PI*M_PI*x*x*sin(M_PI*(sqrt(x*x + y*y) - 0.5))*sin(atan2(y, x) + 1.4)/(x*x + y*y) + x*x*sin(M_PI*(sqrt(x*x + y*y) - 0.5))*sin(atan2(y, x) + 1.4)/pow((x*x + y*y),2) + M_PI*x*x*sin(atan2(y, x) + 1.4)*cos(M_PI*(sqrt(x*x + y*y) - 0.5))/pow((x*x + y*y),(3./2)) + M_PI*M_PI*y*y*sin(M_PI*(sqrt(x*x + y*y) - 0.5))*sin(atan2(y, x) + 1.4)/(x*x + y*y) + y*y*sin(M_PI*(sqrt(x*x + y*y) - 0.5))*sin(atan2(y, x) + 1.4)/pow((x*x + y*y),2) + M_PI*y*y*sin(atan2(y, x) + 1.4)*cos(M_PI*(sqrt(x*x + y*y) - 0.5))/pow((x*x + y*y),(3./2)) - 2*M_PI*sin(atan2(y, x) + 1.4)*cos(M_PI*(sqrt(x*x + y*y) - 0.5))/sqrt(x*x + y*y);
}


double polar_tsin_sol(double x, double y, double z)
{
	double tmpx = sqrt(x*x + y*y) - .5;
	y = atan2(y,x) / (2*M_PI);
	x = tmpx;

	return sin(M_PI*x)*sin(2*M_PI*y + 1.4);
}


double polar_sin_rhs(double x, double y, double z)
{
	return -y*(-6*M_PI*x*x*cos(M_PI*(2*sqrt(x*x + y*y) - 1.0))/pow((x*x + y*y),2) - 4*(M_PI*M_PI)*x*x*sin(M_PI*(2*sqrt(x*x + y*y) - 1.0))/pow((x*x + y*y),(3./2)) + 3*x*x*sin(M_PI*(2*sqrt(x*x + y*y) - 1.0))/pow((x*x + y*y),(5./2)) + 2*M_PI*cos(M_PI*(2*sqrt(x*x + y*y) - 1.0))/(x*x + y*y) - sin(M_PI*(2*sqrt(x*x + y*y) - 1.0))/pow((x*x + y*y),(3./2))) - y*(-6*M_PI*y*y*cos(M_PI*(2*sqrt(x*x + y*y) - 1.0))/pow((x*x + y*y),2) - 4*(M_PI*M_PI)*y*y*sin(M_PI*(2*sqrt(x*x + y*y) - 1.0))/pow((x*x + y*y),(3./2)) + 3*y*y*sin(M_PI*(2*sqrt(x*x + y*y) - 1.0))/pow((x*x + y*y),(5./2)) + 6*M_PI*cos(M_PI*(2*sqrt(x*x + y*y) - 1.0))/(x*x + y*y) - 3*sin(M_PI*(2*sqrt(x*x + y*y) - 1.0))/pow((x*x + y*y),(3./2)));
}


double polar_sin_sol(double x, double y, double z)
{
	double tmpx = sqrt(x*x + y*y) - .5;
	y = atan2(y,x) / (2*M_PI);
	x = tmpx;

	return sin(2*M_PI*x)*sin(2*M_PI*y);
}


double polar_sin_rhs_3d(double r, double s, double t)
{
	return -s*(-6*M_PI*pow(r, 2)*cos(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/pow(pow(r, 2) + pow(s, 2), 2) - 4*pow(M_PI, 2)*pow(r, 2)*sin(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/pow(pow(r, 2) + pow(s, 2), 3.0L/2.0L) + 3*pow(r, 2)*sin(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/pow(pow(r, 2) + pow(s, 2), 5.0L/2.0L) + 2*M_PI*cos(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/(pow(r, 2) + pow(s, 2)) - sin(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/pow(pow(r, 2) + pow(s, 2), 3.0L/2.0L))*sin(2*M_PI*t) - s*(-6*M_PI*pow(s, 2)*cos(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/pow(pow(r, 2) + pow(s, 2), 2) - 4*pow(M_PI, 2)*pow(s, 2)*sin(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/pow(pow(r, 2) + pow(s, 2), 3.0L/2.0L) + 3*pow(s, 2)*sin(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/pow(pow(r, 2) + pow(s, 2), 5.0L/2.0L) + 6*M_PI*cos(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/(pow(r, 2) + pow(s, 2)) - 3*sin(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/pow(pow(r, 2) + pow(s, 2), 3.0L/2.0L))*sin(2*M_PI*t) + 4*pow(M_PI, 2)*s*sin(2*M_PI*t)*sin(M_PI*(2*sqrt(pow(r, 2) + pow(s, 2)) - 1.0))/sqrt(pow(r, 2) + pow(s, 2));
}


double polar_sin_sol_3d(double x, double y, double z)
{
	double tmpx = sqrt(x*x + y*y) - .5;
	y = atan2(y,x) / (2*M_PI);
	x = tmpx;

	return sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
}


double sin_rhs_3d(double x, double y, double z)
{
	return 12*(M_PI*M_PI)*sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
}


double sin_sol_3d(double x, double y, double z)
{
	return sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
}


double tsin_rhs_3d(double x, double y, double z)
{
	return 12*(M_PI*M_PI)*sin(2*M_PI*x + 1.4)*sin(2*M_PI*y+1)*sin(2*M_PI*z+1.6);
}


double tsin_sol_3d(double x, double y, double z)
{
	return sin(2*M_PI*x+1.4)*sin(2*M_PI*y+1)*sin(2*M_PI*z+1.6);
}


double mixed_rhs(double x, double y, double z)
{
	//return 4*(M_PI*M_PI)*(cos(2*M_PI*x) - 1)*cos(2*M_PI*y) + 4*(M_PI*M_PI)*(cos(2*M_PI*y) - 1)*cos(2*M_PI*x);
	return 1.25*(M_PI*M_PI)*sin(0.5*M_PI*y)*cos(M_PI*x);
}


double mixed_sol(double x, double y, double z)
{
	//return (1- cos(2*M_PI*x))*(1- cos(2*M_PI*y));
	return cos(M_PI*x)*sin(.5*M_PI*y);
}


double mixed1_rhs(double x, double y, double z)
{
	return 1.25*(M_PI*M_PI)*sin(.5*M_PI*x)*cos(M_PI*y);
}


double mixed1_sol(double x, double y, double z)
{
	return cos(M_PI*y)*sin(.5*M_PI*x);
}


double mixed2_rhs(double x, double y, double z)
{
	return -2;
}


double mixed2_sol(double x, double y, double z)
{
	return x*x;
}

double mixed_rhs_3d(double x, double y, double z)
{
	return -4*(M_PI*M_PI)*(cos(2*M_PI*x) - 1)*(cos(2*M_PI*y) - 1)*cos(2*M_PI*z) - 4*(M_PI*M_PI)*(cos(2*M_PI*x) - 1)*(cos(2*M_PI*z) - 1)*cos(2*M_PI*y) - 4*(M_PI*M_PI)*(cos(2*M_PI*y) - 1)*(cos(2*M_PI*z) - 1)*cos(2*M_PI*x);
}


double mixed_sol_3d(double x, double y, double z)
{
	return (1- cos(2*M_PI*x))*(1- cos(2*M_PI*y)) * (1-cos(2*M_PI*z));
}


double rot_sol(double x, double y, double z)
{
	return sin(M_PI*(x + y + 1))*sin(M_PI*(y - x + 1));
}


double rot_rhs(double x, double y, double z)
{
	return -2*(M_PI*M_PI)*(-sin((M_PI)*(-x + y + 1))*sin(M_PI*(x + y + 1)) + cos(M_PI*(-x + y + 1))*cos(M_PI*(x + y + 1))) + 2*(M_PI*M_PI)*(sin(M_PI*(-x + y + 1))*sin(M_PI*(x + y + 1)) + cos(M_PI*(-x + y + 1))*cos(M_PI*(x + y + 1)));

}


double jump_sol(double x, double y, double z)
{
	if (x <= 0)
		return (2*(x*x) + (13./4)*x);
	else
		return (-4*(x*x) + 4*x);
}


double jump_rhs(double x, double y, double z)
{
	if (x <= 0)
		return -16;
	else
		return 16;
}

double cboard_sol(double x, double y, double z){
    if (x <= 0 && y < 0)
		return 2*(x*x) + 13.0/4*x - 4*(y*y) + 4*y;
    else if (x > 0 && y >= 0)
        return -4*(x*x) + 4*x + 2*(y*y) + 13.0/4*y;
    else if (x <= 0 && y >= 0)
        return 2*(x*x) + 13.0/4*x + 2*(y*y) + 13.0/4*y;
    else
        return -4*(x*x) + 4*x - 4*(y*y) + 4*y;
}

double cboard_rhs(double x, double y, double z){
    if (x <= 0 && y < 0)
        return 8;
    else if (x > 0 && y >= 0)
        return 8;
    else if (x <= 0 && y >=0)
        return -32;
    else
        return 64;
}

double jump_sin_rhs(double x, double y, double z)
{
    if (x <= 0)
	   return 25*(M_PI*M_PI)*sin(M_PI*x)*sin(2*M_PI*y);
    else
        return 16*(M_PI*M_PI)*sin(2*M_PI*x)*sin(2*M_PI*y);
}


double jump_sin_sol(double x, double y, double z)
{
    if (x <= 0)
	   return sin(M_PI*x)*sin(2*M_PI*y);
    else
        return sin(2*M_PI*x)*sin(2*M_PI*y);
}

problem *problem_create(problem_id id, int nd, int map_id)
{
	problem *pb = (problem*) malloc(sizeof(problem));
	pb->id = id;
	pb->nd = nd;

	pb->nholes = 0;

	switch(id) {
	case (ELECTRODE):
		pb->nholes += 2;
		if (nd == 3)
			pb->nholes += 2;
		pb->holes = (electrode*) malloc(pb->nholes*sizeof(electrode));

		pb->holes[0].rel_size[0] = 1.0 / 16.0;
		pb->holes[0].rel_size[1] = 1.0 / 16.0;
		pb->holes[0].rel_offset[0] = 0.25;
		pb->holes[0].rel_offset[1] = 0.25;

		pb->holes[1].rel_size[0] = 1.0 / 16.0;
		pb->holes[1].rel_size[1] = 1.0 / 16.0;
		pb->holes[1].rel_offset[0] = .75;
		pb->holes[1].rel_offset[1] = .75;

		if (nd == 2) {
			pb->rhs = &electrode_rhs;
			pb->sol = &electrode_sol;
		} else {
			pb->rhs = &electrode_rhs_3d;
			pb->sol = &electrode_sol_3d;
			pb->holes[0].rel_size[2] = 1.0 / 16.0;
			pb->holes[0].rel_offset[2] = 0.25;

			pb->holes[1].rel_size[2] = 1.0 / 16.0;
			pb->holes[1].rel_offset[2] = .25;

		pb->holes[2].rel_size[0] = 1.0 / 16.0;
		pb->holes[2].rel_size[1] = 1.0 / 16.0;
		pb->holes[2].rel_offset[0] = 0.25;
		pb->holes[2].rel_offset[1] = 0.25;
		pb->holes[2].rel_size[2] = 1.0/16.0;
		pb->holes[2].rel_offset[2] = 0.75;

		pb->holes[3].rel_size[0] = 1.0 / 16.0;
		pb->holes[3].rel_size[1] = 1.0 / 16.0;
		pb->holes[3].rel_offset[0] = .75;
		pb->holes[3].rel_offset[1] = .75;
		pb->holes[3].rel_size[2] = 1.0 / 16.0;
		pb->holes[3].rel_offset[2] = 0.75;

			pb->boundary[FRONT] = NEUMANN;
			pb->boundary[BACK] = NEUMANN;
		}

		pb->boundary[NORTH] = NEUMANN;
		pb->boundary[SOUTH] = NEUMANN;
		pb->boundary[EAST] =  NEUMANN;
		pb->boundary[WEST] =  NEUMANN;
		break;
	case (MTUT):
        pb->boundary[NORTH] = DIRICHLET;
        pb->boundary[SOUTH] = DIRICHLET;
        pb->boundary[EAST] =  DIRICHLET;
        pb->boundary[WEST] =  DIRICHLET;
        if (nd == 3) {
            pb->rhs = &mtut_rhs_3d;
    		pb->sol = &mtut_sol_3d;
        } else {
            pb->rhs = &mtut_rhs;
    		pb->sol = &mtut_sol;
        }
		break;
	case (SIN):
		pb->boundary[NORTH] = DIRICHLET;
		pb->boundary[SOUTH] = DIRICHLET;
		pb->boundary[EAST] =  DIRICHLET;
		pb->boundary[WEST] =  DIRICHLET;
		if (nd == 3) {
			if (map_id == MAP_POLAR) {
				pb->rhs = &polar_sin_rhs_3d;
				pb->sol = &polar_sin_sol_3d;
			} else {
				pb->rhs = &sin_rhs_3d;
				pb->sol = &sin_sol_3d;
			}
			pb->boundary[FRONT] =  DIRICHLET;
			pb->boundary[BACK] =  DIRICHLET;
		} else {
			if (map_id == MAP_POLAR) {
				pb->rhs = &polar_sin_rhs;
				pb->sol = &polar_sin_sol;
			} else {
				pb->rhs = &sin_rhs;
				pb->sol = &sin_sol;
			}
		}
		break;
	case(TSINE):
		pb->boundary[NORTH] = DIRICHLET;
		pb->boundary[SOUTH] = DIRICHLET;
		pb->boundary[EAST] =  DIRICHLET;
		pb->boundary[WEST] =  DIRICHLET;
		if (map_id == MAP_POLAR) {
			pb->sol = &polar_tsin_sol;
			pb->rhs = &polar_tsin_rhs;
		} else {
			if (nd == 2) {
				pb->sol = &tsin_sol;
				pb->rhs = &tsin_rhs;
			} else {
				pb->sol = &tsin_sol_3d;
				pb->rhs = &tsin_rhs_3d;
				pb->boundary[FRONT] =  DIRICHLET;
				pb->boundary[BACK] =  DIRICHLET;
			}
		}
		break;
	case (MIXED):
		if (nd == 2) {
			pb->rhs = &mixed_rhs;
			pb->sol = &mixed_sol;
			pb->boundary[NORTH] = NEUMANN;
			pb->boundary[SOUTH] = DIRICHLET;
			pb->boundary[EAST] =  NEUMANN;
			pb->boundary[WEST] =  NEUMANN;
		} else {
			pb->rhs = &mixed_rhs_3d;
			pb->sol = &mixed_sol_3d;
			pb->boundary[NORTH] = DIRICHLET;
			pb->boundary[SOUTH] = NEUMANN;
			pb->boundary[EAST] =  NEUMANN;
			pb->boundary[WEST] =  DIRICHLET;
			pb->boundary[FRONT] = NEUMANN;
			pb->boundary[BACK] = NEUMANN;
		}

		break;
	case(MIXED_1):
		pb->rhs = &mixed1_rhs;
		pb->sol = &mixed1_sol;
		pb->boundary[NORTH] = NEUMANN;
		pb->boundary[SOUTH] = NEUMANN;
		pb->boundary[EAST] =  NEUMANN;
		pb->boundary[WEST] =  DIRICHLET;
		break;
	case(MIXED_2):
		pb->rhs = &mixed2_rhs;
		pb->sol = &mixed2_sol;
		pb->boundary[NORTH] = NEUMANN;
		pb->boundary[SOUTH] = NEUMANN;
		pb->boundary[EAST] =  DIRICHLET;
		pb->boundary[WEST] =  DIRICHLET;
		break;
	case (ROT):
		pb->rhs = &rot_rhs;
		pb->sol = &rot_sol;
		pb->boundary[NORTH] = DIRICHLET;
		pb->boundary[SOUTH] = DIRICHLET;
		pb->boundary[EAST] =  DIRICHLET;
		pb->boundary[WEST] =  DIRICHLET;
		break;
	case (JUMP):
		pb->rhs = &jump_rhs;
		pb->sol = &jump_sol;
		pb->boundary[NORTH] = DIRICHLET;
		pb->boundary[SOUTH] = DIRICHLET;
		pb->boundary[EAST] =  DIRICHLET;
		pb->boundary[WEST] =  DIRICHLET;
		if (nd == 3) {
			pb->boundary[FRONT] =  DIRICHLET;
			pb->boundary[BACK] =  DIRICHLET;
		}
		break;
	case (AXISYMMETRIC):
		pb->rhs = &axisymmetric_rhs;
		pb->sol = &axisymmetric_sol;
		pb->boundary[NORTH] = DIRICHLET;
		pb->boundary[SOUTH] = DIRICHLET;
		pb->boundary[EAST] =  DIRICHLET;
		pb->boundary[WEST] =  DIRICHLET;
		break;
	case (PERIODIC):
		pb->rhs = &polar_tsin_rhs;
		pb->sol = &polar_tsin_sol;
		pb->boundary[EAST] =  DIRICHLET;
		pb->boundary[WEST] =  DIRICHLET;
		pb->boundary[FRONT] = DIRICHLET;
		pb->boundary[BACK] = DIRICHLET;
		break;
	default:
		pb->rhs = &mtut_rhs;
		pb->sol = &mtut_sol;
		break;
    case (CBOARD):
    	pb->rhs = &cboard_rhs;
    	pb->sol = &cboard_sol;
    	pb->boundary[NORTH] = DIRICHLET;
    	pb->boundary[SOUTH] = DIRICHLET;
    	pb->boundary[EAST] =  DIRICHLET;
    	pb->boundary[WEST] =  DIRICHLET;
    	if (nd == 3) {
    		pb->boundary[FRONT] =  DIRICHLET;
    		pb->boundary[BACK] =  DIRICHLET;
    	}
    	break;
    case (JSINE):
    	pb->rhs = &jump_sin_rhs;
    	pb->sol = &jump_sin_sol;
    	pb->boundary[NORTH] = DIRICHLET;
    	pb->boundary[SOUTH] = DIRICHLET;
    	pb->boundary[EAST] =  DIRICHLET;
    	pb->boundary[WEST] =  DIRICHLET;
    	if (nd == 3) {
    		pb->boundary[FRONT] =  DIRICHLET;
    		pb->boundary[BACK] =  DIRICHLET;
    	}
    	break;
	}

	return pb;
}


void problem_destroy(problem *pb)
{
	if (pb->nholes > 0)
		free(pb->holes);
}
