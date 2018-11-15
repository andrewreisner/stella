#ifndef PROBLEM_H
#define PROBLEM_H

#include "base.h"

#ifdef __cplusplus
extern "C" {
#endif

#define NUM_PROBLEMS 13

typedef enum {
	MTUT=0, MIXED=1, MIXED_1=2, MIXED_2=3,
	SIN=4, TSINE=5, ROT=6, ELECTRODE=7,
	JUMP=8, AXISYMMETRIC=9, PERIODIC=10, CBOARD=11, JSINE=12
} problem_id;

extern const char *problem_name[NUM_PROBLEMS];
extern const char *problem_key[NUM_PROBLEMS];

typedef struct {
	double rel_size[3];
	double rel_offset[3];
} electrode;

typedef struct {
	problem_id id;
	func rhs;
	func sol;
	int nd;
	int boundary[6];
	electrode *holes;
	int nholes;
} problem;

problem *problem_create(problem_id, int nd, int map_id);

void problem_destroy(problem*);


#ifdef __cplusplus
}
#endif

#endif
