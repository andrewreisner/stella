#ifndef STELLA_STATE_H
#define STELLA_STATE_H


/**
 * Contains PlasComCM points to state variables
 * see stella_interface.h for descriptions of members
 */
typedef struct {
	double *phi;
	double *jump[3];
	double *dcoef;
	double *bcoef;
	double *rhs;
	double *sol;
} stella_state;


#endif
