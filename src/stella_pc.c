#include "stella_mat.h"
#include "stella_pc.h"
#include "stella_cedar.h"


#undef __FUNCT__
#define __FUNCT__ "stella_pc_create"
PetscErrorCode stella_pc_create(stella_pc **pc)
{
	stella_pc *ctx;

	ctx = (stella_pc*) malloc(sizeof(stella_pc));
	*pc = ctx;

	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "stella_pc_setup"
PetscErrorCode stella_pc_setup(PC pc)
{
	#ifdef WITH_CEDAR
	PetscErrorCode ierr;
	Mat pmat;
	stella_pc *pc_ctx;
	stella_bmg_mat *mat_ctx;

	ierr = PCGetOperators(pc,PETSC_NULL,&pmat);CHKERRQ(ierr);
	ierr = PCShellGetContext(pc, (void**)&pc_ctx);CHKERRQ(ierr);
	ierr = MatShellGetContext(pmat, (void**)&mat_ctx);CHKERRQ(ierr);

	pc_ctx->nd = mat_ctx->nd;
	pc_ctx->solvec = mat_ctx->solvec;
	pc_ctx->rhsvec = mat_ctx->rhsvec;
    
	int cedar_err;
	cedar_config conf;

	/***********************************/
	// HERE WE NEED TO HAVE THE MESH ID
	size_t meshid = 0;
	/***********************************/

	char fname[19];
	sprintf(fname, "config-mesh%zu.json", meshid);

	cedar_err = cedar_config_create(fname, &conf);
	if (cedar_err == CEDAR_ERR_FNAME) {
		char err_str[80];
		sprintf(err_str, "Cedar config not found: %s - Falling back to config.json.", fname);
		stella_io_print(PETSC_COMM_WORLD, err_str);
		cedar_config_create("config.json", &conf);
	}

	cedar_err = cedar_solver_create(mat_ctx->so, conf, &pc_ctx->solver);chkerr(cedar_err);

	if (cedar_err) {
		SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_LIB, "Cedar error code: %d", cedar_err);
	}
	#endif

	return 0;
}


#undef __FUNCT__
#define __FUNCT__ "stella_pc_apply"
PetscErrorCode stella_pc_apply(PC pc, Vec x, Vec y)
{
	#ifdef WITH_CEDAR
	int cedar_err;
	PetscErrorCode ierr;
	stella_pc *pc_ctx;
	double *yarr;
	const double *xarr;

	ierr = PCShellGetContext(pc, (void**)&pc_ctx);CHKERRQ(ierr);
	ierr = VecGetArrayRead(x, &xarr);CHKERRQ(ierr);
	ierr = VecGetArray(y, &yarr);CHKERRQ(ierr);

	cedar_copyto(xarr, pc_ctx->rhsvec);
	cedar_err = cedar_solver_run(pc_ctx->solver, pc_ctx->solvec, pc_ctx->rhsvec);chkerr(cedar_err);
	cedar_copyfrom(pc_ctx->solvec, yarr);

	ierr = VecRestoreArrayRead(x, &xarr);CHKERRQ(ierr);
	ierr = VecRestoreArray(y, &yarr);CHKERRQ(ierr);
	#endif

	return 0;
}



#undef __FUNCT__
#define __FUNCT__ "stella_pc_destroy"
PetscErrorCode stella_pc_destroy(PC pc)
{
	PetscErrorCode ierr;
	stella_pc *pc_ctx;

	ierr = PCShellGetContext(pc, (void**)&pc_ctx);CHKERRQ(ierr);
	free(pc_ctx);

	return 0;
}
