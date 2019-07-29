#include <ceed.h>
#include "IsotropicLinearElasticityCEED.h"




void pylith::fekernels::IsotropicLinearElasticityCEED::CEED_integrate(){

	Ceed ceed;
	CeedQfunction qf_setup, qf_mass, qf_physics;

	///setup problem

	///Setup setup Qfunction 
	CeedQfunctionCreateInterior(ceed, 1, Setup, __FILE__":Setup", &qf_setup);
	CeedQFunctionAddInput(qf_setup, "J", dim, CEED_EVAL_GRAD);
	CeedQFunctionAddInput(qf_setup, "w", 1, CEED_EVAL_WEIGHT);
	CeedQFunctionAddOutput(qf_setup, "qdata", dim*dim, CEED_EVAL_NONE);

	///Setup mass Qfunction
	CeedQfunctionCreateInterior(ceed, 1, ApplyMass, __FILE__":ApplyMass", &qf_mass);
	CeedQFunctionAddInput(qf_mass,"u", dim, CEED_EVAL_INTERP);
	CeedQFunctionAddInput(qf_mass,"w", dim*dim, CEED_EVAL_NONE);
	CeedQFunctionAddOutput(qf_mass,"v", dim, CEED_EVAL_INTERP );

	///Setup physics Qfunction
	CeedQfunctionCreateInterior(ceed, 1, IsotropicLinearElasticity, __FILE__":IsotropicLinearElasticity", &qf_physics);
	CeedQFunctionAddInput(qf_mass, "disp", dim, CEED_EVAL_INTERP);
	CeedQFunctionAddInput(qf_mass, "graddisp", dim*dim, CEED_EVAL_NONE);
	CeedQFunctionAddOutput(qf_mass, "v", dim, CEED_EVAL_INTERP);
	CeedQFunctionAddOutput(qf_mass, "dv", dim*dim, CEED_EVAL_GRAD);

	///operator for setup

	///operator for mass

	///operator for physics

	return 0;
}




