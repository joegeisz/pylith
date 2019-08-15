// -*- C++ -*-
//
// ---------------------------------------------------------------------------------------------------------------------
//
//Joe Geisz
// 
// IntegratorDomain with CEED implementation
//
// ---------------------------------------------------------------------------------------------------------------------
//

#include <portinfo>

#include "IntegratorDomainCEED.hh" // implementation of object methods

#include "pylith/feassemble/UpdateStateVars.hh" // HOLDSA UpdateStateVars
#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/MeshOps.hh" // USES createSubdomainMesh()
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/CoordsVisitor.hh" // USES CoordsVisitor

#include "pylith/fekernels/IsotropicLinearElasticityCEED.hh"//uses q-functions

#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField
#include "petscds.h" // USES PetscDS


#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error
#include <ceed.h>

// Utility function to create local CEED restriction
static int CreateRestrictionFromPlex(Ceed ceed, DM dm, CeedInt P,
                                     CeedElemRestriction *Erestrict, CeedInt dim) {

  PetscSection   section;
  PetscInt       c, cStart, cEnd, Nelem, Ndof, *erestrict, eoffset, nfields;
  PetscInt       ierr, pp;
  Vec Uloc;

  if(dim == 2){pp=P*P;}
  else if (dim == 3){pp=P*P*P;}
  else{printf("dimension error in CreateRestrictionFromPlex\n");}

  PetscFunctionBeginUser;
  ierr = DMGetDefaultSection(dm,&section);CHKERRQ(ierr);
  ierr = PetscSectionGetNumFields(section, &nfields);CHKERRQ(ierr);
  PetscInt ncomp[nfields], fieldoff[nfields+1];
  fieldoff[0] = 0;
  for (PetscInt f=0; f<nfields; f++) {
    ierr = PetscSectionGetFieldComponents(section, f, &ncomp[f]);CHKERRQ(ierr);
    fieldoff[f+1] = fieldoff[f] + ncomp[f];
  }

  ierr = DMPlexGetHeightStratum(dm,0,&cStart,&cEnd);CHKERRQ(ierr);
  Nelem = cEnd - cStart;
  ierr = PetscMalloc1(Nelem*pp, &erestrict);CHKERRQ(ierr);
  for (c=cStart,eoffset=0; c<cEnd; c++) {
    PetscInt numindices, *indices, nnodes;
    ierr = DMPlexGetClosureIndices(dm,section,section,c,&numindices,&indices,NULL);CHKERRQ(ierr);
    if (numindices % fieldoff[nfields]) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Number of closure indices not compatible with Cell %D",c);
    nnodes = numindices / fieldoff[nfields];
    for (PetscInt i=0; i<nnodes; i++) {
      // Check that indices are blocked by node and thus can be coalesced as a single field with
      // fieldoff[nfields] = sum(ncomp) components.
      for (PetscInt f=0; f<nfields; f++) {
        for (PetscInt j=0; j<ncomp[f]; j++) {
          if (indices[fieldoff[f]*nnodes + i*ncomp[f] + j]
              != indices[i*ncomp[0]] + fieldoff[f] + j)
            SETERRQ4(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Cell %D closure indices not interlaced for node %D field %D component %D",c,i,f,j);
        }
      }
      // Essential boundary conditions are encoded as -(loc+1), but we don't care so we decode.
      PetscInt loc = indices[i*ncomp[0]];
      if (loc < 0) loc = -(loc+1);
      erestrict[eoffset++] = loc / fieldoff[nfields];
    }
    ierr = DMPlexRestoreClosureIndices(dm,section,section,c,&numindices,&indices,NULL);CHKERRQ(ierr);
  }
  if (eoffset != Nelem*pp) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_LIB,"ElemRestriction of size (%D,%D) initialized %D nodes",Nelem,pp,eoffset);
  ierr = DMGetLocalVector(dm, &Uloc);CHKERRQ(ierr);
  ierr = VecGetLocalSize(Uloc, &Ndof);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm, &Uloc);CHKERRQ(ierr);
  CeedElemRestrictionCreate(ceed, Nelem, pp, Ndof/fieldoff[nfields], fieldoff[nfields],
                            CEED_MEM_HOST, CEED_COPY_VALUES, erestrict, Erestrict);
  ierr = PetscFree(erestrict);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static int CreateVectorFromPetscVec(Ceed ceed, Vec p, CeedVector *v) {
  PetscErrorCode ierr;
  PetscInt m;

  PetscFunctionBeginUser;
  ierr = VecGetLocalSize(p, &m);CHKERRQ(ierr);
  ierr = CeedVectorCreate(ceed, m, v);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::feassemble::IntegratorDomainCEED::IntegratorDomainCEED(pylith::problems::Physics* const physics) :
   	IntegratorDomain(physics)
    {
    printf("constructing CEED Integrator Domain\n");
}; // constructor

// ---------------------------------------------------------------------------------------------------------------------
// Default destructor.
pylith::feassemble::IntegratorDomainCEED::~IntegratorDomainCEED(void){

   // this->IntegratorDomain::~IntegratorDomain();
    
    ///Clean up libCEED
	CeedOperatorDestroy(&op);

}; // constructor





// ---------------------------------------------------------------------------------------------------------------------
// Setup everything needed to apply operator
void
pylith::feassemble::IntegratorDomainCEED::initialize(const pylith::topology::Field& solution) {

    printf("initialize IntegratorDomainCEED\n");



    //hard-coded- todo: find how to get these
	char ceedresource[4096] = "/cpu/self";
	PetscInt degree = 1;
	PetscInt qextra = 2;
	PetscInt dim = 2;

	
	Ceed ceed;
	CeedQFunction qf_setup, qf_mass, qf_physics;
	CeedElemRestriction restrictq, restrictx,restrictqdi,restrictxi;
	CeedBasis basis;
	CeedOperator op_setup, op_mass;
	CeedVector qdata, localcoordsceed, q0ceed, mceed, xceed, onesvec;
	PetscVec localcoords;
	PetscInt numP = degree + 1;
	PetscInt numQ = numP + qextra;
  PetscInt Qs = numQ;
  PetscInt jacobiansize = 0;

  if(dim == 2){
    Qs=numQ*numQ;
    jacobiansize = 4;
  }
  else if(dim == 3){
    Qs=numQ*numQ*numQ;
    jacobiansize = 9;
  }
  else{
    printf("dimension is wrong\n");
  }


	///setup CEED
	CeedInit(ceedresource, &ceed);
	

	///setup Bases
	CeedBasisCreateTensorH1Lagrange(ceed, dim, 3, numP, numQ, CEED_GAUSS, &basis);
	


	//setup subdomain mesh
	/*
	delete _materialMesh;
    _materialMesh = pylith::topology::MeshOps::createSubdomainMesh(solution.mesh(), "material-id", _materialId, ":UNKOWN:");
    pylith::topology::CoordsVisitor::optimizeClosure(_materialMesh->dmMesh());
	*/
  printf("before init\n");
  this->IntegratorDomain::initialize(solution);
  printf("after init\n");



	PetscDM dmSoln = solution.dmMesh();
  //or =  _materialMesh->dmMesh();???


	DMGetCoordinatesLocal(dmSoln, &localcoords);
	CreateVectorFromPetscVec(ceed,localcoords,&localcoordsceed);


	///setup Restrictions/setup vector for qdata

	CeedInt Nqpts;
	PetscInt localNelem, End, Start;
	CreateRestrictionFromPlex(ceed, dmSoln, degree + 1, &restrictq, dim);
	CreateRestrictionFromPlex(ceed, dmSoln, degree + 1, &restrictx, dim);
	DMPlexGetHeightStratum(dmSoln, 0, &Start, &End);
	localNelem = End - Start;
  CeedBasisGetNumQuadraturePoints(basis, &Nqpts);


	CeedVectorCreate(ceed, jacobiansize*localNelem*Nqpts ,&qdata);
	CeedElemRestrictionCreateVector(restrictq, &q0ceed, NULL);
	CeedElemRestrictionCreateVector(restrictq, &mceed, NULL);
	CeedElemRestrictionCreateVector(restrictq, &onesvec, NULL);
	CeedVectorSetValue(onesvec,1.0);
  CeedElemRestrictionCreateVector(restrictx, &xceed, NULL);
  CeedElemRestrictionCreateIdentity(ceed, localNelem, jacobiansize*Qs,
                                    jacobiansize*localNelem*Qs, 1,
                                    &restrictqdi);
  CeedElemRestrictionCreateIdentity(ceed, localNelem, Qs,
                                    localNelem*Qs, 1,
                                    &restrictxi);



	///Setup setup Qfunction 
	CeedQFunctionCreateInterior(ceed, 1, Setup, __FILE__":Setup", &qf_setup);
	CeedQFunctionAddInput(qf_setup, "J", dim, CEED_EVAL_GRAD);
	CeedQFunctionAddInput(qf_setup, "w", 1, CEED_EVAL_WEIGHT);
	CeedQFunctionAddOutput(qf_setup, "qdata", dim*dim, CEED_EVAL_NONE);

	///Setup mass Qfunction
	CeedQFunctionCreateInterior(ceed, 1, ApplyMass, __FILE__":ApplyMass", &qf_mass);
	CeedQFunctionAddInput(qf_mass,"u", dim, CEED_EVAL_INTERP);
	CeedQFunctionAddInput(qf_mass,"qdata", dim*dim, CEED_EVAL_NONE);
	CeedQFunctionAddOutput(qf_mass,"v", dim, CEED_EVAL_INTERP );

	///Setup physics Qfunction
	CeedQFunctionCreateInterior(ceed, 1, IsotropicLinearElasticityQFunction, __FILE__":IsotropicLinearElasticityQFunction", &qf_physics);
	CeedQFunctionAddInput(qf_physics, "disp", dim, CEED_EVAL_INTERP);
	CeedQFunctionAddInput(qf_physics, "graddisp", dim*dim, CEED_EVAL_NONE);
	CeedQFunctionAddOutput(qf_physics, "v", dim, CEED_EVAL_INTERP);
	CeedQFunctionAddOutput(qf_physics, "dv", dim*dim, CEED_EVAL_GRAD);


	///operator for setup
	CeedOperatorCreate(ceed, qf_setup, NULL, NULL, &op_setup);
	CeedOperatorSetField(op_setup, "J", restrictx, CEED_NOTRANSPOSE, basis, CEED_VECTOR_ACTIVE);
	CeedOperatorSetField(op_setup, "w", restrictxi, CEED_NOTRANSPOSE, basis, CEED_VECTOR_NONE);
	CeedOperatorSetField(op_setup, "qdata", restrictqdi, CEED_NOTRANSPOSE, CEED_BASIS_COLLOCATED, CEED_VECTOR_ACTIVE);

	///operator for mass
	CeedOperatorCreate(ceed, qf_mass, NULL, NULL, &op_mass);
	CeedOperatorSetField(op_mass, "u", restrictq, CEED_TRANSPOSE, basis, CEED_VECTOR_ACTIVE);
	CeedOperatorSetField(op_mass, "qdata", restrictqdi, CEED_NOTRANSPOSE, basis, qdata);
	CeedOperatorSetField(op_mass, "v", restrictq, CEED_TRANSPOSE, basis, CEED_VECTOR_ACTIVE);

	///operator for physics
	CeedOperatorCreate(ceed, qf_physics, NULL, NULL, &op);
	CeedOperatorSetField(op, "disp", restrictq, CEED_TRANSPOSE, basis, CEED_VECTOR_ACTIVE);
	CeedOperatorSetField(op, "graddisp", restrictq, CEED_TRANSPOSE, basis, CEED_VECTOR_ACTIVE);
	CeedOperatorSetField(op, "v", restrictq, CEED_TRANSPOSE, basis, CEED_VECTOR_ACTIVE);
	CeedOperatorSetField(op, "dv", restrictq, CEED_TRANSPOSE, basis, CEED_VECTOR_ACTIVE);


	//setup context todo: don't hardcode these!
	setupcontext.dim = dim;
	physicscontext.dim = dim;
	physicscontext.BM = 1;
	physicscontext;


  CeedQFunctionSetContext(qf_setup, &setupcontext, sizeof setupcontext);
  CeedQFunctionSetContext(qf_mass, &setupcontext, sizeof setupcontext);
  CeedQFunctionSetContext(qf_physics, &physicscontext, sizeof physicscontext);

  printf("after setup\n");
	//todo: apply setup and mass operators
	CeedOperatorApply(op_setup,localcoordsceed, qdata, CEED_REQUEST_IMMEDIATE);
	CeedOperatorApply(op_mass, onesvec, mceed, CEED_REQUEST_IMMEDIATE);


	//clean up
	CeedBasisDestroy(&basis);
  printf("A\n");
	CeedElemRestrictionDestroy(&restrictx);
	CeedElemRestrictionDestroy(&restrictxi);	
	CeedElemRestrictionDestroy(&restrictq);
	CeedElemRestrictionDestroy(&restrictqdi);
  printf("B\n");
	CeedQFunctionDestroy(&qf_setup);
	CeedQFunctionDestroy(&qf_mass);
	CeedQFunctionDestroy(&qf_physics);
	CeedOperatorDestroy(&op_setup);
	CeedOperatorDestroy(&op_mass);
  printf("C\n");
	CeedDestroy(&ceed);
  printf("end of initialize\n");

}; // initialize


void
pylith::feassemble::IntegratorDomainCEED::_computeResidual(pylith::topology::Field* residual,
                                                       const std::vector<ResidualKernels>& kernels,
                                                       const PylithReal t,
                                                       const PylithReal dt,
                                                       const pylith::topology::Field& solution,
                                                       const pylith::topology::Field& solutionDot) 
{
  printf("in compute residual\n");
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_computeResidual(residual="<<residual<<", # kernels="<<kernels.size()<<", t="<<t<<", dt="<<dt<<", solution="<<solution.label()<<", solutionDot="<<solutionDot.label()<<")");

    assert(residual);
    assert(_auxiliaryField);

    PetscDS prob = NULL;
    PetscIS cells = NULL;
    PetscInt cStart = 0, cEnd = 0;
    PetscErrorCode err;

    PetscDM dmSoln = solution.dmMesh();
    PetscDM dmAux = _auxiliaryField->dmMesh();
    PetscDMLabel dmLabel = NULL;

    // :KLUDGE: Potentially we may have multiple PetscDS objects. This assumes that the first one (with a NULL label) is
    // the correct one.

    // Set pointwise function (kernels) in DS
    err = DMGetDS(dmSoln, &prob);PYLITH_CHECK_ERROR(err);
    for (size_t i = 0; i < kernels.size(); ++i) {
        const PetscInt i_field = solution.subfieldInfo(kernels[i].subfield.c_str()).index;
        err = PetscDSSetResidual(prob, i_field, kernels[i].r0, kernels[i].r1);PYLITH_CHECK_ERROR(err);
    } // for

    // Get auxiliary data
    err = PetscObjectCompose((PetscObject) dmSoln, "dmAux", (PetscObject) dmAux);PYLITH_CHECK_ERROR(err);
    err = PetscObjectCompose((PetscObject) dmSoln, "A", (PetscObject) _auxiliaryField->localVector());PYLITH_CHECK_ERROR(err);

    // Compute the local residual
    assert(solution.localVector());
    assert(residual->localVector());
    err = DMGetLabel(dmSoln, "material-id", &dmLabel);PYLITH_CHECK_ERROR(err);
    err = DMLabelGetStratumBounds(dmLabel, _materialId, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
    assert(cEnd > cStart); // Double-check that this material has cells.

    PYLITH_JOURNAL_DEBUG("DMPlexComputeResidual_Internal() with material-id '"<<_materialId<<"' and cells ["<<cStart<<","<<cEnd<<").");
    err = ISCreateStride(PETSC_COMM_SELF, cEnd-cStart, cStart, 1, &cells);PYLITH_CHECK_ERROR(err);
   // err = DMPlexComputeResidual_Internal(dmSoln, cells, PETSC_MIN_REAL, solution.localVector(), solutionDot.localVector(), residual->localVector(), NULL);PYLITH_CHECK_ERROR(err);

    Vec Gloc, Qloc;
    PetscScalar *q, *g;
    CeedVector qceed, gceed;

    //local in and out vectors
    Qloc = solution.localVector();
    Gloc = solution.localVector();
    
    err = VecZeroEntries(Gloc); PYLITH_CHECK_ERROR(err);


    //CEED vectors
 	err = VecGetArrayRead(Qloc, (const PetscScalar**)&q); PYLITH_CHECK_ERROR(err);
 	err = VecGetArray(Gloc, &g); PYLITH_CHECK_ERROR(err);
 	CeedVectorSetArray(qceed, CEED_MEM_HOST, CEED_USE_POINTER, q);
 	CeedVectorSetArray(gceed, CEED_MEM_HOST, CEED_USE_POINTER, g);

 	//apply
    CeedOperatorApply(op, qceed, gceed, CEED_REQUEST_IMMEDIATE);



    err = ISDestroy(&cells);PYLITH_CHECK_ERROR(err);

    //todo: apply physics operator


    PYLITH_METHOD_END;
} // _computeResidual




// End of file
