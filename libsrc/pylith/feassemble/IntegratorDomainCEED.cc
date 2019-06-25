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

#include "spatialdata/spatialdb/GravityField.hh" // HASA GravityField
#include "petscds.h" // USES PetscDS

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL_*

#include <cassert> // USES assert()
#include <stdexcept> // USES std::runtime_error


// ---------------------------------------------------------------------------------------------------------------------
// Initialize integration domain, auxiliary field, and derived field. Update observers.
void
pylith::feassemble::IntegratorDomainCEED::initialize(const pylith::topology::Field& solution) {

    printf("initialize IntegratorDomainCEED");

} // initialize




// End of file
