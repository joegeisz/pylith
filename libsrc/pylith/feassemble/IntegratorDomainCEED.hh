// -*- C++ -*-
//
// ======================================================================
//
//Joe Geisz
// 
// IntegratorDomain with CEED implementation
//
// ======================================================================
//

/**
 * @file libsrc/feassemble/IntegratorDomainCEED.hh
 *
 * @brief Object for finite-element integration over a subset (material) of the simulation domain.
 *        Uses libCEED, Inherits from normal IntegratorDomain
 */

#if !defined(pylith_feassemble_integratordomain_ceed_hh)
#define pylith_feassemble_integratordomain_ceed_hh

#include "pylith/feassemble/IntegratorDomain.hh" //
#include <ceed.h>
#include "pylith/fekernels/IsotropicLinearElasticityCEED.hh"//uses context classes


class pylith::feassemble::IntegratorDomainCEED : public pylith::feassemble::IntegratorDomain {
    friend class TestIntegratorDomain; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

	/// Constructor
    IntegratorDomainCEED(pylith::problems::Physics* const physics);

    /// Destructor
    ~IntegratorDomainCEED(void);

    /// Initialize
    void initialize(const pylith::topology::Field& solution);


    //
    void _computeResidual(pylith::topology::Field* residual,
				    	const std::vector<ResidualKernels>& kernels, 
				    	const PylithReal t,
				    	const PylithReal dt,
				    	const pylith::topology::Field& solution,
				    	const pylith::topology::Field& solutionDot);

    //Ceed objects
    CeedOperator op;

	//contexts
	CoordinatesContext setupcontext;
	ProblemContext physicscontext;

}; // IntegratorDomainCEED

#endif // pylith_feassemble_integratordomain_ceed_hh

// End of file
