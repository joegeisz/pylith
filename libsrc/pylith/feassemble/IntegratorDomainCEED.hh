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

class pylith::feassemble::IntegratorDomainCEED : public pylith::feassemble::IntegratorDomain {
    friend class TestIntegratorDomain; // unit testing

    // PUBLIC MEMBERS //////////////////////////////////////////////////////////////////////////////////////////////////
public:

	/// Constructor
    IntegratorDomainCEED(pylith::problems::Physics* const physics);

    /// Initialize
    void initialize(const pylith::topology::Field& solution);

}; // IntegratorDomainCEED

#endif // pylith_feassemble_integratordomain_ceed_hh

// End of file
