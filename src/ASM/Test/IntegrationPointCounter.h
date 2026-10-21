//==============================================================================
//!
//! \file IntegrationPointCounter.h
//!
//! \date Oct 21 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand counting the interior integration points of a patch.
//!
//==============================================================================

#ifndef _INTEGRATION_POINT_COUNTER_H
#define _INTEGRATION_POINT_COUNTER_H

#include "ASMbase.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "TimeDomain.h"

#include "Catch2Support.h"

#include <cmath>


/*!
  \brief Integrand counting the interior integration points.
*/

class IntegrationPointCounter : public IntegrandBase
{
public:
  //! \brief The constructor sets the number of spatial dimensions.
  explicit IntegrationPointCounter(unsigned short int n) : IntegrandBase(n) {}

  using IntegrandBase::evalInt;
  //! \brief Counts the integration point.
  bool evalInt(LocalIntegral&, const FiniteElement&, const Vec3&) const override
  {
#pragma omp atomic
    ++count;
    return true;
  }

  mutable size_t count = 0; //!< Number of integration points
};


/*!
  \brief Checks that a patch integrates with its current quadrature rule.
  \details The quadrature rule is changed between the integrations, as for
  the solution norms, with each cache policy.
*/

template<class Patch>
void checkQuadratureChange(size_t ndim)
{
  const ASM::CachePolicy policy = ASM::cachePolicy;
  for (ASM::CachePolicy cp : {ASM::NO_CACHE, ASM::PRE_CACHE,
                              ASM::ON_THE_FLY, ASM::FULL_CACHE})
  {
    ASM::cachePolicy = cp;
    ASMbase::resetNumbering();
    Patch pch;
    for (size_t d = 0; d < ndim; d++)
      REQUIRE(pch.uniformRefine(d,1));
    REQUIRE(pch.generateFEMTopology());

    IntegrationPointCounter itg(ndim);
    GlobalIntegral dummy;
    static_cast<ASMbase&>(pch).generateThreadGroups(itg,true,false);
    for (int ng : {3, 5, 3})
    {
      INFO("cache policy " << cp << " nGauss " << ng);
      pch.setGauss(ng);
      itg.count = 0;
      REQUIRE(pch.integrate(itg,dummy,TimeDomain()));
      CHECK(itg.count == pch.getNoElms()*size_t(std::pow(ng,ndim)));
    }
  }
  ASM::cachePolicy = policy;
}

#endif
