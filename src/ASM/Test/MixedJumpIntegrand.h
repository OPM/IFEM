// $Id$
//==============================================================================
//!
//! \file MixedJumpIntegrand.h
//!
//! \date Oct 25 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand for tests of interface integrals of mixed patches.
//!
//==============================================================================

#ifndef _MIXED_JUMP_INTEGRAND_H
#define _MIXED_JUMP_INTEGRAND_H

#include "ASMbase.h"
#include "ElmMats.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "TimeDomain.h"
#include "Vec3.h"

#include "Catch2Support.h"

#include <cmath>
#include <map>


/*!
  \brief Integrates the squared jump of the normal derivative of a field on
  the first basis of a mixed element over its interfaces.
*/

class MixedJumpIntegrand : public IntegrandBase
{
public:
  //! \brief The constructor sets the number of spatial dimensions.
  explicit MixedJumpIntegrand(unsigned short int n) : IntegrandBase(n) {}

  //! \brief Returns the integrand type.
  int getIntegrandType() const override { return INTERFACE_TERMS; }

  using IntegrandBase::getLocalIntegral;
  //! \brief Returns a scalar element vector.
  LocalIntegral* getLocalIntegral(const std::vector<size_t>&,
                                  size_t, bool) const override
  {
    ElmMats* result = new ElmMats();
    result->resize(0,1);
    result->b[0].resize(1);
    return result;
  }

  using IntegrandBase::initElement;
  //! \brief Sets the field values of the element nodes.
  bool initElement(const IntVec& MNPC, const MxFiniteElement&,
                   const std::vector<size_t>&, const std::vector<size_t>&,
                   LocalIntegral& elmInt) override
  {
    Vector v;
    for (int n : MNPC)
      v.push_back(std::sin(1.0 + 0.3*n) + 0.2*std::cos(0.7*n));
    elmInt.vec = { v };
    return true;
  }

  using IntegrandBase::evalIntMx;
  //! \brief Evaluates the integrand at an interface point.
  bool evalIntMx(LocalIntegral& elmInt, const MxFiniteElement& fe,
                 const Vec3&, const Vec3& normal) const override
  {
    // The bases of this element, then those of the neighbor
    const size_t nB = fe.getNoBasis()/2;
    double jump = 0.0;
    for (size_t side = 0; side < 2; side++)
    {
      const Vector& v = elmInt.vec[side];
      const Matrix& dNdX = fe.grad(1 + side*nB);
      for (size_t n = 1; n <= dNdX.rows(); n++)
        for (size_t k = 1; k <= dNdX.cols(); k++)
          jump += (side ? -1.0 : 1.0) * dNdX(n,k)*v(n)*normal[k-1];
    }

    static_cast<ElmMats&>(elmInt).b[0](1) += jump*jump*fe.detJxW;
    return true;
  }
};


/*!
  \brief Keeps the scalar element integrals by element number.
*/

class ElementValues : public GlobalIntegral
{
public:
  //! \brief Stores the element value.
  bool assemble(const LocalIntegral* elmObj, int elmId) override
  {
    values[elmId] = static_cast<const ElmMats*>(elmObj)->b[0](1);
    return true;
  }

  std::map<int,double> values; //!< The element values
};


//! \brief Checks that a continuously differentiable field has no jumps in
//! its normal derivative over the element interfaces of a mixed patch.
inline void checkNoJumps(ASMbase& pch, const ASM::InterfaceChecker& iChk)
{
  MixedJumpIntegrand itg(pch.getNoSpaceDim());
  ElementValues jumps;
  REQUIRE(pch.integrate(itg,jumps,TimeDomain(),iChk));
  REQUIRE(jumps.values.size() == pch.getNoElms());
  for (const std::pair<const int,double>& elm : jumps.values)
  {
    INFO("element " << elm.first);
    CHECK(elm.second < 1e-20);
  }
}

#endif
