//==============================================================================
//!
//! \file LaplacianIntegrand.h
//!
//! \date Oct 22 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand for testing the second derivatives of a spline patch.
//!
//==============================================================================

#ifndef LAPLACIAN_INTEGRAND_H
#define LAPLACIAN_INTEGRAND_H

#include "ASMbase.h"
#include "FiniteElement.h"
#include "IntegrandBase.h"
#include "LocalIntegral.h"
#include "Vec3.h"

#include <cmath>
#include <vector>


/*!
  \brief Integrates the square of the Laplacian of a field over a patch.
  \details The field has the value of a smooth function at the control point
  of each node, so the result does not depend on how the patch numbers its
  nodes, and two patches with the same basis can be compared.
*/

class LaplacianIntegrand : public IntegrandBase
{
  //! \brief The field values at the nodes of an element.
  struct Values : public LocalIntegral
  {
    std::vector<double> u; //!< Field value at each element node
  };

public:
  //! \brief The constructor sets the patch to evaluate the field on.
  explicit LaplacianIntegrand(const ASMbase& patch)
    : IntegrandBase(patch.getNoSpaceDim()), pch(patch) {}

  //! \brief Second derivatives are needed.
  int getIntegrandType() const override { return SECOND_DERIVATIVES; }

  //! \brief Returns a local integral for an element.
  LocalIntegral* getLocalIntegral(size_t, size_t, bool) const override
  { return new Values(); }

  using IntegrandBase::initElement;
  //! \brief Sets up the field values at the nodes of an element.
  bool initElement(const std::vector<int>& MNPC, const FiniteElement&,
                   const Vec3&, size_t, LocalIntegral& elmInt) override
  {
    return this->initElement(MNPC,elmInt);
  }
  //! \brief Sets up the field values at the nodes of an element.
  bool initElement(const std::vector<int>& MNPC,
                   LocalIntegral& elmInt) override
  {
    Values& v = static_cast<Values&>(elmInt);
    v.u.clear();
    for (int n : MNPC)
    {
      const Vec3 X = pch.getCoord(1+n);
      v.u.push_back(std::sin(X.x) + X.x*X.y*X.y + X.y*X.z);
    }
    return true;
  }

  using IntegrandBase::evalInt;
  //! \brief Adds the square of the Laplacian at an integration point.
  bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
               const Vec3&) const override
  {
    const Values& v = static_cast<const Values&>(elmInt);
    double lap = 0.0;
    for (size_t i = 1; i <= fe.N.size(); i++)
      for (size_t k = 1; k <= pch.getNoSpaceDim(); k++)
        lap += v.u[i-1]*fe.d2NdX2(i,k,k);
#pragma omp atomic
    result += lap*lap*fe.detJxW;
    return true;
  }

  mutable double result = 0.0; //!< The integral

private:
  const ASMbase& pch; //!< The patch to evaluate the field on
};

#endif
