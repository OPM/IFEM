//==============================================================================
//!
//! \file FieldValueIntegrand.h
//!
//! \date Oct 22 2026
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Integrand for testing the secondary solution evaluation of a patch.
//!
//==============================================================================

#ifndef FIELD_VALUE_INTEGRAND_H
#define FIELD_VALUE_INTEGRAND_H

#include "FiniteElement.h"
#include "IntegrandBase.h"
#include "Vec3.h"


/*!
  \brief A secondary solution which is the value of a nodal field.
  \details The second component is the element of the evaluation point.
  The first is in the spline space, so it is reproduced by a projection.
*/

class FieldValueIntegrand : public IntegrandBase
{
public:
  //! \brief The constructor sets the nodal field.
  FieldValueIntegrand(unsigned short int n, const Vector& field)
    : IntegrandBase(n), nodal(field) {}

  //! \brief Returns the number of secondary solution components.
  size_t getNoFields(int fld) const override { return fld < 2 ? 1 : 2; }

  using IntegrandBase::evalSol;
  //! \brief Evaluates the secondary solution at a result point.
  bool evalSol(Vector& s, const FiniteElement& fe, const Vec3&,
               const std::vector<int>& MNPC) const override
  {
    s.resize(2);
    s[0] = 0.0;
    for (size_t i = 0; i < MNPC.size(); i++)
      s[0] += fe.N(1+i)*nodal(1+MNPC[i]);
    s[1] = fe.iel;
    return true;
  }

private:
  Vector nodal; //!< The nodal field
};

#endif
