// $Id$
//==============================================================================
//!
//! \file HasGravityBase.h
//!
//! \date May 27 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class representing FEM integrands with gravity terms.
//!
//==============================================================================

#ifndef _HAS_GRAVITY_BASE_H
#define _HAS_GRAVITY_BASE_H

#include "IntegrandBase.h"
#include "Vec3.h"


/*!
  \brief Base class representing a FEM integrand with gravity terms.
*/

class HasGravityBase : public IntegrandBase
{
protected:
  //! \brief The default constructor is protected to allow sub-classes only.
  explicit HasGravityBase(unsigned char n = 0) : IntegrandBase(n) {}

public:
  //! \brief Parses a data section from an XML-element.
  virtual bool parse(const tinyxml2::XMLElement* elem);

  //! \brief Defines the gravitation vector.
  void setGravity(double gx, double gy = 0.0, double gz = 0.0)
  { gravity.x = gx; gravity.y = gy; gravity.z = gz; }
  //! \brief Defines the gravitation vector.
  void setGravity(const Vec3& g) { gravity = g; }

  //! \brief Returns the gravitation vector.
  const Vec3& getGravity() const { return gravity; }

protected:
  Vec3 gravity; //!< Gravitation vector
};

#endif
