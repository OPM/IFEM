// $Id$
//==============================================================================
//!
//! \file ASMLagBase.h
//!
//! \date May 1 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Common base class for %Lagrange FE models.
//!
//==============================================================================

#ifndef _ASM_LAG_BASE_H
#define _ASM_LAG_BASE_H

#include "MatVec.h"
#include "Vec3.h"


/*!
  \brief Common base class for %Lagrange FE models.
  \details This class contains a container for storing the nodal coordinates.
*/

class ASMLagBase
{
protected:
  //! \brief Default constructor.
  ASMLagBase() : coord(myCoord) {}
  //! \brief Copy constructor.
  ASMLagBase(const ASMLagBase& patch, bool cpCoord = true);
  //! \brief Empty destructor.
  virtual ~ASMLagBase() {}

  //! \brief Direct nodal evaluation of a solution field.
  bool nodalField(Matrix& field, const Vector& sol, size_t nno) const;

  //! \brief Returns the geometric center of an element.
  Vec3 getGeometricCenter(const std::vector<int>& MNPC) const;

protected:
  const Vec3Vec& coord; //!< Const reference to the nodal coordinate container
  Vec3Vec      myCoord; //!< Nodal coordinate container
};

#endif
