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

  //! \brief Updates the nodal coordinates for this patch.
  //! \param[in] displ Incremental displacements to update the coordinates with
  //! \param[in] nsd Number of space dimensions
  bool updateCoords(const Vector& displ, unsigned char nsd);

public:
  //! \brief Updates patch origin by adding a constant to all nodes.
  //! \param[in] origin The new origin of the patch
  void updateOrigin(const Vec3& origin);

protected:
  const Vec3Vec& coord; //!< Const reference to the nodal coordinate container
  Vec3Vec      myCoord; //!< Nodal coordinate container
};

#endif
