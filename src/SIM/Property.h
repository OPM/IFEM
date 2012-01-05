// $Id$
//==============================================================================
//!
//! \file Property.h
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of a distributed physical property.
//!
//==============================================================================

#ifndef _PROPERTY_H
#define _PROPERTY_H

#include <vector>


/*!
  \brief Struct for representing a distributed physical property.
*/

struct Property
{
  //! \brief The available property types.
  enum Type
  {
    UNDEFINED,
    MATERIAL,
    BODYLOAD,
    NEUMANN,
    DIRICHLET,
    DIRICHLET_INHOM
  };

  Type   pcode; //!< Physical property code
  size_t pindx; //!< Physical property index
  size_t patch; //!< Patch index [0,nPatch>
  char   lindx; //!< Local entity index which is assigned the property
  char   ldim;  //!< Local entity dimension flag [0,3]

  //! \brief Default constructor.
  Property() : pcode(UNDEFINED), pindx(0), patch(0), lindx(0), ldim(0) {}

  //! \brief Constructor creating an initialized property instance.
  Property(Type t, size_t px, size_t p, char ld, char lx = 0) :
    pcode(t), pindx(px), patch(p), lindx(lx), ldim(ld) {}
};

typedef std::vector<Property> PropertyVec; //!< Vector of properties

#endif
