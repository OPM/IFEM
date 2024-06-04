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
#include <cstddef>


/*!
  \brief Struct for representing a distributed physical property.
*/

struct Property
{
  //! \brief The available property types.
  //! \note The sequence of these enums are of importance, since the less-than
  //! and greater-than operators are used on instances of them. Therefore,
  //! do *not* alter the below order unless you know what you are doing.
  enum Type
  {
    UNDEFINED,
    MATERIAL,
    BODYLOAD,
    NEUMANN,
    NEUMANN_ANASOL,
    NEUMANN_GENERIC,
    ROBIN,
    ROBIN_ANASOL,
    RIGID,
    DIRICHLET,
    DIRICHLET_INHOM,
    DIRICHLET_ANASOL,
    OTHER
  };

  Type   pcode; //!< Physical property code
  int    pindx; //!< Physical property index (0-based)
  size_t patch; //!< Patch index [1,nPatch]
  char   lindx; //!< Local entity index (1-based) which is assigned the property
  char   ldim;  //!< Local entity dimension flag [0,3]
  char   basis; //!< Which basis the property is defined on

  //! \brief Default constructor.
  Property() : pcode(UNDEFINED), pindx(0), patch(0) { lindx=ldim = basis = 0; }

  //! \brief Constructor creating an initialized property instance.
  Property(Type t, int px, size_t p, char ld, char lx = 0, char b = 0) :
    pcode(t), pindx(px), patch(p), lindx(lx), ldim(ld), basis(b) {}
};

typedef std::vector<Property> PropertyVec; //!< Vector of properties

#endif
