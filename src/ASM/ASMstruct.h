// $Id$
//==============================================================================
//!
//! \file ASMstruct.h
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for structured spline-based FE assembly drivers.
//!
//==============================================================================

#ifndef _ASM_STRUCT_H
#define _ASM_STRUCT_H

#include "ASMbase.h"

namespace Go {
  class GeomObject;
}

/*!
  \brief Base class for structured spline-based FE assembly drivers.
  \details This class contains methods common for structured spline patches.
*/

class ASMstruct : public ASMbase
{
protected:
  //! \brief The constructor sets the number of space dimensions and fields.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMstruct(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Copy constructor.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  ASMstruct(const ASMstruct& patch, unsigned char n_f = 0);

public:
  //! \brief The destructor frees the dynamically allocated spline object.
  virtual ~ASMstruct();

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geo == 0; }

  //! \brief Resets the global element and node counters.
  static void resetNumbering(int nnod = 0);

  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integr Object with problem-specific data and methods
  virtual Go::GeomObject* evalSolution(const IntegrandBase& integr) const = 0;

protected:
  Go::GeomObject* geo; //!< Pointer to the actual spline geometry object

  //! Auxilliary node number map used when establishing Dirichlet constraints
  static std::map<int,int> xNode;

  static int gEl;  //!< Global element counter
  static int gNod; //!< Global node counter
};

#endif
