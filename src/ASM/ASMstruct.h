// $Id: ASMstruct.h,v 1.8 2010-12-29 18:02:10 kmo Exp $
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

public:
  //! \brief The destructor frees the dynamically allocated spline object.
  virtual ~ASMstruct();

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geo == 0; }

  //! \brief Defines the numerical integration scheme to use.
  //! \param[in] ng Number of Gauss points in each parameter direction
  void setGauss(int ng) { nGauss = ng; }

  //! \brief Resets the global element and node counters.
  static void resetNumbering() { gEl = gNod = 0; }

protected:
  Go::GeomObject* geo; //!< Pointer to the actual spline geometry object

  int        nGauss;   //!< Numerical integration scheme
  static int gEl;      //!< Global element counter
  static int gNod;     //!< Global node counter
};

#endif
