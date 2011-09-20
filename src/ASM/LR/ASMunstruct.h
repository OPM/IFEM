// $Id$
//==============================================================================
//!
//! \file ASMunstruct.h
//!
//! \date September 2011
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for unstructured spline-based FE assembly drivers.
//!
//==============================================================================

#ifndef _ASM_UNSTRUCT_H
#define _ASM_UNSTRUCT_H

#include "ASMbase.h"

namespace Go {
  class GeomObject;
  class BoundingBox;
}
namespace LR {
  class LRSplineSurface;
}

/*!
  \brief Base class for structured spline-based FE assembly drivers.
  \details This class contains methods common for structured spline patches.
*/

class ASMunstruct : public ASMbase
{
protected:
  //! \brief The constructor sets the space dimensions.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMunstruct(unsigned char n_p, unsigned char n_s, unsigned char n_f);

public:
  //! \brief The destructor frees the dynamically allocated spline object.
  virtual ~ASMunstruct();

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geo == 0; }

  //! \brief Defines the numerical integration scheme to use.
  //! \param[in] ng Number of Gauss points in each parameter direction
  void setGauss(int ng) { nGauss = ng; }

  //! \brief Resets global element and node counters
  static void resetNumbering() { gEl = gNod = 0; }

protected:
  LR::LRSplineSurface* geo; //!< Pointer to the actual spline geometry object

  int        nGauss;                        //!< Numerical integration scheme
  static int gEl;                           //!< Global element counter
  static int gNod;                          //!< Global node counter

};

#endif

