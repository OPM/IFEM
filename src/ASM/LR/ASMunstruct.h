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
}
namespace LR {
  class LRSplineSurface;
}


/*!
  \brief Base class for unstructured spline-based FE assembly drivers.
  \details This class contains methods common for unstructured spline patches.
*/

class ASMunstruct : public ASMbase
{
protected:
  //! \brief The constructor sets the space dimensions.
  //! \param[in] n_p Number of parameter dimensions
  //! \param[in] n_s Number of spatial dimensions
  //! \param[in] n_f Number of primary solution fields
  ASMunstruct(unsigned char n_p, unsigned char n_s, unsigned char n_f);
  //! \brief Copy constructor.
  //! \param[in] patch The patch to use FE data from
  //! \param[in] n_f Number of primary solution fields
  ASMunstruct(const ASMunstruct& patch, unsigned char n_f = 0);

public:
  //! \brief The destructor frees the dynamically allocated spline object.
  virtual ~ASMunstruct();

  //! \brief Checks if the patch is empty.
  virtual bool empty() const { return geo == 0; }

  //! \brief Defines the numerical integration scheme to use.
  //! \param[in] ng Number of Gauss points in each parameter direction
  void setGauss(int ng) { nGauss = ng; }

  //! \brief Resets global element and node counters.
  static void resetNumbering() { gEl = gNod = 0; }

  //! \brief Computes the total number of integration points in this patch.
  virtual void getNoIntPoints(size_t& nPt) { nPt = 0; } // later...
  //! \brief Computes the number of boundary integration points in this patch.
  virtual void getNoBouPoints(size_t& nPt, int, int) { nPt = 0; } // later...

protected:
  LR::LRSplineSurface* geo; //!< Pointer to the actual spline geometry object

  int    nGauss;   //!< Numerical integration scheme
  static int gEl;  //!< Global element counter
  static int gNod; //!< Global node counter

};

#endif

