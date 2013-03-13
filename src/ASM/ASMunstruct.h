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

namespace LR {
  class LRSpline;
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

  //! \brief Refines a specified list of elements.
  //! \param[in] elements 0-based indices of the elements to refine
  //! \param[in] options Additional input parameters to control the refinement,
  //! options[0] is the beta percentage of elements to refine,
  //! options[1] is the knotline multiplicity (default 1),
  //! options[2] is the refinement scheme (default 0),
  //! (FULLSPAN=0, MINSPAN=1, ISOTROPIC ELEMENTS=2, ISOTROPIC FUNCTIONS=3),
  //! options[3] is the symmetry, i.e., always refine a multiple of this value
  //! options[4] is nonzero if testing for linear independence at all iterations
  //! options[5] is the maximum number of T-joints allowed in the model
  //! options[6] is the maximum allowed parametric aspect ratio of an element
  //! options[7] is one if all "gaps" are to be closed
  //! options[8] is one if using true beta
  //! \param[in] fName Optional file name for an image of the resulting mesh
  bool refine(const IntVec& elements, const IntVec& options,
              const char* fName = NULL);
  //! \brief Refines a set of elements based on a list of element errors.
  //! \param[in] elementError Element-wise errors
  //! \param[in] options Additional input parameters to control the refinement
  //! \param[in] fName Optional file name for an image of the resulting mesh
  bool refine(const RealArray& elementError, const IntVec& options,
              const char* fName = NULL);

  //! \brief Resets global element and node counters.
  static void resetNumbering() { gEl = gNod = 0; }

  //! \brief Projects the secondary solution field onto the primary basis.
  //! \param[in] integrand Object with problem-specific data and methods
  virtual LR::LRSpline* evalSolution(const IntegrandBase& integrand) const = 0;

protected:
  LR::LRSpline* geo; //!< Pointer to the actual spline geometry object

  static int gEl;  //!< Global element counter
  static int gNod; //!< Global node counter
};

#endif
