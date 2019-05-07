// $Id$
//==============================================================================
//!
//! \file ASMunstruct.h
//!
//! \date September 2011
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Abstract interface for unstructured FE assembly drivers.
//!
//==============================================================================

#ifndef _ASM_UNSTRUCT_H
#define _ASM_UNSTRUCT_H

#include "MatVec.h"
#include <set>
#include <string>

typedef std::vector<int> IntVec; //!< General integer vector
typedef std::set<int>    IntSet; //!< General integer set

class RealFunc;


namespace LR //! Utilities for LR-splines.
{
  /*!
    \brief A struct of data to control the mesh refinement.
    \details The \a options parameters have the following interpretation:
    - options[0] : the beta percentage of elements to refine
    - options[1] : the knotline multiplicity (default 1)
    - options[2] : the refinement scheme (default 0),
                   see AdaptiveSetup::RefScheme
    - options[3] : nonzero if testing for linear independence at all iterations
    - options[4] : the maximum number of T-joints allowed in the model
    - options[5] : the maximum allowed parametric aspect ratio of an element
    - options[6] : one if all "gaps" are to be closed
    - options[7] : one if using "true beta"
  */
  struct RefineData
  {
    bool      refShare; //!< If \e true, force refinement of shared FE grids
    IntVec    options;  //!< Parameters used to control the refinement
    IntVec    elements; //!< 0-based indices of the elements to refine
    RealArray errors;   //!< List of error indicators for the elements

    //! \brief Default constructor.
    explicit RefineData(bool rs = false) : refShare(rs) {}
    //! \brief Clears the refinement parameters.
    void clear() { options.clear(); elements.clear(); errors.clear(); }
  };
}


/*!
  \brief Abstract interface for unstructured spline patches.
*/

class ASMunstruct
{
protected:
  //! \brief The constructor is protected to allow objects of sub-classes only.
  ASMunstruct() {}

public:
  //! \brief Empty destructor.
  virtual ~ASMunstruct() {}

  //! \brief Refines the mesh adaptively.
  //! \param[in] prm Input data used to control the mesh refinement
  //! \param sol Control point results values that are transferred to new mesh
  virtual bool refine(const LR::RefineData& prm, Vectors& sol) = 0;

  //! \brief Remaps element-wise errors from geometry mesh to refinement mesh.
  //! \param[out] errors The remapped errors
  //! \param[in] orig The element-wise errors on the geometry mesh
  virtual void remapErrors(RealArray& errors, const RealArray& orig,
                           bool = false) const { errors = orig; }

  //! \brief Refines the parametrization based on a mesh density function.
  //! \param[in] refC Mesh refinement criteria function
  //! \param[in] refTol Mesh refinement threshold
  virtual bool refine(const RealFunc& refC, double refTol) = 0;

  //! \brief Returns all boundary functions that are covered by the given nodes.
  virtual IntVec getBoundaryCovered(const IntSet&) const { return IntVec(); }
  //! \brief Extends the refinement domain with information for neighbors.
  virtual void extendRefinementDomain(IntSet&, const IntSet&) const {}

  //! \brief Stores the mesh basis to encapsulated postscript files.
  //! \param[in] fName Prefix for file names
  //! \param[in] fType Flag telling which file type(s) to write (15 means all)
  virtual void storeMesh(const std::string& fName, int fType = 15) const {}
};

#endif
