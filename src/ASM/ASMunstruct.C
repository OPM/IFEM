// $Id$
//==============================================================================
//!
//! \file ASMunstruct.C
//!
//! \date December 2010
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Base class for unstructured spline-based FE assembly drivers.
//!
//==============================================================================

#include "ASMunstruct.h"

#ifndef HAS_LRSPLINE
namespace LR {
  //! \brief Dummy class, referred only when built without LR-Spline support.
  class LRSpline {};
}
#endif


int ASMunstruct::gEl = 0;
int ASMunstruct::gNod = 0;


ASMunstruct::ASMunstruct (unsigned char n_p, unsigned char n_s,
                          unsigned char n_f)
  : ASMbase(n_p,n_s,n_f)
{
  geo = nullptr;
}


ASMunstruct::ASMunstruct (const ASMunstruct& patch, unsigned char n_f)
  : ASMbase(patch,n_f)
{
  nGauss = patch.nGauss;
  geo = patch.geo;
}


#ifndef HAS_LRSPLINE

// Dummy implementations, referred only when compiled without LR-Spline library.
// The real methods are found in the file LR/ASMLRSpline.C

bool ASMunstruct::refine (const LR::RefineData&, Vectors&, const char*)
{
  std::cerr <<" *** ASMunstruct::refine: Not available without LR-Spline"
            << std::endl;
  return false;
}

#endif
