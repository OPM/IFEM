// $Id$
//==============================================================================
//!
//! \file GlbNorm.C
//!
//! \date Dec 09 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of integrated norm quantities over an element.
//!
//==============================================================================

#include "GlbNorm.h"
#include "ElmNorm.h"
#include <cmath>


GlbNorm::~GlbNorm ()
{
  for (size_t i = 0; i < myVals.size(); i++)
    for (size_t j = 0; j < myVals[i].size(); j++)
      this->applyFinalOp(myVals[i][j]);
}


bool GlbNorm::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmNorm* ptr = dynamic_cast<const ElmNorm*>(elmObj);
  if (!ptr) return false;

  // If the element norms are requested (i.e. the internal buffer
  // is not used) the actuall summation of element norms into the
  // global norm is postponed to the very end (when invoked with
  // elmId=0). This must be done like this to allow more than one
  // loop over the elements during the norm integration.
  if (elmId > 0 && ptr->externalStorage()) return true;

  ElmNorm& elVals = *const_cast<ElmNorm*>(ptr);
  size_t i, j, k;
  for (i = j = k = 0; i < elVals.size(); i++)
  {
    if (j >= myVals[k].size())
      k++, j = 0;
    if (k < myVals.size() && j < myVals[k].size())
      myVals[k][j++] += elVals[i];
    this->applyFinalOp(elVals[i]);
  }

  return true;
}


void GlbNorm::applyFinalOp (double& value) const
{
  switch (myOp)
    {
    case ASM::ABS:
      value = fabs(value);
      break;

    case ASM::SQRT:
      if (value < 0.0)
	value = -sqrt(-value);
      else
	value = sqrt(value);
      break;

    default:
      break;
    }
}
