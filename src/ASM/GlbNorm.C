// $Id: GlbNorm.C,v 1.2 2010-01-16 16:02:54 kmo Exp $
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
#include <math.h>


GlbNorm::GlbNorm (std::vector<double>& vec, FinalOp op) : myVals(vec)
{
  myOp = op;
}


GlbNorm::~GlbNorm ()
{
  for (size_t i = 0; i < myVals.size(); i++)
    this->applyFinalOp(myVals[i]);
}


bool GlbNorm::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmNorm* ptr = dynamic_cast<const ElmNorm*>(elmObj);
  if (!ptr) return false;

  ElmNorm& elVals = *const_cast<ElmNorm*>(ptr);
  for (size_t i = 0; i < myVals.size(); i++)
  {
    myVals[i] += elVals[i];
    this->applyFinalOp(elVals[i]);
  }

  return true;
}


void GlbNorm::applyFinalOp (double& value) const
{
  switch (myOp)
    {
    case ABS:
      value = fabs(value);
      break;

    case SQRT:
      if (value < 0.0)
	value = -sqrt(-value);
      else
	value = sqrt(value);
      break;
    }
}
