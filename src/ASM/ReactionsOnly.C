// $Id$
//==============================================================================
//!
//! \file ReactionsOnly.C
//!
//! \date Nov 7 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Global integral sub-class for reaction force integration.
//!
//==============================================================================

#include "ReactionsOnly.h"
#include "ElmMats.h"
#include "SAM.h"


ReactionsOnly::ReactionsOnly (Vector& rf, const SAM* sam) : mySam(sam), R(rf)
{
  mySam->initForAssembly(b,&R);
}


void ReactionsOnly::initialize (bool)
{
  b.init();
  R.fill(0.0);
}


bool ReactionsOnly::finalize (bool)
{
  R *= -1.0;
#if SP_DEBUG > 2
  std::cout <<"\nReaction forces:"<< R;
#endif
  return true;
}


bool ReactionsOnly::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmMats* elMat = dynamic_cast<const ElmMats*>(elmObj);
  if (elMat && mySam->assembleSystem(b,elMat->getRHSVector(),elmId,&R))
    return true;

  std::cerr <<" *** ReactionsOnly::assemble: Failure for element "<< elmId
            << std::endl;
  return false;
}
