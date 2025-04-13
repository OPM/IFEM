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
#include "Property.h"
#include "ElmMats.h"
#include "SAM.h"
#include "ProcessAdm.h"


ReactionsOnly::ReactionsOnly (const SAM* sam, const ProcessAdm& adm,
                              RealArray* rf, Vector* sf)
  : mySam(sam), myAdm(adm), R(rf), S(sf)
{
  b.redim(sam->getNoEquations());
  if (R)
    R->resize(sam->getNoSpecifiedDOFs());
}


void ReactionsOnly::initialize (char)
{
  b.init();
  if (R)
    std::fill(R->begin(),R->end(),0.0);
}


bool ReactionsOnly::finalize (bool)
{
  if (R)
  {
    if (myAdm.dd.isPartitioned())
      myAdm.allReduceAsSum(*R);

    for (double& r : *R) r *= -1.0;
#if SP_DEBUG > 2
    std::cout <<"\nReaction forces:"<< Vector(*R);
#endif
  }

  if (!S)
    return true;
  else if (!b.beginAssembly() || !b.endAssembly())
    return false;
  else if (!mySam->expandSolution(b,*S,0.0))
    return false;

#if SP_DEBUG > 2
  mySam->printVector(std::cout,*S,"Internal forces");
#endif
  return true;
}


bool ReactionsOnly::assemble (const LocalIntegral* elmObj, int elmId)
{
  const ElmMats* elMat = dynamic_cast<const ElmMats*>(elmObj);
  if (elMat && mySam->assembleSystem(b,elMat->getRHSVector(),elmId,R))
    return true;

  std::cerr <<" *** ReactionsOnly::assemble: Failure for element "<< elmId
            << std::endl;
  return false;
}


bool ReactionsOnly::haveContributions (size_t pidx,
                                       const PropertyVec& pvec) const
{
  // Lambda function checking if a property have force contribution on a patch.
  auto&& contributed = [this,pidx](const Property& p)
  {
    if (p.patch != pidx)
      return false;
    else if (S)
      return p.pcode == Property::OTHER;
    else if (!R)
      return false;

    return (p.pcode >= Property::DIRICHLET &&
            p.pcode <= Property::DIRICHLET_ANASOL);
  };

  return std::find_if(pvec.begin(),pvec.end(),contributed) != pvec.end();
}
