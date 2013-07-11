// $Id$
//==============================================================================
//!
//! \file MultiStepSIM.C
//!
//! \date Jul 11 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for multi-step solution drivers.
//!
//==============================================================================

#include "MultiStepSIM.h"
#include "SIMbase.h"
#include "SIMenums.h"
#include "Profiler.h"


MultiStepSIM::MultiStepSIM (SIMbase& sim) : SIMinput(sim), model(sim)
{
#ifndef SP_DEBUG
  msgLevel = 1;   // prints the convergence history only
#elif SP_DEBUG > 2
  msgLevel = 100; // prints the linear solution vector if its size is < 100
#endif

  geoBlk = nBlock = 0;
}


const char** MultiStepSIM::getPrioritizedTags () const
{
  return model.getPrioritizedTags();
}


bool MultiStepSIM::saveModel (char* fileName)
{
  PROFILE1("MultiStepSIM::saveModel");

  geoBlk = nBlock = 0; // initialize the VTF block counters

  // Write VTF-file with model geometry
  if (!model.writeGlvG(geoBlk,fileName))
    return false;

  // Write Dirichlet boundary conditions
  return model.writeGlvBC(nBlock);
}


bool MultiStepSIM::saveStep (int iStep, double time,
                             bool psolOnly, const char* vecName)
{
  PROFILE1("MultiStepSIM::saveStep");

  // Negative iStep means we are saving the initial state only
  if (!model.setMode(iStep < 0 ? SIM::INIT : SIM::RECOVERY))
    return false;
  else if (iStep < 0)
    iStep = -iStep;

  // Write boundary tractions, if any
  if (!psolOnly)
    if (!model.writeGlvT(iStep,nBlock))
      return false;

  // Write residual force vector, but only when no extra visualization points
  if (!psolOnly && opt.nViz[0] == 2 && opt.nViz[1] <= 2 && opt.nViz[2] <= 2)
    if (!model.writeGlvV(residual,"Residual forces",iStep,nBlock))
      return false;

  // Write solution fields
  if (!solution.empty())
    if (!model.writeGlvS(solution.front(),iStep,nBlock,time,psolOnly,vecName))
      return false;

  // Write any problem-specific data (rigid body transformations, etc.)
  if (!model.writeGlvA(nBlock,iStep))
    return false;

  // Write time step information
  return model.writeGlvStep(iStep,time);
}


void MultiStepSIM::dumpStep (int iStep, double time, std::ostream& os,
                             bool withID) const
{
  if (withID)
    os <<"\n\n# Dump of primary solution at Step "<< iStep
       <<", Time="<< time <<"\n";

  model.dumpPrimSol(solution.front(),os,withID);
}


void MultiStepSIM::dumpResults (double time, std::ostream& os,
                                std::streamsize precision, bool formatted) const
{
  model.dumpResults(solution.front(),time,os,formatted,precision);
  model.dumpMoreResults(time,os,precision);
}
