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
#include "SIMoutput.h"
#include "Profiler.h"


MultiStepSIM::MultiStepSIM (SIMbase& sim)
  : SIMinput(sim), model(static_cast<SIMoutput&>(sim))
{
#ifndef SP_DEBUG
  msgLevel = 1;   // prints the convergence history only
#elif SP_DEBUG > 2
  msgLevel = 100; // prints the linear solution vector if its size is < 100
#endif

  refNopt = MAX;
  refNorm = 1.0;

  geoBlk = nBlock = 0;
}


void MultiStepSIM::printProblem (std::ostream& os) const
{
  model.printProblem(os);
}


const char** MultiStepSIM::getPrioritizedTags () const
{
  return model.getPrioritizedTags();
}


bool MultiStepSIM::initEqSystem (bool withRF)
{
  return model.initSystem(opt.solver,1,1,withRF);
}


void MultiStepSIM::setStartGeo (int gID)
{
  model.setStartGeo(gID);
}


bool MultiStepSIM::saveModel (char* fileName)
{
  geoBlk = nBlock = 0; // initialize the VTF block counters
  return this->saveModel(geoBlk,nBlock,fileName);
}


bool MultiStepSIM::saveModel (int& gBlk, int& rBlk, char* fileName)
{
  PROFILE1("MultiStepSIM::saveModel");

  // Write VTF-file with model geometry
  if (!model.writeGlvG(gBlk,fileName,gBlk==0))
    return false;

  // Write Dirichlet boundary conditions
  return model.writeGlvBC(rBlk);
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
    if (!model.writeGlvT(iStep,geoBlk,nBlock))
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


/*!
  This method only writes the primary solution vector as a vector field.
  It is mainly used when this simulator is a component in a coupled simulation,
  and where the secondary solution is of minor interest.
*/

bool MultiStepSIM::saveStep (int iStep, int& nBlock, const char* vecName)
{
  return model.writeGlvV(solution.front(),vecName,iStep,nBlock);
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
