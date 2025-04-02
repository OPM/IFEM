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
#include "HDF5Restart.h"
#include "SIMoutput.h"
#include "TimeStep.h"
#include "Profiler.h"
#include "IFEM.h"


MultiStepSIM::MultiStepSIM (SIMbase& sim)
  : SIMadmin(sim), model(static_cast<SIMoutput&>(sim))
{
#ifndef SP_DEBUG
  msgLevel = 1;   // prints the convergence history only
#elif SP_DEBUG > 2
  msgLevel = 100; // prints the linear solution vector if its size is < 100
#endif

  refNopt = MAX;
  refNorm = 0.0;
  rCond   = -1.0;
  subiter = NONE;
  nRHSvec = 1;
  rotUpd  = false;

  geoBlk = nBlock = lastSt = 0;
}


void MultiStepSIM::printProblem () const
{
  model.printProblem();
}


const char** MultiStepSIM::getPrioritizedTags () const
{
  return model.getPrioritizedTags();
}


bool MultiStepSIM::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  return model.preprocess(ignored,fixDup);
}


void MultiStepSIM::initSol (size_t nSol, size_t nDof)
{
  if (!solution.empty())
    nSol = solution.size();
  else if (nSol > 0 && nSol < model.getNoSolutions(false))
    nSol = model.getNoSolutions(false);

  if (nDof == 0)
    nDof = model.getNoDOFs();

  this->initSolution(nDof,nSol);
}


bool MultiStepSIM::initEqSystem (bool withRF, size_t nScl)
{
  return model.initSystem(opt.solver,1,nRHSvec,nScl,withRF);
}


void MultiStepSIM::setStartGeo (int gID)
{
  model.setStartGeo(gID);
}


bool MultiStepSIM::saveModel (const char* fileName)
{
  geoBlk = nBlock = 0; // initialize the VTF block counters
  return this->saveModel(geoBlk,nBlock,fileName);
}


bool MultiStepSIM::saveModel (int& gBlock, int& rBlock,
                              const char* fileName, bool clearG)
{
  PROFILE1("MultiStepSIM::saveModel");

  // Reset nViz to 1 if less than 3 parametric dimensions
  for (unsigned short int i = model.getNoParamDim(); i < 3; i++)
    model.opt.nViz[i] = 1;

  // Write VTF-file with model geometry
  if (!model.writeGlvG(gBlock,fileName,clearG))
    return false;

  // Write Dirichlet boundary conditions
  return model.writeGlvBC(rBlock);
}


bool MultiStepSIM::saveStep (int iStep, double time, const char* vecName)
{
  if (abs(iStep) <= lastSt) return true; // We have already saved this step

  PROFILE1("MultiStepSIM::saveStep");

  // Negative iStep means we are saving the initial state only
  if (!model.setMode(iStep < 0 ? SIM::INIT : SIM::RECOVERY))
    return false;
  else if (iStep < 0)
    iStep = -iStep;
  lastSt = iStep;

  // Write boundary tractions, if any
  if (!opt.pSolOnly || opt.saveTrac)
    if (!model.writeGlvT(iStep,geoBlk,nBlock))
      return false;

  // Write load and residual force vectors.
  // Only for Lagrange, or for Splines with no extra visualization points.
  if ((opt.discretization < ASM::Spline) ||
      (opt.nViz[0] == 2 && opt.nViz[1] <= 2 && opt.nViz[2] <= 2))
  {
    if (loadVec.size() == model.getNoDOFs())
      if (!model.writeGlvV(loadVec,"Load vector",iStep,nBlock,2))
        return false;

    if (!opt.pSolOnly && residual.size() == model.getNoDOFs())
      if (!model.writeGlvV(residual,"Residual forces",iStep,nBlock,3))
        return false;
  }

  // Write solution fields
  if (this->numSolution())
    if (!model.writeGlvS(this->realSolution(),iStep,nBlock,time,vecName))
      return false;

  // Write any problem-specific data (rigid body transformations, etc.)
  if (!model.writeGlvA(nBlock,iStep))
    return false;

  // Write time step information
  return model.writeGlvStep(iStep,time);
}


/*!
  Use this method when other simulators write results to the same VTF-file.
  The internal result block counter \a nBlock is syncronized with the
  argument \a rBlock before the results of this simulator are written,
  to avoid that multiple result blocks recieve the same result block ID.
*/

bool MultiStepSIM::saveStep (int iStep, int& rBlock, double time,
                             const char* vecName)
{
  if (rBlock > nBlock) nBlock = rBlock;
  bool s = this->saveStep(iStep,time,vecName);
  rBlock = nBlock;
  return s;
}


/*!
  This method only writes the primary solution vector as a vector field.
  It is mainly used when this simulator is a component in a coupled simulation,
  and where the secondary solution is of minor interest.
*/

bool MultiStepSIM::saveStep (int iStep, int& rBlock, const char* vecName)
{
  if (!this->numSolution())
    return true;

  return model.writeGlvV(this->realSolution(),vecName,iStep,rBlock);
}


bool MultiStepSIM::serialize (SerializeMap& data) const
{
  data["MultiStepSIM::refNorm"] = SIMsolution::serialize(&refNorm,1);

  return this->saveSolution(data,model.getName());
}


bool MultiStepSIM::deSerialize (const SerializeMap& data)
{
  SerializeMap::const_iterator sit = data.find("MultiStepSIM::refNorm");
  if (sit != data.end())
    SIMsolution::deSerialize(sit->second,&refNorm,1);

  return this->restoreSolution(data,model.getName());
}


void MultiStepSIM::dumpStep (int iStep, double time, utl::LogStream& os,
                             bool withID) const
{
  if (!this->numSolution())
    return;

  if (withID)
    os <<"\n\n# Dump of primary solution at Step "<< iStep
       <<", Time="<< time <<"\n";

  model.dumpPrimSol(this->realSolution(),os,withID);
}


void MultiStepSIM::dumpResults (double time, utl::LogStream& os,
                                std::streamsize precision, bool formatted) const
{
  if (this->numSolution())
    model.dumpResults(this->realSolution(),time,os,formatted,precision);

  model.dumpMoreResults(time,os,precision);
}


bool MultiStepSIM::hasPointResultFile () const
{
  return model.hasPointResultFile();
}


bool MultiStepSIM::savePoints (double time, int step) const
{
  if (!this->numSolution())
    return true;

  return model.savePoints(this->realSolution(),time,step);
}


void MultiStepSIM::printStep (const TimeStep& param) const
{
  double digits = log10(param.time.t) - log10(param.time.dt);
  if (digits > 6.0)
  {
    // Increase precision on physical time for larger time series
    utl::LogStream& cout = model.getProcessAdm().cout;
    std::streamsize prec = cout.precision(ceil(digits));
    model.printStep(param.step,param.time);
    cout.precision(prec);
  }
  else
    model.printStep(param.step,param.time);
}


bool MultiStepSIM::advanceStep (TimeStep& param, bool updateTime)
{
  if (param.step > 0 && rotUpd == 't')
    // Update nodal rotations of previous time step
    model.updateRotations(Vector());

  return updateTime ? param.increment() : true;
}


bool MultiStepSIM::solutionNorms (const TimeDomain&, double zero_tolerance,
                                  std::streamsize outPrec)
{
  if (msgLevel < 0 || !this->numSolution())
    return true;

  double old_zero_tol = utl::zero_print_tol;
  utl::zero_print_tol = zero_tolerance;
  model.printSolutionSummary(this->realSolution(),0,nullptr,outPrec);
  utl::zero_print_tol = old_zero_tol;

  return true;
}


bool MultiStepSIM::checkForRestart ()
{
  if (opt.restartFile.empty())
    return true; // No restart

  HDF5Restart::SerializeData data;
  HDF5Restart hdf(opt.restartFile,adm);
  int restartStep = hdf.readData(data,opt.restartStep-1);
  if (restartStep < 0 || !this->deSerialize(data))
  {
    std::cerr <<" *** Failed to read restart data."<< std::endl;
    return false;
  }

  IFEM::cout <<"\n === Restarting from a serialized state ==="
             <<"\n     file = "<< opt.restartFile
             <<"\n     step = "<< restartStep+1 << std::endl;
  return true;
}
