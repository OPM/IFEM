// $Id$
//==============================================================================
//!
//! \file EigenModeSIM.C
//!
//! \date Jan 8 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for computing a time history from a set of eigen modes.
//!
//==============================================================================

#include "EigenModeSIM.h"
#include "SIMoutput.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml.h"


EigenModeSIM::EigenModeSIM (SIMbase& sim) : MultiStepSIM(sim)
{
  rotUpd  = 't';
  opt.eig = 4;
  myStart = -1.0;
}


bool EigenModeSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"eigenmodes"))
    return model.parse(elem);

  size_t imode = 0;
  double freq = 0.0;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement()) {
    const char* value = nullptr;
    if ((value = utl::getValue(child,"mode")))
      if (utl::getAttribute(child,"number",imode) && imode > 0)
      {
        if (imode > amplitude.size())
          amplitude.resize(imode,0.0);
        amplitude[imode-1] = atof(value);
        if (utl::getAttribute(child,"frequency",freq) && freq > 0.0)
        {
          if (imode > omega.size())
            omega.resize(imode,0.0);
          omega[imode-1] = 2.0*M_PI*freq;
        }
      }
  }

  opt.nev = amplitude.size();
  if (opt.ncv < opt.nev*2)
    opt.ncv = opt.nev*2;

  return true;
}


void EigenModeSIM::printProblem () const
{
  model.printProblem();

  bool multiModes = false;
  IFEM::cout <<"EigenMode combination: u(x,t) =";
  for (size_t i = 0; i < amplitude.size(); i++)
    if (amplitude[i] != 0.0)
    {
      IFEM::cout << (multiModes ? "\n                              + ":" ")
                 << amplitude[i] <<"*v"<< i+1 <<"(x)*sin(";
      if (i < omega.size() && omega[i] > 0.0)
        IFEM::cout << omega[i];
      else
        IFEM::cout <<"omega"<< i+1;
      IFEM::cout <<"*t)";
      multiModes = true;
    }

  IFEM::cout << std::endl;
}


bool EigenModeSIM::initSol (size_t nSol)
{
  if (!this->MultiStepSIM::initSol(nSol))
    return false;
  else if (nSol > 0)
    return true;

  // Solve the eigenvalue problem giving the natural eigenfrequencies
  model.setMode(SIM::VIBRATION);
  model.initSystem(opt.solver,2,0);
  model.setQuadratureRule(opt.nGauss[0],true);
  if (!model.assembleSystem())
    return false;

  if (!model.systemModes(modes))
    return false;

  // Scale all eigenvectors to have max amplitude equal to one,
  // and convert the eigenvalues back to angular frequency
  for (size_t i = 0; i < modes.size(); i++)
  {
    modes[i].eigVec /= modes[i].eigVec.normInf();
    if (i < omega.size() && omega[i] > 0.0)
      modes[i].eigVal = omega[i];
    else
      modes[i].eigVal *= 2.0*M_PI;
#if SP_DEBUG > 1
    std::cout <<"\nEigenvector #"<< i+1 <<":"<< modes[i].eigVec;
#endif
  }

  return true;
}


bool EigenModeSIM::advanceStep (TimeStep& param, bool updateTime)
{
  this->pushSolution(); // Update solution vectors between time steps
  return this->MultiStepSIM::advanceStep(param,updateTime);
}


SIM::ConvStatus EigenModeSIM::solveStep (TimeStep& param, SIM::SolutionMode,
                                         double zero_tolerance,
                                         std::streamsize outPrec)
{
  if (myStart < 0.0) myStart = param.time.t - param.time.dt;
  double t = param.time.t - myStart; // Time since the start of this simulator

  solution.front().fill(0.0);
  for (size_t i = 0; i < amplitude.size() && i < modes.size(); i++)
    solution.front().add(modes[i].eigVec,amplitude[i]*sin(modes[i].eigVal*t));

  if (msgLevel >= 0) model.getProcessAdm().cout << std::endl;
  return this->solutionNorms(param.time,zero_tolerance,outPrec) ?
    SIM::CONVERGED : SIM::FAILURE;
}
