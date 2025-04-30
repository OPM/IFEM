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
#include "tinyxml2.h"


EigenModeSIM::EigenModeSIM (SIMbase& sim) : MultiStepSIM(sim)
{
  rotUpd  = 't';
  opt.eig = 4;
  myStart = -1.0;
}


bool EigenModeSIM::parse (const tinyxml2::XMLElement* elem)
{
  if (strcasecmp(elem->Value(),"eigenmodes"))
    return model.parse(elem);

  const tinyxml2::XMLElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
  {
    size_t imode = 0;
    double freq = 0.0;
    utl::getAttribute(child,"frequency",freq);
    const char* value = utl::getValue(child,"mode");
    if (value && utl::getAttribute(child,"number",imode) && imode > 0)
    {
      if (imode > amplitude.size())
        amplitude.resize(imode,0.0);
      amplitude[imode-1] = atof(value);

      if (freq > epsZ)
      {
        if (imode > omega.size())
          omega.resize(imode,0.0);
        omega[imode-1] = 2.0*M_PI*freq;
      }
    }
    else if ((value = utl::getValue(child,"modes")))
    {
      double scale = 1.0;
      std::vector<int> modes;
      utl::parseIntegers(modes,value);
      utl::getAttribute(child,"scale",scale);
      std::sort(modes.begin(),modes.end());
      imode = modes.back();
      if (imode > amplitude.size())
        amplitude.resize(imode,0.0);
      if (freq > epsZ && imode > omega.size())
        omega.resize(imode,0.0);
      for (int mode : modes)
      {
        amplitude[mode-1] = scale;
        if (freq > epsZ)
          omega[imode-1] = 2.0*M_PI*freq;
      }
    }
  }

  if (opt.nev < static_cast<int>(amplitude.size()))
    opt.nev = amplitude.size();
  int mincv = opt.nev < 10 ? opt.nev*2 : opt.nev+10;
  if (opt.ncv < mincv)
    opt.ncv = mincv;

  return true;
}


void EigenModeSIM::printProblem (bool stopInputTimer) const
{
  this->MultiStepSIM::printProblem(stopInputTimer);

  bool multiModes = false;
  IFEM::cout <<"Eigenmode combination: u(x,t) =";
  for (size_t i = 0; i < amplitude.size(); i++)
    if (amplitude[i] > epsZ)
    {
      IFEM::cout << (multiModes ? "\n                              + ":" ");
      if (fabs(amplitude[i]-1.0) > epsZ)
        IFEM::cout << amplitude[i] <<"*";
      IFEM::cout <<"v"<< i+1 <<"(x)*sin(";
      if (i < omega.size() && omega[i] > epsZ)
        IFEM::cout << omega[i];
      else
        IFEM::cout <<"omega"<< i+1;
      IFEM::cout <<"*t)";
      multiModes = true;
    }

  IFEM::cout << std::endl;
}


void EigenModeSIM::initSol (size_t nSol, size_t nDof)
{
  this->MultiStepSIM::initSol(nSol,nDof);
  if (nSol > 0) return;

  // Solve the eigenvalue problem giving the natural eigenfrequencies
  model.setMode(SIM::VIBRATION);
  model.initSystem(opt.solver,2,0);
  model.setQuadratureRule(opt.nGauss[0],true);
  if (!model.assembleSystem() || !model.systemModes(modes))
  {
    solution.clear();
    std::cerr <<" *** EigenModeSIM::initSol: No eigenvalue solution."
              << std::endl;
    return;
  }

  // Scale all eigenvectors to have max amplitude equal to one,
  // and convert the eigenvalues back to angular frequency
  for (size_t i = 0; i < modes.size(); i++)
  {
    modes[i].eigVec /= modes[i].eigVec.normInf();
    if (i < omega.size() && omega[i] > epsZ)
      modes[i].eigVal = omega[i];
    else
      modes[i].eigVal *= 2.0*M_PI;
#if SP_DEBUG > 1
    std::cout <<"\nEigenvector #"<< i+1 <<":"<< modes[i].eigVec;
#endif
  }
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
  if (solution.empty())
    return SIM::FAILURE;

  if (msgLevel >= 0)
    this->printStep(param);

  if (myStart < 0.0)
    myStart = param.time.t - param.time.dt;

  double t = param.time.t - myStart; // Time since the start of this simulator

  solution.front().fill(0.0);
  for (size_t i = 0; i < amplitude.size() && i < modes.size(); i++)
    if (amplitude[i] > epsZ)
      solution.front().add(modes[i].eigVec,amplitude[i]*sin(modes[i].eigVal*t));

  return this->solutionNorms(param.time,zero_tolerance,outPrec) ?
    SIM::CONVERGED : SIM::FAILURE;
}
