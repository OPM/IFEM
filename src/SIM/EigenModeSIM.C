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
  const char* value = NULL;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if ((value = utl::getValue(child,"mode")))
      if (utl::getAttribute(child,"number",imode) && imode > 0)
      {
        if (imode > amplitude.size())
          amplitude.resize(imode,0.0);
        amplitude[imode-1] = atof(value);
      }

  opt.nev = amplitude.size();
  if (opt.ncv < opt.nev*2)
    opt.ncv = opt.nev*2;

  return true;
}


void EigenModeSIM::printProblem (std::ostream& os) const
{
  model.printProblem(os);

  bool multiModes = false;
  os <<"EigenMode combination: u(x,t) =";
  for (size_t i = 0; i < amplitude.size(); i++)
    if (amplitude[i] != 0.0)
    {
      os << (multiModes ? "\n                              + ":" ")
         << amplitude[i] <<"*v"<< i+1 <<"(x)*sin(omega"<< i+1 <<"*t)";
      multiModes = true;
    }

  os << std::endl;
}


bool EigenModeSIM::initSol (size_t nSol)
{
  if (!this->MultiStepSIM::initSol(nSol))
    return false;
  else if (nSol > 0)
    return true;

  // Solve the eigenvalue problem giving the natural eigenfrequencies
  model.setMode(SIM::VIBRATION);
  model.initSystem(opt.solver,2,0,false);
  model.setQuadratureRule(opt.nGauss[0]);
  if (!model.assembleSystem())
    return false;

  if (!model.systemModes(modes))
    return false;

  // Scale all eigenvectors to have max amplitude equal to one,
  // and convert the eigenvalues back to angular frequency
  for (size_t i = 0; i < modes.size(); i++)
  {
    modes[i].eigVec /= modes[i].eigVec.normInf();
    modes[i].eigVal *= 2.0*M_PI;
  }

  return true;
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

  if (msgLevel >= 0 && myPid == 0) std::cout << std::endl;
  return this->solutionNorms(param.time,zero_tolerance,outPrec) ?
    SIM::CONVERGED : SIM::FAILURE;
}
