// $Id$
//==============================================================================
//!
//! \file GenAlphaSIM.C
//!
//! \date Aug 21 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Generalized-alpha driver for isogeometric dynamic FEM simulators.
//!
//==============================================================================

#include "GenAlphaSIM.h"
#include "SIMoutput.h"
#include "TimeStep.h"
#include "IFEM.h"
#include "tinyxml.h"


GenAlphaSIM::GenAlphaSIM (SIMbase& s) : NewmarkSIM(s), prevSol(3), tempSol(3)
{
  // Default Newmark parameters (alpha = -0.1)
  alpha_m = 1.0;
  alpha_f = 0.9;
  beta = 0.3025;
  gamma = 0.6;
}


bool GenAlphaSIM::parse (const TiXmlElement* elem)
{
  bool ok = this->NewmarkSIM::parse(elem);

  if (!strcasecmp(elem->Value(),inputContext))
  {
    double alpha = -0.1;
    const char* attr = elem->Attribute("alpha");
    if (attr)
    {
      alpha = atof(attr);
      alpha_f = alpha + 1.0;
      alpha_m = 1.0;
    }
    else if ((attr = elem->Attribute("rho_inf")))
    {
      double rho = atof(attr);
      alpha_f = 1.0/(1.0+rho);
      alpha_m = (2.0-rho)/(1.0+rho);
      alpha = alpha_f - alpha_m;
    }
    beta = 0.25*(1.0-alpha)*(1.0-alpha);
    gamma = 0.5 - alpha;
  }

  return ok;
}


void GenAlphaSIM::initPrm ()
{
  model.setIntegrationPrm(0,alpha1);
  model.setIntegrationPrm(1,fabs(alpha2));
  model.setIntegrationPrm(2,solveDisp ? -alpha_m : alpha_m);
  model.setIntegrationPrm(3,alpha_f);
  model.setIntegrationPrm(4,2.0);
}


bool GenAlphaSIM::initSol (size_t nSol)
{
  this->MultiStepSIM::initSol(nSol);

  size_t nDOFs = model.getNoDOFs();
  for (size_t i = 0; i < 3; i++)
  {
    prevSol[i].resize(nDOFs,true);
    tempSol[i].resize(nDOFs,true);
  }

  return true;
}


bool GenAlphaSIM::advanceStep (TimeStep& param, bool updateTime)
{
  // Update displacement solutions between time steps
  for (int n = solution.size()-3; n > 0; n--)
    std::copy(solution[n-1].begin(),solution[n-1].end(),solution[n].begin());

  // Update the previous solution
  std::copy(solution.front().begin(),solution.front().end(),prevSol[0].begin());
  if (solution.size() > 2)
  {
    size_t iA = solution.size() - 1;
    size_t iV = solution.size() - 2;
    std::copy(solution[iV].begin(),solution[iV].end(),prevSol[1].begin());
    std::copy(solution[iA].begin(),solution[iA].end(),prevSol[2].begin());
  }

  return this->NewmarkSIM::advanceStep(param,updateTime);
}


bool GenAlphaSIM::predictStep (TimeStep& param)
{
  if (!this->NewmarkSIM::predictStep(param))
    return false;

  size_t ipD = 0;
  size_t ipA = solution.size() - 1;
  size_t ipV = solution.size() - 2;
  tempSol[0] = solution[ipD];
  tempSol[1] = solution[ipV];
  tempSol[2] = solution[ipA];

  // Compute the intermediate solution
  // {U}_(n+alpha) = (1-alpha)*{U}_n + alpha*{U}_(n+1)
  solution[ipD].relax(alpha_f,prevSol[0]);
  solution[ipV].relax(alpha_f,prevSol[1]);
  solution[ipA].relax(alpha_m,prevSol[2]);
  return true;
}


bool GenAlphaSIM::correctStep (TimeStep& param, bool converged)
{
  if (solution.size() < 3)
  {
    std::cerr <<" *** GenAlphaSIM::correctStep: Too few solution vectors "
              << solution.size() << std::endl;
    return false;
  }

#ifdef SP_DEBUG
  std::cout <<"\nGenAlphaSIM::correctStep(converged="
            << std::boolalpha << converged <<")";
#endif

  size_t ipD = 0;
  size_t ipA = solution.size() - 1;
  size_t ipV = solution.size() - 2;

  // Update current displacement, velocity and acceleration solutions
  double bdt = beta*param.time.dt;
  double bdt2 = bdt*param.time.dt;
  tempSol[0].add(linsol, solveDisp ? 1.0 : bdt2);
  tempSol[1].add(linsol, solveDisp ? gamma/bdt : gamma*param.time.dt);
  tempSol[2].add(linsol, solveDisp ? 1.0/bdt2 : 1.0);

  if (converged)
  {
    // Converged solution
    solution[ipD] = tempSol[0];
    solution[ipV] = tempSol[1];
    solution[ipA] = tempSol[2];
  }
  else
  {
    // Compute new intermediate solution
    // {U}_(n+alpha) = (1-alpha)*{U}_n + alpha*{U}_(n+1)
    solution[ipD].relax(alpha_f,prevSol[0],tempSol[0]);
    solution[ipV].relax(alpha_f,prevSol[1],tempSol[1]);
    solution[ipA].relax(alpha_m,prevSol[2],tempSol[2]);
  }

#if SP_DEBUG > 1
  std::cout <<"\nCorrected displacement:"<< solution[ipD]
            <<"Corrected velocity:"<< solution[ipV]
            <<"Corrected acceleration:"<< solution[ipA];
#elif defined(SP_DEBUG)
  if (converged)
    std::cout <<"\nConverged displacement:"<< solution[ipD].max()
              <<"\nConverged velocity:"<< solution[ipV].max()
              <<"\nConverged acceleration:"<< solution[ipA].max() << std::endl;
#endif

  model.updateRotations(linsol,alpha_f); //TODO look at this
  return model.updateConfiguration(tempSol[0]);
}


void GenAlphaSIM::setSolution (const Vector& newSol, int idx)
{
  solution[idx] = newSol;

  if (idx == 0)
    tempSol[0] = newSol;
  else if (solution.size() > 2)
    tempSol[solution.size()-3+idx] = newSol;
}
