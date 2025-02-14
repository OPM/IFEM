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
#include "tinyxml2.h"


GenAlphaSIM::GenAlphaSIM (SIMbase& s) : NewmarkSIM(s), prevSol(3), tempSol(3)
{
  // Default Newmark parameters (alpha = -0.1)
  alpha_m = 1.0;
  alpha_f = 0.9;
  beta = 0.3025;
  gamma = 0.6;
}


bool GenAlphaSIM::parse (const tinyxml2::XMLElement* elem)
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


void GenAlphaSIM::initSol (size_t nSol, size_t nDof)
{
  if (nSol < 3) nSol = 3;
  this->NewmarkSIM::initSol(nSol,nDof);

  nDof = solution.front().size();
  for (size_t i = 0; i < 3; i++)
  {
    prevSol[i].resize(nDof,true);
    tempSol[i].resize(nDof,true);
  }
}


bool GenAlphaSIM::advanceStep (TimeStep& param, bool updateTime)
{
  // Update the previous solution
  for (unsigned short int i = 0; i < 3 && i < solution.size(); i++)
    std::copy(solution[i].begin(),solution[i].end(),prevSol[i].begin());

  return this->NewmarkSIM::advanceStep(param,updateTime);
}


bool GenAlphaSIM::predictStep (TimeStep& param)
{
  if (!this->NewmarkSIM::predictStep(param))
    return false;

  for (unsigned short int i = 0; i < 3 && i < solution.size(); i++)
    std::copy(solution[i].begin(),solution[i].end(),tempSol[i].begin());

  // Compute the intermediate solution
  // {U}_(n+alpha) = (1-alpha)*{U}_n + alpha*{U}_(n+1)
  solution[0].relax(alpha_f,prevSol[0]);
  solution[1].relax(alpha_f,prevSol[1]);
  solution[2].relax(alpha_m,prevSol[2]);

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

  // Update current displacement, velocity and acceleration solutions
  double bdt = beta*param.time.dt;
  double bdt2 = bdt*param.time.dt;
  tempSol[0].add(linsol, solveDisp ? 1.0 : bdt2);
  tempSol[1].add(linsol, solveDisp ? gamma/bdt : gamma*param.time.dt);
  tempSol[2].add(linsol, solveDisp ? 1.0/bdt2 : 1.0);

  if (converged) // Converged solution
    for (unsigned short int i = 0; i < 3 && i < solution.size(); i++)
      std::copy(tempSol[i].begin(),tempSol[i].end(),solution[i].begin());
  else
  {
    // Compute new intermediate solution
    // {U}_(n+alpha) = (1-alpha)*{U}_n + alpha*{U}_(n+1)
    solution[0].relax(alpha_f,prevSol[0],tempSol[0]);
    solution[1].relax(alpha_f,prevSol[1],tempSol[1]);
    solution[2].relax(alpha_m,prevSol[2],tempSol[2]);
  }

#if SP_DEBUG > 1
  std::cout <<"\nCorrected displacement:"<< solution[0]
            <<"Corrected velocity:"<< solution[1]
            <<"Corrected acceleration:"<< solution[2];
#elif defined(SP_DEBUG)
  if (converged)
    std::cout <<"\nConverged displacement:"<< solution[0].max()
              <<"\nConverged velocity:"<< solution[1].max()
              <<"\nConverged acceleration:"<< solution[2].max() << std::endl;
#endif

  model.updateRotations(linsol,alpha_f); // TODO: Look at this
  return model.updateConfiguration(tempSol[0]);
}


void GenAlphaSIM::setSolution (const RealArray& newSol, int idx)
{
  if (idx == 0)
    tempSol[0] = newSol;

  solution[idx] = newSol;
}
