// $Id$
//==============================================================================
//!
//! \file NewmarkNLSIM.C
//!
//! \date Jul 4 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Newmark solution driver for isogeometric dynamic FEM simulators.
//!
//==============================================================================

#include "NewmarkNLSIM.h"
#include "SystemMatrix.h"
#include "SIMbase.h"
#include "TimeStep.h"
#include "Utilities.h"
#include "tinyxml.h"


//! \brief Convenience enum defining total solution vector indices.
enum { iD=0, iV=2, iA=3, nSOL=4 };


NewmarkNLSIM::NewmarkNLSIM (SIMbase& sim) : NewmarkSIM(sim), Finert(NULL)
{
  // Default Newmark parameters (alpha = -0.1)
  beta = 0.3025;
  gamma = 0.6;
}


bool NewmarkNLSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"newmarksolver"))
    return model.parse(elem);

  double alpha = -0.1;
  utl::getAttribute(elem,"alpha",alpha);
  beta = 0.25*(1.0-alpha)*(1.0-alpha);
  gamma = 0.5 - alpha;

  utl::getAttribute(elem,"alpha1",alpha1);
  utl::getAttribute(elem,"alpha2",alpha2);

  const char* value = 0;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if ((value = utl::getValue(child,"maxits")))
      maxit = atoi(value);
    else if ((value = utl::getValue(child,"rtol")))
      convTol = atof(value);
    else if ((value = utl::getValue(child,"dtol")))
      divgLim = atof(value);
    else if ((value = utl::getValue(child,"referencenorm")))
    {
      if (!strcasecmp(value,"all"))
        refNopt = ALL;
      else if (!strcasecmp(value,"max"))
        refNopt = MAX;
    }

  return true;
}


void NewmarkNLSIM::init (size_t)
{
  model.setIntegrationPrm(0,alpha1);
  model.setIntegrationPrm(1,alpha2);
  model.setIntegrationPrm(2,0.5-gamma);

  solution.resize(nSOL);
  for (Vectors::iterator it = solution.begin(); it != solution.end(); ++it)
    it->resize(model.getNoDOFs(),true);
}


bool NewmarkNLSIM::initEqSystem (bool withRF)
{
  return model.initSystem(opt.solver,1,2,withRF);
}


bool NewmarkNLSIM::advanceStep (TimeStep& param, bool updateTime)
{
  // Update displacement solutions between time steps
  for (int n = solution.size()-3; n > 0; n--)
    std::copy(solution[n-1].begin(),solution[n-1].end(),solution[n].begin());

  return updateTime ? param.increment() : true;
}


void NewmarkNLSIM::finalizeRHSvector ()
{
  if (Finert)
    model.addToRHSvector(0,*Finert,gamma-0.5);
}


bool NewmarkNLSIM::predictStep (TimeStep& param)
{
  if (solution.size() < nSOL) return false;

  Vector v(solution[iV]);
  Vector a(solution[iA]);

  // Predicted new velocity
  solution[iV] *= gamma/beta - 1.0;
  solution[iV].add(a,(0.5*gamma/beta-1.0)*param.time.dt);

  // Predicted new acceleration
  solution[iA] *= 0.5/beta - 1.0;
  solution[iA].add(v,1.0/(beta*param.time.dt));

  return true;
}


bool NewmarkNLSIM::correctStep (TimeStep& param, bool converged)
{
  if (solution.size() < nSOL) return false;

  // Update current displacement, velocity and acceleration solutions
  solution[iD].add(linsol,1.0);
  solution[iV].add(linsol,gamma/(beta*param.time.dt));
  solution[iA].add(linsol,1.0/(beta*param.time.dt*param.time.dt));

  if (converged)
  {
    // Save the actual inertia vector (minus the residual) from converged step
    delete Finert;
    Finert = model.getRHSvector(1,true);
  }

  return model.updateConfiguration(solution[iD]);
}
