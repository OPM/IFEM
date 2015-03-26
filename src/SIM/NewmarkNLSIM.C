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
#include "SIMoutput.h"
#include "TimeStep.h"
#include "tinyxml.h"


NewmarkNLSIM::NewmarkNLSIM (SIMbase& sim) : NewmarkSIM(sim), Finert(NULL)
{
  // Default Newmark parameters (alpha = -0.1)
  beta = 0.3025;
  gamma = 0.6;

  predictor = 'd'; // default predictor (constant displacement)
}


bool NewmarkNLSIM::parse (const TiXmlElement* elem)
{
  bool ok = this->NewmarkSIM::parse(elem);

  if (!strcasecmp(elem->Value(),"newmarksolver"))
  {
    double alpha = -0.1;
    const char* attr = elem->Attribute("alpha");
    if (attr)
      alpha = atof(attr);
    else if ((attr = elem->Attribute("rho_inf")))
    {
      double rho = atof(attr);
      double alpha_f = 1.0/(1.0+rho);
      double alpha_m = (2.0-rho)/(1.0+rho);
      alpha = alpha_f - alpha_m;
    }
    beta = 0.25*(1.0-alpha)*(1.0-alpha);
    gamma = 0.5 - alpha;

    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"trueinertia"))
        nRHSvec = 2;
  }

  return ok;
}


void NewmarkNLSIM::printProblem (utl::LogStream& os) const
{
  this->NewmarkSIM::printProblem(os);

  if (alpha2 > 0.0)
    os <<"- based on the tangential stiffness matrix\n";
  else if (alpha2 < 0.0)
    os <<"- based on the material stiffness matrix only\n";

  if (nRHSvec > 1)
    os <<"- including true inertia forces from previous step in residual\n";
  else if (alpha2 == 0.0)
    return;

  os << std::endl;
}


void NewmarkNLSIM::initPrm ()
{
  model.setIntegrationPrm(0,alpha1);
  model.setIntegrationPrm(1,fabs(alpha2));
  model.setIntegrationPrm(2,0.5-gamma);
  if (alpha2 < 0.0) // Flag that stiffness-proportional damping should depend
    model.setIntegrationPrm(3,-1.0); // on the material stiffness matrix only
}


bool NewmarkNLSIM::initSol (size_t nSol)
{
  size_t nDOFs = model.getNoDOFs();
  incDis.resize(nDOFs,true);
  predVel.resize(nDOFs,true);
  predAcc.resize(nDOFs,true);

  return this->MultiStepSIM::initSol(nSol);
}


bool NewmarkNLSIM::advanceStep (TimeStep& param, bool updateTime)
{
  // Update displacement solutions between time steps
  for (int n = solution.size()-3; n > 0; n--)
    std::copy(solution[n-1].begin(),solution[n-1].end(),solution[n].begin());

  return this->NewmarkSIM::advanceStep(param,updateTime);
}


void NewmarkNLSIM::finalizeRHSvector (bool)
{
  if (Finert)
    model.addToRHSvector(0,*Finert,gamma-0.5);
}


bool NewmarkNLSIM::predictStep (TimeStep& param)
{
  if (solution.size() < 3)
  {
    std::cerr <<" *** NewmarkNLSIM::predictStep: Too few solution vectors "
              << solution.size() << std::endl;
    return false;
  }

  if (rotUpd && model.getNoFields(1) == 6)
  {
    // Initalize the total angular rotations for this time step
    Vector& psol = solution.front();
    if (psol.size() != 6*model.getNoNodes(true))
    {
      std::cerr <<" *** NewmarkNLSIM::predictStep: Invalid dimension on"
                <<" the displacement vector "<< psol.size()
                <<" != "<< 6*model.getNoNodes(true) << std::endl;
      return false;
    }
    for (size_t i = 3; i < psol.size(); i += 6)
      psol[i] = psol[i+1] = psol[i+2] = 0.0;
  }

#ifdef SP_DEBUG
  std::cout <<"\nNewmarkNLSIM::predictStep";
#endif

  size_t iA = solution.size() - 1;
  size_t iV = solution.size() - 2;

  // Predicted new velocity
  predVel = solution[iV];
  predVel *= gamma/beta - 1.0;
  predVel.add(solution[iA],(0.5*gamma/beta-1.0)*param.time.dt);

  // Predicted new acceleration
  predAcc = solution[iA];
  predAcc *= 0.5/beta - 1.0;
  predAcc.add(solution[iV],1.0/(beta*param.time.dt));

#ifdef SP_DEBUG
  std::cout <<"\nPredicted velocity:";
#if SP_DEBUG > 1
  std::cout << predVel;
#else
  std::cout << predVel.max() <<"\n";
#endif
  std::cout <<"Predicted acceleration:";
#if SP_DEBUG > 1
  std::cout << predAcc;
#else
  std::cout << predAcc.max() << std::endl;
#endif
#endif

  solution[iV] = predVel;
  solution[iA] = predAcc;

  incDis.fill(0.0);
  predVel *= -1.0;
  predAcc *= -1.0;

  return true;
}


bool NewmarkNLSIM::correctStep (TimeStep& param, bool converged)
{
#ifdef SP_DEBUG
  std::cout <<"\nNewmarkNLSIM::correctStep(converged="
            << std::boolalpha << converged <<")";
#endif

  size_t iD = 0;
  size_t iA = solution.size() - 1;
  size_t iV = solution.size() - 2;

  // Update current displacement, velocity and acceleration solutions
  incDis.add(linsol,1.0);
  solution[iD].add(linsol,1.0);
  solution[iV] = predVel;
  solution[iV].add(incDis,gamma/(beta*param.time.dt));
  solution[iA] = predAcc;
  solution[iA].add(incDis,1.0/(beta*param.time.dt*param.time.dt));

  if (converged)
  {
    // Save the actual inertia vector (minus the residual) from converged step
    delete Finert;
    Finert = model.getRHSvector(1,true);
  }

#if SP_DEBUG > 1
  std::cout <<"\nCorrected displacement:"<< solution[iD]
            <<"Corrected velocity:"<< solution[iV]
            <<"Corrected acceleration:"<< solution[iA];
  if (converged && Finert)
    std::cout <<"Actual inertia force:"<< *Finert;
#elif defined(SP_DEBUG)
  if (converged)
    std::cout <<"\nConverged displacement:"<< solution[iD].max()
              <<"\nConverged velocity:"<< solution[iV].max()
              <<"\nConverged acceleration:"<< solution[iA].max() << std::endl;
#endif

  if (rotUpd == 't')
    model.updateRotations(solution[iD]);
  else if (rotUpd)
    model.updateRotations(linsol,1.0);

  return model.updateConfiguration(solution[iD]);
}


void NewmarkNLSIM::setSolution (const Vector& newSol, int idx)
{
  if (idx == 0)
  {
    // When updating the displacements within sub-iterations, we must also
    // update the incremental displacement accordingly (well spotted, Knut)
    incDis.add(solution.front(),-1.0);
    incDis.add(newSol,1.0);
  }

  solution[idx] = newSol;
}
