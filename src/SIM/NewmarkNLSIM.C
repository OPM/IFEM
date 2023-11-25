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
#include "SIMoutput.h"
#include "SystemMatrix.h"
#include "TimeStep.h"
#include "IFEM.h"
#include "tinyxml2.h"


NewmarkNLSIM::NewmarkNLSIM (SIMbase& sim) : NewmarkSIM(sim), Finert(nullptr)
{
  // Default Newmark parameters (alpha = -0.1)
  beta = 0.3025;
  gamma = 0.6;

  predictor = 'd'; // default predictor (constant displacement)

  iD = iV = iA = 0;
}


NewmarkNLSIM::~NewmarkNLSIM ()
{
  delete Finert;
}


bool NewmarkNLSIM::parse (const tinyxml2::XMLElement* elem)
{
  bool ok = this->NewmarkSIM::parse(elem);

  if (!strcasecmp(elem->Value(),inputContext))
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

    const tinyxml2::XMLElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      if (!strcasecmp(child->Value(),"trueinertia"))
        nRHSvec = 2; // using true inertia forces from previous time step
  }

  return ok;
}


void NewmarkNLSIM::printProblem () const
{
  this->NewmarkSIM::printProblem();

  if (alpha2 > 0.0)
    IFEM::cout <<"- based on the tangential stiffness matrix\n";
  else if (alpha2 < 0.0)
    IFEM::cout <<"- based on the material stiffness matrix only\n";

  if (nRHSvec > 1)
    IFEM::cout <<"- including true inertia forces from previous step"
               <<" in the force residual\n";
  else if (alpha2 == 0.0)
    return;

  IFEM::cout << std::endl;
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
  this->MultiStepSIM::initSol(nSol);
  if (solution.size() < 3)
  {
    std::cerr <<" *** NewmarkNLSIM::initSol: Too few solution vectors "
              << solution.size() << std::endl;
    return false;
  }

  iA = solution.size() - 1;
  iV = solution.size() - 2;

  return true;
}


bool NewmarkNLSIM::advanceStep (TimeStep& param, bool updateTime)
{
  // Update displacement solutions between time steps
  if (solution.size() > 3)
    this->pushSolution(solution.size()-2);

  return this->NewmarkSIM::advanceStep(param,updateTime);
}


void NewmarkNLSIM::finalizeRHSvector (bool)
{
  // Add in the actual inertia force, computed from the equilibrium equation
  if (Finert) // RHS -= alphaH*Finert
    model.addToRHSvector(0,*Finert,gamma-0.5);
}


bool NewmarkNLSIM::predictStep (TimeStep& param)
{
  if (rotUpd && model.getNoFields(1) == 6)
  {
    // Initalize the total angular rotations for this time step
    Vector& psol = solution.front();
    if (psol.size() != 6*model.getNoNodes(1))
    {
      std::cerr <<" *** NewmarkNLSIM::predictStep: Invalid dimension on"
                <<" the displacement vector "<< psol.size()
                <<" != "<< 6*model.getNoNodes(1) << std::endl;
      return false;
    }
    for (size_t i = 3; i < psol.size(); i += 6)
      psol[i] = psol[i+1] = psol[i+2] = 0.0;
  }

#ifdef SP_DEBUG
  std::cout <<"\nNewmarkNLSIM::predictStep";
#endif

  // Predicted velocity, V_n = v_n-1*(gamma/beta-1) + a_n-1*dt*(gamma/beta-2)/2
  predVel = solution[iV];
  predVel *= gamma/beta - 1.0;
  predVel.add(solution[iA],param.time.dt*(gamma/beta-2.0)*0.5);

  // Predicted acceleration, A_n = a_n-1*(1/(2*beta)-1) + v_n-1*1/(dt*beta)
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
  std::cout <<"\nNewmarkNLSIM::correctStep(iter="<< param.iter
            <<",converged="<< std::boolalpha << converged <<")";
#endif

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
  std::cout <<"\nDisplacement increment:"<< incDis
            <<"Corrected displacement:"<< solution[iD]
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


void NewmarkNLSIM::setSolution (const RealArray& newSol, int idx)
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


bool NewmarkNLSIM::serialize (SerializeMap& data) const
{
  if (!this->saveSolution(data,model.getName()))
    return false;

  if (Finert)
    data["HHT::Finert"] = SIMsolution::serialize(Finert->getPtr(),
                                                 Finert->dim());

  return true;
}


bool NewmarkNLSIM::deSerialize (const SerializeMap& data)
{
  if (!this->restoreSolution(data,model.getName()))
    return false;

  SerializeMap::const_iterator sit = data.find("HHT::Finert");
  if (sit == data.end()) return false;

  Finert = new StdVector(model.getNoEquations());
  SIMsolution::deSerialize(sit->second,Finert->getPtr(),Finert->dim());

  return true;
}
