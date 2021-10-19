// $Id$
//==============================================================================
//!
//! \file HHTSIM.C
//!
//! \date Nov 13 2014
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Newmark solution driver for isogeometric dynamic FEM simulators.
//!
//==============================================================================

#include "HHTSIM.h"
#include "SIMoutput.h"
#include "SystemMatrix.h"
#include "TimeStep.h"
#include "IFEM.h"
#include "tinyxml.h"


HHTSIM::HHTSIM (SIMbase& sim) : NewmarkSIM(sim), Finert(nullptr), Fext(nullptr)
{
  nRHSvec = 3; // 0: residual force, 1: inertia force, 2: external forces

  // Default Newmark parameters (alpha = -0.1)
  beta = 0.3025;
  gamma = 0.6;

  predictor = 'd'; // default predictor (constant displacement)

  pV = pA = iD = iV = iA = 0;
}


bool HHTSIM::parse (const TiXmlElement* elem)
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
  }

  return ok;
}


void HHTSIM::printProblem () const
{
  this->NewmarkSIM::printProblem();

  if (alpha2 > 0.0)
    IFEM::cout <<"- based on the tangential stiffness matrix\n";
  else if (alpha2 < 0.0)
    IFEM::cout <<"- based on the material stiffness matrix only\n";
  else
    return;

  IFEM::cout << std::endl;
}


void HHTSIM::initPrm ()
{
  model.setIntegrationPrm(0,alpha1);
  model.setIntegrationPrm(1,fabs(alpha2));
  model.setIntegrationPrm(2,0.5-gamma);
  if (alpha2 < 0.0) // Flag that stiffness-proportional damping should depend
    model.setIntegrationPrm(3,-1.0); // on the material stiffness matrix only
  // Flag to the integrand that the Hilber-Hughes-Taylor algorithm is used
  model.setIntegrationPrm(4,1.0);
}


bool HHTSIM::initSol (size_t nSol)
{
  incDis.resize(model.getNoDOFs(),true);
  this->MultiStepSIM::initSol(nSol > 0 ? nSol+2 : 0);
  if (solution.size() < 5)
  {
    std::cerr <<" *** HHTSIM::initSol: Too few solution vectors "
              << solution.size() << std::endl;
    return false;
  }

  iA = solution.size() - 1;
  iV = solution.size() - 2;
  pA = solution.size() - 3;
  pV = solution.size() - 4;

  return true;
}


bool HHTSIM::advanceStep (TimeStep& param, bool updateTime)
{
  // Update displacement solutions between time steps
  if (solution.size() > 5)
    this->pushSolution(solution.size()-4);

  return this->NewmarkSIM::advanceStep(param,updateTime);
}


/*!
  \param[in] predicting If \e true, this is the first iteration of time step > 1

  This method merges the nodal point loads into the right-hand-side vector of
  the linearized dynamic equilibrium equation, according to the HHT scheme.
  It also adds the actual inertia force \e *Finert which has been computed
  from the dynamic equilibrium equation.
*/

void HHTSIM::finalizeRHSvector (bool predicting)
{
  if (predicting)
  {
    if (Fext && Fext->L1norm() > 0.0)
      model.addToRHSvector(1, *Fext); // Finert += R_ext,n-1

    delete Finert;
    Finert = model.getRHSvector(1,true);
#if SP_DEBUG > 1
    if (Finert)
      std::cout <<"\nActual inertia force:"<< *Finert;
#endif
  }

  // Add in the actual inertia force, computed from the equilibrium equation
  if (Finert) // RHS += Finert (predictor), or RHS -= alphaH*Finert (corrector)
    model.addToRHSvector(0, *Finert, predicting ? 1.0 : gamma-0.5);

  // Add in the external nodal loads, if any
  SystemVector* Rext = model.getRHSvector(2);
  if (Rext && Rext->L1norm() > 0.0)
    model.addToRHSvector(0, *Rext, 1.5-gamma); // RHS += (1+alphaH)*R_ext,n
  if (Fext && predicting && Fext->L1norm() > 0.0)
    model.addToRHSvector(0, *Fext, gamma-1.5); // RHS -= (1+alphaH)*R_ext,n-1

#if SP_DEBUG > 1
  if (Rext && predicting && Rext->L1norm() > 0.0)
    std::cout <<"\nExternal nodal loads, this time step:"<< *Rext;
  if (Fext && predicting && Fext->L1norm() > 0.0)
    std::cout <<"\nExternal nodal loads, last time step:"<< *Fext;
#endif
}


bool HHTSIM::predictStep (TimeStep& param)
{
  if (rotUpd && model.getNoFields(1) == 6)
  {
    // Initalize the total angular rotations for this time step
    Vector& psol = solution.front();
    if (psol.size() != 6*model.getNoNodes(1))
    {
      std::cerr <<" *** HHTSIM::predictStep: Invalid dimension on"
                <<" the displacement vector "<< psol.size()
                <<" != "<< 6*model.getNoNodes(1) << std::endl;
      return false;
    }
    for (size_t i = 3; i < psol.size(); i += 6)
      psol[i] = psol[i+1] = psol[i+2] = 0.0;
  }

#ifdef SP_DEBUG
  std::cout <<"\nHHTSIM::predictStep";
#endif

  // Predicted velocity, V_n = v_n-1*gamma/beta + a_n-1*dt*(gamma/beta-2)/2
  solution[pV] = solution[iV];
  solution[pV] *= gamma/beta;
  solution[pV].add(solution[iA],param.time.dt*(gamma/beta-2.0)*0.5);

  // Predicted acceleration, A_n = a_n-1*1/(2*beta) + v_n-1*1/(dt*beta)
  solution[pA] = solution[iA];
  solution[pA] *= 0.5/beta;
  solution[pA].add(solution[iV],1.0/(beta*param.time.dt));

#ifdef SP_DEBUG
  std::cout <<"\nPredicted velocity:";
#if SP_DEBUG > 1
  std::cout << solution[pV];
#else
  std::cout << solution[pV].max() <<"\n";
#endif
  std::cout <<"Predicted acceleration:";
#if SP_DEBUG > 1
  std::cout << solution[pA];
#else
  std::cout << solution[pA].max() << std::endl;
#endif
#endif

  incDis.fill(0.0);
  return true;
}


bool HHTSIM::correctStep (TimeStep& param, bool converged)
{
#ifdef SP_DEBUG
  std::cout <<"\nHHTSIM::correctStep(iter="<< param.iter
            <<",converged="<< std::boolalpha << converged <<")";
#endif

  if (param.iter == 1 && !converged)
  {
    // V_n = v_n^1 = v_{n-1} - V_n
    // A_n = a_n^1 = a_{n-1} - A_n
    solution[pV] = solution[iV].add(solution[pV],-1.0);
    solution[pA] = solution[iA].add(solution[pA],-1.0);
  }
  else
  {
    // v_n^i = V_n
    // a_n^i = A_n
    solution[iV] = solution[pV];
    solution[iA] = solution[pA];
  }

  // Update current displacement, velocity and acceleration solutions
  incDis.add(linsol,1.0);
  solution[iD].add(linsol,1.0);
  solution[iV].add(incDis,gamma/(beta*param.time.dt));
  solution[iA].add(incDis,1.0/(beta*param.time.dt*param.time.dt));

  if (converged)
  {
    // Save the external nodal loads from converged step
    delete Fext;
    Fext = model.getRHSvector(2,true);
  }

#if SP_DEBUG > 1
  std::cout <<"\nDisplacement increment:"<< incDis
            <<"Corrected displacement:"<< solution[iD]
            <<"Corrected velocity:"<< solution[iV]
            <<"Corrected acceleration:"<< solution[iA];
  if (converged && Fext)
    std::cout <<"External force vector:"<< *Fext;
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


void HHTSIM::setSolution (const RealArray& newSol, int idx)
{
  if (idx == 0)
  {
    // When updating the displacements within sub-iterations, we must also
    // update the incremental displacement accordingly (well spotted, Knut)
    incDis.add(solution.front(),-1.0);
    incDis.add(newSol,1.0);
    if (rotUpd)
    {
      // TODO: Must handle rotations also, check if this is correct for GETBEAM
      model.updateRotations(solution.front(),-1.0);
      model.updateRotations(newSol,1.0);
    }
  }

  solution[idx] = newSol;
}


bool HHTSIM::serialize (SerializeMap& data) const
{
  if (!this->saveSolution(data,model.getName()))
    return false;

  if (Finert)
    data["HHT::Finert"] = SIMsolution::serialize(Finert->getPtr(),
                                                 Finert->dim());

  if (Fext)
    data["HHT::Fext"]   = SIMsolution::serialize(Fext->getPtr(),
                                                 Fext->dim());

  return true;
}


bool HHTSIM::deSerialize (const SerializeMap& data)
{
  if (!this->restoreSolution(data,model.getName()))
    return false;

  SerializeMap::const_iterator sit = data.find("HHT::Finert");
  if (sit == data.end()) return false;

  Finert = new StdVector(model.getNoEquations());
  SIMsolution::deSerialize(sit->second,Finert->getPtr(),Finert->dim());

  sit = data.find("HHT::Fext");
  if (sit == data.end()) return false;

  Fext = new StdVector(model.getNoEquations());
  SIMsolution::deSerialize(sit->second,Fext->getPtr(),Fext->dim());

  return true;
}
