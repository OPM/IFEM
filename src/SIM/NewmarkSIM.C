// $Id$
//==============================================================================
//!
//! \file NewmarkSIM.C
//!
//! \date Jul 4 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Newmark solution driver for isogeometric dynamic FEM simulators.
//!
//==============================================================================

#include "NewmarkSIM.h"
#include "SIMoutput.h"
#include "TimeStep.h"
#include "Profiler.h"
#include "Utilities.h"
#include "tinyxml.h"


NewmarkSIM::NewmarkSIM (SIMbase& sim) : MultiStepSIM(sim)
{
  // Default Newmark parameters
  beta = 0.25;
  gamma = 0.5;
  alpha1 = 0.0;
  alpha2 = 0.0;

  predictor = 'a'; // default predictor (zero acceleration)

  // Default iteration parameters
  maxit   = 20;
  convTol = 0.000001;
  divgLim = 10.0;
}


bool NewmarkSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"newmarksolver"))
    return model.parse(elem);

  utl::getAttribute(elem,"alpha1",alpha1);
  utl::getAttribute(elem,"alpha2",alpha2);
  utl::getAttribute(elem,"beta",beta);
  utl::getAttribute(elem,"gamma",gamma);

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
    else if ((value = utl::getValue(child,"predictor")))
    {
      if (!strncasecmp(value,"constant dis",12))
        predictor = 'd';
      else if (!strncasecmp(value,"constant vel",12))
        predictor = 'v';
      else if (!strncasecmp(value,"zero acc",8))
        predictor = 'a';
    }

  return true;
}


void NewmarkSIM::printProblem (std::ostream& os) const
{
  model.printProblem(os);

  os <<"Newmark predictor/multicorrector: beta = "<< beta <<" gamma = "<< gamma;
  switch (predictor) {
  case 'd': os <<"\n- using constant displacement predictor"; break;
  case 'v': os <<"\n- using constant velocity predictor"; break;
  case 'a': os <<"\n- using zero acceleration predictor"; break;
  }

  if (alpha1 > 0.0)
    os <<"\nMass-proportional damping (alpha1): "<< alpha1;
  if (alpha2 != 0.0)
    os <<"\nStiffness-proportional damping (alpha2): "<< fabs(alpha2);

  os << std::endl;
}


void NewmarkSIM::init (size_t nSol)
{
  model.setIntegrationPrm(0,alpha1);
  model.setIntegrationPrm(1,fabs(alpha2));
  model.setIntegrationPrm(2,beta);
  model.setIntegrationPrm(3,gamma);

  this->MultiStepSIM::init(nSol);
}


bool NewmarkSIM::advanceStep (TimeStep& param, bool updateTime)
{
  return updateTime ? param.increment() : true;
}


bool NewmarkSIM::predictStep (TimeStep& param)
{
  if (solution.size() < 3)
  {
    std::cerr <<" *** NewmarkSIM::predictStep: Too few solution vectors "
              << solution.size() << std::endl;
    return false;
  }

  const double dt = param.time.dt;

  size_t iD = 0;
  size_t iA = solution.size() - 1;
  size_t iV = solution.size() - 2;
  Vector oldSol(solution[iD]);

  switch (predictor) {
  case 'a': // zero acceleration predictor

    // Predicted new displacement
    solution[iD].add(solution[iV],dt);
    solution[iD].add(solution[iA],dt*dt*(0.5-beta));

    // Predicted new velocity
    solution[iV].add(solution[iA],dt*(1.0-gamma));

    // Predicted new acceleration (zero)
    solution[iA].fill(0.0);
    break;

  case 'v': // constant velocity predictor

    // Predicted new displacement
    solution[iD].add(solution[iV],dt);
    solution[iD].add(solution[iA],dt*dt*(0.5-beta/gamma));

    // Predicted new acceleration
    solution[iA] *= 1.0 - 1.0/gamma;
    break;

  case 'd': // constant displacement predictor
    oldSol = solution[iV];

    // Predicted new velocity
    solution[iV] *= 1.0 - gamma/beta;
    solution[iV].add(solution[iA],dt*(1.0 - 0.5*gamma/beta));

    // Predicted new acceleration
    solution[iA] *= 1.0 - 0.5/beta;
    solution[iA].add(oldSol,-1.0/(beta*dt));
    break;

  default:
    return false;
  }

#if SP_DEBUG > 1
  std::cout <<"Predicted displacement:"<< solution[iD];
  std::cout <<"Predicted velocity:"<< solution[iV];
  std::cout <<"Predicted acceleration:"<< solution[iA];
#endif

  if (predictor == 'd') return true;

  model.updateRotations(solution[iD]-oldSol,1.0);

  return model.updateConfiguration(solution[iD]);
}


bool NewmarkSIM::correctStep (TimeStep& param, bool)
{
  const double dt = param.time.dt;

  size_t iD = 0;
  size_t iA = solution.size() - 1;
  size_t iV = solution.size() - 2;

  // Corrected displacement
  solution[iD].add(linsol,beta*dt*dt);
  model.updateRotations(linsol,beta*dt*dt);

  // Corrected velocity
  solution[iV].add(linsol,gamma*dt);

  // Corrected acceleration
  solution[iA].add(linsol,1.0);

#if SP_DEBUG > 1
  std::cout <<"Corrected displacement:"<< solution[iD];
  std::cout <<"Corrected velocity:"<< solution[iV];
  std::cout <<"Corrected acceleration:"<< solution[iA];
#endif

  return model.updateConfiguration(solution[iD]);
}


SIM::ConvStatus NewmarkSIM::solveStep (TimeStep& param, SIM::SolutionMode,
                                       double zero_tolerance,
                                       std::streamsize outPrec)
{
  PROFILE1("NewmarkSIM::solveStep");

  if (msgLevel >= 0 && myPid == 0)
  {
    std::streamsize oldPrec = 0;
    double digits = log10(param.time.t)-log10(param.time.dt);
    if (digits > 6.0) oldPrec = std::cout.precision(ceil(digits));
    std::cout <<"\n  step="<< param.step
              <<"  time="<< param.time.t << std::endl;
    if (digits > 6.0) std::cout.precision(oldPrec);
  }

  if (subiter&FIRST && !model.updateDirichlet(param.time.t,&solution.front()))
    return SIM::FAILURE;

  param.iter = 0;
  if (subiter&FIRST && !this->predictStep(param))
    return SIM::FAILURE;

  if (!model.setMode(SIM::DYNAMIC))
    return SIM::FAILURE;

  model.setQuadratureRule(opt.nGauss[0],true);
  if (!model.assembleSystem(param.time,solution))
    return SIM::FAILURE;

  this->finalizeRHSvector();

  if (!model.extractLoadVec(residual))
    return SIM::FAILURE;

  if (!model.solveSystem(linsol,msgLevel-1))
    return SIM::FAILURE;

  while (param.iter <= maxit)
    switch (this->checkConvergence(param))
      {
      case SIM::CONVERGED:
        if (!this->correctStep(param,subiter&LAST))
          return SIM::FAILURE;

        if (!this->solutionNorms(zero_tolerance,outPrec))
          return SIM::FAILURE;

        if (subiter&LAST) param.time.first = false;
        return SIM::CONVERGED;

      case SIM::DIVERGED:
        return SIM::DIVERGED;

      default:
        param.iter++;
        if (!this->correctStep(param))
          return SIM::FAILURE;

        if (subiter&FIRST && param.iter == 1 && !model.updateDirichlet())
          return SIM::FAILURE;

        if (!model.assembleSystem(param.time,solution))
          return SIM::FAILURE;

        this->finalizeRHSvector();

        if (!model.extractLoadVec(residual))
          return SIM::FAILURE;

        if (!model.solveSystem(linsol,msgLevel-1))
          return SIM::FAILURE;
      }

  return SIM::DIVERGED;
}


SIM::ConvStatus NewmarkSIM::checkConvergence (TimeStep& param)
{
  static double prevNorm  = 0.0;
  static int    nIncrease = 0;

  SIM::ConvStatus status = SIM::OK;
  double enorm, resNorm, linsolNorm;
  model.iterationNorms(linsol,residual,enorm,resNorm,linsolNorm);

  double norm = 1.0;
  if (param.iter > 0)
    norm = resNorm / refNorm;
  else
  {
    if (refNopt == ALL || fabs(resNorm) > refNorm)
      refNorm = fabs(resNorm);
    else
      norm = resNorm / refNorm;
    prevNorm = norm;
    nIncrease = 0;
  }

  if (msgLevel > 0 && myPid == 0)
  {
    // Print convergence history
    std::ios::fmtflags oldFlags = std::cout.flags(std::ios::scientific);
    std::streamsize oldPrec = std::cout.precision(3);
    std::cout <<"  iter="<< param.iter
              <<"  conv="<< fabs(norm)
              <<"  enen="<< enorm
              <<"  resn="<< resNorm
              <<"  incn="<< linsolNorm << std::endl;
    std::cout.flags(oldFlags);
    std::cout.precision(oldPrec);
  }

  // Check for convergence or divergence
  if (fabs(norm) < convTol)
    status = SIM::CONVERGED;
  else if (fabs(norm) <= fabs(prevNorm))
    nIncrease = 0;
  else if (++nIncrease > 2 || fabs(norm) > divgLim)
    status = SIM::DIVERGED;

  prevNorm = norm;
  return status;
}


bool NewmarkSIM::solutionNorms (double zero_tolerance, std::streamsize outPrec)
{
  if (msgLevel < 0 || solution.size() < 3) return true;

  // Cannot use the enums here because this method is inherited
  size_t a = solution.size()-1;
  size_t v = solution.size()-2;

  size_t d, nf = model.getNoFields(1);
  size_t iMax[nf], jMax[nf], kMax[nf];
  double dMax[nf], vMax[nf], aMax[nf];
  double disL2 = model.solutionNorms(solution[0],dMax,iMax,nf);
  double velL2 = model.solutionNorms(solution[v],vMax,jMax,nf);
  double accL2 = model.solutionNorms(solution[a],aMax,kMax,nf);

  if (myPid > 0) return true;

  std::streamsize stdPrec = outPrec > 0 ? std::cout.precision(outPrec) : 0;
  double old_tol = utl::zero_print_tol;
  utl::zero_print_tol = zero_tolerance;
  char D;

  std::cout <<"  Displacement L2-norm            : "<< utl::trunc(disL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(dMax[d]) != 0.0)
      std::cout <<"\n               Max "<< D
                <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];

  std::cout <<"\n  Velocity L2-norm                : "<< utl::trunc(velL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(vMax[d]) != 0.0)
      std::cout <<"\n               Max "<< D
                <<"-velocity     : "<< vMax[d] <<" node "<< jMax[d];

  std::cout <<"\n  Acceleration L2-norm            : "<< utl::trunc(accL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(aMax[d]) != 0.0)
      std::cout <<"\n               Max "<< D
                <<"-acceleration : "<< aMax[d] <<" node "<< kMax[d];

  std::cout << std::endl;
  utl::zero_print_tol = old_tol;
  if (stdPrec > 0) std::cout.precision(stdPrec);

  return true;
}


void NewmarkSIM::dumpResults (double time, std::ostream& os,
                              std::streamsize precision, bool formatted) const
{
  model.dumpResults(solution.front(),time,os,formatted,precision);
  model.dumpMoreResults(time,os,precision);
  if (formatted)
    model.dumpVector(this->getVelocity(),"velocity",os,precision);
}
