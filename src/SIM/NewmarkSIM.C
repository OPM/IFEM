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
#include "AlgEqSystem.h"
#include "TimeStep.h"
#include "IFEM.h"
#include "Profiler.h"
#include "Utilities.h"
#include "tinyxml.h"

const char* NewmarkSIM::inputContext = "newmarksolver";


NewmarkSIM::NewmarkSIM (SIMbase& sim) : MultiStepSIM(sim)
{
  // Default Newmark parameters
  beta = 0.25;
  gamma = 0.5;
  alpha1 = 0.0;
  alpha2 = 0.0;

  solveDisp = false; // default use acceleration as primary variables
  predictor = 'a'; // default predictor (zero acceleration)
  cNorm = 1; // default convergence check, force residual

  // Default iteration parameters
  maxIncr = 2;
  maxit   = 20;
  rTol    = 1.0e-6;
  aTol    = 0.0;
  divgLim = 10.0;
  saveIts = 0;
}


bool NewmarkSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),inputContext))
    return model.parse(elem);

  std::string rotUpdate;
  if (utl::getAttribute(elem,"rotation",rotUpdate,true) && !rotUpdate.empty())
    rotUpd = rotUpdate[0];

  utl::getAttribute(elem,"alpha1",alpha1);
  utl::getAttribute(elem,"alpha2",alpha2);
  utl::getAttribute(elem,"beta",beta);
  utl::getAttribute(elem,"gamma",gamma);

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement()) {
    const char* value;
    if ((value = utl::getValue(child,"maxits")))
      maxit = atoi(value);
    else if ((value = utl::getValue(child,"maxIncr")))
      maxIncr = atoi(value);
    else if ((value = utl::getValue(child,"rtol")))
      rTol = atof(value);
    else if ((value = utl::getValue(child,"atol")))
      aTol = atof(value);
    else if ((value = utl::getValue(child,"dtol")))
      divgLim = atof(value);
    else if ((value = utl::getValue(child,"saveiterations")))
      saveIts = atoi(value);
    else if ((value = utl::getValue(child,"referencenorm")))
    {
      if (!strcasecmp(value,"all"))
        refNopt = ALL;
      else if (!strcasecmp(value,"max"))
        refNopt = MAX;
    }
    else if ((value = utl::getValue(child,"convnorm")))
    {
      if (!strncasecmp(value,"energ",5))
        cNorm = 0;
      else if (!strncasecmp(value,"res",3))
        cNorm = 1;
      else if (!strncasecmp(value,"dis",3))
        cNorm = 2;
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
    else if ((value = utl::getValue(child,"rotation")))
      rotUpd = tolower(value[0]);
    else if (!strncasecmp(child->Value(),"solve_dis",9))
      solveDisp = true; // no need for value here
    else if (!strcasecmp(child->Value(),"printCond"))
      rCond = 0.0;
  }

  return true;
}


void NewmarkSIM::printProblem () const
{
  model.printProblem();

  IFEM::cout <<"Newmark predictor/multicorrector: beta = "<< beta
             <<" gamma = "<< gamma;
  switch (predictor) {
  case 'd': IFEM::cout <<"\n- using constant displacement predictor"; break;
  case 'v': IFEM::cout <<"\n- using constant velocity predictor"; break;
  case 'a': IFEM::cout <<"\n- using zero acceleration predictor"; break;
  }
  if (solveDisp)
    IFEM::cout <<"\n- using displacement increments as primary unknowns";
  if (alpha1 > 0.0)
    IFEM::cout <<"\nMass-proportional damping (alpha1): "<< alpha1;
  if (alpha2 != 0.0)
    IFEM::cout <<"\nStiffness-proportional damping (alpha2): "<< fabs(alpha2);

  IFEM::cout << std::endl;
}


void NewmarkSIM::initPrm ()
{
  model.setIntegrationPrm(0,alpha1);
  model.setIntegrationPrm(1,fabs(alpha2));
  model.setIntegrationPrm(2,solveDisp ? -beta : beta);
  model.setIntegrationPrm(3,gamma);
}


bool NewmarkSIM::initAcc (double zero_tolerance, std::streamsize outPrec)
{
  if (!model.setMode(SIM::MASS_ONLY))
    return false;

  // Assemble mass matrix and external forces
  model.setQuadratureRule(opt.nGauss[0],true);
  if (!model.assembleSystem())
    return false;

  // Add in the external nodal loads, if any
  SystemVector* Rext = model.getRHSvector(2);
  if (Rext && Rext->L1norm() > 0.0)
    model.addToRHSvector(0,*Rext);

  // Solve for the initial accelerations
  size_t iA = solution.size() - 1;
  if (!model.solveSystem(solution[iA],msgLevel-1,"acceleration"))
    return false;
  else if (msgLevel < 1)
    return true;

  size_t d, nf = model.getNoFields(1);
  size_t kMax[nf];
  double aMax[nf];
  double accL2 = model.solutionNorms(solution[iA],aMax,kMax,nf);

  utl::LogStream& cout = model.getProcessAdm().cout;
  std::streamsize stdPrec = outPrec > 0 ? cout.precision(outPrec) : 0;
  double old_tol = utl::zero_print_tol;
  utl::zero_print_tol = zero_tolerance;
  char D;

  cout <<"\n  Initial acceleration L2-norm    : "<< utl::trunc(accL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(aMax[d]) != 0.0)
      cout <<"\n               Max "<< D
           <<"-acceleration : "<< aMax[d] <<" node "<< kMax[d];

  cout << std::endl;
  utl::zero_print_tol = old_tol;
  if (stdPrec > 0) cout.precision(stdPrec);

  return true;
}


bool NewmarkSIM::predictStep (TimeStep& param)
{
  if (solution.size() < 3)
  {
    std::cerr <<" *** NewmarkSIM::predictStep: Too few solution vectors "
              << solution.size() << std::endl;
    return false;
  }

  if (rotUpd && model.getNoFields(1) == 6)
  {
    // Initalize the total angular rotations for this time step
    Vector& psol = solution.front();
    if (psol.size() != 6*model.getNoNodes(1))
    {
      std::cerr <<" *** NewmarkSIM::predictStep: Invalid dimension on"
                <<" the displacement vector "<< psol.size()
                <<" != "<< 6*model.getNoNodes(1) << std::endl;
      return false;
    }
    for (size_t i = 3; i < psol.size(); i += 6)
      psol[i] = psol[i+1] = psol[i+2] = 0.0;
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

  if (rotUpd == 't')
    model.updateRotations(solution[iD]);
  else if (rotUpd)
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
  solution[iD].add(linsol, solveDisp ? 1.0 : beta*dt*dt);

  // Corrected velocity
  solution[iV].add(linsol, solveDisp ? gamma/(beta*dt) : gamma*dt);

  // Corrected acceleration
  solution[iA].add(linsol, solveDisp ? 1.0/(beta*dt*dt) : 1.0);

#if SP_DEBUG > 1
  std::cout <<"Corrected displacement:"<< solution[iD];
  std::cout <<"Corrected velocity:"<< solution[iV];
  std::cout <<"Corrected acceleration:"<< solution[iA];
#endif

  if (rotUpd == 't')
    model.updateRotations(solution[iD]);
  else if (rotUpd)
    model.updateRotations(linsol, solveDisp ? 1.0 : beta*dt*dt);

  return model.updateConfiguration(solution[iD]);
}


SIM::ConvStatus NewmarkSIM::solveStep (TimeStep& param, SIM::SolutionMode,
                                       double zero_tolerance,
                                       std::streamsize outPrec)
{
  PROFILE1("NewmarkSIM::solveStep");

  if (msgLevel >= 0)
  {
    double digits = log10(param.time.t) - log10(param.time.dt);
    if (digits > 6.0)
    {
      utl::LogStream& cout = model.getProcessAdm().cout;
      std::streamsize oldPrec = cout.precision(ceil(digits));
      model.printStep(param.step,param.time);
      cout.precision(oldPrec);
    }
    else
      model.printStep(param.step,param.time);
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

  this->finalizeRHSvector(!param.time.first);

  if (!model.extractLoadVec(residual))
    return SIM::FAILURE;

  double* rCondPtr = rCond < 0.0 ? nullptr : &rCond;
  if (!model.solveSystem(linsol,msgLevel-1,rCondPtr))
    return SIM::FAILURE;

  while (param.iter <= maxit)
    switch (this->checkConvergence(param))
      {
      case SIM::CONVERGED:
        if (!this->correctStep(param,subiter&LAST))
          return SIM::FAILURE;

        if (!this->solutionNorms(param.time,zero_tolerance,outPrec))
          return SIM::FAILURE;

        if (subiter&LAST) param.time.first = false;
        return SIM::CONVERGED;

      case SIM::DIVERGED:
        model.getProcessAdm().cout <<" *** Iterations diverged, terminating..."
                                   << std::endl;
        return SIM::DIVERGED;

      default:
        param.iter++;
        if (!this->correctStep(param))
          return SIM::FAILURE;

        if (param.step == saveIts)
        {
          double time = param.time.t + (param.time.dt*param.iter)/maxit;
          if (!this->saveStep(this->getLastSavedStep()+1,time))
            return SIM::FAILURE;
          else if (!model.setMode(SIM::DYNAMIC))
            return SIM::FAILURE;
        }

        if (subiter&FIRST && param.iter == 1 && !model.updateDirichlet())
          return SIM::FAILURE;

        if (!model.assembleSystem(param.time,solution))
          return SIM::FAILURE;

        this->finalizeRHSvector(false);

        if (!model.extractLoadVec(residual))
          return SIM::FAILURE;

        if (!model.solveSystem(linsol,msgLevel-1,rCondPtr))
          return SIM::FAILURE;
      }

  model.getProcessAdm().cout <<" *** No convergence in "<< maxit
                             <<" iterations, terminating..."<< std::endl;
  return SIM::DIVERGED;
}


SIM::ConvStatus NewmarkSIM::solveIteration (TimeStep& param)
{
  bool ok = true;
  if (param.iter == 0)
    ok = model.updateDirichlet(param.time.t,&solution.front());
  else if (param.iter == 1)
    ok = model.updateDirichlet();

  if (param.iter > 0)
    ok &= this->correctStep(param);
  else
    ok &= this->predictStep(param);
  if (!ok)
    return SIM::FAILURE;

  if (!model.setMode(SIM::DYNAMIC))
    return SIM::FAILURE;

  if (!model.assembleSystem(param.time,solution))
    return SIM::FAILURE;

  this->finalizeRHSvector(!param.time.first && param.iter == 0);

  if (!model.extractLoadVec(residual))
    return SIM::FAILURE;

  double* rCondPtr = rCond < 0.0 ? nullptr : &rCond;
  if (!model.solveSystem(linsol,msgLevel-1,rCondPtr))
    return SIM::FAILURE;

  SIM::ConvStatus result = this->checkConvergence(param);
  if (result == SIM::CONVERGED)
    if (!this->solutionNorms(param.time))
      return SIM::FAILURE;

  return result;
}


SIM::ConvStatus NewmarkSIM::checkConvergence (TimeStep& param)
{
  static double convTol   = 0.0;
  static double prevNorm  = 0.0;
  static int    nIncrease = 0;

  double norms[3];
  model.iterationNorms(linsol,residual,norms[0],norms[1],norms[2]);
  double norm = norms[cNorm];

  if (param.iter == 0)
  {
    if (norms[2] == 0.0)
      return SIM::CONVERGED; // No load on this step

    if ((subiter&FIRST && refNopt == ALL) || fabs(norm) > refNorm)
      refNorm = fabs(norm);

    if (refNorm*rTol > aTol) {
      convTol = rTol;
      norm /= refNorm;
    }
    else {
      convTol = aTol;
      refNorm = 1.0;
    }

    prevNorm = norm;
    nIncrease = 0;
  }
  else
    norm /= refNorm;

  if (msgLevel > 0)
  {
    // Print convergence history
    utl::LogStream& cout = model.getProcessAdm().cout;
    std::ios::fmtflags stdFlags = cout.flags(std::ios::scientific);
    std::streamsize stdPrec = cout.precision(3);
    cout <<"  iter="<< param.iter
         <<"  conv="<< fabs(norm)
         <<"  enen="<< norms[0]
         <<"  resn="<< norms[1]
         <<"  incn="<< norms[2];
    if (rCond > 0.0)
      cout <<"  cond="<< 1.0/rCond;
    cout << std::endl;
    cout.flags(stdFlags);
    cout.precision(stdPrec);
  }

  // Check for convergence or divergence
  SIM::ConvStatus status = SIM::OK;
  if (fabs(norm) < convTol && (param.iter > 0 || refNopt == ALL))
    status = SIM::CONVERGED;
  else if (std::isnan(norms[2]))
    status = SIM::DIVERGED;
  else if (fabs(norm) <= fabs(prevNorm))
    nIncrease = 0;
  else if (++nIncrease > maxIncr || fabs(norm) > divgLim)
    status = SIM::DIVERGED;

  prevNorm = norm;
  return status;
}


bool NewmarkSIM::solutionNorms (const TimeDomain&,
                                double zero_tolerance, std::streamsize outPrec)
{
  // Establish the real FE solution from the time integration solution space
  const Vectors& newSol = this->realSolutions();
  if (msgLevel < 0 || newSol.size() < 3)
    return true;

  // Cannot use the enums here because this method is inherited
  size_t a = newSol.size()-1;
  size_t v = newSol.size()-2;

  size_t d, nf = model.getNoFields(1);
  size_t iMax[nf], jMax[nf], kMax[nf];
  double dMax[nf], vMax[nf], aMax[nf];
  double disL2 = model.solutionNorms(newSol[0],dMax,iMax,nf);
  double velL2 = model.solutionNorms(newSol[v],vMax,jMax,nf);
  double accL2 = model.solutionNorms(newSol[a],aMax,kMax,nf);

  utl::LogStream& cout = model.getProcessAdm().cout;
  std::streamsize stdPrec = outPrec > 0 ? cout.precision(outPrec) : 0;
  double old_tol = utl::zero_print_tol;
  utl::zero_print_tol = zero_tolerance;
  char D;

  cout <<"  Displacement L2-norm            : "<< utl::trunc(disL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(dMax[d]) != 0.0)
      cout <<"\n               Max "<< D
           <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];

  cout <<"\n  Velocity L2-norm                : "<< utl::trunc(velL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(vMax[d]) != 0.0)
      cout <<"\n               Max "<< D
           <<"-velocity     : "<< vMax[d] <<" node "<< jMax[d];

  cout <<"\n  Acceleration L2-norm            : "<< utl::trunc(accL2);
  for (d = 0, D = 'X'; d < nf; d++, D=='Z' ? D='x' : D++)
    if (utl::trunc(aMax[d]) != 0.0)
      cout <<"\n               Max "<< D
           <<"-acceleration : "<< aMax[d] <<" node "<< kMax[d];

  cout << std::endl;
  utl::zero_print_tol = old_tol;
  if (stdPrec > 0) cout.precision(stdPrec);

  return true;
}


void NewmarkSIM::dumpResults (double time, utl::LogStream& os,
                              std::streamsize precision, bool formatted) const
{
  model.dumpResults(this->realSolution(),time,os,formatted,precision);
  model.dumpMoreResults(time,os,precision);
  if (formatted)
  {
    model.dumpVector(this->getVelocity(),"velocity",os,precision);
    model.dumpVector(this->getAcceleration(),"acceleration",os,precision);
  }
}
