// $Id$
//==============================================================================
//!
//! \file NonLinSIM.C
//!
//! \date Jun 1 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear solution driver for isogeometric FEM simulators.
//!
//==============================================================================

#include "NonLinSIM.h"
#include "SIMoutput.h"
#include "IFEM.h"
#include "IntegrandBase.h"
#include "TimeStep.h"
#include "Profiler.h"
#include "Utilities.h"
#include "tinyxml.h"
#include <sstream>
#include <iomanip>

using namespace SIM;

const char* NonLinSIM::inputContext = "nonlinearsolver";


NonLinSIM::NonLinSIM (SIMbase& sim, CNORM n) : MultiStepSIM(sim), iteNorm(n)
{
  // Default solution parameters
  fromIni = iteNorm == NONE;
  maxIncr = 2;
  maxit   = 20;
  nupdat  = 20;
  prnSlow = 0;
  rTol    = 0.000001;
  aTol    = 0.0;
  divgLim = 10.0;
  alpha   = alphaO = 1.0;
  eta     = 0.0;
}


NonLinSIM::~NonLinSIM ()
{
  if (slowNodes.empty()) return;

  utl::LogStream& cout = model.getProcessAdm().cout;

  cout <<"\n *** Here are the nodal points flagged with slow convergence"
       <<"\n     ======================================================="
       <<"\n     Node  Count  Patch  Coordinates\n";
  for (const std::pair<const int,int>& node : slowNodes)
  {
    Vec4 X = model.getNodeCoord(node.first);
    cout << std::setw(9) << node.first
         << std::setw(5) << node.second
         << std::setw(7) << X.idx <<"   ";
    for (int i = 0; i < 3; i++)
      cout <<" "<< utl::trunc(X[i]);
    cout << std::endl;
  }
}


bool NonLinSIM::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"NONLINEAR_SOLVER",16))
  {
    std::istringstream cline(utl::readLine(is));
    cline >> maxit >> rTol;
    if (cline.fail() || cline.bad()) return false;

    double tmp;
    cline >> tmp;
    if (!cline.fail() && !cline.bad())
      divgLim = tmp;

    int itmp;
    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      nupdat = itmp;
    else
      nupdat = maxit;

    cline >> tmp;
    if (!cline.fail() && !cline.bad())
      eta = tmp;
  }
  else
    return model.parse(keyWord,is);

  return true;
}


bool NonLinSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),inputContext))
    return model.parse(elem);

  std::string rotUpdate;
  if (utl::getAttribute(elem,"rotation",rotUpdate,true) && !rotUpdate.empty())
    rotUpd = rotUpdate[0];

  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement()) {
    const char* value;
    if ((value = utl::getValue(child,"maxits")))
      maxit = atoi(value);
    else if ((value = utl::getValue(child,"maxIncr")))
      maxIncr = atoi(value);
    else if ((value = utl::getValue(child,"nupdate")))
      nupdat = atoi(value);
    else if ((value = utl::getValue(child,"rtol")))
      rTol = atof(value);
    else if ((value = utl::getValue(child,"atol")))
      aTol = atof(value);
    else if ((value = utl::getValue(child,"dtol")))
      divgLim = atof(value);
    else if ((value = utl::getValue(child,"eta")))
      eta = atof(value);
    else if ((value = utl::getValue(child,"printSlow")))
      prnSlow = atoi(value);
    else if ((value = utl::getValue(child,"referencenorm")))
    {
      if (!strcasecmp(value,"all"))
        refNopt = ALL;
      else if (!strcasecmp(value,"max"))
        refNopt = MAX;
    }
    else if (iteNorm != NONE && (value = utl::getValue(child,"convnorm")))
    {
      if (!strncasecmp(value,"energ",5))
        iteNorm = ENERGY;
      else if (!strncasecmp(value,"res",3))
        iteNorm = L2;
      else if (!strncasecmp(value,"dis",3))
        iteNorm = L2SOL;
      else if (!strncasecmp(value,"none",4))
      {
        iteNorm = NONE;
        fromIni = true;
      }
    }
    else if (!strcasecmp(child->Value(),"fromZero"))
      fromIni = true;
    else if (!strcasecmp(child->Value(),"printCond"))
      rCond = 0.0; // Compute and report condition number in the iteration log
  }

  return true;
}


void NonLinSIM::initPrm ()
{
  if (iteNorm == NONE) // Flag to integrand that a linear solver is used
    model.setIntegrationPrm(3,1.0);
}


void NonLinSIM::init (size_t nSol, const RealArray& value)
{
  this->MultiStepSIM::initSol(nSol);
  if (value.empty() || solution.empty())
    return;

  size_t ndim = solution.front().size();
  if (value.size() > ndim)
    std::copy(value.begin(),value.begin()+ndim,solution.front().begin());
  else
    std::copy(value.begin(),value.end(),solution.front().begin());
}


bool NonLinSIM::advanceStep (TimeStep& param, bool updateTime)
{
  this->pushSolution(); // Update solution vectors between time steps
  return this->MultiStepSIM::advanceStep(param,updateTime);
}


ConvStatus NonLinSIM::solve (double zero_tolerance, std::streamsize outPrec)
{
  TimeStep singleStep; // Solves the nonlinear equations in one single step
  return this->solveStep(singleStep,STATIC,zero_tolerance,outPrec);
}


ConvStatus NonLinSIM::solveStep (TimeStep& param, SolutionMode mode,
                                 double zero_tolerance, std::streamsize outPrec)
{
  PROFILE1("NonLinSIM::solveStep");

  if (solution.empty())
    return FAILURE;

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

  param.iter = 0;
  alpha = alphaO = 1.0;
  if (fromIni) // Always solve from initial configuration
    solution.front().fill(0.0);

  if (subiter&FIRST && !model.updateDirichlet(param.time.t,&solution.front()))
    return FAILURE;

  bool poorConvg = false;
  bool newTangent = param.time.first || iteNorm != NONE;
  model.setMode(newTangent ? mode : RHS_ONLY, false);
  model.setQuadratureRule(opt.nGauss[0],true);
  if (!this->assembleSystem(param.time,solution,newTangent))
    return model.getProblem()->diverged() ? DIVERGED : FAILURE;

  if (iteNorm != NONE)
    if (!model.extractLoadVec(residual))
      return FAILURE;

  double* rCondPtr = rCond < 0.0 ? nullptr : &rCond;
  if (!model.solveSystem(linsol,msgLevel-1,rCondPtr))
    return FAILURE;

  while (param.iter <= maxit)
    switch (this->checkConvergence(param))
      {
      case CONVERGED:
	if (!this->updateConfiguration(param))
	  return FAILURE;

	if (!this->solutionNorms(param.time,zero_tolerance,outPrec))
	  return FAILURE;

	param.time.first = false;
	return CONVERGED;

      case DIVERGED:
	return DIVERGED;

      case SLOW:
	poorConvg = true;

      default:
	param.iter++;
	if (!this->updateConfiguration(param))
	  return FAILURE;

	if (subiter&FIRST && param.iter == 1 && !model.updateDirichlet())
	  return FAILURE;

	if (param.iter > nupdat) newTangent = false;
	model.setMode(newTangent ? mode : RHS_ONLY, false);
	if (!this->assembleSystem(param.time,solution,newTangent,poorConvg))
	  return model.getProblem()->diverged() ? DIVERGED : FAILURE;

	if (!model.extractLoadVec(residual))
	  return FAILURE;

	if (!model.solveSystem(linsol,msgLevel-1,rCondPtr))
	  return FAILURE;

	if (!this->lineSearch(param))
	  return FAILURE;

	poorConvg = false;
      }

  return DIVERGED;
}


SIM::ConvStatus NonLinSIM::solveIteration (TimeStep& param)
{
  if (param.iter == 0 && !model.updateDirichlet(param.time.t,&solution.front()))
    return FAILURE;
  else if (param.iter == 1 && !model.updateDirichlet())
    return FAILURE;

  model.setMode(SIM::STATIC,false);

  if (!this->assembleSystem(param.time,solution))
    return SIM::FAILURE;

  if (!model.extractLoadVec(residual))
    return SIM::FAILURE;

  double* rCondPtr = rCond < 0.0 ? nullptr : &rCond;
  if (!model.solveSystem(linsol,msgLevel-1,rCondPtr))
    return SIM::FAILURE;

  if (!this->lineSearch(param))
    return SIM::FAILURE;

  if (!this->updateConfiguration(param))
    return SIM::FAILURE;

  SIM::ConvStatus result = this->checkConvergence(param);
  if (result == SIM::CONVERGED)
    if (!this->solutionNorms(param.time))
      return SIM::FAILURE;

  return result;
}


/*!
  This procedure is as described on pages 115,116 in Kjell Magne Mathisen's
  Dr.Ing. thesis: "Large displacement analysis of flexible and rigid systems
  considering displacement-dependent loads and nonlinear constraints", 1990.
*/

bool NonLinSIM::lineSearch (TimeStep& param)
{
  if (eta <= 0.0) return true; // No line search

  model.setMode(SIM::RHS_ONLY,false);

  double s0 = residual.dot(linsol);
  double smin = fabs(s0);
  double cmin = 1.0;
  double ck = 1.0;
  double cp = 0.0;
#ifdef SP_DEBUG
  std::cout << std::setw(8) << 0 <<": ck=0        sk="<< s0 << std::endl;
#endif

  for (int iter = 1; iter <= 10; iter++)
  {
    alpha = ck - cp;
    if (!this->updateConfiguration(param))
      return false;

    if (!this->assembleSystem(param.time,solution,false))
      return false;

    if (!model.extractLoadVec(residual))
      return false;

    double sk = residual.dot(linsol);
#ifdef SP_DEBUG
    std::cout << std::setw(8) << iter
              <<": ck="<< std::left << std::setw(9) << ck
              <<"sk="<< sk << std::right << std::endl;
#endif
    if (fabs(sk) < eta*fabs(s0))
    {
      alpha = 0.0;
      alphaO = ck;
      return true;
    }
    else if (fabs(sk) < smin)
    {
      smin = fabs(sk);
      cmin = ck;
    }

    cp = ck;
    ck *= s0/(s0-sk);
    if (ck > 1.5)
      ck = 1.5;
    else if (ck < 0.5)
      ck = 0.5;

    if (fabs(ck-cp) < 0.5*eta*fabs(ck+cp))
    {
      alpha = 0.0;
      alphaO = ck;
      return true;
    }
  }

  alpha = cmin - cp;
  alphaO = cmin;
  return true;
}


ConvStatus NonLinSIM::checkConvergence (TimeStep& param)
{
  if (iteNorm == NONE)
    return CONVERGED; // No iterations, we are solving a linear problem

  static double convTol   = 0.0;
  static double prevNorm  = 0.0;
  static int    nIncrease = 0;

  ConvStatus status = OK;
  double enorm, resNorm, linsolNorm;
  model.iterationNorms(linsol,residual,enorm,resNorm,linsolNorm);
  double norm = iteNorm == ENERGY ? enorm : resNorm;
  if (iteNorm == L2SOL) norm = linsolNorm;

  if (param.iter == 0)
  {
    if (linsolNorm == 0.0)
      return CONVERGED; // No load on this step

    if ((subiter&FIRST && refNopt == ALL) || fabs(norm) > refNorm)
      refNorm = fabs(norm);

    if (refNorm*rTol > aTol) {
      convTol = rTol;
      prevNorm = (norm /= refNorm);
    }
    else {
      convTol = aTol;
      refNorm = 1.0;
      prevNorm = norm;
    }

    nIncrease = 0;
  }
  else
    norm /= refNorm;

  // Check for slow convergence
  if (param.iter > 1 && prevNorm > 0.0 && fabs(norm) > prevNorm*0.1)
    status = SLOW;

  if (msgLevel > 0)
  {
    // Print convergence history
    utl::LogStream& cout = model.getProcessAdm().cout;
    std::ios::fmtflags oldFlags = cout.flags(std::ios::scientific);
    std::streamsize oldPrec = cout.precision(3);
    cout <<"  iter="<< param.iter
         <<"  conv="<< fabs(norm)
         <<"  enen="<< enorm
         <<"  resn="<< resNorm
         <<"  incn="<< linsolNorm;
    if (rCond > 0.0)
      cout <<"  cond="<< 1.0/rCond;
    if (alphaO != 1.0)
      cout <<"  alpha="<< alphaO;
    cout << std::endl;

    // Find and print out the worst DOF(s) when detecting slow convergence
    if (status == SLOW && prnSlow > 0)
      this->printWorst(cout,convTol*refNorm);

    cout.flags(oldFlags);
    cout.precision(oldPrec);
  }

  // Check for convergence or divergence
  if (fabs(norm) < convTol && (param.iter > 0 || refNopt == ALL))
    status = CONVERGED;
  else if (std::isnan(linsolNorm))
    status = DIVERGED;
  else if (fabs(norm) <= fabs(prevNorm))
    nIncrease = 0;
  else if (++nIncrease > maxIncr || fabs(norm) > divgLim)
    status = DIVERGED;

  prevNorm = norm;
  return status;
}


void NonLinSIM::printWorst (utl::LogStream& os, double eps)
{
  // Find and print out the worst DOF(s) when detecting slow convergence
  std::map<std::pair<int,int>,RealArray> worstDOFs;
  model.getWorstDofs(linsol,residual,prnSlow,eps,iteNorm,worstDOFs);

  os <<"  ** Slow convergence detected";
  if (worstDOFs.size() > 1)
    os <<", here are the "<< worstDOFs.size() <<" worst DOFs:";
  else if (worstDOFs.size() == 1)
    os <<", here is the worst DOF:";
  else
    os <<".";

  for (const std::pair<const std::pair<int,int>,RealArray>& wd : worstDOFs)
  {
    os <<"\n     Node "<< wd.first.first <<" local DOF "<< wd.first.second;
    char nodeType = model.getNodeType(wd.first.first);
    if (nodeType != ' ') os <<" ("<< nodeType <<")";
    os <<" :\tEnergy = "<< wd.second[0]
       <<"\tdu = "<< wd.second[1] <<"\tres = "<< wd.second[2];
    std::map<int,int>::iterator nit = slowNodes.find(wd.first.first);
    if (nit == slowNodes.end())
      slowNodes.insert(std::make_pair(wd.first.first,1));
    else
      ++nit->second;
  }
  os << std::endl;
}


bool NonLinSIM::updateConfiguration (TimeStep& time)
{
  if (solution.front().empty() || iteNorm == NONE)
  {
    solution.front() = linsol;
    model.updateRotations(linsol);
  }
  else if (alpha != 0.0)
  {
    if (time.iter == 1 && rotUpd && model.getNoFields(1) == 6)
    {
      // Initalize the total angular rotations for this load increment
      Vector& psol = solution.front();
      if (psol.size() != 6*model.getNoNodes(1))
      {
        std::cerr <<" *** NonLinSIM::updateConfiguration: Invalid dimension on"
                  <<" the displacement vector "<< psol.size()
                  <<" != "<< 6*model.getNoNodes(1) << std::endl;
        return false;
      }
      for (size_t i = 3; i < psol.size(); i += 6)
        psol[i] = psol[i+1] = psol[i+2] = 0.0;
    }

    solution.front().add(linsol,alpha);
    if (rotUpd == 't')
      model.updateRotations(solution.front());
    else if (rotUpd)
      model.updateRotations(linsol,alpha);
  }

  return model.updateConfiguration(solution.front());
}


bool NonLinSIM::assembleSystem(const TimeDomain& time, const Vectors& pSol,
                               bool newLHSmatrix, bool poorConvg)
{
  return model.assembleSystem(time, pSol, newLHSmatrix, poorConvg);
}
