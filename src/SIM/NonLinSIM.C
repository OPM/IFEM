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


NonLinSIM::NonLinSIM (SIMbase& sim, CNORM n) : MultiStepSIM(sim), iteNorm(n)
{
  // Default solution parameters
  fromIni = iteNorm == NONE;
  maxit   = 20;
  nupdat  = 20;
  prnSlow = 0;
  rTol    = 0.000001;
  aTol    = 0.0;
  divgLim = 10.0;
  alpha   = 1.0;
  eta     = 0.0;
}


NonLinSIM::~NonLinSIM ()
{
  if (slowNodes.empty()) return;

  std::map<int,int>::const_iterator nit;
  std::cout <<"\n *** Here are the nodal points flagged with slow convergence"
            <<"\n     ======================================================="
            <<"\n     Node  Count  Patch  Coordinates\n";
  for (nit = slowNodes.begin(); nit != slowNodes.end(); ++nit)
  {
    Vec4 X = model.getNodeCoord(nit->first);
    std::cout << std::setw(9) << nit->first
              << std::setw(5) << nit->second
              << std::setw(7) << X.idx <<"    ";
    X.Vec3::print(std::cout) << std::endl;
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
  if (strcasecmp(elem->Value(),"nonlinearsolver"))
    return model.parse(elem);

  std::string rotUpdate;
  if (utl::getAttribute(elem,"rotation",rotUpdate,true) && !rotUpdate.empty())
    rotUpd = rotUpdate[0];

  const char* value = 0;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())
    if ((value = utl::getValue(child,"maxits")))
      maxit = atoi(value);
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
    }
    else if (!strcasecmp(child->Value(),"fromZero"))
      fromIni = true;

  return true;
}


void NonLinSIM::getTolerances (double& atol, double& rtol,
                               double& dtol, int& mxit) const
{
  atol = aTol;
  rtol = rTol;
  dtol = divgLim;
  mxit = maxit;
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
  // Update solution vectors between time steps
  for (int n = solution.size()-1; n > 0; n--)
    std::copy(solution[n-1].begin(),solution[n-1].end(),solution[n].begin());

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

  if (msgLevel >= 0 && model.getProcessAdm().getProcId() == 0)
  {
    std::streamsize oldPrec = 0;
    double digits = log10(param.time.t)-log10(param.time.dt);
    if (digits > 6.0) oldPrec = std::cout.precision(ceil(digits));
    std::cout <<"\n  step="<< param.step <<"  time="<< param.time.t;
    if (param.maxCFL > 0.0)
      std::cout <<"  CFL = "<< param.time.CFL << std::endl;
    else
      std::cout << std::endl;
    if (oldPrec > 0) std::cout.precision(oldPrec);
  }

  param.iter = 0;
  alpha = 1.0;
  if (fromIni) // Always solve from initial configuration
    solution.front().fill(0.0);

  if (!model.updateDirichlet(param.time.t,&solution.front()))
    return FAILURE;

  if (!model.setMode(mode))
    return FAILURE;

  bool poorConvg = false;
  bool newTangent = true;
  model.setQuadratureRule(opt.nGauss[0],true);
  if (!model.assembleSystem(param.time,solution,newTangent))
    return model.getProblem()->diverged() ? DIVERGED : FAILURE;

  if (iteNorm != NONE)
    if (!model.extractLoadVec(residual))
      return FAILURE;

  if (!model.solveSystem(linsol,msgLevel-1))
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

	if (param.iter == 1)
	  if (!model.updateDirichlet())
	    return FAILURE;

	if (param.iter > nupdat)
	{
	  newTangent = false;
	  model.setMode(RHS_ONLY);
	}
	else
	  model.setMode(mode);

	if (!model.assembleSystem(param.time,solution,newTangent,poorConvg))
	  return model.getProblem()->diverged() ? DIVERGED : FAILURE;

	if (!model.extractLoadVec(residual))
	  return FAILURE;

	if (!model.solveSystem(linsol,msgLevel-1))
	  return FAILURE;

	if (!this->lineSearch(param))
	  return FAILURE;

	poorConvg = false;
      }

  return DIVERGED;
}


bool NonLinSIM::solveLinearizedSystem (const TimeStep& param, double& norm)
{
  if (!model.assembleSystem(param.time,solution))
    return false;

  Vector inc;
  model.extractLoadVec(inc);
  norm = model.solutionNorms(inc);
  if (!model.solveSystem(inc))
    return false;

  solution.front().add(inc);
  return true;
}


bool NonLinSIM::lineSearch (TimeStep& param)
{
  if (eta <= 0.0) return true; // No line search

  double s0 = residual.dot(linsol);
  double smin = fabs(s0);
  double cmin = 1.0;
  double ck = 1.0;
  double cp = 1.0;
#ifdef SP_DEBUG
  std::cout <<"\t0: ck="<< ck <<" sk="<< s0 << std::endl;
#endif

  alpha = 1.0;
  for (int iter = 1; iter <= 10; iter++)
  {
    if (!this->updateConfiguration(param))
      return false;

    if (!model.assembleSystem(param.time,solution,false))
      return false;

    if (!model.extractLoadVec(residual))
      return false;

    double sk = residual.dot(linsol);
#ifdef SP_DEBUG
    std::cout <<"\t"<< iter <<": ck="<< ck <<" sk="<< sk << std::endl;
#endif
    if (fabs(sk) < eta*fabs(s0))
    {
      alpha = 0.0;
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
      return true;
    }

    alpha = ck - alpha;
  }

  alpha = cmin - alpha;
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

    if (refNopt == ALL || fabs(norm) > refNorm)
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

  if (msgLevel > 0 && model.getProcessAdm().getProcId() == 0)
  {
    // Print convergence history
    std::ios::fmtflags oldFlags = std::cout.flags(std::ios::scientific);
    std::streamsize oldPrec = std::cout.precision(3);
    std::cout <<"  iter="<< param.iter
              <<"  conv="<< fabs(norm)
              <<"  enen="<< enorm
              <<"  resn="<< resNorm
              <<"  incn="<< linsolNorm << std::endl;
    if (status == SLOW && prnSlow > 0)
    {
      // Find and print out the worst DOF(s) when detecting slow convergence
      std::map<std::pair<int,int>,RealArray> worstDOFs;
      model.getWorstDofs(linsol,residual,prnSlow,convTol*refNorm,worstDOFs);
      std::cout <<"  ** Slow convergence detected";
      if (worstDOFs.size() > 1)
        std::cout <<", here are the "<< worstDOFs.size() <<" worst DOFs:";
      else if (worstDOFs.size() == 1)
        std::cout <<", here is the worst DOF:";
      else
        std::cout <<".";
      std::map<std::pair<int,int>,RealArray>::const_iterator it;
      for (it = worstDOFs.begin(); it != worstDOFs.end(); it++)
      {
        std::cout <<"\n     Node "<< it->first.first
                  <<" local DOF "<< it->first.second;
        char nodeType = model.getNodeType(it->first.first);
        if (nodeType != ' ')
          std::cout <<" ("<< nodeType <<")";
        std::cout <<" :\tEnergy = "<< it->second[0]
                  <<"\tdu = "<< it->second[1]
                  <<"\tres = "<< it->second[2];
        std::map<int,int>::iterator nit = slowNodes.find(it->first.first);
        if (nit == slowNodes.end())
          slowNodes.insert(std::make_pair(it->first.first,1));
        else
          ++nit->second;
      }
      std::cout << std::endl;
    }
    std::cout.flags(oldFlags);
    std::cout.precision(oldPrec);
  }

  // Check for convergence or divergence
  if (fabs(norm) < convTol && (param.iter > 0 || refNopt == ALL))
    status = CONVERGED;
  else if (fabs(norm) <= fabs(prevNorm))
    nIncrease = 0;
  else if (++nIncrease > 2 || fabs(norm) > divgLim)
    status = DIVERGED;

  prevNorm = norm;
  return status;
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
      if (psol.size() != 6*model.getNoNodes(true))
      {
        std::cerr <<" *** NonLinSIM::updateConfiguration: Invalid dimension on"
                  <<" the displacement vector "<< psol.size()
                  <<" != "<< 6*model.getNoNodes(true) << std::endl;
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
