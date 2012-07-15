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
#include "SIMbase.h"
#include "Profiler.h"
#include "Utilities.h"
#include "tinyxml.h"
#include <sstream>


NonLinSIM::NonLinSIM (SIMbase* sim, CNORM n) : model(sim), iteNorm(n), nBlock(0)
{
#ifndef SP_DEBUG
  msgLevel = 1;   // prints the convergence history only
#elif SP_DEBUG > 2
  msgLevel = 100; // prints the linear solution vector if its size is < 100
#endif

  // Default solution parameters
  maxit     = 20;
  nupdat    = 20;
  startTime = 0.0;
  stopTime  = 1.0;
  convTol   = 1.0e-6;
  divgLim   = 1.0;
  eta       = 0.0;
}


NonLinSIM::~NonLinSIM ()
{
  if (model) delete model;
}


bool NonLinSIM::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"TIME_STEPPING",13))
  {
    int nstep = atoi(keyWord+13);
    if (nstep < 1) nstep = 1;

    double dt;
    steps.resize(nstep);
    for (int i = 0; i < nstep; i++)
    {
      std::istringstream cline(utl::readLine(is));
      if (i == 0) cline >> startTime;
      cline >> steps[i].second >> dt;
      if (cline.fail() || cline.bad()) return false;
      if (dt > 1.0 && ceil(dt) == dt)
      {
	// The number of time steps are specified
	dt = (steps[i].second - (i == 0 ? startTime : steps[i-1].second))/dt;
	steps[i].first.push_back(dt);
      }
      else while (!cline.fail() && !cline.bad())
      {
	// The time step size(s) is/are specified
	steps[i].first.push_back(dt);
	cline >> dt;
      }
    }
    stopTime = steps.back().second;
  }
  else if (!strncasecmp(keyWord,"NONLINEAR_SOLVER",16))
  {
    std::istringstream cline(utl::readLine(is));
    cline >> maxit >> convTol;
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
    return model->parse(keyWord,is);

  return true;
}


bool NonLinSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"nonlinearsolver"))
    return model->parse(elem);

  const char* value = 0;
  const TiXmlElement* child = elem->FirstChildElement();
  for (; child; child = child->NextSiblingElement())

    if (!strcasecmp(child->Value(),"timestepping")) {
      steps.clear();
      const TiXmlElement* step = child->FirstChildElement("step");
      for (; step; step = step->NextSiblingElement()) {
        double start = 0.0, end = 0.0, dt = 0.0;
        std::pair<std::vector<double>,double> timeStep;
	utl::getAttribute(step,"start",start);
	utl::getAttribute(step,"end",end);
        if (steps.empty()) startTime = start;
        if (step->FirstChild()) {
          std::istringstream cline(step->FirstChild()->Value());
          cline >> dt;
          if (dt > 1.0 && ceil(dt) == dt) {
            // The number of steps are specified
            dt = (end - (steps.empty() ? start : steps.back().second))/dt;
            timeStep.first.push_back(dt);
          }
          else while (!cline.fail() && !cline.bad()) {
            // The time step size(s) is/are specified
            timeStep.first.push_back(dt);
            cline >> dt;
          }
          timeStep.second = end;
          steps.push_back(timeStep);
        }
      }
      stopTime = steps.back().second;
    }

    else if ((value = utl::getValue(child,"maxits")))
      maxit = atoi(value);

    else if ((value = utl::getValue(child,"nupdate")))
      nupdat = atoi(value);

    else if ((value = utl::getValue(child,"rtol")))
      convTol = atof(value);

    else if ((value = utl::getValue(child,"dtol")))
      divgLim = atof(value);

    else if ((value = utl::getValue(child,"eta")))
      eta = atof(value);

  return true;
}


void NonLinSIM::init (SolvePrm& param, const RealArray& initVal)
{
  param.initTime(startTime,stopTime,steps);
  param.maxit   = maxit;
  param.nupdat  = nupdat;
  param.convTol = convTol;
  param.divgLim = divgLim;
  param.eta     = eta;

  size_t nSols = model->getNoSolutions();
  if (nSols < 2) nSols = 2;
  solution.resize(nSols);
  for (size_t n = 0; n < nSols; n++)
    solution[n].resize(model->getNoDOFs(),true);

  // Set initial conditions for time-dependent problems
  this->setInitialGuess(initVal);
}


void NonLinSIM::setInitialGuess (const RealArray& value)
{
  if (value.empty() || solution.empty())
    return;

  size_t ndim = solution.front().size();
  if (value.size() > ndim)
    std::copy(value.begin(),value.begin()+ndim,solution.front().begin());
  else
    std::copy(value.begin(),value.end(),solution.front().begin());
}


bool NonLinSIM::advanceStep (SolvePrm& param, bool updateTime)
{
  // Update solution vectors between time steps
  for (int n = solution.size()-1; n > 0; n--)
    solution[n] = solution[n-1];

  return updateTime ? param.increment() : true;
}


bool NonLinSIM::solveStep (SolvePrm& param, SIM::SolutionMode mode,
			   bool energyNorm, double zero_tolerance,
			   std::streamsize outPrec)
{
  PROFILE1("NonLinSIM::solveStep");

  if (msgLevel >= 0 && myPid == 0)
  {
    std::streamsize oldPrec = 0;
    double digits = log10(param.time.t)-log10(param.time.dt);
    if (digits > 6.0) oldPrec = std::cout.precision(ceil(digits));
    std::cout <<"\n  step="<< param.step
	      <<"  time="<< param.time.t << std::endl;
    if (oldPrec > 0) std::cout.precision(oldPrec);
  }

  param.iter = 0;
  param.alpha = 1.0;

  if (!model->updateDirichlet(param.time.t,&solution.front()))
    return false;

  if (!model->setMode(mode))
    return false;

  bool poorConvg = false;
  bool newTangent = true;
  model->setQuadratureRule(model->opt.nGauss[0],true);
  if (!model->assembleSystem(param.time,solution,newTangent))
    return false;

  if (!model->extractLoadVec(residual))
    return false;

  if (!model->solveSystem(linsol,msgLevel-1))
    return false;

  while (param.iter <= param.maxit)
    switch (this->checkConvergence(param))
      {
      case CONVERGED:
	if (!this->updateConfiguration(param))
	  return false;

	if (!this->solutionNorms(param.time,energyNorm,zero_tolerance,outPrec))
	  return false;

	param.time.first = false;
	return true;

      case DIVERGED:
	return false;

      case SLOW:
	poorConvg = true;
      default:
	param.iter++;
	if (!this->updateConfiguration(param))
	  return false;

	if (param.iter == 1)
	  if (!model->updateDirichlet())
	    return false;

	if (param.iter > param.nupdat)
	{
	  newTangent = false;
	  model->setMode(SIM::RHS_ONLY);
	}
	else
	  model->setMode(mode);

	if (!model->assembleSystem(param.time,solution,newTangent,poorConvg))
	  return false;

	if (!model->extractLoadVec(residual))
	  return false;

	if (!model->solveSystem(linsol,msgLevel-1))
	  return false;

	if (!this->lineSearch(param))
	  return false;

	poorConvg = false;
      }

  return false;
}


bool NonLinSIM::lineSearch (SolvePrm& param)
{
  if (param.eta <= 0.0) return true; // No line search

  double s0 = residual.dot(linsol);
  double smin = fabs(s0);
  double cmin = 1.0;
  double ck = 1.0;
  double cp = 1.0;
  if (myPid == 0)
    std::cout <<"\t0: ck="<< ck <<" sk="<< s0 << std::endl;

  param.alpha = 1.0;
  for (int iter = 1; iter <= 10; iter++)
  {
    if (!this->updateConfiguration(param))
      return false;

    if (!model->assembleSystem(param.time,solution,false))
      return false;

    if (!model->extractLoadVec(residual))
      return false;

    double sk = residual.dot(linsol);
    if (myPid == 0)
      std::cout <<"\t"<< iter <<": ck="<< ck <<" sk="<< sk << std::endl;

    if (fabs(sk) < param.eta*fabs(s0))
    {
      param.alpha = 0.0;
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

    if (fabs(ck-cp) < 0.5*param.eta*fabs(ck+cp))
    {
      param.alpha = 0.0;
      return true;
    }

    param.alpha = ck - param.alpha;
  }

  param.alpha = cmin - param.alpha;
  return true;
}


NonLinSIM::ConvStatus NonLinSIM::checkConvergence (SolvePrm& param)
{
  static double prevNorm  = 0.0;
  static int    nIncrease = 0;

  ConvStatus status = OK;
  double enorm, resNorm, linsolNorm;
  model->iterationNorms(linsol,residual,enorm,resNorm,linsolNorm);
  double norm = iteNorm == ENERGY ? enorm : resNorm;

  if (param.iter == 0)
  {
    if (fabs(norm) > param.refNorm)
      param.refNorm = fabs(norm);
    norm = prevNorm = 1.0;
    nIncrease = 0;
  }
  else
    norm /= param.refNorm;

  // Check for slow convergence
  if (param.iter > 1 && prevNorm > 0.0 && fabs(norm) > prevNorm*0.1)
    status = SLOW;

  if (msgLevel > 0 && myPid == 0 && !solution.empty())
  {
    if (status == SLOW)
    {
      // Find and print out the worst DOFs
      std::map<std::pair<int,int>,RealArray> worstDOFs;
      model->getWorstDofs(linsol,residual,param.convTol,worstDOFs);
      std::cout <<"  ** Slow convergence detected, here are the "
                << worstDOFs.size() <<" worst DOFs:";
      std::map<std::pair<int,int>,RealArray>::const_iterator it;
      for (it = worstDOFs.begin(); it != worstDOFs.end(); it++)
        std::cout <<"\n     Node "<< it->first.first
                  <<" Local DOF "<< it->first.second
                  <<" ("<< model->getNodeType(it->first.first)
                  <<") :\tEnergy = "<< it->second[0]
                  <<"\tdu = "<< it->second[1]
                  <<" res = "<< it->second[2];
      std::cout << std::endl;
    }

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
  if (fabs(norm) < param.convTol)
    status = CONVERGED;
  else if (fabs(norm) <= fabs(prevNorm))
    nIncrease = 0;
  else if (++nIncrease > 2 || fabs(norm) > param.divgLim)
    status = DIVERGED;

  prevNorm = norm;
  return status;
}


bool NonLinSIM::updateConfiguration (SolvePrm& param)
{
  if (solution.empty()) return false;

  if (solution.front().empty())
    solution.front() = linsol;
  else if (param.alpha != 0.0)
    solution.front().add(linsol,param.alpha);

  return model->updateConfiguration(solution.front());
}


bool NonLinSIM::solutionNorms (const TimeDomain& time, bool,
			       double zero_tolerance, std::streamsize outPrec)
{
  if (msgLevel < 0 || solution.empty()) return true;

  const size_t nf = model->getNoFields(1);

  size_t iMax[nf];
  double dMax[nf];
  double normL2 = model->solutionNorms(solution.front(),dMax,iMax);

  if (myPid == 0)
  {
    std::streamsize stdPrec = outPrec > 0 ? std::cout.precision(outPrec) : 0;
    double old_tol = utl::zero_print_tol;
    utl::zero_print_tol = zero_tolerance;

    std::cout <<"  Primary solution summary: L2-norm         : "
	      << utl::trunc(normL2);
    if (nf == 1 && utl::trunc(dMax[0]) != 0.0)
      std::cout <<"\n                            Max value       : "
	        << dMax[0] <<" node "<< iMax[0];
    else for (unsigned char d = 0; d < nf; d++)
      if (utl::trunc(dMax[d]) != 0.0)
        std::cout <<"\n                            Max "<< 'X'+d
                  <<"-component : "<< dMax[d] <<" node "<< iMax[d];

    utl::zero_print_tol = old_tol;
    if (stdPrec > 0) std::cout.precision(stdPrec);
  }

  return true;
}


bool NonLinSIM::saveModel (char* fileName)
{
  PROFILE1("NonLinSIM::saveModel");

  nBlock = 0; // initialize the result block counter

  // Write VTF-file with model geometry
  if (!model->writeGlvG(nBlock,fileName))
    return false;

  // Write Dirichlet boundary conditions
  return model->writeGlvBC(nBlock);
}


bool NonLinSIM::saveStep (int iStep, double time, int& iBlock,
			  bool psolOnly, const char* vecName)
{
  PROFILE1("NonLinSIM::saveStep");

  // Negative iStep means we are saving the initial state only
  if (!model->setMode(iStep < 0 ? SIM::INIT : SIM::RECOVERY))
    return false;
  else if (iStep < 0)
    iStep = -iStep;

  // Write boundary tractions, if any
  if (!psolOnly)
    if (!model->writeGlvT(iStep,iBlock))
      return false;

  // Write residual force vector, but only when no extra visualization points
  if (!psolOnly && model->opt.nViz[0] == 2 &&
      model->opt.nViz[1] <= 2 && model->opt.nViz[2] <= 2)
    if (!model->writeGlvV(residual,"Residual forces",iStep,iBlock))
      return false;

  // Write solution fields
  if (!solution.empty())
    if (!model->writeGlvS(solution.front(),iStep,iBlock,time,psolOnly,vecName))
      return false;

  // Write element norms
  if (!psolOnly)
    if (!model->writeGlvN(eNorm,iStep,iBlock))
      return false;

  // Write time/load step information
  return model->writeGlvStep(iStep,time);
}


void NonLinSIM::dumpStep (int iStep, double time, std::ostream& os,
			  bool withID) const
{
  if (withID)
    os <<"\n\n# Dump of primary solution at Step "<< iStep
       <<", Time="<< time <<"\n";

  model->dumpPrimSol(solution.front(),os,withID);
}


bool NonLinSIM::dumpResults (double time, std::ostream& os, int precision) const
{
  return model->dumpResults(solution.front(),time,os,true,precision);
}
