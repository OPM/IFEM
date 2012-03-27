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


NonLinSIM::NonLinSIM (SIMbase* sim) : model(sim), nBlock(0)
{
#if SP_DEBUG > 2
  msgLevel = 100; // prints the linear solution vector if size < 100
#else
  msgLevel = 1;   // prints the convergence history only
#endif

  // Default solution parameters
  maxit     = 20;
  nupdat    = 20;
  startTime = 0.0;
  stopTime  = 1.0;
  convTol   = 1.0e-6;
  divgLim   = 1.0;
  eta       = 0.0;
  iteNorm   = ENERGY;
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


void NonLinSIM::initSystem (int mType, size_t nGauss)
{
  model->initSystem(mType,1,1);
  model->setAssociatedRHS(0,0);
  model->setQuadratureRule(nGauss);
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
  if (initVal.size() == model->getNoDOFs())
    solution.front() = initVal;
}


bool NonLinSIM::advanceStep (SolvePrm& param, bool updateTime)
{
  // Update solution vectors between time steps
  for (int n = solution.size()-1; n > 0; n--)
    solution[n] = solution[n-1];

  return updateTime ? param.increment() : true;
}


bool NonLinSIM::solveStep (SolvePrm& param, SIM::SolutionMode mode,
			   const char* compName, bool energyNorm,
			   double zero_tolerance, std::streamsize outPrec)
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

  bool newTangent = true;
  model->setQuadratureRule(model->opt.nGauss[0]);
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

	model->setMode(SIM::RECOVERY);
	if (!this->solutionNorms(param.time,compName,energyNorm,
				 zero_tolerance,outPrec))
	  return false;

	param.time.first = false;
	return true;

      case DIVERGED:
	return false;

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

	if (!model->assembleSystem(param.time,solution,newTangent))
	  return false;

	if (!model->extractLoadVec(residual))
	  return false;

	if (!model->solveSystem(linsol,msgLevel-1))
	  return false;

	if (!this->lineSearch(param))
	  return false;
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

  double enorm, resNorm, linsolNorm;
  model->iterationNorms(linsol,residual,enorm,resNorm,linsolNorm);
  double norm = iteNorm == ENERGY ? enorm : resNorm;

  if (param.iter == 0)
  {
    if (norm > param.refNorm)
      param.refNorm = norm;
    norm = prevNorm = 1.0;
    nIncrease = 0;
  }
  else
    norm = fabs(norm)/param.refNorm;

  if (msgLevel > 0 && myPid == 0 && !solution.empty())
  {
    // Print convergence history
    std::ios::fmtflags oldFlags = std::cout.flags(std::ios::scientific);
    std::streamsize oldPrec = std::cout.precision(3);
    std::cout <<"  iter="<< param.iter
	      <<"  conv="<< norm
	      <<"  enen="<< enorm
	      <<"  resn="<< resNorm
	      <<"  incn="<< linsolNorm << std::endl;
    std::cout.flags(oldFlags);
    std::cout.precision(oldPrec);
  }

  // Check for convergence
  if (norm < param.convTol)
    return CONVERGED;

  // Check for divergence (increasing norm in three consequtive iterations)
  if (norm <= prevNorm)
    nIncrease = 0;
  else if (++nIncrease > 1 || norm > param.divgLim)
    return DIVERGED;

  prevNorm = norm;
  return NONE;
}


bool NonLinSIM::updateConfiguration (SolvePrm& param)
{
  if (solution.empty()) return false;

  if (solution.front().empty())
    solution.front() = linsol;
  else if (param.alpha != 0.0)
    solution.front().add(linsol,param.alpha);

  return true;
}


bool NonLinSIM::solutionNorms (const TimeDomain& time, const char* compName,
			       bool energyNorm, double zero_tolerance,
			       std::streamsize outPrec)
{
  if (msgLevel < 0 || solution.empty()) return true;

  const int nsd = model->getNoSpaceDim();

  size_t iMax[nsd];
  double dMax[nsd];
  double normL2 = model->solutionNorms(solution.front(),dMax,iMax);

  RealArray RF;
  bool haveRFs = model->getCurrentReactions(RF,solution.front());

  if (energyNorm)
  {
    model->setQuadratureRule(model->opt.nGauss[1]);
    if (!model->solutionNorms(time,solution,gNorm))
      gNorm.clear();
  }

  if (myPid == 0)
  {
    std::streamsize stdPrec = outPrec > 0 ? std::cout.precision(outPrec) : 0;
    double old_tol = utl::zero_print_tol;
    utl::zero_print_tol = zero_tolerance;
    std::cout <<"  Primary solution summary: L2-norm            : "
	      << utl::trunc(normL2);
    char D = 'X';
    for (int d = 0; d < nsd; d++, D++)
      if (utl::trunc(dMax[d]) != 0.0)
	std::cout <<"\n                            Max "<< D <<'-'<< compName
		  <<" : "<< dMax[d] <<" node "<< iMax[d];
    if (haveRFs)
    {
      std::cout <<"\n  Total reaction forces: Sum(R) =";
      for (size_t i = 1; i < RF.size(); i++)
        std::cout <<" "<< utl::trunc(RF[i]);
      if (utl::trunc(RF.front()) != 0.0)
	std::cout <<"\n  "<< compName <<"*reactions: (R,u) = "<< RF.front();
    }
    if (gNorm.size() > 0)
    {
      std::cout <<"\n  Energy norm:    |u^h| = a(u^h,u^h)^0.5 : "
		<< utl::trunc(gNorm(1));
      std::streamsize oldPrec = std::cout.precision(10);
      std::cout <<"\t a(u^h,u^h) = "<< utl::trunc(gNorm(1)*gNorm(1));
      std::cout.precision(oldPrec);
    }
    if (gNorm.size() > 1 && utl::trunc(gNorm(2)) != 0.0)
    {
      std::cout <<"\n  External energy: ((f,u^h)+(t,u^h))^0.5 : "<< gNorm(2);
      std::streamsize oldPrec = std::cout.precision(10);
      std::cout <<"\t(f,u)+(t,u) = "<< gNorm(2)*gNorm(2);
      std::cout.precision(oldPrec);
    }
    if (gNorm.size() > 2)
      std::cout <<"\n  Stress norm, L2: (sigma^h,sigma^h)^0.5 : "<< gNorm(3);
    if (gNorm.size() > 3)
      std::cout <<"\n  Pressure norm, L2:       (p^h,p^h)^0.5 : "<< gNorm(4)
		<<"\t(p^h = trace(sigma^h)/3)";
    if (gNorm.size() > 4)
      std::cout <<"\n  Deviatoric stress norm:  (s^d,s^d)^0.5 : "<< gNorm(5)
		<<"\t(s^d = sigma^h - p^h*I)";
    if (gNorm.size() > 5)
      std::cout <<"\n  Stress norm, von Mises: vm(sigma^h)    : "<< gNorm(6);
    std::cout << std::endl;

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


bool NonLinSIM::project ()
{
  return model->project(solution.front());
}
