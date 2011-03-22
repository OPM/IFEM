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
    double dt;
    tInc.clear();
    tInc.reserve(5);
    std::istringstream cline(utl::readLine(is));

    cline >> startTime >> stopTime >> dt;
    if (cline.fail() || cline.bad()) return false;
    while (!cline.fail() && !cline.bad())
    {
      tInc.push_back(dt);
      cline >> dt;
    }
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


void NonLinSIM::init (SolvePrm& param, const Vector& initVal)
{
  param.startTime = startTime;
  param.stopTime = stopTime;
  param.tInc    = tInc;
  param.time.t  = startTime;
  param.time.dt = tInc.empty() ? stopTime-startTime : tInc.front();
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


bool NonLinSIM::advanceStep (SolvePrm& param)
{
  // Update solution vectors between time steps
  for (int n = solution.size()-1; n > 0; n--)
    solution[n] = solution[n-1];

  return param.increment();
}


bool NonLinSIM::solveStep (SolvePrm& param, SIM::SolutionMode mode)
{
  PROFILE1("NonLinSIM::solveStep");

  if (msgLevel >= 0 && myPid == 0)
    std::cout <<"\n  step="<< param.step
	      <<"  time="<< param.time.t << std::endl;

  param.iter = 0;
  param.alpha = 1.0;

  if (!model->updateDirichlet(param.time.t,&solution.front()))
    return false;

  if (!model->setMode(mode))
    return false;

  bool newTangent = true;
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

	param.time.first = false;
	model->setMode(SIM::RECOVERY);
	return this->solutionNorms(param.time);

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


bool NonLinSIM::solutionNorms (const TimeDomain& time, const char* compName)
{
  if (msgLevel < 0 || solution.empty()) return true;

  const int nsd = model->getNoSpaceDim();

  size_t iMax[nsd];
  double dMax[nsd];
  double normL2 = model->solutionNorms(solution.front(),dMax,iMax);

  if (!model->solutionNorms(time,solution,eNorm,gNorm))
    gNorm.clear();

  if (myPid == 0)
  {
    std::cout <<"  Primary solution summary: L2-norm            : "<< normL2;
    char D = 'X';
    for (int d = 0; d < nsd; d++, D++)
      std::cout <<"\n                            Max "<< D <<'-'<< compName
		<<" : "<< dMax[d] <<" node "<< iMax[d];
    if (gNorm.size() > 0)
    {
      std::cout <<"\n  Energy norm:    |u^h| = a(u^h,u^h)^0.5 : "<< gNorm(1);
      int oldPrec = std::cout.precision(10);
      std::cout <<"\t a(u^h,u^h) = "<< gNorm(1)*gNorm(1);
      std::cout.precision(oldPrec);
    }
    if (gNorm.size() > 1 && gNorm(2) != 0.0)
    {
      std::cout <<"\n  External energy: ((f,u^h)+(t,u^h))^0.5 : "<< gNorm(2);
      int oldPrec = std::cout.precision(10);
      std::cout <<"\t(f,u)+(t,u) = "<< gNorm(2)*gNorm(2);
      std::cout.precision(oldPrec);
    }
    std::cout << std::endl;
  }

  return true;
}


bool NonLinSIM::saveModel (char* fileName, int format, int* nViz)
{
  PROFILE1("NonLinSIM::saveModel");

  nBlock = 0; // initialize the result block counter

  // Write VTF-file with model geometry
  if (!model->writeGlv(fileName,nViz,format))
    return false;

  // Write Dirichlet boundary conditions
  return model->writeGlvBC(nViz,nBlock);
}


bool NonLinSIM::saveStep (int iStep, double time, int* nViz, bool psolOnly)
{
  PROFILE1("NonLinSIM::saveStep");

  // Negative iStep means we are saving the initial state only
  if (!model->setMode(iStep < 0 ? SIM::INIT : SIM::RECOVERY))
    return false;
  else if (iStep < 0)
    iStep = -iStep;

  // Write boundary tractions, if any
  if (!psolOnly)
    if (!model->writeGlvT(iStep,nBlock))
      return false;

  // Write residual force vector, but only when no extra visualization points
  if (!psolOnly && nViz[0] == 2 && nViz[1] <= 2 && nViz[2] <= 2)
    if (!model->writeGlvV(residual,"Residual forces",nViz,iStep,nBlock))
      return false;

  // Write solution fields
  if (!solution.empty())
    if (!model->writeGlvS(solution.front(),nViz,iStep,nBlock,psolOnly))
      return false;

  // Write element norms (only when no additional visualization points are used)
  if (!psolOnly && nViz[0] == 2 && nViz[1] <= 2 && nViz[2] <= 2)
    if (!model->writeGlvN(eNorm,iStep,nBlock))
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
