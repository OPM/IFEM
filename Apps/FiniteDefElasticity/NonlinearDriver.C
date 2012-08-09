// $Id$
//==============================================================================
//!
//! \file NonlinearDriver.C
//!
//! \date Jul 15 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Nonlinear driver for isogeometric finite deformation FEM analysis.
//!
//==============================================================================

#include "NonlinearDriver.h"
#include "SIMbase.h"
#include "Elasticity.h"
#include "DataExporter.h"
#include "tinyxml.h"


bool NonlinearDriver::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"TIME_STEPPING",13))
    return params.parse(keyWord,is);

  return this->NonLinSIM::parse(keyWord,is);
}


bool NonlinearDriver::parse (const TiXmlElement* elem)
{
  if (!strcasecmp(elem->Value(),"nonlinearsolver"))
  {
    const TiXmlElement* child = elem->FirstChildElement();
    for (; child; child = child->NextSiblingElement())
      params.parse(child);
  }

  return this->NonLinSIM::parse(elem);
}


bool NonlinearDriver::solutionNorms (const TimeDomain& time, bool energyNorm,
                                     double zero_tol, std::streamsize outPrec)
{
  if (msgLevel < 0 || solution.empty()) return true;

  const size_t nsd = model->getNoSpaceDim();

  size_t iMax[nsd];
  double dMax[nsd];
  double normL2 = model->solutionNorms(solution.front(),dMax,iMax);

  RealArray RF;
  bool haveReac = model->getCurrentReactions(RF,solution.front());

  Vectors gNorm;
  if (energyNorm)
  {
    model->setMode(SIM::RECOVERY);
    model->setQuadratureRule(model->opt.nGauss[1]);
    if (!model->solutionNorms(time,solution,gNorm))
      gNorm.clear();
  }

  if (myPid > 0) return true;

  std::streamsize stdPrec = outPrec > 0 ? std::cout.precision(outPrec) : 0;
  double old_tol = utl::zero_print_tol;
  utl::zero_print_tol = zero_tol;

  std::cout <<"  Primary solution summary: L2-norm            : "
            << utl::trunc(normL2);

  for (unsigned char d = 0; d < nsd; d++)
    if (utl::trunc(dMax[d]) != 0.0)
      std::cout <<"\n                            Max "<< char('X'+d)
                <<"-displacement : "<< dMax[d] <<" node "<< iMax[d];

  if (haveReac)
  {
    std::cout <<"\n  Total reaction forces: Sum(R) =";
    for (size_t i = 1; i < RF.size(); i++)
      std::cout <<" "<< utl::trunc(RF[i]);
    if (utl::trunc(RF.front()) != 0.0)
      std::cout <<"\n  displacement*reactions: (R,u) = "<< RF.front();
  }

  if (!gNorm.empty())
  {
    const Vector& norm = gNorm.front();
    if (norm.size() > 0)
    {
      std::cout <<"\n  Energy norm:    |u^h| = a(u^h,u^h)^0.5 : "
                << utl::trunc(norm(1));
      std::streamsize oldPrec = std::cout.precision(10);
      std::cout <<"\t a(u^h,u^h) = "<< utl::trunc(norm(1)*norm(1));
      std::cout.precision(oldPrec);
    }
    if (norm.size() > 1 && utl::trunc(norm(2)) != 0.0)
    {
      std::cout <<"\n  External energy: ((f,u^h)+(t,u^h))^0.5 : "<< norm(2);
      std::streamsize oldPrec = std::cout.precision(10);
      std::cout <<"\t(f,u)+(t,u) = "<< norm(2)*norm(2);
      std::cout.precision(oldPrec);
    }
    if (norm.size() > 2)
      std::cout <<"\n  Stress norm, L2: (sigma^h,sigma^h)^0.5 : "<< norm(3);
    if (norm.size() > 3)
      std::cout <<"\n  Pressure norm, L2:       (p^h,p^h)^0.5 : "<< norm(4)
                <<"\t(p^h = trace(sigma^h)/3)";
    if (norm.size() > 4)
      std::cout <<"\n  Deviatoric stress norm:  (s^d,s^d)^0.5 : "<< norm(5)
                <<"\t(s^d = sigma^h - p^h*I)";
    if (norm.size() > 5)
      std::cout <<"\n  Stress norm, von Mises: vm(sigma^h)    : "<< norm(6);
    std::cout << std::endl;
  }

  utl::zero_print_tol = old_tol;
  if (stdPrec > 0) std::cout.precision(stdPrec);
  return true;
}


/*!
  This method controls the load incrementation loop of the finite deformation
  simulation. It uses the automatic increment size adjustment of the TimeStep
  class and supports iteration cut-back in case of divergence.
*/

int NonlinearDriver::solveProblem (bool skip2nd, bool energyNorm,
                                   DataExporter* writer,
                                   std::ostream* oss, double dtDump,
                                   double zero_tol, std::streamsize outPrec)
{
  std::streamsize normPrec = outPrec > 3 ? outPrec : 0;

  if (dtDump <= 0.0) dtDump = params.stopTime + 1.0;
  double nextDump = params.time.t + dtDump;
  double nextSave = params.time.t + model->opt.dtSave;
  const Elasticity* elp = dynamic_cast<const Elasticity*>(model->getProblem());

  int iStep = 0; // Save initial state to VTF
  if (model->opt.format >= 0 && params.multiSteps() && params.time.dt > 0.0)
    if (!this->saveStep(-(++iStep),params.time.t,skip2nd))
      return 4;

  // Initialize the linear solver
  model->initSystem(model->opt.solver,1,1);

  // Invoke the time-step loop
  NonLinSIM::ConvStatus stat = NonLinSIM::OK;
  while (this->advanceStep(params))
  {
    do
    {
      if (stat == NonLinSIM::DIVERGED)
      {
        // Try cut-back with a smaller time step when diverging
        if (!params.cutback()) break;

	std::copy(solution[1].begin(),solution[1].end(),solution[0].begin());
	model->updateConfiguration(solution.front());
	refNorm = 1.0; // Reset the reference norm
      }

      // Solve the nonlinear FE problem at this load step
      stat = this->solveStep(params,SIM::STATIC,energyNorm,zero_tol,normPrec);
    }
    while (stat == NonLinSIM::DIVERGED);

    if (stat != NonLinSIM::CONVERGED)
      return 5;

    // Print solution components at the user-defined points
    this->dumpResults(params.time.t,std::cout,outPrec);

    if (params.hasReached(nextDump))
    {
      // Dump primary solution for inspection or external processing
      if (oss)
        this->dumpStep(params.step,params.time.t,*oss,false);
      else
        this->dumpStep(params.step,params.time.t,std::cout);

      nextDump = params.time.t + dtDump;
    }

    if (params.hasReached(nextSave))
    {
      // Save solution variables to VTF for visualization
      if (model->opt.format >= 0)
      {
        if (!this->saveStep(++iStep,params.time.t,skip2nd))
          return 6;

	// Print out the maximum values
	if (elp)
	{
	  elp->printMaxVals(std::cout,outPrec,7);
	  elp->printMaxVals(std::cout,outPrec,8);
	}
      }

      // Save solution variables to HDF5
      if (writer)
        if (!writer->dumpTimeLevel())
          return 7;

      nextSave = params.time.t + model->opt.dtSave;
      if (nextSave > params.stopTime)
        nextSave = params.stopTime; // Always save the final step
    }
  }

  return 0;
}
