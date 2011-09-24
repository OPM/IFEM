// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.h
//!
//! \date Sep 22 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution driver for linear isogeometric FEM simulators.
//!
//==============================================================================

#include "AdaptiveSIM.h"
#ifdef HAS_LRSPLINE
#include "ASMunstruct.h"
#endif
#include "SIMbase.h"
#include "SIMenums.h"
#include "Utilities.h"
#include <sstream>
#include <cstdio>


AdaptiveSIM::AdaptiveSIM (SIMbase* sim) : model(sim)
{
  // Default adaptation parameters
  beta = 10.0;
  errTol = 1.0;
  maxStep = 10;
  maxDOFs = 1000000;
}


AdaptiveSIM::~AdaptiveSIM ()
{
  if (model) delete model;
}


bool AdaptiveSIM::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ADAPTIVE",8))
  {
    options.clear();
    std::istringstream cline(utl::readLine(is));

    cline >> beta >> errTol;
    if (cline.fail() || cline.bad()) return false;

    double itmp;
    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      maxStep = itmp;

    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      maxDOFs = itmp;

    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      options.push_back(itmp);

    std::string ctmp;
    cline >> ctmp;
    if (!cline.fail() && !cline.bad())
      options.push_back(1);
  }
  else
    return model->parse(keyWord,is);

  return true;
}


bool AdaptiveSIM::solveStep (const char* inputfile, SystemMatrix::Type solver,
			     int iStep)
{
  std::cout <<"\nAdaptive step "<< iStep << std::endl;
  if (iStep > 1)
  {
    // Re-generate the FE model after the refinement
    model->clearProperties();
#ifdef HAS_LRSPLINE
    ASMunstruct::resetNumbering();
#endif
    if (!model->read(inputfile) || !model->preprocess())
      return false;
  }

  // Assemble the linear FE equation system
  model->setMode(SIM::STATIC,true);
  model->initSystem(iStep == 1 ? SystemMatrix::DENSE : solver, 1,1);
  model->setAssociatedRHS(0,0);
  if (!model->assembleSystem())
    return false;

  // Solve the linear system of equations
  if (!model->solveSystem(linsol,1))
    return false;

  // Evaluate solution norms
  model->setMode(SIM::RECOVERY);
  return (model->solutionNorms(Vectors(1,linsol),eNorm,gNorm) &&
	  model->dumpResults(linsol,0.0,std::cout,true,6));
}


bool AdaptiveSIM::adaptMesh (int iStep)
{
  printNorms(gNorm,std::cout);

  // Check if further refinement is required
  if (iStep > maxStep || model->getNoDOFs() > (size_t)maxDOFs) return false;
  if (gNorm.size() > 3 && 100.0*gNorm(4) < errTol*gNorm(3)) return false;

  // Find the list of elements to refine (the beta % with the highest error)
  std::vector<int> elements;
  Vector errors(eNorm.getRow(4));
  size_t ipivot = ceil(errors.size()*beta/100.0);
  if (ipivot < 1 || ipivot > errors.size()) return false;
  std::partial_sort(errors.begin(),errors.begin()+ipivot,
		    errors.end(),std::greater<double>());

  double pivot = errors(ipivot);
  std::cout <<"\nRefining "<< ipivot <<" elements with errors in range ["
	    << pivot <<","<< errors.front() <<"]"<< std::endl;

  elements.reserve(ipivot);
  for (size_t e = 1; e <= errors.size(); e++)
    if (eNorm(4,e) >= pivot)
      elements.push_back(e-1);

  // Now refine the mesh
  char fname[12];
  sprintf(fname,"mesh_%02d.eps",iStep);
  return model->refine(elements,options,fname);
}


void AdaptiveSIM::printNorms (const Vector& norms, std::ostream& os)
{
  os <<"Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< norms(1)
     <<"\nExternal energy ((h,u^h)+(t,u^h)^0.5 : "<< norms(2);
  if (norms.size() > 2)
    os <<"\nExact norm  |u|   = a(u,u)^0.5     : "<< norms(3);
  if (norms.size() > 3)
    os <<"\nExact error a(e,e)^0.5, e=u-u^h    : "<< norms(4)
       <<"\nExact relative error (%) : "<< 100.0*norms(4)/norms(3);
  os << std::endl;
}
