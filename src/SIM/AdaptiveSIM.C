// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.h
//!
//! \date Sep 22 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution driver for isogeometric FEM simulators.
//!
//==============================================================================

#include "AdaptiveSIM.h"
#include "ASMunstruct.h"
#include "SIMbase.h"
#include "SIMenums.h"
#include "Utilities.h"
#include <sstream>


AdaptiveSIM::AdaptiveSIM (SIMbase* sim) : model(sim)
{
  // Default adaptation parameters
  nStep = 10;
  stopTol = 1.0;
  beta = 25.0;
}


AdaptiveSIM::~AdaptiveSIM ()
{
  if (model) delete model;
}


bool AdaptiveSIM::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ADAPTIVE",8))
  {
    std::istringstream cline(utl::readLine(is));
    cline >> nStep >> stopTol >> beta;
    if (cline.fail() || cline.bad()) return false;
  }
  else
    return model->parse(keyWord,is);

  return true;
}


bool AdaptiveSIM::solveStep (const char* inputfile, SystemMatrix::Type solver,
			     int iStep)
{
  if (iStep > 1)
  {
    // Re-generate the FE model after the refinement
    ASMunstruct::resetNumbering();
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
  return model->solutionNorms(Vectors(1,linsol),eNorm,gNorm);
}


static bool larger (double a, double b) { return a > b; }


bool AdaptiveSIM::adaptMesh (int iStep)
{
  printNorms(gNorm,std::cout);

  double eta = gNorm.size() > 3 ? 100.0*gNorm(4)/gNorm(3) : 0.0;
  if (eta < stopTol || iStep > nStep) return false;

  // Find the list of elements to refine (the beta % with the highest error)
  std::vector<int> elements;
  Vector errors(eNorm.getRow(4));
  size_t ipivot = ceil(errors.size()*beta/100.0);
  if (ipivot < 1 || ipivot > errors.size()) return false;
  std::partial_sort(errors.begin(),errors.begin()+ipivot,errors.end(),larger);

  double pivot = errors(ipivot);
  std::cout <<"\nRefining "<< ipivot <<" elements with errors in range ["
	    << pivot <<","<< errors.front() <<"]"<< std::endl;

  elements.reserve(ipivot);
  for (size_t e = 1; e <= errors.size(); e++)
    if (eNorm(4,e) >= pivot)
      elements.push_back(e-1);

  char fname[10] = "mesh_.eps";
  fname[4] = '0' + iStep;
  return model->refine(elements,fname);
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
  std::cout << std::endl;
}
