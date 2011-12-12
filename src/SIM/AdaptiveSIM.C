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
#include "ASMu2D.h"
#endif
#include "ASMbase.h"
#include "SIMbase.h"
#include "SIMenums.h"
#include "Utilities.h"
#include <sstream>
#include <cstdio>


AdaptiveSIM::AdaptiveSIM (SIMbase* sim) : model(sim)
{
  // Default grid adaptation parameters
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
    options.push_back(beta);

    double itmp;
    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      maxStep = itmp;

    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      maxDOFs = itmp;

    cline >> itmp; // read knotline multiplicity
    if (!cline.fail() && !cline.bad())
      options.push_back(itmp);

    cline >> itmp; // read refinement scheme (0=refine all, 1=minimum span, 2=isotropic)
    if (!cline.fail() && !cline.bad())
      options.push_back(itmp);

    cline >> itmp; // read symmetry (0=none, <n> always requests a multiplum of <n> elements for refinement)
    if (!cline.fail() && !cline.bad())
      options.push_back(itmp);
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
  else
    // Output the initial grid to eps-file
    model->refine(std::vector<int>(),options,"mesh_001.eps");

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
  double RMS_norm = 0;
  double avg_norm = 0;
  double min_err = 10000000;
  double max_err  = 0;
  for(size_t i=1; i<=eNorm.cols(); i++) {
    avg_norm += eNorm(4,i);
    min_err = (min_err < eNorm(4,i)) ? min_err : eNorm(4,i);
    max_err = (max_err > eNorm(4,i)) ? max_err : eNorm(4,i);
  }
  avg_norm  /= eNorm.cols();
  for(size_t i=1; i<=eNorm.cols(); i++)
    RMS_norm += pow(eNorm(4,i)-avg_norm,2);
  RMS_norm = sqrt(RMS_norm/eNorm.cols())/avg_norm;
  gNorm.push_back(RMS_norm);
  gNorm.push_back(min_err);
  gNorm.push_back(max_err);
  gNorm.push_back(avg_norm);
  printNorms(gNorm,std::cout);

  // Check if further refinement is required
  if (iStep > maxStep || model->getNoDOFs() > (size_t)maxDOFs) return false;
  if (gNorm.size() > 3 && 100.0*gNorm(4) < errTol*gNorm(3)) return false;

  std::vector<int> toBeRefined;
  // Vector errors;
  std::vector<IndexDouble> errors;

  if(options.size() > 2 && options[2]==3)
  {
    // sum up the total function error over all supported elements for that function
    ASMbase *patch = model->getFEModel()[0];
    IntMat::const_iterator eit;
    IntVec::const_iterator nit;
    for(size_t i=0; i<patch->getNoNodes(); i++)
      errors.push_back(IndexDouble(i, 0));
    int i = 1;
    for(eit=patch->begin_elm(); eit<patch->end_elm(); eit++, i++)
      for(nit=eit->begin(); nit<eit->end(); nit++)
        errors[*nit].v += eNorm(4,i);
  }
  else
  {
    for(size_t i=1; i<=eNorm.cols(); i++)
      errors.push_back(IndexDouble(i-1, eNorm(4,i)));
  }

  // Find the list of elements/functions to refine (the beta % with the highest error)
  size_t ipivot = ceil(errors.size()*beta/100.0);
  // make ipivot a multiplum of options[3] in case of symmetric problems
  if(options.size()>3 && options[3]!=0)
    ipivot += (options[3]-ipivot%options[3]);

  if (ipivot < 1 || ipivot > errors.size()) return false;
  std::sort(errors.begin(),errors.end(),
	    std::greater<IndexDouble>());

/*
  double pivot = errors[ipivot].v;
  std::cout <<"\nRefining "<< ipivot <<" elements with errors in range ["
	    << pivot <<","<< errors.front().v <<"]"<< std::endl;
*/

  std::cout <<"\nRefining up to "<< beta <<" \% new DOFs " << std::endl;

  toBeRefined.reserve(ipivot);
  for (size_t i = 0; i < errors.size(); i++)
    toBeRefined.push_back(errors[i].i);

  // Now refine the mesh
  char fname[13];
  sprintf(fname,"mesh_%03d.eps",iStep);
  return model->refine(toBeRefined,options,fname);
}


void AdaptiveSIM::printNorms (const Vector& norms, std::ostream& os)
{
  os <<    "Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< norms(1)
     <<  "\nExternal energy ((h,u^h)+(t,u^h)^0.5 : "<< norms(2);
  if (norms.size() > 2)
    os <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< norms(3);
  if (norms.size() > 3)
    os <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< norms(4)
       <<"\nExact relative error (%) : "<< 100.0*norms(4)/norms(3);
  if (norms.size() > 4)
    os <<"\nRoot mean square (RMS) of est.err    : "<< norms(5);
  if (norms.size() > 5) {
    os <<"\nMin element error                    : "<< norms(6);
    os <<"\nMax element error                    : "<< norms(7);
    os <<"\nAvg element error                    : "<< norms(8);
  }
  os << std::endl;
}


bool AdaptiveSIM::writeGlv (const char* infile, int format, const int* nViz,
			    int iStep, int& nBlock)
{
  if (format < 0) return true;

  // Write VTF-file with model geometry
  if (!model->writeGlvG(nViz, nBlock, iStep == 1 ? infile : 0, format))
    return false;

  // Write boundary tractions, if any
  if (!model->writeGlvT(iStep,nBlock))
    return false;

  // Write Dirichlet boundary conditions
  if (!model->writeGlvBC(nViz,nBlock,iStep))
    return false;

  // Write solution fields
  if (!model->writeGlvS(linsol,nViz,iStep,nBlock))
    return false;

  // Write element norms
  if (!model->writeGlvN(eNorm,iStep,nBlock))
    return false;

  // Write state information
  return model->writeGlvStep(iStep,iStep,1);
}
