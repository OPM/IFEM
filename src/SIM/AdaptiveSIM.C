// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.C
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
#else
#include "ASMbase.h"
#endif
#include "SIMbase.h"
#include "SIMenums.h"
#include "Utilities.h"
#include <sstream>
#include <cstdio>

#include "tinyxml.h"


AdaptiveSIM::AdaptiveSIM (SIMbase* sim) : model(sim)
{
  // Default grid adaptation parameters
  storeMesh = false;
  beta      = 10.0;
  errTol    = 1.0;
  maxStep   = 10;
  maxDOFs   = 1000000;
  scheme    = 0; // fullspan
  symmetry  = 1; // no symmetry
  knot_mult = 1; // maximum regularity (continuity)
}


AdaptiveSIM::~AdaptiveSIM ()
{
  if (model) delete model;
}

bool AdaptiveSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"adaptive"))
    return model->parse(elem);

  const TiXmlElement* child = elem->FirstChildElement();
  std::string str;
  while (child) {
    if (!strcasecmp(child->Value(), "maxstep")) {
      utl::getAttribute(child, "value", maxStep);
    } else if (!strcasecmp(child->Value(), "beta")) {
      utl::getAttribute(child, "value", beta);
    } else if (!strcasecmp(child->Value(), "maxdof")) {
      utl::getAttribute(child, "value", maxDOFs);
    } else if (!strcasecmp(child->Value(), "errtol")) {
      utl::getAttribute(child, "value", errTol);
    } else if (!strcasecmp(child->Value(), "symmetry")) {
      utl::getAttribute(child, "value", symmetry);
    } else if (!strcasecmp(child->Value(), "knot_mult")) {
      utl::getAttribute(child, "value", knot_mult);
    } else if (!strcasecmp(child->Value(), "store_eps_mesh")) {
      utl::getAttribute(child, "value", storeMesh);
    } else if (!strcasecmp(child->Value(), "scheme")) {
      utl::getAttribute(child, "value", str, true);
      if( !strcasecmp(str.c_str(), "fullspan") )
        scheme = 0;
      else if( !strcasecmp(str.c_str(), "minspan") )
        scheme = 1;
      else if( !strcasecmp(str.c_str(), "isotropic_element") )
        scheme = 2;
      else if( !strcasecmp(str.c_str(), "isotropic_function") )
        scheme = 3;
      else
      {
      	std::cerr << "Error parsing adaptive refinement scheme: unknown value\n";
      	scheme = 0;
        return false;
      }
    }

    child = child->NextSiblingElement();
  }

  return true;
}

bool AdaptiveSIM::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"ADAPTIVE",8))
  {
    std::istringstream cline(utl::readLine(is));

    cline >> beta >> errTol;
    if (cline.fail() || cline.bad()) return false;

    int itmp;
    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      maxStep = itmp;

    cline >> itmp;
    if (!cline.fail() && !cline.bad())
      maxDOFs = itmp;

    cline >> itmp; // read knotline multiplicity
    if (!cline.fail() && !cline.bad())
      knot_mult = itmp;

    cline >> itmp; // read refinement scheme
    // (0=refine all, 1=minimum span, 2=isotropic)
    if (!cline.fail() && !cline.bad())
      scheme = itmp;

    cline >> itmp; // read symmetry
    // (0=none, <n> always requests a multiplum of <n> elements for refinement)
    if (!cline.fail() && !cline.bad())
      symmetry = itmp;
  }
  else
    return model->parse(keyWord,is);

  return true;
}


bool AdaptiveSIM::solveStep (const char* inputfile, SystemMatrix::Type solver,
			     const std::map<SIMbase::ProjectionMethod,std::string>& pOpt,
			     size_t adaptor, int iStep)
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
  else if (storeMesh)
    // Output the initial grid to eps-file
    model->refine(std::vector<int>(),std::vector<int>(),"mesh_001.eps");

  // Assemble the linear FE equation system
  model->setMode(SIM::STATIC,true);
  model->initSystem(iStep == 1 ? SystemMatrix::DENSE : solver, 1,1);
  model->setAssociatedRHS(0,0);
  if (!model->assembleSystem())
    return false;

  // Solve the linear system of equations
  if (!model->solveSystem(linsol,1))
    return false;

  Matrix  ssol;
  Vectors projs;
  std::map<SIMbase::ProjectionMethod,std::string>::const_iterator pit;

  // Project the secondary solution onto the splines basis
  model->setMode(SIM::RECOVERY);
  for (pit = pOpt.begin(); pit != pOpt.end(); pit++)
    if (!model->project(ssol,linsol,pit->first))
      return false;
    else
      projs.push_back(ssol);

  // Evaluate solution norms
  return (model->solutionNorms(Vectors(1,linsol),projs,eNorm,gNorm) &&
	  model->dumpResults(linsol,0.0,std::cout,true,6));
}


//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> IndexDouble;


bool AdaptiveSIM::adaptMesh (size_t adaptor, int iStep)
{
  printNorms(gNorm,eNorm,std::cout,adaptor);

  // Define the reference norm
  if (adaptor > gNorm.size() || adaptor > eNorm.rows()) return false;

  double uNorm;
  if (adaptor == 4 && model->haveAnaSol())
    uNorm = gNorm(3); // Using the analytical solution, |u|_ref = |u|
  else // |u|_ref = sqrt( |u^h|^2 + |e^*|^2 )
    uNorm = sqrt(gNorm(1)*gNorm(1) + gNorm(adaptor)*gNorm(adaptor));

  // Check if further refinement is required
  if (iStep > maxStep || model->getNoDOFs() > (size_t)maxDOFs) return false;
  if (eNorm.cols() < 1 || 100.0*gNorm(adaptor) < errTol*uNorm) return false;

  std::vector<int> toBeRefined, options;
  std::vector<IndexDouble> errors;

  options.reserve(4);
  options.push_back(beta);
  options.push_back(knot_mult);
  options.push_back(scheme);
  options.push_back(symmetry);

  size_t i;
  if (scheme == 3)
  {
    // sum up the total function error over all supported elements for that function
    ASMbase* patch = model->getFEModel().front();
    IntMat::const_iterator eit;
    IntVec::const_iterator nit;
    for (i = 0; i < patch->getNoNodes(); i++)
      errors.push_back(IndexDouble(0.0,i));
    for (i = 1, eit = patch->begin_elm(); eit < patch->end_elm(); eit++, i++)
      for (nit = eit->begin(); nit < eit->end(); nit++)
        errors[*nit].first += eNorm(adaptor,i);
  }
  else
    for (i = 0; i < eNorm.cols(); i++)
      errors.push_back(IndexDouble(eNorm(adaptor,1+i),i));

  // Find the list of elements/functions to refine (the beta % with the highest error)
  size_t ipivot = ceil(errors.size()*beta/100.0);
  // make ipivot a multiplum of 'symmetry' in case of symmetric problems
  if (symmetry > 0)
    ipivot += (symmetry-ipivot%symmetry);

  if (ipivot < 1 || ipivot > errors.size()) return false;
  std::sort(errors.begin(),errors.end(),std::greater<IndexDouble>());

  std::cout <<"\nRefining "<< ipivot <<" elements ("<< beta
	    <<"\%) with errors in range ["<< errors[ipivot-1].first
	    <<","<< errors.front().first <<"]"<< std::endl;

  toBeRefined.reserve(ipivot);
  for (i = 0; i < errors.size(); i++)
    toBeRefined.push_back(errors[i].second);

  // Now refine the mesh
  char fname[13];
  sprintf(fname,"mesh_%03d.eps",iStep);
  if(storeMesh)
    return model->refine(toBeRefined,options,fname);
  else
    return model->refine(toBeRefined,options,NULL);
}


std::ostream& AdaptiveSIM::printNorms (const Vector& norms, const Matrix& eNorm,
				       std::ostream& os, size_t adaptor)
{
  // TODO: This needs further work to enable print for all recovered solutions.
  // As for now we only print the norm used for the mesh adaption.
  // Will be completed after the kmo/xmlfixes branch is merged with trunk.
  os <<    "Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< norms(1)
     <<  "\nExternal energy ((h,u^h)+(t,u^h)^0.5 : "<< norms(2);
  if (adaptor == 4 && norms.size() >= 4)
    os <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< norms(3)
       <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< norms(4)
       <<"\nExact relative error (%) : "<< 100.0*norms(4)/norms(3);
  else if (adaptor > 4 && adaptor <= norms.size())
    os <<"\nError estimate a(e,e)^0.5, e=u^r-u^h : "<< norms(adaptor)
       <<"\nRelative error (%) : "<< 100.0*norms(adaptor)/
      sqrt(norms(1)*norms(1) + norms(adaptor)*norms(adaptor));

  if (eNorm.rows() < adaptor || eNorm.cols() < 1)
    return os << std::endl;

  // Compute some additional error measures

  size_t i;
  double avg_norm = eNorm(adaptor,1);
  double min_err  = avg_norm;
  double max_err  = avg_norm;
  for (i = 2; i <= eNorm.cols(); i++)
  {
    avg_norm += eNorm(adaptor,i);
    if (min_err > eNorm(adaptor,i))
      min_err = eNorm(adaptor,i);
    else if (max_err < eNorm(adaptor,i))
      max_err = eNorm(adaptor,i);
  }
  avg_norm /= eNorm.cols();

  double RMS_norm = 0.0;
  for (i = 1; i <= eNorm.cols(); i++)
    RMS_norm += pow(eNorm(adaptor,i)-avg_norm,2.0);
  RMS_norm = sqrt(RMS_norm/eNorm.cols())/avg_norm;

  os <<"\nRoot mean square (RMS) of error      : "<< RMS_norm
     <<"\nMin element error                    : "<< min_err
     <<"\nMax element error                    : "<< max_err
     <<"\nAverage element error                : "<< avg_norm;
  return os << std::endl;
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
