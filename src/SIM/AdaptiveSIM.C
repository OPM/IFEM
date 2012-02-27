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
  {
  	// if no adaptive word is found, set default refine options and return 
    options.clear();
    options.push_back(beta);
    options.push_back(knot_mult);
    options.push_back(scheme);
    options.push_back(symmetry);
    return model->parse(elem);
  }
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
  options.clear();
  options.push_back(beta);
  options.push_back(knot_mult);
  options.push_back(scheme);
  options.push_back(symmetry);

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

    options.clear();
    options.push_back(beta);
    options.push_back(knot_mult);
    options.push_back(scheme);
    options.push_back(symmetry);
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
  else if(storeMesh)
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


//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> IndexDouble;


bool AdaptiveSIM::adaptMesh (int iStep)
{
  printNorms(gNorm,eNorm,std::cout);

  // Check if further refinement is required
  if (iStep > maxStep || model->getNoDOFs() > (size_t)maxDOFs) return false;
  if (gNorm.size() > 3 && 100.0*gNorm(4) < errTol*gNorm(3)) return false;
  if (eNorm.cols() < 1 || eNorm.rows() < 4) return false;

  std::vector<int> toBeRefined;
  std::vector<IndexDouble> errors;

  size_t i;
  if (options.size() > 2 && options[2] == 3)
  {
    // sum up the total function error over all supported elements for that function
    ASMbase* patch = model->getFEModel().front();
    IntMat::const_iterator eit;
    IntVec::const_iterator nit;
    for (i = 0; i < patch->getNoNodes(); i++)
      errors.push_back(IndexDouble(0.0,i));
    for (i = 1, eit = patch->begin_elm(); eit < patch->end_elm(); eit++, i++)
      for (nit = eit->begin(); nit < eit->end(); nit++)
        errors[*nit].first += eNorm(4,i);
  }
  else
    for (i = 0; i < eNorm.cols(); i++)
      errors.push_back(IndexDouble(eNorm(4,1+i),i));

  // Find the list of elements/functions to refine (the beta % with the highest error)
  size_t ipivot = ceil(errors.size()*beta/100.0);
  // make ipivot a multiplum of options[3] in case of symmetric problems
  if (options.size() > 3 && options[3] > 0)
    ipivot += (options[3]-ipivot%options[3]);

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
				       std::ostream& os)
{
  os <<    "Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< norms(1)
     <<  "\nExternal energy ((h,u^h)+(t,u^h)^0.5 : "<< norms(2);
  if (norms.size() > 2)
    os <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< norms(3);
  if (norms.size() > 3)
    os <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< norms(4)
       <<"\nExact relative error (%) : "<< 100.0*norms(4)/norms(3);

  if (eNorm.rows() < 4 || eNorm.cols() < 1)
    return os << std::endl;

  // Compute some additional error measures

  size_t i;
  double avg_norm = eNorm(4,1);
  double min_err  = avg_norm;
  double max_err  = avg_norm;
  for (i = 2; i <= eNorm.cols(); i++)
  {
    avg_norm += eNorm(4,i);
    if (min_err > eNorm(4,i))
      min_err = eNorm(4,i);
    else if (max_err < eNorm(4,i))
      max_err = eNorm(4,i);
  }
  avg_norm /= eNorm.cols();

  double RMS_norm = 0.0;
  for (i = 1; i <= eNorm.cols(); i++)
    RMS_norm += pow(eNorm(4,i)-avg_norm,2.0);
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
