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
#include "SystemMatrix.h"
#include "Utilities.h"
#include "tinyxml.h"
#include <sstream>
#include <cstdio>


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
  adaptor   = 0;
}


AdaptiveSIM::~AdaptiveSIM ()
{
  if (model) delete model;
}


bool AdaptiveSIM::parse (const TiXmlElement* elem)
{
  if (strcasecmp(elem->Value(),"adaptive"))
    return model->parse(elem);

  const char* value = 0;
  const TiXmlElement* child = elem->FirstChildElement();
  while (child) {
    if ((value = utl::getValue(child,"maxstep")))
      maxStep = atoi(value);
    else if ((value = utl::getValue(child,"beta")))
      beta = atof(value);
    else if ((value = utl::getValue(child,"maxdof")))
      maxDOFs = atoi(value);
    else if ((value = utl::getValue(child,"errtol")))
      errTol = atof(value);
    else if ((value = utl::getValue(child,"symmetry")))
      symmetry = atoi(value);
    else if ((value = utl::getValue(child,"knot_mult")))
      knot_mult = atoi(value);
    else if (!strcasecmp(child->Value(), "store_eps_mesh"))
      storeMesh = true; // no need for value here
    else if ((value = utl::getValue(child,"scheme"))) {
      if (!strcasecmp(value,"fullspan"))
        scheme = 0;
      else if (!strcasecmp(value,"minspan"))
        scheme = 1;
      else if (!strcasecmp(value,"isotropic_element"))
        scheme = 2;
      else if (!strcasecmp(value,"isotropic_function"))
        scheme = 3;
      else
      	std::cerr <<"  ** AdaptiveSIM::parse: Unknown refinement scheme \""
                  << value <<"\" (ignored)"<< std::endl;
    }
    else if ((value = utl::getValue(child,"use_norm")))
      adaptor = atoi(value);

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


bool AdaptiveSIM::initAdaptor (size_t indxProj, size_t nNormProj)
{
  if (indxProj == 0) indxProj = adaptor; // use value from XML input

  SIMoptions::ProjectionMap::const_iterator pit = model->opt.project.begin();
  for (size_t j = 0; pit != model->opt.project.end(); pit++, j++)
    if (j+1 == indxProj)
    {
      // Compute the index into eNorm for the error indicator to adapt on
      adaptor = model->haveAnaSol() ? 6+(nNormProj+2)*j : 4+nNormProj*j;
      break;
    }

  std::cout <<"\n\n >>> Starting adaptive simulation based on";
  if (pit != model->opt.project.end())
    std::cout <<"\n     "<< pit->second <<" error estimates (index="
                << adaptor <<") <<<\n";
  else if (model->haveAnaSol())
  {
    std::cout <<" exact errors <<<\n";
    adaptor = 4;
  }
  else
  {
    std::cout <<" - nothing, bailing out ...\n";
    return false;
  }

  if (model->opt.format >= 0)
    prefix.reserve(model->opt.project.size()+1);

  return true;
}


bool AdaptiveSIM::solveStep (const char* inputfile, int iStep)
{
  std::cout <<"\nAdaptive step "<< iStep << std::endl;
  if (iStep > 1)
  {
    // Re-generate the FE model after the refinement
    model->clearProperties();
#ifdef HAS_LRSPLINE
    ASMunstruct::resetNumbering();
#endif
    // Caution: If we are using XML-input and have specified old command-line
    // options in order to override simulation options read from the input file,
    // those options will not be overridden here, so please don't do that..
    if (!model->read(inputfile) || !model->preprocess())
      return false;
  }
  else if (storeMesh)
    // Output the initial grid to eps-file
    model->refine(std::vector<int>(),std::vector<int>(),"mesh_001.eps");

  // Assemble the linear FE equation system
  model->setMode(SIM::STATIC,true);
  model->initSystem(iStep == 1 ? SystemMatrix::DENSE : model->opt.solver, 1, 1);
  model->setAssociatedRHS(0,0);
  model->setQuadratureRule(model->opt.nGauss[0],true);
  if (!model->assembleSystem())
    return false;

  // Solve the linear system of equations
  if (!model->solveSystem(linsol,1))
    return false;

  projs.clear();
  projs.reserve(model->opt.project.size());

  // Project the secondary solution onto the splines basis
  model->setMode(SIM::RECOVERY);
  SIMoptions::ProjectionMap::const_iterator pit;
  for (pit = model->opt.project.begin(); pit != model->opt.project.end(); pit++)
  {
    Matrix ssol;
    if (!model->project(ssol,linsol,pit->first))
      return false;

    projs.push_back(ssol);
    if (iStep == 1 && model->opt.format >= 0)
      prefix.push_back(pit->second.c_str());
  }
  if (iStep == 1 && model->opt.format >= 0)
    prefix.push_back(NULL);

  if (msgLevel > 1 && !projs.empty())
    std::cout << std::endl;

  // Evaluate solution norms
  model->setQuadratureRule(model->opt.nGauss[1]);
  return (model->solutionNorms(Vectors(1,linsol),projs,eNorm,gNorm) &&
	  model->dumpResults(linsol,0.0,std::cout,true,6));
}


//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> IndexDouble;


bool AdaptiveSIM::adaptMesh (int iStep)
{
  printNorms(gNorm,eNorm,std::cout,adaptor,model->haveAnaSol());

  if (adaptor > gNorm.size() || adaptor > eNorm.rows())
    return false;

  // Define the reference norm
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
    // Sum up the total error over all supported elements for each function
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

  // Find the list of elements to refine (the beta % with the highest error)
  size_t ipivot = ceil(errors.size()*beta/100.0);
  // Make ipivot a multiplum of 'symmetry' in case of symmetric problems
  if (symmetry > 0)
    ipivot += (symmetry-ipivot%symmetry);

  if (ipivot < 1 || ipivot > errors.size()) return false;
  std::sort(errors.begin(),errors.end(),std::greater<IndexDouble>());

  std::cout <<"\nRefining "<< ipivot <<" elements ("<< beta
	    <<"\%) with errors in range ["<< errors[ipivot-1].first
	    <<","<< errors.front().first <<"]"<< std::endl;

  toBeRefined.reserve(ipivot);
  for (i = 0; i < ipivot; i++)
    toBeRefined.push_back(errors[i].second);

  // Now refine the mesh
  if (!storeMesh)
    return model->refine(toBeRefined,options);

  char fname[13];
  sprintf(fname,"mesh_%03d.eps",iStep);
  return model->refine(toBeRefined,options,fname);
}


std::ostream& AdaptiveSIM::printNorms (const Vector& norms, const Matrix& eNorm,
				       std::ostream& os, size_t adaptor,
				       bool withExact)
{
  // TODO: This needs further work to enable print for all recovered solutions.
  // As for now we only print the norm used for the mesh adaption.
  os <<    "Energy norm |u^h| = a(u^h,u^h)^0.5   : "<< norms(1)
     <<  "\nExternal energy ((h,u^h)+(t,u^h)^0.5 : "<< norms(2);
  if ((adaptor == 4 || withExact) && norms.size() >= 4)
    os <<"\nExact norm  |u|   = a(u,u)^0.5       : "<< norms(3)
       <<"\nExact error a(e,e)^0.5, e=u-u^h      : "<< norms(4)
       <<"\nExact relative error (%) : "<< 100.0*norms(4)/norms(3);
  if (adaptor != 4 && adaptor <= norms.size())
  {
    os <<"\nError estimate a(e,e)^0.5, e=u^r-u^h : "<< norms(adaptor)
       <<"\nRelative error (%) : "<< 100.0*norms(adaptor)/
      sqrt(norms(1)*norms(1) + norms(adaptor)*norms(adaptor));
    if (withExact && norms.size() >= 4)
      os <<"\nEffectivity index  : "<< norms(adaptor)/norms(4);
  }
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


bool AdaptiveSIM::writeGlv (const char* infile, int iStep, int& nBlock,
			    size_t nNormProj)
{
  if (model->opt.format < 0) return true;

  // Write VTF-file with model geometry
  if (!model->writeGlvG(nBlock, iStep == 1 ? infile : 0))
    return false;

  // Write boundary tractions, if any
  if (!model->writeGlvT(iStep,nBlock))
    return false;

  // Write Dirichlet boundary conditions
  if (!model->writeGlvBC(nBlock,iStep))
    return false;

  // Write solution fields
  if (!model->writeGlvS(linsol,iStep,nBlock))
    return false;

  // Write projected solution fields
  SIMoptions::ProjectionMap::const_iterator pit = model->opt.project.begin();
  for (size_t i = 0; i < projs.size(); i++, pit++)
    if (!model->writeGlvP(projs[i],iStep,nBlock,100+10*i,pit->second.c_str()))
      return false;

  // Write element norms
  if (model->haveAnaSol()) nNormProj += 2;
  if (!model->writeGlvN(eNorm,iStep,nBlock,&prefix.front(),nNormProj))
    return false;

  // Write state information
  return model->writeGlvStep(iStep,iStep,1);
}
