// $Id$
//==============================================================================
//!
//! \file AdaptiveSIM.C
//!
//! \date Sep 22 2011
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Adaptive solution driver for linear static FEM simulators.
//!
//==============================================================================

#include "AdaptiveSIM.h"
#include "SIMoutput.h"
#include "SIMenums.h"
#include "ASMunstruct.h"
#include "IFEM.h"


AdaptiveSIM::AdaptiveSIM (SIMoutput& sim, bool sa)
  : SIMadmin(sim), AdaptiveSetup(sim,sa)
{
  geoBlk = nBlock = 0;
  solution.resize(1);
}


bool AdaptiveSIM::parse (const TiXmlElement* elem)
{
  return this->AdaptiveSetup::parse(elem);
}


bool AdaptiveSIM::parse (char* keyWord, std::istream& is)
{
  return this->AdaptiveSetup::parse(keyWord,is);
}


bool AdaptiveSIM::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  return model.preprocess(ignored,fixDup);
}


bool AdaptiveSIM::initAdaptor (size_t normGroup)
{
  if (!this->initPrm(normGroup))
    return false;

  projs.resize(opt.project.size());
  if (model.haveDualSol())
    projd.resize(opt.project.size());

  if (opt.format >= 0)
  {
    prefix.reserve(opt.project.size());
    for (const SIMoptions::ProjectionMap::value_type& prj : opt.project)
      prefix.push_back(prj.second);
  }

  return true;
}


bool AdaptiveSIM::assembleAndSolveSystem ()
{
  // Assemble the linear FE equation system
  if (!model.setMode(SIM::STATIC,true,true))
    return false;

  model.setQuadratureRule(opt.nGauss[0],true);
  if (!model.assembleSystem())
    return false;

  // Solve the linear system of equations
  int printSol = 1;
  solution.resize(model.getNoRHS());
  for (size_t i = 0; i < solution.size(); i++)
    if (!model.solveSystem(solution[i],printSol,&rCond,"displacement",i))
      return false;
    else if (i == 0)
    {
      for (double value : solution[i])
        if (std::isnan(value))
        {
          std::cerr <<" *** Solution contains NaN, aborting..."<< std::endl;
          return false;
        }
    }
    else if (solution.size() > 2)
      printSol = 0; // Print summary only for the first two solutions

  return true;
}


bool AdaptiveSIM::solveStep (const char* inputfile, int iStep, bool withRF,
                             std::streamsize precision)
{
  // Lambda function to write refined mesh to VTF on failure
  auto&& failure = [this,iStep]()
  {
    if (opt.format >= 0)
      if (model.writeGlvG(geoBlk,nullptr))
        model.writeGlvStep(iStep,iStep,1);

    return false;
  };

  gNorm.clear();
  dNorm.clear();
  eNorm.clear();
  fNorm.clear();

  model.getProcessAdm().cout <<"\nAdaptive step "<< iStep << std::endl;
  if (iStep > 1)
  {
    SIMoptions oldOpt(opt);
    // Re-generate the FE model after the refinement
    model.clearProperties();
    // Caution: If we are using XML-input and have specified old command-line
    // options in order to override simulation options read from the input file,
    // those options will not be overridden here, so please don't do that..
    if (!model.read(inputfile) || !model.preprocess())
      return failure();
    opt = oldOpt;
  }
  else
    this->writeMesh(1); // Output initial grid to eps-file(s)

  // Initialize the linear equation system for the current grid
  if (!model.initSystem(opt.solver,1,model.getNoRHS(),0,withRF))
    return failure();

  // Assemble and solve the FE equation system
  if (!this->assembleAndSolveSystem())
    return failure();

  // Project the secondary solution onto the splines basis
  size_t idx = 0;
  model.setMode(SIM::RECOVERY);
  for (const SIMoptions::ProjectionMap::value_type& prj : opt.project)
    if (prj.first <= SIMoptions::NONE)
      idx++; // No projection for this norm group
    else if (!model.project(projs[idx++],solution.front(),prj.first))
      return failure();
    else if (idx == adaptor && idx <= projd.size() && solution.size() > 1)
      if (!model.project(projd[idx-1],solution[1],prj.first))
        return failure();

  if (msgLevel > 1 && !projs.empty())
    model.getProcessAdm().cout << std::endl;

  // Evaluate solution norms
  model.setMode(SIM::NORMS);
  model.setQuadratureRule(opt.nGauss[1]);
  if (!model.solutionNorms(solution.front(),projs,eNorm,gNorm))
    return failure();

  if (!projd.empty() && solution.size() > 1)
  {
    if (!model.solutionNorms(solution[1],projd,fNorm,dNorm))
      return failure();

    if (eRow <= fNorm.rows() && eRow <= eNorm.rows())
    {
      // Calculate refinement indicators including the dual error estimates.
      // Store in the 2nd row in fNorm (normally holding the external energy).
      // Store the associated global norm in the dNorm array.
      dNorm.front()(2) = 0.0;
      for (size_t j = 1; j <= eNorm.cols() && j <= fNorm.cols(); j++)
        dNorm.front()(2) += fNorm(2,j) = eNorm(eRow,j)*fNorm(eRow,j);
    }
  }

  model.setMode(SIM::RECOVERY);
  if (!model.dumpResults(solution.front(),0.0,
                         model.getProcessAdm().cout,true,precision))
    return failure();

  if (!this->savePoints(0.0, iStep))
    return failure();

  return true;
}


//! \brief Element error and associated index.
//! \note The error value must be first and the index second, such that the
//! internally defined greater-than operator can be used when sorting the
//! error+index pairs in decreasing error order.
typedef std::pair<double,int> DblIdx;


bool AdaptiveSIM::adaptMesh (int iStep, std::streamsize outPrec)
{
  if (iStep < 2)
    return true; // No refinement in the first adaptive cycle

  if (outPrec > 0)
  {
    std::streamsize oldPrec = IFEM::cout.precision(outPrec);
    this->printNorms(gNorm,dNorm,eNorm);
    IFEM::cout.precision(oldPrec);
  }
  else
    this->printNorms(gNorm,dNorm,eNorm);

  // Set up refinement parameters
  LR::RefineData prm;
  Vector refIn = fNorm.empty() ? eNorm.getRow(eRow) : fNorm.getRow(2);
  if (this->calcRefinement(prm,iStep,gNorm,refIn) <= 0)
    return false;

  // Now refine the mesh and write out resulting grid
  return model.refine(prm) & this->writeMesh(iStep);
}


bool AdaptiveSIM::writeGlv (const char* infile, int iStep)
{
  if (opt.format < 0)
    return true;

  // Write VTF-file with model geometry
  if (!model.writeGlvG(geoBlk, iStep == 1 ? infile : nullptr))
    return false;

  // Write boundary tractions, if any
  if (!model.writeGlvT(iStep,geoBlk,nBlock))
    return false;

  // Write Dirichlet boundary conditions
  if (!model.writeGlvBC(nBlock,iStep))
    return false;

  // Write solution fields
  model.setMode(SIM::RECOVERY);
  if (!model.writeGlvS(solution.front(),iStep,nBlock))
    return false;

  if (solution.size() > 1)
    if (!model.writeGlvS1(solution[1],iStep,nBlock,0.0,"Dual solution",90,-1))
      return false;

  // Write projected solution fields
  SIMoptions::ProjectionMap::const_iterator pit = opt.project.begin();
  for (size_t i = 0; i < projs.size(); i++, ++pit)
    if (!model.writeGlvP(projs[i],iStep,nBlock,100+10*i,pit->second.c_str()))
      return false;

  // Write element norms
  if (!model.writeGlvN(eNorm,iStep,nBlock,prefix))
    return false;

  if (!fNorm.empty())
  {
    std::vector<std::string> prefix = { "Dual projected" };
    if (!model.writeGlvN(fNorm,iStep,nBlock,prefix,300,"Dual"))
      return false;
  }

  // Write state information
  return model.writeGlvStep(iStep,iStep,1);
}


bool AdaptiveSIM::savePoints(double time, int iStep) const
{
  return model.savePoints(solution.front(), time, iStep);
}
