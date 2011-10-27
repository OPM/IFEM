// $Id$
//==============================================================================
//!
//! \file SIMbase.C
//!
//! \date Dec 08 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for NURBS-based FEM simulators.
//!
//==============================================================================

#include "SIMbase.h"
#include "ASMbase.h"
#ifdef PARALLEL_PETSC
#include "SAMpatchPara.h"
#include "petscksp.h"
#else
#include "SAMpatch.h"
#endif
#include "IntegrandBase.h"
#include "AlgEqSystem.h"
#include "LinSolParams.h"
#include "GlbNorm.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "Tensor.h"
#include "Vec3Oper.h"
#include "EigSolver.h"
#include "Profiler.h"
#include "ElementBlock.h"
#include "VTF.h"
#include "Utilities.h"
#include <iomanip>
#include <ctype.h>
#include <stdio.h>


ASM::Discretization SIMbase::discretization  = ASM::Spline;
bool                SIMbase::preserveNOrder  = false;
bool                SIMbase::ignoreDirichlet = false;
int                 SIMbase::num_threads_SLU = 1;


SIMbase::SIMbase () : g2l(&myGlb2Loc)
{
  myProblem = 0;
  mySol = 0;
  myVtf = 0;
  myEqSys = 0;
  mySam = 0;
  mySolParams = 0;
  nGlPatches = 0;
  vizIncr = 1;
  format = 1;

  MPCLess::compareSlaveDofOnly = true; // to avoid multiple slave definitions

#ifdef PARALLEL_PETSC
  // In parallel simulations, we need to retain all DOFs in the equation system.
  // The fixed DOFs (if any) will receive a homogeneous constraint instead.
  ASMbase::fixHomogeneousDirichlet = false;
#endif
}


SIMbase::~SIMbase ()
{
#ifdef SP_DEBUG
  std::cout <<"\nEntering SIMbase destructor"<< std::endl;
#endif

  if (myProblem)   delete myProblem;
  if (mySol)       delete mySol;
  if (myVtf)       delete myVtf;
  if (myEqSys)     delete myEqSys;
  if (mySam)       delete mySam;
  if (mySolParams) delete mySolParams;

  for (FEModelVec::iterator i1 = myModel.begin(); i1 != myModel.end(); i1++)
    delete *i1;

  myModel.clear();
  this->clearProperties();

#ifdef SP_DEBUG
  std::cout <<"Leaving SIMbase destructor"<< std::endl;
#endif
}


void SIMbase::clearProperties ()
{
  for (FEModelVec::iterator i1 = myModel.begin(); i1 != myModel.end(); i1++)
    (*i1)->clear(true); // retain the geometry only
  for (SclFuncMap::iterator i2 = myScalars.begin(); i2 != myScalars.end(); i2++)
    delete i2->second;
  for (VecFuncMap::iterator i3 = myVectors.begin(); i3 != myVectors.end(); i3++)
    delete i3->second;
  for (TracFuncMap::iterator i4 = myTracs.begin(); i4 != myTracs.end(); i4++)
    delete i4->second;

  myPatches.clear();
  myScalars.clear();
  myVectors.clear();
  myTracs.clear();
  myProps.clear();
}


bool SIMbase::parse (char* keyWord, std::istream& is)
{
  if (!strncasecmp(keyWord,"LINEARSOLVER",12))
    this->readLinSolParams(is,atoi(keyWord+12));

  else if (!strncasecmp(keyWord,"PARTITIONING",12))
  {
    int nproc = atoi(keyWord+12);
    if (myPid == 0)
      std::cout <<"\nNumber of partitions: "<< nproc << std::endl;

    char* cline = 0;
    nGlPatches = 0;
    for (int i = 0; i < nproc && (cline = utl::readLine(is)); i++)
    {
      int proc  = atoi(strtok(cline," "));
      int first = atoi(strtok(NULL," "));
      int last  = atoi(strtok(NULL," "));

      if (last > nGlPatches) nGlPatches = last;

      if (proc == myPid && last >= first)
      {
	myPatches.reserve(last-first+1);
	for (int j = first; j <= last; j++)
	  myPatches.push_back(j);
      }
    }

    if (myPatches.empty())
    {
      std::cerr <<" *** SIMbase::parse: No partitioning input for processor: "
		<< myPid << std::endl;
      return false;
    }
  }

  else if (!strncasecmp(keyWord,"RESULTPOINTS",12))
  {
    int nres = atoi(keyWord+12);
    if (myPid == 0 && nres > 0)
      std::cout <<"\nNumber of result points: "<< nres;

    char* cline = 0;
    myPoints.resize(nres);
    for (int i = 0; i < nres && (cline = utl::readLine(is)); i++)
    {
      myPoints[i].patch = atoi(strtok(cline," "));
      if (myPid == 0)
	std::cout <<"\n\tPoint "<< i+1 <<": P"<< myPoints[i].patch <<" xi =";
      for (int j = 0; j < 3 && (cline = strtok(NULL," ")); j++)
      {
        myPoints[i].par[j] = atof(cline);
	if (myPid == 0)
	  std::cout <<' '<< myPoints[i].par[j];
      }
    }
    if (myPid == 0 && nres > 0)
      std::cout << std::endl;
  }

  else if (!strncasecmp(keyWord,"VISUALIZATION",13))
  {
    char* cline = utl::readLine(is);
    if ((cline = strtok(cline," "))) format = atoi(cline);
    if ((cline = strtok(NULL," "))) vizIncr = atoi(cline);
  }

#ifdef SP_DEBUG
  // Since the same input file might be parsed by several substep solvers,
  // warnings on ignored keywords are issued when compiled in debug mode only.
  else if (isalpha(keyWord[0]))
    std::cerr <<" *** SIMbase::parse: Unknown keyword: "<< keyWord << std::endl;
#endif

  return true;
}


void SIMbase::readLinSolParams (std::istream& is, int npar)
{
  if (!mySolParams)
    mySolParams = new LinSolParams();

  mySolParams->read(is,npar);
}


bool SIMbase::createFEMmodel ()
{
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel[i]->generateFEMTopology())
      return false;

  if (nGlPatches == 0 && nProc == 1)
    nGlPatches = myModel.size();

  return true;
}


int SIMbase::getLocalPatchIndex (int patchNo) const
{
  if (patchNo < 1 || (patchNo > nGlPatches && nGlPatches > 0))
  {
    std::cerr <<" *** SIMbase::getLocalPatchIndex: Patch number "<< patchNo
	      <<" out of range [1,"<< nGlPatches <<"]"<< std::endl;
    return -1;
  }
  else if (myPatches.empty() || nProc == 1)
    return patchNo;

  for (size_t i = 0; i < myPatches.size(); i++)
    if (myPatches[i] == patchNo)
      return 1+i;

  return 0;
}


bool SIMbase::preprocess (const std::vector<int>& ignored, bool fixDup)
{
  if (!this->createFEMmodel()) return false;

  FEModelVec::const_iterator mit;

  // Erase all patches that should be ignored in the analysis
  std::vector<int>::const_iterator it;
  for (it = ignored.begin(); it != ignored.end(); it++)
    if (*it > 0 && (size_t)*it <= myModel.size())
      myModel[*it-1]->clear();

  // If material properties are specified for at least one patch, assign the
  // property code 999999 to all patches with no material property code yet
  std::vector<int> pchWthMat;
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); p++)
    if (p->pcode == Property::MATERIAL && !myModel[p->patch-1]->empty())
      pchWthMat.push_back(p->patch-1);

  if (!pchWthMat.empty())
    for (size_t i = 0; i < myModel.size(); i++)
      if (!myModel[i]->empty())
	if (std::find(pchWthMat.begin(),pchWthMat.end(),i) == pchWthMat.end())
	  myProps.push_back(Property(Property::MATERIAL,999999,i+1,
				     myModel[i]->getNoSpaceDim()));

  if (fixDup)
  {
    // Check for duplicated nodes (missing topology)
    int patch, nDupl = 0;
    std::map<Vec3,int> globalNodes;
    for (mit = myModel.begin(), patch = 0; mit != myModel.end(); mit++, patch++)
      if (!(*mit)->empty())
      {
	std::cout <<"   * Checking Patch "<< patch << std::endl;
	for (size_t node = 1; node <= (*mit)->getNoNodes(); node++)
	{
	  Vec3 X((*mit)->getCoord(node));
	  std::map<Vec3,int>::const_iterator xit = globalNodes.find(X);
	  if (xit == globalNodes.end())
	    globalNodes.insert(std::make_pair(X,(*mit)->getNodeID(node)));
	  else if ((*mit)->mergeNodes(node,xit->second))
	    nDupl++;
	}
      }
    if (nDupl > 0)
      std::cout <<"   * "<< nDupl <<" duplicated nodes merged."<< std::endl;
  }

  // Renumber the nodes to account for overlapping nodes and erased patches.
  // In parallel simulations, the resulting global-to-local node number mapping
  // will map the global node numbers to local node numbers on the current
  // processor. In serial simulations, the global-to-local mapping will be unity
  // unless the original global node number sequence had "holes" due to
  // duplicated nodes and/or erased patches.
  int ngnod = 0;
  int renum = 0;
  if (preserveNOrder)
  {
    renum = ASMbase::renumberNodes(myModel,myGlb2Loc);
    ngnod = g2l->size();
  }
  else for (mit = myModel.begin(); mit != myModel.end(); mit++)
    renum += (*mit)->renumberNodes(myGlb2Loc,ngnod);

  if (renum > 0)
    std::cout <<"\nRenumbered "<< renum <<" nodes"<< std::endl;

  for (mit = myModel.begin(); mit != myModel.end(); mit++)
    (*mit)->renumberNodes(*g2l);

  // Process the specified Dirichlet boundary conditions
  std::cout <<"\nResolving Dirichlet boundary conditions"<< std::endl;

  bool ok = true;
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); p++)
    switch (p->pcode) {

    case Property::DIRICHLET:
      if (this->addConstraint(p->patch,p->lindx,p->ldim,p->pindx))
	std::cout << std::endl;
      else
	ok = false;
      break;

    case Property::DIRICHLET_INHOM:
      if (this->addConstraint(p->patch,p->lindx,p->ldim,p->pindx%1000,p->pindx))
	std::cout << std::endl;
      else
	ok = false;
      break;

    default:
      break;
    }

  // Set initial values for the inhomogeneous diriclet conditions, if any
  this->initDirichlet();

  // Resolve possibly chained multi-point constraints
  ASMbase::resolveMPCchains(myModel);

  // Preprocess the result points
  for (ResPointVec::iterator p = myPoints.begin(); p != myPoints.end();)
  {
    int pid = this->getLocalPatchIndex(p->patch);
    if (pid < 1 || myModel[pid-1]->empty())
      p = myPoints.erase(p);
    else if ((p->inod = myModel[pid-1]->evalPoint(p->par,p->par,p->X)) < 0)
      p = myPoints.erase(p);
    else
    {
      p->npar = myModel[pid-1]->getNoParamDim();
      int ipt = 1 + (int)(p-myPoints.begin());
      if (ipt == 1) std::cout <<'\n';
      std::cout <<"Result point #"<< ipt <<": patch #"<< p->patch;
      switch (p->npar) {
      case 1: std::cout <<" u="; break;
      case 2: std::cout <<" (u,v)=("; break;
      case 3: std::cout <<" (u,v,w)=("; break;
      }
      std::cout << p->par[0];
      for (unsigned char c = 1; c < p->npar; c++)
	std::cout <<','<< p->par[c];
      if (p->npar > 1) std::cout <<')';
      if (p->inod > 0) std::cout <<", node #"<< p->inod;
      std::cout <<", X = "<< p->X << std::endl;
      p++;
    }
  }

  // Initialize data structures for the algebraic system
  if (mySam) delete mySam;
#ifdef PARALLEL_PETSC
  mySam = new SAMpatchPara(*g2l);
#else
  mySam = new SAMpatch();
#endif

  return mySam->init(myModel,ngnod) && ok;
}


bool SIMbase::setPropertyType (int code, Property::Type ptype, int pindex)
{
  if (code < 0)
  {
    std::cerr <<"  ** SIMbase::setPropertyType: Negative property code "
	      << code <<" (ignored)"<< std::endl;
    return false;
  }

  for (size_t j = 0; j < myProps.size(); j++)
    if (myProps[j].pindx == (size_t)code &&
	myProps[j].pcode == Property::UNDEFINED)
    {
      myProps[j].pcode = ptype;
      if (pindex >= 0) myProps[j].pindx = pindex;
    }

  return true;
}


bool SIMbase::setVecProperty (int code, Property::Type ptype, VecFunc* field)
{
  if (field) myVectors[code] = field;
  return this->setPropertyType(code,ptype);
}


bool SIMbase::initSystem (SystemMatrix::Type mType, size_t nMats, size_t nVec)
{
  if (!mySam) return false;
#if SP_DEBUG > 1
  mySam->print(std::cout);
  std::string heading("\n\nNodal coordinates for Patch 1");
  size_t i, j = heading.size()-1;
  for (i = 0; i < myModel.size() && i < 9; i++, heading[j]++)
    myModel[i]->printNodes(std::cout,heading.c_str());
#endif

  if (myEqSys) delete myEqSys;
  myEqSys = new AlgEqSystem(*mySam);
  myEqSys->init(mType,mySolParams,nMats,nVec,num_threads_SLU);
  myEqSys->initAssembly();
  return true;
}


bool SIMbase::setAssociatedRHS (size_t iMat, size_t iVec)
{
  if (!myEqSys) return false;

  return myEqSys->setAssociatedVector(iMat,iVec);
}


bool SIMbase::setMode (int mode, bool resetSol)
{
  if (!myProblem) return false;

  myProblem->setMode((SIM::SolutionMode)mode);
  if (resetSol) myProblem->resetSolution();

  return true;
}


void SIMbase::printProblem (std::ostream& os) const
{
  if (myProblem)
  {
    std::cout <<"\nProblem definition:"<< std::endl;
    myProblem->print(os);
  }

#if SP_DEBUG > 1
  std::cout <<"\nProperty mapping:";
  for (PropertyVec::const_iterator i = myProps.begin(); i != myProps.end(); i++)
    std::cout <<"\n"<< i->pcode <<" "<< i->pindx <<" "<< i->patch
	      << (int)i->lindx <<" "<< (int)i->ldim;
  std::cout << std::endl;
#endif
}


size_t SIMbase::getNoFields (int basis) const
{
  return myModel.empty() ? 0 : myModel.front()->getNoFields(basis);
}


size_t SIMbase::getNoSpaceDim () const
{
  return myModel.empty() ? 0 : myModel.front()->getNoSpaceDim();
}


size_t SIMbase::getNoDOFs () const
{
  return mySam ? mySam->getNoDOFs() : 0;
}


size_t SIMbase::getNoNodes (bool unique) const
{
  size_t nnod = 0;
  if (unique && mySam)
    nnod = mySam->getNoNodes();
  else for (size_t i = 0; i < myModel.size(); i++)
    nnod += myModel[i]->getNoNodes();

  return nnod;
}


size_t SIMbase::getNoSolutions () const
{
  return myProblem ? myProblem->getNoSolutions() : 0;
}


bool SIMbase::initDirichlet (double time)
{
  Vector dummy;
  return this->updateDirichlet(time,&dummy);
}


bool SIMbase::updateDirichlet (double time, const Vector* prevSol)
{
  if (prevSol)
    for (size_t i = 0; i < myModel.size(); i++)
      if (!myModel[i]->updateDirichlet(myScalars,myVectors,time))
	return false;

  if (mySam)
    return mySam->updateConstraintEqs(myModel,prevSol);
  else
    return true;
}


bool SIMbase::updateGrid (const Vector& displ)
{
  bool ok = true;
  Vector locdisp;
  for (size_t i = 0; i < myModel.size() && ok; i++)
  {
    myModel[i]->extractNodeVec(displ,locdisp,myModel[i]->getNoSpaceDim());
    ok = myModel[i]->updateCoords(locdisp);
  }

  return true;
}


bool SIMbase::assembleSystem (const TimeDomain& time, const Vectors& prevSol,
			      bool newLHSmatrix)
{
  PROFILE1("Element assembly");

  if (!myProblem) return false;

  myProblem->initIntegration(time);
  myEqSys->init(newLHSmatrix);
  bool ok = true;

  // Loop over the different material regions, integrating interior
  // coefficient matrix terms for the patch associated with each material
  size_t i, j = 0, lp = 0;
  for (i = 0; i < myProps.size() && ok; i++)
    if (myProps[i].pcode == Property::MATERIAL)
      if ((j = myProps[i].patch) < 1 || j > myModel.size())
      {
	std::cerr <<" *** SIMbase::assembleSystem: Patch index "<< j
		  <<" out of range [1,"<< myModel.size() <<"]"<< std::endl;
	ok = false;
      }
      else if (this->initMaterial(myProps[i].pindx))
      {
	if (msgLevel > 1)
	  std::cout <<"\nAssembling interior matrix terms for P"<< j
		    << std::endl;
	this->initBodyLoad(j);
	this->extractPatchSolution(myModel[j-1],prevSol);
	ok = myModel[j-1]->integrate(*myProblem,*myEqSys,time);
	lp = j;
      }
      else
	ok = false;

  if (j == 0 && ok)
    // All patches are referring to the same material, and we assume it has
    // been initialized during input processing (thus no initMaterial call here)
    for (i = 0; i < myModel.size() && ok; i++)
    {
      if (msgLevel > 1)
	std::cout <<"\nAssembling interior matrix terms for P"<< i+1
		  << std::endl;
      this->initBodyLoad(i+1);
      this->extractPatchSolution(myModel[i],prevSol);
      ok = myModel[i]->integrate(*myProblem,*myEqSys,time);
      lp = i+1;
    }

  // Assemble right-hand-side contributions from the Neumann boundary conditions
  if (myEqSys->getVector())
    for (i = 0; i < myProps.size() && ok; i++)
      if (myProps[i].pcode == Property::NEUMANN)
	if ((j = myProps[i].patch) < 1 || j > myModel.size())
	{
	  std::cerr <<" *** SIMbase::assembleSystem: Patch index "<< j
		    <<" out of range [1,"<< myModel.size() <<"]"<< std::endl;
	  ok = false;
	}

	else if (myProps[i].ldim+1 == myModel[j-1]->getNoSpaceDim())
	  if (this->initNeumann(myProps[i].pindx))
	  {
	    int bIndex = myProps[i].lindx;
	    if (msgLevel > 1)
	      std::cout <<"\nAssembling Neumann matrix terms for boundary "
			<< bIndex <<" on P"<< j << std::endl;
	    if (j != lp) this->extractPatchSolution(myModel[j-1],prevSol);
	    ok = myModel[j-1]->integrate(*myProblem,bIndex,*myEqSys,time);
	    lp = j;
	  }
	  else
	    ok = false;

	else if (myProps[i].ldim+2 == myModel[j-1]->getNoSpaceDim())
	  if (this->initNeumann(myProps[i].pindx))
	  {
	    int bIndex = myProps[i].lindx;
	    if (msgLevel > 1)
	      std::cout <<"\nAssembling Neumann matrix terms for edge "
			<< bIndex <<" on P"<< j << std::endl;
	    if (j != lp) this->extractPatchSolution(myModel[j-1],prevSol);
	    ok = myModel[j-1]->integrateEdge(*myProblem,bIndex,*myEqSys,time);
	    lp = j;
	  }
	  else
	    ok = false;

  return ok && this->finalizeAssembly(newLHSmatrix);
}


bool SIMbase::finalizeAssembly (bool newLHSmatrix)
{
  // Communication of matrix and vector assembly (for PETSc only)
  SystemMatrix* A = myEqSys->getMatrix();
  if (A && newLHSmatrix)
  {
    if (!A->beginAssembly()) return false;
    if (!A->endAssembly())   return false;
#if SP_DEBUG > 3
    std::cout <<"\nSystem coefficient matrix:"<< *A;
#endif
  }

  SystemVector* b = myEqSys->getVector();
  if (b)
  {
    if (!b->beginAssembly()) return false;
    if (!b->endAssembly())   return false;
#if SP_DEBUG > 2
    std::cout <<"\nSystem right-hand-side vector:"<< *b;
#endif
  }

  return true;
}


bool SIMbase::extractLoadVec (Vector& loadVec) const
{
  SystemVector* b = myEqSys->getVector();
  if (!b) return false;

  // Expand load vector from equation ordering to DOF-ordering
  return mySam->expandSolution(*b,loadVec,0.0);
}


void SIMbase::extractPatchSolution (const ASMbase* patch, const Vectors& sol)
{
  for (size_t i = 0; i < sol.size() && i < myProblem->getNoSolutions(); i++)
    if (!sol[i].empty())
      patch->extractNodeVec(sol[i],myProblem->getSolution(i));
}


bool SIMbase::solveSystem (Vector& solution, int printSol,
			   const char* compName, bool newLHS)
{
  SystemMatrix* A = myEqSys->getMatrix();
  SystemVector* b = myEqSys->getVector();
  if (!A) std::cerr <<" *** SIMbase::solveSystem: No LHS matrix"<< std::endl;
  if (!b) std::cerr <<" *** SIMbase::solveSystem: No RHS vector"<< std::endl;
  if (!A || !b) return false;

  // Solve the linear system of equations
  if (msgLevel > 1)
    std::cout <<"\nSolving the equation system ..."<< std::endl;
  {
    PROFILE1("Equation solving");
    if (!A->solve(*b,newLHS)) return false;
  }

  // Expand solution vector from equation ordering to DOF-ordering
  if (!mySam->expandSolution(*b,solution)) return false;
  if (printSol < 1) return true;

  // Compute and print solution norms
  const size_t nf = this->getNoFields(1);
  size_t iMax[nf];
  double dMax[nf];
  double dNorm = this->solutionNorms(solution,dMax,iMax,nf);

  if (myPid == 0)
  {
    std::cout <<"\n >>> Solution summary <<<\n"
	      <<"\nL2-norm            : "<< utl::trunc(dNorm);
    if (nf == 1 && utl::trunc(dMax[0]) != 0.0)
      std::cout <<"\nMax "<< compName <<"   : "<< dMax[0] <<" node "<< iMax[0];
    else if (nf > 1)
    {
      char D = 'X';
      for (size_t d = 0; d < nf; d++, D++)
	if (utl::trunc(dMax[d]) != 0.0)
	  std::cout <<"\nMax "<< D <<'-'<< compName <<" : "
		    << dMax[d] <<" node "<< iMax[d];
    }
    std::cout << std::endl;
  }

  // Print entire solution vector if it is small enough
  if (mySam->getNoEquations() < printSol)
  {
    std::cout <<"\nSolution vector:";
    for (int inod = 1; inod <= mySam->getNoNodes(); inod++)
    {
      std::cout <<"\nNode "<< inod <<":";
      std::pair<int,int> dofs = mySam->getNodeDOFs(inod);
      for (int d = dofs.first-1; d < dofs.second; d++)
	std::cout <<" "<< utl::trunc(solution[d]);
    }
    std::cout << std::endl;
  }
#if SP_DEBUG > 2
  else
    std::cout <<"\nSolution vector:"<< *b;
#endif

  return true;
}


void SIMbase::iterationNorms (const Vector& u, const Vector& r,
			      double& eNorm, double& rNorm, double& dNorm) const
{
  eNorm = mySam->dot(r,u,'A');
  rNorm = mySam->norm2(r,'D');
  dNorm = mySam->norm2(u,'D');
}


double SIMbase::solutionNorms (const Vector& x, double* inf,
			       size_t* ind, size_t nf) const
{
  if (nf == 0) nf = this->getNoSpaceDim();

  for (size_t d = 0; d < nf; d++)
  {
    ind[d] = d+1;
    inf[d] = mySam->normInf(x,ind[d]);
  }

  return mySam->normL2(x);
}


NormBase* SIMbase::getNormIntegrand () const
{
  return myProblem->getNormIntegrand(mySol);
}


bool SIMbase::solutionNorms (const TimeDomain& time,
			     const Vectors& psol, const Vectors& ssol,
			     Vector& gNorm, Matrix* eNorm)
{
  NormBase* norm = myProblem->getNormIntegrand(mySol);
  if (!norm)
  {
#ifdef SP_DEBUG
    std::cerr <<" *** SIMbase::solutionNorms: No integrand."<< std::endl;
#endif
    return false;
  }

  PROFILE1("Norm integration");

  myProblem->initIntegration(time);
  const Vector& primsol = psol.front();

  size_t i, j, k;
  for (i = 0; i < ssol.size(); i++)
    if (!ssol[i].empty())
      norm->getProjection(i+1);

  size_t nCmp = ssol.empty() ? 0 : ssol.front().size() / mySam->getNoNodes();
  size_t nNorms = norm->getNoFields();
  gNorm.resize(nNorms,true);

  // Initialize norm integral classes
  GlbNorm globalNorm(gNorm,GlbNorm::SQRT);
  LintegralVec elementNorms;
  if (eNorm)
  {
    eNorm->resize(nNorms,mySam->getNoElms(),true);
    elementNorms.reserve(eNorm->cols());
    for (i = 0; i < eNorm->cols(); i++)
      elementNorms.push_back(new ElmNorm(eNorm->ptr(i),nNorms));
  }

  // Loop over the different material regions, integrating solution norm terms
  // for the patch domain associated with each material
  bool ok = true;
  size_t lp = 0;
  for (i = j = 0; i < myProps.size() && ok; i++)
    if (myProps[i].pcode == Property::MATERIAL)
      if ((j = myProps[i].patch) < 1 || j > myModel.size())
	ok = false;
      else if (this->initMaterial(myProps[i].pindx))
      {
	myModel[j-1]->extractNodeVec(primsol,myProblem->getSolution());
	for (k = 0; k < ssol.size(); k++)
	  if (!ssol[k].empty())
	    myModel[j-1]->extractNodeVec(ssol[k],norm->getProjection(k+1),nCmp);
	ok = myModel[j-1]->integrate(*norm,globalNorm,time,elementNorms);
	lp = j;
      }
      else
	ok = false;

  if (j == 0 && ok)
    // All patches are referring to the same material, and we assume it has
    // been initialized during input processing (thus no initMaterial call here)
    for (i = 0; i < myModel.size() && ok; i++)
    {
      myModel[i]->extractNodeVec(primsol,myProblem->getSolution());
      for (k = 0; k < ssol.size(); k++)
	if (!ssol[k].empty())
	  myModel[i]->extractNodeVec(ssol[k],norm->getProjection(k+1),nCmp);
      ok = myModel[i]->integrate(*norm,globalNorm,time,elementNorms);
      lp = i+1;
    }

  // Integrate norm contributions due to Neumann boundary conditions, if any.
  // Note: Currently, only the global norms are considered here.
  // The corresponding element-level norms are not stored. This is mainly
  // because the current design only allows one loop over the elements
  // in the element-norm calculations. Consider rivising this later.
  if (norm->hasBoundaryTerms())
    for (i = 0; i < myProps.size() && ok; i++)
      if (myProps[i].pcode == Property::NEUMANN)
	if ((j = myProps[i].patch) < 1 || j > myModel.size())
	  ok = false;

	else if (myProps[i].ldim+1 == myModel[j-1]->getNoSpaceDim())
	  if (this->initNeumann(myProps[i].pindx))
	  {
	    if (j != lp)
	      myModel[j-1]->extractNodeVec(primsol,myProblem->getSolution());
	    ok = myModel[j-1]->integrate(*norm,myProps[i].lindx,
					 globalNorm,time);
	    lp = j;
	  }
	  else
	    ok = false;

	else if (myProps[i].ldim+2 == myModel[j-1]->getNoSpaceDim())
	  if (this->initNeumann(myProps[i].pindx))
	  {
	    if (j != lp)
	      myModel[j-1]->extractNodeVec(primsol,myProblem->getSolution());
	    ok = myModel[j-1]->integrateEdge(*norm,myProps[i].lindx,
					     globalNorm,time);
	    lp = j;
	  }
	  else
	    ok = false;

  // Add problem-dependent external norm contributions
  gNorm(2) += this->externalEnergy(psol);

#ifdef PARALLEL_PETSC
  if (nProc > 1)
  {
    double* tmp = new double[nNorms];
    MPI_Allreduce(gNorm.ptr(),tmp,nNorms,MPI_DOUBLE,MPI_SUM,PETSC_COMM_WORLD);
    memcpy(gNorm.ptr(),tmp,nNorms*sizeof(double));
    delete[] tmp;
  }
#endif

  delete norm;
  for (i = 0; i < elementNorms.size(); i++)
    delete elementNorms[i];

  return ok;
}


double SIMbase::externalEnergy (const Vectors& psol) const
{
  const Vector* reactionForces = myEqSys->getReactions();
  if (!reactionForces || psol.empty()) return 0.0;

  // Add norm contributions due to inhomogeneous Dirichlet boundary conditions.
  // That is, the path integral of the total solution vector times the
  // reaction forces at the prescribed DOFs.
  if (psol.size() == 1)
    return mySam->normReact(psol.front(),*reactionForces);

  static double extEnergy = 0.0;
  static Vector prevForces(reactionForces->size());
  extEnergy += mySam->normReact(psol[0]-psol[1],*reactionForces+prevForces);
  prevForces = *reactionForces;
  return extEnergy;
}


/*!
  The content of the output array \a RF is as follows:
  \f[
  RF[0] = \sum_{i=n}^{\rm nnod} \sum_{i=1}^{\rm nsd} f_n^i u_n^i
  \quad\quad\mbox{(energy)}\f]
  \f[
  RF[i] = \sum_{n=1}^{\rm nnod} f_n^i \quad\forall\quad i=1,\ldots,{\rm nsd}
  \f]
*/

bool SIMbase::getCurrentReactions (RealArray& RF, const Vector& psol) const
{
  const Vector* reactionForces = myEqSys->getReactions();
  if (!reactionForces) return false;

  RF.resize(1+this->getNoSpaceDim());
  RF.front() = 2.0*mySam->normReact(psol,*reactionForces);
  for (size_t dir = 1; dir < RF.size(); dir++)
    RF[dir] = mySam->getReaction(dir,*reactionForces);

  return true;
}


bool SIMbase::systemModes (std::vector<Mode>& solution,
			   int nev, int ncv, int iop, double shift,
			   size_t iA, size_t iB)
{
  if (nev < 1 || ncv <= nev) return false;

  PROFILE1("Eigenvalue analysis");

  Vector eigVal;
  Matrix eigVec;
  if (nev > mySam->getNoEquations()) nev = mySam->getNoEquations();
  if (ncv > mySam->getNoEquations()) ncv = mySam->getNoEquations();

  // Solve the eigenvalue problem
  std::cout <<"\nSolving the eigenvalue problem ..."<< std::endl;
  SystemMatrix* A = myEqSys->getMatrix(iA);
  SystemMatrix* B = myEqSys->getMatrix(iB);
#ifdef HAS_SELPC
  // To interface SLEPC another interface is used
  bool ok = eig::solve(A,B,eigVal,eigVec,nev);
#else
  bool ok = eig::solve(A,B,eigVal,eigVec,nev,ncv,iop,shift);
#endif

  // Expand eigenvectors to DOF-ordering and print out eigenvalues
  bool freq = iop == 3 || iop == 4 || iop == 6;
  std::cout <<"\n >>> Computed Eigenvalues <<<\n     Mode\t"
	    << (freq ? "Frequency [Hz]" : "Eigenvalue");
  solution.resize(nev);
  for (int i = 1; i <= nev && ok; i++)
  {
    solution[i-1].eigNo = i;
    if (!mySam->expandVector(eigVec.getColumn(i),solution[i-1].eigVec))
      ok = false;
    else if (!freq)
      solution[i-1].eigVal = eigVal(i);
    else if (eigVal(i) < 0.0)
      solution[i-1].eigVal = -sqrt(-eigVal(i))*0.5/M_PI;
    else
      solution[i-1].eigVal =  sqrt( eigVal(i))*0.5/M_PI;

    std::cout <<"\n     "<< i <<"\t\t"<< utl::trunc(solution[i-1].eigVal);
  }
  std::cout << std::endl;
  return ok;
}


bool SIMbase::writeGlv (const char* inpFile, const int* nViz, int format)
{
  int nBlock = 0;
  return this->writeGlvG(nViz,nBlock,inpFile,format);
}


bool SIMbase::writeGlvG (const int* nViz, int& nBlock,
			 const char* inpFile, int format)
{
  if (inpFile)
  {
    if (myVtf) return false;

#if HAS_VTFAPI == 2
    const char* ext = ".vtfx";
#else
    const char* ext = ".vtf";
#endif

    // Open a new VTF-file
    char* vtfName = new char[strlen(inpFile)+10];
    strtok(strcpy(vtfName,inpFile),".");
    if (nProc > 1)
      sprintf(vtfName+strlen(vtfName),"_p%04d%s",myPid,ext);
    else
      strcat(vtfName,ext);

    std::cout <<"\nWriting VTF-file "<< vtfName << std::endl;
    myVtf = new VTF(vtfName,format);
    delete[] vtfName;
  }
  else
    myVtf->clearGeometryBlocks();

  // Convert and write model geometry
  char pname[16];
  for (size_t i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 0)
      std::cout <<"Writing geometry for patch "<< i+1 << std::endl;
    size_t nd = myModel[i]->getNoParamDim();
    ElementBlock* lvb = new ElementBlock(nd == 3 ? 8 : (nd == 2 ? 4 : 2));
    if (!myModel[i]->tesselate(*lvb,nViz))
      return false;

    sprintf(pname,"Patch %ld",i+1);
    if (!myVtf->writeGrid(lvb,pname,++nBlock))
      return false;
  }

  // Do not write the geometry blocks to file yet, writeVectors might create
  // an additional block for the point vectors, see method VTF::writeVectors
  return true;
}


bool SIMbase::writeGlvBC (const int* nViz, int& nBlock, int iStep) const
{
  if (!myVtf) return false;

  Matrix field;
  size_t i, j;
  int node, geomID = 0;
  std::vector<int> dID[3];

  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    geomID++;
    size_t nbc = myModel[i]->getNoFields(1);
    Matrix bc(nbc,myModel[i]->getNoNodes());
    RealArray flag(3,0.0);
    ASMbase::BCVec::const_iterator bit;
    for (bit = myModel[i]->begin_BC(); bit != myModel[i]->end_BC(); bit++)
      if ((node = myModel[i]->getNodeIndex(bit->node)))
      {
	if (!bit->CX && nbc > 0) bc(1,node) = flag[0] = 1.0;
	if (!bit->CY && nbc > 1) bc(2,node) = flag[1] = 1.0;
	if (!bit->CZ && nbc > 2) bc(3,node) = flag[2] = 1.0;
      }

    if (flag[0]+flag[1]+flag[2] == 0.0)
      continue; // nothing on this patch

    if (msgLevel > 1)
      std::cout <<"Writing boundary conditions for patch "<< i+1 << std::endl;

    if (!myModel[i]->evalSolution(field,bc,nViz))
      return false;

    // The BC fields should either be 0.0 or 1.0
    if (nViz[0] > 2 || nViz[1] > 2 || nViz[2] > 2)
      for (j = 1; j <= 3; j++)
	if (flag[j-1] == 1.0)
	  for (size_t n = 1; n <= field.cols(); n++)
	    if (field(j,n) < 0.9999) field(j,n) = 0.0;

    for (j = 0; j < 3; j++)
      if (flag[j] == 1.0)
	if (myVtf->writeNres(field.getRow(1+j),++nBlock,geomID))
	  dID[j].push_back(nBlock);
	else
	  return false;
  }

  const char* label[3] = {
    "fix_X", "fix_Y", "fix_Z"
  };

  for (j = 0; j < 3; j++)
    if (!dID[j].empty())
      if (!myVtf->writeSblk(dID[j],label[j],1+j,iStep))
	return false;

  return true;
}


bool SIMbase::writeGlvT (int iStep, int& nBlock) const
{
  if (myProblem->hasTractionValues())
  {
    if (msgLevel > 1)
      std::cout <<"Writing boundary tractions"<< std::endl;
    return myProblem->writeGlvT(myVtf,iStep,nBlock);
  }

  return true;
}


bool SIMbase::writeGlvV (const Vector& vec, const char* fieldName,
			 const int* nViz, int iStep, int& nBlock) const
{
  if (vec.empty())
    return true;
  else if (!myVtf)
    return false;

  Matrix field;
  Vector lovec;
  int geomID = 0;
  std::vector<int> vID;

  for (size_t i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing vector field for patch "<< i+1 << std::endl;

    myModel[i]->extractNodeVec(vec,lovec);
    if (!myModel[i]->evalSolution(field,lovec,nViz))
      return false;

    if (!myVtf->writeVres(field,++nBlock,++geomID,this->getNoSpaceDim()))
      return false;
    else
      vID.push_back(nBlock);
  }

  int idBlock = 2;
  return myVtf->writeVblk(vID,fieldName,idBlock,iStep);
}


bool SIMbase::writeGlvS (const Vector& psol, const int* nViz,
			 int iStep, int& nBlock, double time,
			 char psolOnly, const char* pvecName,
			 int idBlock, int psolComps)
{
  if (psol.empty())
    return true;
  else if (!myVtf)
    return false;

  if (!psolOnly && myProblem)
    myProblem->initResultPoints(time);

  Matrix field;
  Vector lovec;
  size_t i, j;
  int geomID = 0;
  const size_t pMAX = 6;
  const size_t sMAX = 21;
  std::vector<int> vID[2], dID[pMAX], sID[sMAX];
  bool scalarEq = this->getNoFields() == 1;
  size_t nVcomp = scalarEq ? 1 : this->getNoSpaceDim();
  bool haveAsol = false;
  if (mySol)
    if (scalarEq)
      haveAsol = mySol->hasScalarSol() > 1;
    else
      haveAsol = mySol->hasVectorSol() > 1;

  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing FE solution for patch "<< i+1 << std::endl;

    // 1. Evaluate primary solution variables

    Vector& locvec = myProblem ? myProblem->getSolution() : lovec;
    myModel[i]->extractNodeVec(psol,locvec,psolComps);
    if (!myModel[i]->evalSolution(field,locvec,nViz))
      return false;

    if (psolOnly > 1)
      geomID++;
    else if (!myVtf->writeVres(field,++nBlock,++geomID,nVcomp))
      return false;
    else
      vID[0].push_back(nBlock);

    for (j = 0; j < field.rows() && j < pMAX; j++)
      if (!myVtf->writeNres(field.getRow(1+j),++nBlock,geomID))
	return false;
      else
	dID[j].push_back(nBlock);

    if (psolOnly || !myProblem) continue; // skip secondary solution

    // 2. Direct evaluation of secondary solution variables

    LocalSystem::patch = i;
    if (!myModel[i]->evalSolution(field,*myProblem,nViz))
      return false;

    if (scalarEq && field.rows() == this->getNoSpaceDim())
      if (!myVtf->writeVres(field,++nBlock,geomID))
	return false;
      else
	vID[1].push_back(nBlock);

    size_t k = 0;
    for (j = 1; j <= field.rows() && k < sMAX; j++)
      if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
	return false;
      else
	sID[k++].push_back(nBlock);

    if (discretization == ASM::Spline || discretization == ASM::SplineC1)
    {
      // 3. Projection of secondary solution variables (tensorial splines only)

      if (!myModel[i]->evalSolution(field,*myProblem,nViz,true))
	return false;

      for (j = 1; j <= field.rows() && k < sMAX; j++)
	if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
	  return false;
	else
	  sID[k++].push_back(nBlock);
    }

    if (haveAsol)
    {
      // 4. Evaluate analytical solution variables

      if (msgLevel > 1)
	std::cout <<"Writing exact solution for patch "<< i+1 << std::endl;

      const ElementBlock* grid = myVtf->getBlock(geomID);
      std::vector<Vec3>::const_iterator cit = grid->begin_XYZ();
      Vector solPt; field.fill(0.0);
      bool ok = true;
      for (j = 1; cit != grid->end_XYZ() && ok; j++, cit++)
      {
	if (mySol->hasScalarSol() == 3 || mySol->hasVectorSol() == 3)
	  ok = myProblem->evalSol(solPt,*mySol->getStressSol(),*cit);
	else if (scalarEq)
	  ok = myProblem->evalSol(solPt,*mySol->getScalarSecSol(),*cit);
	else
	  ok = myProblem->evalSol(solPt,*mySol->getVectorSecSol(),*cit);
	if (ok)
	  field.fillColumn(j,solPt);
      }

      for (j = 1; j <= field.rows() && k < sMAX && ok; j++)
	if (!myVtf->writeNres(field.getRow(j),++nBlock,geomID))
	  return false;
	else
	  sID[k++].push_back(nBlock);
    }
  }

  // Write result block identifications

  bool ok = true;
  if (!vID[0].empty())
    if (pvecName)
      ok = myVtf->writeVblk(vID[0],pvecName,idBlock,iStep);
    else
      ok = myVtf->writeDblk(vID[0],"Solution",idBlock,iStep);

  if (ok && !vID[1].empty())
    ok = myVtf->writeVblk(vID[1],"Flux",idBlock+1,iStep);

  std::string pname;
  if (!myProblem)
  {
    pname = pvecName ? pvecName : "Solution";
    if (psolOnly < 2) pname += "_w";
  }

  for (j = 0; j < pMAX && !dID[j].empty() && ok; j++)
  {
    if (myProblem)
      pname = myProblem->getField1Name(j);
    else if (psolOnly < 2)
      (*pname.rbegin()) ++;
    ok = myVtf->writeSblk(dID[j],pname.c_str(),++idBlock,iStep);
  }

  if (psolOnly || !myProblem || !ok) return ok;

  size_t nf = myProblem->getNoFields(2);
  for (i = j = 0; i < nf && !sID[j].empty(); i++, j++)
    if (!myVtf->writeSblk(sID[j],myProblem->getField2Name(i,haveAsol?"FE":0),
			  ++idBlock,iStep)) return false;

  if (discretization == ASM::Spline || discretization == ASM::SplineC1)
    for (i = 0; i < nf && !sID[j].empty(); i++, j++)
      if (!myVtf->writeSblk(sID[j],myProblem->getField2Name(i,"Projected"),
			    ++idBlock,iStep)) return false;

  if (haveAsol)
    for (i = 0; i < nf && !sID[j].empty(); i++, j++)
      if (!myVtf->writeSblk(sID[j],myProblem->getField2Name(i,"Exact"),
			    ++idBlock,iStep)) return false;

  return true;
}


bool SIMbase::writeGlvP (const Vector& ssol, const int* nViz,
			 int iStep, int& nBlock, double time,
			 int idBlock, const char* prefix)
{
  if (ssol.empty())
    return true;
  else if (!myVtf)
    return false;

  Matrix field;
  Vector lovec;
  const size_t nf = myProblem->getNoFields(2);
  std::vector<int> sID[nf];

  size_t i, j;
  int geomID = 0;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing projected solution for patch "<< i+1 << std::endl;

    // Evaluate the solution variables at the visualization points
    myModel[i]->extractNodeVec(ssol,lovec,nf);
    if (!myModel[i]->evalSolution(field,lovec,nViz))
      return false;

    // Write out to VTF-file as scalar fields
    geomID++;
    for (j = 0; j < field.rows(); j++)
      if (!myVtf->writeNres(field.getRow(1+j),++nBlock,geomID))
	return false;
      else
	sID[j].push_back(nBlock);
  }

  // Write result block identifications
  for (j = 0; j < nf && !sID[j].empty(); j++)
    if (!myVtf->writeSblk(sID[j],myProblem->getField2Name(j,prefix),
			  ++idBlock,iStep)) return false;

  return true;
}


bool SIMbase::writeGlvStep (int iStep, double value, int itype)
{
  myVtf->writeGeometryBlocks(iStep);

  if (itype == 0)
    return myVtf->writeState(iStep,"Time %g",value,itype);
  else
    return myVtf->writeState(iStep,"Step %g",value,itype);
}


bool SIMbase::writeGlvM (const Mode& mode, bool freq,
			 const int* nViz, int& nBlock)
{
  if (mode.eigVec.empty())
    return true;
  else if (!myVtf)
    return false;

  if (msgLevel > 1)
    std::cout <<"Writing eigenvector for Mode "<< mode.eigNo << std::endl;

  Vector displ;
  Matrix field;
  int geomID = 0;
  std::vector<int> vID;

  for (size_t i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches
    if (myModel.size() > 1 && msgLevel > 1) std::cout <<"."<< std::flush;

    geomID++;
    myModel[i]->extractNodeVec(mode.eigVec,displ);
    if (!myModel[i]->evalSolution(field,displ,nViz))
      return false;

    if (!myVtf->writeVres(field,++nBlock,geomID))
      return false;
    else
      vID.push_back(nBlock);
  }
  if (myModel.size() > 1 && msgLevel > 1) std::cout << std::endl;

  int idBlock = 10;
  if (!myVtf->writeDblk(vID,"Mode Shape",idBlock,mode.eigNo))
    return false;

  return myVtf->writeState(mode.eigNo, freq ? "Frequency %g" : "Eigenvalue %g",
			   mode.eigVal, 1);
}


bool SIMbase::writeGlvN (const Matrix& norms, int iStep, int& nBlock)
{
  if (norms.empty())
    return true;
  else if (!myVtf)
    return false;

  Matrix field;
  int geomID = 0;
  std::vector<int> sID[10];

  size_t i, j, k;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (msgLevel > 1)
      std::cout <<"Writing element norms for patch "<< i+1 << std::endl;

    const ElementBlock* grid = myVtf->getBlock(++geomID);
    myModel[i]->extractElmRes(norms,field);

    if (grid->getNoElms() > field.cols())
    {
      // Expand the element result array
      Matrix efield(field);
      field.resize(field.rows(),grid->getNoElms());
      for (j = 1; j <= field.cols(); j++)
	field.fillColumn(j,efield.getColumn(grid->getElmId(j)));
    }

    for (j = k = 0; j < field.rows() && j < 10; j++)
      if (j != 1) // Norm #1 (external norm) does not have element contributions
	if (!myVtf->writeEres(field.getRow(1+j),++nBlock,geomID))
	  return false;
	else
	  sID[k++].push_back(nBlock);
  }

  int idBlock = 200;
  for (k = 0; !sID[k].empty(); k++)
    if (!myVtf->writeSblk(sID[k],NormBase::getName(k),++idBlock,iStep,true))
      return false;

  return true;
}


void SIMbase::closeGlv ()
{
  if (myVtf) delete myVtf;
  myVtf = 0;
}


void SIMbase::dumpPrimSol (const Vector& psol, std::ostream& os,
			   bool withID) const
{
  if (psol.empty()) return;

  size_t i, j, ip;
  unsigned char k, n;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    Vector patchSol;
    myModel[i]->extractNodeVec(psol,patchSol);

    if (withID)
    {
      if (myModel.size() > 1)
	os <<"\n# Patch: "<< i+1;
      os <<"\n# inod/gnod\tNodal Coordinates\tSolution\n";
    }
    for (ip = 0, j = 1; j <= myModel[i]->getNoNodes(); j++)
    {
      if ((n = myModel[i]->getNodalDOFs(j)) == 0)
	continue;
      else if (withID)
	os << j <<' '<< myModel[i]->getNodeID(j)
	   <<"\t\t"<< myModel[i]->getCoord(j) <<"\t\t";
      os << utl::trunc(patchSol[ip++]);
      for (k = 1; k < n; k++)
	os <<' '<< utl::trunc(patchSol[ip++]);
      os <<'\n';
    }
  }

  os.flush();
}


bool SIMbase::dumpGeometry (std::ostream& os) const
{
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel[i]->empty())
      if (!myModel[i]->write(os))
	return false;

  return true;
}


bool SIMbase::dumpBasis (std::ostream& os, int basis, size_t patch) const
{
  size_t start = patch ? patch-1 : 0;
  size_t end = patch ? start+1 : myModel.size();
  for (size_t i = start; i < end && i < myModel.size(); i++)
    if (!myModel[i]->empty())
      if (!myModel[i]->write(os,basis))
	return false;

  return true;
}


bool SIMbase::dumpSolution (const Vector& psol, std::ostream& os) const
{
  if (psol.empty())
    return true;

  Matrix field;
  size_t i, j, k;

  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    if (myModel.size() > 1)
      os <<"\n# Patch: "<< i+1;

    // Extract and write primary solution
    size_t nf = myModel[i]->getNoFields(1);
    Vector& patchSol = myProblem->getSolution();
    myModel[i]->extractNodeVec(psol,patchSol);
    for (k = 1; k <= nf; k++)
    {
      os << myProblem->getField1Name(k,"# FE");
      for (j = 1; j <= myModel[i]->getNoNodes(); j++)
      {
	std::pair<int,int> dofs = mySam->getNodeDOFs(j);
	int idof = dofs.first+k-1;
	if (idof <= dofs.second)
	  os <<"\n"<< utl::trunc(patchSol[idof-1]);
      }
      os << std::endl;
    }

    // Evaluate secondary solution variables
    LocalSystem::patch = i;
    if (!myModel[i]->evalSolution(field,*myProblem))
      return false;

    // Write the secondary solution
    for (j = 1; j <= field.rows(); j++)
    {
      os << myProblem->getField2Name(j-1,"# FE");
      for (k = 1; k <= field.cols(); k++)
	os <<"\n"<< utl::trunc(field(j,k));
      os << std::endl;
    }
  }

  return true;
}


bool SIMbase::dumpResults (const Vector& psol, double time, std::ostream& os,
			   bool formatted, std::streamsize precision) const
{
  if (psol.empty() || myPoints.empty())
    return true;

  myProblem->initResultPoints(time);

  size_t i, j, k;
  Matrix sol1, sol2;
  Vector reactionFS;
  const Vector* reactionForces = myEqSys->getReactions();

  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    ResPointVec::const_iterator p;
    std::vector<int> points;
    RealArray params[3];

    // Find all evaluation points within this patch, if any
    for (j = 0, p = myPoints.begin(); p != myPoints.end(); j++, p++)
      if (this->getLocalPatchIndex(p->patch) == (int)(i+1))
	if (discretization >= ASM::Spline)
	{
	  points.push_back(p->inod > 0 ? p->inod : -(j+1));
	  for (k = 0; k < myModel[i]->getNoParamDim(); k++)
	    params[k].push_back(p->par[k]);
	}
	else if (p->inod > 0)
	  points.push_back(p->inod);

    if (points.empty()) continue; // no points in this patch

    myModel[i]->extractNodeVec(psol,myProblem->getSolution());
    if (discretization >= ASM::Spline)
    {
      // Evaluate the primary solution variables
      if (!myModel[i]->evalSolution(sol1,myProblem->getSolution(),params,false))
	return false;

      // Evaluate the secondary solution variables
      LocalSystem::patch = i;
      if (myProblem->getNoFields(2) > 0)
        if (!myModel[i]->evalSolution(sol2,*myProblem,params,false))
          return false;
    }
    else
      // Extract nodal primary solution variables
      if (!myModel[i]->getSolution(sol1,myProblem->getSolution(),points))
	return false;

    // Formatted output, use scientific notation with fixed field width
    std::streamsize flWidth = 8 + precision;
    std::streamsize oldPrec = os.precision(precision);
    std::ios::fmtflags oldF = os.flags(std::ios::scientific | std::ios::right);
    for (j = 0; j < points.size(); j++)
    {
      if (!formatted)
        os << time <<" ";
      else if (points[j] < 0)
	os <<"  Point #"<< -points[j] <<":\tsol1 =";
      else
      {
	points[j] = myModel[i]->getNodeID(points[j]);
	os <<"  Node #"<< points[j] <<":\tsol1 =";
      }

      for (k = 1; k <= sol1.rows(); k++)
	os << std::setw(flWidth) << utl::trunc(sol1(k,j+1));

      if (discretization >= ASM::Spline)
      {
        if (formatted && sol2.rows() > 0)
          os <<"\n\t\tsol2 =";
	for (k = 1; k <= sol2.rows(); k++)
	  os << std::setw(flWidth) << utl::trunc(sol2(k,j+1));
      }

      if (reactionForces && points[j] > 0)
	// Print nodal reaction forces for nodes with prescribed DOFs
	if (mySam->getNodalReactions(points[j],*reactionForces,reactionFS))
	{
          if (formatted)
            os <<"\n\t\treac =";
	  for (k = 0; k < reactionFS.size(); k++)
	    os << std::setw(flWidth) << utl::trunc(reactionFS[k]);
	}

      os << std::endl;
    }
    os.precision(oldPrec);
    os.flags(oldF);
  }

  return true;
}


bool SIMbase::project (Matrix& ssol, const Vector& psol,
		       ProjectionMethod pMethod) const
{
  PROFILE1("Solution projection");

  if (msgLevel > 1)
    std::cout <<"\nProjecting secondary solution ...\n"<< std::endl;

  ssol.resize(0,0);

  size_t i, j, n;
  size_t ngNodes = mySam->getNoNodes();

  Matrix values;
  Vector count(myModel.size() > 1 ? ngNodes : 0);

  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    // Extract the primary solution control point values for this patch
    myModel[i]->extractNodeVec(psol,myProblem->getSolution());

    // Project the secondary solution and retrieve control point values
    switch (pMethod) {
    case GLOBAL:
      if (!myModel[i]->evalSolution(values,*myProblem))
	return false;
      break;

    case LOCAL:
      // Annette, add your local projection stuff here...

    default:
      std::cerr <<" *** SIMbase::project: Projection method "<< pMethod
		<<" not implemented."<< std::endl;
      return false;
    }

    size_t nComps = values.rows();
    size_t nNodes = myModel[i]->getNoNodes();
    if (ssol.empty())
      ssol.resize(nComps,ngNodes);

    // Nodal averaging for nodes referred to by two or more patches
    // (these are typically the interface nodes)
    for (n = 1; n <= nNodes; n++)
      if (count.empty())
	ssol.fillColumn(myModel[i]->getNodeID(n),values.getColumn(n));
      else
      {
	int inod = myModel[i]->getNodeID(n);
	for (j = 1; j <= nComps; j++)
	  ssol(j,inod) += values(j,n);
	count(inod) ++;
      }
  }

  // Divide through by count(n) to get the nodal average at the interface nodes
  for (n = 1; n <= count.size(); n++)
    if (count(n) > 1.0)
      for (j = 1; j <= ssol.rows(); j++)
	ssol(j,n) /= count(n);

  return true;
}


bool SIMbase::project (Vector& sol) const
{
  Matrix secsol;
  if (!this->project(secsol,sol))
    return false;

  sol = secsol;
  return true;
}


size_t SIMbase::extractPatchSolution (const Vector& sol, int pindx)
{
  if (pindx < 0 || (size_t)pindx >= myModel.size() || sol.empty())
    return 0;

  myModel[pindx]->extractNodeVec(sol,myProblem->getSolution());
  return myModel[pindx]->getNoFields(1)*myModel[pindx]->getNoNodes(1);
}


bool SIMbase::injectPatchSolution (Vector& sol, const Vector& locSol,
				   int pindx, unsigned char nndof)
{
  if (pindx < 0 || (size_t)pindx >= myModel.size())
    return false;

  return myModel[pindx]->injectNodeVec(locSol,sol,nndof);
}


bool SIMbase::evalSecondarySolution (Matrix& field, int pindx) const
{
  if (pindx < 0 || (size_t)pindx >= myModel.size())
    return false;

  return myModel[pindx]->evalSolution(field,*myProblem);
}


bool SIMbase::extractPatchElmRes (const Matrix& globRes, Matrix& elmRes,
				  int pindx)
{
  if (pindx < 0 || (size_t)pindx >= myModel.size() || globRes.empty())
    return false;

  myModel[pindx]->extractElmRes(globRes,elmRes);
  return true;
}
