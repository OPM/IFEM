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
#include "SIMoptions.h"
#include "ASMs2DC1.h"
#ifdef HAS_PETSC
#include "SAMpatchPETSc.h"
#else
#include "SAMpatch.h"
#endif
#include "IntegrandBase.h"
#include "AlgEqSystem.h"
#include "LinSolParams.h"
#include "EigSolver.h"
#include "GlbNorm.h"
#include "ElmNorm.h"
#include "AnaSol.h"
#include "TensorFunction.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include "Profiler.h"
#include "IFEM.h"
#include <fstream>
#ifdef SP_DEBUG
#include <cassert>
#endif


bool SIMbase::preserveNOrder  = false;
bool SIMbase::ignoreDirichlet = false;


SIMbase::SIMbase (IntegrandBase* itg) : g2l(&myGlb2Loc)
{
  isRefined = false;
  nsd = 3;
  myProblem = itg;
  mySol = nullptr;
  myEqSys = nullptr;
  mySam = nullptr;
  mySolParams = nullptr;
  nGlPatches = 0;
  nIntGP = nBouGP = 0;
  lagMTOK = false;
  extEnergy = 0.0;

  MPCLess::compareSlaveDofOnly = true; // to avoid multiple slave definitions
}


SIMbase::~SIMbase ()
{
#ifdef SP_DEBUG
  IFEM::cout <<"\nEntering SIMbase destructor"<< std::endl;
#endif

  for (auto& itg : myInts)
    if (itg.second != myProblem)
      delete itg.second;

  if (myProblem)   delete myProblem;
  if (mySol)       delete mySol;
  if (myEqSys)     delete myEqSys;
  if (mySam)       delete mySam;
  if (mySolParams) delete mySolParams;

  for (ASMbase* patch : myModel)
    delete patch;

  myModel.clear();
  this->SIMbase::clearProperties();

#ifdef SP_DEBUG
  IFEM::cout <<"Leaving SIMbase destructor"<< std::endl;
#endif
}


void SIMbase::clearProperties ()
{
  for (ASMbase* patch : myModel)
    patch->clear(true); // retain the geometry only

  for (auto& i2 : myScalars)
    delete i2.second;
  for (auto& i3 : myVectors)
    delete i3.second;
  for (auto& i4 : myTracs)
    delete i4.second;

  myPatches.clear();
  myGlb2Loc.clear();
  myScalars.clear();
  myVectors.clear();
  myTracs.clear();
  myProps.clear();
  myInts.clear();
  mixedMADOFs.clear();
  adm.dd.setElms({},"");
}


int SIMbase::getLocalPatchIndex (int patchNo) const
{
  if (patchNo < 1 || (patchNo > nGlPatches && nGlPatches > 0))
  {
    std::cerr <<" *** SIMbase::getLocalPatchIndex: Patch number "<< patchNo
              <<" is out of range [1,"<< nGlPatches <<"]."<< std::endl;
    return -1;
  }
  else if (myPatches.empty() || nProc == 1)
    return patchNo;

  for (size_t i = 0; i < myPatches.size(); i++)
    if (myPatches[i] == patchNo)
      return 1+i;

  return 0;
}


ASMbase* SIMbase::getPatch (int idx, bool glbIndex) const
{
  int pid = glbIndex ? this->getLocalPatchIndex(idx) : idx;
  return pid > 0 && (size_t)pid <= myModel.size() ? myModel[pid-1] : nullptr;
}


#if SP_DEBUG > 2
static void printNodalConnectivity (const ASMVec& model, std::ostream& os)
{
  typedef std::pair<int,int> Ipair;
  typedef std::vector<Ipair> Ipairs;
  std::map<int,Ipairs> nodeInfo;
  for (ASMbase* pch : model)
    if (!pch->empty())
      for (size_t n = 1; n <= pch->getNoNodes(); n++)
        nodeInfo[pch->getNodeID(n)].push_back(std::make_pair(pch->idx,n));

  for (const std::pair<int,Ipairs>& node : nodeInfo)
    if (node.second.size() > 1)
    {
      os <<"\nConnectivity for node "<< node.first <<":";
      for (const Ipair& n : node.second)
        os <<" P"<< n.first+1 <<","<< n.second;
    }

  std::cout << std::endl;
}
#endif


bool SIMbase::preprocess (const IntVec& ignored, bool fixDup)
{
  if (myModel.empty())
    return true; // Empty simulator, nothing to preprocess

  if (mySam && !isRefined)
  {
    std::cerr <<" *** SIMbase::preprocess: Logic error, invoked more than once"
              <<" for "<< (myHeading.empty() ? this->getName() : myHeading)
              <<"."<< std::endl;
    return false;
  }

  PROFILE1("Model preprocessing");

  static int substep = 10;
  this->printHeading(substep);

  // Perform some sub-class specific pre-preprocessing, if any
  this->preprocessA();

  // Create the classical FE data structures
  if (!this->createFEMmodel('Y')) return false;

  PatchVec::const_iterator mit;
  ASMbase* pch;
  size_t patch;

  // Erase all patches that should be ignored in the analysis
  for (int idx : ignored)
    if ((pch = this->getPatch(idx)))
      pch->clear();

  // If material properties are specified for at least one patch, assign the
  // property code 999999 to all patches with no material property code yet
  PatchVec pchWthMat;
  for (const Property& p : myProps)
    if (p.pcode == Property::MATERIAL && (pch = this->getPatch(p.patch)))
      if (!pch->empty()) pchWthMat.push_back(pch);

  if (!pchWthMat.empty())
    for (mit = myModel.begin(), patch = 1; mit != myModel.end(); ++mit, patch++)
      if (std::find(pchWthMat.begin(),pchWthMat.end(),*mit) == pchWthMat.end())
	myProps.push_back(Property(Property::MATERIAL,999999,patch,
				   (*mit)->getNoParamDim()));

  if (fixDup)
  {
    // Check for duplicated nodes (missing topology)
    int nDupl = 0;
    std::map<Vec3,int> globalNodes;
    for (mit = myModel.begin(), patch = 1; mit != myModel.end(); ++mit, patch++)
      if (!(*mit)->empty())
      {
	IFEM::cout <<"   * Checking Patch "<< patch << std::endl;
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
      IFEM::cout <<"   * "<< nDupl <<" duplicated nodes merged."<< std::endl;
  }

#if SP_DEBUG > 2
  printNodalConnectivity(myModel,std::cout);
#endif

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
  else for (mit = myModel.begin(); mit != myModel.end(); ++mit)
    renum += (*mit)->renumberNodes(myGlb2Loc,ngnod);

  if (renum > 0)
    IFEM::cout <<"\nRenumbered "<< renum <<" nodes."<< std::endl;

  // Apply the old-to-new node number mapping to all node tables in the model
  if (!this->renumberNodes(*g2l))
    return false;

  // Perform specialized preprocessing before the assembly initialization.
  // This typically involves the system-level Lagrange multipliers, etc.
  if (!this->preprocessBeforeAsmInit(ngnod))
    return false;

  IFEM::cout <<"\nResolving Dirichlet boundary conditions"<< std::endl;
  ASMbase::resetNumbering(ngnod); // to account for possibly added nodes

  // Process the Dirichlet boundary conditions in the order of increasing
  // dimension, such that vertex definitions override definitions on edges,
  // and edge definitions override definitions on faces
  size_t nprop = 0;
  int code, dofs, ierr = 0, iwar = 0;
  for (unsigned char dim = 0; nprop < myProps.size(); dim++)
    for (const Property& p : myProps)
      if (abs(p.ldim) == dim || (dim == 2 && abs(p.ldim) > 3))
      {
        nprop++;
        code = p.pindx;
        dofs = abs(code%1000000);
        switch (p.pcode) {
        case Property::DIRICHLET:
          code = 0;
        case Property::DIRICHLET_INHOM:
          break;

        case Property::UNDEFINED:
          ++iwar;
#ifdef SP_DEBUG
          std::cout <<"  ** SIMbase::preprocess: Undefined property set, code="
                    << p.pindx <<" Patch="<< p.patch <<" Item="
                    << (int)p.lindx <<" "<< (int)p.ldim <<"D"<< std::endl;
#endif
        default:
          dofs = 0;
          break;
        }

        if (dofs > 0)
          if (this->addConstraint(p.patch,p.lindx,p.ldim,dofs,code,ngnod,
                                  p.basis))
            IFEM::cout << std::endl;
          else
            ++ierr;
      }

  if (iwar > 0)
    std::cerr <<"\n  ** SIMbase::preprocess: Warning: "<< iwar
              <<" undefined property sets were detected.\n";
  if (ierr > 0)
    std::cerr <<"\n *** SIMbase::preprocess: Error: "<< ierr
              <<" invalid Dirichlet properties were detected.\n"
              <<"     Please check your model, execution aborts..."
              << std::endl;
  else if (iwar > 0)
    std::cerr <<"     Please verify your model, execution continues..."
              << std::endl;

  // Compute the set of all MPCs over the whole model. This will also merge
  // multiple constraint equations defined on interfaces between patches.
  MPCSet allMPCs;
  ASMbase::mergeAndGetAllMPCs(myModel,allMPCs);

  // Set initial values for the inhomogeneous dirichlet conditions, if any,
  // and compute coupling coefficients for the C1-continuity constraints
  if (!this->initDirichlet()) return false;

  // Resolve possibly chaining of the MPC equations
  if (!allMPCs.empty())
    ASMbase::resolveMPCchains(allMPCs,this->hasTimeDependentDirichlet());

  // Generate element groups for multi-threading
  bool silence = msgLevel < 1 || (msgLevel < 3 && nGlPatches > 1);
  for (mit = myModel.begin(); mit != myModel.end() && myProblem; ++mit)
    if (!(*mit)->empty())
      (*mit)->generateThreadGroups(*myProblem,silence,lagMTOK);

  for (const Property& p : myProps)
    if (p.pcode == Property::NEUMANN ||
        p.pcode == Property::NEUMANN_GENERIC ||
        p.pcode == Property::ROBIN)
      this->generateThreadGroups(p,silence);

  // Preprocess the result points
  this->preprocessResultPoints();

  // Check if the integrand is monolithic with different nodal DOF types
  std::vector<char> dofTypes;
  if (myProblem)
    myProblem->getNodalDofTypes(dofTypes);

  // Initialize data structures for the algebraic system
  if (mySam) delete mySam;
#ifdef HAS_PETSC
  if (opt.solver == SystemMatrix::PETSC)
    mySam = new SAMpatchPETSc(*g2l,adm);
  else
    mySam = new SAMpatch();
#else
  mySam = new SAMpatch();
#endif
  if (!static_cast<SAMpatch*>(mySam)->init(myModel,ngnod,dofTypes))
    return false;

  if (nProc > 1 && myPatches.empty() && adm.isParallel() && adm.dd.getElms().empty())
  {
    IFEM::cout <<" *** SIMbase::preprocess: No partitioning information for "
               << nProc <<" processors found. Using graph partitioning\n";
  }

  if (!adm.dd.setup(adm,*this))
  {
    std::cerr <<"\n *** SIMbase::preprocess(): Failed to establish "
              <<" domain decomposition data."<< std::endl;
    return false;
  }
  else if (adm.dd.isPartitioned()) // reuse thread groups for element list
  {
    for (ASMbase* pch : this->getFEModel())
      pch->generateThreadGroupsFromElms(adm.dd.getElms());
  }

  if (!myProblem)
  {
    std::cerr <<"\n *** SIMbase::preprocess(): No problem integrand for the "
              << this->getName() <<"-simulator."<< std::endl;
    return false;
  }

  // Now perform the sub-class specific final preprocessing, if any
  return this->preprocessB() && ierr == 0;
}


bool SIMbase::renumberNodes (const std::map<int,int>& nodeMap)
{
  bool ok = true;
  for (ASMbase* pch : myModel)
    ok &= pch->renumberNodes(nodeMap);

  ASMs2DC1::renumberNodes(nodeMap);
  return ok;
}


void SIMbase::generateThreadGroups (const Property& p, bool silence)
{
  ASMbase* pch = this->getPatch(p.patch);
  if (pch && !pch->empty() && abs(p.ldim)+1 == pch->getNoParamDim())
    pch->generateThreadGroups(p.lindx%10,silence,lagMTOK);
}


VecFunc* SIMbase::getVecFunc (size_t patch, Property::Type ptype) const
{
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (p->patch == patch)
      if (p->pcode == ptype || ptype == Property::UNDEFINED)
      {
        VecFuncMap::const_iterator it = myVectors.find(p->pindx);
        if (it != myVectors.end()) return it->second;
      }

  return nullptr;
}


RealFunc* SIMbase::getSclFunc (int code) const
{
  SclFuncMap::const_iterator it = myScalars.find(code);
  return it == myScalars.end() ? nullptr : it->second;
}


bool SIMbase::initSystem (int mType, size_t nMats, size_t nVec, size_t nScl,
                          bool withRF)
{
  if (!mySam) return true; // Silently ignore when no algebraic system

#if SP_DEBUG > 2
  mySam->print(std::cout);
  std::string heading("\n\nNodal coordinates for Patch 1");
  char& patchNo = heading[heading.size()-1];
  for (ASMbase* pch : myModel)
  {
    if (!pch->empty())
      pch->printNodes(std::cout,heading.c_str());
    ++patchNo;
  }
  printNodalConnectivity(myModel,std::cout);
#endif

  if (myEqSys) delete myEqSys;
  myEqSys = new AlgEqSystem(*mySam,adm);

  // Workaround SuperLU bug for tiny systems
  if (mType == SystemMatrix::SPARSE && this->getNoElms(true) < 3)
  {
    std::cerr <<"  ** System too small for SuperLU, falling back to Dense."
              << std::endl;
    mType = SystemMatrix::DENSE;
  }

  return myEqSys->init(static_cast<SystemMatrix::Type>(mType),
                       mySolParams, nMats, nVec, nScl, withRF,
                       myProblem->getLinearSystemType(), opt.num_threads_SLU);
}


bool SIMbase::setAssociatedRHS (size_t iMat, size_t iVec)
{
  return myEqSys ? myEqSys->setAssociatedVector(iMat,iVec) : false;
}


bool SIMbase::setMode (int mode, bool resetSol)
{
  if (myInts.empty())
    myInts.insert(std::make_pair(0,myProblem));

  for (IntegrandMap::iterator it = myInts.begin(); it != myInts.end(); ++it)
    if (it->second)
    {
      it->second->setMode((SIM::SolutionMode)mode);
      if (resetSol) it->second->resetSolution();
    }
    else
    {
      std::cerr <<" *** SIMbase::setMode: No integrand yet";
      if (it->first > 0) std::cerr <<", code="<< it->first;
      std::cerr <<"."<< std::endl;
      return false;
    }

  return true;
}


void SIMbase::setIntegrationPrm (unsigned short int i, double prm)
{
  if (!myInts.empty())
    for (IntegrandMap::iterator it = myInts.begin(); it != myInts.end(); ++it)
      it->second->setIntegrationPrm(i,prm);
  else if (myProblem)
    myProblem->setIntegrationPrm(i,prm);
  else
    std::cerr <<"  ** SIMbase::setIntegrationPrm: myProblem not set yet."
              << std::endl;
}


void SIMbase::setQuadratureRule (size_t ng, bool redimBuffers, bool printQP)
{
  size_t nInterfaceGP = nIntGP = nBouGP = 0;
  for (size_t i = 0; i < myModel.size(); i++)
    if (!myModel[i]->empty())
    {
      myModel[i]->setGauss(ng);
      // Count the interior integration points
      myModel[i]->getNoIntPoints(nIntGP,nInterfaceGP);
    }

  if (!myProblem) return;

  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (p->pcode == Property::NEUMANN && myProblem->hasBoundaryTerms())
    {
      // Account for possibly more than one Neumann property on a boundary
      bool notCounted = true;
      for (PropertyVec::const_iterator q = myProps.begin(); q != p; ++q)
        if (q->patch == p->patch && q->lindx%10 == p->lindx%10 &&
            q->pcode == p->pcode)
          notCounted = false;

      if (notCounted) // Count the boundary integration points
        myModel[p->patch-1]->getNoBouPoints(nBouGP,abs(p->ldim),p->lindx%10);
    }

  // Let the integrands know how many integration points in total we have
  if (redimBuffers)
    for (IntegrandMap::iterator it = myInts.begin(); it != myInts.end(); ++it)
      it->second->initIntegration(nIntGP,nBouGP);

  if (printQP)
  {
    IFEM::cout <<"Number of quadrature points "<< nIntGP;
    if (nInterfaceGP > 0) IFEM::cout <<" "<< nInterfaceGP;
    if (nBouGP > 0) IFEM::cout <<" "<< nBouGP;
    IFEM::cout << std::endl;
  }
}


void SIMbase::printProblem () const
{
  if (myProblem)
  {
    IFEM::cout <<"\nProblem definition:"<< std::endl;
    myProblem->printLog();
  }

#if SP_DEBUG > 1
  std::cout <<"\nProperty mapping:";
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    std::cout <<"\n"<< p->pcode <<" "<< p->pindx <<" "<< p->patch
              <<" "<< (int)p->lindx <<" "<< (int)p->ldim;
  std::cout << std::endl;
#endif
}


size_t SIMbase::getNoFields (int basis) const
{
  return myModel.empty() ? 0 : myModel.front()->getNoFields(basis);
}


size_t SIMbase::getNoDOFs () const
{
  return mySam ? mySam->getNoDOFs() : 0;
}


size_t SIMbase::getNoNodes (int basis) const
{
  if (mySam)
  {
    if (basis == 1)
      return (mySam->getNoNodes('D') + // Default node type for 1st-basis nodes
              mySam->getNoNodes('X') + // Associate X-nodes with 1st basis
              mySam->getNoNodes('0')); // Associate 0-dof nodes with 1st basis
    else if (basis > 1)
      return mySam->getNoNodes('N'+basis);
    else
      return mySam->getNoNodes();
  }
  else if (myModel.size() == 1)
    return myModel.front()->getNoNodes(basis);

  std::cerr <<" *** SIMbase::getNoNodes: Number of nodes in a multi-patch model"
            <<" is not known at this point."<< std::endl;
  return 0;
}


size_t SIMbase::getNoElms (bool includeXelms) const
{
  if (mySam && includeXelms)
    return mySam->getNoElms();

  size_t noElms = 0;
  for (size_t i = 0; i < myModel.size(); i++)
    noElms += myModel[i]->getNoElms(false,includeXelms);

  return noElms;
}


bool SIMbase::getElmNodes (IntVec& mnpc, int iel) const
{
  if (mySam)
    return mySam->getElmNodes(mnpc,iel);
  else if (myModel.size() != 1)
  {
    std::cerr <<" *** SIMbase::getElmNodes: The nodal point correspondance is"
              <<" unknown at this point in a multi-patch model."<< std::endl;
    return false;
  }

  mnpc = myModel.front()->getElementNodes(iel);
  return !mnpc.empty();
}


size_t SIMbase::getNoSolutions () const
{
  return myProblem ? myProblem->getNoSolutions() : 0;
}


size_t SIMbase::getNoEquations () const
{
  return mySam ? mySam->getNoEquations() : 0;
}


size_t SIMbase::getNoRHS () const
{
  return myEqSys ? myEqSys->getNoRHS() : 1;
}


unsigned char SIMbase::getNoBasis () const
{
  size_t result = myModel.empty() ? 0 : myModel.front()->getNoBasis();
#ifdef SP_DEBUG
  for (size_t i = 1; i < myModel.size(); i++)
    assert(myModel[i]->getNoBasis() == result);
#endif

  return result;
}


bool SIMbase::hasTimeDependentDirichlet () const
{
  for (size_t i = 0; i < myModel.size(); i++)
    if (myModel[i]->hasTimeDependentDirichlet(myScalars,myVectors))
      return true;

  return false;
}


bool SIMbase::initDirichlet (double time)
{
  if (time == 0.0)
    for (size_t i = 0; i < myModel.size(); i++)
      if (!myModel[i]->initConstraints())
	return false;

  Vector dummy;
  return this->updateDirichlet(time,&dummy);
}


bool SIMbase::updateDirichlet (double time, const Vector* prevSol)
{
  if (prevSol)
    for (size_t i = 0; i < myModel.size(); i++)
      if (!myModel[i]->updateDirichlet(myScalars,myVectors,time))
	return false;

  SAMpatch* pSam = dynamic_cast<SAMpatch*>(mySam);
  return pSam ? pSam->updateConstraintEqs(myModel,prevSol) : true;
}


bool SIMbase::updateGrid (const Vector& displ)
{
  if (displ.empty()) return true; // No displacements (yet), totally fine

  bool ok = true;
  Vector locdisp;
  for (size_t i = 0; i < myModel.size() && ok; i++)
  {
    myModel[i]->extractNodeVec(displ,locdisp,myModel[i]->getNoSpaceDim(),-1);
    ok = myModel[i]->updateCoords(locdisp);
  }

  return ok;
}


bool SIMbase::updateGrid (const std::string& field)
{
  const Vector* displ = this->getDependentField(field);
  if (displ) return this->updateGrid(*displ);

  std::cerr <<" *** SIMbase::updateGrid: No such field \""<< field
	    <<"\" registered for \""<< this->getName() <<"\"."<< std::endl;
  return false;
}


void SIMbase::getBoundaryNodes (int pcode, IntVec& glbNodes, Vec3Vec* XYZ) const
{
  glbNodes.clear();
  if (XYZ) XYZ->clear();

  ASMbase* pch;
  size_t node;
  for (PropertyVec::const_iterator p = myProps.begin(); p != myProps.end(); ++p)
    if (abs(p->pindx) == pcode && (pch = this->getPatch(p->patch))) {
      if (abs(p->ldim)+1 == pch->getNoParamDim() || p->ldim == 4) {
        // The boundary is of one dimension lower than the patch
        IntVec nodes;
        if (p->ldim == 4) // The boundary nodes are stored explicitly
          nodes = pch->getNodeSet(p->lindx);
        else
          pch->getBoundaryNodes(abs(p->lindx%10),nodes);
        for (int n : nodes)
          if (std::find(glbNodes.begin(),glbNodes.end(),n) == glbNodes.end()) {
            glbNodes.push_back(n);
            if (XYZ) {
              if ((node = pch->getNodeIndex(n,true)))
                XYZ->push_back(pch->getCoord(node));
              else
                XYZ->push_back(Vec3());
            }
          }
      }
      else if (pch->getNoParamDim() == abs(p->ldim)) {
        // The boundary and the patch are of same dimension
        for (node = 1; node <= pch->getNoNodes(); node++)
        {
          glbNodes.push_back(pch->getNodeID(node));
          if (XYZ) XYZ->push_back(pch->getCoord(node));
        }
      }
    }
}


int SIMbase::findClosestNode (const Vec3& X) const
{
  if (myModel.empty()) return -1;

  ASMbase* closestPch = nullptr;
  std::pair<size_t,double> closest(0,1.0e99);
  for (size_t i = 0; i < myModel.size(); i++)
  {
    std::pair<size_t,double> node = myModel[i]->findClosestNode(X);
    if (node.first > 0 && node.second < closest.second)
    {
      closest = node;
      closestPch = myModel[i];
    }
  }
  if (!closestPch) return -2;

#ifdef SP_DEBUG
  std::cout <<"SIMbase::findClosestNode("<< X <<") -> Node "<< closest.first
            <<" in Patch "<< closestPch->idx+1 <<" distance="<< closest.second
            << std::endl;
#endif

  return closestPch->getNodeID(closest.first);
}


bool SIMbase::assembleSystem (const TimeDomain& time, const Vectors& prevSol,
			      bool newLHSmatrix, bool poorConvg)
{
  PROFILE1("Element assembly");

  bool ok = true;
  bool isAssembling = (myProblem->getMode() > SIM::INIT &&
                       myProblem->getMode() < SIM::RECOVERY);
  if (isAssembling && myEqSys)
    myEqSys->initialize(newLHSmatrix);

  // Loop over the integrands
  IntegrandMap::const_iterator it;
  for (it = myInts.begin(); it != myInts.end() && ok; ++it)
  {
    if (msgLevel > 1)
      IFEM::cout <<"\n\nProcessing integrand associated with code "<< it->first
                << std::endl;

    GlobalIntegral& sysQ = it->second->getGlobalInt(myEqSys);
    if (&sysQ != myEqSys && isAssembling)
      sysQ.initialize(newLHSmatrix);

    if (!prevSol.empty())
      it->second->initIntegration(time,prevSol.front(),poorConvg);

    // Loop over the different material regions, integrating interior
    // coefficient matrix terms for the patch associated with each material
    size_t lp = 0;
    ASMbase* pch = nullptr;
    PropertyVec::const_iterator p, p2;
    if (it->second->hasInteriorTerms())
    {
      for (p = myProps.begin(); p != myProps.end() && ok; ++p)
        if (p->pcode == Property::MATERIAL &&
            (it->first == 0 || it->first == p->pindx))
          if (!(pch = this->getPatch(p->patch)))
          {
            std::cerr <<" *** SIMbase::assembleSystem: Patch index "<< p->patch
                      <<" out of range [1,"<< myModel.size() <<"]."<< std::endl;
            ok = false;
          }
          else if (this->initMaterial(p->pindx))
          {
            lp = p->patch;
            if (msgLevel > 1)
              IFEM::cout <<"\nAssembling interior matrix terms for P"<< lp
                         << std::endl;
            ok &= this->initBodyLoad(lp);
            ok &= this->extractPatchSolution(it->second,prevSol,lp-1);
            if (mySol)
              mySol->initPatch(pch->idx);
#if SP_DEBUG > 2
            it->second->printSolution(std::cout,pch->idx+1);
#endif
            ok &= pch->integrate(*it->second,sysQ,time);
          }
          else
            ok = false;

      if (lp == 0 && it->first == 0)
        // All patches refer to the same material, and we assume it has been
        // initialized during input processing (thus no initMaterial call here)
        for (size_t k = 0; k < myModel.size() && ok; k++)
        {
          lp = k+1;
          if (msgLevel > 1)
            IFEM::cout <<"\nAssembling interior matrix terms for P"<< lp
                       << std::endl;
          ok &= this->initBodyLoad(lp);
          ok &= this->extractPatchSolution(it->second,prevSol,k);
          if (mySol)
            mySol->initPatch(myModel[k]->idx);
#if SP_DEBUG > 2
          it->second->printSolution(std::cout,myModel[k]->idx+1);
#endif
          ok &= myModel[k]->integrate(*it->second,sysQ,time);
        }
    }

    // Assemble contributions from the Neumann boundary conditions
    // and other boundary integrals (Robin properties, contact, etc.)
    if (it->second->hasBoundaryTerms() && myEqSys && myEqSys->getVector())
      for (p = myProps.begin(); p != myProps.end() && ok; ++p)
        if ((p->pcode == Property::NEUMANN && it->first == 0) ||
            ((p->pcode == Property::NEUMANN_GENERIC ||
              p->pcode == Property::ROBIN) && it->first == p->pindx))
        {
          if (!(pch = this->getPatch(p->patch)))
          {
            std::cerr <<" *** SIMbase::assembleSystem: Patch index "<< p->patch
                      <<" out of range [1,"<< myModel.size() <<"]."<< std::endl;
            ok = false;
            break;
          }

          if (mySol)
            mySol->initPatch(pch->idx);

          for (p2 = myProps.begin(); p2 != myProps.end() && ok; ++p2)
            if (p2->pcode == Property::MATERIAL && p->patch == p2->patch)
              if (!(ok = this->initMaterial(p2->pindx)))
                std::cerr <<" *** SIMbase::assembleSystem: Failed to initialize"
                          <<" material for patch "<< p2->patch << std::endl;

          if (abs(p->ldim)+1 == pch->getNoParamDim())
          {
            if (p->pcode == Property::NEUMANN_GENERIC ||
                this->initNeumann(p->pindx))
            {
              if (msgLevel > 1)
                IFEM::cout <<"\nAssembling Neumann matrix terms for boundary "
                           << p->lindx%10 <<" on P"<< p->patch << std::endl;
              if (p->patch != lp)
                ok &= this->extractPatchSolution(it->second,prevSol,p->patch-1);
              ok &= pch->integrate(*it->second,p->lindx,sysQ,time);
              lp = p->patch;
            }
            else
              ok = false;
          }
          else if (abs(p->ldim) == 1 && pch->getNoParamDim() == 3)
          {
            if (p->pcode == Property::NEUMANN_GENERIC ||
                this->initNeumann(p->pindx))
            {
              if (msgLevel > 1)
                IFEM::cout <<"\nAssembling Neumann matrix terms for edge "
                           << p->lindx%10 <<" on P"<< p->patch << std::endl;
              if (p->patch != lp)
                ok &= this->extractPatchSolution(it->second,prevSol,p->patch-1);
              ok &= pch->integrateEdge(*it->second,p->lindx,sysQ,time);
              lp = p->patch;
            }
            else
              ok = false;
          }
        }

    if (ok) ok = this->assembleDiscreteTerms(it->second,time);
    if (ok && &sysQ != myEqSys && isAssembling)
      ok = sysQ.finalize(newLHSmatrix);
  }
  if (ok && isAssembling && myEqSys)
    ok = myEqSys->finalize(newLHSmatrix);

  if (!ok)
    std::cerr <<" *** SIMbase::assembleSystem: Failure.\n"<< std::endl;

  return ok;
}


bool SIMbase::extractLoadVec (Vector& loadVec, size_t idx) const
{
  if (!myEqSys || !mySam)
    return false;

  // Expand load vector from equation ordering to DOF-ordering
  SystemVector* b = myEqSys->getVector(idx);
  if (!b || !mySam->expandSolution(*b,loadVec,0.0))
    return false;

#if SP_DEBUG > 1
  std::cout <<"\nLoad vector:"<< loadVec;
#endif
  return true;
}


double SIMbase::extractScalar (size_t idx) const
{
  return myEqSys ? myEqSys->getScalar(idx) : 0.0;
}


bool SIMbase::applyDirichlet (Vector& glbVec) const
{
  return mySam ? mySam->applyDirichlet(glbVec) : false;
}


bool SIMbase::solveSystem (Vector& solution, int printSol, double* rCond,
                           const char* compName, bool newLHS, size_t idxRHS)
{
  if (!myEqSys) return false;

  SystemMatrix* A = myEqSys->getMatrix();
  SystemVector* b = myEqSys->getVector(idxRHS);
  if (!A) std::cerr <<" *** SIMbase::solveSystem: No LHS matrix."<< std::endl;
  if (!b) std::cerr <<" *** SIMbase::solveSystem: No RHS vector."<< std::endl;
  if (!A || !b) return false;

  // Dump system matrix to file, if requested
  for (DumpData& dmp : lhsDump)
    if (dmp.doDump()) {
      IFEM::cout <<"\nDumping system matrix to file "<< dmp.fname << std::endl;
      std::ofstream os(dmp.fname.c_str(),
                       dmp.step.size() == 1 ? std::ios::out : std::ios::app);
      os << std::setprecision(17);
      double old_tol = utl::zero_print_tol;
      utl::zero_print_tol = dmp.eps;
      SystemMatrix* M = myEqSys->getMatrix(0);
      char matName[16];
      if (dmp.step.size() == 1)
        strcpy(matName,"A");
      else
        sprintf(matName,"A%d",dmp.count);
      for (int i = 0; M; M = myEqSys->getMatrix(++i), ++matName[0])
        M->dump(os,dmp.format,matName); // label matrices as A,B,C,...
      utl::zero_print_tol = old_tol;
    }

  // Dump right-hand-side vector to file, if requested
  for (DumpData& dmp : rhsDump)
    if (dmp.doDump()) {
      IFEM::cout <<"\nDumping RHS vector to file "<< dmp.fname << std::endl;
      std::ofstream os(dmp.fname.c_str(),
                       dmp.step.size() == 1 ? std::ios::out : std::ios::app);
      os << std::setprecision(17);
      double old_tol = utl::zero_print_tol;
      utl::zero_print_tol = dmp.eps;
      SystemVector* c = myEqSys->getVector(0);
      char vecName[16];
      if (dmp.step.size() == 1)
        strcpy(vecName,"b");
      else
        sprintf(vecName,"b%d",dmp.count);
      for (int i = 0; c; c = myEqSys->getVector(++i), ++vecName[0])
        c->dump(os,dmp.format,vecName); // label vectors as b,c,d,...
      utl::zero_print_tol = old_tol;
    }

  // Solve the linear system of equations
  if (msgLevel > 1)
    IFEM::cout <<"\nSolving the equation system ..."<< std::endl;

  double rcn = 1.0;
  double* rp = msgLevel > 1 ? &rcn : rCond;

  utl::profiler->start("Equation solving");
  bool status = A->solve(*b,newLHS,rp);
  utl::profiler->stop("Equation solving");

  if (msgLevel > 1)
  {
    if (rcn < 1.0)
      IFEM::cout <<"\tCondition number: "<< 1.0/rcn << std::endl;
    if (rCond) *rCond = rcn;
  }

  // Dump solution vector to file, if requested
  for (DumpData& dmp : solDump)
    if (dmp.doDump()) {
      IFEM::cout <<"Dumping solution vector to file "<< dmp.fname << std::endl;
      std::ofstream os(dmp.fname.c_str(),
                       dmp.step.size() == 1 ? std::ios::out : std::ios::app);
      os << std::setprecision(17);
      double old_tol = utl::zero_print_tol;
      utl::zero_print_tol = dmp.eps;
      char vecName[16];
      if (dmp.step.size() == 1)
        strcpy(vecName,"x");
      else
        sprintf(vecName,"x%d",dmp.count);
      b->dump(os,dmp.format,vecName);
      utl::zero_print_tol = old_tol;
    }

  // Expand solution vector from equation ordering to DOF-ordering
  if (status)
    status = mySam ? mySam->expandSolution(*b,solution) : false;

#if SP_DEBUG > 2
  if (printSol < 1000) printSol = 1000;
#endif
  if (printSol > 0 && status)
    this->printSolutionSummary(solution,printSol,compName);

  return status;
}


bool SIMbase::solveMatrixSystem (Vectors& solution, int printSol,
                                 const char* compName)
{
  solution.resize(this->getNoRHS());
  for (size_t i = 0; i < solution.size(); i++)
    if (!this->solveSystem(solution[i],printSol,nullptr,compName,i==0,i))
      return false;
    else if (solution.size() > 2)
      printSol = 0; // Print summary only for the first two solutions

  return true;
}


void SIMbase::printStep (int istep, const TimeDomain& time) const
{
  adm.cout <<"\n  step="<< istep <<"  time="<< time.t << std::endl;
}


void SIMbase::printSolutionSummary (const Vector& solution, int printSol,
                                    const char* compName,
                                    std::streamsize outPrec)
{
  // Compute and print solution norms
  const size_t nf = this->getNoFields(1);
  size_t iMax[nf > 0 ? nf : nsd];
  double dMax[nf > 0 ? nf : nsd];
  double dNorm = this->solutionNorms(solution,dMax,iMax,nf);

  int oldPrec = adm.cout.precision();
  if (outPrec > 0)
    adm.cout << std::setprecision(outPrec);

  if (compName)
    adm.cout <<"\n >>> Solution summary <<<\n\nL2-norm            : ";
  else if (nf > 1)
    adm.cout <<"  Primary solution summary: L2-norm         : ";
  else
    adm.cout <<"  Primary solution summary: L2-norm      : ";
  adm.cout << utl::trunc(dNorm);

  if (nf == 1 && utl::trunc(dMax[0]) != 0.0)
  {
    if (compName)
      adm.cout <<"\nMax "<< compName <<"   : ";
    else
      adm.cout <<"\n                            Max value    : ";
    adm.cout << dMax[0] <<" node "<< iMax[0];
  }
  else if (nf > 1)
  {
    char D = 'X';
    for (size_t d = 0; d < nf; d++, D=='Z' ? D='x' : D++)
      if (utl::trunc(dMax[d]) != 0.0)
      {
        if (compName)
          adm.cout <<"\nMax "<< D <<'-'<< compName <<" : ";
        else
          adm.cout <<"\n                            Max "<< D <<"-component : ";
        adm.cout << dMax[d] <<" node "<< iMax[d];
      }
  }
  adm.cout << std::endl;
  adm.cout << std::setprecision(oldPrec);

  // Print entire solution vector if it is small enough
  if (mySam && mySam->getNoEquations() < printSol)
  {
    adm.cout <<"\nSolution vector:";
    for (int inod = 1; inod <= mySam->getNoNodes(); inod++)
    {
      adm.cout <<"\nNode "<< inod <<":";
      std::pair<int,int> dofs = mySam->getNodeDOFs(inod);
      for (int d = dofs.first-1; d < dofs.second; d++)
        adm.cout <<" "<< utl::trunc(solution[d]);
    }
    adm.cout << std::endl;
  }
}


void SIMbase::getWorstDofs (const Vector& u, const Vector& r,
                            size_t nWorst, double eps, int iteNorm,
                            std::map<std::pair<int,int>,RealArray>& worst) const
{
  size_t i;
  RealArray data(3);
  std::multimap<double,size_t> energy;

  // Compute the energy at each DOF and insert into a map sorted on the energy
  for (i = 0; i < u.size() && i < r.size(); i++)
    if (iteNorm == 1) // L2-norm of residual
      energy.insert(std::make_pair(fabs(r[i]),i+1));
    else if (iteNorm == 2) // L2-norm of solution correction
      energy.insert(std::make_pair(fabs(u[i]),i+1));
    else // Energy norm
      energy.insert(std::make_pair(fabs(u[i]*r[i]),i+1));

  // Pick the nWorst highest energies from the back of the map
  std::multimap<double,size_t>::reverse_iterator rit = energy.rbegin();
  for (i = 0; i < nWorst && rit != energy.rend(); i++, ++rit)
    if (rit->first > eps)
    {
      data[0] = rit->first;
      data[1] = u(rit->second);
      data[2] = r(rit->second);
      worst[mySam->getNodeAndLocalDof(rit->second)] = data;
    }
}


char SIMbase::getNodeType (int inod) const
{
  return mySam ? mySam->getNodeType(inod) : ' ';
}


Vec4 SIMbase::getNodeCoord (int inod) const
{
  Vec4 Xnod;
  for (ASMbase* pch : myModel)
  {
    size_t node = pch->getNodeIndex(inod,true);
    if (node > 0)
    {
      Xnod = pch->getCoord(node);
      if (nGlPatches > 1)
        Xnod.idx = pch->idx; // Store patch index, if multi-patch model
      break;
    }
  }

  return Xnod;
}


bool SIMbase::isFixed (int inod, int dof) const
{
  for (ASMbase* pch : myModel)
  {
    size_t node = pch->getNodeIndex(inod,true);
    if (node > 0)
      return pch->isFixed(node,dof,true);
  }

  return true;
}


int SIMbase::getGlobalNode (int node) const
{
  std::map<int,int>::const_iterator it = utl::findValue(*g2l,node);
  return it != g2l->end() ? it->first : -1;
}


int SIMbase::getLocalNode (int node) const
{
  std::map<int,int>::const_iterator it = g2l->find(node);
  return it != g2l->end() ? it->second : -1;
}


SystemVector* SIMbase::getRHSvector (size_t idx, bool copy) const
{
  if (!myEqSys) return nullptr;

  SystemVector* rhs = myEqSys->getVector(idx);
  return rhs && copy ? rhs->copy() : rhs;
}


void SIMbase::addToRHSvector (size_t idx, const SystemVector& vec, double scale)
{
  if (!myEqSys) return;

  SystemVector* rhs = myEqSys->getVector(idx);
  if (!rhs || scale == 0.0) return;

  rhs->add(vec,scale);
}


void SIMbase::iterationNorms (const Vector& u, const Vector& r,
			      double& eNorm, double& rNorm, double& dNorm) const
{
  eNorm = mySam->dot(r,u,'A');
  rNorm = mySam->norm2(r,'D');
  dNorm = mySam->norm2(u,'D');
}


double SIMbase::solutionNorms (const Vector& x, double* inf,
			       size_t* ind, size_t nf, char type) const
{
  if (inf && ind && nf == 0) nf = nsd;

  for (size_t d = 0; d < nf; d++)
  {
    ind[d] = d+1;
    inf[d] = mySam->normInf(x,ind[d],type);
  }

  return mySam->normL2(x,type);
}


NormBase* SIMbase::getNormIntegrand () const
{
  return myProblem->getNormIntegrand(mySol);
}


ForceBase* SIMbase::getBoundaryForceIntegrand (const Vec3* X0) const
{
  return myProblem->getForceIntegrand(X0,mySol);
}


ForceBase* SIMbase::getNodalForceIntegrand () const
{
  return myProblem->getForceIntegrand();
}


bool SIMbase::solutionNorms (const TimeDomain& time,
			     const Vectors& psol, const Vectors& ssol,
			     Vectors& gNorm, Matrix* eNorm, const char* name)
{
  if (!mySam) return true; // Silently ignore when uninitialized system

  PROFILE1("Norm integration");

  NormBase* norm = myProblem->getNormIntegrand(mySol);
  if (!norm)
  {
#ifdef SP_DEBUG
    std::cerr <<"  ** SIMbase::solutionNorms: No integrand."<< std::endl;
#endif
    return true; // Silently ignore when no norm integrand is provided
  }

  if (msgLevel > 1 && name)
    IFEM::cout <<"\nIntegrating solution norms ("<< name <<") ..."<< std::endl;

  myProblem->initIntegration(time,psol.front());
  norm->initProjection(ssol.size());
  norm->initIntegration(nIntGP,nBouGP);

  // Number of recovered solution components
  size_t i, nCmp = 0;
  for (i = 0; i < ssol.size() && nCmp == 0; i++)
    nCmp = ssol[i].size() / this->getNoNodes(1);

#ifdef USE_OPENMP
  // When assembling in parallel, we must always do the norm summation
  // at the end in a serial loop, to avoid that the threads try to update
  // the same memory address simultaneously.
  Matrix dummy;
  if (!eNorm) eNorm = &dummy;
#endif

  // Initialize norm integral classes
  gNorm.resize(norm->getNoFields(0));
  size_t nNorms = 0;
  auto proj_idx = opt.project.begin();
  // count norms if they are:
  // 1) associated with the primary solution
  // 2) associated with a projected secondary solution
  // 3) associated with a residual secondary solution
  // 4) associated with a secondary solution provided through external means
  // 5) associated with an additional group defined in the integrand
  for (i = 0; i < gNorm.size(); i++)
    if (i == 0 || i > ssol.size() || !ssol[i-1].empty() ||
        proj_idx->first == SIMoptions::NONE || norm->hasExternalProjections())
    {
      size_t nNrm = norm->getNoFields(1+i);
      gNorm[i].resize(nNrm,true);
      nNorms += nNrm;
      if (i > 0 && proj_idx != opt.project.end())
        ++proj_idx;
    }

  GlbNorm globalNorm(gNorm,norm->getFinalOperation());
  LintegralVec elementNorms;
  if (eNorm)
  {
    globalNorm.delayAssembly();
    eNorm->resize(nNorms,mySam->getNoElms(),true);
    elementNorms.reserve(eNorm->cols());
    for (i = 0; i < eNorm->cols(); i++)
      elementNorms.push_back(new ElmNorm(eNorm->ptr(i),nNorms));
    norm->setLocalIntegrals(&elementNorms);
  }

  // Loop over the different material regions, integrating solution norm terms
  // for the patch domain associated with each material
  bool ok = true;
  size_t k, lp = 0;
  ASMbase* pch = nullptr;
  PropertyVec::const_iterator p;
  size_t projOfs = 0;
  for (p = myProps.begin(); p != myProps.end() && ok; ++p)
    if (p->pcode == Property::MATERIAL)
      if (!(pch = this->getPatch(p->patch)))
	ok = false;
      else if (this->initMaterial(p->pindx))
      {
	lp = p->patch;
	ok = this->extractPatchSolution(psol,lp-1);
	size_t nval = pch->getNoProjectionNodes()*myProblem->getNoFields(2);
	for (k = 0; k < ssol.size(); k++)
          if (ssol[k].empty())
            norm->getProjection(k).clear();
          else if (this->fieldProjections()) {
            Vector c(nval);
            std::copy(ssol[k].begin()+projOfs,
                      ssol[k].begin()+projOfs+nval, c.begin());
            Fields* f = pch->getProjectedFields(c, myProblem->getNoFields(2));
            norm->setProjectedFields(f, k);
            projOfs += nval;
          } else
            this->extractPatchSolution(ssol[k],norm->getProjection(k),pch,nCmp,1);

        if (mySol)
          mySol->initPatch(pch->idx);

	ok &= pch->integrate(*norm,globalNorm,time);
        if (norm->getIntegrandType() & IntegrandBase::INTERFACE_TERMS) {
          ASM::InterfaceChecker* iChk = this->getInterfaceChecker(pch->idx);
          if (iChk) {
            ok &= !pch->integrate(*norm,globalNorm,time,*iChk);
            delete iChk;
          }
        }
      }
      else
	ok = false;

  if (lp == 0)
    // All patches are referring to the same material, and we assume it has
    // been initialized during input processing (thus no initMaterial call here)
    for (i = 0; i < myModel.size() && ok; i++)
    {
      ok = this->extractPatchSolution(psol,i);
      size_t nval = myModel[i]->getNoProjectionNodes()*myProblem->getNoFields(2);
      for (k = 0; k < ssol.size(); k++)
        if (ssol[k].empty())
          norm->getProjection(k).clear();
        else if (this->fieldProjections()) {
          Vector c(nval);
          std::copy(ssol[k].begin()+projOfs,
                    ssol[k].begin()+projOfs+nval, c.begin());
          Fields* f = myModel[i]->getProjectedFields(c, myProblem->getNoFields(2));
          norm->setProjectedFields(f, k);
          projOfs += nval;
        } else
          this->extractPatchSolution(ssol[k],norm->getProjection(k),myModel[i],nCmp,1);

      if (mySol)
        mySol->initPatch(myModel[i]->idx);

      ok &= myModel[i]->integrate(*norm,globalNorm,time);
      if (norm->getIntegrandType() & IntegrandBase::INTERFACE_TERMS) {
        ASM::InterfaceChecker* iChk = this->getInterfaceChecker(myModel[i]->idx);
        if (iChk) {
          ok &= myModel[i]->integrate(*norm,globalNorm,time,*iChk);
          delete iChk;
        }
      }
      lp = i+1;
    }

  if (norm->hasBoundaryTerms())
    for (p = myProps.begin(); p != myProps.end() && ok; ++p)
      if (p->pcode == Property::NEUMANN)
	if (!(pch = this->getPatch(p->patch)))
	  ok = false;

	else if (abs(p->ldim)+1 == pch->getNoParamDim())
	  if (this->initNeumann(p->pindx))
	  {
	    if (p->patch != lp)
	      ok = this->extractPatchSolution(psol,p->patch-1);
            if (mySol)
              mySol->initPatch(pch->idx);
	    ok &= pch->integrate(*norm,p->lindx,globalNorm,time);
	    lp = p->patch;
	  }
	  else
	    ok = false;

	else if (abs(p->ldim)+2 == pch->getNoParamDim())
	  if (this->initNeumann(p->pindx))
	  {
	    if (p->patch != lp)
	      ok = this->extractPatchSolution(psol,p->patch-1);
            if (mySol)
              mySol->initPatch(pch->idx);
	    ok &= pch->integrateEdge(*norm,p->lindx,globalNorm,time);
	    lp = p->patch;
	  }
	  else
	    ok = false;

  if (!ok) std::cerr <<" *** SIMbase::solutionNorms: Failure.\n"<< std::endl;

  // Clean up the dynamically allocated norm objects. This will also perform
  // the actual global norm assembly, in case the element norms are stored,
  // and always when multi-threading is used.
  for (k = 0; k < elementNorms.size(); k++)
  {
    globalNorm.assemble(elementNorms[k]);
    delete elementNorms[k];
  }

  // Add problem-dependent external norm contributions
  norm->addBoundaryTerms(gNorm,this->externalEnergy(psol,time));

  delete norm;

  if (!adm.dd.isPartitioned())
    for (k = 0; k < gNorm.size(); k++)
      adm.allReduceAsSum(gNorm[k]);

  if (ok)
    ok = this->postProcessNorms(gNorm, eNorm);

  return ok;
}


double SIMbase::externalEnergy (const Vectors& psol, const TimeDomain&) const
{
  const Vector* reactionForces = this->getReactionForces();
  if (!reactionForces || !mySam || psol.empty()) return 0.0;

  // Add norm contributions due to inhomogeneous Dirichlet boundary conditions.
  // That is, the path integral of the total solution vector times the
  // reaction forces at the prescribed DOFs.
  if (psol.size() == 1)
    return mySam->normReact(psol.front(),*reactionForces);

  if (prevForces.empty())
    extEnergy += mySam->normReact(psol[0]-psol[1],*reactionForces);
  else
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
  const Vector* reactionForces = this->getReactionForces();
  if (!reactionForces || !mySam) return false;

  RF.resize(1+nsd);
  RF.front() = 2.0*mySam->normReact(psol,*reactionForces);
  for (unsigned char dir = 1; dir <= nsd; dir++)
    RF[dir] = mySam->getReaction(dir,*reactionForces);

  return true;
}


bool SIMbase::getCurrentReactions (RealArray& RF, int pcode) const
{
  const Vector* reactionForces = this->getReactionForces();
  if (!reactionForces || !mySam) return false;

  IntVec glbNodes;
  this->getBoundaryNodes(pcode,glbNodes);

  RF.resize(nsd);
  for (unsigned char dir = 1; dir <= nsd; dir++)
    RF[dir-1] = mySam->getReaction(dir,*reactionForces,&glbNodes);

  return true;
}


const Vector* SIMbase::getReactionForces () const
{
  return myEqSys ? myEqSys->getReactions() : nullptr;
}


bool SIMbase::systemModes (std::vector<Mode>& solution,
			   int nev, int ncv, int iop, double shift,
			   size_t iA, size_t iB)
{
  if (nev < 1 || ncv <= nev || !myEqSys || !mySam) return false;

  PROFILE1("Eigenvalue analysis");

  Vector eigVal;
  Matrix eigVec;
  if (nev > mySam->getNoEquations()) nev = mySam->getNoEquations();
  if (ncv > mySam->getNoEquations()) ncv = mySam->getNoEquations();

  // Solve the eigenvalue problem
  IFEM::cout <<"\nSolving the eigenvalue problem ..."<< std::endl;
  SystemMatrix* A = myEqSys->getMatrix(iA);
  SystemMatrix* B = myEqSys->getMatrix(iB);
#ifdef HAS_SLEPC
  // To interface SLEPC another interface is used
  bool ok = eig::solve(A,B,eigVal,eigVec,nev);
#else
  bool ok = eig::solve(A,B,eigVal,eigVec,nev,ncv,iop,shift);
#endif

  // Expand eigenvectors to DOF-ordering and print out eigenvalues
  bool freq = iop == 3 || iop == 4 || iop == 6;
  IFEM::cout <<"\n >>> Computed Eigenvalues <<<\n     Mode\t"
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

    IFEM::cout <<"\n     "<< i <<"\t\t"<< utl::trunc(solution[i-1].eigVal);
  }
  IFEM::cout << std::endl;
  return ok;
}


bool SIMbase::project (Matrix& ssol, const Vector& psol,
		       SIMoptions::ProjectionMethod pMethod,
		       const TimeDomain& time) const
{
  PROFILE1("Solution projection");

  if (msgLevel > 1)
    IFEM::cout <<"\nProjecting secondary solution ..."<< std::endl;

  ssol.clear();

  size_t ngNodes = this->fieldProjections() ? 0 : this->getNoNodes(1);
  Vector count(myModel.size() > 1 ? ngNodes : 0);
  Matrix values;

  if (this->fieldProjections())
    // No nodal averaging - full (potentially discontinuous) representation
    for (ASMbase* pch : myModel)
      if (!pch->empty())
        ngNodes += pch->getNoProjectionNodes();

  if (pMethod == SIMoptions::DGL2 || pMethod == SIMoptions::CGL2)
    const_cast<SIMbase*>(this)->setQuadratureRule(opt.nGauss[1]);
  else if (pMethod == SIMoptions::CGL2_INT)
  {
    // Reinitialize the integration point buffers within the integrands (if any)
    // Note: We here use the same integration scheme as for the main problem
    const_cast<SIMbase*>(this)->setQuadratureRule(opt.nGauss[0]);
    myProblem->initIntegration(time,psol);
  }

  size_t i, ofs = 0;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    // Extract the primary solution control point values for this patch
    myProblem->initResultPoints(time.t);
    if (!this->extractPatchSolution(myProblem,Vectors(1,psol),i))
      return false;

    // Initialize material properties for this patch in case of multiple regions
    const_cast<SIMbase*>(this)->setPatchMaterial(i+1);

    // Project the secondary solution and retrieve control point values
    bool ok = false;
    switch (pMethod) {
    case SIMoptions::GLOBAL:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tGreville point projection"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem);
      break;

    case SIMoptions::DGL2:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tDiscrete global L2-projection"<< std::endl;
      ok = myModel[i]->globalL2projection(values,*myProblem);
      break;

    case SIMoptions::CGL2:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tContinuous global L2-projection"<< std::endl;
      ok = myModel[i]->globalL2projection(values,*myProblem,true);
      break;

    case SIMoptions::CGL2_INT:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tContinuous global L2-projection"<< std::endl;
      ok = myModel[i]->L2projection(values,myProblem,time);
      break;

    case SIMoptions::SCR:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tSuperconvergent recovery"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem,nullptr,'S');
      break;

    case SIMoptions::VDSA:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tVariation diminishing projection"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem,nullptr,'A');
      break;

    case SIMoptions::QUASI:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tQuasi interpolation"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem,nullptr,'L');
      break;

    case SIMoptions::LEASTSQ:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tLeast squares projection"<< std::endl;
      ok = myModel[i]->evalSolution(values,*myProblem,nullptr,'W');
      break;

    default:
      std::cerr <<" *** SIMbase::project: Projection method "<< pMethod
                <<" not implemented."<< std::endl;
    }

    if (!ok)
    {
      std::cerr <<" *** SIMbase::project: Failure when projecting patch "
                << myModel[i]->idx+1 <<"."<< std::endl;
      return false;
    }

    // If true, we cannot assume we have a multi-patch numbering for
    // patch projections, so simply append each vector successively.
    if (this->fieldProjections()) {
      if (ssol.empty())
        ssol.resize(myProblem->getNoFields(2),ngNodes);

      ssol.fillBlock(values, 1, ofs+1);
      ofs += myModel[i]->getNoProjectionNodes()*myProblem->getNoFields(2);
    } else {
      size_t nComps = values.rows();
      size_t nNodes = values.cols();
      if (ssol.empty())
        ssol.resize(nComps,ngNodes);

      // Nodal averaging for nodes referred to by two or more patches
      // (these are typically the interface nodes)
      for (size_t n = 1; n <= nNodes; n++)
        if (count.empty())
          ssol.fillColumn(myModel[i]->getNodeID(n),values.getColumn(n));
        else
        {
          int inod = myModel[i]->getNodeID(n);
          for (size_t j = 1; j <= nComps; j++)
            ssol(j,inod) += values(j,n);
          count(inod) ++;
        }
    }
  }

  // Divide through by count(n) to get the nodal average at the interface nodes
  for (size_t n = 1; n <= count.size(); n++)
    if (count(n) > 1.0)
      for (size_t j = 1; j <= ssol.rows(); j++)
        ssol(j,n) /= count(n);

  return true;
}


bool SIMbase::project (Vector& values, const FunctionBase* f,
                       int basis, int iField, int nFields,
                       SIMoptions::ProjectionMethod pMethod, double time) const
{
  bool ok = true;
  for (size_t j = 0; j < myModel.size() && ok; j++)
  {
    if (myModel[j]->empty()) continue; // skip empty patches

    Vector loc_values;
    Matrix f_values;
    switch (pMethod) {
    case SIMoptions::GLOBAL:
      // Greville point projection
      ok = myModel[j]->evaluate(f,loc_values,basis,time);
      break;

    case SIMoptions::CGL2:
    case SIMoptions::CGL2_INT:
      // Continuous global L2-projection
      if (myModel[j]->separateProjectionBasis())
      {
        // Not implemented yet, silently ignore unless debug build
#ifdef SP_DEBUG
        std::cerr <<"  ** L2 projection of explicit functions onto a"
                  <<" separate basis is not available."<< std::endl;
#endif
        return false;
      }
      ok = myModel[j]->L2projection(f_values,const_cast<FunctionBase*>(f),time);
      loc_values = f_values;
      break;

    default:
      std::cerr <<" *** SIMbase::project: Projection method "<< pMethod
                <<" not implemented for functions."<< std::endl;
      return false;
    }

    if (nFields <= (int)f->dim())
      ok &= this->injectPatchSolution(values,loc_values,
                                      myModel[j],f->dim(),basis);
    else if (f->dim() > 1)
    {
      std::cerr <<" *** SIMbase::project: Cannot interleave non-scalar function"
                <<" into a multi-dimensional field."<< std::endl;
      return false;
    }
    else if (iField < nFields)
    {
      // Interleave
      size_t i, k = iField;
      Vector loc_vector(loc_values.size()*nFields);
      myModel[j]->extractNodeVec(values,loc_vector,0,basis);
      for (i = 0; i < loc_values.size(); i++, k += nFields)
        loc_vector[k] = loc_values[i];
      ok &= myModel[j]->injectNodeVec(loc_vector,values,0,basis);
    }
    else
    {
      std::cerr <<" *** SIMbase::project: Field component "<< iField+1
                <<" is out of range [1,"<< nFields <<"]."<< std::endl;
      return false;
    }
  }

  return ok;
}


bool SIMbase::extractPatchSolution (IntegrandBase* problem,
                                    const Vectors& sol, size_t pindx) const
{
  ASMbase* pch = this->getPatch(pindx+1);
  if (!pch || !mySam) return false;

  problem->initNodeMap(pch->getGlobalNodeNums());
  for (size_t i = 0; i < sol.size() && i < problem->getNoSolutions(); i++)
    if (!sol[i].empty())
      pch->extractNodalVec(sol[i],problem->getSolution(i),mySam->getMADOF());

  return this->extractPatchDependencies(problem,myModel,pindx);
}


bool SIMbase::project (Vector& ssol, const Vector& psol,
		       SIMoptions::ProjectionMethod pMethod, size_t iComp) const
{
  Matrix stmp;
  if (!this->project(stmp,psol,pMethod))
    return false;

  if (iComp > 0)
    ssol = stmp.getRow(iComp);
  else
    ssol = stmp;

  return true;
}


bool SIMbase::projectAnaSol (Vector& ssol,
                             SIMoptions::ProjectionMethod pMethod) const
{
  if (!mySol) return false;

  FunctionBase* f = mySol->getScalarSecSol();
  if (f)
  {
    ssol.resize(f->dim()*this->getNoNodes());
    return this->project(ssol,f,0,0,0,pMethod);
  }

  f = mySol->getVectorSecSol();
  if (f)
  {
    ssol.resize(f->dim()*this->getNoNodes());
    return this->project(ssol,f,0,0,0,pMethod);
  }

  f = mySol->getStressSol();
  if (f)
  {
    ssol.resize(f->dim()*this->getNoNodes());
    return this->project(ssol,f,0,0,0,pMethod);
  }

  return false;
}


size_t SIMbase::extractPatchSolution (const Vector& sol, Vector& vec,
                                      const ASMbase* pch, unsigned char nndof,
                                      unsigned char basis) const
{
  if (!pch || sol.empty()) return 0;

  if (basis && nndof != pch->getNoFields(basis) && pch->getNoFields(2) > 0)
  {
    // Need an additional MADOF
    const IntVec& madof = this->getMADOF(basis,nndof);
    pch->extractNodalVec(sol,vec,madof.data(),madof.size());
  }
  else
    pch->extractNodeVec(sol,vec,nndof,basis);

  return vec.size();
}


bool SIMbase::injectPatchSolution (Vector& sol, const Vector& vec,
                                   const ASMbase* pch, unsigned char nndof,
                                   unsigned char basis) const
{
  if (!pch) return false;

  if (basis && nndof != pch->getNoFields(basis) && pch->getNoFields(2) > 0)
    // Need an additional MADOF
    return pch->injectNodalVec(vec,sol,this->getMADOF(basis,nndof),basis);
  else
    return pch->injectNodeVec(vec,sol,nndof,basis);
}


bool SIMbase::evalSecondarySolution (Matrix& field, int pindx) const
{
  ASMbase* pch = pindx >= 0 ? this->getPatch(pindx+1) : nullptr;
  if (!pch) return false;

  const_cast<SIMbase*>(this)->setPatchMaterial(pindx+1);
  return pch->evalSolution(field,*myProblem);
}


bool SIMbase::fieldProjections () const
{
  for (const ASMbase* pch : myModel)
    if (pch->separateProjectionBasis())
      return true;

  return false;
}


bool SIMbase::extractPatchElmRes (const Matrix& glbRes, Matrix& elRes,
				  int pindx) const
{
  ASMbase* pch = pindx >= 0 ? this->getPatch(pindx+1) : nullptr;
  if (!pch || glbRes.empty()) return false;

  pch->extractElmRes(glbRes,elRes);
  return true;
}


bool SIMbase::setPatchMaterial (size_t patch)
{
  for (const Property& p : myProps)
    if (p.pcode == Property::MATERIAL && p.patch == patch)
      return this->initMaterial(p.pindx);

  return false;
}


bool SIMbase::addMADOF (unsigned char basis, unsigned char nndof, bool other)
{
  int key = basis << 16 + nndof;
  if (mixedMADOFs.find(key) != mixedMADOFs.end())
    return false; // This MADOF already calculated

  IntVec& madof = mixedMADOFs[key];
  madof.resize(this->getNoNodes()+1,0);
  if (madof.size() < 2) return false;

  int ierr = 0;
  char nType = basis <= 1 ? 'D' : 'P' + basis-2;
  for (const ASMbase* pch : myModel)
    for (size_t inod = 1; inod <= pch->getNoNodes(); inod++)
    {
      int node = pch->getNodeID(inod);
      if (node > 0)
      {
        int nndofs = 0;
        if (pch->getNodeType(inod) == nType)
          nndofs = nndof;
        else if (other)
          nndofs = pch->getNodalDOFs(inod);
        if (nndofs > 0 && madof[node] == 0)
          madof[node] = nndofs;
        else if (madof[node] != nndofs)
          ierr++;
      }
    }

  madof[0] = 1;
  for (size_t n = 1; n < madof.size(); n++)
    madof[n] += madof[n-1];

  if (ierr == 0)
    return true;

  std::cerr <<" *** SIMbase::addMADOF: Detected "<< ierr <<" nodes with"
            <<" conflicting number of DOFs in adjacent patches."<< std::endl;
  return false;
}


const IntVec& SIMbase::getMADOF (unsigned char basis,
                                 unsigned char nndof) const
{
  int key = basis << 16 + nndof;
  auto it = mixedMADOFs.find(key);
  if (it != mixedMADOFs.end())
    return it->second;

  static IntVec empty;
  return empty;
}


void SIMbase::registerDependency (const std::string& name, SIMdependency* sim,
                                  short int nvc, unsigned char basis)
{
  this->registerDependency(sim,name,nvc,myModel,
                           this->getMADOF(basis,nvc).data());
}
