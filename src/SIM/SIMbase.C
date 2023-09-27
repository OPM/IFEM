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
#include "ASMmxBase.h"
#include "ASMs2DC1.h"
#ifdef HAS_PETSC
#include "SAMpatchPETSc.h"
#else
#include "SAMpatch.h"
#endif
#include "IntegrandBase.h"
#include "AlgEqSystem.h"
#include "ReactionsOnly.h"
#include "SystemMatrix.h"
#include "LinSolParams.h"
#include "EigSolver.h"
#include "GlbNorm.h"
#include "GlbL2projector.h"
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
  nsd = 3;
  myProblem = itg;
  mySol = nullptr;
  myEqSys = nullptr;
  mySam = nullptr;
  mySolParams = nullptr;
  myGl2Params = nullptr;
  dualField = nullptr;
  isRefined = lagMTOK = false;
  nGlbNodes = nGlPatches = 0;
  nIntGP = nBouGP = 0;
  nDofS = 0;
  mdFlag = 0;
  extEnergy = 0.0;

  MPCLess::compareSlaveDofOnly = true; // to avoid multiple slave definitions
}


SIMbase::~SIMbase ()
{
#ifdef SP_DEBUG
  IFEM::cout <<"\nEntering SIMbase destructor"<< std::endl;
#endif

  for (IntegrandMap::value_type& itg : myInts)
    if (itg.second != myProblem)
      delete itg.second;

  delete myProblem;
  delete mySol;
  if (mdFlag <= 1)
  {
    delete myEqSys;
    delete mySam;
  }
  delete mySolParams;
  delete myGl2Params;

  for (ASMbase* patch : myModel)
    delete patch;

  myModel.clear();
  this->SIMbase::clearProperties();

#ifdef SP_DEBUG
  IFEM::cout <<"Leaving SIMbase destructor"<< std::endl;
#endif
}


/*!
  Adds profiling to the virtual read() method.
*/

bool SIMbase::readModel (const char* fileName)
{
  PROFILE1("Model input");
  return this->read(fileName);
}


void SIMbase::clearProperties ()
{
  delete dualField;
  dualField = nullptr;

  for (ASMbase* patch : myModel)
    patch->clear(true); // retain the geometry only

  for (auto& i2 : myScalars)
    delete i2.second;
  for (auto& i3 : myVectors)
    delete i3.second;
  for (auto& i4 : myTracs)
    delete i4.second;
  for (FunctionBase* f : extrFunc)
    delete f;

  nGlbNodes = 0;
  myPatches.clear();
  myGlb2Loc.clear();
  myScalars.clear();
  myVectors.clear();
  myTracs.clear();
  myProps.clear();
  myInts.clear();
  extrFunc.clear();
  extraMADOFs.clear();
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

  for (const std::pair<const int,Ipairs>& node : nodeInfo)
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

  if (mySol)
    mySol->setupSecondarySolutions();

  // Perform some sub-class specific pre-preprocessing, if any
  this->preprocessA();

  // Create the classical FE data structures
  if (!this->createFEMmodel('Y'))
    return false;

  // Erase all patches that should be ignored in the analysis
  for (int idx : ignored)
  {
    ASMbase* pch = this->getPatch(idx);
    if (pch) pch->clear();
  }

  // If material properties are specified for at least one patch, assign the
  // property code 999999 to all patches with no material property code yet
  PatchVec patches;
  for (const Property& p : myProps)
    if (p.pcode == Property::MATERIAL)
    {
      ASMbase* pch = this->getPatch(p.patch);
      if (pch && !pch->empty())
        patches.push_back(pch);
    }

  if (!patches.empty())
    for (size_t i = 0; i < myModel.size(); i++)
      if (std::find(patches.begin(),patches.end(),myModel[i]) == patches.end())
        myProps.push_back(Property(Property::MATERIAL,999999,i+1,
                                   myModel[i]->getNoParamDim()));

  if (fixDup)
  {
    // Check for duplicated nodes (missing topology)
    int nDupl = 0;
    std::map<Vec3,int> globalNodes;
    for (ASMbase* pch : myModel)
      if (!pch->empty())
      {
	IFEM::cout <<"   * Checking Patch "<< pch->idx+1 << std::endl;
	for (size_t node = 1; node <= pch->getNoNodes(); node++)
	{
	  Vec3 X(pch->getCoord(node));
	  std::map<Vec3,int>::const_iterator xit = globalNodes.find(X);
	  if (xit == globalNodes.end())
	    globalNodes.insert(std::make_pair(X,pch->getNodeID(node)));
	  else if (pch->mergeNodes(node,xit->second))
	    nDupl++;
	}
      }
    if (nDupl > 0)
      IFEM::cout <<"   * "<< nDupl <<" duplicated nodes merged."<< std::endl;
  }

#if SP_DEBUG > 2
  printNodalConnectivity(myModel,std::cout);
#endif

  // Renumber the nodes to account for resolved patch topology
  if (!nGlbNodes) nGlbNodes = this->renumberNodes();
  if (nGlbNodes < 0) return false;

  // Perform specialized preprocessing before the assembly initialization.
  // This typically involves the system-level Lagrange multipliers, etc.
  ASMbase::resetNumbering(nGlbNodes); // to account for possibly added nodes
  if (!this->preprocessBeforeAsmInit(nGlbNodes))
    return false;

  IFEM::cout <<"\nResolving Dirichlet boundary conditions"<< std::endl;

  // Process the Dirichlet boundary conditions in the order of increasing
  // dimension, such that vertex definitions override definitions on edges,
  // and edge definitions override definitions on faces
  size_t nprop = 0;
  int code, dofs, ierr = 0, iwar = 0;
  for (unsigned char dim = 0; nprop < myProps.size(); dim++)
    for (const Property& p : myProps)
      if (abs(p.ldim)%4 == dim)
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
          if (this->addConstraint(p.patch,p.lindx,p.ldim,dofs,code,nGlbNodes,
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

  // Compute coupling coefficients for C1-continuity constraints, if any
  for (ASMbase* pch : myModel)
    if (!pch->initConstraints())
      return false;

  // Set the inhomogeneous dirichlet conditions, in stationary case
  bool timeDependent = this->hasTimeDependentDirichlet();
  if (!timeDependent && !this->initDirichlet())
    return false;

  // Resolve possibly chaining of the MPC equations
  if (!allMPCs.empty())
    ASMbase::resolveMPCchains(allMPCs,myModel,timeDependent);

  // Set initial values for the inhomogeneous dirichlet conditions, if any
  if (timeDependent && !this->initDirichlet())
    return false;

  // Generate element groups for multi-threading
  bool silence = msgLevel < 1 || (msgLevel < 3 && nGlPatches > 1);
  for (ASMbase* pch : myModel)
    if (!pch->empty())
      pch->generateThreadGroups(*myProblem,silence,lagMTOK);

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
  delete mySam;
#ifdef HAS_PETSC
  if (opt.solver == LinAlg::PETSC)
    mySam = new SAMpatchPETSc(*g2l,adm);
  else
    mySam = new SAMpatch();
#else
  mySam = new SAMpatch();
#endif

  if (!static_cast<SAMpatch*>(mySam)->init(myModel,nGlbNodes,dofTypes))
  {
#ifdef SP_DEBUG
    for (ASMbase* pch : myModel)
      pch->printElements(std::cout);
#endif
    return false;
  }

  if (mdFlag > 0)
    nDofS = mySam->getNoDOFs();
  else if (!adm.dd.setup(adm,*this))
  {
    std::cerr <<"\n *** SIMbase::preprocess(): Failed to establish "
              <<" domain decomposition data."<< std::endl;
    return false;
  }

  // Now perform the sub-class specific final preprocessing, if any
  return this->preprocessB() && ierr == 0;
}


bool SIMbase::merge (SIMbase* that, const std::map<int,int>* old2new, int poff)
{
  if (this == that)
    return true;
  else if (!mySam || this->mdFlag > 1 || that->mdFlag < 2)
  {
    std::cerr <<" *** SIMbase::merge: Logic error, this->mdFlag = "
              << (short int)this->mdFlag <<" that->mdFlag = "
              << (short int)that->mdFlag << std::endl;
    return false;
  }

  for (ASMbase* pch : that->myModel)
    pch->idx += poff;

  that->shiftGlobalNums(mySam->getNoNodes(),mySam->getNoElms());
  if (!mySam->merge(that->mySam,old2new))
    return false;

  delete that->mySam;
  that->mySam = this->mySam;

  // Set up domain decomposition data when we have merged all models
  if (that->mdFlag == 3 && !adm.dd.setup(adm,*this))
  {
    std::cerr <<"\n *** SIMbase::merge(): Failed to establish "
              <<" domain decomposition data."<< std::endl;
    return false;
  }

  return true;
}


/*!
  This method renumbers the nodes to account for overlapping nodes
  for multi-patch models, erased patches, etc. In parallel simulations,
  the resulting global-to-local node number mapping \ref myGlb2Loc will map
  the global node numbers to local node numbers on the current processor.
  In serial simulations, this mapping will be unity unless the original global
  node number sequence had "holes" due to duplicated nodes or erased patches.
*/

int SIMbase::renumberNodes ()
{
  int ngnod = 0;
  int renum = 0;
  if (preserveNOrder)
  {
    renum = ASMbase::renumberNodes(myModel,myGlb2Loc);
    ngnod = g2l->size();
  }
  else for (ASMbase* pch : myModel)
    renum += pch->renumberNodes(myGlb2Loc,ngnod);

  if (renum > 0)
    IFEM::cout <<"\nRenumbered "<< renum <<" nodes."<< std::endl;

  // Apply the old-to-new node number mapping to all node tables in the model
  return this->renumberNodes(*g2l) ? ngnod : -ngnod;
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


bool SIMbase::initSystem (LinAlg::MatrixType mType,
                          size_t nMats, size_t nVec, size_t nScl, bool withRF)
{
  if (!mySam) return true; // Silently ignore when no algebraic system

#if SP_DEBUG > 2
  mySam->print(std::cout);
  for (ASMbase* pch : myModel)
    if (!pch->empty())
    {
      pch->printNodes(std::cout);
      pch->printElements(std::cout);
    }
  printNodalConnectivity(myModel,std::cout);
#endif

  delete myEqSys;
  myEqSys = new AlgEqSystem(*mySam,&adm);

  // Workaround SuperLU bug for tiny systems
  if (mType == LinAlg::SPARSE && this->getNoElms(true) < 3)
  {
    std::cerr <<"  ** System too small for SuperLU, falling back to Dense."
              << std::endl;
    mType = LinAlg::DENSE;
  }

  return myEqSys->init(mType, mySolParams, nMats, nVec, nScl,
                       withRF, opt.num_threads_SLU);
}


bool SIMbase::initSystem (const SIMbase* that)
{
  if (this == that)
    return false;
  else if (myEqSys || this->mdFlag < 2 || that->mdFlag > 1)
  {
    std::cerr <<" *** SIMbase::initSystem: Logic error, this->mdFlag = "
              << (short int)this->mdFlag <<" that->mdFlag = "
              << (short int)that->mdFlag << std::endl;
    return false;
  }

  this->myEqSys = that->myEqSys;
  return true;
}


void SIMbase::initLHSbuffers ()
{
  if (myProblem)
    myProblem->initLHSbuffers(this->getNoElms());
}


bool SIMbase::setAssociatedRHS (size_t iMat, size_t iVec)
{
  return myEqSys ? myEqSys->setAssociatedVector(iMat,iVec) : false;
}


bool SIMbase::setMode (int mode, bool needIntegr, bool resetSol)
{
  if (myInts.empty() && (myProblem || needIntegr))
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
  for (ASMbase* pch : myModel)
    if (!pch->empty())
    {
      pch->setGauss(ng);
      pch->setCachePolicy(opt.policy);
      // Count the interior integration points
      pch->getNoIntPoints(nIntGP,nInterfaceGP);
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
        myModel[p->patch-1]->getNoBouPoints(nBouGP,abs(p->ldim),p->lindx);
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
  std::cout <<"\nProperty mapping:\n";
  for (const Property& p : myProps)
  {
    switch (p.pcode) {
    case Property::MATERIAL:
      std::cout <<"MATERIAL        : ";
      break;
    case Property::BODYLOAD:
      std::cout <<"BODYLOAD        : ";
      break;
    case Property::NEUMANN:
      std::cout <<"NEUMANN         : ";
      break;
    case Property::NEUMANN_ANASOL:
      std::cout <<"NEUMANN_ANASOL  : ";
      break;
    case Property::NEUMANN_GENERIC:
      std::cout <<"NEUMANN_GENERIC : ";
      break;
    case Property::ROBIN:
      std::cout <<"ROBIN           : ";
      break;
    case Property::RIGID:
      std::cout <<"RIGID           : ";
      break;
    case Property::DIRICHLET:
      std::cout <<"DIRICHLET       : ";
      break;
    case Property::DIRICHLET_INHOM:
      std::cout <<"DIRICHLET_INHOM : ";
      break;
    case Property::DIRICHLET_ANASOL:
      std::cout <<"DIRICHLET_ANASOL: ";
      break;
    case Property::OTHER:
      std::cout <<"OTHER           : ";
      break;
    default:
      continue;
    }
    std::cout <<"code="<< p.pindx <<"\tP"<< p.patch
              <<"\ttop. item "<< (int)p.lindx
              <<" "<< (int)p.ldim <<"D"<< std::endl;
  }
#endif
}


size_t SIMbase::getNoFields (int basis) const
{
  return myModel.empty() ? 0 : myModel.front()->getNoFields(basis);
}


size_t SIMbase::getNoDOFs (bool subSim) const
{
  return subSim ? nDofS : (mySam ? mySam->getNoDOFs() : 0);
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
  else if (basis < 1 && nGlbNodes > 0)
    return nGlbNodes;

  std::cerr <<" *** SIMbase::getNoNodes: Number of nodes in a multi-patch model"
            <<" is not known at this point."<< std::endl;
  return 0;
}


size_t SIMbase::getNoElms (bool includeXelms) const
{
  if (mySam && includeXelms)
    return mySam->getNoElms();

  size_t noElms = 0;
  for (ASMbase* pch : myModel)
    noElms += pch->getNoElms(false,includeXelms);

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


size_t SIMbase::getNoConstraints () const
{
  return mySam ? mySam->getNoConstraints() : 0;
}


size_t SIMbase::getNoRHS () const
{
  return myEqSys ? myEqSys->getNoRHS() : (this->haveDualSol() ? 2 : 1);
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
  for (ASMbase* pch : myModel)
    if (pch->hasTimeDependentDirichlet(myScalars,myVectors))
      return true;

  return false;
}


bool SIMbase::initDirichlet (double time)
{
  Vector dummy;
  return this->updateDirichlet(time,&dummy);
}


bool SIMbase::updateDirichlet (double time, const Vector* prevSol)
{
  if (prevSol)
    for (ASMbase* pch : myModel)
      if (!pch->updateDirichlet(myScalars,myVectors,time))
        return false;

  SAMpatch* pSam = dynamic_cast<SAMpatch*>(mySam);
  return pSam ? pSam->updateConstraintEqs(prevSol) : true;
}


bool SIMbase::updateGrid (const RealArray& displ)
{
  if (displ.empty()) return true; // No displacements (yet), totally fine

  for (ASMbase* pch : myModel)
  {
    Vector locdisp;
    if (this->mixedProblem())
      this->extractPatchSolution(displ,locdisp,pch,
                                 pch->getNoSpaceDim(),ASMmxBase::itgBasis);
    else
      pch->extractNodeVec(displ,locdisp,pch->getNoSpaceDim(),-1);

    if (!pch->updateCoords(locdisp))
      return false;
  }

  return true;
}


bool SIMbase::updateGrid (const std::string& field)
{
  const RealArray* displ = this->getDependentField(field);
  if (displ) return this->updateGrid(*displ);

  std::cerr <<" *** SIMbase::updateGrid: No such field \""<< field
	    <<"\" registered for \""<< this->getName() <<"\"."<< std::endl;
  return false;
}


void SIMbase::getBoundaryNodes (int pcode, IntVec& glbNodes, Vec3Vec* XYZ) const
{
  glbNodes.clear();
  if (XYZ) XYZ->clear();

  for (const Property& prop : myProps)
    if (abs(prop.pindx) == pcode)
    {
      ASMbase* pch = this->getPatch(prop.patch);
      if (!pch) continue; // Silently ignore invalid patch index

      IntVec nodes;
      if (abs(prop.ldim)+1 == pch->getNoParamDim())
        // The boundary is of one dimension lower than the patch
        pch->getBoundaryNodes(abs(prop.lindx%10),nodes);
      else if (abs(prop.ldim)+2 == pch->getNoParamDim())
        // The boundary is of two dimensions lower than the patch
        pch->getBoundary1Nodes(prop.lindx,nodes);
      else if (abs(prop.ldim) == pch->getNoParamDim())
        // The boundary and the patch are of same dimension
        for (size_t node = 1; node <= pch->getNoNodes(); node++)
          glbNodes.push_back(pch->getNodeID(node));
      else if (prop.ldim == 4)
        // The boundary nodes are stored explicitly
        nodes = pch->getNodeSet(prop.lindx);
      else
        continue; // Silently ignore other/invalid dimension specification

      // Ensure the global node numbers are unique if more than one patch
      for (int n : nodes)
        if (std::find(glbNodes.begin(),glbNodes.end(),n) == glbNodes.end())
        {
          glbNodes.push_back(n);
          if (XYZ)
          {
            size_t node = pch->getNodeIndex(n,true);
            XYZ->push_back(node ? pch->getCoord(node) : Vec3());
          }
        }
    }

#ifdef SP_DEBUG
  std::cout <<"Boundary nodes associated with property code "<< pcode <<":";
  if (glbNodes.empty())
    std::cout <<" (none)";
  else for (size_t i = 0; i < glbNodes.size(); i++)
    std::cout << (i%10 ? " " : "\n") << glbNodes[i];
  std::cout << std::endl;
#endif
}


int SIMbase::findClosestNode (const Vec3& X) const
{
  if (myModel.empty()) return -1;

  ASMbase* closestPch = nullptr;
  std::pair<size_t,double> closest(0,1.0e99);
  for (ASMbase* pch : myModel)
  {
    std::pair<size_t,double> node = pch->findClosestNode(X);
    if (node.first > 0 && node.second < closest.second)
    {
      closest = node;
      closestPch = pch;
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

  // Lambda function for assembling the interior terms for a given patch
  auto&& assembleInterior = [this,time,prevSol](IntegrandBase* integrand,
                                                GlobalIntegral& integral,
                                                ASMbase* pch, int pidx)
  {
    if (!integral.haveContributions(pidx,myProps))
      return true;

    if (msgLevel > 1)
      IFEM::cout <<"\nAssembling interior matrix terms for P"<< pidx
                 << std::endl;

    if (!this->initBodyLoad(pidx))
      return false;

    if (!this->extractPatchSolution(integrand,prevSol,pidx-1))
      return false;

    if (myProblem->getExtractionField())
    {
      if (dualField && dualField->initPatch(pch->idx))
      {
        Matrix extrField(*myProblem->getExtractionField());
        if (!pch->L2projection(extrField,dualField))
          return false;
      }
      else
        myProblem->getExtractionField()->clear();
    }

    integrand->initPatch(pch->idx);
    if (mySol)
      mySol->initPatch(pch->idx);

#if SP_DEBUG > 2
    integrand->printSolution(std::cout,pch->idx+1);
#endif
    if (!pch->integrate(*integrand,integral,time))
      return false;

    if (!(integrand->getIntegrandType() & IntegrandBase::INTERFACE_TERMS))
      return true;

    ASM::InterfaceChecker* iChk = this->getInterfaceChecker(pch->idx);
    if (!iChk) return true;

    bool ok = pch->integrate(*integrand,integral,time,*iChk);
    delete iChk;
    return ok;
  };

  bool ok = true;
  bool isAssembling = (myProblem->getMode() > SIM::INIT &&
                       myProblem->getMode() < SIM::RECOVERY);
  if (isAssembling && myEqSys && mdFlag <= 1)
    myEqSys->initialize(newLHSmatrix);

  // Loop over the integrands
  IntegrandMap::const_iterator it;
  for (it = myInts.begin(); it != myInts.end() && ok; ++it)
  {
    if (msgLevel > 1)
      IFEM::cout <<"\n\nProcessing integrand associated with code "<< it->first
                << std::endl;

    GlobalIntegral& sysQ = it->second->getGlobalInt(myEqSys);
    if (&sysQ != myEqSys && isAssembling && mdFlag <= 1)
      sysQ.initialize(newLHSmatrix);

    if (isAssembling && mdFlag <= 1)
      it->second->initLHSbuffers(newLHSmatrix);

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
        {
          if (!(pch = this->getPatch(p->patch)))
          {
            std::cerr <<" *** SIMbase::assembleSystem: Patch index "<< p->patch
                      <<" out of range [1,"<< myModel.size() <<"]."<< std::endl;
            ok = false;
          }
          else if (this->initMaterial(p->pindx))
          {
            if (p->ldim == 5)
              it->second->activateElmGroup(pch->getElementSet(p->lindx));
            else
              it->second->activateElmGroup();
            lp = p->patch;
            ok = assembleInterior(it->second,sysQ,pch,lp);
          }
          else
            ok = false;
        }

      if (lp == 0 && it->first == 0)
        // All patches refer to the same material, and we assume it has been
        // initialized during input processing (thus no initMaterial call here)
        for (lp = 0; lp < myModel.size() && ok; lp++)
          ok = assembleInterior(it->second,sysQ,myModel[lp],1+lp);
    }

    // Assemble contributions from the Neumann boundary conditions
    // and other boundary integrals (Robin properties, contact, etc.)
    it->second->activateElmGroup();
    if (it->second->hasBoundaryTerms())
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

          it->second->initPatch(pch->idx);
          if (mySol)
            mySol->initPatch(pch->idx);

          for (p2 = myProps.begin(); p2 != myProps.end() && ok; ++p2)
            if (p2->pcode == Property::MATERIAL && p->patch == p2->patch)
            {
              ok = this->initMaterial(p2->pindx);
              break;
            }

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
                           << (int)p->lindx <<" on P"<< p->patch << std::endl;
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
    if (ok && msgLevel > 1) IFEM::cout <<"\nDone."<< std::endl;
    if (ok && &sysQ != myEqSys && isAssembling && mdFlag%2 == 0)
      ok = sysQ.finalize(newLHSmatrix);
  }
  if (ok && isAssembling && myEqSys && mdFlag%2 == 0)
    ok = myEqSys->finalize(newLHSmatrix);

  if (!ok)
    std::cerr <<" *** SIMbase::assembleSystem: Failure.\n"<< std::endl;

  return ok;
}


bool SIMbase::extractLoadVec (Vector& loadVec, size_t idx, const char* hd) const
{
  if (!myEqSys || !mySam)
    return false;

  // Expand load vector from equation ordering to DOF-ordering
  SystemVector* b = myEqSys->getVector(idx);
  if (!b || !mySam->expandSolution(*b,loadVec,0.0))
    return false;

#if SP_DEBUG > 1
  mySam->printVector(std::cout,loadVec,"\nLoad vector");
#endif
  if (hd)
  {
    Vec3 sumLoad;
    for (int inod = 1; inod <= mySam->getNoNodes(); inod++)
    {
      std::pair<int,int> dofs = mySam->getNodeDOFs(inod);
      int idof = dofs.first-1;
      for (int i = 0; i < 3 && idof < dofs.second; i++, idof++)
        sumLoad[i] += loadVec[idof];
    }
    IFEM::cout <<"\nSum "<< hd <<" : "<< sumLoad << std::endl;
  }

  return true;
}


bool SIMbase::extractScalars (RealArray& values) const
{
  if (!myEqSys) return false;

  values = myEqSys->getScalars();
  for (double v : values)
    if (fabs(v) > utl::zero_print_tol)
      return true;

  return false; // All scalars are zero
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
                           const char* compName, size_t idxRHS)
{
  if (!myEqSys) return false;

  SystemMatrix* A = myEqSys->getMatrix();
  SystemVector* b = myEqSys->getVector(idxRHS);
  if (!A) std::cerr <<" *** SIMbase::solveSystem: No LHS matrix."<< std::endl;
  if (!b) std::cerr <<" *** SIMbase::solveSystem: No RHS vector."<< std::endl;
  if (!A || !b) return false;

  // Dump equation system to file(s) if requested
  this->dumpEqSys();

  // Solve the linear system of equations
  if (msgLevel > 1)
    IFEM::cout <<"\nSolving the equation system ..."<< std::endl;

  double rcn = 1.0;
  utl::profiler->start("Equation solving");
  bool status = A->solve(*b, msgLevel > 1 ? &rcn : rCond);
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
  if (status && mySam)
    status = mySam->expandSolution(*b, solution, idxRHS == 0 ? 1.0 : 0.0);
  else
    status = false;

#if SP_DEBUG > 2
  if (printSol < 1000) printSol = 1000;
#endif
  if (printSol > 0 && status)
    this->printSolutionSummary(solution,printSol,compName);

  return status;
}


bool SIMbase::solveSystem (Vectors& solution, int printSol, const char* cmpName)
{
  size_t nSol = myEqSys ? myEqSys->getNoRHS() : 0;
  if (solution.size() < nSol)
    solution.resize(nSol);

  bool status = nSol > 0;
  for (size_t i = 0; i < nSol && status; i++)
    status = this->solveSystem(solution[i],printSol,nullptr,cmpName,i);

  return status;
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
    mySam->printVector(adm.cout,solution,"Solution vector");
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


SystemMatrix* SIMbase::getRayleighDampingMatrix (size_t iM, size_t iK) const
{
  SystemMatrix* Cmat = nullptr;
  if (!myEqSys || !myProblem)
    return Cmat;

  const double eps = 1.0e-8;

  double alpha1 = myProblem->getIntegrationPrm(0);
  double alpha2 = myProblem->getIntegrationPrm(1);

  if (alpha2 > eps) // We have stiffness-proportional damping
  {
    Cmat = this->getLHSmatrix(iK,true);
    if (!Cmat)
      std::cerr <<" *** SIMbase::getRayleighDampingMatrix: No stiffness matrix "
                << iK << std::endl;
    else
    {
      Cmat->mult(alpha2);
      if (alpha1 > eps) // We have mass-proportional damping also
      {
        SystemMatrix* Mmat = this->getLHSmatrix(iM,false);
        if (!Mmat)
        {
          std::cerr <<" *** SIMbase::getRayleighDampingMatrix: No mass matrix "
                    << iM << std::endl;
          delete Cmat;
          Cmat = nullptr;
        }
        else if (!Cmat->add(*Mmat,alpha1))
        {
          std::cerr <<" *** SIMbase::getRayleighDampingMatrix: Incompatible"
                    <<" mass- and stiffness matrices"<< std::endl;
          delete Cmat;
          Cmat = nullptr;
        }
      }
    }
  }
  else if (alpha1 > eps) // Only mass-proportional damping
  {
    Cmat = this->getLHSmatrix(iM,true);
    if (!Cmat)
      std::cerr <<" *** SIMbase::getRayleighDampingMatrix: No mass matrix "
                << iM << std::endl;
    else
      Cmat->mult(alpha1);
  }

  return Cmat;
}


SystemMatrix* SIMbase::getLHSmatrix (size_t idx, bool copy) const
{
  if (!myEqSys) return nullptr;

  SystemMatrix* lhs = myEqSys->getMatrix(idx);
  return lhs && copy ? lhs->copy() : lhs;
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

  size_t nCmp = 0;

  // Lambda function for assembling the interior norm terms for a given patch
  auto&& assembleNorms = [this,time,psol,ssol,&nCmp](NormBase* norm,
                                                     GlbNorm& integral,
                                                     ASMbase* pch, int pidx)
  {
    if (!integral.haveContributions(pidx,myProps))
      return true;

    if (!this->extractPatchSolution(psol,pidx-1))
      return false;

    bool ok = true;
    Vector* extrVec = myProblem->getExtractionField();
    if (extrVec && extrFunc.size() == 1)
    {
      if (extrFunc.front()->initPatch(pch->idx))
      {
        Matrix extrField(*extrVec);
        ok = pch->L2projection(extrField,extrFunc.front());
      }
      else
        extrVec->clear();
    }
    else if (extrFunc.size() > 1)
    {
      std::vector<FunctionBase*> activeFunc;
      std::vector<Matrix*>       extrFields;
      activeFunc.reserve(extrFunc.size());
      extrFields.reserve(extrFunc.size());

      for (size_t k = 1; (extrVec = myProblem->getExtractionField(k)); k++)
        if (extrFunc[k-1]->initPatch(pch->idx))
        {
          activeFunc.push_back(extrFunc[k-1]);
          extrFields.push_back(new Matrix(*extrVec));
        }
        else
          extrVec->clear();

      if (!activeFunc.empty())
      {
        ok = pch->L2projection(extrFields,activeFunc);
        for (Matrix* m : extrFields) delete m;
      }
    }

    size_t nfld = myProblem->getNoFields(2);
    size_t nval = pch->getNoProjectionNodes()*nfld;
    for (size_t k = 0; k < ssol.size(); k++)
      if (ssol[k].empty())
        norm->getProjection(k).clear();
      else if (this->fieldProjections())
      {
        Vector c(nval);
        std::copy(ssol[k].begin()+nCmp,ssol[k].begin()+nCmp+nval,c.begin());
        norm->setProjectedFields(pch->getProjectedFields(c,nfld),k);
        nCmp += nval;
      }
      else if (nCmp != pch->getNoFields(1) && pch->getNoFields(2) > 0)
      {
        // Mixed problem using first basis for projection,
        // need a separate MADOF array with nCmp dofs per node
        const IntVec& madof = this->getMADOF(1,nCmp);
        pch->extractNodalVec(ssol[k],norm->getProjection(k),
                             madof.data(),madof.size());
      }
      else
        pch->extractNodeVec(ssol[k],norm->getProjection(k),nCmp,1);

    norm->initPatch(pch->idx);
    if (mySol)
      mySol->initPatch(pch->idx);

    if (!pch->integrate(*norm,integral,time))
      return false;

    if (!(norm->getIntegrandType() & IntegrandBase::INTERFACE_TERMS))
      return ok;

    ASM::InterfaceChecker* iChk = this->getInterfaceChecker(pch->idx);
    if (!iChk) return ok;

    ok &= pch->integrate(*norm,integral,time,*iChk);
    delete iChk;
    return ok;
  };

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
  size_t nNodes = this->getNoNodes(1);
  if (nNodes > 0 && !this->fieldProjections())
    for (const Vector& s : ssol)
      if ((nCmp = s.size() / nNodes) > 0)
        break;

#ifdef USE_OPENMP
  // When assembling in parallel, we must always do the norm summation
  // at the end in a serial loop, to avoid that the threads try to update
  // the same memory address simultaneously.
  Matrix dummy;
  if (!eNorm) eNorm = &dummy;
#endif

  // Initialize norm integral classes
  if (!extrFunc.empty())
    gNorm.reserve(norm->getNoFields(0)+1);
  gNorm.resize(norm->getNoFields(0));
  size_t nNorms = 0;
  SIMoptions::ProjectionMap::const_iterator proj_idx = opt.project.begin();
  // count norms if they are:
  // 1) associated with the primary solution
  // 2) associated with a projected secondary solution
  // 3) associated with a residual secondary solution
  // 4) associated with a secondary solution provided through external means
  // 5) associated with an additional group defined in the integrand
  for (size_t i = 0; i < gNorm.size(); i++)
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
    for (size_t c = 0; c < eNorm->cols(); c++)
      elementNorms.push_back(new ElmNorm(eNorm->ptr(c),nNorms));
    norm->setLocalIntegrals(&elementNorms);
  }

  // Loop over the different material regions, integrating solution norm terms
  // for the patch domain associated with each material
  bool ok = true;
  size_t lp = 0;
  ASMbase* pch = nullptr;
  PropertyVec::const_iterator p;
  for (p = myProps.begin(); p != myProps.end() && ok; ++p)
    if (p->pcode == Property::MATERIAL)
      if (!(pch = this->getPatch(p->patch)))
        ok = false;
      else if (this->initMaterial(p->pindx))
      {
        if (p->ldim == 5)
          norm->activateElmGroup(pch->getElementSet(p->lindx));
        else
          norm->activateElmGroup();
        lp = p->patch;
        ok = assembleNorms(norm,globalNorm,pch,lp);
      }
      else
        ok = false;

  if (lp == 0)
    // All patches are referring to the same material, and we assume it has
    // been initialized during input processing (thus no initMaterial call here)
    for (lp = 0; lp < myModel.size() && ok; lp++)
      ok = assembleNorms(norm,globalNorm,myModel[lp],1+lp);

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
            norm->initPatch(pch->idx);
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
            norm->initPatch(pch->idx);
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
  for (LocalIntegral* elmNorm : elementNorms)
  {
    globalNorm.assemble(elmNorm);
    delete elmNorm;
  }

  // Add problem-dependent external norm contributions
  norm->addBoundaryTerms(gNorm,this->externalEnergy(psol,time));

  if (!extrFunc.empty()) // Make room for extraction function domain volumes
    gNorm.resize(norm->getNoFields(0)+1,Vector(extrFunc.size()));

  delete norm;

  if (!adm.dd.isPartitioned())
    for (Vector& glbNorm : gNorm)
      adm.allReduceAsSum(glbNorm);

  return ok && this->postProcessNorms(gNorm,eNorm);
}


double SIMbase::externalEnergy (const Vectors& psol, const TimeDomain&) const
{
  const RealArray* reactionForces = this->getReactionForces();
  if (!reactionForces || !mySam || psol.empty()) return 0.0;

  // Add norm contributions due to inhomogeneous Dirichlet boundary conditions.
  // That is, the path integral of the total solution vector times the
  // reaction forces at the prescribed DOFs.
  if (psol.size() == 1)
    return mySam->normReact(psol.front(),*reactionForces);

  if (prevForces.empty())
    extEnergy += 0.5*mySam->normReact(psol[0]-psol[1],*reactionForces);
  else
    extEnergy += 0.5*mySam->normReact(psol[0]-psol[1],
                                      prevForces.add(*reactionForces));
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
  const RealArray* reactionForces = this->getReactionForces();
  if (!reactionForces || !mySam) return false;

  RF.resize(1+nsd);
  RF.front() = mySam->normReact(psol,*reactionForces);
  for (unsigned char dir = 1; dir <= nsd; dir++)
    RF[dir] = mySam->getReaction(dir,*reactionForces);

#if SP_DEBUG > 1
  std::cout <<"\nHere are the nodal reaction forces:";
  Vector nrf;
  for (int inod = 1; inod <= mySam->getNoNodes(); inod++)
    if (mySam->getNodalReactions(inod,*reactionForces,nrf))
    {
      std::cout <<"\nNode "<< inod <<":";
      for (double f : nrf) IFEM::cout <<" "<< f;
    }
  std::cout << std::endl;
#endif
  return true;
}


bool SIMbase::getCurrentReactions (RealArray& RF, int pcode) const
{
  const RealArray* reactionForces = this->getReactionForces();
  if (!reactionForces || !mySam) return false;

  IntVec glbNodes;
  this->getBoundaryNodes(pcode,glbNodes);

  RF.resize(nsd);
  for (unsigned char dir = 1; dir <= nsd; dir++)
    RF[dir-1] = mySam->getReaction(dir,*reactionForces,&glbNodes);

#if SP_DEBUG > 1
  std::cout <<"\nReaction forces associated with property code "<< pcode;
  Vector nrf;
  for (int inod : glbNodes)
    if (mySam->getNodalReactions(inod,*reactionForces,nrf))
    {
      std::cout <<"\nNode "<< inod <<":";
      for (double f : nrf) IFEM::cout <<" "<< f;
    }
  std::cout << std::endl;
#endif
  return true;
}


bool SIMbase::haveReactions (int pcode) const
{
  IntVec glbNodes;
  if (pcode > 0)
    this->getBoundaryNodes(pcode,glbNodes);

  for (unsigned char dir = 1; dir <= nsd; dir++)
    if (mySam->haveReaction(dir, pcode > 0 ? &glbNodes : nullptr))
      return true;

  return false;
}


const RealArray* SIMbase::getReactionForces () const
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

  // Dump equation system to file(s) if requested
  this->dumpEqSys();

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
    solution[i-1].eqnVec = eigVec.getColumn(i);
    if (!mySam->expandVector(solution[i-1].eqnVec,solution[i-1].eigVec))
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


bool Mode::orthonormalize (const SystemMatrix& mat)
{
  StdVector tmp;
  if (!mat.multiply(StdVector(eqnVec),tmp))
    return false;

  double scale = tmp.dot(eqnVec);
  if (scale < 1.0e-16)
  {
    std::cerr <<" *** Mode::orthonormalize: Mode shape "<< eigNo
              <<" is zero."<< std::endl;
    return false;
  }

  IFEM::cout <<"  * Mode shape "<< eigNo;
  if (fabs(scale-1.0) < 1.0e-16)
    IFEM::cout <<" is already orthonormal."<< std::endl;
  else
  {
    eqnVec /= scale;
    eigVec /= scale;
    IFEM::cout <<": scale = "<< 1.0/scale << std::endl;
  }

  return true;
}


bool Mode::computeDamping (const SystemMatrix& mat)
{
  StdVector tmp;
  if (!mat.multiply(StdVector(eqnVec),tmp))
    return false;

  damping = tmp.dot(eqnVec);
  IFEM::cout <<"  * Mode "<< eigNo <<": damping = "<< damping << std::endl;

  return true;
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

  // Initialize result point buffers within the integrand (if any).
  // The large negative time argument is here used to flag that we are going to
  // do numerical integration, possibly using integration point buffers instead
  // of the result evaluation buffers.
  myProblem->initResultPoints(-9999.0);

  size_t i, ofs = 0;
  for (i = 0; i < myModel.size(); i++)
  {
    if (myModel[i]->empty()) continue; // skip empty patches

    // Extract the primary solution control point values for this patch
    if (!this->extractPatchSolution(myProblem,Vectors(1,psol),i))
      return false;

    // Initialize material properties for this patch in case of multiple regions
    this->setPatchMaterial(i+1);

    LocalSystem::patch = i; // Hack: Used for patch-wise max-value calculation

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
      ok = myModel[i]->globalL2projection(values,L2ProbIntegrand(*myModel[i], *myProblem));
      break;

    case SIMoptions::CGL2:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tContinuous global L2-projection"<< std::endl;
      GlbL2::SolverParams = myGl2Params;
      ok = myModel[i]->globalL2projection(values,L2ProbIntegrand(*myModel[i], *myProblem),true);
      break;

    case SIMoptions::CGL2_INT:
      if (msgLevel > 1 && i == 0)
        IFEM::cout <<"\tContinuous global L2-projection"<< std::endl;
      GlbL2::SolverParams = myGl2Params;
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


bool SIMbase::project (RealArray& values, const FunctionBase* f,
                       int basis, int iField, int nFields,
                       SIMoptions::ProjectionMethod pMethod, double time) const
{
  bool ok = true;
  for (size_t j = 0; j < myModel.size() && ok; j++)
  {
    if (myModel[j]->empty()) continue; // skip empty patches

    Vector loc_values;
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
        if (myModel.size() > 1) {
          std::cerr <<"  ** L2 projection of explicit functions onto a"
                    <<" separate basis is not available for multi-patch models."
                    << std::endl;
          return false;
        }
        Matrix ftmp(loc_values);
        ok = myModel[j]->globalL2projection(ftmp,L2FuncIntegrand(*myModel[j], *f),true);
        values = loc_values;
        return ok;
      }
      else
      {
        Matrix ftmp(loc_values);
        ok = myModel[j]->L2projection(ftmp,const_cast<FunctionBase*>(f),time);
      }
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
  for (size_t i = 0; i < problem->getNoSolutions(); i++)
    if (i < sol.size() && !sol[i].empty())
      pch->extractNodalVec(sol[i],problem->getSolution(i),mySam->getMADOF());
    else
      problem->getSolution(i).clear();

  return this->extractPatchDependencies(problem,myModel,pindx);
}


bool SIMbase::project (Vector& ssol, const Vector& psol,
                       SIMoptions::ProjectionMethod pMethod, size_t iComp) const
{
  if (iComp > 0)
  {
    Matrix stmp;
    bool ok = this->project(stmp,psol,pMethod);
    ssol = stmp.getRow(iComp);
    return ok;
  }

  Matrix stmp(ssol);
  return this->project(stmp,psol,pMethod);
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


size_t SIMbase::extractPatchSolution (const RealArray& sol, RealArray& vec,
                                      const ASMbase* pch, unsigned char nndof,
                                      unsigned char basis) const
{
  if (!pch || sol.empty()) return 0;

  unsigned char ncmp = pch->getNoFields(basis);
  if (basis && nndof > 0 && nndof != ncmp && pch->getNoFields(2) > 0)
  {
    // Mixed problem, and the incoming global vector field has more (or less)
    // components per node than used on this basis for the primary solution.
    // ==> Need to use a separate MADOF array (pre-computed) on this basis.
    const IntVec& madof = this->getMADOF(basis,nndof);
    pch->extractNodalVec(sol,vec,madof.data(),madof.size());
  }
  else if (mySam && (mySam->getNoNodes('X') > 0 || this->getMDflag() > 0))
  {
    // This is a model with extraordinary nodes,
    // and/or with monolithic coupled simulators
    const int* madof = mySam->getMADOF();
    if (nndof > 0 && nndof != ncmp)
    {
      // Construct a separate MADOF array with nndof DOFs per node
      IntVec& MADOF = extraMADOFs[-nndof];
      if (MADOF.empty())
      {
        MADOF.resize(mySam->getNoNodes()+1,0);
        for (const ASMbase* pch : myModel)
          for (size_t inod = 1; inod <= pch->getNoNodes(); inod++)
          {
            int node = pch->getNodeID(inod);
            if (node > 0 && pch->getNodeType(inod) == 'D')
              MADOF[node] = nndof;
          }

        MADOF[0] = 1;
        for (size_t n = 1; n < MADOF.size(); n++)
          MADOF[n] += MADOF[n-1];

      }
      madof = MADOF.data();
    }

    // Excluding the extraordinary DOFs
    pch->extractNodalVec(sol,vec,madof,-2);
  }
  else
    pch->extractNodeVec(sol,vec,nndof,basis);

  return vec.size();
}


bool SIMbase::injectPatchSolution (RealArray& sol, const RealArray& vec,
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

  this->setPatchMaterial(pindx+1);
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


bool SIMbase::setPatchMaterial (size_t patch) const
{
  for (const Property& p : myProps)
    if (p.pcode == Property::MATERIAL && p.patch == patch)
      return const_cast<SIMbase*>(this)->initMaterial(p.pindx);

  return false;
}


bool SIMbase::addMADOF (unsigned char basis, unsigned char nndof, bool other)
{
  int key = basis*65536 + nndof;
  if (extraMADOFs.find(key) != extraMADOFs.end())
    return false; // This MADOF already calculated

  IntVec& madof = extraMADOFs[key];
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
  int key = basis*65536 + nndof;
  auto it = extraMADOFs.find(key);
  if (it != extraMADOFs.end())
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


void SIMbase::dumpEqSys ()
{
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

  if (myEqSys->getNoRHS() == 0)
    return;

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
}


bool SIMbase::assembleForces (const Vector& solution, double t0,
                              Vector* R, Vector* S)
{
  if (!myProblem) return false;

  // Assign secondary integral to the integrand
  myProblem->setSecondaryInt(new ReactionsOnly(mySam,adm,R,S));

  // Temporarily nullify myEqSys such that the secondary integral will be used
  AlgEqSystem* tmpEqSys = myEqSys;
  myEqSys = nullptr;

  // Integrate the reaction and/or interface forces
  bool ok = this->setMode(SIM::RHS_ONLY) && this->assembleSystem(t0,{solution});

  // Release the secondary integral and restore the system matrix pointer
  myProblem->setSecondaryInt();
  myEqSys = tmpEqSys;

  return ok;
}
