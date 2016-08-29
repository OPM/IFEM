// $Id$
//==============================================================================
//!
//! \file ASMbase.C
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for spline-based finite element (FE) assembly drivers.
//!
//==============================================================================

#include "ASMbase.h"
#include "ASM2D.h"
#include "ASM3D.h"
#include "IFEM.h"
#include "MPC.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include <algorithm>

#ifdef USE_OPENMP
#include <omp.h>
#endif


bool ASMbase::fixHomogeneousDirichlet = true;


/*!
  \brief Convenience function writing error message for non-implemented methods.
*/

static bool Aerror (const char* name)
{
  std::cerr <<" *** ASMbase::"<< name
	    <<": Must be implemented in sub-class."<< std::endl;
  return false;
}


ASMbase::ASMbase (unsigned char n_p, unsigned char n_s, unsigned char n_f)
  : MLGE(myMLGE), MLGN(myMLGN), MNPC(myMNPC), shareFE(0)
{
  nf = n_f;
  nsd = n_s > 3 ? 3 : n_s;
  ndim = n_p > nsd ? nsd : n_p;
  nLag = 0;
  nGauss = 0;
  nel = nnod = 0;
  idx = 0;
  firstIp = 0;
  myLMs.first = myLMs.second = 0;
}


ASMbase::ASMbase (const ASMbase& patch, unsigned char n_f)
  : MLGE(patch.MLGE), MLGN(patch.MLGN), MNPC(patch.MNPC), shareFE('F')
{
  nf = n_f > 0 ? n_f : patch.nf;
  nsd = patch.nsd;
  ndim = patch.ndim;
  nLag = patch.nLag;
  nGauss = patch.nGauss;
  nel = patch.nel;
  nnod = patch.nnod;
  idx = patch.idx;
  firstIp = patch.firstIp;
  firstBp = patch.firstBp;
  myLMs = patch.myLMs;
  myLMTypes = patch.myLMTypes;
  // Note: Properties are _not_ copied
}


ASMbase::ASMbase (const ASMbase& patch)
  : MLGE(myMLGE), MLGN(myMLGN), MNPC(myMNPC), shareFE('S')
{
  nf = patch.nf;
  nsd = patch.nsd;
  ndim = patch.ndim;
  nGauss = patch.nGauss;
  nel = patch.nel;
  nnod = patch.nnod;
  idx = patch.idx;
  firstIp = patch.firstIp;
  firstBp = patch.firstBp;
  BCode = patch.BCode;

  // Only copy the regular part of the FE data, leave out any extraordinaries

  if (patch.MLGE.size() > nel)
    myMLGE.insert(myMLGE.begin(),patch.MLGE.begin(),patch.MLGE.begin()+nel);
  else
    myMLGE = patch.MLGE;

  if (patch.MLGN.size() > nnod)
    myMLGN.insert(myMLGN.begin(),patch.MLGN.begin(),patch.MLGN.begin()+nnod);
  else
    myMLGN = patch.MLGN;

  if (patch.MNPC.size() > nel)
    myMNPC.insert(myMNPC.begin(),patch.MNPC.begin(),patch.MNPC.begin()+nel);
  else
    myMNPC = patch.MNPC;

  // Can not copy pointers as it might cause problems on destruction
  if (!patch.dCode.empty() || !patch.mpcs.empty())
    std::cerr <<"  ** ASMbase copy constructor: The copied patch has"
              <<" multi-point constraints, these are not copied.\n";

  nLag = 0; // Lagrange multipliers are not copied
  myLMs.first = myLMs.second = 0;
}


ASMbase::~ASMbase ()
{
  for (MPCIter it = mpcs.begin(); it != mpcs.end(); it++)
    delete *it;
}


ASMbase* ASMbase::cloneUnShared() const
{
  const ASM2D* patch2 = dynamic_cast<const ASM2D*>(this);
  if (patch2) return patch2->clone();

  const ASM3D* patch3 = dynamic_cast<const ASM3D*>(this);
  if (patch3) return patch3->clone();

  return nullptr;
}


void ASMbase::clear (bool retainGeometry)
{
  if (retainGeometry)
  {
    // Clear all FE structures, including the elements
    myMLGE.clear();
    myMNPC.clear();
  }
  else // Don't erase the elements, but set them to have zero nodes
    for (size_t i = 0; i < myMNPC.size(); i++) myMNPC[i].clear();

  // Erase the nodes, boundary conditions and multi-point constraints
  for (MPCIter it = mpcs.begin(); it != mpcs.end(); it++)
    delete *it;

  myLMs.first = myLMs.second = 0;
  myLMTypes.clear();

  myMLGN.clear();
  BCode.clear();
  dCode.clear();
  mpcs.clear();
}


bool ASMbase::addXElms (short int, short int, size_t, std::vector<int>&)
{
  return Aerror("addXElms(short int,short int,size_t,std::vector<int>&)");
}


bool ASMbase::addLagrangeMultipliers (size_t iel, const IntVec& mGLag,
                                      unsigned char nnLag)
{
  if (iel > MNPC.size())
  {
    std::cerr <<" *** ASMbase::addLagrangeMultipliers: Element index "<< iel
              <<" is out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
  else if (shareFE == 'F')
    return false;

  if (nLag == 0 || iel == 0)
    nLag = nnLag;
  else if (nnLag != nLag)
    return false;

  for (size_t i = 0; i < mGLag.size(); i++)
  {
    size_t node = MLGN.size();
    IntVec::const_iterator it = std::find(MLGN.begin(),MLGN.end(),mGLag[i]);
    if (it == MLGN.end())
      myMLGN.push_back(mGLag[i]); // Add a new Lagrange multiplier node
    else
      node = it - MLGN.begin(); // Existing Lagrange multiplier node

    // Update the nodal range (1-based indices) of the Lagrange multipliers
    if (myLMs.first == 0)
      myLMs.first = myLMs.second = node+1;
    else if (node+1 < myLMs.first)
    {
      std::cerr <<" *** ASMbase::addLagrangeMultipliers: Node "<< node+1
                <<" is out of range ["<< myLMs.first <<","<< myLMs.second
                <<"]."<< std::endl;
      return false;
    }
    else if (node >= myLMs.second)
      myLMs.second = node+1;

    if (myLMTypes.size() < node-myLMs.first+2)
      myLMTypes.resize(node-myLMs.first+2);

    myLMTypes[node+1-myLMs.first] = iel == 0 ? 'G' : 'L';

    // Extend the element connectivity table
    if (iel > 0)
      myMNPC[iel-1].push_back(node);
    else for (auto& it : myMNPC)
      it.push_back(node);
  }

  return true;
}


bool ASMbase::addGlobalLagrangeMultipliers (const IntVec& mGLag,
                                            unsigned char nnLag)
{
#ifdef USE_OPENMP
  if (omp_get_max_threads() > 1) {
    std::cerr << "** Cannot do multi-threaded assembly with global multipliers. **" << std::endl
              << "\t Setting OMP_NUM_THREADS = 1" << std::endl;
    omp_set_max_threads(1);
  }
#endif

  return this->addLagrangeMultipliers(0,mGLag,nnLag);
}


size_t ASMbase::getNodeIndex (int globalNum, bool) const
{
  IntVec::const_iterator it = std::find(MLGN.begin(),MLGN.end(),globalNum);
  if (it == MLGN.end()) return 0;

  return 1 + (it-MLGN.begin());
}


int ASMbase::getNodeID (size_t inod, bool) const
{
  if (inod < 1 || inod > MLGN.size())
    return 0;

  return MLGN[inod-1];
}


char ASMbase::getLMType(size_t n) const
{
  if (n >= myLMs.first && n <= myLMs.second)
    return myLMTypes[n-myLMs.first];

  return 'L';
}


int ASMbase::getElmID (size_t iel) const
{
  if (iel < 1 || iel > MLGE.size())
    return 0;

  return abs(MLGE[iel-1]);
}


unsigned char ASMbase::getNodalDOFs (size_t inod) const
{
  return this->isLMn(inod) ? nLag : nf;
}


char ASMbase::getNodeType (size_t inod) const
{
  return this->isLMn(inod) ? this->getLMType(inod) : (inod > nnod ? 'X' : 'D');
}


size_t ASMbase::getNoNodes (int basis) const
{
  if (basis > 0)
    return nnod;
  else if (basis < 0 && myLMs.first > 0)
    return myLMs.first - 1;
  else
    return MLGN.size();
}


size_t ASMbase::getNoElms (bool includeZeroVolElms, bool includeXElms) const
{
  if (includeZeroVolElms)
    return nel;

  size_t numels = 0;
  for (size_t i = 0; i < MLGE.size(); i++)
    if (MLGE[i] > 0 || (includeXElms && MLGE[i] < 0))
      numels++;

  return numels;
}


void ASMbase::getNoIntPoints (size_t& nPt, size_t& nIPt)
{
  size_t nGp = 1;
  for (unsigned char d = 0; d < ndim; d++)
    nGp *= nGauss;

  firstIp = nPt;
  nPt += nel*nGp; // Note: Includes also the 0-span elements

  // Count additional interface quadrature points
  size_t nInterface = MLGE.size() - nel;
  if (nInterface > 0 && nInterface != nel && nGauss > 0)
    nIPt += nInterface*nGp/nGauss;
}


void ASMbase::getNoBouPoints (size_t& nPt, char ldim, char lindx)
{
  size_t nGp = 1;
  for (char d = 0; d < ldim; d++)
    nGp *= nGauss;

  firstBp[lindx] = nPt;

  nPt += this->getNoBoundaryElms(lindx,ldim)*nGp; // Includes 0-span elements
}


void ASMbase::printNodes (std::ostream& os, const char* heading) const
{
  Matrix X;
  this->getNodalCoordinates(X);
  if (heading) os << heading;
  for (size_t inod = 1; inod <= X.cols(); inod++)
  {
    os <<'\n'<< inod <<' '<< MLGN[inod-1] <<':';
    for (size_t i = 1; i <= X.rows(); i++)
      os <<' '<< X(i,inod);
  }
  os << std::endl;
}


/*!
  \brief A helper class used by ASMbase::isFixed.
  \details The class is just an unary function that checks whether a DOF object
  matches the fixed status of a given BC object.
*/

class fixed
{
  int myNode; //!< The internal node number to compare with
  int myDofs; //!< The local DOFs to compare with

public:
  //! \brief Constructor initializing the node and local dof index.
  fixed(int node, int dof) : myNode(node), myDofs(dof) {}
  //! \brief Returns \e true if the DOF has the same fixed status as \a bc.
  bool operator()(const ASMbase::BC& bc)
  {
    if (bc.node == myNode)
      for (int dof = myDofs; dof > 0; dof /= 10)
	switch (dof%10)
	  {
	  case 1: return bc.CX == 0;
	  case 2: return bc.CY == 0;
	  case 3: return bc.CZ == 0;
	  case 4: return bc.RX == 0;
	  case 5: return bc.RY == 0;
	  case 6: return bc.RZ == 0;
	  }
    return false;
  }
};


bool ASMbase::isFixed (int node, int dof, bool all) const
{
  BCVec::const_iterator bit = BCode.begin();
  if (dof < 10 || !all)
    bit = std::find_if(BCode.begin(),BCode.end(),fixed(node,dof));
  else for (int d = dof; d > 0 && bit != BCode.end(); d /= 10)
    if (d <= nf)
      bit = std::find_if(BCode.begin(),BCode.end(),fixed(node,d%10));

  return bit != BCode.end();
}


bool ASMbase::addMPC (MPC*& mpc, int code, bool silence)
{
  if (!mpc) return true;

  if (this->isFixed(mpc->getSlave().node,mpc->getSlave().dof))
  {
    // Silently ignore MPC's on dofs that already are marked as FIXED
    delete mpc;
    mpc = 0;
    return true;
  }

  std::pair<MPCIter,bool> mit = mpcs.insert(mpc);
  if (mit.second)
  {
#if SP_DEBUG > 1
    if (!silence) std::cout <<"Added constraint: "<< *mpc;
#endif
    if (code > 0) dCode[mpc] = code;
    return true;
  }

#ifdef SP_DEBUG
  if (!silence) std::cout <<"Ignored constraint (duplicated slave): "<< *mpc;
#endif
  delete mpc;
  mpc = *mit.first; // This dof is already a slave in another MPC
  return false;
}


bool ASMbase::add2PC (int slave, int dir, int master, int code)
{
  if (dir < 1 || dir > nf) return true;

  MPC* cons = new MPC(slave,dir);
  bool stat = this->addMPC(cons,code,true);
  if (!cons) return stat;

  cons->addMaster(master,dir);
#if SP_DEBUG > 1
  std::cout <<"Added constraint: "<< *cons;
#endif
  return stat;
}


bool ASMbase::add3PC (int slave, int dir, int master1, int master2, int code)
{
  if (dir < 1 || dir > nf) return true;

  MPC* cons = new MPC(slave,dir);
  bool stat = this->addMPC(cons,code,true);
  if (!cons) return stat;

  cons->addMaster(master1,dir);
  cons->addMaster(master2,dir);
#if SP_DEBUG > 1
  std::cout <<"Added constraint: "<< *cons;
#endif
  return stat;
}


MPC* ASMbase::findMPC (int node, int dof) const
{
  MPC slave(node,dof);
  MPCIter cit = mpcs.find(&slave);
  return cit == mpcs.end() ? nullptr : *cit;
}


bool ASMbase::addPeriodicity (size_t master, size_t slave, int dir)
{
  int slaveNode  = this->getNodeID(slave);
  int masterNode = this->getNodeID(master);
  if (slaveNode < 1 || masterNode < 1)
  {
    std::cerr <<" *** ASMbase::addPeriodicity: Invalid node indices "
	      << master <<", "<< slave << std::endl;
    return false;
  }

  if (this->add2PC(masterNode,slaveNode,dir) ||
      this->add2PC(slaveNode,masterNode,dir))
    return true;

  std::cerr <<" *** ASMbase::addPeriodicity: Failed to connect nodes "
	    << masterNode <<" and "<< slaveNode <<" in direction "
	    << dir << std::endl;
  return false;
}


void ASMbase::makePeriodic (size_t master, size_t slave, int dirs)
{
  std::set<int> dofs(utl::getDigits(dirs));
  if (dofs.size() == nf && *dofs.rbegin() == nf)
    // If all DOFs are going to be coupled, assign a common global node number
    ASMbase::collapseNodes(*this,master,*this,slave);
  else for (std::set<int>::iterator it = dofs.begin(); it != dofs.end(); ++it)
    this->addPeriodicity(master,slave,*it);
}


void ASMbase::constrainPatch (int dof, int code)
{
  if (code > 0)
    std::cerr <<"  ** ASMbase::constrainPatch: Projection onto the spline basis"
              <<" not yet implemented!"<< std::endl;
  else if (code < 0)
    code = -code;

  for (size_t node = 1; node <= this->getNoNodes(1); node++)
    this->prescribe(node,dof,code);
}


int ASMbase::prescribe (size_t inod, int dirs, int code)
{
  if (code == 0 && fixHomogeneousDirichlet)
    return this->fix(inod,dirs);

  int node = this->getNodeID(inod);
  if (node < 1 || dirs < 1) return dirs;

  int ignoredDirs = 0;
  std::set<int> dofs(utl::getDigits(dirs));
  for (std::set<int>::const_iterator it = dofs.begin(); it != dofs.end(); ++it)
    if (*it <= nf)
    {
      MPC* mpc = new MPC(node,*it);
      if (!this->addMPC(mpc,code))
        ignoredDirs = 10*ignoredDirs + *it;
    }
    else
    {
      ignoredDirs = 10*ignoredDirs + *it;
      std::cerr <<"  ** ASMbase::prescribe: Ignoring invalid DOF code "<< *it
                << std::endl;
    }

  return ignoredDirs;
}


/*!
  \brief Equality operator for BC objects comparing node numbers and DOF codes.
*/

bool operator== (const ASMbase::BC& rhs, const ASMbase::BC& lhs)
{
  return (rhs.node == lhs.node &&
          rhs.CX == lhs.CX && rhs.CY == lhs.CY && rhs.CZ == lhs.CZ &&
          rhs.RX == lhs.RX && rhs.RY == lhs.RY && rhs.RZ == lhs.RZ);
}


/*!
  \brief Equality operator for BC objects comparing node numbers only.
*/

bool operator== (const ASMbase::BC& rhs, const int& lhs)
{
  return rhs.node == lhs;
}


int ASMbase::fix (size_t inod, int dirs)
{
  int node = this->getNodeID(inod);
  if (node < 1 || dirs < 1) return dirs;

  BCVec::iterator bit = std::find(BCode.begin(),BCode.end(),node);
  if (bit == BCode.end())
  {
    BCode.push_back(BC(node));
    bit = BCode.end()-1;
  }

#if SP_DEBUG > 1
  BC old = *bit;
#endif

  int invalidDOFs = 0;
  std::set<int> dofs(utl::getDigits(dirs));
  for (std::set<int>::const_iterator it = dofs.begin(); it != dofs.end(); ++it)
    if (*it <= nf)
      switch (*it) {
      case 1: bit->CX = 0; break;
      case 2: bit->CY = 0; break;
      case 3: bit->CZ = 0; break;
      case 4: bit->RX = 0; break;
      case 5: bit->RY = 0; break;
      case 6: bit->RZ = 0; break;
      }
    else
    {
      invalidDOFs = 10*invalidDOFs + *it;
      std::cerr <<"  ** ASMbase::fix: Ignoring invalid DOF code "<< *it
                << std::endl;
    }

#if SP_DEBUG > 1
  if (!(old == *bit))
    std::cout <<"\tFixed node: "<< node <<" "<< dirs << std::endl;
#endif
  return invalidDOFs;
}


void ASMbase::mergeAndGetAllMPCs (const ASMVec& model, MPCSet& allMPCs)
{
  // Build the set of constraint equations over all patches in the model
  if (model.size() == 1)
  {
    // Trivial for single-patch models
    allMPCs.insert(model.front()->begin_MPC(),model.front()->end_MPC());
    return;
  }

  // In multi-patch models interface nodes may be constrained twice or more.
  // Resolve this such that allMPCs only contains the set of unique MPCs.
  int nmerged = 0;
  int ndeleted = 0;
  for (ASMVec::const_iterator it = model.begin(); it != model.end(); it++)
  {
    std::pair<MPCIter,bool> ret;
    std::vector<MPC*> uniqueMPC;
    uniqueMPC.reserve((*it)->getNoMPCs());
    for (MPCIter cit = (*it)->begin_MPC(); cit != (*it)->end_MPC(); cit++)
      if ((ret = allMPCs.insert(*cit)).second)
	uniqueMPC.push_back(*cit);
      else
      {
	// Merge multiple constraint equations with common slave definition
	if ((*ret.first)->getSlave().coeff == 0.0 &&
	    (*ret.first)->getNoMaster() > 1 &&
	    (*ret.first)->merge(*cit))
	{
#if SP_DEBUG > 1
	  std::cout <<"Merging constraint "<< **cit;
	  std::cout <<"Resulting constraint "<< **ret.first;
#endif
	  nmerged++;
	}
	else
	{
	  // The found constraint *ret.first is either a prescribed movement
	  // or a single-master constraint. Such MPCs are not to be merged but
	  // should superseed any multi-master constraints with matching slave.
#if SP_DEBUG > 1
	  std::cout <<"Deleted constraint "<< **cit;
#endif
	  ndeleted++;
	}
	(*it)->dCode.erase(*cit);
	delete *cit;
      }

    if (uniqueMPC.size() < (*it)->getNoMPCs())
    {
      // Compress the MPC set for this patch removing duplicated entries
      (*it)->mpcs.clear();
      (*it)->mpcs.insert(uniqueMPC.begin(),uniqueMPC.end());
    }
  }

  if (nmerged > 0)
    IFEM::cout <<"Merged "<< nmerged <<" MPC equations."<< std::endl;
  if (ndeleted > 0)
    IFEM::cout <<"Deleted "<< ndeleted <<" MPC equations."<< std::endl;

#if SP_DEBUG > 1
  if (allMPCs.empty()) return;
  std::cout <<"\nMulti-point constraints:\n";
  for (MPCIter c = allMPCs.begin(); c != allMPCs.end(); c++) std::cout << **c;
#endif
}


/*!
  If \a setPtrOnly is \e true, the MPC equations are not modified. Instead
  the pointers to the next MPC in the chain is assigned for the master DOFs
  which are slaves in other MPCs.

  This is used in time-dependent/nonlinear simulators where the constraint
  coefficients might be updated due to time variation. The resolving is then
  done directly in the MMCEQ/TTCC arrays of the SAM data object.
  \sa SAMpatch::updateConstraintEqs.
*/

void ASMbase::resolveMPCchains (const MPCSet& allMPCs, bool setPtrOnly)
{
#if SP_DEBUG > 1
  std::cout <<"\nResolving MPC chains"<< std::endl;
  for (MPCIter c = allMPCs.begin(); c != allMPCs.end(); c++) std::cout << **c;
#endif

  int nresolved = 0;
  for (MPCIter cit = allMPCs.begin(); cit != allMPCs.end(); cit++)
    if (setPtrOnly)
      for (size_t i = 0; i < (*cit)->getNoMaster(); i++)
      {
        MPC master((*cit)->getMaster(i).node,(*cit)->getMaster(i).dof);
        MPCIter chain = allMPCs.find(&master);
        if (chain != allMPCs.end())
          (*cit)->updateMaster(i,*chain); // Set pointer to next MPC in chain
      }
    else if (ASMbase::resolveMPCchain(allMPCs,*cit))
      nresolved++;

  if (nresolved > 0)
    IFEM::cout <<"Resolved "<< nresolved <<" MPC chains."<< std::endl;
}


/*!
  Recursive resolving of (possibly multi-level) chaining in multi-point
  constraint equations (MPCs). If a master dof in one MPC is specified as a
  slave by another MPC, it is replaced by the master(s) of that other equation.
*/

bool ASMbase::resolveMPCchain (const MPCSet& allMPCs, MPC* mpc)
{
  if (!mpc) return false;

  bool resolved = false;
  for (size_t i = 0; i < mpc->getNoMaster();)
  {
    MPC master(mpc->getMaster(i).node,mpc->getMaster(i).dof);
    MPCIter cit = allMPCs.find(&master);
    if (cit != allMPCs.end())
    {
      // We have a master dof which is a slave in another constraint equation.
      // Invoke resolveMPCchain recursively to ensure that all master dofs
      // of that equation are not slaves themselves.
      ASMbase::resolveMPCchain(allMPCs,*cit);

      // Remove current master specification
      double coeff = mpc->getMaster(i).coeff;
      mpc->removeMaster(i);

      // Add constant offset from the other equation
      mpc->addOffset(coeff*(*cit)->getSlave().coeff);

      // Add masters from the other equations
      for (size_t j = 0; j < (*cit)->getNoMaster(); j++)
        mpc->addMaster((*cit)->getMaster(j).node,
                       (*cit)->getMaster(j).dof,
                       (*cit)->getMaster(j).coeff*coeff);
      resolved = true;
    }
    else
      i++;
  }

#if SP_DEBUG > 1
  if (resolved) std::cout <<"Resolved constraint: "<< *mpc;
#endif
  return resolved;
}


void ASMbase::setNodeNumbers (const std::vector<int>& nodes)
{
  for (size_t i = 0; i < nodes.size() && i < myMLGN.size(); i++)
    myMLGN[i] = 1 + nodes[i]; // Make node numbers 1-based
}


bool ASMbase::hasTimeDependentDirichlet (const std::map<int,RealFunc*>& func,
                                         const std::map<int,VecFunc*>& vfunc)
{
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;
  for (MPCMap::iterator cit = dCode.begin(); cit != dCode.end(); cit++)
    if ((fit = func.find(cit->second)) != func.end())
    {
      if (!fit->second->isConstant())
        return true;
    }
    else if ((vfit = vfunc.find(cit->second)) != vfunc.end())
    {
      if (!vfit->second->isConstant())
        return true;
    }

  return false;
}


bool ASMbase::updateDirichlet (const std::map<int,RealFunc*>& func,
                               const std::map<int,VecFunc*>& vfunc, double time,
                               const std::map<int,int>* g2l)
{
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;
  for (MPCMap::iterator cit = dCode.begin(); cit != dCode.end(); cit++)
  {
    size_t inod = this->getNodeIndex(cit->first->getSlave().node);
    if (inod < 1)
    {
      std::cerr <<" *** ASMbase::updateDirichlet: Invalid slave node in MPC, "
                << *cit->first;
      return false;
    }

    int node = cit->first->getSlave().node;
    if (g2l) node = utl::findKey(*g2l,node);
    Vec4 X(this->getCoord(inod),time,node);
    if ((fit = func.find(cit->second)) != func.end())
    {
      RealFunc& g = *fit->second;
      if (g.isZero())
        cit->first->setSlaveCoeff(0.0);
      else
        cit->first->setSlaveCoeff(g(X));
    }
    else if ((vfit = vfunc.find(cit->second)) != vfunc.end())
    {
      int idof = cit->first->getSlave().dof;
      VecFunc& g = *vfit->second;
      if (g.isZero())
        cit->first->setSlaveCoeff(0.0);
      else
        cit->first->setSlaveCoeff(g(X)[idof-1]);
    }
    else
    {
      std::cerr <<" *** ASMbase::updateDirichlet: Code "<< cit->second
                <<" is not associated with any function."<< std::endl;
      return false;
    }
  }

  return true;
}


void ASMbase::addNeighbor (ASMbase* pch)
{
  if (std::find(neighbors.begin(),neighbors.end(),pch) != neighbors.end())
    return;

  neighbors.push_back(pch);
  pch->neighbors.push_back(this);
}


bool ASMbase::collapseNodes (ASMbase& pch1, int node1, ASMbase& pch2, int node2)
{
  if (node1 < 1 || (size_t)node1 > pch1.myMLGN.size()) return false;
  if (node2 < 1 || (size_t)node2 > pch2.myMLGN.size()) return false;

  if (pch1.myMLGN[node1-1] > pch2.myMLGN[node2-1])
    pch1.mergeNodes(node1,pch2.myMLGN[node2-1],true);
  else if (pch2.myMLGN[node2-1] > pch1.myMLGN[node1-1])
    pch2.mergeNodes(node2,pch1.myMLGN[node1-1],true);
  else
    return false;

#if SP_DEBUG > 1
  std::cout <<"Node "<< node1 <<" in P"<< pch1.idx+1
            <<" and Node "<< node2 <<" in P"<< pch2.idx+1
            <<" assigned common global node number "
            << pch1.MLGN[node1-1] << std::endl;
#endif

  return true;
}


bool ASMbase::mergeNodes (size_t inod, int globalNum, bool silence)
{
  if (inod < 1 || inod > myMLGN.size())
    return false;

  int oldNum = myMLGN[inod-1];
  if (oldNum == globalNum)
    return false;

  if (!silence)
    IFEM::cout <<"  ** Merging duplicated nodes "<< globalNum <<" and "<< oldNum
               <<" at X="<< this->getCoord(inod) << std::endl;

  std::map<int,int> old2New;
  myMLGN[inod-1] = old2New[oldNum] = globalNum;
  for (ASMVec::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    (*it)->renumberNodes(old2New,true);

  return this->renumberNodes(old2New);
}


int ASMbase::renumberNodes (const ASMVec& model, std::map<int,int>& old2new)
{
  ASMVec::const_iterator it;
  std::map<int,int>::iterator nit;

  for (it = model.begin(); it != model.end(); it++)
    if (!(*it)->shareFE)
      for (size_t i = 0; i < (*it)->myMLGN.size(); i++)
        old2new[(*it)->myMLGN[i]] = (*it)->myMLGN[i];

  int n, renum = 0;
  for (n = 1, nit = old2new.begin(); nit != old2new.end(); nit++, n++)
    if (nit->second > n)
    {
      nit->second = n;
      renum++;
    }

  if (renum > 0)
    for (it = model.begin(); it != model.end(); it++)
      for (size_t i = 0; i < (*it)->myMLGN.size(); i++)
	utl::renumber((*it)->myMLGN[i],old2new);

  return renum;
}


int ASMbase::renumberNodes (std::map<int,int>& old2new, int& nNod)
{
  if (shareFE) return 0;

  int renum = 0;
  for (size_t j = 0; j < myMLGN.size(); j++)
    if (utl::renumber(myMLGN[j],nNod,old2new))
      renum++;

  if (renum == 0)
    nNod = std::max(nNod,*std::max_element(MLGN.begin(),MLGN.end()));

  return renum;
}


bool ASMbase::renumberNodes (const std::map<int,int>& old2new, bool renumNodes)
{
#ifdef SP_DEBUG
  bool printInvalidNodes = old2new.size() > 1;
#else
  bool printInvalidNodes = false;
#endif

  int invalid = 0;
  if (renumNodes)
    for (size_t j = 0; j < myMLGN.size(); j++)
      if (!utl::renumber(myMLGN[j],old2new,printInvalidNodes))
	if (old2new.size() > 1)
	  invalid++;

  for (BCVec::iterator bit = BCode.begin(); bit != BCode.end();)
    if (utl::renumber(bit->node,old2new,printInvalidNodes))
      bit++;
    else if (old2new.size() > 1)
    {
      bit = BCode.erase(bit);
      invalid++;
    }
    else
      bit++;

  for (MPCIter mit = mpcs.begin(); mit != mpcs.end(); mit++)
    invalid += (*mit)->renumberNodes(old2new,printInvalidNodes);

  if (invalid == 0 || old2new.size() == 1) return true;

  std::cerr <<" *** "<< invalid <<" invalid nodes found while renumbering\n";
  return false;
}


void ASMbase::extractElmRes (const Matrix& globRes, Matrix& elmRes) const
{
  elmRes.resize(globRes.rows(),MLGE.size(),true);

  size_t iel, ivel = 0;
  for (iel = 0; iel < MLGE.size(); iel++)
    if (MLGE[iel] > 0)
      elmRes.fillColumn(++ivel,globRes.getColumn(MLGE[iel]));

  elmRes.resize(globRes.rows(),ivel);
}


void ASMbase::extractNodeVec (const Vector& globRes, Vector& nodeVec,
			      const int* madof) const
{
  nodeVec.clear();
  nodeVec.reserve(nf*MLGN.size());
  for (size_t i = 0; i < MLGN.size(); i++)
  {
    int inod = MLGN[i];
    int idof = madof[inod-1] - 1;
    int jdof = madof[inod] - 1;
#ifdef INDEX_CHECK
    if (inod < 1 || jdof > (int)globRes.size())
      std::cerr <<" *** ASMbase::extractNodeVec: Global DOF "<< jdof
                <<" is out of range [1,"<< globRes.size() <<"]."<< std::endl;
#endif
    nodeVec.insert(nodeVec.end(),globRes.ptr()+idof,globRes.ptr()+jdof);
  }
}

void ASMbase::injectNodeVec (const std::vector<int>& madof,
                             const Vector& nodeVec, Vector& globVec,
                             int basis) const
{
  size_t ldof = 0;
  char bType = basis == 1 ? 'D' : 'P'+basis-2;
  for (size_t i = 0; i < MLGN.size(); i++)
  {
    if (basis == 0 || getNodeType(i+1) == bType) {
      int inod = MLGN[i];
      int idof = madof[inod-1] - 1;
      int jdof = madof[inod] - 1;
#ifdef INDEX_CHECK
      if (inod < 1 || jdof > (int)globVec.size())
        std::cerr <<" *** ASMbase::injectNodeVec: Global DOF "<< jdof
                  <<" is out of range [1,"<< globVec.size() <<"]"<< std::endl;
#endif
      std::copy(nodeVec.begin()+ldof, nodeVec.begin()+ldof+(jdof-idof), globVec.begin()+idof);
      ldof += jdof-idof;
    }
  }
}


void ASMbase::extractNodeVec (const Vector& globRes, Vector& nodeVec,
			      unsigned char nndof, int basis) const
{
  if (nndof == 0) nndof = nf;

  // Don't extract the Lagrange multipliers, if any
  size_t nNod = myLMs.first > 0 ? myLMs.first-1 : MLGN.size();

  nodeVec.resize(nndof*nNod);
  double* nodeP = nodeVec.ptr();
  for (size_t i = 1; i <= nNod; i++, nodeP += nndof)
  {
    // Note: If nndof==1 (scalar field) we should ignore possibly added nodes
    // due to Dirichlet conditions defined in local axes because these nodes are
    // assumed not to have entries in the globRes vector. In that case, we use
    // the "real" node at that location instead (see, e.g., ASMs2D::getNodeID).
    // The same is assumed if basis < 0 on input (for vector fields) HACK!
    int n = this->getNodeID(i, nndof == 1 || basis < 0) - 1;
#ifdef INDEX_CHECK
    if (n < 0 || nndof*(size_t)(n+1) > globRes.size())
      std::cerr <<" *** ASMbase::extractNodeVec: Global DOF "<< nndof*(n+1)
                <<" is out of range [1,"<< globRes.size() <<"]."<< std::endl;
#endif
    memcpy(nodeP,globRes.ptr()+nndof*n,nndof*sizeof(double));
  }
}


bool ASMbase::injectNodeVec (const Vector& nodeVec, Vector& globRes,
			     unsigned char nndof, int) const
{
  if (nndof == 0) nndof = nf;

  // Don't inject the Lagrange multipliers, if any
  size_t nNod = myLMs.first > 0 ? myLMs.first-1 : MLGN.size();

  if (nodeVec.size() != nNod*nndof)
  {
    std::cerr <<" *** ASMbase::injectNodeVec:: Invalid patch vector, size = "
              << nodeVec.size() <<" != "<< nNod*nndof << std::endl;
    return false;
  }

  const double* nodeP = nodeVec.ptr();
  for (size_t i = 0; i < nNod; i++, nodeP += nndof)
  {
    int n = MLGN[i];
    if (n > 0 && nndof*(size_t)n <= globRes.size())
      memcpy(globRes.ptr()+nndof*(n-1),nodeP,nndof*sizeof(double));
#ifdef SP_DEBUG
    else // This is most likely OK, print message only in debug mode
      std::cerr <<" *** ASMbase::injectNodeVec: Global DOF "<< nndof*n
                <<" is out of range [1,"<< globRes.size() <<"]."<< std::endl;
#endif
  }

  return true;
}


bool ASMbase::getSolution (Matrix& sField, const Vector& locSol,
			   const IntVec& nodes) const
{
  sField.resize(nf,nodes.size());
  for (size_t i = 0; i < nodes.size(); i++)
    if (nodes[i] < 1 || (size_t)nodes[i] > MLGN.size())
    {
      std::cerr <<" *** ASMbase::getSolution: Node #"<< nodes[i]
		<<" is out of range [1,"<< MLGN.size() <<"]."<< std::endl;
      return false;
    }
    else if (this->isLMn(nodes[i]))
    {
      std::cerr <<"  ** ASMbase::getSolution: Node #"<< nodes[i]
		<<" is a Lagrange multiplier, returning 0.0."<< std::endl;
      sField.fillColumn(i+1,RealArray(nf,0.0));
    }
    else
      sField.fillColumn(i+1,locSol.ptr()+nf*nodes[i]-nf);

  return true;
}


int ASMbase::searchCtrlPt (RealArray::const_iterator cit,
                           RealArray::const_iterator end,
                           const Vec3& X, int dimension, double tol) const
{
  Vec3 Xnod;
  double distance = 0.0;
  size_t inod = 0, iclose = 0;
  for (int i = 0; cit != end; ++cit, i++)
  {
    if (i < nsd) Xnod[i] = *cit;
    if (i+1 == dimension)
    {
      i = -1;
      inod++;
      if (X.equal(Xnod,tol))
        if (iclose == 0 || (X-Xnod).length() < distance)
        {
          // In case of a very fine grid where several control points might
          // be withing the tolerance, choose the closest point
          iclose = inod;
          distance = (X-Xnod).length();
        }
    }
  }

  return iclose;
}


bool ASMbase::evalSolution (Matrix&, const Vector&, const int*) const
{
  return Aerror("evalSolution(Matrix&,const Vector&,const int*)");
}


bool ASMbase::evalSolution (Matrix&, const Vector&,
			    const RealArray*, bool, int) const
{
  return Aerror("evalSolution(Matrix&,const Vector&,const RealArray*,bool,int)");
}


bool ASMbase::evalSolution (Matrix&, const IntegrandBase&,
			    const int*, char) const
{
  return Aerror("evalSolution(Matrix&,const IntegrandBase&,const int*,char)");
}


bool ASMbase::evalSolution (Matrix&, const IntegrandBase&,
			    const RealArray*, bool) const
{
  return Aerror("evalSolution(Matrix&,const IntegrandBase&,const RealArray*)");
}


bool ASMbase::globalL2projection (Matrix&, const IntegrandBase&, bool) const
{
  return Aerror("globalL2projection(Matrix&,const IntegrandBase&,bool)");
}


bool ASMbase::evaluate (const ASMbase*, const Vector&, RealArray&, int) const
{
  return Aerror("evaluate(const ASMbase*,const Vector&,RealArray&,int)");
}


bool ASMbase::evaluate (const Field*, RealArray&, int) const
{
  return Aerror("evaluate(const Field*,RealArray&,int)");
}


bool ASMbase::evaluate (const RealFunc*, RealArray&, int, double) const
{
  return Aerror("evaluate(const RealFunc*,RealArray&,int,double)");
}
