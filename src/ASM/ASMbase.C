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
#include "Tensor.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "Function.h"
#include "Utilities.h"
#include <algorithm>
#include <functional>
#include <iomanip>


bool ASMbase::fixHomogeneousDirichlet = true;
bool ASMbase::refineGeometry = false;
int  ASMbase::dbgElm = 0;
ASM::CachePolicy ASMbase::cachePolicy = ASM::PRE_CACHE;

//! This quantitiy is used to scale the characteristic element sizes which
//! are used by residual error estimates, etc., such that they always are in
//! the range [0,1.0]. The applications have to set an appropriate value,
//! when needed.
double ASMbase::modelSize = 1.0;

int ASMbase::gEl = 0;
int ASMbase::gNod = 0;
std::map<int,int> ASMbase::xNode;


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
}


ASMbase::ASMbase (const ASMbase& patch, unsigned char n_f)
  : MLGE(patch.MLGE), MLGN(patch.MLGN), MNPC(patch.MNPC), shareFE('F'),
    firstBp(patch.firstBp), myLMTypes(patch.myLMTypes), myLMs(patch.myLMs),
    myRmaster(patch.myRmaster)
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
  // Note: Properties are _not_ copied
}


ASMbase::ASMbase (const ASMbase& patch)
  : MLGE(myMLGE), MLGN(myMLGN), MNPC(myMNPC), shareFE('S'),
    BCode(patch.BCode), firstBp(patch.firstBp)
{
  nf = patch.nf;
  nsd = patch.nsd;
  ndim = patch.ndim;
  nGauss = patch.nGauss;
  nel = patch.nel;
  nnod = patch.nnod;
  idx = patch.idx;
  firstIp = patch.firstIp;

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
}


ASMbase::~ASMbase ()
{
  for (MPC* mpc : mpcs)
    delete mpc;
}


ASMbase* ASMbase::cloneUnShared () const
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
    if (shareFE == 'F')
    {
      const_cast<IntVec&>(MLGE).clear();
      const_cast<IntMat&>(MNPC).clear();
    }
  }
  else // Don't erase the elements, but set them to have zero nodes
    for (IntVec& mnpc : myMNPC)
      mnpc.clear();

  // Erase the nodes, boundary conditions and multi-point constraints
  for (MPC* mpc : mpcs)
    delete mpc;

  myLMs.clear();
  myLMTypes.clear();
  myRmaster.clear();

  myMLGN.clear();
  if (shareFE == 'F')
    const_cast<IntVec&>(MLGN).clear();

  BCode.clear();
  dCode.clear();
  mpcs.clear();
}


bool ASMbase::addXElms (short int, short int, size_t, IntVec&)
{
  return Aerror("addXElms(short int,short int,size_t,IntVec&)");
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

  for (int iLag : mGLag)
  {
    size_t node = 1 + utl::findIndex(MLGN,iLag);
    if (node == 0)
    {
      // Add a new Lagrange multiplier node
      myMLGN.push_back(iLag);
      node = myMLGN.size();
    }

    // Update the nodal (1-based) indices of the Lagrange multipliers
    if (myLMs.empty() || node >= *myLMs.begin())
      myLMs.insert(node);
    else
    {
      std::cerr <<" *** ASMbase::addLagrangeMultipliers: Node "<< node
                <<" is out of range ["<< *myLMs.begin() <<","<< *myLMs.rbegin()
                <<"]."<< std::endl;
      return false;
    }

    size_t idxLag = node - *myLMs.begin();
    if (myLMTypes.size() < idxLag+1)
      myLMTypes.resize(idxLag+1,0);

    myLMTypes[idxLag] = iel == 0 ? 'G' : 'L';

    // Extend the element connectivity table
    if (iel > 0)
      myMNPC[iel-1].push_back(node-1);
    else for (IntVec& mnpc : myMNPC)
      mnpc.push_back(node-1);
  }

  return true;
}


bool ASMbase::addGlobalLagrangeMultipliers (const IntVec& mGLag,
                                            unsigned char nnLag)
{
  return this->addLagrangeMultipliers(0,mGLag,nnLag);
}


void ASMbase::resetNumbering (int n)
{
  gEl = 0;
  gNod = n;
  xNode.clear();
}


size_t ASMbase::getNodeIndex (int globalNum, bool) const
{
  return 1 + utl::findIndex(MLGN,globalNum);
}


int ASMbase::getNodeID (size_t inod, bool) const
{
  return inod < 1 || inod > MLGN.size() ? 0 : MLGN[inod-1];
}


char ASMbase::getLMType (size_t inod) const
{
  std::set<size_t>::const_iterator firstLM = myLMs.begin();
  return this->isLMn(inod) ? myLMTypes[inod-(*firstLM)] : 0;
}


size_t ASMbase::getElmIndex (int globalNum) const
{
  return 1 + utl::findIndex(MLGE,globalNum);
}


int ASMbase::getElmID (size_t iel) const
{
  return iel < 1 || iel > MLGE.size() ? 0 : abs(MLGE[iel-1]);
}


const IntVec& ASMbase::getElementNodes (int iel) const
{
  if (iel > 0 && (size_t)iel <= MNPC.size())
    return MNPC[iel-1];

  static IntVec empty;
  return empty;
}


unsigned char ASMbase::getNodalDOFs (size_t inod) const
{
  if (this->isLMn(inod))
    return nLag;
  else if (this->isRMn(inod))
    return nsd < 3 ? 3 : 6; // Including rotational DOFs
  else
    return nf;
}


char ASMbase::getNodeType (size_t inod) const
{
  return this->isLMn(inod) ? this->getLMType(inod) : (inod > nnod ? 'X' : 'D');
}


size_t ASMbase::getNoNodes (int basis) const
{
  if (basis > 0)
    return nnod;
  else if (basis < 0 && !myLMs.empty())
    return *myLMs.begin() - 1;
  else
    return MLGN.size();
}


size_t ASMbase::getNoElms (bool includeZeroVolElms, bool includeXElms) const
{
  if (includeZeroVolElms)
    return nel;

  size_t numels = 0;
  for (int iel : MLGE)
    if (iel > 0 || (includeXElms && iel < 0))
      numels++;

  return numels;
}


void ASMbase::getNoIntPoints (size_t& nPt, size_t& nIPt)
{
  size_t nGp = 1;
  if (nGauss > 0 && nGauss <= 10)
    for (unsigned char d = 0; d < ndim; d++)
      nGp *= nGauss;
  else
  {
    // Use polynomial order to define number of quadrature points
    int ng[3] = { 0, 0, 0 };
    this->getOrder(ng[0],ng[1],ng[2]);
    for (unsigned char d = 0; d < ndim && d < 3; d++)
      nGp *= ng[d] + nGauss%10;
  }

  firstIp = nPt;
  nPt += nel*nGp; // Note: Includes also the 0-span elements

  // Count additional interface quadrature points
  size_t nInterface = MLGE.size() - nel;
  if (nInterface > 0 && nInterface != nel && nGauss > 0 && nGauss <= 10)
    nIPt += nInterface*nGp/nGauss;
}


void ASMbase::getNoBouPoints (size_t& nPt, char ldim, char lindx)
{
  if (ldim+1 == ndim)
    lindx %= 10; // Mask off Neumann order flag

  size_t nGp = 1;
  if (nGauss > 0 && nGauss <= 10)
    for (char d = 0; d < ldim; d++)
      nGp *= nGauss;
  else
  {
    // Use polynomial order to define number of quadrature points
    int ng[3] = { 0, 0, 0 };
    this->getOrder(ng[0],ng[1],ng[2]);
    ng[(lindx-1)/2] = 1;
    for (unsigned char d = 0; d < ndim; d++)
      nGp *= ng[d];
  }

  firstBp[lindx] = nPt;

  nPt += this->getNoBoundaryElms(lindx,ldim)*nGp; // Includes 0-span elements
}


void ASMbase::printNodes (std::ostream& os) const
{
  Matrix X;
  this->getNodalCoordinates(X);
  os <<"\n\nNodal coordinates for Patch "<< idx+1;
  for (size_t inod = 1; inod <= X.cols(); inod++)
  {
    os <<'\n'<< std::setw(4) << inod
       << ' '<< std::setw(4) << MLGN[inod-1] <<':';
    for (size_t i = 1; i <= X.rows(); i++)
      os <<' '<< X(i,inod);
  }
  os << std::endl;
}


void ASMbase::printElements (std::ostream& os) const
{
  if (MNPC.empty()) return;

  os <<"\n\nElement connectivities for Patch "<< idx+1;
  for (size_t iel = 0; iel < MLGE.size(); iel++)
  {
    os <<'\n'<< std::setw(4) << iel+1
       << ' '<< std::setw(4) << MLGE[iel] <<':';
    if (iel < MNPC.size())
      for (int node : MNPC[iel]) os <<" "<< node;
  }
  for (size_t ielx = MLGE.size(); ielx < MNPC.size(); ielx++)
  {
    os <<'\n'<< std::setw(4) << ielx+1 <<" ----:";
    for (int node : MNPC[ielx]) os <<" "<< node;
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


bool ASMbase::addMPC (MPC*& mpc, int code, bool verbose)
{
  if (!mpc) return true;

  if (this->isFixed(mpc->getSlave().node,mpc->getSlave().dof))
  {
    // Silently ignore MPC's on dofs that already are marked as FIXED
    delete mpc;
    mpc = nullptr;
    return true;
  }

  std::pair<MPCIter,bool> mit = mpcs.insert(mpc);
  if (mit.second)
  {
#if SP_DEBUG > 1
    if (verbose) std::cout <<"Added constraint: "<< *mpc;
#endif
    if (code > 0) dCode[mpc] = code;
    return true;
  }

#ifdef SP_DEBUG
  if (verbose) std::cout <<"Ignored constraint (duplicated slave): "<< *mpc;
#endif
  delete mpc;
  mpc = *mit.first; // This dof is already a slave in another MPC
  return false;
}


bool ASMbase::add2PC (int slave, int dir, int master, int code)
{
  if (dir < 1 || dir > nf) return true;
  if (slave == master) return true;

  MPC* cons = new MPC(slave,dir);
  bool stat = this->addMPC(cons,code);
  if (!cons) return stat;

  cons->addMaster(master,dir);
#if SP_DEBUG > 1
  std::cout <<"Added constraint: "<< *cons;
#endif
  return stat;
}


bool ASMbase::add3PC (int slave, int dir, int master1, int master2, int code)
{
  if (master1 == master2)
    return this->add2PC(slave,dir,master1,code);

  if (dir < 1 || dir > nf) return true;

  MPC* cons = new MPC(slave,dir);
  bool stat = this->addMPC(cons,code);
  if (!cons) return stat;

  if (master1 != slave) cons->addMaster(master1,dir);
  if (master2 != slave) cons->addMaster(master2,dir);
#if SP_DEBUG > 1
  std::cout <<"Added constraint: "<< *cons;
#endif
  return stat;
}


void ASMbase::addLocal2GlobalCpl (int iSlave, int master, const Tensor& Tlg)
{
  // Establish constraint equations relating the global DOFs
  // of the slave node to the local DOFs of the master node.
  // We here assume there are (at least) nsd unknowns per node,
  // and that only the first nsd DOFs are subjected to transformation.
  int fixDirs = 0;
  for (unsigned char d = 1; d <= nf; d++)
  {
    MPC* cons = new MPC(MLGN[iSlave],d);
    if (this->addMPC(cons) && cons)
    {
      if (d > nsd)
      {
        if (!this->isFixed(master,d))
          cons->addMaster(master,d);
      }
      else for (unsigned char c = 1; c <= nsd; c++)
      {
        if (!this->isFixed(master,c))
          cons->addMaster(master,c,Tlg(d,c));
      }

      if (cons->getNoMaster() == 0)
      {
        // All master DOFs are fixed.
        // Then the MPC is not needed, fix the slave DOF instead.
        mpcs.erase(cons);
        delete cons;
        fixDirs = d + 10*fixDirs;
      }
#if SP_DEBUG > 1
      else
        std::cout <<"Added constraint: "<< *cons;
#endif
    }
  }

  if (fixDirs > 0)
    this->fix(iSlave+1,fixDirs);
}


bool ASMbase::createRgdMasterNode (int& gMaster, const Vec3& Xpt)
{
  bool newNode = gMaster == 0;
  if (newNode)
    gMaster = ++gNod;
  else if (std::find(MLGN.begin()+nnod,MLGN.end(),gMaster) != MLGN.end())
    return newNode; // This node has already been created

#if SP_DEBUG > 1
  std::cout <<"Adding extra-ordinary node "<< gMaster
            <<" for rigid coupling in Patch "<< idx+1 << std::endl;
#endif
  myMLGN.push_back(gMaster);
  myRmaster[myMLGN.size()] = { Xpt.x, Xpt.y, Xpt.z };
  return newNode;
}


void ASMbase::addRigidMPC (int gSlave, int gMaster, const Vec3& dX)
{
  for (unsigned short int dof = 1; dof <= nf; dof++)
  {
    MPC* cons = new MPC(gSlave,dof);
    if (this->addMPC(cons) && cons)
    {
      // Add one-to-one translation coupling
      cons->addMaster(gMaster,dof,1.0);
      // Add translation-to-rotation coupling
      switch (dof) {
      case 1: // u_x = dZ*theta_y - dY*theta_z
        cons->addMaster(gMaster,5, dX.z);
        cons->addMaster(gMaster,6,-dX.y);
        break;
      case 2: // u_y = dX*theta_z - dZ*theta_x
        cons->addMaster(gMaster,4,-dX.z);
        cons->addMaster(gMaster,6, dX.x);
        break;
      case 3: // u_z = dY*theta_x - dX*theta_y
        cons->addMaster(gMaster,4, dX.y);
        cons->addMaster(gMaster,5,-dX.x);
        break;
      }
#if SP_DEBUG > 1
      std::cout <<"Added constraint: "<< *cons;
#endif
    }
  }
}


bool ASMbase::addRigidCpl (int lindx, int ldim, int basis,
                           int& gMaster, const Vec3& Xmaster, bool extraPt)
{
  if (ldim+1 != ndim)
  {
    IFEM::cout <<"  ** ASMbase::addRigidCpl: Not implemented for "
               << ldim <<"-dimensional boundaries (ignored)."<< std::endl;
    return false;
  }

  if (extraPt) // The master point is not a patch node, create an extra node
    extraPt = this->createRgdMasterNode(gMaster,Xmaster);

  IntVec nodes;
  this->getBoundaryNodes(lindx,nodes,basis,1,0,true);
  for (int node : nodes)
  {
    Vec3 dX = this->getCoord(node) - Xmaster;
    if (nsd == 3)
      this->addRigidMPC(MLGN[node-1],gMaster,dX);
    else if (nsd == 2 && nf > 1)
    {
      // Special for 2D problems with (at least) 2 nodal DOFs
      MPC* cons = new MPC(MLGN[node-1],1);
      if (this->addMPC(cons) && cons)
      {
        cons->addMaster(gMaster,1, 1.0);
        cons->addMaster(gMaster,3,-dX.y);
      }
#if SP_DEBUG > 1
      std::cout <<"Added constraint: "<< *cons;
#endif
      cons = new MPC(MLGN[node-1],2);
      if (this->addMPC(cons) && cons)
      {
        cons->addMaster(gMaster,2, 1.0);
        cons->addMaster(gMaster,3, dX.x);
#if SP_DEBUG > 1
        std::cout <<"Added constraint: "<< *cons;
#endif
      }
    }
  }

  return extraPt;
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

  if (this->add2PC(masterNode,dir,slaveNode) ||
      this->add2PC(slaveNode,dir,masterNode))
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
  else for (int dof : dofs)
    this->addPeriodicity(master,slave,dof);
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


void ASMbase::constrainNodes (const IntVec& nodes, int dof, int code)
{
  if (code < 0) code = -code;

  for (int node : nodes)
    if (node > 0 && node <= (int)this->getNoNodes(1))
      this->prescribe(node,dof,code);
    else
      std::cerr <<"  ** ASMbase::constrainNodes: Node "<< node
                <<" is out of range [1,"<< this->getNoNodes(1)
                <<"]."<< std::endl;
}


bool ASMbase::constrainXnode (int node, int dof, int code)
{
  for (size_t inod = nnod+1; inod <= MLGN.size(); inod++)
    if (this->isRMn(inod) && MLGN[inod-1] == node)
    {
      this->prescribe(inod,dof,code);
      return true;
    }

  return false;
}


int ASMbase::prescribe (size_t inod, int dirs, int code)
{
  if (code == 0 && fixHomogeneousDirichlet)
    return this->fix(inod,dirs);

  int node = this->getNodeID(inod);
  if (node < 1 || dirs < 1) return dirs;

  int ignoredDirs = 0;
  for (int dof : utl::getDigits(dirs))
    if (dof <= nf)
    {
      MPC* mpc = new MPC(node,dof,1.0);
      if (!this->addMPC(mpc,code,true))
        ignoredDirs = 10*ignoredDirs + dof;
    }
    else
    {
      ignoredDirs = 10*ignoredDirs + dof;
      std::cerr <<"  ** ASMbase::prescribe: Ignoring invalid DOF code "<< dof
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
  for (int dof : utl::getDigits(dirs))
    if (dof <= nf)
      switch (dof) {
      case 1: bit->CX = 0; break;
      case 2: bit->CY = 0; break;
      case 3: bit->CZ = 0; break;
      case 4: bit->RX = 0; break;
      case 5: bit->RY = 0; break;
      case 6: bit->RZ = 0; break;
      }
    else
    {
      invalidDOFs = 10*invalidDOFs + dof;
      std::cerr <<"  ** ASMbase::fix: Ignoring invalid DOF code "<< dof
                << std::endl;
    }

#if SP_DEBUG > 1
  if (!(old == *bit))
    std::cout <<"\tFixed node: "<< node <<" "<< dirs << std::endl;
#endif
  return invalidDOFs;
}


bool ASMbase::allDofs (int dirs) const
{
  int myDof = 0;
  for (int dof : utl::getDigits(dirs))
    if (++myDof != dof)
      return false;

  return myDof == nf;
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
  for (ASMbase* pch : model)
  {
    std::pair<MPCIter,bool> ret;
    std::vector<MPC*> uniqueMPC;
    uniqueMPC.reserve(pch->getNoMPCs());
    for (MPC* mpc : pch->mpcs)
      if ((ret = allMPCs.insert(mpc)).second)
	uniqueMPC.push_back(mpc);
      else
      {
	// Merge multiple constraint equations with common slave definition
	if ((*ret.first)->getSlave().coeff == 0.0 &&
	    (*ret.first)->getNoMaster() > 1 &&
	    (*ret.first)->merge(mpc))
	{
#if SP_DEBUG > 1
	  std::cout <<"Merging constraint "<< *mpc;
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
	  std::cout <<"Deleted constraint "<< *mpc;
#endif
	  ndeleted++;
	}
	pch->dCode.erase(mpc);
	delete mpc;
      }

    if (uniqueMPC.size() < pch->getNoMPCs())
    {
      // Compress the MPC set for this patch removing duplicated entries
      pch->mpcs.clear();
      pch->mpcs.insert(uniqueMPC.begin(),uniqueMPC.end());
    }
  }

  if (nmerged > 0)
    IFEM::cout <<"Merged "<< nmerged <<" MPC equations."<< std::endl;
  if (ndeleted > 0)
    IFEM::cout <<"Deleted "<< ndeleted <<" MPC equations."<< std::endl;

#if SP_DEBUG > 1
  if (allMPCs.empty()) return;
  std::cout <<"\nMulti-point constraints:\n";
  for (const MPC* mpc : allMPCs) std::cout << *mpc;
#endif
}


/*!
  Recursive resolving of (possibly multi-level) chaining in the multi-point
  constraint equations (MPCs). If a master dof in one MPC is specified as a
  slave by another MPC, it is replaced by the master(s) of that other equation.
  Since an MPC-equation may couple nodes belonging to different patches,
  this method must have access to all patches in the model.

  If \a setPtrOnly is \e true, the MPC equations are not modified. Instead
  the pointers to the next MPC in the chain is assigned for the master DOFs
  which are slaves in other MPCs.

  This is used in time-dependent/nonlinear simulators where the constraint
  coefficients might be updated due to time variation. The resolving is then
  done directly in the MMCEQ/TTCC arrays of the SAM data object.
  \sa SAMpatch::updateConstraintEqs.
*/

void ASMbase::resolveMPCchains (const MPCSet& allMPCs,
                                const ASMVec& model, bool setPtrOnly)
{
#if SP_DEBUG > 1
  int count = 0;
  std::cout <<"\nResolving MPC chains"<< std::endl;
  for (MPC* mpc : allMPCs) std::cout << std::setw(4) << ++count <<": "<< *mpc;
  if (setPtrOnly) std::cout <<"\nResolved MPCs"<< std::endl;
  count = 0;
#endif

  // Recursive lambda function for resolving the chained MPC-equations.
  std::function<bool(MPC*)> resolve = [&allMPCs,&resolve](MPC* mpc)
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
        // Invoke resolve() recursively to ensure that all master dofs of that
        // equation are not slaves themselves.
        resolve(*cit);

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
  };

  // Recursive lambda function for linking the chained MPC-equations.
  std::function<bool(MPC*)> resolveLnk = [&allMPCs,&model,&resolveLnk](MPC* mpc)
  {
    if (!mpc) return false;

    for (size_t i = 0; i < mpc->getNoMaster();)
    {
      int node = mpc->getMaster(i).node;
      int ldof = mpc->getMaster(i).dof;
      for (ASMbase* pch : model)
        if (pch->isFixed(node,ldof))
        {
          ldof = -1;
          break;
        }

      if (ldof > 0)
      {
        MPC master(node,ldof);
        MPCIter chain = allMPCs.find(&master);
        if (chain == allMPCs.end())
          i++; // Free master DOF
        else if (!resolveLnk(*chain))
          mpc->removeMaster(i); // This master DOF is fixed, remove from link
        else
          mpc->updateMaster(i++,*chain); // Set pointer to next MPC in chain
      }
      else
        mpc->removeMaster(i);
    }

    return mpc->getNoMaster() > 0 || mpc->getSlave().coeff != 0.0;
  };

  int nresolved = 0;
  for (MPC* mpc : allMPCs)
    if (setPtrOnly)
    {
      resolveLnk(mpc);
#if SP_DEBUG > 1
      count++;
      if (mpc->isChained() || mpc->getNoMaster() == 0)
        std::cout << std::setw(4) << count <<": "<< *mpc;
#endif
    }
    else if (resolve(mpc))
      nresolved++;

  if (nresolved > 0)
    IFEM::cout <<"Resolved "<< nresolved <<" MPC chains."<< std::endl;
}


void ASMbase::setNodeNumbers (const IntVec& nodes)
{
  for (size_t i = 0; i < nodes.size() && i < myMLGN.size(); i++)
    myMLGN[i] = 1 + nodes[i]; // Make node numbers 1-based
}


bool ASMbase::hasTimeDependentDirichlet (const std::map<int,RealFunc*>& func,
                                         const std::map<int,VecFunc*>& vfunc)
{
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;
  for (const MPCMap::value_type& cit : dCode)
    if ((fit = func.find(cit.second)) != func.end())
    {
      if (!fit->second->isConstant())
        return true;
    }
    else if ((vfit = vfunc.find(cit.second)) != vfunc.end())
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
  for (const MPCMap::value_type& cit : dCode)
  {
    size_t inod = this->getNodeIndex(cit.first->getSlave().node);
    if (inod < 1)
    {
      std::cerr <<" *** ASMbase::updateDirichlet: Invalid slave node in MPC, "
                << *cit.first;
      return false;
    }

    int node = cit.first->getSlave().node;
    if (g2l) node = utl::findKey(*g2l,node);
    Vec4 X(this->getCoord(inod),time,node);
    if ((fit = func.find(cit.second)) != func.end())
    {
      RealFunc& g = *fit->second;
      if (g.isZero())
        cit.first->setSlaveCoeff(0.0);
      else
        cit.first->setSlaveCoeff(g(X));
    }
    else if ((vfit = vfunc.find(cit.second)) != vfunc.end())
    {
      int idof = cit.first->getSlave().dof;
      VecFunc& g = *vfit->second;
      if (g.isZero())
        cit.first->setSlaveCoeff(0.0);
      else
        cit.first->setSlaveCoeff(g(X)[idof-1]);
    }
    else
    {
      std::cerr <<" *** ASMbase::updateDirichlet: Code "<< cit.second
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
    pch1.mergeNodes(node1,pch2.myMLGN[node2-1],false);
  else if (pch2.myMLGN[node2-1] > pch1.myMLGN[node1-1])
    pch2.mergeNodes(node2,pch1.myMLGN[node1-1],false);
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


bool ASMbase::mergeNodes (size_t inod, int globalNum, bool verbose)
{
  if (inod < 1 || inod > myMLGN.size())
    return false;

  int oldNum = myMLGN[inod-1];
  if (oldNum == globalNum)
    return false;

  if (verbose)
    IFEM::cout <<"  ** Merging duplicated nodes "<< globalNum <<" and "<< oldNum
               <<" at X="<< this->getCoord(inod) << std::endl;

  std::map<int,int> old2New;
  myMLGN[inod-1] = old2New[oldNum] = globalNum;
  for (ASMbase* pch : neighbors)
    pch->renumberNodes(old2New,true);

  return this->renumberNodes(old2New);
}


int ASMbase::renumberNodes (const ASMVec& model, std::map<int,int>& old2new)
{
  for (const ASMbase* pch : model)
    if (!pch->shareFE)
      for (int n : pch->myMLGN)
        old2new[n] = n;

  int n = 0, renum = 0;
  for (std::pair<const int,int>& node : old2new)
    if (node.second > ++n)
    {
      node.second = n;
      renum++;
    }

  if (renum > 0)
    for (ASMbase* pch : model)
      for (int& node : pch->myMLGN)
        utl::renumber(node,old2new);

  return renum;
}


int ASMbase::renumberNodes (std::map<int,int>& old2new, int& nNod)
{
  int renum = 0;
  if (!shareFE)
    for (int& node : myMLGN)
      if (utl::renumber(node,nNod,old2new))
        renum++;

  if (renum == 0)
    nNod = std::max(nNod,*std::max_element(MLGN.begin(),MLGN.end()));

  return renum;
}


/*!
  This method renumbers the global node numbers referred by the boundary
  condition- and multi-point constraint equation objects in the patch,
  according to the provided mapping \a old2new.

  The nodes themselves (in \a MLGN) are assumed already to be up to date,
  unless \a renumGN is non-zero. If \a renumGN equals one, all node numbers
  are assumed present in the \a old2new mapping and an error is flagged if
  that is not the case.
*/

bool ASMbase::renumberNodes (const std::map<int,int>& old2new, char renumGN)
{
  bool renumAll = old2new.size() > 1 && renumGN < 2;
#ifdef SP_DEBUG
  bool printInvalidNodes = renumAll;
#else
  bool printInvalidNodes = false;
#endif

  int invalid = 0;
  if (renumGN > 0)
    for (int& node : myMLGN)
      if (!utl::renumber(node,old2new,printInvalidNodes) && renumAll)
        ++invalid;

  for (BCVec::iterator bit = BCode.begin(); bit != BCode.end();)
    if (utl::renumber(bit->node,old2new,printInvalidNodes))
      ++bit;
    else
    {
      if (renumAll)
        bit = BCode.erase(bit); // Remove fixation codes on non-existing node
      else
        ++bit;
      ++invalid;
    }

  for (MPC* mpc : mpcs)
    invalid += mpc->renumberNodes(old2new,printInvalidNodes);

  if (invalid == 0 || !renumAll) return true;

  std::cerr <<" *** "<< invalid <<" invalid nodes found while renumbering\n";
  return false;
}


void ASMbase::shiftGlobalNodeNums (int nshift)
{
  for (int& node : myMLGN)
    node += nshift;

  for (BC& bc : BCode)
    bc.node += nshift;

  for (MPC* mpc : mpcs)
    mpc->shiftNodes(nshift);
}


void ASMbase::shiftGlobalElmNums (int eshift)
{
  for (int& iel : myMLGE)
    iel += eshift;
}


void ASMbase::extractElmRes (const Matrix& globRes, Matrix& elmRes) const
{
  elmRes.resize(globRes.rows(),MLGE.size(),true);

  size_t ivel = 0;
  for (int iel : MLGE)
    if (iel > 0)
      elmRes.fillColumn(++ivel,globRes.getColumn(iel));

  elmRes.resize(globRes.rows(),ivel);
}


bool ASMbase::extractNodalVec (const RealArray& globRes, RealArray& nodeVec,
                               const int* madof, int ngnod) const
{
  nodeVec.clear();
  if (ngnod == 0)
  {
    std::cerr <<" *** ASMbase::extractNodalVec: Empty MADOF array."<< std::endl;
    return false;
  }

  size_t nPchNod = ngnod == -2 ? nnod : MLGN.size();
  nodeVec.reserve(nf*nPchNod);
  for (size_t i = 0; i < nPchNod; i++)
  {
    int inod = MLGN[i];
#ifdef INDEX_CHECK
    if (inod < 1 || (inod > ngnod && ngnod > 0))
    {
      std::cerr <<" *** ASMbase::extractNodalVec: Global node "<< inod;
      if (ngnod > 0)
        std::cerr <<" is out of range [1,"<< ngnod <<"]."<< std::endl;
      else
        std::cerr <<" is out of range."<< std::endl;
      return false;
    }
#endif
    int idof = madof[inod-1] - 1;
    int jdof = madof[inod] - 1;
    if (idof == jdof)
      continue; // DOF-less node
#ifdef INDEX_CHECK
    else if (idof < 0 || idof > jdof || jdof > (int)globRes.size())
    {
      std::cerr <<" *** ASMbase::extractNodalVec: Global DOFs "
                << idof+1 <<" "<< jdof
                <<" out of range [1,"<< globRes.size() <<"]."<< std::endl;
      return false;
    }
#endif
    nodeVec.insert(nodeVec.end(),globRes.data()+idof,globRes.data()+jdof);
  }

  return true;
}


bool ASMbase::injectNodalVec (const RealArray& nodeVec, RealArray& globVec,
                              const IntVec& madof, int basis) const
{
  if (madof.empty())
  {
    std::cerr <<" *** ASMbase::injectNodalVec: Empty MADOF array."<< std::endl;
    return false;
  }

  size_t i = 0, ldof = 0;
  char bType = basis == 1 ? 'D' : 'P'+basis-2;
  for (int inod : MLGN)
    if (basis == 0 || this->getNodeType(++i) == bType)
    {
#ifdef INDEX_CHECK
      if (inod < 1 || inod > (int)madof.size())
      {
        std::cerr <<" *** ASMbase::injectNodalVec: Node "<< inod
                  <<" is out of range [1,"<< madof.size() <<"]."<< std::endl;
        return false;
      }
#endif
      int idof = madof[inod-1] - 1;
      int ndof = madof[inod] - 1 - idof;
      if (ndof == 0)
        continue; // DOF-less node
#ifdef INDEX_CHECK
      else if (ndof < 0 || ldof+ndof > nodeVec.size())
      {
        std::cerr <<" *** ASMbase::injectNodalVec: Local DOF "<< ldof+ndof
                  <<" is out of range [1,"<< nodeVec.size() <<"]"<< std::endl;
        return false;
      }
      else if (idof+ndof > (int)globVec.size())
      {
        std::cerr <<" *** ASMbase::injectNodalVec: Global DOF "<< idof+ndof
                  <<" is out of range [1,"<< globVec.size() <<"]"<< std::endl;
        return false;
      }
#endif
      std::copy(nodeVec.begin()+ldof, nodeVec.begin()+ldof+ndof,
                globVec.begin()+idof);
      ldof += ndof;
    }

  return true;
}


void ASMbase::extractNodeVec (const RealArray& globRes, RealArray& nodeVec,
                              unsigned char nndof, int basis) const
{
  if (nndof == 0) nndof = nf;

  // Don't extract the Lagrange multipliers, if any
  size_t nNod = myLMs.empty() ? MLGN.size() : *myLMs.begin()-1;

  nodeVec.resize(nndof*nNod);
  double* nodeP = nodeVec.data();
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
    memcpy(nodeP,globRes.data()+nndof*n,nndof*sizeof(double));
  }
}


bool ASMbase::injectNodeVec (const RealArray& nodeVec, RealArray& globRes,
                             unsigned char nndof, int) const
{
  if (nndof == 0) nndof = nf;

  // Don't inject the Lagrange multipliers, if any
  size_t nNod = myLMs.empty() ? MLGN.size() : *myLMs.begin()-1;

  if (nodeVec.size() != nNod*nndof)
  {
    std::cerr <<" *** ASMbase::injectNodeVec:: Invalid patch vector, size = "
              << nodeVec.size() <<" != "<< nNod*nndof << std::endl;
    return false;
  }

  const double* nodeP = nodeVec.data();
  for (size_t i = 0; i < nNod; i++, nodeP += nndof)
  {
    int n = MLGN[i];
    if (n > 0 && nndof*(size_t)n <= globRes.size())
      memcpy(globRes.data()+nndof*(n-1),nodeP,nndof*sizeof(double));
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
      sField.fillColumn(i+1,locSol.data()+nf*nodes[i]-nf);

  return true;
}


/*!
  This method calculates the updated (control point) coordinates of current
  element by adding the current displacement vector stored as the first vector
  in \a eVec to the nodal coordinates of the reference configuration.
  The result is stored as the last vector in \a eVec.
*/

bool ASMbase::deformedConfig (const RealArray& Xnod,
                              Vectors& eVec, bool force2nd) const
{
  if (eVec.empty() || eVec.front().empty())
    return true; // No primary solution yet, silently OK

  const Vector& eV = eVec.front();
  size_t nen = Xnod.size()/nsd;
  size_t npv = nen > 0 ? eV.size()/nen : 0;
  if (npv < nsd || npv*nen != eV.size() || nsd*nen != Xnod.size())
  {
    std::cerr <<" *** ASMbase::deformedConfig: Inconsistent size of "
              <<" nodal coordinate and displacement arrays "
              << eV.size() <<", "<< Xnod.size() << std::endl;
    return false;
  }

  Vector Xdef(Xnod.data(),Xnod.size());
  if (npv == nsd)
    Xdef.add(eV);
  else
  {
    size_t i, j;
    for (i = j = 0; i < Xdef.size(); i++)
    {
      Xdef[i] += eV[j+i%nsd];
      if ((i+1)%nsd == 0) j += npv;
    }
  }

  if (force2nd && eVec.size() >= 2)
    eVec[1].swap(Xdef);
  else
    eVec.push_back(Xdef);
#if SP_DEBUG > 2
  std::cout <<"Element solution vector "<< eVec.size() << eVec.back();
#endif
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


int ASMbase::getNoGaussPt (int p, bool neumann) const
{
  if (nGauss > 0 && nGauss < 11)
    return nGauss;
  else if (neumann)
    return p;

  return p + nGauss%10;
}


bool ASMbase::evalSolution (Matrix&, const Vector&, const int*, int) const
{
  return Aerror("evalSolution(Matrix&,const Vector&,const int*,int)");
}


bool ASMbase::evalSolution (Matrix&, const Vector&,
                            const RealArray*, bool, int, int) const
{
  return Aerror("evalSolution(Matrix&,const Vector&,const RealArray*,bool,int,int)");
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


bool ASMbase::evaluate (const ASMbase*, const Vector&, RealArray&, int) const
{
  return Aerror("evaluate(const ASMbase*,const Vector&,RealArray&,int)");
}


bool ASMbase::evaluate (const Field*, RealArray&, int) const
{
  return Aerror("evaluate(const Field*,RealArray&,int)");
}


bool ASMbase::evaluate (const FunctionBase*, RealArray&, int, double) const
{
  return Aerror("evaluate(const FunctionBase*,RealArray&,int,double)");
}


bool ASMbase::writeLagBasis (std::ostream& os, const char* type) const
{
  os << "# LAGRANGIAN nodes=" << this->getNoNodes()
     << " elements=" << this->getNoElms()
     << " type=" << type << "\n";
  for (size_t i = 1; i <= this->getNoNodes(); ++i)
    os << this->getCoord(i) << "\n";
  for (const IntVec& mnpc : MNPC)
    for (size_t i = 0; i < mnpc.size(); ++i)
      os << mnpc[i] << (i+1 == mnpc.size() ? '\n' : ' ');

  return true;
}
