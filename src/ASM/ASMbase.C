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
#include "MPC.h"
#include "Vec3.h"
#include "Vec3Oper.h"
#include "Utilities.h"
#include <algorithm>


bool ASMbase::fixHomogeneousDirichlet = true;


ASMbase::ASMbase (unsigned char n_p, unsigned char n_s, unsigned char n_f)
  : MLGE(myMLGE), MLGN(myMLGN), MNPC(myMNPC), shareFE(false)
{
  nf = n_f;
  nsd = n_s > 3 ? 3 : n_s;
  ndim = n_p > nsd ? nsd : n_p;
}


ASMbase::ASMbase (const ASMbase& patch, unsigned char n_f)
  : MLGE(patch.myMLGE), MLGN(patch.myMLGN), MNPC(patch.myMNPC), shareFE(true)
{
  nf = n_f > 0 ? n_f : patch.nf;
  nsd = patch.nsd;
  ndim = patch.ndim;
  // Note: Properties are _not_ copied
}


ASMbase::~ASMbase ()
{
  for (MPCIter it = mpcs.begin(); it != mpcs.end(); it++)
    delete *it;
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

  myMLGN.clear();
  BCode.clear();
  dCode.clear();
  mpcs.clear();
}


size_t ASMbase::getNodeIndex (int globalNum) const
{
  IntVec::const_iterator it = find(MLGN.begin(),MLGN.end(),globalNum);
  if (it == MLGN.end()) return 0;

  return 1 + (it-MLGN.begin());
}


int ASMbase::getNodeID (size_t inod) const
{
  if (inod < 1 || inod > MLGN.size())
    return 0;

  return MLGN[inod-1];
}


int ASMbase::getElmID (size_t iel) const
{
  if (iel < 1 || iel > MLGE.size())
    return 0;

  return MLGE[iel-1];
}


size_t ASMbase::getNoElms (bool includeZeroVolumeElms) const
{
  if (includeZeroVolumeElms) return MLGE.size();

  size_t nel = 0;
  for (size_t i = 0; i < MLGE.size(); i++)
    if (MLGE[i] > 0) nel++;

  return nel;
}


void ASMbase::printNodes (std::ostream& os, const char* heading) const
{
  Matrix X;
  this->getNodalCoordinates(X);
  if (heading) os << heading;
  for (size_t inod = 1; inod <= X.cols(); inod++)
  {
    os <<'\n'<< MLGN[inod-1] <<':';
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

class fixed : public std::unary_function<const ASMbase::BC&,bool>
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
	  }
    return false;
  }
};


bool ASMbase::isFixed (int node, int ldof) const
{
  return std::find_if(BCode.begin(),BCode.end(),fixed(node,ldof)) != BCode.end();
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


bool ASMbase::addSPC (int node, int dir, int code)
{
  if (dir < 1 || dir > nf) return true;

  MPC* mpc = new MPC(node,dir);
  return this->addMPC(mpc,code);
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
  switch (dirs)
    {
    case 1:
    case 2:
    case 3:
      this->addPeriodicity(master,slave,dirs);
      break;
    case 12:
    case 21:
      for (int dir = 1; dir <= 2; dir++)
	this->addPeriodicity(master,slave,dir);
      break;
    case 13:
    case 31:
      for (int dir = 1; dir <= 3; dir += 2)
	this->addPeriodicity(master,slave,dir);
      break;
    case 23:
    case 32:
      for (int dir = 2; dir <= 3; dir++)
	this->addPeriodicity(master,slave,dir);
      break;
    default:
      // If all DOFs are going to be coupled, assign a common global node number
      ASMbase::collapseNodes(*this,master,*this,slave);
    }
}


void ASMbase::prescribe (size_t inod, int dirs, int code)
{
  if (code == 0 && fixHomogeneousDirichlet)
    return this->fix(inod,dirs);

  int node = this->getNodeID(inod);
  if (node < 1) return;

  switch (dirs)
    {
    case 1:
    case 2:
    case 3:
      this->addSPC(node,dirs,code);
      break;
    case 12:
    case 21:
      for (int dir = 1; dir <= 2; dir++)
	this->addSPC(node,dir,code);
      break;
    case 13:
    case 31:
      for (int dir = 1; dir <= 3; dir += 2)
	this->addSPC(node,dir,code);
      break;
    case 23:
    case 32:
      for (int dir = 2; dir <= 3; dir++)
	this->addSPC(node,dir,code);
      break;
    default:
      for (int dir = 1; dir <= 3; dir++)
	this->addSPC(node,dir,code);
    }
}


/*!
  \brief Equality operator for BC objects comparing node numbers only.
*/

bool operator== (const ASMbase::BC& rhs, const int& lhs)
{
  return rhs.node == lhs;
}


void ASMbase::fix (size_t inod, int dirs)
{
  int node = this->getNodeID(inod);
  if (node < 1) return;

  BCVec::iterator bit = find(BCode.begin(),BCode.end(),node);
  if (bit == BCode.end())
  {
    BCode.push_back(BC(node));
    bit = BCode.end()-1;
  }

  switch (dirs)
    {
    case 1:
      if (bit->CX == 0) return;
      bit->CX = 0;
      break;
    case 2:
      if (bit->CY == 0) return;
      bit->CY = 0;
      break;
    case 3:
      if (bit->CZ == 0) return;
      bit->CZ = 0;
      break;
    case 12:
    case 21:
      if (bit->CX + bit->CY == 0) return;
      bit->CX = bit->CY = 0;
      break;
    case 13:
    case 31:
      if (bit->CX + bit->CZ == 0) return;
      bit->CX = bit->CZ = 0;
      break;
    case 23:
    case 32:
      if (bit->CY + bit->CZ == 0) return;
      bit->CY = bit->CZ = 0;
      break;
    default:
      if (bit->CX + bit->CY + bit->CZ == 0) return;
      bit->CX = bit->CY = bit->CZ = 0;
    }

#if SP_DEBUG > 1
  std::cout <<"\tFixed node: "<< node <<" "<< dirs << std::endl;
#endif
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
    std::cout <<"Merged "<< nmerged <<" MPC equations."<< std::endl;
  if (ndeleted > 0)
    std::cout <<"Deleted "<< ndeleted <<" MPC equations."<< std::endl;

#if SP_DEBUG > 1
  if (allMPCs.empty()) return;
  std::cout <<"\nMulti-point constraints:\n";
  for (MPCIter c = allMPCs.begin(); c != allMPCs.end(); c++) std::cout << **c;
#endif
}


void ASMbase::resolveMPCchains (const MPCSet& allMPCs)
{
#if SP_DEBUG > 1
  std::cout <<"\nResolving MPC chains"<< std::endl;
  for (MPCIter c = allMPCs.begin(); c != allMPCs.end(); c++) std::cout << **c;
#endif

  int nresolved = 0;
  for (MPCIter cit = allMPCs.begin(); cit != allMPCs.end(); cit++)
    if (ASMbase::resolveMPCchain(allMPCs,*cit)) nresolved++;

  if (nresolved > 0)
    std::cout <<"Resolved "<< nresolved <<" MPC chains."<< std::endl;
}


// Recursive method to resolve (possibly multi-level) chaining in multi-point
// constraint equations (MPCs). If a master dof in one MPC is specified as a
// slave by another MPC, it is replaced by the master(s) of that other equation.

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


bool ASMbase::updateDirichlet (const std::map<int,RealFunc*>& func,
			       const std::map<int,VecFunc*>& vfunc, double time)
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

    Vec4 X(this->getCoord(inod),time);
    if ((fit = func.find(cit->second)) != func.end())
    {
      RealFunc& g = *fit->second;
      cit->first->setSlaveCoeff(g(X));
    }
    else if ((vfit = vfunc.find(cit->second)) != vfunc.end())
    {
      int idof = cit->first->getSlave().dof;
      VecFunc& g = *vfit->second;
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
    std::cout <<"  ** Merging duplicated nodes "<< globalNum <<" and "<< oldNum
	      <<" at X="<< this->getCoord(inod) << std::endl;

  std::map<int,int> old2New;
  myMLGN[inod-1] = old2New[oldNum] = globalNum;
  return this->renumberNodes(old2New);
}


int ASMbase::renumberNodes (const ASMVec& model, std::map<int,int>& old2new)
{
  ASMVec::const_iterator it;
  std::map<int,int>::iterator nit;

  for (it = model.begin(); it != model.end(); it++)
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


int ASMbase::renumberNodes (std::map<int,int>& old2new, int& nnod)
{
  int renum = 0;
  for (size_t j = 0; j < myMLGN.size(); j++)
    if (utl::renumber(myMLGN[j],nnod,old2new))
      renum++;

  if (renum == 0)
    nnod = std::max(nnod,*std::max_element(MLGN.begin(),MLGN.end()));

  return renum;
}


bool ASMbase::renumberNodes (const std::map<int,int>& old2new)
{
#ifdef SP_DEBUG
  bool printInvalidNodes = old2new.size() > 1;
#else
  bool printInvalidNodes = false;
#endif

  int invalid = 0;
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
			      unsigned char nndof, int) const
{
  if (nndof == 0) nndof = nf;

  nodeVec.resize(nndof*MLGN.size());
  for (size_t i = 0; i < MLGN.size(); i++)
  {
    int n = MLGN[i]-1;
    for (unsigned char j = 0; j < nndof; j++)
      nodeVec[nndof*i+j] = globRes[nndof*n+j];
  }
}


bool ASMbase::injectNodeVec (const Vector& nodeVec, Vector& globRes,
			     unsigned char nndof) const
{
  if (nndof == 0) nndof = nf;

  if (nodeVec.size() != MLGN.size()*nndof)
  {
    std::cerr <<" *** ASMbase::injectNodeVec:: Invalid patch vector, size = "
	      << nodeVec.size() <<" != "<< MLGN.size()*nndof << std::endl;
    return false;
  }

  for (size_t i = 0; i < MLGN.size(); i++)
  {
    int n = MLGN[i]-1;
    for (unsigned char j = 0; j < nndof; j++)
      globRes[nndof*n+j] = nodeVec[nndof*i+j];
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
    else
      sField.fillColumn(i+1,locSol.ptr()+nf*nodes[i]-nf);

  return true;
}


bool ASMbase::evalSolution (Matrix&, const Vector&, const int*) const
{
  std::cerr <<" *** ASMBase::evalSolution: Must be implemented in sub-class."
	    << std::endl;
  return false;
}


bool ASMbase::evalSolution (Matrix&, const Vector&,
			    const RealArray*, bool) const
{
  std::cerr <<" *** ASMBase::evalSolution: Must be implemented in sub-class."
	    << std::endl;
  return false;
}


bool ASMbase::evalSolution (Matrix&, const IntegrandBase&,
			    const int*, bool) const
{
  std::cerr <<" *** ASMBase::evalSolution: Must be implemented in sub-class."
	    << std::endl;
  return false;
}


bool ASMbase::evalSolution (Matrix&, const IntegrandBase&,
			    const RealArray*, bool) const
{
  std::cerr <<" *** ASMBase::evalSolution: Must be implemented in sub-class."
	    << std::endl;
  return false;
}
