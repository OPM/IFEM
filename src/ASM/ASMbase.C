// $Id$
//==============================================================================
//!
//! \file ASMbase.C
//!
//! \date Sep 20 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Base class for spline-based FE assembly drivers.
//!
//==============================================================================

#include "ASMbase.h"
#include "MPC.h"
#include "Vec3.h"
#include "Utilities.h"
#include <algorithm>


bool ASMbase::fixHomogeneousDirichlet = true;


ASMbase::ASMbase (unsigned char n_p, unsigned char n_s, unsigned char n_f)
{
  nf = n_f;
  nsd = n_s > 3 ? 3 : n_s;
  ndim = n_p > nsd ? nsd : n_p;
}


void ASMbase::clear ()
{
  // Don't erase the elements, but set them to have zero nodes
  for (size_t i = 0; i < MNPC.size(); i++) MNPC[i].clear();

  // Erase the nodes, boundary conditions and multi-point constraints
  MLGN.clear();
  BCode.clear();
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
  \brief An unary function that checks whether a DOF object matches the fixed
  status of a BC object.
*/

class fixed : public std::unary_function<const ASMbase::BC&,bool>
{
  const MPC::DOF& myDof; //!< The DOF object to compare with

public:
  //! \brief Constructor initializing the myDof reference.
  fixed(const MPC::DOF& slaveDof) : myDof(slaveDof) {}
  //! \brief Returns \e true if \a myDof has the same fixed status as \a bc.
  bool operator()(const ASMbase::BC& bc)
  {
    if (bc.node == myDof.node)
      switch (myDof.dof)
	{
	case 1: return bc.CX == 0;
	case 2: return bc.CY == 0;
	case 3: return bc.CZ == 0;
	}
    return false;
  }
};


bool ASMbase::addMPC (MPC* mpc, int code)
{
  if (!mpc) return true;

  // Silently ignore MPC's on dofs that already are marked as FIXED
  bool retVal = true;
  if (find_if(BCode.begin(),BCode.end(),fixed(mpc->getSlave())) == BCode.end())
    if (mpcs.insert(mpc).second)
    {
#if SP_DEBUG > 1
      std::cout <<"Added constraint: "<< *mpc;
#endif
      if (code > 0) dCode[mpc] = code;
      return retVal;
    }
    else
    {
#ifdef SP_DEBUG
      std::cout <<"Ignored constraint (duplicated slave): "<< *mpc;
#endif
      retVal = false; // This dof is already a slave in another MPC
    }

  delete mpc;
  return retVal;
}


bool ASMbase::addSPC (int node, int dir, int code)
{
  if (dir < 1 || dir > nsd) return true;
  return this->addMPC(new MPC(node,dir), code);
}


bool ASMbase::addPeriodicity (size_t master, size_t slave, int dir)
{
  if (dir < 1 || dir > nsd) return true;

  int slaveNode  = this->getNodeID(slave);
  int masterNode = this->getNodeID(master);
  if (slaveNode < 1 || masterNode < 1)
  {
    std::cerr <<" *** ASMbase::addPeriodicity: Invalid node indices "
	      << master <<", "<< slave << std::endl;
    return false;
  }

  MPC* mpc = new MPC(slaveNode,dir);
  mpc->addMaster(masterNode,dir);
  if (this->addMPC(mpc))
    return true;

  // Try to swap master and slave
  mpc = new MPC(masterNode,dir);
  mpc->addMaster(slaveNode,dir);
  if (this->addMPC(mpc))
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
      for (int dir = 1; dir <= 3; dir++)
	this->addPeriodicity(master,slave,dir);
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


void ASMbase::fix (size_t inod, int dirs)
{
  int node = this->getNodeID(inod);
  if (node < 1) return;

  switch (dirs)
    {
    case 1:
      BCode.push_back(BC(node,0,1,1));
      break;
    case 2:
      BCode.push_back(BC(node,1,0,1));
      break;
    case 3:
      BCode.push_back(BC(node,1,1,0));
      break;
    case 12:
    case 21:
      BCode.push_back(BC(node,0,0,1));
      break;
    case 13:
    case 31:
      BCode.push_back(BC(node,0,1,0));
      break;
    case 23:
    case 32:
      BCode.push_back(BC(node,1,0,0));
      break;
    default:
      BCode.push_back(BC(node,0,0,0));
    }

#if SP_DEBUG > 1
  std::cout <<"\tFixed node: "<< node <<" "<< dirs << std::endl;
#endif
}


void ASMbase::resolveMPCchains (const ASMVec& model)
{
  MPCSet allMPCs;
  for (size_t i = 0; i < model.size(); i++)
    allMPCs.insert(model[i]->begin_MPC(),model[i]->end_MPC());

  int nresolved = 0;
  for (MPCSet::iterator cit = allMPCs.begin(); cit != allMPCs.end(); cit++)
    if (ASMbase::resolveMPCchain(allMPCs,*cit)) nresolved++;

  if (nresolved > 0)
    std::cout <<"Resolved "<< nresolved <<" MPC chains."<< std::endl;
}


// Recursive method to resolve (possibly multi-level) chaining in multi-point
// constraint equations (MPCs). If a master dof in one MPC is specified as a
// slave by another MPC, it is replaced by the master(s) of that other equation.

bool ASMbase::resolveMPCchain (const MPCSet& allMPCs, MPC* mpc)
{
  bool resolved = false;
  for (size_t i = 0; i < mpc->getNoMaster();)
  {
    MPC master(mpc->getMaster(i).node,mpc->getMaster(i).dof);
    MPCSet::iterator cit = allMPCs.find(&master);
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


bool ASMbase::updateDirichlet (const std::map<int,RealFunc*>& func, double time)
{
  std::map<int,RealFunc*>::const_iterator fit;
  for (MPCMap::iterator cit = dCode.begin(); cit != dCode.end(); cit++)
    if ((fit = func.find(cit->second)) == func.end())
    {
      std::cerr <<" *** ASMbase::updateDirichlet: Code "<< cit->second
		<<" is not associated with a scalar function"<< std::endl;
      return false;
    }
    else
    {
      size_t inod = this->getNodeIndex(cit->first->getSlave().node);
      if (inod < 1) return false;

      RealFunc& g = *fit->second;
      Vec4 X(this->getCoord(inod),time);
      cit->first->setSlaveCoeff(g(X));
    }

  return true;
}


void ASMbase::collapseNodes (int& node1, int& node2)
{
  if (node1 > node2)
    node1 = node2;
  else if (node2 > node1)
    node2 = node1;
}


bool ASMbase::mergeNodes (size_t inod, int globalNum)
{
  if (inod < 1 || inod > MLGN.size())
    return false;
  else if (MLGN[--inod] <= globalNum)
    return false;

  MLGN[inod] = globalNum;
  return true;
}


int ASMbase::renumberNodes (const ASMVec& model, IntVec* l2gn)
{
  int nnod = 0;
  int renum = 0;
  size_t i, j;
  std::map<int,int> old2new;
  for (i = 0; i < model.size(); i++)
    for (j = 0; j < model[i]->MLGN.size(); j++)
      if (utl::renumber(model[i]->MLGN[j],nnod,old2new))
	renum++;

  if (renum > 0)
  {
    for (i = 0; i < model.size(); i++)
      model[i]->renumberNodes(old2new,false);
    std::cout <<"\nRenumbered "<< renum <<" nodes"<< std::endl;
  }

  if (l2gn)
  {
    l2gn->resize(old2new.size(),0);
    std::map<int,int>::const_iterator it;
    for (it = old2new.begin(); it != old2new.end(); it++)
      (*l2gn)[it->second-1] = it->first;
  }

  return nnod;
}


bool ASMbase::renumberNodes (const std::map<int,int>& old2new, bool silent)
{
#ifdef SP_DEBUG
  bool printInvalidNodes = true;
#else
  bool printInvalidNodes = false;
#endif

  int invalid = 0;
  for (size_t j = 0; j < BCode.size(); j++)
    if (!utl::renumber(BCode[j].node,old2new,printInvalidNodes))
      invalid++;

  for (MPCSet::iterator mit = mpcs.begin(); mit != mpcs.end(); mit++)
    invalid += (*mit)->renumberNodes(old2new,printInvalidNodes);

  if (invalid == 0 || silent) return true;

  std::cerr <<" *** "<< invalid <<" invalid nodes found while renumbering\n";
  return false;
}


void ASMbase::extractElmRes (const Matrix& globRes, Matrix& elmRes) const
{
  elmRes.resize(globRes.rows(),MLGE.size(),true);

  size_t i, iel, ivel = 0;
  for (iel = 0; iel < MLGE.size(); iel++)
    if (MLGE[iel] > 0)
      for (++ivel, i = 1; i <= globRes.rows(); i++)
	elmRes(i,ivel) = globRes(i,MLGE[iel]);

  elmRes.resize(globRes.rows(),ivel);
}


void ASMbase::extractNodeVec (const Vector& globRes, Vector& nodeVec,
			      unsigned char nndof) const
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


void ASMbase::injectNodeVec (const Vector& nodeVec, Vector& globRes,
			     unsigned char nndof) const
{
  if (nndof == 0) nndof = nf;

  for (size_t i = 0; i < MLGN.size(); i++)
  {
    int n = MLGN[i]-1;
    for (unsigned char j = 0; j < nndof; j++)
      globRes[nndof*n+j] = nodeVec[nndof*i+j];
  }
}


bool ASMbase::tesselate (ElementBlock&, const int*) const
{
  std::cerr <<" *** ASMBase::tesselate: Must be implemented in sub-class."
	    << std::endl;
  return false;
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


bool ASMbase::evalSolution (Matrix&, const Integrand&, const int*) const
{
  std::cerr <<" *** ASMBase::evalSolution: Must be implemented in sub-class."
	    << std::endl;
  return false;
}


bool ASMbase::evalSolution (Matrix&, const Integrand&,
			    const RealArray*, bool) const
{
  std::cerr <<" *** ASMBase::evalSolution: Must be implemented in sub-class."
	    << std::endl;
  return false;
}
