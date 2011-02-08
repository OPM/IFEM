// $Id: VolumePatch.C,v 1.25 2010-12-07 12:56:25 kmo Exp $
//==============================================================================
//!
//! \file VolumePatch.C
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of a topological cube as a SplineVolume.
//!
//==============================================================================

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"

#include "VolumePatch.h"
#include "GaussQuadrature.h"
#include "LinEqSystem.h"
#include "ElementBlock.h"
#include "Vec3Oper.h"
#include "Tensor.h"
#include "MPC.h"
#include "SAM.h"
#include "Utilities.h"
#include <ctype.h>
#include <fstream>
#include <algorithm>


typedef Go::SplineVolume::Dmatrix DoubleMat;  //!< 2D double array
typedef Go::SplineVolume::Dvector DoubleVec;  //!< 1D double array
typedef DoubleVec::const_iterator DoubleIter; //!< Iterator over DoubleVec

int  VolumePatch::gEl = 0;
int  VolumePatch::gNod = 0;
int  VolumePatch::splineEvalMethod = 2;
bool VolumePatch::mergeDuplNodes = true;
bool VolumePatch::swapJac = false;


VolumePatch::VolumePatch (const char* fileName, bool checkRHS)
{
  E = 2.1e7;
  nu = 0.3;
  rho = 7.85e3;
  svol = 0;
  swapW = false;

  std::cout <<"\nReading patch file "<< fileName << std::endl;
  std::ifstream is(fileName);
  if (!is.good())
    std::cerr <<" *** VolumePatch: Failure opening patch file"<< std::endl;
  else if (this->read(is) && checkRHS)
    this->checkRightHandSystem();
}


VolumePatch::VolumePatch (std::istream& is, bool checkRHS)
{
  E = 2.1e7;
  nu = 0.3;
  rho = 7.85e3;
  svol = 0;
  swapW = false;

  if (this->read(is) && checkRHS)
    this->checkRightHandSystem();
}


bool VolumePatch::read (std::istream& is)
{
  Go::ObjectHeader head;
  svol = new Go::SplineVolume;
  is >> head >> *svol;

  // Eat white-space characters to see if there is more data to read
  char c;
  while (is.get(c))
    if (!isspace(c))
    {
      is.putback(c);
      break;
    }

  if (is.good() || is.eof())
    return true;

  std::cerr <<" *** VolumePatch::read: Failure reading spline data"<< std::endl;
  delete svol;
  svol = 0;
  return false;
}


bool VolumePatch::write (std::ostream& os) const
{
  if (!svol) return false;

  os <<"700 1 0 0\n";
  os << *svol;

  return os.good();
}


VolumePatch::~VolumePatch ()
{
  if (svol) delete svol;
}


void VolumePatch::clear ()
{
  // Erase spline data
  delete svol;
  svol = 0;

  // Don't erase the elements, but set them to have zero nodes
  for (size_t i = 0; i < MNPC.size(); i++) MNPC[i].clear();

  // Erase the nodes, boundary conditions and multi-point constraints
  nodeInd.clear();
  BCode.clear();
  mpcs.clear();
}


bool VolumePatch::checkRightHandSystem ()
{
  if (!svol) return false;

  // Evaluate the spline volume at its center
  DoubleVec u(1,0.5*(svol->startparam(0) + svol->endparam(0)));
  DoubleVec v(1,0.5*(svol->startparam(1) + svol->endparam(1)));
  DoubleVec w(1,0.5*(svol->startparam(2) + svol->endparam(2)));
  DoubleVec X(3), dXdu(3), dXdv(3), dXdw(3);
  svol->gridEvaluator(u,v,w,X,dXdu,dXdv,dXdw);

  // Check that |J| = (dXdu x dXdv) * dXdw > 0.0
  if (Vec3(dXdu,dXdv) * dXdw > 0.0) return false;

  // This patch has a negative Jacobian determinant. Probably it is modelled
  // in a left-hand-system. Swap the w-parameter direction to correct for this.
  svol->reverseParameterDirection(2);
  std::cout <<"\tSwapped."<< std::endl;
  return swapW = true;
}


bool VolumePatch::uniformRefine (int dir, int nInsert)
{
  if (!svol || dir < 0 || dir > 2 || nInsert < 1) return false;

  DoubleVec  extraKnots;
  DoubleIter uit = svol->basis(dir).begin();
  double ucurr, uprev = *(uit++);
  while (uit != svol->basis(dir).end())
  {
    ucurr = *(uit++);
    if (ucurr > uprev)
      for (int i = 0; i < nInsert; i++)
      {
	double xi = (double)(i+1)/(double)(nInsert+1);
	extraKnots.push_back(ucurr*xi + uprev*(1.0-xi));
      }
    uprev = ucurr;
  }

  svol->insertKnot(dir,extraKnots);
  return true;
}


bool VolumePatch::raiseOrder (int ru, int rv, int rw)
{
  if (!svol) return false;

  svol->raiseOrder(ru,rv,rw);
  return true;
}


bool VolumePatch::generateFEMTopology ()
{
  if (!svol) return false;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  if (!nodeInd.empty())
    return nodeInd.size() == n1*n2*n3;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1 <<" "<< n2 <<" "<< n3 << std::endl;
  std::cout <<"order: "<< p1 <<" "<< p2 <<" "<< p3 << std::endl;
  for (int d = 0; d < 3; d++)
  {
    std::cout <<'d'<< char('u'+d) <<':';
    for (int i = 0; i < svol->numCoefs(d); i++)
      std::cout <<' '<< svol->knotSpan(d,i);
    std::cout << std::endl;
  }
#endif
  // Consistency checks, just to be fool-proof
  if (n1 <  2 || n2 <  2 || n3 <  2) return false;
  if (p1 <  1 || p1 <  1 || p3 <  1) return false;
  if (p1 > n1 || p2 > n2 || p3 > n3) return false;

  nodeInd.resize(n1*n2*n3);
  MNPC.resize((n1-p1+1)*(n2-p2+1)*(n3-p3+1));
  MLGE.resize(MNPC.size());

  int iel = 0;
  int inod = 0;
  for (int i3 = 1; i3 <= n3; i3++)
    for (int i2 = 1; i2 <= n2; i2++)
      for (int i1 = 1; i1 <= n1; i1++)
      {
	nodeInd[inod].I = i1-1;
	nodeInd[inod].J = i2-1;
	nodeInd[inod].K = i3-1;
	if (i1 >= p1 && i2 >= p2 && i3 >= p3)
	{
	  MNPC[iel].resize(p1*p2*p3,0);
	  for (int j3 = 0; j3 < p3; j3++)
	    for (int j2 = 0; j2 < p2; j2++)
	      for (int j1 = 0; j1 < p1; j1++)
	      {
		int gnod = inod - n1*n2*j3 - n1*j2 - j1;
		int lnod = p1*p2*j3 + p1*j2 + j1;
		MNPC[iel][lnod] = gnod;
	      }
	  MLGE[iel++] = ++gEl; // global element number over all patches
	}
	nodeInd[inod++].global = ++gNod; // global node number over all patches
      }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< iel <<" NNOD = "<< inod << std::endl;
#endif
  return true;
}


int VolumePatch::getNodeID (int inod) const
{
  if (inod < 0)
    return -inod;
  else if (inod < 1 || inod > nodeInd.size())
    return 0;

  return nodeInd[inod-1].global;
}


int VolumePatch::getNodeIndex (int globalNum) const
{
  for (int inod = 0; inod < nodeInd.size(); inod++)
    if (nodeInd[inod].global == globalNum)
      return inod+1;

  return 0;
}


/*!
  \brief An unary function that checks whether a DOF object matches the fixed
  status of a BC object.
*/

class fixed : public std::unary_function<const VolumePatch::BC&,bool>
{
  const MPC::DOF& myDof; //!< The DOF object to compare with
public:
  //! \brief Constructor initializing the myDof reference.
  fixed (const MPC::DOF& slaveDof) : myDof(slaveDof) {}
  //! \brief Returns \e true if \a myDof has the same fixed status as \a bc.
  bool operator() (const VolumePatch::BC& bc)
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


bool VolumePatch::addMPC (MPC* mpc)
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


bool VolumePatch::addSPC (int node, int dir, double value)
{
  return this->addMPC(new MPC(node,dir,value));
}


bool VolumePatch::addPeriodicity (int master, int slave, int dir)
{
  slave = this->getNodeID(slave);
  master = this->getNodeID(master);
  if (slave < 1 || master < 1) return false;

  MPC* mpc = new MPC(slave,dir);
  mpc->addMaster(master,dir);
  if (this->addMPC(mpc))
    return true;

  // Try to swap master and slave
  mpc = new MPC(master,dir);
  mpc->addMaster(slave,dir);
  if (this->addMPC(mpc))
    return true;

  std::cerr <<" *** VolumePatch::addPeriodicity: Failed to connect nodes "
	    << master <<" and "<< slave <<" in direction "<< dir << std::endl;
  return false;
}


void VolumePatch::makePeriodic (int master, int slave, int code)
{
  switch (code)
    {
    case 1:
    case 2:
    case 3:
      this->addPeriodicity(master,slave,code);
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


void VolumePatch::prescribe (int node, int code, double value)
{
  node = this->getNodeID(node);
  if (node < 1) return;

  switch (code)
    {
    case 1:
    case 2:
    case 3:
      this->addSPC(node,code,value);
      break;
    case 12:
    case 21:
      for (int dir = 1; dir <= 2; dir++)
	this->addSPC(node,dir,value);
      break;
    case 13:
    case 31:
      for (int dir = 1; dir <= 3; dir += 2)
	this->addSPC(node,dir,value);
      break;
    case 23:
    case 32:
      for (int dir = 2; dir <= 3; dir++)
	this->addSPC(node,dir,value);
      break;
    default:
      for (int dir = 1; dir <= 3; dir++)
	this->addSPC(node,dir,value);
    }
}


void VolumePatch::fix (int node, int code)
{
  node = this->getNodeID(node);
  if (node < 1) return;

  switch (code)
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
  std::cout <<"\tFixed node: "<< node <<" "<< code << std::endl;
#endif
}


bool VolumePatch::connectPatch (int face, VolumePatch& neighbor,
				int nface, int norient)
{
  if (!svol) return false;

  if (swapW) // Account for swapped parameter direction
    if (face == 5)
      face = 6;
    else if (face == 6)
      face = 5;

  if (neighbor.swapW) // Account for swapped parameter direction
    if (nface == 5)
      nface = 6;
    else if (nface == 6)
      nface = 5;

  // Set up the slave node numbers for this volume patch

  int n1 = svol->numCoefs(0);
  int n2 = svol->numCoefs(1);
  int n3 = svol->numCoefs(2);
  int node = 1, i1 = 0, i2 = 0;

  switch (face)
    {
    case 2: // Positive I-direction
      node = n1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      n2 = n3;
      break;

    case 4: // Positive J-direction
      node = n1*(n2-1)+1;
    case 3: // Negative J-direction
      i2 = n1*(n2-1);
      i1 = 1;
      n2 = n3;
      break;

    case 6: // Positive K-direction
      node = n1*n2*(n3-1)+1;
    case 5: // Negative K-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** VolumePatch::connectPatch: Invalid slave face "
		<< face << std::endl;
      return false;
    }

  int i, j;
  IntMat slaveNodes(n1,IntVec(n2,0));
  for (j = 0; j < n2; j++, node += i2)
    for (i = 0; i < n1; i++, node += i1)
      slaveNodes[i][j] = node;

  // Set up the master node numbers for the neighboring volume patch

  n1 = neighbor.svol->numCoefs(0);
  n2 = neighbor.svol->numCoefs(1);
  n3 = neighbor.svol->numCoefs(2);
  node = 1; i1 = i2 = 0;

  switch (nface)
    {
    case 2: // Positive I-direction
      node = n1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      n2 = n3;
      break;

    case 4: // Positive J-direction
      node = n1*(n2-1)+1;
    case 3: // Negative J-direction
      i2 = n1*(n2-1);
      i1 = 1;
      n2 = n3;
      break;

    case 6: // Positive K-direction
      node = n1*n2*(n3-1)+1;
    case 5: // Negative K-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** VolumePatch::connectPatch: Invalid master face "
		<< nface << std::endl;
      return false;
    }

  bool matching = false;
  if (norient < 0 || norient > 7)
  {
    std::cerr <<" *** VolumePatch::connectPatch: Orientation flag "
	      << norient <<" is out of range [0,7]"<< std::endl;
    return false;
  }
  else if (norient < 4)
    matching = (n1 == slaveNodes.size() && n2 == slaveNodes.front().size());
  else
    matching = (n2 == slaveNodes.size() && n1 == slaveNodes.front().size());
  if (!matching)
  {
    std::cerr <<" *** VolumePatch::connectPatch: Non-matching faces, sizes "
	      << n1 <<","<< n2 <<" and "
	      << slaveNodes.size() <<","<< slaveNodes.front().size()
	      << std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  for (j = 0; j < n2; j++, node += i2)
    for (i = 0; i < n1; i++, node += i1)
    {
      int k = i, l = j;
      switch (norient)
	{
	case  1: k =    i  ; l = n2-j-1; break;
	case  2: k = n1-i-1; l =    j  ; break;
	case  3: k = n1-i-1; l = n2-j-1; break;
	case  4: k =    j  ; l =    i  ; break;
	case  5: k =    j  ; l = n1-i-1; break;
	case  6: k = n2-j-1; l =    i  ; break;
	case  7: k = n2-j-1; l = n1-i-1; break;
	default: k =    i  ; l = j     ;
	}
      if (!neighbor.getCoord(node).equal(this->getCoord(slaveNodes[k][l]),xtol))
      {
	std::cerr <<" *** VolumePatch::connectPatch: Non-matching nodes "
		  << node <<": "<< neighbor.getCoord(node)
		  <<"\n                                               and "
		  << slaveNodes[k][l] <<": "<< this->getCoord(slaveNodes[k][l])
		  << std::endl;
	return false;
      }
      else if (mergeDuplNodes)
	VolumePatch::mergeNodes(neighbor.nodeInd[node-1],
				nodeInd[slaveNodes[k][l]-1]);
      else
	this->makePeriodic(-neighbor.getNodeID(node),slaveNodes[k][l]);
    }

  return true;
}


void VolumePatch::closeFaces (int dir)
{
  if (!svol) return;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  int master = 1;

  switch (dir)
    {
    case 1: // Faces are closed in I-direction
      for (int i3 = 1; i3 <= n3; i3++)
	for (int i2 = 1; i2 <= n2; i2++, master += n1)
	  this->makePeriodic(master,master+n1-1);
      break;

    case 2: // Faces are closed in J-direction
      for (int i3 = 1; i3 <= n3; i3++, master += n1*(n2-1))
	for (int i1 = 1; i1 <= n1; i1++, master++)
	  this->makePeriodic(master,master+n1*(n2-1));
      break;

    case 3: // Faces are closed in K-direction
      for (int i2 = 1; i2 <= n2; i2++)
	for (int i1 = 1; i1 <= n1; i1++, master++)
	  this->makePeriodic(master,master+n1*n2*(n3-1));
      break;
    }
}


void VolumePatch::constrainFace (int dir, int dof, double value)
{
  if (!svol) return;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  int node = 1;

  if (swapW) // Account for swapped parameter direction
    if (dir == 3 || dir == -3) dir = -dir;

  switch (dir)
    {
    case  1: // Right face (positive I-direction)
      node += n1-1;
    case -1: // Left face (negative I-direction)
      for (int i3 = 1; i3 <= n3; i3++)
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      break;

    case  2: // Back face (positive J-direction)
      node += n1*(n2-1);
    case -2: // Front face (negative J-direction)
      for (int i3 = 1; i3 <= n3; i3++, node += n1*(n2-1))
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      break;

    case  3: // Top face (positive K-direction)
      node += n1*n2*(n3-1);
    case -3: // Bottom face (negative K-direction)
      for (int i2 = 1; i2 <= n2; i2++)
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      break;
    }
}


void VolumePatch::constrainLine (int fdir, int ldir, double xi,
				 int dof, double value)
{
  if (!svol) return;
  if (xi < 0.0 || xi > 1.0) return;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  int node = 1;

  if (swapW) // Account for swapped parameter direction
    if (fdir == 3 || fdir == -3) fdir = -fdir;

  switch (fdir)
    {
    case  1: // Right face (positive I-direction)
      node += n1-1;
    case -1: // Left face (negative I-direction)
      if (ldir == 2)
      {
	// Line goes in J-direction
	node += n1*n2*int(0.5+(n3-1)*(swapW ? 1.0-xi : xi));
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      }
      else if (ldir == 3)
      {
	// Line goes in K-direction
	node += n1*int(0.5+(n2-1)*xi);
	for (int i3 = 1; i3 <= n3; i3++, node += n1*n2)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      }
      break;

    case  2: // Back face (positive J-direction)
      node += n1*(n2-1);
    case -2: // Front face (negative J-direction)
      if (ldir == 1)
      {
	// Line goes in I-direction
	node += n1*n2*int(0.5+(n3-1)*(swapW ? 1.0-xi : xi));
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      }
      else if (ldir == 3)
      {
	// Line goes in K-direction
	node += int(0.5+(n1-1)*xi);
	for (int i3 = 1; i3 <= n3; i3++, node += n1*n2)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      }
      break;

    case  3: // Top face (positive K-direction)
      node += n1*n2*(n3-1);
    case -3: // Bottom face (negative K-direction)
      if (ldir == 1)
      {
	// Line goes in I-direction
	node += n1*int(0.5+(n2-1)*xi);
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      }
      else if (ldir == 2)
      {
	// Line goes in J-direction
	node += int(0.5+(n1-1)*xi);
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
	  if (value == 0.0)
	    this->fix(node,dof);
	  else
	    this->prescribe(node,dof,value);
      }
      break;
    }
}


void VolumePatch::constrainCorner (int I, int J, int K, int dof, double value)
{
  if (!svol) return;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  if (swapW) // Account for swapped parameter direction
    K = -K;

  int node = 1;
  if (I > 0) node += n1-1;
  if (J > 0) node += n1*(n2-1);
  if (K > 0) node += n1*n2*(n3-1);

  if (value == 0.0)
    this->fix(node,dof);
  else
    this->prescribe(node,dof,value);
}


void VolumePatch::constrainNode (double xi, double eta, double zeta,
				 int dof, double value)
{
  if (!svol) return;
  if (xi   < 0.0 || xi   > 1.0) return;
  if (eta  < 0.0 || eta  > 1.0) return;
  if (zeta < 0.0 || zeta > 1.0) return;

  if (swapW) // Account for swapped parameter direction
    zeta = 1.0-zeta;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  int node = 1;
  if (xi   > 0.0) node += int(0.5+(n1-1)*xi);
  if (eta  > 0.0) node += n1*int(0.5+(n2-1)*eta);
  if (zeta > 0.0) node += n1*n2*int(0.5+(n3-1)*zeta);

  if (value == 0.0)
    this->fix(node,dof);
  else
    this->prescribe(node,dof,value);
}


void VolumePatch::resolveMPCchains (const std::vector<VolumePatch*>& model)
{
  MPCSet allMPCs;
  for (size_t i = 0; i < model.size(); i++)
    allMPCs.insert(model[i]->begin_MPC(),model[i]->end_MPC());

  int nresolved = 0;
  for (MPCSet::iterator cit = allMPCs.begin(); cit != allMPCs.end(); cit++)
    if (VolumePatch::resolveMPCchain(allMPCs,*cit)) nresolved++;

  if (nresolved > 0)
    std::cout <<"Resolved "<< nresolved <<" MPC chains."<< std::endl;
}


// Recursive method to resolve (possibly multi-level) chaining in multi-point
// constraint equations (MPCs). If a master dof in one MPC is specified as a
// slave by another MPC, it is replaced by the master(s) of that other equation.

bool VolumePatch::resolveMPCchain (const MPCSet& allMPCs, MPC* mpc)
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
      VolumePatch::resolveMPCchain(allMPCs,*cit);

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


void VolumePatch::mergeNodes (IJK& node1, IJK& node2)
{
  if (node1.global > node2.global)
    node1.global = node2.global;
  else if (node2.global > node1.global)
    node2.global = node1.global;
}


bool VolumePatch::mergeNodes (int node, int globalNum)
{
  if (node < 1 || node > nodeInd.size())
    return false;
  else if (nodeInd[node-1].global <= globalNum)
    return false;

  nodeInd[node-1].global = globalNum;
  return true;
}


int VolumePatch::renumberNodes (const std::vector<VolumePatch*>& model)
{
  int nnod = 0;
  int renum = 0;
  size_t i, j;
  std::map<int,int> old2new;
  for (i = 0; i < model.size(); i++)
    for (j = 0; j < model[i]->nodeInd.size(); j++)
      if (utl::renumber(model[i]->nodeInd[j].global,nnod,old2new))
	renum++;

  if (renum > 0)
  {
    for (i = 0; i < model.size(); i++)
      model[i]->renumberNodes(old2new,false);
    std::cout <<"\nRenumbered "<< renum <<" nodes"<< std::endl;
  }

  return nnod;
}


bool VolumePatch::renumberNodes (const std::map<int,int>& old2new, bool silent)
{
  int invalid = 0;
  for (size_t j = 0; j < BCode.size(); j++)
    if (!utl::renumber(BCode[j].node,old2new))
      invalid++;

  for (MPCSet::iterator mit = mpcs.begin(); mit != mpcs.end(); mit++)
    invalid += (*mit)->renumberNodes(old2new);

  if (invalid == 0 || silent) return true;

  std::cerr <<" *** "<< invalid <<" invalid nodes found while renumbering\n";
  return false;
}


#define DERR -999.99

double VolumePatch::getParametricVolume (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > MNPC.size())
  {
    std::cerr <<" *** VolumePatch::getParametricVolume: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif

  int inod1 = MNPC[iel-1][0];
#ifdef INDEX_CHECK
  if (inod1 < 0 || inod1 >= nodeInd.size())
  {
    std::cerr <<" *** VolumePatch::getParametricVolume: Node index "<< inod1
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
    return DERR;
  }
#endif

  double du = svol->knotSpan(0,nodeInd[inod1].I);
  double dv = svol->knotSpan(1,nodeInd[inod1].J);
  double dw = svol->knotSpan(2,nodeInd[inod1].K);
  return du*dv*dw;
}


double VolumePatch::getParametricArea (int iel, int dir) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > MNPC.size())
  {
    std::cerr <<" *** VolumePatch::getParametricArea: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif

  int inod1 = MNPC[iel-1][0];
#ifdef INDEX_CHECK
  if (inod1 < 0 || inod1 >= nodeInd.size())
  {
    std::cerr <<" *** VolumePatch::getParametricArea: Node index "<< inod1
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
    return DERR;
  }
#endif

  const int ni = nodeInd[inod1].I;
  const int nj = nodeInd[inod1].J;
  const int nk = nodeInd[inod1].K;
  switch (dir)
    {
    case 1: return svol->knotSpan(1,nj)*svol->knotSpan(2,nk);
    case 2: return svol->knotSpan(0,ni)*svol->knotSpan(2,nk);
    case 3: return svol->knotSpan(0,ni)*svol->knotSpan(1,nj);
    }

  std::cerr <<" *** VolumePatch::getParametricArea: Invalid face direction "
	    << dir << std::endl;
  return DERR;
}


bool VolumePatch::getMaterialMatrix (Matrix& C, bool inverse) const
{
  C.resize(6,6,true);

  if (nu < 0.0 || nu >= 0.5)
  {
    std::cerr <<" *** VolumePatch::getMaterialMatrix: Poisson's ratio "<< nu
	      <<" out of range [0,0.5>."<< std::endl;
    return false;
  }

  const double one  = 1.0;
  const double half = 0.5;
  const double two  = 2.0;
  const double fact = E / ((one + nu) * (one - nu - nu));

  C(1,1) = inverse ? one / E : (one - nu) * fact;
  C(2,1) = inverse ? -nu / E : nu * fact;
  C(3,1) = C(2,1);

  C(1,2) = C(2,1);
  C(2,2) = C(1,1);
  C(3,2) = C(2,1);

  C(1,3) = C(2,1);
  C(2,3) = C(2,1);
  C(3,3) = C(1,1);

  C(4,4) = inverse ? (two + nu + nu) / E : E / (two + nu + nu);
  C(5,5) = C(4,4);
  C(6,6) = C(4,4);

  return true;
}


int VolumePatch::coeffInd (int inod) const
{
#ifdef INDEX_CHECK
  if (inod < 0 || inod >= nodeInd.size())
  {
    std::cerr <<" *** VolumePatch::coeffInd: Node index "<< inod
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
    return -1;
  }
#endif

  const int ni = nodeInd[inod].I;
  const int nj = nodeInd[inod].J;
  const int nk = nodeInd[inod].K;
  return (nk*svol->numCoefs(1) + nj)*svol->numCoefs(0) + ni;
}


bool VolumePatch::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > MNPC.size())
  {
    std::cerr <<" *** VolumePatch::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  const IntVec& mnpc = MNPC[iel-1];
  X.resize(3,mnpc.size());

  DoubleIter cit = svol->coefs_begin();
  for (int n = 0; n < mnpc.size(); n++)
  {
    int ip = this->coeffInd(mnpc[n])*3;
    if (ip < 0) return false;

    for (int i = 0; i < 3; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


void VolumePatch::getNodalCoordinates (Matrix& X) const
{
  X.resize(3,nodeInd.size());

  DoubleIter cit = svol->coefs_begin();
  for (int inod = 0; inod < nodeInd.size(); inod++)
  {
    int ip = this->coeffInd(inod)*3;
    for (int i = 0; i < 3; i++)
      X(i+1,inod+1) = *(cit+(ip+i));
  }
}


Vec3 VolumePatch::getCoord (int inod) const
{
  int ip = this->coeffInd(inod-1)*3;
  if (ip < 0) return Vec3();

  DoubleIter cit = svol->coefs_begin() + ip;
  return Vec3(*cit,*(cit+1),*(cit+2));
}



/*!
  \brief Auxilliary function to compute GoTools basis function indices.
*/

static void scatterInd (int n1, int n2, int n3,
			int p1, int p2, int p3,
			int* start, IntVec& index)
{
  index.reserve(p1*p2*p3);
  int ip = ((start[2]-p3+1)*n2 + (start[1]-p2+1))*n1 + (start[0]-p1+1);
  for (int i3 = 0; i3 < p3; i3++, ip += n1*(n2-p2))
    for (int i2 = 0; i2 < p2; i2++, ip += n1-p1)
      for (int i1 = 0; i1 < p1; i1++, ip++)
	index.push_back(ip);
}


void VolumePatch::formBmatrix (Matrix& B, const Matrix& dNdX)
{
  const int nne = dNdX.rows();
  B.resize(18,nne,true);

  // Strain-displacement matrix for volume elements:
  //
  //         | d/dx   0     0   |
  //         |  0    d/dy   0   |
  //   [B] = |  0     0    d/dz | * [N]
  //         | d/dy  d/dx   0   |
  //         | d/dz   0    d/dx |
  //         |  0    d/dz  d/dy |

#define index(i,j) i+6*(j-1)
  for (int i = 1; i <= nne; i++)
  {
    // Normal strain part
    B(index(1,1),i) = dNdX(i,1);
    B(index(2,2),i) = dNdX(i,2);
    B(index(3,3),i) = dNdX(i,3);
    // Shear strain part
    B(index(4,1),i) = dNdX(i,2);
    B(index(4,2),i) = dNdX(i,1);
    B(index(5,1),i) = dNdX(i,3);
    B(index(5,3),i) = dNdX(i,1);
    B(index(6,2),i) = dNdX(i,3);
    B(index(6,3),i) = dNdX(i,2);
  }

  B.resize(6,3*nne);
#undef index
}


#if SP_DEBUG > 3
std::ostream& operator<<(std::ostream& os, const Go::BasisPts& basis)
{
  os <<" : (u,v,w) = "
     << basis.param[0] <<" "
     << basis.param[1] <<" "
     << basis.param[2] <<"  left_idx = "
     << basis.left_idx[0] <<" "
     << basis.left_idx[1] <<" "
     << basis.left_idx[2] << std::endl;
  for (unsigned int i = 0; i < basis.basisValues.size(); i++)
    os << 1+i <<'\t'<< basis.basisValues[i] << std::endl;
  return os;
}

std::ostream& operator<<(std::ostream& os, const Go::BasisDerivs& bder)
{
  os <<" : (u,v,w) = "
     << bder.param[0] <<" "
     << bder.param[1] <<" "
     << bder.param[2] <<"  left_idx = "
     << bder.left_idx[0] <<" "
     << bder.left_idx[1] <<" "
     << bder.left_idx[2] << std::endl;
  for (unsigned int i = 0; i < bder.basisValues.size(); i++)
    os << 1+i <<'\t'
       << bder.basisValues[i] <<'\t'
       << bder.basisDerivs_u[i] <<'\t'
       << bder.basisDerivs_v[i] <<'\t'
       << bder.basisDerivs_w[i] << std::endl;
  return os;
}
#endif


bool VolumePatch::assembleSystem (LinEqSystem& sys, const SAM& sam,
				  const Vec3& gravity, int nGauss,
				  const Vector& displ)
{
  if (!svol) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  int dir;
  Matrix gpar[3];
  for (dir = 0; dir < 3; dir++)
  {
    int pm1 = svol->order(dir) - 1;
    DoubleIter uit = svol->basis(dir).begin() + pm1;
    double ucurr, uprev = *(uit++);
    int nCol = svol->numCoefs(dir) - pm1;
    gpar[dir].resize(nGauss,nCol);
    for (int j = 1; j <= nCol; uit++, j++)
    {
      ucurr = *uit;
      for (int i = 1; i <= nGauss; i++)
	gpar[dir](i,j) = 0.5*((ucurr-uprev)*xg[i-1] + ucurr+uprev);
      uprev = ucurr;
    }
  }

  // Evaluate the shape function derivatives at all integration points
  std::cout <<"Spline evaluation "
	    << gpar[0].size() <<" "<< gpar[1].size() <<" "<< gpar[2].size()
	    <<" ... "<< std::flush;
  DoubleMat N0, dN1, dN2, dN3;
  std::vector<Go::BasisDerivs> spline;
  if (splineEvalMethod == 1)
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],N0,dN1,dN2,dN3);
  else
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
  std::cout << std::endl;

#if SP_DEBUG > 3
  unsigned int i, j;
  for (i = 0; i < spline.size(); i++)
    std::cout <<"\nShape functions for point "<< 1+i << spline[i];
  for (i = 0; i < N0.size(); i++)
  {
    std::cout <<"\nShape functions for point "<< 1+i << std::endl;
    for (j = 0; j < N0[i].size(); j++)
      std::cout << j+1 <<'\t'<< N0[i][j] <<'\t'
		<< dN1[i][j] <<'\t'<< dN2[i][j] <<'\t'<< dN3[i][j] << std::endl;
  }
#endif

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  Vector N(p1*p2*p3), stress;
  Matrix dNdu(p1*p2*p3,3);
  Matrix dNdX(p1*p2*p3,3);
  Matrix Xnod(3,p1*p2*p3);
  Matrix Jac(3,3);
  Matrix EK, EM, ES, B, CB, C, Dtmp;

  // Set up the consitutive matrix
  this->getMaterialMatrix(C);

  // Also assemble geometric stiffness?
  bool assembleKg = sys.M && !displ.empty();

  // Also assemble mass matrix?
  bool assembleMass = sys.M && rho > 0.0;
  if (assembleKg) assembleMass = false;

  // Also assemble load vector?
  bool assembleLoad = gravity.x != 0.0 || gravity.y != 0.0 || gravity.z != 0.0;
  if (sys.RHS.empty()) assembleLoad = false;

  // Set up the body force vector due to gravity loading
  Vector fb(3);
  if (assembleLoad && rho > 0.0)
    for (dir = 0; dir < 3; dir++)
      fb[dir] = rho*gravity[dir];


  // === Assembly loop over all elements in the patch ==========================

  int iel = 1, counter = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
	// Check that the current element has nonzero volume
	double dV = this->getParametricVolume(iel);
	if (dV == DERR) return false; // topology error (probably logic error)
	if (dV <= 0.0)  continue;     // zero volume in the parametric domain

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	const IntVec& mnpc = MNPC[iel-1];

	bool addTo = false;
	if (assembleMass || assembleKg) EM.resize(3*p1*p2*p3,3*p1*p2*p3,true);

	// --- Integration loop over all Gauss points in each direction --------

	int ip = (((i3-p3)*nGauss*nel2 + i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
	for (int k = 0; k < nGauss; k++, ip += nGauss*(nel2-1)*nGauss*nel1)
	  for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
	    for (int i = 0; i < nGauss; i++, ip++)
	    {
	      // Weight of current integration point
	      double weight = 0.125*dV*wg[i]*wg[j]*wg[k];

	      // Fetch shape function derivatives at current integration point
	      for (int n = 0; n < mnpc.size(); n++)
		if (splineEvalMethod == 1)
		{
		  int jp = this->coeffInd(mnpc[n]);
		  if (jp < 0) return false;

		     N(1+n)   =  N0[ip][jp];
		  dNdu(1+n,1) = dN1[ip][jp];
		  dNdu(1+n,2) = dN2[ip][jp];
		  dNdu(1+n,3) = dN3[ip][jp];
		}
		else
		{
		  int jp = mnpc.size()-n-1;

		     N(1+n)   = spline[ip].basisValues[jp];
		  dNdu(1+n,1) = spline[ip].basisDerivs_u[jp];
		  dNdu(1+n,2) = spline[ip].basisDerivs_v[jp];
		  dNdu(1+n,3) = spline[ip].basisDerivs_w[jp];
		}
#if SP_DEBUG > 4
	      std::cout <<"\niel, ip = "<< iel <<" "<< ip << N << std::endl;
#endif

	      // Compute Jacobian determinant and inverse
	      Jac.multiply(Xnod,dNdu); // Jac = Xnod * dNdu
	      double detJ = swapJac ? -Jac.inverse() : Jac.inverse();

	      // Compute strain-displacement matrix B from dNdX = dNdu * Jac^-1
	      formBmatrix(B,dNdX.multiply(dNdu,Jac));

	      // Integrate the element stiffness matrix
	      CB.multiply(C,B).multiply(detJ*weight); // CB = C*B*|J|*w
	      EK.multiply(B,CB,true,false,addTo);     // EK += B^T * CB

	      // Integrate geometric stiffness contributions
	      if (assembleKg)
	      {
		utl::gather(mnpc,3,displ,Dtmp);
		CB.multiply(Dtmp,stress); // Sigma = CB * Dtmp
		SymmTensor Sigma(stress);
		for (int a = 1; a <= N.size(); a++)
		  for (int b = 1; b <= N.size(); b++)
		  {
		    double kg = 0.0;
		    for (int c = 1; c <= 3; c++)
		      for (int d = 1; d <= 3; d++)
			kg += dNdX(a,c)*Sigma(c,d)*dNdX(b,d);
		    for (int d = 1; d <= 3; d++)
		      EM(3*(a-1)+d,3*(b-1)+d) += kg;
		  }
	      }

	      // Integrate the element mass matrix
	      else if (assembleMass)
	      {
		double rhow = rho*detJ*weight;
		for (int a = 1; a <= N.size(); a++)
		  for (int b = 1; b <= N.size(); b++)
		    for (int d = 1; d <= 3; d++)
		      EM(3*(a-1)+d,3*(b-1)+d) += rhow*N(a)*N(b);
	      }

	      // Integrate body force vector
	      if (assembleLoad && rho > 0.0)
		ES.outer_product(fb*detJ*weight,N,addTo);

	      addTo = true;
	    }

	// Assembly of system stiffness, mass and load
#if SP_DEBUG > 2
	std::cout <<"Stiffness matrix for element "<< iel << EK << std::endl;
	if (assembleKg)
	  std::cout <<"Geometric stiffness for "<< iel << EM << std::endl;
	else if (assembleMass)
	  std::cout <<"Mass matrix for element "<< iel << EM << std::endl;
#elif defined(SP_DEBUG)
	std::cout <<"Assembling stiffness for element "<< iel << std::endl;
#else
	std::cout << iel << ((counter++)%10 ? " " : "\n") << std::flush;
#endif
	bool status = false;
	if (!sys.RHS.empty())
	  status = sam.assembleSystem(*sys.K,sys.RHS,EK,MLGE[iel-1]);
	else
	  status = sam.assembleSystem(*sys.K,EK,MLGE[iel-1]);

	if ((assembleMass || assembleKg) && status)
	  status = sam.assembleSystem(*sys.M,EM,MLGE[iel-1]);

	if (assembleLoad && status && rho > 0.0)
	  status = sam.assembleSystem(sys.RHS,(const Vector&)ES,MLGE[iel-1]);

	if (!status) return false;
      }

  if (counter%10 != 1) std::cout << std::endl;
  return true;
}


bool VolumePatch::assembleForces (SystemVector& S, const SAM& sam,
				  const TractionFunc& t, int dir, int nGauss,
				  std::map<Vec3,Vec3>* trac)
{
  if (!svol) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch face
  Matrix gpar[3];
  for (int d = 0; d < 3; d++)
    if (-1-d == dir)
    {
      gpar[d].resize(1,1);
      gpar[d](1,1) = svol->startparam(d);
    }
    else if (1+d == dir)
    {
      gpar[d].resize(1,1);
      gpar[d](1,1) = svol->endparam(d);
    }
    else
    {
      int pm1 = svol->order(d) - 1;
      DoubleIter uit = svol->basis(d).begin() + pm1;
      double ucurr, uprev = *(uit++);
      int nCol = svol->numCoefs(d) - pm1;
      gpar[d].resize(nGauss,nCol);
      for (int j = 1; j <= nCol; uit++, j++)
      {
	ucurr = *uit;
	for (int i = 1; i <= nGauss; i++)
	  gpar[d](i,j) = 0.5*((ucurr-uprev)*xg[i-1] + ucurr+uprev);
	uprev = ucurr;
      }
    }

  // Evaluate the shape function derivatives at all integration points
  std::cout <<"Spline evaluation "
	    << gpar[0].size() <<" "<< gpar[1].size() <<" "<< gpar[2].size()
	    <<" ... "<< std::flush;
  DoubleMat N0, dN1, dN2, dN3;
  std::vector<Go::BasisDerivs> spline;
  if (splineEvalMethod == 1)
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],N0,dN1,dN2,dN3);
  else
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
  std::cout << std::endl;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  Vector N(p1*p2*p3);
  Matrix dNdu(p1*p2*p3,3);
  Matrix dNdX(p1*p2*p3,3);
  Matrix Xnod(3,p1*p2*p3);
  Matrix Jac(3,3);
  Vector ES;
  Vec3   T, X, normal;


  // === Assembly loop over all elements on the patch face =====================

  int iel = 1, counter = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
	// Skip elements that are not on current boundary face
	bool skipMe = false;
	switch (dir)
	  {
	  case -1: if (i1 > p1) skipMe = true; break;
	  case  1: if (i1 < n1) skipMe = true; break;
	  case -2: if (i2 > p2) skipMe = true; break;
	  case  2: if (i2 < n2) skipMe = true; break;
	  case -3: if (i3 > p3) skipMe = true; break;
	  case  3: if (i3 < n3) skipMe = true; break;
	  }
	if (skipMe) continue;

	// Check that the current element has nonzero face area
	double dA = this->getParametricArea(iel,abs(dir));
	if (dA == DERR) return false; // topology error (probably logic error)
	if (dA <= 0.0)  continue;     // zero area in the parametric domain

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	const IntVec& mnpc = MNPC[iel-1];

	// Define some loop control variables depending on which face we are on
	int nf1, j1, j2;
	switch (abs(dir))
	  {
	  case 1: nf1 = nel2; j2 = i3-p3; j1 = i2-p2; break;
	  case 2: nf1 = nel1; j2 = i3-p3; j1 = i1-p1; break;
	  case 3: nf1 = nel1; j2 = i2-p2; j1 = i1-p1; break;
	  }
	int t1 = 1 + abs(dir)%3; // first tangent direction
	int t2 = 1 + t1%3;       // second tangent direction

	// --- Integration loop over all Gauss points in each direction --------

	ES.resize(3*p1*p2*p3,true);
	int n, ip = (j2*nGauss*nf1 + j1)*nGauss;
	for (int j = 0; j < nGauss; j++, ip += nGauss*(nf1-1))
	  for (int i = 0; i < nGauss; i++, ip++)
	  {
	    // Weight of current integration point
	    double weight = 0.25*dA*wg[i]*wg[j];

	    // Fetch shape function values at current integration point
	    for (n = 0; n < mnpc.size(); n++)
	      if (splineEvalMethod == 1)
	      {
		int jp = this->coeffInd(mnpc[n]);
		if (jp < 0) return false;

		   N(1+n)   =  N0[ip][jp];
		dNdu(1+n,1) = dN1[ip][jp];
		dNdu(1+n,2) = dN2[ip][jp];
		dNdu(1+n,3) = dN3[ip][jp];
	      }
	      else
	      {
		int jp = mnpc.size()-n-1;

		   N(1+n)   = spline[ip].basisValues[jp];
		dNdu(1+n,1) = spline[ip].basisDerivs_u[jp];
		dNdu(1+n,2) = spline[ip].basisDerivs_v[jp];
		dNdu(1+n,3) = spline[ip].basisDerivs_w[jp];
	      }

	    // Cartesian coordinates of current integration point
	    X = Xnod * N;

	    // Compute Jacobian matrix
	    Jac.multiply(Xnod,dNdu); // Jac = Xnod * dNdu

	    // Compute the face normal
	    normal.cross(Jac.getColumn(t1),Jac.getColumn(t2));
	    double dS = normal.normalize();
	    if (dir < 0) normal *= -1.0;

	    // Evaluate the surface traction
	    T = t(X,normal);
#if SP_DEBUG > 3
	    std::cout <<"iel, ip "<< iel <<" "<< ip
		      <<" : T = "<< T << std::endl;
#endif
	    // Store the traction value for vizualization
	    if (trac && !T.isZero())
	      trac->insert(std::make_pair(X,T));

	    // Integrate the force vector
	    int idof, d;
	    T *= weight*dS;
	    for (idof = n = 1; n <= mnpc.size(); n++)
	      for (d = 0; d < 3; d++, idof++)
		ES(idof) += T[d]*N(n);
	  }

	// Assembly of system force vector
#if SP_DEBUG > 2
	std::cout <<"Force vector for element "<< iel << ES << std::endl;
#elif defined(SP_DEBUG)
	std::cout <<"Assembling pressure force for element "<< iel << std::endl;
#else
	std::cout << iel << ((counter++)%10 ? " " : "\n") << std::flush;
#endif
	if (!sam.assembleSystem(S,ES,MLGE[iel-1])) return false;
      }

  if (counter%10 != 1) std::cout << std::endl;
  return true;
}


bool VolumePatch::solutionNorms (Vector& gNorm, Matrix& eNorm, int nGauss,
				 const Vector& displ, const TensorFunc* sol)
{
  if (!svol) return true; // silently ignore empty patches
  if (displ.empty()) return false;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  int dir;
  Matrix gpar[3];
  for (dir = 0; dir < 3; dir++)
  {
    int pm1 = svol->order(dir) - 1;
    DoubleIter uit = svol->basis(dir).begin() + pm1;
    double ucurr, uprev = *(uit++);
    int nCol = svol->numCoefs(dir) - pm1;
    gpar[dir].resize(nGauss,nCol);
    for (int j = 1; j <= nCol; uit++, j++)
    {
      ucurr = *uit;
      for (int i = 1; i <= nGauss; i++)
	gpar[dir](i,j) = 0.5*((ucurr-uprev)*xg[i-1] + ucurr+uprev);
      uprev = ucurr;
    }
  }

  // Evaluate the shape function derivatives at all integration points
  DoubleMat N0, dN1, dN2, dN3;
  std::vector<Go::BasisDerivs> spline;
  if (splineEvalMethod == 1)
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],N0,dN1,dN2,dN3);
  else
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  Vector N(p1*p2*p3), sigma, sigmah;
  Matrix dNdu(p1*p2*p3,3);
  Matrix dNdX(p1*p2*p3,3);
  Matrix Xnod(3,p1*p2*p3);
  Matrix Jac(3,3);
  Matrix B, CB, C, Cinv, Dtmp;
  double pnorm[3];

  // Set up the consitutive matrix and its inverse
  this->getMaterialMatrix(C,false);
  this->getMaterialMatrix(Cinv,true);

  // === Integration loop over all elements in the patch =======================

  int iel = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
	// Check that the current element has nonzero volume
	double dV = this->getParametricVolume(iel);
	if (dV == DERR) return false; // topology error (probably logic error)
	if (dV <= 0.0)  continue;     // zero volume in the parametric domain

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	const IntVec& mnpc = MNPC[iel-1];

	// --- Integration loop over all Gauss points in each direction --------

	pnorm[0] = pnorm[1] = pnorm[2] = 0.0;
	int ip = (((i3-p3)*nGauss*nel2 + i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
	for (int k = 0; k < nGauss; k++, ip += nGauss*(nel2-1)*nGauss*nel1)
	  for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
	    for (int i = 0; i < nGauss; i++, ip++)
	    {
	      // Weight of current integration point
	      double weight = 0.125*dV*wg[i]*wg[j]*wg[k];

	      // Fetch shape function derivatives at current integration point
	      for (int n = 0; n < mnpc.size(); n++)
		if (splineEvalMethod == 1)
		{
		  int jp = this->coeffInd(mnpc[n]);
		  if (jp < 0) return false;

		     N(1+n)   =  N0[ip][jp];
		  dNdu(1+n,1) = dN1[ip][jp];
		  dNdu(1+n,2) = dN2[ip][jp];
		  dNdu(1+n,3) = dN3[ip][jp];
		}
		else
		{
		  int jp = mnpc.size()-n-1;

		     N(1+n)   = spline[ip].basisValues[jp];
		  dNdu(1+n,1) = spline[ip].basisDerivs_u[jp];
		  dNdu(1+n,2) = spline[ip].basisDerivs_v[jp];
		  dNdu(1+n,3) = spline[ip].basisDerivs_w[jp];
		}

	      // Compute Jacobian determinant and inverse
	      Jac.multiply(Xnod,dNdu); // Jac = Xnod * dNdu
	      double detJ = swapJac ? -Jac.inverse() : Jac.inverse();

	      // Compute strain-displacement matrix B from dNdX = dNdu * Jac^-1
	      formBmatrix(B,dNdX.multiply(dNdu,Jac));

	      // Evaluate the FE stress field
	      utl::gather(mnpc,3,displ,Dtmp);
	      CB.multiply(C,B).multiply(Dtmp,sigmah); // sigmah = C * B * Dtmp

	      // Integrate the energy norm a(u^h,u^h)
	      pnorm[0] += sigmah.dot(Cinv*sigmah)*detJ*weight;
	      if (sol)
	      {
		// Evaluate the analytical stress field
		sigma = (*sol)(Xnod*N);
		// Integrate the energy norm a(u,u)
		pnorm[1] += sigma.dot(Cinv*sigma)*detJ*weight;
		// Integrate the error in energy norm a(u-u^h,u-u^h)
		sigma -= sigmah;
		pnorm[2] += sigma.dot(Cinv*sigma)*detJ*weight;
	      }
	    }

	// Accumulate element and global norms
	int i;
	for (i = 1; i <= eNorm.rows() && i <= 3; i++)
	  eNorm(i,MLGE[iel-1]) = sqrt(pnorm[i-1]);
	for (i = 1; i <= gNorm.size() && i <= 3; i++)
	  gNorm(i) += pnorm[i-1];
      }

  return true;
}


bool VolumePatch::getGridParameters (RealArray& prm, int dir,
				     int nSegPerSpan) const
{
  if (!svol) return false;

  if (nSegPerSpan < 1)
  {
    std::cerr <<" *** VolumePatch::getGridParameters: Too few knot-span points "
	      << nSegPerSpan+1 <<" in direction "<< dir << std::endl;
    return false;
  }

  DoubleIter uit = svol->basis(dir).begin();
  double ucurr, uprev = *(uit++);
  while (uit != svol->basis(dir).end())
  {
    ucurr = *(uit++);
    if (ucurr > uprev)
      for (int i = 0; i < nSegPerSpan; i++)
      {
	double xg = (double)(2*i-nSegPerSpan)/(double)nSegPerSpan;
	prm.push_back(0.5*(ucurr*(1.0+xg) + uprev*(1.0-xg)));
      }
    uprev = ucurr;
  }

  prm.push_back(svol->basis(dir).endparam());
  return true;
}


bool VolumePatch::convertToElementBlock (ElementBlock& grid,
					 const int* npe) const
{
  // Compute parameter values of the nodal points
  RealArray gpar[3];
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the spline volume at all points
  size_t nx = gpar[0].size();
  size_t ny = gpar[1].size();
  size_t nz = gpar[2].size();
  DoubleVec XYZ(3*nx*ny*nz);
  svol->gridEvaluator(gpar[0],gpar[1],gpar[2],XYZ);

  // Establish the block grid coordinates
  size_t i, j, k, l;
  grid.resize(nx,ny,nz);
  for (i = j = 0; i < grid.getNoNodes(); i++, j += 3)
    grid.setCoor(i,XYZ[j],XYZ[j+1],XYZ[j+2]);

  // Establish the block grid topology
  if (grid.getNoElmNodes() == 8) // Linear elements
  {
    int n[8], ip = 0;
    for (k = 1, n[2] = 0; k < nz; k++)
      for (j = 1, n[1] = n[2]; j < ny; j++)
      {
	n[0] = n[1];
	n[1] = n[0] + 1;
	n[2] = n[1] + nx;
	n[3] = n[1] + nx-1;
	n[4] = n[0] + nx*ny;
	n[5] = n[4] + 1;
	n[6] = n[5] + nx;
	n[7] = n[5] + nx-1;
	for (i = 1; i < nx; i++)
	  for (l = 0; l < 8; l++)
	    grid.setNode(ip++,n[l]++);
      }

    return true;
  }
  else if (grid.getNoElmNodes() == 27) // parabolic elements
    if (nx%2 && ny%2 && nz%2 && nx > 2 && ny > 2 && nz > 2)
    {
      int n[27], ip = 0;
      for (k = 1, n[23] = 0; k < nz; k += 2)
	for (j = 1, n[20] = n[23]; j < ny; j += 2)
	{
	  // Corner nodes
	  n[0] = n[20];
	  n[1] = n[0] + 2;
	  n[2] = n[1] + 2*nx;
	  n[3] = n[1] + 2*nx-2;
	  n[4] = n[0] + 2*nx*ny;
	  n[5] = n[4] + 2;
	  n[6] = n[5] + 2*nx;
	  n[7] = n[5] + 2*nx-2;
	  // Mid-edge nodes
	  n[ 8] = n[ 0] + 1;
	  n[ 9] = n[ 1] + nx;
	  n[10] = n[ 3] + 1;
	  n[11] = n[ 0] + nx;
	  n[12] = n[ 8] + 2*nx*ny;
	  n[13] = n[ 9] + 2*nx*ny;
	  n[14] = n[10] + 2*nx*ny;
	  n[15] = n[11] + 2*nx*ny;
	  n[16] = n[ 0] + nx*ny;
	  n[17] = n[ 1] + nx*ny;
	  n[18] = n[ 2] + nx*ny;
	  n[19] = n[ 3] + nx*ny;
	  // Mid-face nodes
	  n[20] = n[11] + 1;
	  n[21] = n[16] + 1;
	  n[22] = n[17] + nx;
	  n[23] = n[19] + 1;
	  n[24] = n[16] + nx;
	  n[25] = n[15] + 1;
	  // Interior node
	  n[26] = n[24] + 1;
	  // Generate elements in I-direction
	  for (i = 1; i < nx; i += 2)
	    for (l = 0; l < 27; l++)
	    {
	      grid.setNode(ip++,n[l]);
	      n[l] += 2;
	    }
	}

      return true;
    }

  std::cerr <<" *** VolumePatch::convertToElementBlock: Can not convert "
	    << nx <<"x"<< ny <<"x"<< nz <<" points into a\n"
	    <<"                                         "
	    << grid.getNoElmNodes() <<"-noded hexahedron grid."<< std::endl;
  return false;
}


void VolumePatch::extractElmRes (const Matrix& globRes, Matrix& elmRes) const
{
  elmRes.resize(globRes.rows(),MLGE.size(),true);
  int i, iel, ivel = 0;
  for (iel = 1; iel <= MLGE.size(); iel++)
    if (this->getParametricVolume(iel) > 0.0)
      for (++ivel, i = 1; i <= globRes.rows(); i++)
	elmRes(i,ivel) = globRes(i,MLGE[iel-1]);

  elmRes.resize(globRes.rows(),ivel);
}


void VolumePatch::extractSolution (const Vector& solution, Vector& displ) const
{
  displ.resize(3*nodeInd.size());
  for (size_t i = 0; i < nodeInd.size(); i++)
  {
    int n = nodeInd[i].global-1;
    for (int j = 0; j < 3; j++)
      displ[3*i+j] = solution[3*n+j];
  }
}


bool VolumePatch::evalDisplField (Matrix& dField, const Vector& displ,
				  const int* npe, const LocalSystem* cs) const
{
  // Compute parameter values of the nodal points
  RealArray gpar[3];
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the displacement field
  if (splineEvalMethod == 1) // slow and memory intensive
    return evalDisplField1 (dField,displ,gpar,cs);
  else if (splineEvalMethod == 2) // faster and less memory intensive
    return evalDisplField2 (dField,displ,gpar,cs);
  else
    return false;
}


bool VolumePatch::evalDisplField1 (Matrix& dField,
				   const Vector& displ,
				   const RealArray* gpar,
				   const LocalSystem* cs) const
{
  // Evaluate the shape functions at all points
#ifdef SP_DEBUG
  std::cout <<"Spline evaluation "
	    << gpar[0].size() <<" "<< gpar[1].size() <<" "<< gpar[2].size()
	    <<" ... "<< std::flush;
#endif
  DoubleMat N;
  svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],N);

#if SP_DEBUG > 3
  for (unsigned int k = 0; k < N.size(); k++)
  {
    std::cout <<"\nShape functions for point "<< 1+k << std::endl;
    for (unsigned int l = 0; l < N[k].size(); l++)
      std::cout << l+1 <<'\t'<< N[k][l] << std::endl;
  }
#elif defined(SP_DEBUG)
  std::cout << std::endl;
#endif

  size_t nComp = displ.size() / this->getNoNodes();
  if (nComp*this->getNoNodes() != displ.size())
    return false;

  // Fetch nodal (control point) coordinates if local coordinates are used
  Matrix Xnod;
  if (cs && nComp == 3) this->getNodalCoordinates(Xnod);

  // Evaluate the displacement field at each point
  dField.resize(nComp,N.size());
  for (size_t i = 0; i < N.size(); i++)
  {
    for (size_t j = 0; j < nComp; j++)
      dField(1+j,1+i) = displ.dot(N[i],j,nComp);
#if SP_DEBUG > 2
    std::cout <<"\nGlobal displacements at point "<< 1+i
	      <<" "<< dField.getColumn(1+i) << std::endl;
#endif

    if (cs && nComp == 3)
    {
      // Transformation to local coordinate system at current point
      Vec3 d = cs->getTmat(Xnod*N[i]) * Vec3(dField.getColumn(1+i));
      dField(1,1+i) = d.x;
      dField(2,1+i) = d.y;
      dField(3,1+i) = d.z;
    }
  }

  return true;
}


bool VolumePatch::evalDisplField2 (Matrix& dField,
				   const Vector& displ,
				   const RealArray* gpar,
				   const LocalSystem* cs) const
{
  // Evaluate the shape functions at all points
#ifdef SP_DEBUG
  std::cout <<"Spline evaluation "
	    << gpar[0].size() <<" "<< gpar[1].size() <<" "<< gpar[2].size()
	    <<" ... "<< std::flush;
#endif
  std::vector<Go::BasisPts> spline;
  svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);

#if SP_DEBUG > 3
  for (unsigned int k = 0; k < spline.size(); k++)
    std::cout <<"\nShape functions for point "<< 1+k << spline[k];
#elif defined(SP_DEBUG)
  std::cout << std::endl;
#endif

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  size_t nComp = displ.size() / this->getNoNodes();
  if (nComp*this->getNoNodes() != displ.size())
    return false;

  // Fetch nodal (control point) coordinates if local coordinates are used
  Vector Ytmp;
  Matrix Xnod, Xtmp;
  if (cs && nComp == 3) this->getNodalCoordinates(Xnod);

  // Evaluate the displacement field at each point
  dField.resize(nComp,spline.size());
  for (size_t i = 0; i < spline.size(); i++)
  {
    IntVec ip;
    scatterInd(n1,n2,n3,p1,p2,p3,spline[i].left_idx,ip);
    utl::gather(ip,nComp,displ,Xtmp);
    Xtmp.multiply(spline[i].basisValues,Ytmp);
    dField.fillColumn(1+i,Ytmp);
#if SP_DEBUG > 2
    std::cout <<"\nGlobal displacements at point "<< 1+i
	      <<" "<< dField.getColumn(1+i) << std::endl;
#endif

    if (cs && nComp == 3)
    {
      // Transformation to local coordinate system at current point
      utl::gather(ip,3,Xnod,Xtmp);
      Xtmp.multiply(spline[i].basisValues,Ytmp);

      Vec3 d = cs->getTmat(Ytmp) * Vec3(dField.getColumn(1+i));
      dField(1,1+i) = d.x;
      dField(2,1+i) = d.y;
      dField(3,1+i) = d.z;
    }
  }

  return true;
}


bool VolumePatch::evalStressField (Matrix& sField, const Vector& displ,
				   const int* npe, const LocalSystem* cs) const
{
  // Compute parameter values of the nodal points
  RealArray gpar[3];
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the stress field
  if (splineEvalMethod == 1) // slow and memory intensive
    return evalStressField1 (sField,displ,gpar,cs);
  else if (splineEvalMethod == 2) // faster and less memory intensive
    return evalStressField2 (sField,displ,gpar,cs);
  else
    return false;
}


#ifndef epsZ
//! \brief Zero tolerance for the Jacobian determinant in stress evaluation.
#define epsZ 1.0e-16
#endif

bool VolumePatch::evalStressField1 (Matrix& sField,
				    const Vector& displ,
				    const RealArray* gpar,
				    const LocalSystem* cs) const
{
  // Evaluate the shape function derivatives at all points
  std::cout <<"Spline evaluation "
	    << gpar[0].size() <<" "<< gpar[1].size() <<" "<< gpar[2].size()
	    <<" ... "<< std::flush;
  DoubleMat N, dN1, dN2, dN3;
  svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],N,dN1,dN2,dN3);
  std::cout << std::endl;

  // Fetch nodal (control point) coordinates and material matrix
  Matrix B, C, CB, Xnod, Jac(3,3);
  this->getNodalCoordinates(Xnod);
  this->getMaterialMatrix(C);

  size_t nCoeffs = Xnod.cols();
  size_t nPoints = dN1.size();

  // Evaluate the stress field at each point
  sField.resize(6,nPoints);
  Matrix dNdu(nCoeffs,3), dNdX(nCoeffs,3);
  SymmTensor Sigma(3);
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch shape function derivates at current point
    for (size_t j = 0; j < nCoeffs; j++)
#ifdef INDEX_CHECK
      if (j >= dN1[i].size())
      {
	std::cerr <<" *** VolumePatch::evalStressField: Array size mismatch "
		  << nCoeffs <<", "<< dN1[i].size() << std::endl;
	return false;
      }
      else
#endif
      {
	dNdu(1+j,1) = dN1[i][j];
	dNdu(1+j,2) = dN2[i][j];
	dNdu(1+j,3) = dN3[i][j];
      }

    // Evaluate the stress tensor
    if (Jac.multiply(Xnod,dNdu).inverse() <= epsZ) // Jac = (Xnod * dNdu)^-1
      Sigma.zero(); // Set zero stress when Jacobian is singular
    else
    {
      formBmatrix(B,dNdX.multiply(dNdu,Jac)); // B = strain-displacement matrix
      CB.multiply(C,B).multiply(displ,Sigma); // Sigma = C * B * displ

      // Congruence transformation to local coordinate system at current point
      if (cs) Sigma.transform(cs->getTmat(Xnod*N[i]));
    }
    sField.fillColumn(1+i,Sigma);
  }

  return true;
}


bool VolumePatch::evalStressField2 (Matrix& sField,
				    const Vector& displ,
				    const RealArray* gpar,
				    const LocalSystem* cs) const
{
  // Evaluate the shape function derivatives at all points
  std::cout <<"Spline evaluation "
	    << gpar[0].size() <<" "<< gpar[1].size() <<" "<< gpar[2].size()
	    <<" ... "<< std::flush;
  std::vector<Go::BasisDerivs> spline;
  svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
  std::cout << std::endl;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  // Fetch nodal (control point) coordinates and material matrix
  Matrix B, C, CB, Xnod, Xtmp, Dtmp, Jac(3,3);
  this->getNodalCoordinates(Xnod);
  this->getMaterialMatrix(C);

  size_t nCoeffs = p1*p2*p3;
  size_t nPoints = spline.size();

  // Evaluate the stress field at each point
  sField.resize(6,nPoints);
  Matrix dNdu(nCoeffs,3), dNdX(nCoeffs,3);
  SymmTensor Sigma(3);
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch shape function derivates at current point
    IntVec ip;
    scatterInd(n1,n2,n3,p1,p2,p3,spline[i].left_idx,ip);
    utl::gather(ip,3,Xnod,Xtmp);
    utl::gather(ip,3,displ,Dtmp);
    for (size_t j = 0; j < ip.size(); j++)
    {
      dNdu(1+j,1) = spline[i].basisDerivs_u[j];
      dNdu(1+j,2) = spline[i].basisDerivs_v[j];
      dNdu(1+j,3) = spline[i].basisDerivs_w[j];
    }

    // Evaluate the stress tensor
    if (Jac.multiply(Xtmp,dNdu).inverse() <= epsZ) // Jac = (Xtmp * dNdu)^-1
      Sigma.zero(); // Set zero stress when Jacobian is singular
    else
    {
      formBmatrix(B,dNdX.multiply(dNdu,Jac)); // B = strain-displacement matrix
      CB.multiply(C,B).multiply(Dtmp,Sigma);  // Sigma = C * B * Dtmp

      // Congruence transformation to local coordinate system at current point
      if (cs) Sigma.transform(cs->getTmat(Xtmp*spline[i].basisValues));
    }
    sField.fillColumn(1+i,Sigma);
  }

  return true;
}


void VolumePatch::vonMises (const Matrix& s, Vector& vms)
{
  vms.resize(s.cols());
  for (size_t i = 1; i <= vms.size(); i++)
    vms(i) = sqrt(s(1,i)*s(1,i) + s(2,i)*s(2,i) + s(3,i)*s(3,i) -
		  s(1,i)*s(2,i) - s(2,i)*s(3,i) - s(3,i)*s(1,i) +
	     3.0*(s(4,i)*s(4,i) + s(5,i)*s(5,i) + s(6,i)*s(6,i)));
}
