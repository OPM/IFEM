// $Id$
//==============================================================================
//!
//! \file ASMs3D.C
//!
//! \date Dec 10 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 3D spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/VolumeInterpolator.h"

#include "ASMs3D.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Vec3Oper.h"
#include <ctype.h>
#include <fstream>


ASMs3D::ASMs3D (const char* fName, bool checkRHS, unsigned char n_f)
  : ASMstruct(3,3,n_f), svol(0), swapW(false)
{
  if (fName)
  {
    std::cout <<"\nReading patch file "<< fName << std::endl;
    std::ifstream is(fName);
    if (!is.good())
      std::cerr <<" *** ASMs3D: Failure opening patch file"<< std::endl;
    else if (this->read(is) && checkRHS)
      this->checkRightHandSystem();
  }
}


ASMs3D::ASMs3D (std::istream& is, bool checkRHS, unsigned char n_f)
  : ASMstruct(3,3,n_f), svol(0), swapW(false)
{
  if (this->read(is) && checkRHS)
    this->checkRightHandSystem();
}


bool ASMs3D::read (std::istream& is)
{
  if (svol) delete svol;

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

  if (!is.good() && !is.eof())
  {
    std::cerr <<" *** ASMs3D::read: Failure reading spline data"<< std::endl;
    delete svol;
    svol = 0;
    return false;
  }
  else if (svol->dimension() < 3)
  {
    std::cerr <<" *** ASMs3D::read: Invalid spline volume patch, dim="
	      << svol->dimension() << std::endl;
    delete svol;
    svol = 0;
    return false;
  }

  geo = svol;
  return true;
}


bool ASMs3D::write (std::ostream& os, int) const
{
  if (!svol) return false;

  os <<"700 1 0 0\n";
  os << *svol;

  return os.good();
}


void ASMs3D::clear ()
{
  // Erase spline data
  if (svol) delete svol;
  svol = 0;
  geo = 0;

  // Erase the FE data
  nodeInd.clear();
  ASMbase::clear();
}


bool ASMs3D::checkRightHandSystem ()
{
  if (!svol) return false;

  // Evaluate the spline volume at its center
  RealArray u(1,0.5*(svol->startparam(0) + svol->endparam(0)));
  RealArray v(1,0.5*(svol->startparam(1) + svol->endparam(1)));
  RealArray w(1,0.5*(svol->startparam(2) + svol->endparam(2)));
  RealArray X(3), dXdu(3), dXdv(3), dXdw(3);
  svol->gridEvaluator(u,v,w,X,dXdu,dXdv,dXdw);

  // Check that |J| = (dXdu x dXdv) * dXdw > 0.0
  if (Vec3(dXdu,dXdv) * dXdw > 0.0) return false;

  // This patch has a negative Jacobian determinant. Probably it is modelled
  // in a left-hand-system. Swap the w-parameter direction to correct for this.
  svol->reverseParameterDirection(2);
  std::cout <<"\tSwapped."<< std::endl;
  return swapW = true;
}


bool ASMs3D::refine (int dir, const RealArray& xi)
{
  if (!svol || dir < 0 || dir > 2 || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;

  RealArray extraKnots;
  RealArray::const_iterator uit = svol->basis(dir).begin();
  double ucurr, uprev = *(uit++);
  while (uit != svol->basis(dir).end())
  {
    ucurr = *(uit++);
    if (ucurr > uprev)
      for (size_t i = 0; i < xi.size(); i++)
	if (i > 0 && xi[i] < xi[i-1])
	  return false;
	else
	  extraKnots.push_back(ucurr*xi[i] + uprev*(1.0-xi[i]));

    uprev = ucurr;
  }

  svol->insertKnot(dir,extraKnots);
  return true;
}


bool ASMs3D::uniformRefine (int dir, int nInsert)
{
  if (!svol || dir < 0 || dir > 2 || nInsert < 1) return false;

  RealArray extraKnots;
  RealArray::const_iterator uit = svol->basis(dir).begin();
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


bool ASMs3D::raiseOrder (int ru, int rv, int rw)
{
  if (!svol) return false;

  svol->raiseOrder(ru,rv,rw);
  return true;
}


bool ASMs3D::generateFEMTopology ()
{
  if (!svol) return false;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  if (!nodeInd.empty())
  {
    if (nodeInd.size() == (size_t)n1*n2*n3) return true;
    std::cerr <<" *** ASMs3D::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "<< n1*n2*n3
	      <<" in the patch."<< std::endl;
    return false;
  }

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1 <<" "<< n2 <<" "<< n3;
  std::cout <<"\norder: "<< p1 <<" "<< p2 <<" "<< p3;
  for (int d = 0; d < 3; d++)
  {
    std::cout <<"\nd"<< char('u'+d) <<':';
    for (int i = 0; i < svol->numCoefs(d); i++)
      std::cout <<' '<< svol->knotSpan(d,i);
  }
  std::cout << std::endl;
#endif
  // Consistency checks, just to be fool-proof
  if (n1 <  2 || n2 <  2 || n3 <  2) return false;
  if (p1 <  1 || p1 <  1 || p3 <  1) return false;
  if (p1 > n1 || p2 > n2 || p3 > n3) return false;

  MLGE.resize((n1-p1+1)*(n2-p2+1)*(n3-p3+1),0);
  MLGN.resize(n1*n2*n3);
  MNPC.resize(MLGE.size());
  nodeInd.resize(MLGN.size());

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
	  if (svol->knotSpan(0,i1-1) > 0.0)
	    if (svol->knotSpan(1,i2-1) > 0.0)
	      if (svol->knotSpan(2,i3-1) > 0.0)
	      {
		MLGE[iel] = ++gEl; // global element number over all patches
		MNPC[iel].resize(p1*p2*p3,0);

		int lnod = 0;
		for (int j3 = p3-1; j3 >= 0; j3--)
		  for (int j2 = p2-1; j2 >= 0; j2--)
		    for (int j1 = p1-1; j1 >= 0; j1--)
		      MNPC[iel][lnod++] = inod - n1*n2*j3 - n1*j2 - j1;
	      }

	  iel++;
	}
	MLGN[inod++] = ++gNod; // global node number over all patches
      }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< iel <<" NNOD = "<< inod << std::endl;
#endif
  return true;
}


int ASMs3D::Edge::next ()
{
  int ret = icnod;
  icnod += incr;

  return ret;
}


int ASMs3D::Face::next ()
{
  int ret = isnod;
  isnod += incrI;

  if (++indxI >= nnodI-1)
  {
    indxI = 1;
    isnod += incrJ - incrI*(nnodI-2);
  }

  return ret;
}


int ASMs3D::BlockNodes::next ()
{
  int ret = iinod;
  iinod += inc[0];

  if (++indxI >= nnodI-1)
  {
    indxI = 1;
    iinod += inc[1] - inc[0]*(nnodI-2);
    if (++indxJ >= nnodJ-1)
    {
      indxJ = 1;
      iinod += inc[2] - inc[1]*(nnodJ-2);
    }
  }

  return ret;
}


bool ASMs3D::assignNodeNumbers (BlockNodes& nodes, int basis)
{
  int n1, n2, n3;
  if (!this->getSize(n1,n2,n3,basis))
    return false;

  int m1 = 0, m2 = 0, m3 = 0;
  if (basis > 0)
    if (!this->getSize(m1,m2,m3,3-basis))
      return false;

  if (MLGN.size() != (size_t)(n1*n2*n3+m1*m2*m3)) return false;

  nodes.faces[0].nnodI = nodes.faces[1].nnodI = n2;
  nodes.faces[2].nnodI = nodes.faces[3].nnodI = n1;
  nodes.faces[4].nnodI = nodes.faces[5].nnodI = n1;
  nodes.nnodI = n1;
  nodes.nnodJ = n2;

  if (nodes.inc[0] == 0 || nodes.inc[1] == 0 || nodes.inc[2] == 0)
  {
    nodes.inc[0] = 1;
    nodes.inc[1] = n1-2;
    nodes.inc[2] = (n1-2)*(n2-2);
  }

  int inod = basis > 1 ? m1*m2*m3 : 0;
  for (int k = 1; k <= n3; k++)
    for (int j = 1; j <= n2; j++)
      for (int i = 1; i <= n1; i++, inod++)
	if (k == 1)
	{
	  if (j == 1)
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.ibnod[0];
	    else if (i == n1)
	      MLGN[inod] = nodes.ibnod[1];
	    else
	      MLGN[inod] = nodes.edges[0].next();
	  }
	  else if (j == n2)
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.ibnod[2];
	    else if (i == n1)
	      MLGN[inod] = nodes.ibnod[3];
	    else
	      MLGN[inod] = nodes.edges[1].next();
	  }
	  else
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.edges[4].next();
	    else if (i == n1)
	      MLGN[inod] = nodes.edges[5].next();
	    else
	      MLGN[inod] = nodes.faces[4].next();
	  }
	}
	else if (k == n3)
	{
	  if (j == 1)
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.ibnod[4];
	    else if (i == n1)
	      MLGN[inod] = nodes.ibnod[5];
	    else
	      MLGN[inod] = nodes.edges[2].next();
	  }
	  else if (j == n2)
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.ibnod[6];
	    else if (i == n1)
	      MLGN[inod] = nodes.ibnod[7];
	    else
	      MLGN[inod] = nodes.edges[3].next();
	  }
	  else
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.edges[6].next();
	    else if (i == n1)
	      MLGN[inod] = nodes.edges[7].next();
	    else
	      MLGN[inod] = nodes.faces[5].next();
	  }
	}
	else
	{
	  if (j == 1)
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.edges[8].next();
	    else if (i == n1)
	      MLGN[inod] = nodes.edges[9].next();
	    else
	      MLGN[inod] = nodes.faces[2].next();
	  }
	  else if (j == n2)
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.edges[10].next();
	    else if (i == n1)
	      MLGN[inod] = nodes.edges[11].next();
	    else
	      MLGN[inod] = nodes.faces[3].next();
	  }
	  else
	  {
	    if (i == 1)
	      MLGN[inod] = nodes.faces[0].next();
	    else if (i == n1)
	      MLGN[inod] = nodes.faces[1].next();
	    else
	      MLGN[inod] = nodes.next();
	  }
	}

#if SP_DEBUG > 1
  if (basis > 0) std::cout <<"\nBasis "<< basis <<":";
  for (int i = inod-n1*n2*n3; i < inod; i++)
    std::cout <<"\nNode "<< i+1 <<"\t: "
	      << nodeInd[i].I <<" "<< nodeInd[i].J <<" "<< nodeInd[i].K
	      <<"\tglobal no. "<< MLGN[i];
  std::cout << std::endl;
#endif
  return true;
}


bool ASMs3D::connectPatch (int face, ASMs3D& neighbor, int nface, int norient)
{
  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  if (neighbor.swapW && face > 4) // Account for swapped parameter direction
    nface = 11-nface;

  return this->connectBasis(face,neighbor,nface,norient);
}


bool ASMs3D::connectBasis (int face, ASMs3D& neighbor, int nface, int norient,
			   int basis, int slave, int master)
{
  // Set up the slave node numbers for this volume patch

  int n1, n2, n3;
  if (!this->getSize(n1,n2,n3,basis)) return false;
  int node = slave+1, i1 = 0, i2 = 0;

  switch (face)
    {
    case 2: // Positive I-direction
      node += n1-1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      n2 = n3;
      break;

    case 4: // Positive J-direction
      node += n1*(n2-1);
    case 3: // Negative J-direction
      i2 = n1*(n2-1);
      i1 = 1;
      n2 = n3;
      break;

    case 6: // Positive K-direction
      node += n1*n2*(n3-1);
    case 5: // Negative K-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** ASMs3D::connectPatch: Invalid slave face "
		<< face << std::endl;
      return false;
    }

  int i, j;
  IntMat slaveNodes(n1,IntVec(n2,0));
  for (j = 0; j < n2; j++, node += i2)
    for (i = 0; i < n1; i++, node += i1)
      slaveNodes[i][j] = node;

  // Set up the master node numbers for the neighboring volume patch

  if (!neighbor.getSize(n1,n2,n3,basis)) return false;
  node = master+1; i1 = i2 = 0;

  switch (nface)
    {
    case 2: // Positive I-direction
      node += n1-1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      n2 = n3;
      break;

    case 4: // Positive J-direction
      node += n1*(n2-1);
    case 3: // Negative J-direction
      i2 = n1*(n2-1);
      i1 = 1;
      n2 = n3;
      break;

    case 6: // Positive K-direction
      node += n1*n2*(n3-1);
    case 5: // Negative K-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** ASMs3D::connectPatch: Invalid master face "
		<< nface << std::endl;
      return false;
    }

  if (norient < 0 || norient > 7)
  {
    std::cerr <<" *** ASMs3D::connectPatch: Orientation flag "
	      << norient <<" is out of range [0,7]"<< std::endl;
    return false;
  }

  int m1 = slaveNodes.size();
  int m2 = slaveNodes.front().size();
  if (norient < 4 ? (n1 != m1 || n2 != m2) : (n2 != m1 || n1 != m2))
  {
    std::cerr <<" *** ASMs3D::connectPatch: Non-matching faces, sizes "
	      << n1 <<","<< n2 <<" and "<< m1 <<","<< m2 << std::endl;
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
	std::cerr <<" *** ASMs3D::connectPatch: Non-matching nodes "
		  << node <<": "<< neighbor.getCoord(node)
		  <<"\n                                          and "
		  << slaveNodes[k][l] <<": "<< this->getCoord(slaveNodes[k][l])
		  << std::endl;
	return false;
      }
      else
	ASMbase::collapseNodes(neighbor.MLGN[node-1],MLGN[slaveNodes[k][l]-1]);
    }

  return true;
}


void ASMs3D::closeFaces (int dir, int basis, int master)
{
  int n1, n2, n3;
  if (basis < 1) basis = 1;
  if (!this->getSize(n1,n2,n3,basis)) return;

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


void ASMs3D::constrainFace (int dir, int dof, int code)
{
  int n1, n2, n3, node = 1;
  if (!this->getSize(n1,n2,n3,1)) return;

  if (swapW) // Account for swapped parameter direction
    if (dir == 3 || dir == -3) dir = -dir;

  switch (dir)
    {
    case  1: // Right face (positive I-direction)
      node += n1-1;
    case -1: // Left face (negative I-direction)
      for (int i3 = 1; i3 <= n3; i3++)
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
	  this->prescribe(node,dof,code);
      break;

    case  2: // Back face (positive J-direction)
      node += n1*(n2-1);
    case -2: // Front face (negative J-direction)
      for (int i3 = 1; i3 <= n3; i3++, node += n1*(n2-1))
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  this->prescribe(node,dof,code);
      break;

    case  3: // Top face (positive K-direction)
      node += n1*n2*(n3-1);
    case -3: // Bottom face (negative K-direction)
      for (int i2 = 1; i2 <= n2; i2++)
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  this->prescribe(node,dof,code);
      break;
    }
}


void ASMs3D::constrainEdge (int lEdge, int dof, int code)
{
  int n1, n2, n3, n, node = 1, inc = 1;
  if (!this->getSize(n1,n2,n3,1)) return;

  if (swapW && lEdge <= 8) // Account for swapped parameter direction
    lEdge += (lEdge-1)%4 < 2 ? 2 : -2;

  if (lEdge > 8)
  {
    inc = n1*n2;
    n = n3;
  }
  else if (lEdge > 4)
  {
    inc = n1;
    n = n2;
  }
  else
    n = n1;

  switch (lEdge)
    {
    case  6:
    case 10:
      node = n1;
      break;
    case  2:
    case 12:
      node = n1*(n2-1) + 1;
      break;
    case 11:
      node = n1*n2;
      break;
    case  3:
    case  7:
      node = n1*n2*(n3-1) + 1;
      break;
    case  8:
      node = n1*(n2*(n3-1) + 1);
      break;
    case  4:
      node = n1*(n2*n3-1) + 1;
      break;
    }

  for (int i = 0; i < n; i++, node += inc)
    this->prescribe(node,dof,code);
}


void ASMs3D::constrainLine (int fdir, int ldir, double xi, int dof, int code)
{
  if (xi < 0.0 || xi > 1.0) return;

  int n1, n2, n3, node = 1;
  if (!this->getSize(n1,n2,n3,1)) return;

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
	  this->prescribe(node,dof,code);
      }
      else if (ldir == 3)
      {
	// Line goes in K-direction
	node += n1*int(0.5+(n2-1)*xi);
	for (int i3 = 1; i3 <= n3; i3++, node += n1*n2)
	  this->prescribe(node,dof,code);
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
	  this->prescribe(node,dof,code);
      }
      else if (ldir == 3)
      {
	// Line goes in K-direction
	node += int(0.5+(n1-1)*xi);
	for (int i3 = 1; i3 <= n3; i3++, node += n1*n2)
	  this->prescribe(node,dof,code);
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
	  this->prescribe(node,dof,code);
      }
      else if (ldir == 2)
      {
	// Line goes in J-direction
	node += int(0.5+(n1-1)*xi);
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
	  this->prescribe(node,dof,code);
      }
      break;
    }
}


void ASMs3D::constrainCorner (int I, int J, int K, int dof, int code)
{
  int n1, n2, n3;
  if (!this->getSize(n1,n2,n3,1)) return;

  if (swapW) // Account for swapped parameter direction
    K = -K;

  int node = 1;
  if (I > 0) node += n1-1;
  if (J > 0) node += n1*(n2-1);
  if (K > 0) node += n1*n2*(n3-1);

  this->prescribe(node,dof,code);
}


void ASMs3D::constrainNode (double xi, double eta, double zeta,
			    int dof, int code)
{
  if (xi   < 0.0 || xi   > 1.0) return;
  if (eta  < 0.0 || eta  > 1.0) return;
  if (zeta < 0.0 || zeta > 1.0) return;

  if (swapW) // Account for swapped parameter direction
    zeta = 1.0-zeta;

  int n1, n2, n3;
  if (!this->getSize(n1,n2,n3,1)) return;

  int node = 1;
  if (xi   > 0.0) node += int(0.5+(n1-1)*xi);
  if (eta  > 0.0) node += n1*int(0.5+(n2-1)*eta);
  if (zeta > 0.0) node += n1*n2*int(0.5+(n3-1)*zeta);

  this->prescribe(node,dof,code);
}


#define DERR -999.99

double ASMs3D::getParametricVolume (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3D::getParametricVolume: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  int inod1 = MNPC[iel-1].back();
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nodeInd.size())
  {
    std::cerr <<" *** ASMs3D::getParametricVolume: Node index "<< inod1
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
    return DERR;
  }
#endif

  double du = svol->knotSpan(0,nodeInd[inod1].I);
  double dv = svol->knotSpan(1,nodeInd[inod1].J);
  double dw = svol->knotSpan(2,nodeInd[inod1].K);
  return du*dv*dw;
}


double ASMs3D::getParametricArea (int iel, int dir) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3D::getParametricArea: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  int inod1 = MNPC[iel-1].back();
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nodeInd.size())
  {
    std::cerr <<" *** ASMs3D::getParametricArea: Node index "<< inod1
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

  std::cerr <<" *** ASMs3D::getParametricArea: Invalid face direction "
	    << dir << std::endl;
  return DERR;
}


int ASMs3D::coeffInd (size_t inod) const
{
#ifdef INDEX_CHECK
  if (inod >= nodeInd.size())
  {
    std::cerr <<" *** ASMs3D::coeffInd: Node index "<< inod
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
    return -1;
  }
#endif

  const int ni = nodeInd[inod].I;
  const int nj = nodeInd[inod].J;
  const int nk = nodeInd[inod].K;
  return (nk*svol->numCoefs(1) + nj)*svol->numCoefs(0) + ni;
}


bool ASMs3D::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3D::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  const IntVec& mnpc = MNPC[iel-1];
  X.resize(3,mnpc.size());

  RealArray::const_iterator cit = svol->coefs_begin();
  for (size_t n = 0; n < mnpc.size(); n++)
  {
    int ip = this->coeffInd(mnpc[n])*svol->dimension();
    if (ip < 0) return false;

    for (size_t i = 0; i < 3; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


void ASMs3D::getNodalCoordinates (Matrix& X) const
{
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  X.resize(3,n1*n2*n3);

  RealArray::const_iterator cit = svol->coefs_begin();
  size_t inod = 1;
  for (int i3 = 0; i3 < n3; i3++)
    for (int i2 = 0; i2 < n2; i2++)
      for (int i1 = 0; i1 < n1; i1++, inod++)
      {
	int ip = ((i3*n2 + i2)*n1 + i1)*svol->dimension();
	for (size_t i = 0; i < 3; i++)
	  X(i+1,inod) = *(cit+(ip+i));
      }
}


Vec3 ASMs3D::getCoord (size_t inod) const
{
  if (inod == 0) return Vec3();
  int ip = this->coeffInd(inod-1)*svol->dimension();
  if (ip < 0) return Vec3();

  RealArray::const_iterator cit = svol->coefs_begin() + ip;
  return Vec3(*cit,*(cit+1),*(cit+2));
}


bool ASMs3D::getSize (int& n1, int& n2, int& n3, int) const
{
  if (!svol) return false;

  n1 = svol->numCoefs(0);
  n2 = svol->numCoefs(1);
  n3 = svol->numCoefs(2);
  return true;
}


void ASMs3D::getGaussPointParameters (Matrix& uGP, int dir, int nGauss,
				      const double* xi) const
{
  int pm1 = svol->order(dir) - 1;
  RealArray::const_iterator uit = svol->basis(dir).begin() + pm1;

  int nCol = svol->numCoefs(dir) - pm1;
  uGP.resize(nGauss,nCol);

  double ucurr, uprev = *(uit++);
  for (int j = 1; j <= nCol; uit++, j++)
  {
    ucurr = *uit;
    for (int i = 1; i <= nGauss; i++)
      uGP(i,j) = 0.5*((ucurr-uprev)*xi[i-1] + ucurr+uprev);
    uprev = ucurr;
  }
}


/*!
  \brief Computes the characteristic element length from nodal coordinates.
*/

static double getElmSize (int p1, int p2, int p3, const Matrix& X)
{
  int i, j, k, id1, id2;
  double value, v1, h = 1.0e12;

  // Z-direction
  for (i = 1; i <= p1; i++)
    for (j = 0; j < p2; j++)
    {
      id1 = j*p1 + i;
      id2 = id1 + (p3-1)*p2*p1;
      value = 0.0;
      for (k = 1; k <= 3; k++)
      {
	v1 = X(k,id2) - X(k,id1);
	value += v1*v1;
      }
      if (value < h) h = value;
    }

  // Y-direction
  for (i = 1; i <= p1; i++)
    for (k = 0; k < p2; k++)
    {
      id1 = k*p2*p1 + i;
      id2 = id1 + (p2-1)*p1;
      value = 0.0;
      for (j = 1; j <= 3; j++)
      {
	v1 = X(j,id2) - X(j,id1);
	value += v1*v1;
      }
      if (value < h) h = value;
    }

  // X-direction
  for (j = 0; j < p2; j++)
    for (k = 0; k < p3; k++)
    {
      id1 = k*p1*p2 + j*p1 + 1;
      id2 = id1 + p1 - 1;
      value = 0.0;
      for (i = 1; i <= 3; i++)
      {
	v1 = X(i,id2) - X(i,id1);
	value += v1*v1;
      }
      if (value < h) h = value;
    }

  return sqrt(h);
}


void ASMs3D::extractBasis (const Go::BasisDerivs& spline,
			   Vector& N, Matrix& dNdu)
{
  dNdu.resize(N.size(),3);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
     N  (n)   = spline.basisValues[jp];
    dNdu(n,1) = spline.basisDerivs_u[jp];
    dNdu(n,2) = spline.basisDerivs_v[jp];
    dNdu(n,3) = spline.basisDerivs_w[jp];
  }
}


void ASMs3D::extractBasis (const Go::BasisDerivs2& spline,
			   Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2)
{
   dNdu .resize(N.size(),3);
  d2Ndu2.resize(N.size(),3,3);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
      N   (n)     = spline.basisValues[jp];
     dNdu (n,1)   = spline.basisDerivs_u[jp];
     dNdu (n,2)   = spline.basisDerivs_v[jp];
     dNdu (n,3)   = spline.basisDerivs_w[jp];
    d2Ndu2(n,1,1) = spline.basisDerivs_uu[jp];
    d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline.basisDerivs_uv[jp];
    d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = spline.basisDerivs_uw[jp];
    d2Ndu2(n,2,2) = spline.basisDerivs_vv[jp];
    d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = spline.basisDerivs_vw[jp];
    d2Ndu2(n,3,3) = spline.basisDerivs_ww[jp];
  }
}


bool ASMs3D::integrate (Integrand& integrand,
			GlobalIntegral& glInt,
			const TimeDomain& time,
			const LintegralVec& locInt)
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3D::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Get the reduced integration quadrature points, if needed
  const double* xr = 0;
  int nRed = integrand.getIntegrandType() - 10;
  if (nRed < 1)
    nRed = nRed < 0 ? nGauss : 0;
  else if (!(xr = GaussQuadrature::getCoord(nRed)))
    return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gpar[3], redpar[3];
  for (int d = 0; d < 3; d++)
  {
    this->getGaussPointParameters(gpar[d],d,nGauss,xg);
    if (integrand.getIntegrandType() > 10)
      this->getGaussPointParameters(redpar[d],d,nRed,xr);
  }

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivs>  spline;
  std::vector<Go::BasisDerivs2> spline2;
  std::vector<Go::BasisDerivs>  splineRed;
  {
    PROFILE2("Spline evaluation");
    if (integrand.getIntegrandType() == 2)
      svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
    else
      svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
    if (integrand.getIntegrandType() > 10)
      svol->computeBasisGrid(redpar[0],redpar[1],redpar[2],splineRed);
  }

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  FiniteElement fe(p1*p2*p3);
  Matrix   dNdu, Xnod, Jac;
  Matrix3D d2Ndu2, Hess;
  Vec4     X;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
	fe.iel = MLGE[iel-1];
	if (fe.iel < 1) continue; // zero-volume element

	// Get element volume in the parameter space
	double dV = this->getParametricVolume(iel);
	if (dV < 0.0) return false; // topology error (probably logic error)

	// Set up control point (nodal) coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Compute characteristic element length, if needed
	if (integrand.getIntegrandType() == 2)
	  fe.h = getElmSize(p1,p2,p3,Xnod);

	else if (integrand.getIntegrandType() == 3)
	{
	  // --- Compute average value of basis functions over the element -----

	  fe.Navg.resize(p1*p2*p3,true);
	  double vol = 0.0;
	  int ip = (((i3-p3)*nGauss*nel2 + i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
	  for (int k = 0; k < nGauss; k++, ip += nGauss*(nel2-1)*nGauss*nel1)
	    for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
	      for (int i = 0; i < nGauss; i++, ip++)
	      {
		// Fetch basis function derivatives at current integration point
		extractBasis(spline[ip],fe.N,dNdu);

		// Compute Jacobian determinant of coordinate mapping
		// and multiply by weight of current integration point
		double detJac = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu,false);
		double weight = 0.125*dV*wg[i]*wg[j]*wg[k];

		// Numerical quadrature
		fe.Navg.add(fe.N,detJac*weight);
		vol += detJac*weight;
	      }

	  // Divide by element volume
	  fe.Navg /= vol;
	}

	else if (integrand.getIntegrandType() == 4)
	{
	  // Compute the element center
	  Go::Point X0;
	  double u0 = 0.5*(gpar[0](1,i1-p1+1) + gpar[0](nGauss,i1-p1+1));
	  double v0 = 0.5*(gpar[1](1,i2-p2+1) + gpar[1](nGauss,i2-p2+1));
	  double w0 = 0.5*(gpar[2](1,i3-p3+1) + gpar[2](nGauss,i3-p3+1));
	  svol->point(X0,u0,v0,w0);
	  X.x = X0[0];
	  X.y = X0[1];
	  X.z = X0[2];
	}

	// Initialize element quantities
	if (!integrand.initElement(MNPC[iel-1],X,nRed*nRed*nRed))
	  return false;

	// Caution: Unless locInt is empty, we assume it points to an array of
	// LocalIntegral pointers, of length at least the number of elements in
	// the model (as defined by the highest number in the MLGE array).
	// If the array is shorter than this, expect a segmentation fault.
	LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[fe.iel-1];


	if (integrand.getIntegrandType() > 10)
	{
	  // --- Selective reduced integration loop ----------------------------

	  int ip = (((i3-p3)*nRed*nel2 + i2-p2)*nRed*nel1 + i1-p1)*nRed;
	  for (int k = 0; k < nRed; k++, ip += nRed*(nel2-1)*nRed*nel1)
	    for (int j = 0; j < nRed; j++, ip += nRed*(nel1-1))
	      for (int i = 0; i < nRed; i++, ip++)
	      {
		// Local element coordinates of current integration point
		fe.xi   = xr[i];
		fe.eta  = xr[j];
		fe.zeta = xr[k];

		// Parameter values of current integration point
		fe.u = redpar[0](i+1,i1-p1+1);
		fe.v = redpar[1](j+1,i2-p2+1);
		fe.w = redpar[2](k+1,i3-p3+1);

		// Fetch basis function derivatives at current point
		extractBasis(splineRed[ip],fe.N,dNdu);

		// Compute Jacobian inverse and derivatives
		fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

		// Compute the reduced integration terms of the integrand
		if (!integrand.reducedInt(fe))
		  return false;
	      }
	}


	// --- Integration loop over all Gauss points in each direction --------

	int ip = (((i3-p3)*nGauss*nel2 + i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
	for (int k = 0; k < nGauss; k++, ip += nGauss*(nel2-1)*nGauss*nel1)
	  for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
	    for (int i = 0; i < nGauss; i++, ip++)
	    {
	      // Local element coordinates of current integration point
	      fe.xi   = xg[i];
	      fe.eta  = xg[j];
	      fe.zeta = xg[k];

	      // Parameter values of current integration point
	      fe.u = gpar[0](i+1,i1-p1+1);
	      fe.v = gpar[1](j+1,i2-p2+1);
	      fe.w = gpar[2](k+1,i3-p3+1);

	      // Fetch basis function derivatives at current integration point
	      if (integrand.getIntegrandType() == 2)
		extractBasis(spline2[ip],fe.N,dNdu,d2Ndu2);
	      else
		extractBasis(spline[ip],fe.N,dNdu);

 	      // Compute Jacobian inverse of coordinate mapping and derivatives
	      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
	      if (fe.detJxW == 0.0) continue; // skip singular points

 	      // Compute Hessian of coordinate mapping and 2nd order derivatives
	      if (integrand.getIntegrandType() == 2)
		if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
		  return false;

#if SP_DEBUG > 4
	      std::cout <<"\niel, ip = "<< iel <<" "<< ip
			<<"\nN ="<< fe.N <<"dNdX ="<< fe.dNdX << std::endl;
#endif

	      // Cartesian coordinates of current integration point
	      X = Xnod * fe.N;
	      X.t = time.t;

	      // Evaluate the integrand and accumulate element contributions
	      fe.detJxW *= 0.125*dV*wg[i]*wg[j]*wg[k];
	      if (!integrand.evalInt(elmInt,fe,time,X))
		return false;
	    }

	// Finalize the element quantities
	if (!integrand.finalizeElement(elmInt,time))
	  return false;

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,fe.iel))
	  return false;
      }

  return true;
}


bool ASMs3D::integrate (Integrand& integrand, int lIndex,
			GlobalIntegral& glInt,
			const TimeDomain& time,
			const LintegralVec& locInt)
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3D::integrate(B)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  // Compute parameter values of the Gauss points over the whole patch face
  Matrix gpar[3];
  for (int d = 0; d < 3; d++)
    if (-1-d == faceDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(svol->startparam(d));
    }
    else if (1+d == faceDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(svol->endparam(d));
    }
    else
      this->getGaussPointParameters(gpar[d],d,nGauss,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivs> spline;
  {
    PROFILE2("Spline evaluation");
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
  }

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  FiniteElement fe(p1*p2*p3);
  fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);
  fe.w = gpar[2](1,1);

  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   normal;


  // === Assembly loop over all elements on the patch face =====================

  int iel = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
	fe.iel = MLGE[iel-1];
	if (fe.iel < 1) continue; // zero-volume element

	// Skip elements that are not on current boundary face
	bool skipMe = false;
	switch (faceDir)
	  {
	  case -1: if (i1 > p1) skipMe = true; break;
	  case  1: if (i1 < n1) skipMe = true; break;
	  case -2: if (i2 > p2) skipMe = true; break;
	  case  2: if (i2 < n2) skipMe = true; break;
	  case -3: if (i3 > p3) skipMe = true; break;
	  case  3: if (i3 < n3) skipMe = true; break;
	  }
	if (skipMe) continue;

	// Get element face area in the parameter space
	double dA = this->getParametricArea(iel,abs(faceDir));
	if (dA < 0.0) return false; // topology error (probably logic error)

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	if (!integrand.initElementBou(MNPC[iel-1])) return false;

	// Define some loop control variables depending on which face we are on
	int nf1, j1, j2;
	switch (abs(faceDir))
	  {
	  case 1: nf1 = nel2; j2 = i3-p3; j1 = i2-p2; break;
	  case 2: nf1 = nel1; j2 = i3-p3; j1 = i1-p1; break;
	  case 3: nf1 = nel1; j2 = i2-p2; j1 = i1-p1; break;
	  default: nf1 = j1 = j2 = 0;
	  }

	// Caution: Unless locInt is empty, we assume it points to an array of
	// LocalIntegral pointers, of length at least the number of elements in
	// the model (as defined by the highest number in the MLGE array).
	// If the array is shorter than this, expect a segmentation fault.
	LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[fe.iel-1];


	// --- Integration loop over all Gauss points in each direction --------

	int k1, k2, k3;
	int ip = (j2*nGauss*nf1 + j1)*nGauss;
	for (int j = 0; j < nGauss; j++, ip += nGauss*(nf1-1))
	  for (int i = 0; i < nGauss; i++, ip++)
	  {
	    // Local element coordinates and parameter values
	    // of current integration point
	    switch (abs(faceDir)) {
	    case 1: k2 = i+1; k3 = j+1; k1 = 0; break;
	    case 2: k1 = i+1; k3 = j+1; k2 = 0; break;
	    case 3: k1 = i+1; k2 = j+1; k3 = 0; break;
	    default: k1 = k2 = k3 = 0;
	    }
	    if (gpar[0].size() > 1)
	    {
	      fe.xi = xg[k1];
	      fe.u = gpar[0](k1,i1-p1+1);
	    }
	    if (gpar[1].size() > 1)
	    {
	      fe.eta = xg[k2];
	      fe.v = gpar[1](k2,i2-p2+1);
	    }
	    if (gpar[2].size() > 1)
	    {
	      fe.zeta = xg[k3];
	      fe.w = gpar[2](k3,i3-p3+1);
	    }

	    // Fetch basis function derivatives at current integration point
	    extractBasis(spline[ip],fe.N,dNdu);

	    // Compute basis function derivatives and the face normal
	    fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
	    if (fe.detJxW == 0.0) continue; // skip singular points

	    if (faceDir < 0) normal *= -1.0;

	    // Cartesian coordinates of current integration point
	    X = Xnod * fe.N;
	    X.t = time.t;

	    // Evaluate the integrand and accumulate element contributions
	    fe.detJxW *= 0.25*dA*wg[i]*wg[j];
	    if (!integrand.evalBou(elmInt,fe,time,X,normal))
	      return false;
	  }

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,fe.iel))
	  return false;
      }

  return true;
}


bool ASMs3D::integrateEdge (Integrand& integrand, int lEdge,
			    GlobalIntegral& glInt,
			    const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3D::integrate(E)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points along the whole patch edge
  Matrix gpar[3];
  for (int d = 0; d < 3; d++)
    if (lEdge < d*4+1 || lEdge >= d*4+5)
    {
      gpar[d].resize(1,1);
      if (lEdge%4 == 1)
	gpar[d].fill(svol->startparam(d));
      else if (lEdge%4 == 0)
	gpar[d].fill(svol->endparam(d));
      else if (lEdge == 6 || lEdge == 10)
	gpar[d].fill(d == 0 ? svol->endparam(d) : svol->startparam(d));
      else if (lEdge == 2 || lEdge == 11)
	gpar[d].fill(d == 1 ? svol->endparam(d) : svol->startparam(d));
      else if (lEdge == 3 || lEdge == 7)
	gpar[d].fill(d == 2 ? svol->endparam(d) : svol->startparam(d));
    }
    else
    {
      int pm1 = svol->order(d) - 1;
      RealArray::const_iterator uit = svol->basis(d).begin() + pm1;
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

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivs> spline;
  {
    PROFILE2("Spline evaluation");
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
  }

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  FiniteElement fe(p1*p2*p3);
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);
  fe.w = gpar[2](1,1);
  if (gpar[0].size() == 1) fe.xi = fe.u == svol->startparam(0) ? -1.0 : 1.0;
  if (gpar[1].size() == 1) fe.eta = fe.v == svol->startparam(1) ? -1.0 : 1.0;
  if (gpar[2].size() == 1) fe.zeta = fe.w == svol->startparam(2) ? -1.0 : 1.0;

  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   tang;


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
	fe.iel = MLGE[iel-1];
	if (fe.iel < 1) continue; // zero-volume element

	// Skip elements that are not on current boundary edge
	bool skipMe = false;
	switch (lEdge)
	  {
	  case  1: if (i2 > p2 || i3 > p3) skipMe = true; break;
	  case  2: if (i2 < n2 || i3 > p3) skipMe = true; break;
	  case  3: if (i2 > p2 || i3 < n3) skipMe = true; break;
	  case  4: if (i2 < n2 || i3 < n3) skipMe = true; break;
	  case  5: if (i1 > p1 || i3 > p3) skipMe = true; break;
	  case  6: if (i1 < n1 || i3 > p3) skipMe = true; break;
	  case  7: if (i1 > p1 || i3 < n3) skipMe = true; break;
	  case  8: if (i1 < n1 || i3 < n3) skipMe = true; break;
	  case  9: if (i1 > p1 || i2 > p2) skipMe = true; break;
	  case 10: if (i1 < n1 || i2 > p2) skipMe = true; break;
	  case 11: if (i1 > p1 || i2 < n2) skipMe = true; break;
	  case 12: if (i1 < n1 || i2 < n2) skipMe = true; break;
	  }
	if (skipMe) continue;

	// Get element edge length in the parameter space
	double dS = 0.0;
	int ip = MNPC[iel-1].back();
#ifdef INDEX_CHECK
	if (ip < 0 || (size_t)ip >= nodeInd.size()) return false;
#endif
	if (lEdge < 5)
	{
	  dS = svol->knotSpan(0,nodeInd[ip].I);
	  ip = (i1-p1)*nGauss;
	}
	else if (lEdge < 9)
	{
	  dS = svol->knotSpan(1,nodeInd[ip].J);
	  ip = (i2-p2)*nGauss;
	}
	else if (lEdge < 13)
	{
	  dS = svol->knotSpan(2,nodeInd[ip].K);
	  ip = (i3-p3)*nGauss;
	}

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
	if (!integrand.initElementBou(MNPC[iel-1])) return false;


	// --- Integration loop over all Gauss points along the edge -----------

	LocalIntegral* elmInt = 0;
	for (int i = 0; i < nGauss; i++, ip++)
	{
	  // Parameter values of current integration point
	  if (gpar[0].size() > 1) fe.u = gpar[0](i+1,i1-p1+1);
	  if (gpar[1].size() > 1) fe.v = gpar[1](i+1,i2-p2+1);
	  if (gpar[2].size() > 1) fe.w = gpar[2](i+1,i3-p3+1);

	  // Fetch basis function derivatives at current integration point
	  extractBasis(spline[ip],fe.N,dNdu);

	  // Compute basis function derivatives and the edge tang
	  fe.detJxW = utl::Jacobian(Jac,tang,fe.dNdX,Xnod,dNdu,1+(lEdge-1)/4);
	  if (fe.detJxW == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X = Xnod * fe.N;
	  X.t = time.t;

	  // Evaluate the integrand and accumulate element contributions
	  fe.detJxW *= 0.5*dS*wg[i];
	  if (!integrand.evalBou(elmInt,fe,time,X,tang))
	    return false;
	}

	// Assembly of global system integral
	if (!glInt.assemble(elmInt,fe.iel))
	  return false;
      }

  return true;
}


int ASMs3D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!svol) return -3;

  int i;
  for (i = 0; i < 3; i++)
    param[i] = (1.0-xi[i])*svol->startparam(i) + xi[i]*svol->endparam(i);

  Go::Point X0;
  svol->point(X0,param[0],param[1],param[2]);
  for (i = 0; i < 3 && i < svol->dimension(); i++)
    X[i] = X0[i];

  // Check if this point matches any of the control points (nodes)
  Vec3 Xnod;
  size_t inod = 1;
  RealArray::const_iterator cit = svol->coefs_begin();
  for (i = 0; cit != svol->coefs_end(); cit++, i++)
  {
    if (i < 3) Xnod[i] = *cit;
    if (i+1 == svol->dimension())
      if (X.equal(Xnod,0.001))
	return inod;
      else
      {
	inod++;
	i = -1;
      }
  }

  return 0;
}


bool ASMs3D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
  if (!svol) return false;

  if (nSegPerSpan < 1)
  {
    std::cerr <<" *** ASMs3D::getGridParameters: Too few knot-span points "
	      << nSegPerSpan+1 <<" in direction "<< dir << std::endl;
    return false;
  }

  RealArray::const_iterator uit = svol->basis(dir).begin();
  double ucurr = 0.0, uprev = *(uit++);
  while (uit != svol->basis(dir).end())
  {
    ucurr = *(uit++);
    if (ucurr > uprev)
      if (nSegPerSpan == 1)
	prm.push_back(uprev);
      else for (int i = 0; i < nSegPerSpan; i++)
      {
	double xg = (double)(2*i-nSegPerSpan)/(double)nSegPerSpan;
	prm.push_back(0.5*(ucurr*(1.0+xg) + uprev*(1.0-xg)));
      }
    uprev = ucurr;
  }

  if (ucurr > prm.back())
    prm.push_back(ucurr);
  return true;
}


bool ASMs3D::getGrevilleParameters (RealArray& prm, int dir) const
{
  if (!svol) return false;

  const Go::BsplineBasis& basis = svol->basis(dir);

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}


bool ASMs3D::tesselate (ElementBlock& grid, const int* npe) const
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
  RealArray XYZ(svol->dimension()*nx*ny*nz);
  svol->gridEvaluator(gpar[0],gpar[1],gpar[2],XYZ);

  // Establish the block grid coordinates
  size_t i, j, k, l;
  grid.resize(nx,ny,nz);
  for (i = j = 0; i < grid.getNoNodes(); i++, j += svol->dimension())
    grid.setCoor(i,XYZ[j],XYZ[j+1],XYZ[j+2]);

  // Establish the block grid topology
  int nel1 = svol->numCoefs(0) - svol->order(0) + 1;
  int nel2 = svol->numCoefs(1) - svol->order(1) + 1;
  int ie, nse1 = npe[0] - 1;
  int je, nse2 = npe[1] - 1;
  int ke, nse3 = npe[2] - 1;
  int n[8], ip = 0;
  for (k = ke = 1, n[2] = 0; k < nz; k++)
  {
    for (j = je = 1, n[1] = n[2]; j < ny; j++)
    {
      n[0] = n[1];
      n[1] = n[0] + 1;
      n[2] = n[1] + nx;
      n[3] = n[1] + nx-1;
      n[4] = n[0] + nx*ny;
      n[5] = n[4] + 1;
      n[6] = n[5] + nx;
      n[7] = n[5] + nx-1;
      for (i = ie = 1; i < nx; i++)
      {
	for (l = 0; l < 8; l++)
	  grid.setNode(ip++,n[l]++);
	grid.setElmId(((k-1)*(ny-1)+j-1)*(nx-1)+i,((ke-1)*nel2+je-1)*nel1+ie);
	if (i%nse1 == 0) ie++;
      }
      if (j%nse2 == 0) je++;
    }
    if (k%nse3 == 0) ke++;
  }

  return true;
}


void ASMs3D::scatterInd (int n1, int n2, int n3, int p1, int p2, int p3,
			 const int* start, IntVec& index)
{
  index.reserve(p1*p2*p3);
  int ip = ((start[2]-p3+1)*n2 + (start[1]-p2+1))*n1 + (start[0]-p1+1);
  for (int i3 = 0; i3 < p3; i3++, ip += n1*(n2-p2))
    for (int i2 = 0; i2 < p2; i2++, ip += n1-p1)
      for (int i1 = 0; i1 < p1; i1++, ip++)
        index.push_back(ip);
}


bool ASMs3D::evalSolution (Matrix& sField, const Vector& locSol,
			   const int* npe) const
{
  // Compute parameter values of the result sampling points
  RealArray gpar[3];
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  return this->evalSolution(sField,locSol,gpar);
}


bool ASMs3D::evalSolution (Matrix& sField, const Vector& locSol,
			   const RealArray* gpar, bool regular) const
{
  // Evaluate the basis functions at all points
  std::vector<Go::BasisPts> spline;
  if (regular)
  {
    PROFILE2("Spline evaluation");
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    PROFILE2("Spline evaluation");
    spline.resize(gpar[0].size());
    for (size_t i = 0; i < spline.size(); i++)
      svol->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline[i]);
  }
  else
    return false;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  size_t nComp = locSol.size() / this->getNoNodes();
  if (nComp*this->getNoNodes() != locSol.size())
    return false;

  Matrix Xtmp;
  Vector Ytmp;

  // Evaluate the primary solution field at each point
  size_t nPoints = spline.size();
  sField.resize(nComp,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    IntVec ip;
    scatterInd(n1,n2,n3,p1,p2,p3,spline[i].left_idx,ip);

    utl::gather(ip,nComp,locSol,Xtmp);
    Xtmp.multiply(spline[i].basisValues,Ytmp);
    sField.fillColumn(1+i,Ytmp);
  }

  return true;
}


bool ASMs3D::evalSolution (Matrix& sField, const Integrand& integrand,
			   const int* npe, bool project) const
{
  if (npe)
  {
    // Compute parameter values of the result sampling points
    RealArray gpar[3];
    if (this->getGridParameters(gpar[0],0,npe[0]-1) &&
	this->getGridParameters(gpar[1],1,npe[1]-1) &&
	this->getGridParameters(gpar[2],2,npe[2]-1))
      if (project)
      {
	// Project the secondary solution onto the spline basis
	Go::SplineVolume* v = this->projectSolution(integrand);
	if (v)
	{
	  // Evaluate the projected field at the result sampling points
	  const Vector& svec = sField; // using utl::matrix cast operator
	  sField.resize(v->dimension(),
			gpar[0].size()*gpar[1].size()*gpar[2].size());
	  v->gridEvaluator(gpar[0],gpar[1],gpar[2],const_cast<Vector&>(svec));
	  delete v;
	  return true;
	}
      }
      else
	// Evaluate the secondary solution directly at all sampling points
	return this->evalSolution(sField,integrand,gpar);
  }
  else
  {
    // Project the secondary solution onto the spline basis
    Go::SplineVolume* v = this->projectSolution(integrand);
    if (v)
    {
      // Extract control point values from the spline object
      sField.resize(v->dimension(),
		    v->numCoefs(0)*v->numCoefs(1)*v->numCoefs(2));
      sField.fill(&(*v->coefs_begin()));
      delete v;
      return true;
    }
  }

  std::cerr <<" *** ASMs3D::evalSolution: Failure!"<< std::endl;
  return false;
}


Go::GeomObject* ASMs3D::evalSolution (const Integrand& integrand) const
{
  return this->projectSolution(integrand);
}


Go::SplineVolume* ASMs3D::projectSolution (const Integrand& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[3];
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGrevilleParameters(gpar[dir],dir))
      return 0;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,gpar))
    return 0;

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume the the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (svol->rational())
    svol->getWeights(weights);

  const Vector& vec = sValues;
  return Go::VolumeInterpolator::regularInterpolation(svol->basis(0),
						      svol->basis(1),
						      svol->basis(2),
						      gpar[0], gpar[1], gpar[2],
						      const_cast<Vector&>(vec),
						      sValues.rows(),
						      svol->rational(),
						      weights);
}


bool ASMs3D::evalSolution (Matrix& sField, const Integrand& integrand,
			   const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  // Evaluate the basis functions and their derivatives at all points
  std::vector<Go::BasisDerivs> spline(regular ? 0 : gpar[0].size());
  if (regular)
  {
    PROFILE2("Spline evaluation");
    svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
  }
  else if (gpar[0].size() == gpar[1].size() && gpar[0].size() == gpar[2].size())
  {
    PROFILE2("Spline evaluation");
    std::vector<Go::BasisDerivs> tmpSpline(1);
    for (size_t i = 0; i < spline.size(); i++)
    {
      svol->computeBasisGrid(RealArray(1,gpar[0][i]),
                             RealArray(1,gpar[1][i]),
                             RealArray(1,gpar[2][i]),
                             tmpSpline);
      spline[i] = tmpSpline.front();
    }
    // TODO: Request a GoTools method replacing the above:
    // void SplineVolume::computeBasisGrid(double param_u,
    //                                     double param_v,
    //                                     double param_w,
    //                                     BasisDerivs& result) const
    /*
    for (size_t i = 0; i < spline.size(); i++)
      svol->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline[i]);
    */
  }
  else
    return false;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  Vector N(p1*p2*p3), solPt;
  Matrix dNdu, dNdX, Jac;

  // Evaluate the secondary solution field at each point
  size_t nPoints = spline.size();
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntVec ip;
    scatterInd(n1,n2,n3,p1,p2,p3,spline[i].left_idx,ip);

    // Fetch associated control point coordinates
    utl::gather(ip,3,Xnod,Xtmp);

    // Fetch basis function derivatives at current integration point
    extractBasis(spline[i],N,dNdu);

    // Compute the Jacobian inverse
    if (utl::Jacobian(Jac,dNdX,Xtmp,dNdu) == 0.0) // Jac = (Xtmp * dNdu)^-1
      continue; // skip singular points

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,N,dNdX,Xtmp*N,ip))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}
