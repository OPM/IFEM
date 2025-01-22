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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SurfaceInterpolator.h"

#include "ASMs3D.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "SplineField3D.h"
#include "SplineFields3D.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "Point.h"
#include "Tensor.h"
#include "MPC.h"
#include "IFEM.h"
#include <array>
#include <utility>


ASMs3D::ASMs3D (unsigned char n_f) : ASMstruct(3,3,n_f), nodeInd(myNodeInd)
{
  svol = nullptr;
  swapW = false;
}


ASMs3D::ASMs3D (const ASMs3D& patch, unsigned char n_f)
  : ASMstruct(patch,n_f), nodeInd(patch.myNodeInd)
{
  svol = patch.svol;

  swapW = patch.swapW;

  // Need to set nnod here,
  // as hasXNodes might be invoked before the FE data is generated
  if (nnod == 0 && svol)
    nnod = svol->numCoefs(0)*svol->numCoefs(1)*svol->numCoefs(2);
}


ASMs3D::ASMs3D (const ASMs3D& patch)
  : ASMstruct(patch), svol(patch.svol),
    nodeInd(myNodeInd), myNodeInd(patch.nodeInd), dirich(patch.dirich)
{
  swapW = patch.swapW;
}


Go::SplineSurface* ASMs3D::getBoundary (int dir, int)
{
  if (dir < -3 || dir == 0 || dir > 3)
    return nullptr;

  // The boundary surfaces are stored internally in the SplineVolume object
  int iface = dir > 0 ? 2*dir-1 : -2*dir-2;
  return svol->getBoundarySurface(iface).get();
}


const Go::SplineVolume* ASMs3D::getBasis (int basis) const
{
  switch (basis) {
    case ASM::GEOMETRY_BASIS:
      return static_cast<const Go::SplineVolume*>(geomB.get());
    case ASM::PROJECTION_BASIS:
      return static_cast<const Go::SplineVolume*>(projB.get());
    case ASM::PROJECTION_BASIS_2:
      return static_cast<const Go::SplineVolume*>(projB2.get());
    case ASM::REFINEMENT_BASIS:
      return nullptr;
    default:
      return svol.get();
  }
}


Go::SplineVolume* ASMs3D::getBasis (int basis)
{
  return const_cast<Go::SplineVolume*>(std::as_const(*this).getBasis(basis));
}


void ASMs3D::copyParameterDomain (const ASMbase* other)
{
  const ASMs3D* o = dynamic_cast<const ASMs3D*>(other);
  if (!o) return;

  Go::Array<double,6> pd = o->getBasis()->parameterSpan();
  this->getBasis()->setParameterDomain(pd[0],pd[1],pd[2],pd[3],pd[4],pd[5]);
}


bool ASMs3D::read (std::istream& is)
{
  if (shareFE) return true;

  Go::ObjectHeader head;
  svol = std::make_shared<Go::SplineVolume>();
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
    svol.reset();
    return false;
  }
  else if (svol->dimension() < 3)
  {
    std::cerr <<" *** ASMs3D::read: Invalid spline volume patch, dim="
              << svol->dimension() << std::endl;
    svol.reset();
    return false;
  }

  geomB = svol;
  return true;
}


bool ASMs3D::write (std::ostream& os, int basis) const
{
  if (!svol) return false;
  if (basis > static_cast<int>(this->getNoBasis())) return false;
  const Go::SplineVolume* spline = this->getBasis(basis);
  if (!spline) return false;

  os <<"700 1 0 0\n" << *spline;

  return os.good();
}


void ASMs3D::clear (bool retainGeometry)
{
  if (!retainGeometry)
  {
    // Erase spline data
    projB.reset();
    svol.reset();
    geomB.reset();
  }

  // Erase the FE data
  this->ASMbase::clear(retainGeometry);
  myNodeInd.clear();
  xnMap.clear();
  nxMap.clear();

  // Erase threading group data
  threadGroupsVol[0].clear();
  threadGroupsVol[1].clear();
  threadGroupsFace.clear();

  myCache.clear();
}


bool ASMs3D::addXElms (short int dim, short int item, size_t nXn, IntVec& nodes)
{
  if (!this->addXNodes(dim,nXn,nodes))
    return false;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  int iel = 0;
  bool skipMe = false;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
        if (MLGE[iel] < 1) continue; // Skip zero-volume element

        // Skip elements that are not on current boundary face
        switch (item)
          {
          case 1: skipMe = i1 > p1; break;
          case 2: skipMe = i1 < n1; break;
          case 3: skipMe = i2 > p2; break;
          case 4: skipMe = i2 < n2; break;
          case 5: skipMe = i3 > p3; break;
          case 6: skipMe = i3 < n3; break;
          }
        if (skipMe) continue;

        IntVec& mnpc = myMNPC[nel+iel];
        if (!mnpc.empty())
        {
          std::cerr <<" *** ASMs3D::addXElms: Only one X-face allowed."
                    << std::endl;
          return false;
        }

        mnpc = MNPC[iel]; // Copy the ordinary element nodes

        // Negate node numbers that are not on the boundary face, to flag that
        // they shall not receive any tangent and/or residual contributions
        int lnod = 0;
        for (int j3 = 0; j3 < p3; j3++)
          for (int j2 = 0; j2 < p2; j2++)
            for (int j1 = 0; j1 < p1; j1++, lnod++)
            {
              switch (item)
                {
                case 1: skipMe = j1 > 0;    break;
                case 2: skipMe = j1 < p1-1; break;
                case 3: skipMe = j2 > 0;    break;
                case 4: skipMe = j2 < p2-1; break;
                case 5: skipMe = j3 > 0;    break;
                case 6: skipMe = j3 < p3-1; break;
                }
              if (skipMe) // Hack for node 0: Using -maxint as flag instead
                mnpc[lnod] = mnpc[lnod] == 0 ? -2147483648 : -mnpc[lnod];
            }

        // Add connectivity to the extra-ordinary nodes
        for (size_t i = 0; i < nXn; i++)
          mnpc.push_back(MLGN.size()-nXn+i);

        myMLGE[nel+iel] = -(++gEl); // Extra-ordinary element => negative sign
      }

  return true;
}


size_t ASMs3D::getNodeIndex (int globalNum, bool noAddedNodes) const
{
  size_t inod = 1 + utl::findIndex(MLGN,globalNum);
  if (noAddedNodes && !xnMap.empty() && inod > nnod)
  {
    std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
    if (it != xnMap.end()) return it->second;
  }

  return inod;
}


int ASMs3D::getNodeID (size_t inod, bool noAddedNodes) const
{
  if (noAddedNodes && !nxMap.empty())
  {
    std::map<size_t,size_t>::const_iterator it = nxMap.find(inod);
    if (it != nxMap.end()) inod = it->second;
  }

  if (inod < 1 || inod > MLGN.size())
    return 0;

  return MLGN[inod-1];
}


size_t ASMs3D::getNoNodes (int basis) const
{
  size_t n = this->ASMbase::getNoNodes(basis);
  if (n > 0 || basis < 1 || !svol) return n;

  // We request the number of nodes before the FE topology has been generated
  return svol->numCoefs(0) * svol->numCoefs(1) * svol->numCoefs(2);
}


bool ASMs3D::checkRightHandSystem ()
{
  if (!svol || shareFE) return false;

  // Evaluate the spline volume at its center
  RealArray u(1,0.5*(svol->startparam(0) + svol->endparam(0)));
  RealArray v(1,0.5*(svol->startparam(1) + svol->endparam(1)));
  RealArray w(1,0.5*(svol->startparam(2) + svol->endparam(2)));
  RealArray X(3), dXdu(3), dXdv(3), dXdw(3);
  svol->gridEvaluator(u,v,w,X,dXdu,dXdv,dXdw);

  // Check that |J| = (dXdu x dXdv) * dXdw > 0.0
  if (Vec3(dXdu,dXdv) * Vec3(dXdw) > 0.0) return false;

  // This patch has a negative Jacobian determinant. Probably it is modelled
  // in a left-hand-system. Swap the w-parameter direction to correct for this.
  svol->reverseParameterDirection(2);
  return swapW = true;
}


bool ASMs3D::refine (int dir, const RealArray& xi)
{
  if (!svol || dir < 0 || dir > 2 || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = svol->basis(dir).begin();
  double uprev = *(uit++);
  while (uit != svol->basis(dir).end())
  {
    double ucurr = *(uit++);
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
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = svol->basis(dir).begin();
  double uprev = *(uit++);
  while (uit != svol->basis(dir).end())
  {
    double ucurr = *(uit++);
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


bool ASMs3D::raiseOrder (int ru, int rv, int rw, bool setOrder)
{
  if (!svol) return false;
  if (shareFE) return true;

  if (setOrder)
  {
    ru -= svol->order(0);
    rv -= svol->order(1);
    rw -= svol->order(2);
  }
  svol->raiseOrder(ru,rv,rw);
  return true;
}


bool ASMs3D::separateProjectionBasis () const
{
  return projB && projB != svol;
}


/*!
  This method is supposed to be invoked twice during the model generation.
  In the first call, with \a init = \e true, the spline volume object \a *svol
  is cloned into \a *proj and the two pointers are then swapped, such that
  the subsequent refine and raiseOrder operations will apply to the projection
  basis and not on the geometry basis.
  In the second call, the pointers are swapped back.

  The method can also be invoked twice with \a init = \e false in case the
  projection basis is to be read from a file.
*/

bool ASMs3D::createProjectionBasis (bool init)
{
  if (!svol)
    return false;
  else if (init && !projB)
    projB.reset( svol->clone());

  std::swap(geomB,projB);
  svol = std::static_pointer_cast<Go::SplineVolume>(geomB);
  return true;
}


bool ASMs3D::generateFEMTopology ()
{
  if (!svol) return false;
  if (!projB) projB = svol;

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  if (!nodeInd.empty())
  {
    nnod = n1*n2*n3;
    if (nodeInd.size() != (size_t)nnod)
    {
      std::cerr <<" *** ASMs3D::generateFEMTopology: Inconsistency between the"
                <<" number of FE nodes "<< nodeInd.size()
                <<"\n     and the number of spline coefficients "<< nnod
                <<" in the patch."<< std::endl;
      return false;
    }
    else if (shareFE == 'F')
    {
      // Must store the global node numbers anyway, in case
      // the patch sharing from gets extraordinary nodes later
      myMLGN = MLGN;
      gNod += nnod;
    }
    nel = (n1-p1+1)*(n2-p2+1)*(n3-p3+1);
    return true;
  }
  else if (shareFE == 'F')
    return true;

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
  if (p1 <  1 || p2 <  1 || p3 <  1) return false;
  if (p1 > n1 || p2 > n2 || p3 > n3) return false;

  myMLGE.resize((n1-p1+1)*(n2-p2+1)*(n3-p3+1),0);
  myMLGN.resize(n1*n2*n3);
  myMNPC.resize(myMLGE.size());
  myNodeInd.resize(myMLGN.size());

  nnod = nel = 0;
  for (int i3 = 1; i3 <= n3; i3++)
    for (int i2 = 1; i2 <= n2; i2++)
      for (int i1 = 1; i1 <= n1; i1++)
      {
        myNodeInd[nnod].I = i1-1;
        myNodeInd[nnod].J = i2-1;
        myNodeInd[nnod].K = i3-1;
        if (i1 >= p1 && i2 >= p2 && i3 >= p3)
        {
          if (svol->knotSpan(0,i1-1) > 0.0)
            if (svol->knotSpan(1,i2-1) > 0.0)
              if (svol->knotSpan(2,i3-1) > 0.0)
              {
                myMLGE[nel] = ++gEl; // global element number over all patches
                myMNPC[nel].resize(p1*p2*p3,0);

                int lnod = 0;
                for (int j3 = p3-1; j3 >= 0; j3--)
                  for (int j2 = p2-1; j2 >= 0; j2--)
                    for (int j1 = p1-1; j1 >= 0; j1--)
                      myMNPC[nel][lnod++] = nnod - n1*n2*j3 - n1*j2 - j1;
              }

          nel++;
        }
        myMLGN[nnod++] = ++gNod; // global node number over all patches
      }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
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
  if (shareFE == 'F') return true;

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
	      myMLGN[inod] = nodes.ibnod[0];
	    else if (i == n1)
	      myMLGN[inod] = nodes.ibnod[1];
	    else
	      myMLGN[inod] = nodes.edges[0].next();
	  }
	  else if (j == n2)
	  {
	    if (i == 1)
	      myMLGN[inod] = nodes.ibnod[2];
	    else if (i == n1)
	      myMLGN[inod] = nodes.ibnod[3];
	    else
	      myMLGN[inod] = nodes.edges[1].next();
	  }
	  else
	  {
	    if (i == 1)
	      myMLGN[inod] = nodes.edges[4].next();
	    else if (i == n1)
	      myMLGN[inod] = nodes.edges[5].next();
	    else
	      myMLGN[inod] = nodes.faces[4].next();
	  }
	}
	else if (k == n3)
	{
	  if (j == 1)
	  {
	    if (i == 1)
	      myMLGN[inod] = nodes.ibnod[4];
	    else if (i == n1)
	      myMLGN[inod] = nodes.ibnod[5];
	    else
	      myMLGN[inod] = nodes.edges[2].next();
	  }
	  else if (j == n2)
	  {
	    if (i == 1)
	      myMLGN[inod] = nodes.ibnod[6];
	    else if (i == n1)
	      myMLGN[inod] = nodes.ibnod[7];
	    else
	      myMLGN[inod] = nodes.edges[3].next();
	  }
	  else
	  {
	    if (i == 1)
	      myMLGN[inod] = nodes.edges[6].next();
	    else if (i == n1)
	      myMLGN[inod] = nodes.edges[7].next();
	    else
	      myMLGN[inod] = nodes.faces[5].next();
	  }
	}
	else
	{
	  if (j == 1)
	  {
	    if (i == 1)
	      myMLGN[inod] = nodes.edges[8].next();
	    else if (i == n1)
	      myMLGN[inod] = nodes.edges[9].next();
	    else
	      myMLGN[inod] = nodes.faces[2].next();
	  }
	  else if (j == n2)
	  {
	    if (i == 1)
	      myMLGN[inod] = nodes.edges[10].next();
	    else if (i == n1)
	      myMLGN[inod] = nodes.edges[11].next();
	    else
	      myMLGN[inod] = nodes.faces[3].next();
	  }
	  else
	  {
	    if (i == 1)
	      myMLGN[inod] = nodes.faces[0].next();
	    else if (i == n1)
	      myMLGN[inod] = nodes.faces[1].next();
	    else
	      myMLGN[inod] = nodes.next();
	  }
	}

#if SP_DEBUG > 2
  if (basis > 0) std::cout <<"\nBasis "<< basis <<":";
  for (int i = inod-n1*n2*n3; i < inod; i++)
  {
    std::cout <<"\nNode "<< i+1 <<"\t: ";
    if (!nodeInd.empty())
      std::cout << nodeInd[i].I <<" "<< nodeInd[i].J <<" "<< nodeInd[i].K;
    std::cout <<"\tglobal no. "<< MLGN[i];
  }
  std::cout << std::endl;
#endif
  return true;
}


bool ASMs3D::connectPatch (int face, ASM3D& neighbor, int nface,
                           int norient, int, bool coordCheck, int thick)
{
  ASMs3D* neighS = dynamic_cast<ASMs3D*>(&neighbor);
  if (!neighS)
    return false;

  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  if (neighS->swapW && nface > 4) // Account for swapped parameter direction
    nface = 11-nface;

  if (!this->connectBasis(face,*neighS,nface,norient,1,0,0,coordCheck,thick))
    return false;

  this->addNeighbor(neighS);
  return true;
}


bool ASMs3D::connectBasis (int face, ASMs3D& neighbor, int nface, int norient,
                           int basis, int slave, int master,
                           bool coordCheck, int thick)
{
  if (this->isShared() && neighbor.isShared())
    return true;
  else if (this->isShared() || neighbor.isShared())
  {
    std::cerr <<" *** ASMs3D::connectPatch: Logic error, cannot"
	      <<" connect a shared patch with an unshared one"<< std::endl;
    return false;
  }
  else if (norient < 0 || norient > 7)
  {
    std::cerr <<" *** ASMs3D::connectPatch: Orientation flag "<< norient
              <<" is out of range [0,7]"<< std::endl;
    return false;
  }

  int m1, m2;
  if (!this->getFaceSize(m1,m2,basis,face))
    return false;

  int n1, n2;
  if (!neighbor.getFaceSize(n1,n2,basis,nface))
    return false;
  else if (norient > 3)
    std::swap(n1,n2);

  // Set up the slave node numbers for this volume patch
  IntVec slaveNodes;
  this->getBoundaryNodes(face, slaveNodes, basis, thick, 0, true);
  for (int& it : slaveNodes)
    it += slave;

  // Set up the master node numbers for the neighboring volume patch
  IntVec masterNodes;
  neighbor.getBoundaryNodes(nface, masterNodes, basis, thick, 0, true);
  for (int& it : masterNodes)
    it += master;

  if (n1 == 0 || n2 == 0 || n1 != m1 || n2 != m2)
  {
    std::cerr <<" *** ASMs3D::connectPatch: Non-matching faces, sizes "
              << n1 <<","<< n2 <<" and "<< m1 <<","<< m2 << std::endl;
    return false;
  }

  if (norient > 3)
    std::swap(n1,n2);

  const double xtol = 1.0e-4;
  int node = 1, failures = 0;
  for (int j = 0; j < n2; j++)
    for (int i = 0; i < n1; i++, node++)
    {
      int k, l;
      switch (norient)
        {
        case 0: k =    i  ; l =    j  ; break;
        case 1: k =    i  ; l = n2-j-1; break;
        case 2: k = n1-i-1; l =    j  ; break;
        case 3: k = n1-i-1; l = n2-j-1; break;
        case 4: k =    j  ; l =    i  ; break;
        case 5: k =    j  ; l = m2-i-1; break;
        case 6: k = m1-j-1; l =    i  ; break;
        case 7: k = m1-j-1; l = m2-i-1; break;
        }

      for (int t = 0; t < thick; t++)
      {
        int snode = slaveNodes[(l*m1+k)*thick+t];
        int mnode = masterNodes[(node-1)*thick+t];
        if (!coordCheck)
          ASMbase::collapseNodes(neighbor,mnode,*this,snode);
        else if (neighbor.getCoord(mnode).equal(this->getCoord(snode),xtol))
          ASMbase::collapseNodes(neighbor,mnode,*this,snode);
        else
        {
          std::cerr <<" *** ASMs3D::connectPatch: Non-matching nodes "
                    << mnode <<": "<< neighbor.getCoord(mnode) <<"\n"
                    << std::string(42,' ') <<"and "
                    << snode <<": "<< this->getCoord(snode) << std::endl;
          failures++;
        }
      }
    }

  return failures == 0;
}


void ASMs3D::closeBoundaries (int dir, int basis, int master)
{
  if (basis < 1) basis = 1;

  int n1, n2, n3;
  if (!this->getSize(n1,n2,n3,basis)) return;

  int dirs = utl::getDirs(this->getNoFields(basis));

  switch (dir)
    {
    case 1: // Faces are closed in I-direction
      for (int i3 = 1; i3 <= n3; i3++)
        for (int i2 = 1; i2 <= n2; i2++, master += n1)
          this->makePeriodic(master,master+n1-1,dirs);
      threadGroupsVol.stripDir = ThreadGroups::U;
      break;

    case 2: // Faces are closed in J-direction
      for (int i3 = 1; i3 <= n3; i3++, master += n1*(n2-1))
        for (int i1 = 1; i1 <= n1; i1++, master++)
          this->makePeriodic(master,master+n1*(n2-1),dirs);
      threadGroupsVol.stripDir = ThreadGroups::V;
      break;

    case 3: // Faces are closed in K-direction
      for (int i2 = 1; i2 <= n2; i2++)
        for (int i1 = 1; i1 <= n1; i1++, master++)
          this->makePeriodic(master,master+n1*n2*(n3-1),dirs);
      threadGroupsVol.stripDir = ThreadGroups::W;
      break;
    }
}


bool ASMs3D::collapseFace (int face, int edge, int basis)
{
  if (basis < 1) basis = 1;

  int n1, n2, n3;
  int node = this->findStartNode(n1,n2,n3,basis);
  if (node < 1) return false;

  if (swapW && face > 4) // Account for swapped parameter direction
    face = 11-face;

  // Lambda function to verify co-location of nodes and collapse them
  auto&& collapse = [this](int master, int slave)
  {
    const double xtol = 1.0e-4;

    if (master == slave)
      return true;
    else if (this->getCoord(slave).equal(this->getCoord(master),xtol))
    {
      ASMbase::collapseNodes(*this,master,*this,slave);
      return true;
    }

    std::cerr <<" *** ASMs3D::collapseFace: Not all nodes on the face"
              <<" are co-located\n"<< std::string(27,' ')
              << master <<": " << this->getCoord(master) << std::endl
              << std::string(23,' ') <<"and "
              << slave <<": "<< this->getCoord(slave) << std::endl;
    return false;
  };

  int iedge, master;
  switch (face)
    {
    case 2: // Right face (positive I-direction)
      node += n1-1;
    case 1: // Left face (negative I-direction)
      iedge = (edge-face-2)/2;
      master = node;
      if (iedge == 2)
        master += n1*n2*(n3-1);
      else if (iedge == 4)
        master += n1*(n2-1);
      for (int i3 = 1; i3 <= n3; i3++)
        for (int i2 = 1; i2 <= n2; i2++, node += n1)
        {
          if (!collapse(master,node))
            return false;
          switch (iedge) {
          case 1:
          case 2:
            master += n1;
            if (i2 == n2)
              master -= n1*n2;
            break;
          case 3:
          case 4:
            if (i2 == n2)
              master += n1*n2;
            break;
          case -1:
          case -2:
            break; // edge=0, collapsing to a vertex
          default:
            std::cerr <<" *** ASMs3D::collapseFace: Invalid face/edge "
                      << face <<"/"<< edge << std::endl;
            return false;
          }
        }

      threadGroupsVol.stripDir = ThreadGroups::U;
      break;

    case 4: // Back face (positive J-direction)
      node += n1*(n2-1);
    case 3: // Front face (negative J-direction)
      iedge = edge < 9 ? (edge-face+4)/2 : edge-2*face;
      master = node;
      if (iedge == 2)
        master += n1*n2*(n3-1);
      else if (iedge == 4)
        master += n1-1;
      for (int i3 = 1; i3 <= n3; i3++, node += n1*(n2-1))
        for (int i1 = 1; i1 <= n1; i1++, node++)
        {
          if (!collapse(master,node))
            return false;
          switch (iedge) {
          case 1:
          case 2:
            master++;
            if (i1 == n1)
              master -= n1;
            break;
          case 3:
          case 4:
            if (i1 == n1)
              master += n1*n2;
            break;
          case 0:
            break; // edge=0, collapsing to a vertex
          default:
            std::cerr <<" *** ASMs3D::collapseFace: Invalid face/edge "
                      << face <<"/"<< edge << std::endl;
            return false;
          }
        }

      threadGroupsVol.stripDir = ThreadGroups::V;
      break;

    case 6: // Top face (positive K-direction)
      node += n1*n2*(n3-1);
    case 5: // Bottom face (negative K-direction)
      iedge = edge-2*face+10;
      master = node;
      for (int i2 = 1; i2 <= n2; i2++)
        for (int i1 = 1; i1 <= n1; i1++, node++)
        {
          if (!collapse(master,node))
            return false;
          switch (iedge) {
          case 1:
          case 2:
            master++;
            if (i1 == n1)
              master -= n1;
            break;
          case 5:
          case 6:
            if (i1 == n1)
              master += n1;
            break;
          case 0:
          case -2:
            break; // edge=0, collapsing to a vertex
          default:
            std::cerr <<" *** ASMs3D::collapseFace: Invalid face/edge "
                      << face <<"/"<< edge << std::endl;
            return false;
          }
        }

      threadGroupsVol.stripDir = ThreadGroups::W;
      break;

    default:
      return false;
    }

  return true;
}


int ASMs3D::findStartNode (int& n1, int& n2, int& n3, char basis) const
{
  int node = 1;
  for (char i = 1; i <= basis; i++)
    if (!this->getSize(n1,n2,n3,i))
      return 0;
    else if (i < basis)
      node += n1*n2*n3;

  return node;
}


/*!
  A negative \a code value implies direct evaluation of the Dirichlet condition
  function at the control point. Positive \a code implies projection onto the
  spline basis representing the boundary surface (needed for curved faces and/or
  non-constant functions).
*/

void ASMs3D::constrainFace (int dir, bool open, int dof,
                            int code, char basis)
{
  if (basis < 1) basis = 1;

  int n1, n2, n3;
  int node = this->findStartNode(n1,n2,n3,basis);
  if (node < 1) return;

  if (swapW && abs(dir) == 3) // Account for swapped parameter direction
    dir = -dir;

  int bcode = code;
  if (code > 0) // Dirichlet projection will be performed
    dirich.push_back(DirichletFace(this->getBoundary(dir,basis),dof,code));
  else if (code < 0)
    bcode = -code;

  switch (dir)
    {
    case  1: // Right face (positive I-direction)
      node += n1-1;
    case -1: // Left face (negative I-direction)
      for (int i3 = 1; i3 <= n3; i3++)
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
	  if (open && (i2 == 1 || i2 == n2 || i3 == 1 || i3 == n3))
	    continue; // skip all edge nodes if an open boundary is requested
	  else if ((i2 == 1 || i2 == n2) && (i3 == 1 || i3 == n3))
	    this->prescribe(node,dof,bcode); // corner node
	  else
	  {
	    // If the Dirichlet condition is to be projected, add this node to
	    // the set of nodes to receive prescribed value from the projection
	    // **unless this node already has a homogeneous constraint**
	    if (this->prescribe(node,dof,-code) == 0 && code > 0)
	      dirich.back().nodes.push_back(std::make_pair(n2*(i3-1)+i2,node));
	  }
      break;

    case  2: // Back face (positive J-direction)
      node += n1*(n2-1);
    case -2: // Front face (negative J-direction)
      for (int i3 = 1; i3 <= n3; i3++, node += n1*(n2-1))
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  if (open && (i1 == 1 || i1 == n1 || i3 == 1 || i3 == n3))
	    continue; // skip all edge nodes if an open boundary is requested
	  else if ((i1 == 1 || i1 == n1) && (i3 == 1 || i3 == n3))
	    this->prescribe(node,dof,bcode); // corner node
	  else
	  {
	    // If the Dirichlet condition is to be projected, add this node to
	    // the set of nodes to receive prescribed value from the projection
	    // **unless this node already has a homogeneous constraint**
	    if (this->prescribe(node,dof,-code) == 0 && code > 0)
	      dirich.back().nodes.push_back(std::make_pair(n1*(i3-1)+i1,node));
	  }
      break;

    case  3: // Top face (positive K-direction)
      node += n1*n2*(n3-1);
    case -3: // Bottom face (negative K-direction)
      for (int i2 = 1; i2 <= n2; i2++)
	for (int i1 = 1; i1 <= n1; i1++, node++)
	  if (open && (i1 == 1 || i1 == n1 || i2 == 1 || i2 == n2))
	    continue; // skip all edge nodes if an open boundary is requested
	  else if ((i1 == 1 || i1 == n1) && (i2 == 1 || i2 == n2))
	    this->prescribe(node,dof,bcode); // corner node
	  else
	  {
	    // If the Dirichlet condition is to be projected, add this node to
	    // the set of nodes to receive prescribed value from the projection
	    // **unless this node already has a homogeneous constraint**
	    if (this->prescribe(node,dof,-code) == 0 && code > 0)
	      dirich.back().nodes.push_back(std::make_pair(n1*(i2-1)+i1,node));
	  }
      break;
    }
}


/*!
  The local coordinate systems in which the constraints are applied,
  are constructed from the tangent directions of the boundary surface,
  evaluated at the Greville points. The local Z-direction is then the
  outward-directed normal, computed as the cross product of the two tangents.
  If \a project is \e true, the normal vector is projected onto the surface
  basis of the face in order to obtain corresponding control point values.
  Otherwise, it is used directly.

  A negative \a code value implies direct evaluation of the Dirichlet condition
  function at the control point. Positive \a code implies projection onto the
  spline basis representing the boundary surface (needed for curved faces and/or
  non-constant functions).
*/

size_t ASMs3D::constrainFaceLocal (int dir, bool open, int dof, int code,
				   bool project, char T1)
{
  if (shareFE == 'F')
  {
    std::cerr <<"\n *** ASMs3D::constrainFaceLocal: Logic error, can not have"
             <<" constraints in local CSs for shared patches."<< std::endl;
    return 0;
  }

  int t1 = abs(dir)%3; // first tangent direction [0,2]
  int t2 = (1+t1)%3;   // second tangent direction [0,2]
  if (swapW && t1 == 0) dir = -dir; // Account for swapped parameter direction
  if (t1 == 2) std::swap(t1,t2);

  // Get parameter values of the Greville points on the face
  RealArray upar, vpar;
  if (!this->getGrevilleParameters(upar,t1) ||
      !this->getGrevilleParameters(vpar,t2))
    return 0;

  // Find the surface representing the face geometry (for tangent evaluation)
  Go::SplineSurface* face = this->getBoundary(dir);
  if (!face) return 0;

  // We need to add extra nodes, check that the global node counter is good.
  // If not, we cannot do anything here
  if (gNod < *std::max_element(MLGN.begin(),MLGN.end()))
  {
    std::cerr <<"\n *** ASMs3D::constrainFaceLocal: Logic error, gNod = "<< gNod
              <<" is too small!"<< std::endl;
    return 0;
  }
  if (nf < 3)
  {
    std::cerr <<"\n *** ASMs3D::constrainFaceLocal: Not for scalar problems!"
              << std::endl;
    return 0;
  }

  // Loop over the Greville points along over the face
  size_t i, j, k = 0;
  unsigned char c, d;
  std::vector<Go::Point> pts(3);
  RealArray gdata(3*upar.size()*vpar.size());
  for (j = 0; j < vpar.size(); j++)
    for (i = 0; i < upar.size(); i++)
    {
      // Compute the outward-directed face normal at this point.
      // That will be the local z-axis of the local coordinate system.
      face->point(pts,upar[i],vpar[j],1);
      Vec3 Zaxis(SplineUtils::toVec3(pts[1]),SplineUtils::toVec3(pts[2]));
      gdata[k++] = dir > 0 ? Zaxis.x : -Zaxis.x;
      gdata[k++] = dir > 0 ? Zaxis.y : -Zaxis.y;
      gdata[k++] = dir > 0 ? Zaxis.z : -Zaxis.z;
    }

  Go::SplineSurface* locs = nullptr;
  if (project)
  {
    // Project the Greville point values onto the spline basis
    // to obtain the corresponding control point values

    RealArray weights;
    if (face->rational())
      face->getWeights(weights);

    locs = Go::SurfaceInterpolator::regularInterpolation(face->basis(0),
							 face->basis(1),
							 upar,vpar,gdata,3,
							 face->rational(),
							 weights);
  }

  // Find start index and increment of the (slave) node with the global DOFs
  int n1, n2, n3, incI = 1, incJ = 0, iSnod = 0;
  this->getSize(n1,n2,n3,1);
  switch (dir)
    {
    case  1: // Right face (positive I-direction)
      iSnod = n1-1;
    case -1: // Left face (negative I-direction)
      incI = n1;
      break;

    case  2: // Back face (positive J-direction)
      iSnod = n1*(n2-1);
    case -2: // Front face (negative J-direction)
      incJ = n1*(n2-1);
      break;

    case  3: // Top face (positive K-direction)
      iSnod = n1*n2*(n3-1);
    case -3: // Bottom face (negative K-direction)
      break;
    }

  int bcode = code;
  if (code > 0) // Dirichlet projection will be performed
    dirich.push_back(DirichletFace(this->getBoundary(dir),dof,code));
  else if (code < 0)
    bcode = -code;

  size_t nxNode = 0;
  RealArray::const_iterator it = locs ? locs->coefs_begin() : gdata.begin();
  for (j = k = 0; j < vpar.size(); j++, iSnod += incJ)
    for (i = 0; i < upar.size(); i++, iSnod += incI, it += 3)
    {
      // Skip the edge points if this should be regarded an open boundary
      if (open && (i == 0 || i+1 == upar.size())) continue;
      if (open && (j == 0 || j+1 == vpar.size())) continue;
      // Check if this node already has been constrained or fixed
      if (this->isFixed(MLGN[iSnod],dof)) continue;

      // We need an extra node representing the local DOFs at this point
      int iMnod = myMLGN.size();
      std::map<int,int>::const_iterator xit = xNode.end();
      // Create an extra node for the local DOFs. The new node, for which
      // the Dirichlet boundary conditions will be defined, then inherits
      // the global node number of the original node. The original node, which
      // do not enter the equation system, receives a new global node number.
      if (i > 0 && i+1 < upar.size() && j > 0 && j+1 < vpar.size())
        // This is not an edge node
        myMLGN.push_back(++gNod);
      else if ((xit = xNode.find(MLGN[iSnod])) != xNode.end())
        // This is an edge node already processed by another patch
        myMLGN.push_back(xit->second);
      else
      {
        // This is an edge node, store its original-to-extra node number
        // mapping in ASMstruct::xNode
        myMLGN.push_back(++gNod);
        xNode[MLGN[iSnod]] = gNod;
      }
      std::swap(myMLGN[iMnod],myMLGN[iSnod]);

      xnMap[1+iMnod] = 1+iSnod; // Store nodal connection needed by getCoord
      nxMap[1+iSnod] = 1+iMnod; // Store nodal connection needed by getNodeID
      if (xit != xNode.end()) continue; // This node has already been processed

      ++nxNode; // Increment number of added nodes

      // Global node numbers of the nodes to be coupled
      int masterNode = MLGN[iMnod];
      int slaveNode  = MLGN[iSnod];

      // Add Dirichlet condition on the local DOF(s) of the added node
      if ((i == 0 || i+1 == upar.size()) && (j == 0 || j+1 == vpar.size()))
        this->prescribe(1+iMnod,dof,bcode); // corner node
      else
      {
        this->prescribe(1+iMnod,dof,-code);
        if (code > 0)
          dirich.back().nodes.push_back(std::make_pair(1+k,1+iMnod));
      }

      // Local-to-global transformation matrix at this point, created by
      // projecting the global X- or Y-axis onto the tangent plane defined by
      // the normal vector (see the Tensor::Tensor(const Vec3&) constructor).
      // If T1 is defined, it is used to determine the first tangent direction.
      Tensor Tlg(3);
      if (T1 < 'x' || T1 > 'z')
        Tlg = Tensor(Vec3(it[0],it[1],it[2]));
      else
      {
        Vec3 v1;
        v1[T1-'x'] = 1.0;
        Tlg = Tensor(Vec3(it[0],it[1],it[2]),v1,true);
      }

      // Now establish constraint equations relating the global and local DOFs.
      // We here assume there are (at least) 3 unknowns per node,
      // and only the first 3 DOFs are subjected to transformation.
      for (d = 1; d <= nf; d++)
      {
        MPC* cons = new MPC(slaveNode,d);
        if (this->addMPC(cons) && cons)
        {
          if (d > 3)
            cons->addMaster(masterNode,d);
          else for (c = 1; c <= 3; c++)
            if (!this->isFixed(masterNode,c))
              cons->addMaster(masterNode,c,Tlg(d,c));
#if SP_DEBUG > 1
          std::cout <<"Added constraint: "<< *cons;
#endif
        }
      }
    }

  if (locs) delete locs;
  return nxNode; // Number of added nodes
}


void ASMs3D::constrainEdge (int lEdge, bool open, int dof,
                            int code, char basis)
{
  IntVec edgeNodes;
  this->getBoundary1Nodes(lEdge,edgeNodes,basis,0,true,open);
  for (int node : edgeNodes)
    this->prescribe(node,dof,code);
}


void ASMs3D::constrainLine (int fdir, int ldir, double xi, int dof,
                            int code, char basis)
{
  if (xi < 0.0 || xi > 1.0) return;
  if (basis < 1) basis = 1;

  int n1, n2, n3;
  int node = this->findStartNode(n1,n2,n3,basis);
  if (node < 1) return;

  if (swapW && abs(fdir) == 3) // Account for swapped parameter direction
    fdir = -fdir;

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


void ASMs3D::constrainCorner (int I, int J, int K, int dof,
                              int code, char basis)
{
  if (basis < 1) basis = 1;

  int node = this->getCorner(I,J,K,basis);
  if (node > 0)
    this->prescribe(node,dof,code);
}


void ASMs3D::constrainNode (double xi, double eta, double zeta, int dof,
                            int code)
{
  if (xi   < 0.0 || xi   > 1.0) return;
  if (eta  < 0.0 || eta  > 1.0) return;
  if (zeta < 0.0 || zeta > 1.0) return;

  int n1 = 0, n2 = 0, n3 = 0;
  if (!this->getSize(n1,n2,n3,1))
    return;

  if (swapW) // Account for swapped parameter direction
    zeta = 1.0-zeta;

  int node = 1;
  if (xi   > 0.0) node += int(0.5+(n1-1)*xi);
  if (eta  > 0.0) node += n1*int(0.5+(n2-1)*eta);
  if (zeta > 0.0) node += n1*n2*int(0.5+(n3-1)*zeta);

  this->prescribe(node,dof,code);
}


void ASMs3D::setNodeNumbers (const IntVec& nodes)
{
  this->ASMbase::setNodeNumbers(nodes);
  if (!swapW) return;

  // Account for swapped parameter direction
  const int n1 = svol->numCoefs(0)*svol->numCoefs(1);
  const int n2 = svol->numCoefs(2);
  for (int j = 0; j < n2/2; j++)
    for (int i = 0; i < n1; i++)
      std::swap(myMLGN[i+n1*j],myMLGN[i+n1*(n2-j-1)]);
}


/*!
  This method projects the function describing the in-homogeneous Dirichlet
  boundary condition onto the spline basis defining the boundary surface,
  in order to find the control point values which are used as the prescribed
  values of the boundary DOFs.
*/

bool ASMs3D::updateDirichlet (const std::map<int,RealFunc*>& func,
			      const std::map<int,VecFunc*>& vfunc, double time,
                              const std::map<int,int>* g2l)
{
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;

  for (size_t i = 0; i < dirich.size(); i++)
  {
    // Project the function onto the spline surface basis
    Go::SplineSurface* dsurf = nullptr;
    if ((fit = func.find(dirich[i].code)) != func.end())
      dsurf = SplineUtils::project(dirich[i].surf,*fit->second,1,time);
    else if ((vfit = vfunc.find(dirich[i].code)) != vfunc.end())
      dsurf = SplineUtils::project(dirich[i].surf,*vfit->second,
                                   vfit->second->dim(),time);
    else
    {
      std::cerr <<" *** ASMs3D::updateDirichlet: Code "<< dirich[i].code
		<<" is not associated with any function."<< std::endl;
      return false;
    }
    if (!dsurf)
    {
      std::cerr <<" *** ASMs3D::updateDirichlet: Projection failure."
		<< std::endl;
      return false;
    }

    // Loop over the (interior) nodes (control points) of this boundary surface
    for (const Ipair& node : dirich[i].nodes)
      for (int dofs = dirich[i].dof; dofs > 0; dofs /= 10)
      {
        int dof = dofs%10;
        // Find the constraint equation for current (node, local DOF)
        MPC pDOF(MLGN[node.second-1],dof);
        MPCIter mit = mpcs.find(&pDOF);
        if (mit == mpcs.end()) continue; // probably a deleted constraint

        // Find the control point index for this (node,local DOF) in dsurf
        RealArray::const_iterator cit = dsurf->coefs_begin();
        if (dsurf->dimension() > 1) // A vector field is specified
          cit += (node.first-1)*dsurf->dimension() + (dof-1);
        else // A scalar field is specified at this DOF
          cit += (node.first-1);

        // Now update the prescribed value in the constraint equation
        (*mit)->setSlaveCoeff(*cit);
#if SP_DEBUG > 1
        std::cout <<"Updated constraint: "<< **mit;
#endif
      }
    delete dsurf;
  }

  // The parent class method takes care of the corner nodes with direct
  // evaluation of the Dirichlet functions (since they are interpolatory)
  return this->ASMbase::updateDirichlet(func,vfunc,time,g2l);
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

  int inod1 = MNPC[iel-1][this->getLastItgElmNode()];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs3D::getParametricVolume: Node index "<< inod1
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
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

  int inod1 = MNPC[iel-1][this->getLastItgElmNode()];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs3D::getParametricArea: Node index "<< inod1
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
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
  if (inod >= nnod)
  {
    std::cerr <<" *** ASMs3D::coeffInd: Node index "<< inod
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
    return -1;
  }
#endif

  const int ni = nodeInd[inod].I;
  const int nj = nodeInd[inod].J;
  const int nk = nodeInd[inod].K;
  return (nk*svol->numCoefs(1) + nj)*svol->numCoefs(0) + ni;
}


Vec3 ASMs3D::getCoord (size_t inod) const
{
  if (inod > nnod && inod <= MLGN.size())
  {
    // This is a node added due to constraints in local directions.
    // Find the corresponding original node (see constrainEdgeLocal)
    std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
    if (it != xnMap.end()) inod = it->second;
  }
  if (inod == 0) return Vec3();

  std::map<size_t,XYZ>::const_iterator it = myRmaster.find(inod);
  if (it != myRmaster.end()) return Vec3(it->second.data());

  int ip = this->coeffInd(inod-1)*svol->dimension();
  if (ip < 0) return Vec3();

  RealArray::const_iterator cit = svol->coefs_begin() + ip;
  return Vec3(*cit,*(cit+1),*(cit+2));
}


bool ASMs3D::getElementCoordinates (Matrix& X, int iel, bool forceItg) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs3D::getElementCoordinates: Element index "<< iel
              <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  const Go::SplineVolume* geo = forceItg ? svol.get() : this->getBasis(ASM::GEOMETRY_BASIS);

  int lnod0 = this->getFirstItgElmNode();
  if (geo != svol.get()) {
    const IJK& nIdx = nodeInd[MNPC[iel-1][lnod0]];
    double u = *(svol->basis(0).begin() + nIdx.I + svol->order(0) - 1);
    double v = *(svol->basis(1).begin() + nIdx.J + svol->order(1) - 1);
    double w = *(svol->basis(2).begin() + nIdx.K + svol->order(2) - 1);
    return this->getElementCoordinatesPrm(X,u,v,w);
  }

  X.resize(3,geo->order(0)*geo->order(1)*geo->order(2));

  RealArray::const_iterator cit = geo->coefs_begin();
  for (size_t n = 0; n < X.cols(); n++)
  {
    int ip = this->coeffInd(MNPC[iel-1][n + lnod0])*svol->dimension();
    if (ip < 0) return false;

    for (size_t i = 0; i < 3; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


bool ASMs3D::getElementCoordinatesPrm (Matrix& X, double u,
                                       double v, double w) const
{
  const Go::SplineVolume* geo = this->getBasis(ASM::GEOMETRY_BASIS);

  if (u < geo->startparam(0) || u > geo->endparam(0) ||
      v < geo->startparam(1) || v > geo->endparam(1) ||
      w < geo->startparam(2) || w > geo->endparam(2))
    return false;

  int ni, nj, nk;
#pragma omp critical
  {
    ni = geo->basis(0).knotInterval(u) - geo->order(0) + 1;
    nj = geo->basis(1).knotInterval(v) - geo->order(1) + 1;
    nk = geo->basis(2).knotInterval(w) - geo->order(2) + 1;
  }

  X.resize(3,geo->order(0)*geo->order(1)*geo->order(2));
  RealArray::const_iterator cit = geo->coefs_begin();
  for (size_t n = 0; n < X.cols(); n++)
  {
    const int iu = n % geo->order(0);
    const int iv = (n / geo->order(0)) % geo->order(1);
    const int iw = n / (geo->order(0) * geo->order(1));
    const int ip = (ni + iu + ((nk + iw)*geo->numCoefs(1) + (nj + iv))*geo->numCoefs(0))*geo->dimension();

    for (size_t i = 0; i < 3; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element containing parameters ("
            << u <<"," << v << "," << w <<"):" << X << std::endl;
#endif
  return true;
}


void ASMs3D::getNodalCoordinates (Matrix& X, bool geo) const
{
  const Go::SplineVolume* spline = geo ? this->getBasis(ASM::GEOMETRY_BASIS) : svol.get();
  const int n1 = spline->numCoefs(0);
  const int n2 = spline->numCoefs(1);
  const int n3 = spline->numCoefs(2);
  X.resize(3,n1*n2*n3);

  RealArray::const_iterator cit = spline->coefs_begin();
  size_t inod = 1;
  for (int i3 = 0; i3 < n3; i3++)
    for (int i2 = 0; i2 < n2; i2++)
      for (int i1 = 0; i1 < n1; i1++, inod++)
      {
        int ip = ((i3*n2 + i2)*n1 + i1)*spline->dimension();
        for (size_t i = 0; i < 3; i++)
          X(i+1,inod) = *(cit+(ip+i));
      }
}


bool ASMs3D::updateCoords (const Vector& displ)
{
  if (!svol) return true; // silently ignore empty patches
  if (shareFE) return true;

  size_t nno = svol->numCoefs(0)*svol->numCoefs(1)*svol->numCoefs(2);
  if (displ.size() != 3*nno && displ.size() != 3*MLGN.size())
  {
    std::cerr <<" *** ASMs3D::updateCoords: Invalid dimension "
              << displ.size() <<" on displacement vector, should be ";
    if (nno != MLGN.size())
      std::cerr <<"either "<< 3*MLGN.size() <<" or ";
    std::cerr << 3*nno << std::endl;
    return false;
  }

  svol->deform(displ,3);
  return true;
}


Vec3 ASMs3D::getElementCenter (int iel) const
{
  std::cerr <<" *** ASMs3D::getElementCenter("<< iel
            <<"): Not implemented yet."<< std::endl;
  return Vec3(); //TODO: Maybe later, if needed
}


void ASMs3D::getBoundaryNodes (int lIndex, IntVec& nodes,
                               int basis, int thick, int, bool local) const
{
  if (basis < 1) basis = 1;

#if SP_DEBUG > 1
  size_t last = nodes.size();
  std::cout <<"Boundary nodes in patch "<< idx+1 <<" face "<< lIndex <<":";
#endif

  int n1, n2, n3, inod, t;
  int node = this->findStartNode(n1,n2,n3,basis);
  if (node < 1)
    return;
  else if (local)
    node = 1;

  if (swapW && lIndex > 4) // Account for swapped parameter direction
    lIndex = 11-lIndex;

  switch (lIndex)
  {
    case 2: // Right face (positive I-direction)
      node += n1-thick;
    case 1: // Left face (negative I-direction)
      for (int i3 = 1; i3 <= n3; i3++)
	for (int i2 = 1; i2 <= n2; i2++, node += n1)
          for (t = 0, inod = node; t < thick; t++, inod++)
            nodes.push_back(local ? inod : this->getNodeID(inod));
      break;

    case 4: // Back face (positive J-direction)
      node += n1*(n2-thick);
    case 3: // Front face (negative J-direction)
      for (int i3 = 1; i3 <= n3; i3++, node += n1*(n2-1))
	for (int i1 = 1; i1 <= n1; i1++, node++)
          for (t = 0, inod = node; t < thick; t++, inod += n1)
            nodes.push_back(local ? inod : this->getNodeID(inod));
      break;

    case 6: // Top face (positive K-direction)
      node += n1*n2*(n3-thick);
    case 5: // Bottom face (negative K-direction)
      for (int i2 = 1; i2 <= n2; i2++)
	for (int i1 = 1; i1 <= n1; i1++, node++)
          for (t = 0, inod = node; t < thick; t++, inod += n1*n2)
            nodes.push_back(local ? inod : this->getNodeID(inod));
      break;
  }

#if SP_DEBUG > 1
  if (nodes.size() == last)
    std::cout <<" (none)";
  else for (size_t i = last; i < nodes.size(); i++)
    std::cout <<" "<< nodes[i];
  std::cout << std::endl;
#endif
}


void ASMs3D::getBoundary1Nodes (int lEdge, IntVec& nodes,
                                int basis, int, bool local, bool open) const
{
  if (basis < 1) basis = 1;

#if SP_DEBUG > 1
  size_t last = nodes.size();
  std::cout <<"Boundary nodes in patch "<< idx+1 <<" edge "<< lEdge <<":";
#endif

  int n1, n2, n3;
  int node = this->findStartNode(n1,n2,n3,basis);
  if (node < 1)
    return;
  else if (local)
    node = 1;

  if (swapW && lEdge <= 8) // Account for swapped parameter direction
    lEdge += (lEdge-1)%4 < 2 ? 2 : -2;

  std::array<int,2> range;
  if (lEdge > 8)
    range = { n1*n2, n3 };
  else if (lEdge > 4)
    range = { n1, n2 };
  else
    range = { 1, n1 };

  switch (lEdge)
    {
    case  6:
    case 10:
      node += n1 - 1;
      break;
    case  2:
    case 11:
      node += n1*(n2-1);
      break;
    case 12:
      node += n1*n2 - 1;
      break;
    case  3:
    case  7:
      node += n1*n2*(n3-1);
      break;
    case  8:
      node += n1*(n2*(n3-1) + 1) - 1;
      break;
    case  4:
      node += n1*(n2*n3-1);
      break;
    }

  // Skip the first and last node if we are requesting an open boundary
  for (int i = 1; i <= range[1]; i++, node += range[0])
    if (!open || (i > 1 && i < range[1]))
      nodes.push_back(local ? node : this->getNodeID(node));

#if SP_DEBUG > 1
  if (nodes.size() == last)
    std::cout <<" (none)";
  else for (size_t i = last; i < nodes.size(); i++)
    std::cout <<" "<< nodes[i];
  std::cout << std::endl;
#endif
}


bool ASMs3D::getOrder (int& p1, int& p2, int& p3) const
{
  if (!svol) {
    p1 = p2 = p3 = 0;
    return false;
  }

  p1 = svol->order(0);
  p2 = svol->order(1);
  p3 = svol->order(2);
  return true;
}


bool ASMs3D::getSize (int& n1, int& n2, int& n3, int) const
{
  if (!svol) {
    n1 = n2 = n3 = 0;
    return false;
  }

  n1 = svol->numCoefs(0);
  n2 = svol->numCoefs(1);
  n3 = svol->numCoefs(2);
  return true;
}


size_t ASMs3D::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (!svol) return 0;

  if (ldim < 1 && lIndex > 0)
    return 1;
  else if (ldim < 2 && lIndex > 0 && lIndex <= 12)
    return svol->numCoefs((lIndex-1)/4) - svol->order((lIndex-1)/4) + 1;

  int n1 = svol->numCoefs(0) - svol->order(0) + 1;
  int n2 = svol->numCoefs(1) - svol->order(1) + 1;
  int n3 = svol->numCoefs(2) - svol->order(2) + 1;

  switch (lIndex)
    {
    case 1:
    case 2:
      return n2*n3;
    case 3:
    case 4:
      return n1*n3;
    case 5:
    case 6:
      return n1*n2;
    }

  return 0;
}


bool ASMs3D::getParameterDomain (Real2DMat& u, IntVec* corners) const
{
  u.resize(3,RealArray(2));
  for (int d = 0; d < 3; d++)
  {
    u[d].front() = svol->basis(d).startparam();
    u[d].back() = svol->basis(d).endparam();
  }

  if (corners)
  {
    corners->resize(8);
    corners->at(0) = 1;
    corners->at(1) = svol->numCoefs(0);
    corners->at(2) = svol->numCoefs(0)*(svol->numCoefs(1)-1) + 1;
    corners->at(3) = svol->numCoefs(0)*svol->numCoefs(1);
    int wnjump = corners->at(3)*(svol->numCoefs(2)-1);
    for (int i = 0; i < 4; i++)
      corners->at(4+i) = corners->at(i) + wnjump;
  }

  return true;
}


const RealArray& ASMs3D::getGaussPointParameters (Matrix& uGP, int dir,
                                                  int nGauss, const double* xi,
                                                  const Go::SplineVolume* spline) const
{
  if (!spline)
    spline = svol.get();

  int pm1 = spline->order(dir) - 1;
  RealArray::const_iterator uit = spline->basis(dir).begin() + pm1;

  int nCol = spline->numCoefs(dir) - pm1;
  uGP.resize(nGauss,nCol);

  double uprev = *(uit++);
  for (int j = 1; j <= nCol; ++uit, j++)
  {
    double ucurr = *uit;
    for (int i = 1; i <= nGauss; i++)
      uGP(i,j) = 0.5*((ucurr-uprev)*xi[i-1] + ucurr+uprev);
    uprev = ucurr;
  }

  return uGP;
}


void ASMs3D::getElementBorders (int iel, double* ub) const
{
  int p1 = svol->order(0);
  int p2 = svol->order(1);
  int p3 = svol->order(2);
  int n1 = svol->numCoefs(0) - p1 + 1;
  int n2 = svol->numCoefs(1) - p2 + 1;
  int i1 = p1 +  (iel-1) % n1;
  int i2 = p2 + ((iel-1) / n1) % n2;
  int i3 = p3 +  (iel-1) / (n1*n2);
  this->getElementBorders(i1-1,i2-1,i3-1,ub,ub+2,ub+4);
}


void ASMs3D::getElementBorders (int i1, int i2, int i3,
                                double* u, double* v, double* w) const
{
  RealArray::const_iterator uit = svol->basis(0).begin();
  RealArray::const_iterator vit = svol->basis(1).begin();
  RealArray::const_iterator wit = svol->basis(2).begin();

  for (int i = 0; i < 2; i++)
  {
    u[i] = uit[i1+i];
    v[i] = vit[i2+i];
    w[i] = wit[i3+i];
  }
}


double ASMs3D::getElementCorners (int i1, int i2, int i3, Vec3Vec& XC,
                                  RealArray* uC) const
{
  // Fetch parameter values of the element (knot-span) corners
  RealArray u(2), v(2), w(2);
  this->getElementBorders(i1,i2,i3,u.data(),v.data(),w.data());

  // Evaluate the spline volume at the corners to find physical coordinates
  int dim = svol->dimension();
  RealArray XYZ(dim*8);
#pragma omp critical
  svol->gridEvaluator(u,v,w,XYZ);

  XC.clear();
  XC.reserve(8);
  if (uC)
  {
    uC->clear();
    uC->reserve(24);
  }

  const double* pt = &XYZ.front();
  for (int k = 0; k < 2; k++)
    for (int j = 0; j < 2; j++)
      for (int i = 0; i < 2; i++, pt += dim)
      {
        XC.push_back(Vec3(pt,nsd));
        if (uC)
        {
          uC->push_back(u[i]);
          uC->push_back(v[j]);
          uC->push_back(w[k]);
        }
      }

  return getElementSize(XC);
}


void ASMs3D::getCornerPoints (int i1, int i2, int i3, PointVec& XC) const
{
  RealArray uC;
  Vec3Vec  XYZ;
  this->getElementCorners(i1,i2,i3,XYZ,&uC);

  XC.clear();
  XC.reserve(8);
  for (int i = 0; i < 4; i++)
    XC.push_back(utl::Point(XYZ[i], { uC[3*i], uC[3*i+1], uC[3*i+2] }));
}


bool ASMs3D::integrate (Integrand& integrand,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3D::integrate(I)");

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  if (myCache.empty())
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this));

  BasisFunctionCache& cache = *myCache.front();
  cache.setIntegrand(&integrand);
  if (!cache.init(use2ndDer ? 2 : 1))
    return false;

  // Get Gaussian quadrature points and weights
  const std::array<int,3>& ng = cache.nGauss();
  const std::array<const double*,3>& xg = cache.coord();
  const std::array<const double*,3>& wg = cache.weight();

  // Get the reduced integration quadrature points, if needed
  const double* xr = cache.coord(true)[0];
  const double* wr = cache.weight(true)[0];

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  ThreadGroups oneGroup;
  if (glInt.threadSafe()) oneGroup.oneStripe(nel);
  const ThreadGroups& groups = glInt.threadSafe() ? oneGroup : threadGroupsVol;


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < groups.size() && ok; g++)
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < groups[g].size(); t++)
    {
      FiniteElement fe(p1*p2*p3);
      Matrix   dNdu, Xnod, Jac;
      Matrix3D d2Ndu2, Hess;
      double   dXidu[3];
      double   param[3];
      Vec4     X(param,time.t);
      for (size_t l = 0; l < groups[g][t].size() && ok; l++)
      {
        int iel = groups[g][t][l];
        fe.iel = MLGE[iel];
        if (!this->isElementActive(fe.iel)) continue; // zero-volume element

#ifdef SP_DEBUG
        if (dbgElm < 0 && 1+iel != -dbgElm)
          continue; // Skipping all elements, except for -dbgElm
#endif

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel2;
        int i3 = p3 + iel / (nel1*nel2);

        // Get element volume in the parameter space
        double dV = 0.125*this->getParametricVolume(++iel);
        if (dV < 0.0)
        {
          ok = false; // topology error (probably logic error)
          break;
        }

        // Set up control point (nodal) coordinates for current element
        if (!this->getElementCoordinates(Xnod,iel))
        {
          ok = false;
          break;
        }

        if (useElmVtx)
          fe.h = this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        {
          // Element size in parametric space
          dXidu[0] = svol->knotSpan(0,i1-1);
          dXidu[1] = svol->knotSpan(1,i2-1);
          dXidu[2] = svol->knotSpan(2,i3-1);
        }

        if (integrand.getIntegrandType() & Integrand::AVERAGE)
        {
          // --- Compute average value of basis functions over the element -----

          fe.Navg.resize(p1*p2*p3,true);
          double vol = 0.0;
          size_t ip = 0;
          for (int k = 0; k < ng[2]; k++)
            for (int j = 0; j < ng[1]; j++)
              for (int i = 0; i < ng[0]; i++, ip++)
              {
                // Fetch basis function derivatives at current integration point
                const BasisFunctionVals& bfs = cache.getVals(iel-1,ip);

                // Compute Jacobian determinant of coordinate mapping
                // and multiply by weight of current integration point
                double detJac = utl::Jacobian(Jac,fe.dNdX,Xnod,bfs.dNdu,false);
                double weight = dV*wg[0][i]*wg[1][j]*wg[2][k];

                // Numerical quadrature
                fe.Navg.add(bfs.N,detJac*weight);
                vol += detJac*weight;
              }

          // Divide by element volume
          fe.Navg /= vol;
        }

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
        {
          // Compute the element center
          param[0] = 0.5*(cache.getParam(0,i1-p1,0) + cache.getParam(0,i1-p1,ng[0]-1));
          param[1] = 0.5*(cache.getParam(1,i2-p2,0) + cache.getParam(1,i2-p2,ng[1]-1));
          param[2] = 0.5*(cache.getParam(2,i3-p3,0) + cache.getParam(2,i3-p3,ng[2]-1));
          SplineUtils::point(X,param[0],param[1],param[2],svol.get());
          if (!useElmVtx)
          {
            // When element corner coordinates are not needed, store coordinates
            // and parameters of the element center in XC, for material usage
            fe.XC.resize(2);
            fe.XC.front() = X;
            fe.XC.back() = Vec3(param);
          }
        }

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
        int nRed = cache.nGauss(true)[0];
        if (!integrand.initElement(MNPC[iel-1],fe,X,nRed*nRed*nRed,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        if (xr)
        {
          // --- Selective reduced integration loop ----------------------------

          size_t ip = 0;
          for (int k = 0; k < nRed; k++)
            for (int j = 0; j < nRed; j++)
              for (int i = 0; i < nRed; i++, ip++)
              {
                // Local element coordinates of current integration point
                fe.xi   = xr[i];
                fe.eta  = xr[j];
                fe.zeta = xr[k];

                // Parameter values of current integration point
                fe.u = param[0] = cache.getParam(0,i1-p1,i,true);
                fe.v = param[1] = cache.getParam(1,i2-p2,j,true);
                fe.w = param[2] = cache.getParam(2,i3-p3,k,true);

                const BasisFunctionVals& bfs = cache.getVals(iel-1,ip,true);
                fe.N = bfs.N;

                // Compute Jacobian inverse and derivatives
                fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,bfs.dNdu);

                // Cartesian coordinates of current integration point
                X.assign(Xnod * fe.N);

                // Compute the reduced integration terms of the integrand
                fe.detJxW *= dV*wr[i]*wr[j]*wr[k];
                if (!integrand.reducedInt(*A,fe,X))
                  ok = false;
              }
        }


        // --- Integration loop over all Gauss points in each direction --------

        size_t ip = 0;
        int jp = (((i3-p3)*nel2 + i2-p2)*nel1 + i1-p1)*ng[0]*ng[1]*ng[2];
        fe.iGP = firstIp + jp; // Global integration point counter
        for (int k = 0; k < ng[2]; k++)
          for (int j = 0; j < ng[1]; j++)
            for (int i = 0; i < ng[0]; i++, ++ip, ++fe.iGP)
            {
              // Local element coordinates of current integration point
              fe.xi   = xg[0][i];
              fe.eta  = xg[1][j];
              fe.zeta = xg[2][k];

              // Parameter values of current integration point
              fe.u = param[0] = cache.getParam(0,i1-p1,i);
              fe.v = param[1] = cache.getParam(1,i2-p2,j);
              fe.w = param[2] = cache.getParam(2,i3-p3,k);

              // Fetch basis function derivatives at current integration point
              const BasisFunctionVals& bfs = cache.getVals(iel-1,ip);
              fe.N = bfs.N;

              // Compute Jacobian inverse of coordinate mapping and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,bfs.dNdu);
              if (fe.detJxW == 0.0) continue; // skip singular points

              // Compute Hessian of coordinate mapping and 2nd order derivatives
              if (use2ndDer)
                if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,bfs.d2Ndu2,fe.dNdX))
                  ok = false;

              // Compute G-matrix
              if (integrand.getIntegrandType() & Integrand::G_MATRIX)
                utl::getGmat(Jac,dXidu,fe.G);

#if SP_DEBUG > 4
              if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
                std::cout <<"\n"<< fe;
#endif

              // Cartesian coordinates of current integration point
              X.assign(Xnod * fe.N);

              // Evaluate the integrand and accumulate element contributions
              fe.detJxW *= dV*wg[0][i]*wg[1][j]*wg[2][k];
#ifndef USE_OPENMP
              PROFILE3("Integrand::evalInt");
#endif
              if (!integrand.evalInt(*A,fe,time,X))
                ok = false;
            }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,fe,time,firstIp+jp))
          ok = false;

        // Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();

#ifdef SP_DEBUG
        if (iel == -dbgElm)
          break; // Skipping all elements, except for -dbgElm
#endif
      }
    }

  cache.finalizeAssembly();
  return ok;
}


bool ASMs3D::integrate (Integrand& integrand,
			GlobalIntegral& glInt,
			const TimeDomain& time,
                        const Real3DMat& itgPts)
{
  if (!svol) return true; // silently ignore empty patches

  if (integrand.getReducedIntegration(2) != 0)
  {
    std::cerr <<" *** ASMs3D::integrate(Integrand&,GlobalIntegral&,"
              <<"const TimeDomain&,const Real3DMat&): Available for standard"
              <<" integrands only."<< std::endl;
    return false;
  }

  PROFILE2("ASMs3D::integrate(I)");

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  // Evaluate basis function derivatives at all integration points
  size_t i, j, k;
  std::vector<size_t> MPitg(itgPts.size()+1,0);
  for (i = MPitg.front() = 0; i < itgPts.size(); i++)
    MPitg[i+1] = MPitg[i] + itgPts[i].size();
  size_t nPoints = MPitg.back();
  std::vector<Go::BasisDerivs>  spline(use2ndDer ? 0 : nPoints);
  std::vector<Go::BasisDerivs2> spline2(!use2ndDer ? 0 : nPoints);
  for (i = k = 0; i < itgPts.size(); i++)
    for (j = 0; j < itgPts[i].size(); j++, k++)
    {
      const RealArray& itgPt = itgPts[i][j];
      if (use2ndDer)
        svol->computeBasis(itgPt[0],itgPt[1],itgPt[2],spline2[k]);
      else
        svol->computeBasis(itgPt[0],itgPt[1],itgPt[2],spline[k]);
    }

  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);

  const int nel1 = n1 - p1 + 1;
  const int nel2 = n2 - p2 + 1;

  ThreadGroups oneGroup;
  if (glInt.threadSafe()) oneGroup.oneStripe(nel);
  const ThreadGroups& groups = glInt.threadSafe() ? oneGroup : threadGroupsVol;


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < groups.size() && ok; g++)
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < groups[g].size(); t++)
    {
      FiniteElement fe(p1*p2*p3);
      Matrix   dNdu, Xnod, Jac;
      Matrix3D d2Ndu2, Hess;
      double   dXidu[3];
      Vec4     X(nullptr,time.t);
      for (size_t e = 0; e < groups[g][t].size() && ok; e++)
      {
        int iel = groups[g][t][e];
        if (itgPts[iel].empty()) continue; // no points in this element

        fe.iel = MLGE[iel];
        if (!this->isElementActive(fe.iel)) continue; // zero-volume element

#ifdef SP_DEBUG
        if (dbgElm < 0 && 1+iel != -dbgElm)
          continue; // Skipping all elements, except for -dbgElm
#endif

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel2;
        int i3 = p3 + iel / (nel1*nel2);

        // Get element volume in the parameter space
        double dV = 0.125*this->getParametricVolume(++iel);
        if (dV < 0.0)
        {
          ok = false; // topology error (probably logic error)
          break;
        }

        // Set up control point (nodal) coordinates for current element
        if (!this->getElementCoordinates(Xnod,iel))
        {
          ok = false;
          break;
        }

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
        {
          // Compute the element center
          fe.h = this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);
	  X = 0.125*(fe.XC[0]+fe.XC[1]+fe.XC[2]+fe.XC[3]+
		     fe.XC[4]+fe.XC[5]+fe.XC[6]+fe.XC[7]);
        }
        else if (useElmVtx)
          fe.h = this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        {
          // Element size in parametric space
          dXidu[0] = svol->knotSpan(0,i1-1);
          dXidu[1] = svol->knotSpan(1,i2-1);
          dXidu[2] = svol->knotSpan(2,i3-1);
        }

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
        if (!integrand.initElement(MNPC[iel-1],fe,X,0,*A))
        {
          A->destruct();
          ok = false;
          break;
        }


        // --- Integration loop over all quadrature points in this element -----

        size_t jp = MPitg[iel]; // Patch-wise integration point counter
        fe.iGP = firstIp + jp;  // Global integration point counter

        for (size_t ip = 0; ip < itgPts[iel].size(); ip++, jp++, fe.iGP++)
        {
          // Parameter values of current integration point
          fe.u = itgPts[iel][ip][0];
          fe.v = itgPts[iel][ip][1];
          fe.w = itgPts[iel][ip][2];

          // Fetch basis function derivatives at current integration point
          if (use2ndDer)
            SplineUtils::extractBasis(spline2[jp],fe.N,dNdu,d2Ndu2);
          else
            SplineUtils::extractBasis(spline[jp],fe.N,dNdu);

          // Compute Jacobian inverse of coordinate mapping and derivatives
          fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
          if (fe.detJxW == 0.0) continue; // skip singular points

          // Compute Hessian of coordinate mapping and 2nd order derivatives
          if (use2ndDer)
            if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,fe.dNdX))
              ok = false;

          // Compute G-matrix
          if (integrand.getIntegrandType() & Integrand::G_MATRIX)
            utl::getGmat(Jac,dXidu,fe.G);

#if SP_DEBUG > 4
          if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
            std::cout <<"\n"<< fe;
#endif

          // Cartesian coordinates of current integration point
          X.assign(Xnod * fe.N);
          X.u = itgPts[iel][ip].data();

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= dV*itgPts[iel][ip][3];
#ifndef USE_OPENMP
          PROFILE3("Integrand::evalInt");
#endif
          if (!integrand.evalInt(*A,fe,time,X))
            ok = false;
        }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,fe,time,firstIp+MPitg[iel]))
          ok = false;

        // Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();

#ifdef SP_DEBUG
        if (iel == -dbgElm)
          break; // Skipping all elements, except for -dbgElm
#endif
      }
    }

  return ok;
}


bool ASMs3D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time,
                        const ASM::InterfaceChecker& iChk)
{
  if (!svol) return true; // silently ignore empty patches
  if (!(integrand.getIntegrandType() & Integrand::INTERFACE_TERMS)) return true;

  std::cerr << __PRETTY_FUNCTION__ <<": Not implemented yet."<< std::endl;
  return false;
}


bool ASMs3D::integrate (Integrand& integrand, int lIndex,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3D::integrate(B)");

  std::map<char,ThreadGroups>::const_iterator tit;
  if ((tit = threadGroupsFace.find(lIndex%10)) == threadGroupsFace.end())
  {
    std::cerr <<" *** ASMs3D::integrate: No thread groups for face "
              << lIndex%10 << std::endl;
    return false;
  }
  const ThreadGroups& threadGrp = tit->second;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  // Get Gaussian quadrature points and weights
  // and compute parameter values of the Gauss points over the whole patch face
  std::array<int,3> ng;
  std::array<const double*,3> xg, wg;
  std::array<Matrix,3> gpar;
  for (int d = 0; d < 3; d++)
    if (-1-d == faceDir)
    {
      ng[d] = 1;
      xg[d] = nullptr;
      wg[d] = nullptr;
      gpar[d].resize(1,1);
      gpar[d].fill(svol->startparam(d));
    }
    else if (1+d == faceDir)
    {
      ng[d] = 1;
      xg[d] = nullptr;
      wg[d] = nullptr;
      gpar[d].resize(1,1);
      gpar[d].fill(svol->endparam(d));
    }
    else
    {
      int n = this->getNoGaussPt(svol->order(d),true);
      ng[d] = integrand.getBouIntegrationPoints(n);
      xg[d] = GaussQuadrature::getCoord(ng[d]);
      wg[d] = GaussQuadrature::getWeight(ng[d]);
      if (xg[d] && wg[d])
        this->getGaussPointParameters(gpar[d],d,ng[d],xg[d]);
      else
        return false;
    }

  const int tt1 = t1 > t2 ? t2-1 : t1-1;
  const int tt2 = t1 > t2 ? t1-1 : t2-1;
  const int nG1 = ng[tt1];
  const int nG2 = ng[tt2];

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

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

  // Integrate the extraordinary elements?
  size_t doXelms = 0;
  if (integrand.getIntegrandType() & Integrand::XO_ELEMENTS)
    if ((doXelms = (n1-p1+1)*(n2-p2+1)*(n3-p3+1))*2 > MNPC.size())
    {
      std::cerr <<" *** ASMs3D::integrate: Too few XO-elements "
                << MNPC.size() - doXelms << std::endl;
      return false;
    }

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (size_t g = 0; g < threadGrp.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGrp[g].size(); t++)
    {
      FiniteElement fe(p1*p2*p3);
      fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
      fe.u = gpar[0](1,1);
      fe.v = gpar[1](1,1);
      fe.w = gpar[2](1,1);

      Matrix dNdu, Xnod, Jac;
      double param[3] = { fe.u, fe.v, fe.w };
      Vec4   X(param,time.t);
      Vec3   normal;
      double dXidu[3];

      for (size_t l = 0; l < threadGrp[g][t].size() && ok; l++)
      {
        int iel = threadGrp[g][t][l];
        fe.iel = abs(MLGE[doXelms+iel]);
        if (!this->isElementActive(fe.iel)) continue; // zero-volume element

#ifdef SP_DEBUG
        if (dbgElm < 0 && 1+iel != -dbgElm)
          continue; // Skipping all elements, except for -dbgElm
#endif

        int i1 = p1 + iel % nel1;
        int i2 = p2 + (iel / nel1) % nel2;
        int i3 = p3 + iel / (nel1*nel2);

        // Get element face area in the parameter space
        double dA = 0.25*this->getParametricArea(++iel,abs(faceDir));
        if (dA < 0.0) // topology error (probably logic error)
        {
          ok = false;
          break;
        }

        // Set up control point coordinates for current element
        if (!this->getElementCoordinates(Xnod,iel))
        {
          ok = false;
          break;
        }

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
          fe.h = this->getElementCorners(i1-1,i2-1,i3-1,fe.XC);

        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        {
          // Element size in parametric space
          dXidu[0] = svol->knotSpan(0,i1-1);
          dXidu[1] = svol->knotSpan(1,i2-1);
          dXidu[2] = svol->knotSpan(2,i3-1);
        }

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
        if (!integrand.initElementBou(MNPC[doXelms+iel-1],*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        // Define some loop control variables depending on which face we are on
        int nf1 = 0, j1 = 0, j2 = 0;
        switch (abs(faceDir))
        {
          case 1: nf1 = nel2; j2 = i3-p3; j1 = i2-p2; break;
          case 2: nf1 = nel1; j2 = i3-p3; j1 = i1-p1; break;
          case 3: nf1 = nel1; j2 = i2-p2; j1 = i1-p1; break;
        }


        // --- Integration loop over all Gauss points in each direction --------

        int ip = (j2*nf1*nG2 + j1)*nG1;
        int jp = (j2*nf1 + j1)*nG1*nG2;
        fe.iGP = firstp + jp; // Global integration point counter

        for (int j = 0; j < nG2; j++, ip += nG1*(nf1-1))
          for (int i = 0; i < nG1; i++, ip++, fe.iGP++)
          {
#if SP_DEBUG > 4
            std::cout <<"Elem "<< iel <<": "<< i1-p1 <<" "<< i2-p2 <<" "<< i3-p3
                      <<", Face "<< faceDir <<": "<< j1 <<" "<< j2
                      <<", Point "<< ip <<": "<< i <<" "<< j <<" "
                      << spline.size() << std::endl;
#endif

            // Local element coordinates and parameter values
            // of current integration point
            int k1 = 0, k2 = 0, k3 = 0;
            switch (abs(faceDir))
            {
              case 1: k2 = i; k3 = j; break;
              case 2: k1 = i; k3 = j; break;
              case 3: k1 = i; k2 = j; break;
            }

            if (xg[0])
            {
              fe.xi = xg[0][k1];
              fe.u = param[0] = gpar[0](k1+1,i1-p1+1);
            }
            if (xg[1])
            {
              fe.eta = xg[1][k2];
              fe.v = param[1] = gpar[1](k2+1,i2-p2+1);
            }
            if (xg[2])
            {
              fe.zeta = xg[2][k3];
              fe.w = param[2] = gpar[2](k3+1,i3-p3+1);
            }

            // Fetch basis function derivatives at current integration point
            SplineUtils::extractBasis(spline[ip],fe.N,dNdu);

            // Compute basis function derivatives and the face normal
            fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
            if (fe.detJxW == 0.0) continue; // skip singular points

            if (faceDir < 0) normal *= -1.0;

            // Compute G-matrix
            if (integrand.getIntegrandType() & Integrand::G_MATRIX)
              utl::getGmat(Jac,dXidu,fe.G);

#if SP_DEBUG > 4
            if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
              std::cout <<"\n"<< fe;
#endif

            // Cartesian coordinates of current integration point
            X.assign(Xnod * fe.N);

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= dA*wg[tt1][i]*wg[tt2][j];
            if (!integrand.evalBou(*A,fe,time,X,normal))
              ok = false;
          }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElementBou(*A,fe,time))
          ok = false;

        // Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();

#ifdef SP_DEBUG
        if (dbgElm < 0 && iel == -dbgElm)
          break; // Skipping all elements, except for -dbgElm
#endif
      }
    }
  }

  return ok;
}


bool ASMs3D::integrateEdge (Integrand& integrand, int lEdge,
			    GlobalIntegral& glInt,
			    const TimeDomain& time)
{
  if (!svol) return true; // silently ignore empty patches

  PROFILE2("ASMs3D::integrate(E)");

  // Get Gaussian quadrature points and weights
  int ng = this->getNoGaussPt(svol->order((lEdge-1)/4),true);
  const double* xg = GaussQuadrature::getCoord(ng);
  const double* wg = GaussQuadrature::getWeight(ng);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points along the whole patch edge
  std::array<Matrix,3> gpar;
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
      double uprev = *(uit++);
      int nCol = svol->numCoefs(d) - pm1;
      gpar[d].resize(ng,nCol);
      for (int j = 1; j <= nCol; ++uit, j++)
      {
	double ucurr = *uit;
	for (int i = 1; i <= ng; i++)
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

  std::map<char,size_t>::const_iterator iit = firstBp.find(lEdge);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  FiniteElement fe(p1*p2*p3);
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);
  fe.w = gpar[2](1,1);
  if (gpar[0].size() == 1) fe.xi = fe.u == svol->startparam(0) ? -1.0 : 1.0;
  if (gpar[1].size() == 1) fe.eta = fe.v == svol->startparam(1) ? -1.0 : 1.0;
  if (gpar[2].size() == 1) fe.zeta = fe.w == svol->startparam(2) ? -1.0 : 1.0;

  Matrix dNdu, Xnod, Jac;
  double param[3] = { 0.0, 0.0, 0.0 };
  Vec4   X(param,time.t);
  Vec3   tang;


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
        fe.iel = MLGE[iel-1];
        if (!this->isElementActive(fe.iel)) continue; // zero-volume element

        if (!myElms.empty() &&
            std::find(myElms.begin(), myElms.end(), fe.iel-1) == myElms.end())
          continue;

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
	int ip = MNPC[iel-1][svol->order(0)*svol->order(1)*svol->order(2)-1];
#ifdef INDEX_CHECK
	if (ip < 0 || (size_t)ip >= nnod) return false;
#endif
	if (lEdge < 5)
	{
	  dS = svol->knotSpan(0,nodeInd[ip].I);
	  ip = (i1-p1)*ng;
	}
	else if (lEdge < 9)
	{
	  dS = svol->knotSpan(1,nodeInd[ip].J);
	  ip = (i2-p2)*ng;
	}
	else if (lEdge < 13)
	{
	  dS = svol->knotSpan(2,nodeInd[ip].K);
	  ip = (i3-p3)*ng;
	}

	// Set up control point coordinates for current element
	if (!this->getElementCoordinates(Xnod,iel)) return false;

	// Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
        bool ok = integrand.initElementBou(MNPC[iel-1],*A);


	// --- Integration loop over all Gauss points along the edge -----------

	fe.iGP = firstp + ip; // Global integration point counter

	for (int i = 0; i < ng && ok; i++, ip++, fe.iGP++)
	{
	  // Parameter values of current integration point
	  if (gpar[0].size() > 1) fe.u = param[0] = gpar[0](i+1,i1-p1+1);
	  if (gpar[1].size() > 1) fe.v = param[1] = gpar[1](i+1,i2-p2+1);
	  if (gpar[2].size() > 1) fe.w = param[2] = gpar[2](i+1,i3-p3+1);

	  // Fetch basis function derivatives at current integration point
	  SplineUtils::extractBasis(spline[ip],fe.N,dNdu);

	  // Compute basis function derivatives and the edge tangent
	  fe.detJxW = utl::Jacobian(Jac,tang,fe.dNdX,Xnod,dNdu,1+(lEdge-1)/4);
	  if (fe.detJxW == 0.0) continue; // skip singular points

	  // Cartesian coordinates of current integration point
	  X.assign(Xnod * fe.N);

	  // Evaluate the integrand and accumulate element contributions
	  fe.detJxW *= 0.5*dS*wg[i];
          ok = integrand.evalBou(*A,fe,time,X,tang);
	}

        // Finalize the element quantities
        if (ok && !integrand.finalizeElementBou(*A,fe,time))
          ok = false;

	// Assembly of global system integral
	if (ok && !glInt.assemble(A->ref(),fe.iel))
	  ok = false;

	A->destruct();

	if (!ok) return false;
      }

  return true;
}


int ASMs3D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  const Go::SplineVolume* geo = this->getBasis(ASM::GEOMETRY_BASIS);
  if (!geo) return -3;

  for (int i = 0; i < 3; i++)
    param[i] = (1.0-xi[i])*geo->startparam(i) + xi[i]*geo->endparam(i);

  SplineUtils::point(X,param[0],param[1],param[2],geo);

  // Check if this point matches any of the control points (nodes)
  return this->searchCtrlPt(svol->coefs_begin(),svol->coefs_end(),
                            X,svol->dimension());
}


int ASMs3D::findElementContaining (const double* param) const
{
  if (!svol) return -2;

  int p1   = svol->order(0) - 1;
  int p2   = svol->order(1) - 1;
  int p3   = svol->order(2) - 1;
  int nel1 = svol->numCoefs(0) - p1;
  int nel2 = svol->numCoefs(1) - p2;
  int uEl  = svol->basis(0).knotInterval(param[0]) - p1;
  int vEl  = svol->basis(1).knotInterval(param[1]) - p2;
  int wEl  = svol->basis(2).knotInterval(param[2]) - p3;
  return 1 + uEl + (vEl + wEl*nel2)*nel1;
}


double ASMs3D::findPoint (Vec3& X, double* param) const
{
  if (!svol) return -1.0;

  // Use with caution, very slow!
  double dist;
  Go::Point Xpt(X.x,X.y,X.z), Xfound(3);
#pragma omp critical
  svol->closestPoint(Xpt, param[0], param[1], param[2], Xfound, dist, 1.0e-5);
  for (int i = 0; i < svol->dimension(); i++) X[i] = Xfound[i];
  return dist;
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

  RealArray::const_iterator uit = svol->basis(dir).begin() + svol->basis(dir).order()-1;
  RealArray::const_iterator uend = svol->basis(dir).begin() + svol->basis(dir).numCoefs()+1;

  double ucurr = 0.0, uprev = *(uit++);
  while (uit != uend)
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


bool ASMs3D::tesselate (ElementBlock& grid, const int* npe) const
{
  // Compute parameter values of the nodal points
  std::array<RealArray,3> gpar;
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
  for (i = j = 0; i < grid.getNoNodes(); i++, j += svol->dimension()) {
    grid.setParams(i,gpar[0][i % nx], gpar[1][i/nx % ny], gpar[2][i / (nx*ny)]);
    grid.setCoor(i,XYZ[j],XYZ[j+1],XYZ[j+2]);
  }

  // Establish the block grid topology
  int ie, nse1 = npe[0] - 1;
  int je, nse2 = npe[1] - 1;
  int ke, nse3 = npe[2] - 1;
  int nel1 = (nx-1)/nse1;
  int nel2 = (ny-1)/nse2;
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
                           const int* npe, int n_f, bool piola) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  if (piola)
    return this->evalSolutionPiola(sField,locSol,gpar.data(),true);
  else
    return this->evalSolution(sField,locSol,gpar.data(),true,0,n_f);
}


bool ASMs3D::evalSolution (Matrix& sField, const Vector& locSol,
                           const RealArray* gpar,
                           bool regular, int deriv, int) const
{
  PROFILE2("ASMs3D::evalSol(P)");

  // Evaluate the basis functions at all points
  size_t nPoints = gpar[0].size();
  std::vector<Go::BasisPts>     spline0(regular || deriv != 0 ? 0 : nPoints);
  std::vector<Go::BasisDerivs>  spline1(regular || deriv != 1 ? 0 : nPoints);
  std::vector<Go::BasisDerivs2> spline2(regular || deriv != 2 ? 0 : nPoints);
  if (regular)
  {
    PROFILE2("Spline evaluation");
    nPoints *= gpar[1].size()*gpar[2].size();
    switch (deriv) {
    case 0:
      svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline0);
      break;
    case 1:
      svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
      break;
    case 2:
      svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
      break;
    default:
      return false;
    }
  }
  else if (nPoints == gpar[1].size() && nPoints == gpar[2].size())
  {
    PROFILE2("Spline evaluation");
    for (size_t i = 0; i < nPoints; i++)
      switch (deriv) {
      case 0:
        svol->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline0[i]);
        break;
      case 1:
        svol->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline1[i]);
        break;
      case 2:
        svol->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline2[i]);
        break;
      default:
        return false;
      }
  }
  else
    return false;

  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  size_t nComp = locSol.size() / (n1*n2*n3);

  Vector   ptSol;
  Matrix   dNdu, dNdX, Xnod, Xtmp, Jac, eSol, ptDer;
  Matrix3D d2Ndu2, d2NdX2, Hess, ptDer2;

  // Fetch nodal (control point) coordinates
  this->getNodalCoordinates(Xnod);

  // Evaluate the primary solution field at each point
  sField.resize(nComp*int(pow(3.0,deriv)),nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    IntVec ip;
    switch (deriv) {

    case 0: // Evaluate the solution
      scatterInd(n1,n2,n3,p1,p2,p3,spline0[i].left_idx,ip);
      utl::gather(ip,nComp,locSol,Xtmp);
      Xtmp.multiply(spline0[i].basisValues,ptSol);
      sField.fillColumn(1+i,ptSol);
      break;

    case 1: // Evaluate first derivatives of the solution
      scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip);
      SplineUtils::extractBasis(spline1[i],ptSol,dNdu);
      utl::gather(ip,3,Xnod,Xtmp);
      utl::Jacobian(Jac,dNdX,Xtmp,dNdu);
      utl::gather(ip,nComp,locSol,Xtmp);
      ptDer.multiply(Xtmp,dNdX);
      sField.fillColumn(1+i,ptDer);
      break;

    case 2: // Evaluate second derivatives of the solution
      scatterInd(n1,n2,n3,p1,p2,p3,spline2[i].left_idx,ip);
      SplineUtils::extractBasis(spline2[i],ptSol,dNdu,d2Ndu2);
      utl::gather(ip,3,Xnod,Xtmp);
      utl::Jacobian(Jac,dNdX,Xtmp,dNdu);
      utl::Hessian(Hess,d2NdX2,Jac,Xtmp,d2Ndu2,dNdX);
      utl::gather(ip,nComp,locSol,Xtmp);
      ptDer2.multiply(Xtmp,d2NdX2);
      sField.fillColumn(1+i,ptDer2);
      break;
    }
  }

  return true;
}


bool ASMs3D::evalProjSolution (Matrix& sField, const Vector& locSol,
                               const int* npe, int n_f) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the projected solution at all sampling points
  if (!this->separateProjectionBasis())
    return this->evalSolution(sField,locSol,gpar.data(),true,0,n_f);

  // The projection uses a separate basis, need to interpolate
  Fields* f = this->getProjectedFields(locSol);
  if (!f) return false;

  // Evaluate the projected solution field at each point
  Vector vals;
  sField.resize(f->getNoFields(),gpar[0].size()*gpar[1].size()*gpar[2].size());

  size_t ipt = 0;
  for (double w : gpar[1])
    for (double v : gpar[1])
      for (double u : gpar[0])
      {
        f->valueFE(ItgPoint(u,v,w),vals);
        sField.fillColumn(++ipt,vals);
      }

  delete f;
  return true;
}


bool ASMs3D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			   const int* npe, char project) const
{
  // Project the secondary solution onto the spline basis
  Go::SplineVolume* v = nullptr;
  if (project == 'A')
    v = this->projectSolutionLocalApprox(integrand);
  else if (project == 'L')
    v = this->projectSolutionLocal(integrand);
  else if (project == 'W')
    v = this->projectSolutionLeastSquare(integrand);
  else if (project == 'D' || !npe)
    v = this->projectSolution(integrand);

  if (npe)
  {
    // Compute parameter values of the result sampling points
    std::array<RealArray,3> gpar;
    if (this->getGridParameters(gpar[0],0,npe[0]-1) &&
	this->getGridParameters(gpar[1],1,npe[1]-1) &&
	this->getGridParameters(gpar[2],2,npe[2]-1))
    {
      if (!project)
        // Evaluate the secondary solution directly at all sampling points
        return this->evalSolution(sField,integrand,gpar.data());
      else if (v)
      {
        // Evaluate the projected field at the result sampling points
        sField.resize(v->dimension(),
                      gpar[0].size()*gpar[1].size()*gpar[2].size());
        v->gridEvaluator(gpar[0],gpar[1],gpar[2],sField);
        delete v;
        return true;
      }
    }
    else if (v)
      delete v;
  }
  else if (v)
  {
    // Extract control point values from the spline object
    sField.resize(v->dimension(),
                  v->numCoefs(0)*v->numCoefs(1)*v->numCoefs(2));
    sField.fill(&(*v->coefs_begin()));
    delete v;
    return true;
  }

  std::cerr <<" *** ASMs3D::evalSolution: Failure!";
  if (project) std::cerr <<" project="<< project;
  std::cerr << std::endl;
  return false;
}


bool ASMs3D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			   const RealArray* gpar, bool regular) const
{
  PROFILE2("ASMs3D::evalSol(S)");

  sField.resize(0,0);

  // Evaluate the basis functions and their derivatives at all points
  size_t nPoints = gpar[0].size();
  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  std::vector<Go::BasisDerivs>  spline1(regular ||  use2ndDer ? 0 : nPoints);
  std::vector<Go::BasisDerivs2> spline2(regular || !use2ndDer ? 0 : nPoints);
  if (regular)
  {
    PROFILE2("Spline evaluation");
    nPoints *= gpar[1].size()*gpar[2].size();
    if (use2ndDer)
      svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline2);
    else
      svol->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline1);
  }
  else if (nPoints == gpar[1].size() && nPoints == gpar[2].size())
  {
    PROFILE2("Spline evaluation");
    for (size_t i = 0; i < nPoints; i++)
      if (use2ndDer)
        svol->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline2[i]);
      else
        svol->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline1[i]);
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

  FiniteElement fe(p1*p2*p3,firstIp);
  Vector        solPt;
  Matrix        dNdu, Jac;
  Matrix3D      d2Ndu2, Hess;

  // Evaluate the secondary solution field at each point
  for (size_t i = 0; i < nPoints; i++, fe.iGP++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntVec ip;
    if (use2ndDer)
    {
      scatterInd(n1,n2,n3,p1,p2,p3,spline2[i].left_idx,ip);
      fe.u = spline2[i].param[0];
      fe.v = spline2[i].param[1];
      fe.w = spline2[i].param[2];
    }
    else
    {
      scatterInd(n1,n2,n3,p1,p2,p3,spline1[i].left_idx,ip);
      fe.u = spline1[i].param[0];
      fe.v = spline1[i].param[1];
      fe.w = spline1[i].param[2];
    }

    // Fetch associated control point coordinates
    utl::gather(ip,3,Xnod,Xtmp);

    // Fetch basis function derivatives at current integration point
    if (use2ndDer)
      SplineUtils::extractBasis(spline2[i],fe.N,dNdu,d2Ndu2);
    else
      SplineUtils::extractBasis(spline1[i],fe.N,dNdu);

    // Compute the Jacobian inverse and derivatives
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xtmp,dNdu);
    if (fe.detJxW == 0.0) continue; // skip singular points

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (use2ndDer)
      if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xtmp,d2Ndu2,fe.dNdX))
        continue;

#if SP_DEBUG > 4
    std::cout <<"\n"<< fe;
#endif

    // Now evaluate the solution field
    utl::Point X4(Xtmp*fe.N,{fe.u,fe.v,fe.w});
    if (!integrand.evalSol(solPt,fe,X4,ip))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMs3D::generateThreadGroups (const Integrand& integrand, bool silence,
                                   bool ignoreGlobalLM)
{
  if (threadGroupsVol.stripDir == ThreadGroups::NONE)
    threadGroupsVol.oneGroup(nel);
  else
    this->generateThreadGroups(svol->order(0)-1, svol->order(1)-1,
                               svol->order(2)-1, silence, ignoreGlobalLM);
}


void ASMs3D::generateThreadGroups (size_t strip1, size_t strip2, size_t strip3,
                                   bool silence, bool ignoreGlobalLM)
{
  const int p1 = svol->order(0) - 1;
  const int p2 = svol->order(1) - 1;
  const int p3 = svol->order(2) - 1;
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  std::vector<bool> el1, el2, el3;
  el1.reserve(n1 - p1);
  el2.reserve(n2 - p2);
  el3.reserve(n3 - p3);

  int ii;
  for (ii = p1; ii < n1; ii++)
    el1.push_back(svol->knotSpan(0,ii) > 0.0);
  for (ii = p2; ii < n2; ii++)
    el2.push_back(svol->knotSpan(1,ii) > 0.0);
  for (ii = p3; ii < n3; ii++)
    el3.push_back(svol->knotSpan(2,ii) > 0.0);

  threadGroupsVol.calcGroups(el1,el2,el3,strip1,strip2,strip3);
  if (silence || threadGroupsVol.size() < 2) return;

  IFEM::cout <<"\nMultiple threads are utilized during element assembly.";
  for (size_t i = 0; i < threadGroupsVol.size(); i++)
  {
    std::vector< std::set<int> > nodes(threadGroupsVol[i].size());

    IFEM::cout <<"\n Thread group "<< i+1;
    for (size_t j = 0; j < threadGroupsVol[i].size(); j++)
    {
      IFEM::cout <<"\n\tthread "<< j+1
                 << ": "<< threadGroupsVol[i][j].size() <<" elements";
      size_t k, l, nzerovol = 0;
      for (k = 0; k < threadGroupsVol[i][j].size(); k++)
      {
        int iel = threadGroupsVol[i][j][k];
        if (MLGE[iel] > 0)
          for (l = 0; l < MNPC[iel].size(); l++)
            nodes[j].insert(MNPC[iel][l]);
        else
          nzerovol++;
      }
      if (nzerovol)
        IFEM::cout <<" ("<< threadGroupsVol[i][j].size() - nzerovol <<" real)";

      // Verify that the nodes on this thread are not present on the others
      this->checkThreadGroups(nodes, j, ignoreGlobalLM);
    }
  }
  IFEM::cout << std::endl;
}


void ASMs3D::generateThreadGroups (char lIndex, bool silence, bool)
{
  std::map<char,ThreadGroups>::iterator tit = threadGroupsFace.find(lIndex);
  if (tit != threadGroupsFace.end())
  {
    if (tit->second.stripDir == ThreadGroups::NONE)
      tit->second.oneGroup(nel);
    return;
  }

  const int p1 = svol->order(0) - 1;
  const int p2 = svol->order(1) - 1;
  const int p3 = svol->order(2) - 1;
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);

  // Flag the non-zero knot-spans
  std::vector<bool> el1, el2, el3;
  if (lIndex > 2) {
    el1.reserve(n1-p1);
    for (int i = p1; i < n1; i++)
      el1.push_back(svol->knotSpan(0,i) > 0.0);
  }
  if (lIndex < 3 || lIndex > 4) {
    el2.reserve(n2-p2);
    for (int i = p2; i < n2; i++)
      el2.push_back(svol->knotSpan(1,i) > 0.0);
  }
  if (lIndex < 6) {
    el3.reserve(n3-p3);
    for (int i = p3; i < n3; i++)
      el3.push_back(svol->knotSpan(2,i) > 0.0);
  }

  ThreadGroups& fGrp = threadGroupsFace[lIndex];
  switch (lIndex)
    {
    case 1:
    case 2:
      fGrp.calcGroups(el2,el3,p2,p3);
      break;
    case 3:
    case 4:
      fGrp.calcGroups(el1,el3,p1,p3);
      break;
    default:
      fGrp.calcGroups(el1,el2,p1,p2);
    }

  // Find elements that are on the boundary face 'lIndex'
  IntVec map;
  this->findBoundaryElms(map,lIndex);

  fGrp.applyMap(map);

  if (silence || fGrp.size() < 2) return;

  for (size_t i = 0; i < fGrp.size(); i++)
  {
    IFEM::cout <<"\n Thread group "<< i+1 <<" for boundary face "<< (int)lIndex;
    for (size_t j = 0; j < fGrp[i].size(); j++)
      IFEM::cout <<"\n\tthread "<< j+1
                 << ": "<< fGrp[i][j].size() <<" elements";
  }
  IFEM::cout << std::endl;
}


int ASMs3D::getLastItgElmNode () const
{
  return svol->order(0)*svol->order(1)*svol->order(2)-1;
}


bool ASMs3D::getNoStructElms (int& n1, int& n2, int& n3) const
{
  n1 = svol->numCoefs(0) - svol->order(0) + 1;
  n2 = svol->numCoefs(1) - svol->order(1) + 1;
  n3 = svol->numCoefs(2) - svol->order(2) + 1;

  return true;
}


void ASMs3D::evaluateBasis (double u, double v, double w, Vector& N) const
{
  Go::BasisPts spline;
  svol->computeBasis(u,v,w,spline);
  N = spline.basisValues;
}


void ASMs3D::extractBasis (double u, double v, double w,
                           Vector& N, Matrix& dNdu, bool fromRight) const
{
  Go::BasisDerivs spline;
  svol->computeBasis(u,v,w,spline,fromRight);
  SplineUtils::extractBasis(spline,N,dNdu);
}


void ASMs3D::extractBasis (double u, double v, double w, Vector& N,
                           Matrix& dNdu, Matrix3D& d2Ndu2, bool fromRight) const
{
  Go::BasisDerivs2 spline;
  svol->computeBasis(u,v,w,spline,fromRight);
  SplineUtils::extractBasis(spline,N,dNdu,d2Ndu2);
}


void ASMs3D::extractBasis (double u, double v, double w, int dir, int p,
                           Vector& dN, bool fromRight) const
{
  /* TODO: Missing GoTools function
  Go::BasisDerivsU spline;
  svol->computeBasis(u,v,w,p,spline,fromRight);
  dN.resize(spline.values.size());
  dir += 3*p-3;
  for (size_t i = 0; i < spline.values.size(); i++)
    dN[i] = spline.values[i][dir];
  */
}


short int ASMs3D::InterfaceChecker::hasContribution (int, int I, int J, int K) const
{
  bool neighbor[6];
  neighbor[0] = I > myPatch.svol->order(0);    // West neighbor
  neighbor[1] = I < myPatch.svol->numCoefs(0); // East neighbor
  neighbor[2] = J > myPatch.svol->order(1);    // South neighbor
  neighbor[3] = J < myPatch.svol->numCoefs(1); // North neighbor
  neighbor[4] = K > myPatch.svol->order(2);    // Back neighbor
  neighbor[5] = K < myPatch.svol->numCoefs(2); // Front neighbor

  // Check for existing neighbors
  short int status = 0, s = 1;
  for (short int i = 0; i < 6; i++, s *= 2)
    if (neighbor[i]) status += s;

  return status;
}


bool ASMs3D::evaluate (const FunctionBase* func, RealArray& vec,
                       int basisNum, double time) const
{
  const Go::SplineVolume* oldVol = this->getBasis(basisNum);
  Go::SplineVolume* newVol = SplineUtils::project(oldVol,*func,
                                                  func->dim(),time);
  if (!newVol)
  {
    std::cerr <<" *** ASMs3D::evaluate: Projection failure."<< std::endl;
    return false;
  }

  vec.assign(newVol->coefs_begin(),newVol->coefs_end());
  delete newVol;

  return true;
}


int ASMs3D::getCorner (int I, int J, int K, int basis) const
{
  int n1, n2, n3;
  int node = this->findStartNode(n1,n2,n3,basis);
  if (node < 1) return -1;

  if (swapW) // Account for swapped parameter direction
    K = -K;

  if (I > 0) node += n1-1;
  if (J > 0) node += n1*(n2-1);
  if (K > 0) node += n1*n2*(n3-1);

  return node;
}


bool ASMs3D::getFaceSize (int& n1, int& n2, int basis, int face) const
{
  int n3;
  if (!this->getSize(n1,n2,n3,basis))
    return false;

  if (face == 1 || face == 2)
    n1 = n2, n2 = n3;
  else if (face == 3 || face == 4)
    n2 = n3;
  else if (face != 5 && face != 6)
    n1 = n2 = 0;

  return true;
}


Field* ASMs3D::getProjectedField (const Vector& coefs) const
{
  if (this->getNoProjectionNodes() == coefs.size())
    return new SplineField3D(static_cast<const Go::SplineVolume*>(projB.get()),coefs);

  std::cerr <<" *** ASMs3D::getProjectedField: Non-matching coefficent array,"
            <<" size="<< coefs.size() <<" nnod="<< this->getNoProjectionNodes()
            << std::endl;
  return nullptr;
}


Fields* ASMs3D::getProjectedFields (const Vector& coefs, size_t) const
{
  if (projB.get() == this->getBasis(1) || this->getNoProjectionNodes() == 0)
    return nullptr;

  size_t ncmp = coefs.size() / this->getNoProjectionNodes();
  if (ncmp*this->getNoProjectionNodes() == coefs.size())
    return new SplineFields3D(static_cast<const Go::SplineVolume*>(projB.get()),coefs,ncmp);

  std::cerr <<" *** ASMs3D::getProjectedFields: Non-matching coefficent array,"
            <<" size="<< coefs.size() <<" nnod="<< this->getNoProjectionNodes()
            << std::endl;
  return nullptr;
}


size_t ASMs3D::getNoProjectionNodes () const
{
  if (!projB) return 0;

  const Go::SplineVolume* proj = static_cast<const Go::SplineVolume*>(projB.get());

  return proj->numCoefs(0) * proj->numCoefs(1) * proj->numCoefs(2);
}


void ASMs3D::getElmConnectivities (IntMat& neigh) const
{
  const int n1 = svol->numCoefs(0);
  const int n2 = svol->numCoefs(1);
  const int n3 = svol->numCoefs(2);
  const int p1 = svol->order(0);
  const int p2 = svol->order(1);
  const int p3 = svol->order(2);
  const int N1 = n1 - p1 + 1;
  const int N2 = n2 - p2 + 1;

  size_t iel = 0;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
        if (MLGE[iel] > 0)
        {
          int idx = MLGE[iel]-1;
          neigh[idx].resize(6,-1);
          if (i1 > p1)
            neigh[idx][0] = MLGE[iel-1]-1;
          if (i1 < n1)
            neigh[idx][1] = MLGE[iel+1]-1;
          if (i2 > p2)
            neigh[idx][2] = MLGE[iel-N1]-1;
          if (i2 < n2)
            neigh[idx][3] = MLGE[iel+N1]-1;
          if (i3 > p3)
            neigh[idx][4] = MLGE[iel-N1*N2]-1;
          if (i3 < n3)
            neigh[idx][5] = MLGE[iel+N1*N2]-1;
        }
}


void ASMs3D::findBoundaryElms (IntVec& elms, int lIndex, int) const
{
  int N1m = svol->numCoefs(0) - svol->order(0) + 1;
  int N2m = svol->numCoefs(1) - svol->order(1) + 1;
  int N3m = svol->numCoefs(2) - svol->order(2) + 1;

  elms.clear();
  switch (lIndex) {
  case 1:
  case 2:
    elms.reserve(N2m*N3m);
    for (int k = 0; k < N3m; ++k)
      for (int j = 0; j < N2m; ++j)
        elms.push_back(j*N1m + k*N1m*N2m + (lIndex-1)*(N1m-1));
    break;
  case 3:
  case 4:
    elms.reserve(N1m*N3m);
    for (int k = 0; k < N3m; ++k)
      for (int i = 0; i < N1m; ++i)
        elms.push_back(i + k*N1m*N2m + (lIndex-3)*(N1m*(N2m-1)));
    break;
  case 5:
  case 6:
    elms.reserve(N1m*N2m);
    for (int j = 0; j < N2m; ++j)
      for (int i = 0; i < N1m; ++i)
        elms.push_back(i + j*N1m + (lIndex-5)*(N1m*N2m*(N3m-1)));
  }
}


void ASMs3D::generateThreadGroupsFromElms (const IntVec& elms)
{
  myElms.clear();
  for (int elm : elms)
    if (this->getElmIndex(elm+1) > 0)
      myElms.push_back(this->getElmIndex(elm+1)-1);

  // We need myElms to be non-empty to flag that partitioning is enabled
  if (myElms.empty())
    myElms.push_back(-1);

  threadGroupsVol = threadGroupsVol.filter(myElms);

  for (std::pair<const char,ThreadGroups>& group : threadGroupsFace)
    group.second = group.second.filter(myElms);
}


bool ASMs3D::addRigidCpl (int lindx, int ldim, int basis,
                          int& gMaster, const Vec3& Xmaster, bool extraPt)
{
  if (threadGroupsVol.stripDir != ThreadGroups::NONE)
    IFEM::cout <<"  ** ASMs3D::addRigidCpl: Multi-threading deactivated"
               <<" for Patch "<< idx+1 << std::endl;
  threadGroupsVol.stripDir = ThreadGroups::NONE;
  threadGroupsFace[lindx].stripDir = ThreadGroups::NONE;

  return this->ASMstruct::addRigidCpl(lindx,ldim,basis,gMaster,Xmaster,extraPt);
}


ASMs3D::BasisFunctionCache::BasisFunctionCache (const ASMs3D& pch) :
  ::BasisFunctionCache<3>(), patch(pch)
{
  for (size_t d = 0; d < 3 && patch.svol; ++d)
    nel[d] = patch.svol->numCoefs(d) - patch.svol->order(d) + 1;
}


ASMs3D::BasisFunctionCache::BasisFunctionCache (const BasisFunctionCache& cache,
                                                int b) :
  ::BasisFunctionCache<3>(cache), patch(cache.patch), nel(cache.nel)
{
  basis = b;
}


bool ASMs3D::BasisFunctionCache::internalInit ()
{
  if (!mainQ->xg[0])
    this->setupQuadrature();

  nTotal = patch.nel * mainQ->ng[0]*mainQ->ng[1]*mainQ->ng[2];
  if (reducedQ->xg[0])
    nTotalRed = patch.nel * reducedQ->ng[0]*reducedQ->ng[1]*reducedQ->ng[2];

  return true;
}


void ASMs3D::BasisFunctionCache::internalCleanup ()
{
  if (basis == 1) {
    mainQ->reset();
    reducedQ->reset();
  }
}


bool ASMs3D::BasisFunctionCache::setupQuadrature ()
{
   int p[3];
   patch.getOrder(p[0],p[1],p[2]);

  // Get Gaussian quadrature points and weights
  for (int d = 0; d < 3; d++)
  {
    mainQ->ng[d] = patch.getNoGaussPt(p[d]);
    mainQ->xg[d] = GaussQuadrature::getCoord(mainQ->ng[d]);
    mainQ->wg[d] = GaussQuadrature::getWeight(mainQ->ng[d]);
    if (!mainQ->xg[d] || !mainQ->wg[d]) return false;
  }

  // Get the reduced integration quadrature points, if needed
  int nRed = integrand ? integrand->getReducedIntegration(mainQ->ng[0]) : 0;
  if (nRed > 0)
  {
    reducedQ->xg[0] = reducedQ->xg[1] = reducedQ->xg[2] = GaussQuadrature::getCoord(nRed);
    reducedQ->wg[0] = reducedQ->wg[1] = reducedQ->wg[2] = GaussQuadrature::getWeight(nRed);
    if (!reducedQ->xg[0] || !reducedQ->wg[0]) return false;
  }
   else if (nRed < 0)
    nRed = mainQ->ng[0];

  reducedQ->ng[0] = reducedQ->ng[1] = reducedQ->ng[2] = nRed;

  this->setupParameters();
  return true;
}


void ASMs3D::BasisFunctionCache::setupParameters ()
{
  // Compute parameter values of the Gauss points over the whole patch
  for (int d = 0; d < 3; d++) {
    patch.getGaussPointParameters(mainQ->gpar[d],d,mainQ->ng[d],mainQ->xg[d]);
    if (reducedQ->xg[0])
      patch.getGaussPointParameters(reducedQ->gpar[d],d,reducedQ->ng[d],reducedQ->xg[d]);
  }
}


BasisFunctionVals ASMs3D::BasisFunctionCache::calculatePt (size_t el,
                                                           size_t gp,
                                                           bool reduced) const
{
  PROFILE2("Spline evaluation");
  std::array<size_t,3> gpIdx = this->gpIndex(gp,reduced);
  std::array<size_t,3> elIdx = this->elmIndex(el);

  BasisFunctionVals result;
  if (nderiv == 1) {
    Go::BasisDerivs spline;
#pragma omp critical
    patch.getBasis(basis)->computeBasis(this->getParam(0,elIdx[0],gpIdx[0],reduced),
                                        this->getParam(1,elIdx[1],gpIdx[1],reduced),
                                        this->getParam(2,elIdx[2],gpIdx[2],reduced), spline);
    SplineUtils::extractBasis(spline,result.N,result.dNdu);
  } else if (nderiv == 2) {
    Go::BasisDerivs2 spline;
#pragma omp critical
    patch.getBasis(basis)->computeBasis(this->getParam(0,elIdx[0],gpIdx[0],reduced),
                                        this->getParam(1,elIdx[1],gpIdx[1],reduced),
                                        this->getParam(2,elIdx[2],gpIdx[2],reduced), spline);
    SplineUtils::extractBasis(spline,result.N,result.dNdu,result.d2Ndu2);
  }

  return result;
}


size_t ASMs3D::BasisFunctionCache::index (size_t el, size_t gp, bool reduced) const
{
  const Quadrature& q = reduced ? *reducedQ : *mainQ;
  std::array<size_t,3> elIdx = this->elmIndex(el);
  std::array<size_t,3> gpIdx = this->gpIndex(gp,reduced);

  return ((elIdx[2]*q.ng[2]*nel[1] + elIdx[1])*q.ng[1]*nel[0] + elIdx[0])*q.ng[0]
          + (gpIdx[2]*q.ng[1]*nel[1] + gpIdx[1])*q.ng[0]*nel[0] + gpIdx[0];
}


std::array<size_t, 3> ASMs3D::BasisFunctionCache::elmIndex (size_t el) const
{
  return { el % nel[0], (el/nel[0]) % nel[1], el / (nel[0]*nel[1]) };
}


void ASMs3D::BasisFunctionCache::calculateAll ()
{
  PROFILE2("Spline evaluation");
  auto&& extract1 = [this](const Matrix& par1,
                           const Matrix& par2,
                           const Matrix& par3,
                           std::vector<BasisFunctionVals>& values)
  {
    std::vector<Go::BasisDerivs>  spline;
    patch.getBasis(basis)->computeBasisGrid(par1,par2,par3,spline);
    size_t idx = 0;
    for (const Go::BasisDerivs& spl : spline) {
      SplineUtils::extractBasis(spl,values[idx].N,values[idx].dNdu);
      ++idx;
    }
  };
  if (nderiv == 1)
    extract1(mainQ->gpar[0],mainQ->gpar[1],mainQ->gpar[2],values);
  else if (nderiv == 2) {
    std::vector<Go::BasisDerivs2> spline;
    patch.getBasis(basis)->computeBasisGrid(mainQ->gpar[0],mainQ->gpar[1],mainQ->gpar[2],spline);
    size_t idx = 0;
    for (const Go::BasisDerivs2& spl : spline) {
      SplineUtils::extractBasis(spl,values[idx].N,values[idx].dNdu,values[idx].d2Ndu2);
      ++idx;
    }
  }
  if (!reducedQ->gpar[0].empty())
    extract1(reducedQ->gpar[0],reducedQ->gpar[1],reducedQ->gpar[2],valuesRed);
}
