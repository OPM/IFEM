// $Id$
//==============================================================================
//!
//! \file ASMs2D.C
//!
//! \date Jan 19 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Driver for assembly of structured 2D spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/CurveInterpolator.h"
#include "GoTools/geometry/RectDomain.h"

#include "ASMs2D.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Vec3Oper.h"
#include "Tensor.h"
#include "MPC.h"
#include <array>
#ifdef USE_OPENMP
#include <omp.h>
#endif


ASMs2D::ASMs2D (unsigned char n_s, unsigned char n_f)
  : ASMstruct(2,n_s,n_f), surf(nullptr), nodeInd(myNodeInd)
{
  bou[0] = bou[1] = bou[2] = bou[3] = nullptr;
  swapV = false;
}


ASMs2D::ASMs2D (const ASMs2D& patch, unsigned char n_f)
  : ASMstruct(patch,n_f), surf(patch.surf), nodeInd(patch.myNodeInd)
{
  for (int i = 0; i < 4; i++)
    bou[i] = patch.bou[i];
  swapV = patch.swapV;

  // Need to set nnod here,
  // as hasXNodes might be invoked before the FE data is generated
  if (nnod == 0 && surf)
    nnod = surf->numCoefs_u()*surf->numCoefs_v();
}


ASMs2D::ASMs2D (const ASMs2D& patch)
  : ASMstruct(patch), surf(patch.surf), nodeInd(myNodeInd)
{
  for (int i = 0; i < 4; i++)
    bou[i] = patch.bou[i];

  swapV = patch.swapV;
  myNodeInd = patch.nodeInd;
  dirich = patch.dirich;
}


ASMs2D::~ASMs2D ()
{
  if (!shareFE)
    for (int i = 0; i < 4; i++)
      delete bou[i];
}


Go::SplineCurve* ASMs2D::getBoundary (int dir, int)
{
  if (dir < -2 || dir == 0 || dir > 2)
    return nullptr;

  int iedge = dir > 0 ? dir : 3*dir+6;
  if (!bou[iedge])
    bou[iedge] = surf->edgeCurve(iedge);

  return bou[iedge];
}


void ASMs2D::copyParameterDomain (const ASMbase* other)
{
  const ASMs2D* o = dynamic_cast<const ASMs2D*>(other);
  if (!o) return;

  Go::RectDomain pd = o->getBasis()->parameterDomain();
  this->getBasis()->setParameterDomain(pd.umin(),pd.umax(),pd.vmin(),pd.vmax());
}


bool ASMs2D::read (std::istream& is)
{
  if (shareFE) return true;
  if (surf) delete surf;

  Go::ObjectHeader head;
  surf = new Go::SplineSurface;
  is >> head >> *surf;

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
    std::cerr <<" *** ASMs2D::read: Failure reading spline data"<< std::endl;
    delete surf;
    surf = 0;
    return false;
  }
  else if (surf->dimension() < 2)
  {
    std::cerr <<" *** ASMs2D::read: Invalid spline surface patch, dim="
	      << surf->dimension() << std::endl;
    delete surf;
    surf = 0;
    return false;
  }
  else if (surf->dimension() < nsd)
  {
    std::cerr <<"  ** ASMs2D::read: The dimension of this surface patch "
	      << surf->dimension() <<" is less than nsd="<< (int)nsd
	      <<".\n                   Resetting nsd to "<< surf->dimension()
	      <<" for this patch."<< std::endl;
    nsd = surf->dimension();
  }

  geo = surf;
  return true;
}


bool ASMs2D::write (std::ostream& os, int) const
{
  if (!surf) return false;

  os <<"200 1 0 0\n";
  os << *surf;

  return os.good();
}


void ASMs2D::clear (bool retainGeometry)
{
  if (!retainGeometry)
  {
    // Erase spline data
    if (surf && !shareFE) delete surf;
    surf = 0;
    geo = 0;
  }

  // Erase the FE data
  this->ASMbase::clear(retainGeometry);
  myNodeInd.clear();
  xnMap.clear();
  nxMap.clear();
}


bool ASMs2D::addXElms (short int dim, short int item, size_t nXn, IntVec& nodes)
{
  if (!this->addXNodes(dim,nXn,nodes))
    return false;

  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  int iel = 0;
  bool skipMe = false;
  for (int i2 = p2; i2 <= n2; i2++)
    for (int i1 = p1; i1 <= n1; i1++, iel++)
    {
      if (MLGE[iel] < 1) continue; // Skip zero-area element

      // Skip elements that are not on current boundary edge
      switch (item)
        {
        case 1: skipMe = i1 > p1; break;
        case 2: skipMe = i1 < n1; break;
        case 3: skipMe = i2 > p2; break;
        case 4: skipMe = i2 < n2; break;
        }
      if (skipMe) continue;

      IntVec& mnpc = myMNPC[nel+iel];
      if (!mnpc.empty())
      {
        std::cerr <<" *** ASMs2D::addXElms: Only one X-edge allowed."
                  << std::endl;
        return false;
      }

      mnpc = MNPC[iel]; // Copy the ordinary element nodes

      // Negate node numbers that are not on the boundary edge, to flag that
      // they shall not receive any tangent and/or residual contributions
      int lnod = 0;
      for (int j2 = 0; j2 < p2; j2++)
        for (int j1 = 0; j1 < p1; j1++, lnod++)
        {
          switch (item)
            {
            case 1: skipMe = j1 > 0;    break;
            case 2: skipMe = j1 < p1-1; break;
            case 3: skipMe = j2 > 0;    break;
            case 4: skipMe = j2 < p2-1; break;
            }
          if (skipMe) // Hack for node 0: Using -maxint as flag instead
            mnpc[lnod] = mnpc[lnod] == 0 ? -2147483648 : -mnpc[lnod];
        }

      // Add connectivity to the extra-ordinary nodes
      for (size_t i = 0; i < nXn; i++)
        mnpc.push_back(MLGN.size()-nXn+i);

      myMLGE[nel+iel] = -(++gEl); // Flag extraordinary element by negative sign
    }

  return true;
}


bool ASMs2D::addInterfaceElms (const InterfaceChecker& iChk)
{
  if (!surf || shareFE == 'F') return false;

  if (MNPC.size() != nel || MLGE.size() != nel)
  {
    // Already added extra elements, currently not allowed
    std::cerr <<" *** ASMs2D::addInterfaceElms: Already have extra elements."
              << std::endl;
    return false;
  }

  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  int iel = 0;
  for (int i2 = p2; i2 <= n2; i2++)
    for (int i1 = p1; i1 <= n1; i1++, iel++)
    {
      if (MLGE[iel] < 1) continue; // Skip zero-area element

      // Loop over the (north and east only) element edges with contributions
      short int status = iChk.hasContribution(i1,i2);
      for (int iedge = 1; iedge <= 4 && status > 0; iedge++, status /= 2)
        if (iedge%2 == 0 && status%2 == 1)
        {
          // Find index of the neighboring element
          int jel = iel + (iedge == 2 ? 1 : n1-p1+1);
          if (MLGE[jel] < 1) continue; // Skip zero-area element

          // Set up connectivity for the interface element
          IntVec mnpc(MNPC[iel]);
          utl::merge(mnpc,MNPC[jel]);
          myMNPC.push_back(mnpc);
          myMLGE.push_back(-(++gEl)); // Flag interface element by negative sign
        }
    }

  return true;
}


size_t ASMs2D::getNodeIndex (int globalNum, bool noAddedNodes) const
{
  IntVec::const_iterator it = std::find(MLGN.begin(),MLGN.end(),globalNum);
  if (it == MLGN.end()) return 0;

  size_t inod = 1 + (it-MLGN.begin());
  if (noAddedNodes && !xnMap.empty() && inod > nnod)
  {
    std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
    if (it != xnMap.end()) return it->second;
  }

  return inod;
}


int ASMs2D::getNodeID (size_t inod, bool noAddedNodes) const
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


bool ASMs2D::checkRightHandSystem ()
{
  if (!surf || shareFE) return false;

  // Evaluate the spline surface at its center
  RealArray u(1,0.5*(surf->startparam_u() + surf->endparam_u()));
  RealArray v(1,0.5*(surf->startparam_v() + surf->endparam_v()));
  RealArray X(3), dXdu(3), dXdv(3);
  surf->gridEvaluator(u,v,X,dXdu,dXdv);

  // Check that |J| = (dXdu x dXdv) * {0,0,1} > 0.0
  if (Vec3(dXdu,dXdv).z > 0.0) return false;

  // This patch has a negative Jacobian determinant. Probably it is modelled
  // in a left-hand-system. Swap the v-parameter direction to correct for this.
  surf->reverseParameterDirection(false);
  return swapV = true;
}


bool ASMs2D::refine (int dir, const RealArray& xi)
{
  if (!surf || dir < 0 || dir > 1 || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = surf->basis(dir).begin();
  double ucurr, uprev = *(uit++);
  while (uit != surf->basis(dir).end())
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

  if (dir == 0)
    surf->insertKnot_u(extraKnots);
  else
    surf->insertKnot_v(extraKnots);

  return true;
}


bool ASMs2D::uniformRefine (int dir, int nInsert)
{
  if (!surf || dir < 0 || dir > 1 || nInsert < 1) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = surf->basis(dir).begin();
  double ucurr, uprev = *(uit++);
  while (uit != surf->basis(dir).end())
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

  if (dir == 0)
    surf->insertKnot_u(extraKnots);
  else
    surf->insertKnot_v(extraKnots);

  return true;
}


bool ASMs2D::raiseOrder (int ru, int rv)
{
  if (!surf) return false;
  if (shareFE) return true;

  surf->raiseOrder(ru,rv);
  return true;
}


bool ASMs2D::generateFEMTopology ()
{
  if (!surf) return false;

  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  if (!nodeInd.empty())
  {
    nnod = n1*n2;
    if (nodeInd.size() != (size_t)nnod)
    {
      std::cerr <<" *** ASMs2D::generateFEMTopology: Inconsistency between the"
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
    nel = (n1-p1+1)*(n2-p2+1);
    return true;
  }
  else if (shareFE == 'F')
    return true;

#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1 <<" "<< n2;
  std::cout <<"\norder: "<< p1 <<" "<< p2;
  std::cout <<"\ndu:";
  for (int i = 0; i < n1; i++)
    std::cout <<' '<< surf->knotSpan(0,i);
  std::cout <<"\ndv:";
  for (int j = 0; j < n2; j++)
    std::cout <<' '<< surf->knotSpan(1,j);
  std::cout << std::endl;
#endif
  // Consistency checks, just to be fool-proof
  if (n1 <  2 || n2 <  2) return false;
  if (p1 <  1 || p1 <  1) return false;
  if (p1 > n1 || p2 > n2) return false;

  myMLGE.resize((n1-p1+1)*(n2-p2+1),0);
  myMLGN.resize(n1*n2);
  myMNPC.resize(myMLGE.size());
  myNodeInd.resize(myMLGN.size());

  nnod = nel = 0;
  for (int i2 = 1; i2 <= n2; i2++)
    for (int i1 = 1; i1 <= n1; i1++)
    {
      myNodeInd[nnod].I = i1-1;
      myNodeInd[nnod].J = i2-1;
      if (i1 >= p1 && i2 >= p2)
      {
        if (surf->knotSpan(0,i1-1) > 0.0)
          if (surf->knotSpan(1,i2-1) > 0.0)
          {
            myMLGE[nel] = ++gEl; // global element number over all patches
            myMNPC[nel].resize(p1*p2,0);

            int lnod = 0;
            for (int j2 = p2-1; j2 >= 0; j2--)
              for (int j1 = p1-1; j1 >= 0; j1--)
                myMNPC[nel][lnod++] = nnod - n1*j2 - j1;
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


int ASMs2D::Edge::next ()
{
  int ret = icnod;
  icnod += incr;

  return ret;
}


int ASMs2D::BlockNodes::next ()
{
  int ret = iinod;
  iinod += inc[0];

  if (++indxI >= nnodI-1)
  {
    indxI = 1;
    iinod += inc[1] - inc[0]*(nnodI-2);
  }

  return ret;
}


bool ASMs2D::assignNodeNumbers (BlockNodes& nodes, int basis)
{
  if (shareFE == 'F') return true;

  int n1, n2;
  if (!this->getSize(n1,n2,basis))
    return false;

  int m1 = 0, m2 = 0;
  if (basis > 0)
    if (!this->getSize(m1,m2,3-basis))
      return false;

  if (MLGN.size() != (size_t)(n1*n2+m1*m2)) return false;

  nodes.nnodI = n1;

  if (nodes.inc[0] == 0 || nodes.inc[1] == 0)
  {
    nodes.inc[0] = 1;
    nodes.inc[1] = n1-2;
  }

  int inod = basis > 1 ? m1*m2 : 0;
  for (int j = 1; j <= n2; j++)
    for (int i = 1; i <= n1; i++, inod++)
      if (j == 1)
      {
        if (i == 1)
          myMLGN[inod] = nodes.ibnod[0];
        else if (i == n1)
          myMLGN[inod] = nodes.ibnod[1];
        else
          myMLGN[inod] = nodes.edges[2].next();
      }
      else if (j == n2)
      {
        if (i == 1)
          myMLGN[inod] = nodes.ibnod[2];
        else if (i == n1)
          myMLGN[inod] = nodes.ibnod[3];
        else
          myMLGN[inod] = nodes.edges[3].next();
      }
      else
      {
        if (i == 1)
          myMLGN[inod] = nodes.edges[0].next();
        else if (i == n1)
          myMLGN[inod] = nodes.edges[1].next();
        else
          myMLGN[inod] = nodes.next();
      }

#if SP_DEBUG > 2
  if (basis > 0) std::cout <<"\nBasis "<< basis <<":";
  for (int i = inod-n1*n2; i < inod; i++)
  {
    std::cout <<"\nNode "<< i+1 <<"\t: ";
    if (!nodeInd.empty())
      std::cout << nodeInd[i].I <<" "<< nodeInd[i].J;
    std::cout <<"\tglobal no. "<< MLGN[i];
  }
  std::cout << std::endl;
#endif
  return true;
}


bool ASMs2D::connectPatch (int edge, ASMs2D& neighbor, int nedge,
                           bool revers, int, bool coordCheck, int thick)
{
  if (swapV && edge > 2) // Account for swapped parameter direction
    edge = 7-edge;

  if (neighbor.swapV && nedge > 2) // Account for swapped parameter direction
    nedge = 7-nedge;

  if (!this->connectBasis(edge,neighbor,nedge,revers,1,0,0,coordCheck,thick))
    return false;

  this->addNeighbor(&neighbor);
  return true;
}


bool ASMs2D::connectBasis (int edge, ASMs2D& neighbor, int nedge, bool revers,
                           int basis, int slave, int master,
                           bool coordCheck, int thick)
{
  if (this->isShared() && neighbor.isShared())
    return true;
  else if (this->isShared() || neighbor.isShared())
  {
    std::cerr <<" *** ASMs2D::connectBasis: Logic error, cannot"
	      <<" connect a shared patch with an unshared one"<< std::endl;
    return false;
  }

  // Set up the slave node numbers for this surface patch
  IntVec slaveNodes;
  this->getBoundaryNodes(edge, slaveNodes, basis, thick, true);
  for (int& it : slaveNodes)
    it += slave;

  // Set up the master node numbers for the neighboring surface patch
  IntVec masterNodes;
  neighbor.getBoundaryNodes(nedge, masterNodes, basis, thick, true);
  for (int& it : masterNodes)
    it += master;

  int n1, n2;
  if (!neighbor.getSize(n1,n2,basis)) return false;
  std::cout << "\tmaster coords:";
  for (int i=1; i<=(nedge<3 ? n1 : n1*n2); i+=(nedge<3 ? 1 : n1))
    std::cout << " " << neighbor.getCoord(i)[nedge<3 ? 0 : 1];
  std::cout << std::endl;

  if (!this->getSize(n1,n2,basis)) return false;
  std::cout << "\tslave coords:";
  for (int i=1; i<=(nedge<3 ? n1 : n1*n2); i+=(nedge<3 ? 1 : n1))
    std::cout << " " << this->getCoord(i)[nedge<3 ? 0 : 1];
  std::cout << std::endl;

  if (masterNodes.size() != slaveNodes.size())
  {
    std::cerr <<" *** ASMs2D::connectBasis: Non-matching edges, sizes "
              << masterNodes.size() <<" and "<< slaveNodes.size() << std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  for (size_t i = 0; i < masterNodes.size(); ++i)
  {
    int node = masterNodes[i];
    int slave = slaveNodes[revers ? slaveNodes.size()-i-1 : i];
    if (!coordCheck)
      ASMbase::collapseNodes(neighbor,node,*this,slave);
    else if (neighbor.getCoord(node).equal(this->getCoord(slave),xtol))
      ASMbase::collapseNodes(neighbor,node,*this,slave);
    else
    {
      std::cerr <<" *** ASMs2D::connectBasis: Non-matching nodes "
                << node <<": "<< neighbor.getCoord(node)
                <<"\n                                          and "
                << slave <<": "<< this->getCoord(slave) << std::endl;
      return false;
    }
  }

  return true;
}


void ASMs2D::closeEdges (int dir, int basis, int master)
{
  int n1, n2;
  if (basis < 1) basis = 1;
  if (!this->getSize(n1,n2,basis)) return;

  switch (dir)
    {
    case 1: // Edges are closed in I-direction
      for (int i2 = 1; i2 <= n2; i2++, master += n1)
	this->makePeriodic(master,master+n1-1);
      break;

    case 2: // Edges are closed in J-direction
      for (int i1 = 1; i1 <= n1; i1++, master++)
	this->makePeriodic(master,master+n1*(n2-1));
      break;
    }
}


/*!
  A negative \a code value implies direct evaluation of the Dirichlet condition
  function at the control point. Positive \a code implies projection onto the
  spline basis representing the boundary curve (needed for curved edges and/or
  non-constant functions).
*/

void ASMs2D::constrainEdge (int dir, bool open, int dof, int code, char basis)
{
  int n1, n2, node = 1;
  for (char i = 1; i <= basis; i++)
    if (!this->getSize(n1,n2,i))
      return;
    else if (i < basis)
      node += n1*n2;

  if (swapV) // Account for swapped parameter direction
    if (dir == 2 || dir == -2) dir = -dir;

  int bcode = code;
  if (code > 0) // Dirichlet projection will be performed
    dirich.push_back(DirichletEdge(this->getBoundary(dir,basis),dof,code));
  else if (code < 0)
    bcode = -code;

  switch (dir)
    {
    case  1: // Right edge (positive I-direction)
      node += n1-1;
    case -1: // Left edge (negative I-direction)
      if (!open)
	this->prescribe(node,dof,bcode);
      node += n1;
      for (int i2 = 2; i2 < n2; i2++, node += n1)
      {
	// If the Dirichlet condition is to be projected, add this node to
	// the set of nodes to receive prescribed value from the projection
	// **unless this node already has a homogeneous constraint**
	if (this->prescribe(node,dof,-code) == 0 && code > 0)
	  dirich.back().nodes.push_back(std::make_pair(i2,node));
      }
      if (!open)
	this->prescribe(node,dof,bcode);
      break;

    case  2: // Back edge (positive J-direction)
      node += n1*(n2-1);
    case -2: // Front edge (negative J-direction)
      if (!open)
	this->prescribe(node,dof,bcode);
      node++;
      for (int i1 = 2; i1 < n1; i1++, node++)
      {
	// If the Dirichlet condition is to be projected, add this node to
	// the set of nodes to receive prescribed value from the projection
	// **unless this node already has a homogeneous constraint**
	if (this->prescribe(node,dof,-code) == 0 && code > 0)
	  dirich.back().nodes.push_back(std::make_pair(i1,node));
      }
      if (!open)
	this->prescribe(node,dof,bcode);
      break;
    }

  if (code > 0)
    if (dirich.back().nodes.empty())
      dirich.pop_back(); // In the unlikely event of a 2-point boundary
#if SP_DEBUG > 1
    else
    {
      std::cout <<"Non-corner boundary nodes:";
      for (size_t i = 0; i < dirich.back().nodes.size(); i++)
	std::cout <<" ("<< dirich.back().nodes[i].first
		  <<","<< dirich.back().nodes[i].second
		  <<")";
      std::cout <<"\nThese nodes will be subjected to Dirichlet projection"
		<< std::endl;
    }
#endif
}


/*!
  The local coordinate systems in which the constraints are applied, are
  constructed from the tangent direction of the edge curve, evaluated at the
  Greville points (giving the local Y-direction). The local X-direction is then
  the outward-directed normal, defined from the cross product between the
  surface normal vector and the edge tangent. If \a project is \e true, the
  tangent- and normal direction vectors are projected onto the curve basis of
  the edge in order to obtain corresponding control point values. Otherwise,
  they are used directly.

  A negative \a code value implies direct evaluation of the Dirichlet condition
  function at the control point. Positive \a code implies projection onto the
  spline basis representing the boundary curve (needed for curved edges and/or
  non-constant functions).
*/

size_t ASMs2D::constrainEdgeLocal (int dir, bool open, int dof, int code,
				   bool project)
{
  if (shareFE == 'F')
  {
    std::cerr <<"\n *** ASMs2D::constrainEdgeLocal: Logic error, can not have"
	      <<" constraints in local CSs for shared patches."<< std::endl;
    return 0;
  }

  int ndir = abs(dir); // normal parameter direction (1,2) for the edge
  int tdir = 2 - ndir; // tangent parameter direction (0,1) for the edge
  if (swapV && tdir == 0) dir = -dir; // Account for swapped parameter direction
  double u[2] = { 0.0, 0.0 }; // running surface parameters along the edge
  if (ndir == 1)
    u[0] = dir < 0 ? surf->startparam_u() : surf->endparam_u();
  else if (ndir == 2)
    u[1] = dir < 0 ? surf->startparam_v() : surf->endparam_v();

  // Get parameter values of the Greville points along the edge
  RealArray gpar;
  if (!this->getGrevilleParameters(gpar,tdir))
    return 0;

  // Find the curve representing the edge geometry (for tangent evaluation)
  Go::SplineCurve* edge = this->getBoundary(dir);
  if (!edge) return 0;

  // We need to add extra nodes, check that the global node counter is good.
  // If not, we cannot do anything here
  if (gNod < *std::max_element(MLGN.begin(),MLGN.end()))
  {
    std::cerr <<"\n *** ASMs2D::constrainEdgeLocal: Logic error, gNod = "<< gNod
	      <<" is too small!"<< std::endl;
    return 0;
  }
  if (nf < nsd)
  {
    std::cerr <<"\n *** ASMs2D::constrainEdgeLocal: Not for scalar problems!"
	      << std::endl;
    return 0;
  }

  // Loop over the Greville points along the edge
  size_t i, k = 0;
  unsigned char c, d;
  std::vector<Go::Point> pts(3);
  RealArray gdata(6*gpar.size(),0.0);
  for (i = 0; i < gpar.size(); i++, k += 6)
  {
    // Find the tangent direction of the edge at this point.
    // That will be the y-axis of the local coordinate system along the edge.
    edge->point(pts,gpar[i],1);
    for (d = 0; d < nsd; d++)
      gdata[k+d] = dir == -1 || dir == 2 ? -pts[1][d] : pts[1][d];

    // Compute the surface normal at this edge point.
    // That will be the local z-axis of the local coordinate system.
    u[tdir] = gpar[i];
    surf->point(pts,u[0],u[1],1);
    Vec3 Zaxis(SplineUtils::toVec3(pts[1],nsd),SplineUtils::toVec3(pts[2],nsd));
    gdata[k+3] = Zaxis.x;
    gdata[k+4] = Zaxis.y;
    gdata[k+5] = Zaxis.z;
  }

  Go::SplineCurve* locc = nullptr;
  if (project)
  {
    // Project the Greville point values onto the spline basis
    // to obtain the corresponding control point values

    RealArray weights;
    if (edge->rational())
      edge->getWeights(weights);

    locc = Go::CurveInterpolator::regularInterpolation(edge->basis(),
                                                       gpar,gdata,6,
                                                       edge->rational(),
                                                       weights);
  }

  // Find start index and increment of the (slave) node with the global DOFs
  int n1, n2, incNod = 0, iSnod = 0;
  this->getSize(n1,n2,1);
  switch (dir)
    {
    case  1: // Right edge (positive I-direction)
      iSnod = n1-1;
    case -1: // Left edge (negative I-direction)
      incNod = n1;
      break;
    case  2: // Back edge (positive J-direction)
      iSnod = n1*(n2-1);
    case -2: // Front edge (negative J-direction)
      incNod = 1;
      break;
    }

  int bcode = code;
  if (code > 0) // Dirichlet projection will be performed
    dirich.push_back(DirichletEdge(edge,dof,code));
  else if (code < 0)
    bcode = -code;

  size_t nxNode = 0;
  RealArray::const_iterator it = locc ? locc->coefs_begin() : gdata.begin();
  for (i = 0; i < gpar.size(); i++, iSnod += incNod, it += 6)
  {
    // Skip the end points if this should be regarded an open boundary
    if (open && (i == 0 || i+1 == gpar.size())) continue;
    // Check if this node already has been constrained or fixed
    if (this->isFixed(MLGN[iSnod],dof)) continue;

    // We need an extra node representing the local (master) DOFs at this point
    int iMnod = myMLGN.size();
    std::map<int,int>::const_iterator xit = xNode.end();
    // Create an extra node for the local DOFs. The new node, for which
    // the Dirichlet boundary conditions will be defined, then inherits
    // the global node number of the original node. The original node, which
    // do not enter the equation system, receives a new global node number.
    if (i > 0 && i+1 < gpar.size())
      // This is not a corner node
      myMLGN.push_back(++gNod);
    else if ((xit = xNode.find(MLGN[iSnod])) != xNode.end())
      // This is a corner node already processed by another patch
      myMLGN.push_back(xit->second);
    else
    {
      // This is a corner node, store its original-to-extra node number
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
    if (i == 0 || i+1 == gpar.size())
      this->prescribe(1+iMnod,dof,bcode);
    else
    {
      this->prescribe(1+iMnod,dof,-code);
      if (code > 0)
        dirich.back().nodes.push_back(std::make_pair(1+i,1+iMnod));
    }

    // Find the local axis directions of the edge at this point
    Vec3 Yaxis(it[0],it[1],it[2]);
    Vec3 Zaxis(it[3],it[4],it[5]);
    Vec3 Xaxis(Yaxis,Zaxis);
    Yaxis.normalize();
    Xaxis.normalize();

    // Local-to-global transformation matrix at this point
    Tensor Tlg(Xaxis,Yaxis,Zaxis.cross(Xaxis,Yaxis));

    // Now establish constraint equations relating the global and local DOFs.
    // We here assume there are (at least) nsd unknowns per node,
    // and only the first nsd DOFs are subjected to transformation.
    for (d = 1; d <= nf; d++)
    {
      MPC* cons = new MPC(slaveNode,d);
      if (this->addMPC(cons,0,true) && cons)
      {
        if (d > nsd)
          cons->addMaster(masterNode,d);
        else for (c = 1; c <= nsd; c++)
          if (!this->isFixed(masterNode,c))
            cons->addMaster(masterNode,c,Tlg(d,c));
#if SP_DEBUG > 1
        std::cout <<"Added constraint: "<< *cons;
#endif
      }
    }
  }

  if (locc) delete locc;
  return nxNode; // Number of added nodes
}


void ASMs2D::constrainCorner (int I, int J, int dof, int code, char basis)
{
  int node = this->getCorner(I, J, basis);
  if (node > 0)
    this->prescribe(node,dof,code);
}


void ASMs2D::constrainNode (double xi, double eta, int dof, int code)
{
  if (xi  < 0.0 || xi  > 1.0) return;
  if (eta < 0.0 || eta > 1.0) return;

  if (swapV) // Account for swapped parameter direction
    eta = 1.0-eta;

  int n1, n2;
  if (!this->getSize(n1,n2,1)) return;

  int I = int(0.5+(n1-1)*xi);
  int J = int(0.5+(n2-1)*eta);

  this->prescribe(n1*J+I+1,dof,code);
}


void ASMs2D::setNodeNumbers (const std::vector<int>& nodes)
{
  this->ASMbase::setNodeNumbers(nodes);
  if (!swapV) return;

  // Account for swapped parameter direction
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  for (int j = 0; j < n2/2; j++)
    for (int i = 0; i < n1; i++)
      std::swap(myMLGN[i+n1*j],myMLGN[i+n1*(n2-j-1)]);
}


/*!
  This method projects the function describing the in-homogeneous Dirichlet
  boundary condition onto the spline basis defining the boundary curve,
  in order to find the control point values which are used as the prescribed
  values of the boundary DOFs.
*/

bool ASMs2D::updateDirichlet (const std::map<int,RealFunc*>& func,
                              const std::map<int,VecFunc*>& vfunc, double time,
                              const std::map<int,int>* g2l)
{
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;
  std::vector<DirichletEdge>::const_iterator dit;
  std::vector<Ipair>::const_iterator nit;

  for (size_t i = 0; i < dirich.size(); i++)
  {
    // Project the function onto the spline curve basis
    Go::SplineCurve* dcrv = 0;
    if ((fit = func.find(dirich[i].code)) != func.end())
      dcrv = SplineUtils::project(dirich[i].curve,*fit->second,time);
    else if ((vfit = vfunc.find(dirich[i].code)) != vfunc.end())
      dcrv = SplineUtils::project(dirich[i].curve,*vfit->second,nf,time);
    else
    {
      std::cerr <<" *** ASMs2D::updateDirichlet: Code "<< dirich[i].code
		<<" is not associated with any function."<< std::endl;
      return false;
    }
    if (!dcrv)
    {
      std::cerr <<" *** ASMs2D::updateDirichlet: Projection failure."
		<< std::endl;
      return false;
    }

    // Loop over the (interior) nodes (control points) of this boundary curve
    for (nit = dirich[i].nodes.begin(); nit != dirich[i].nodes.end(); ++nit)
      for (int dofs = dirich[i].dof; dofs > 0; dofs /= 10)
      {
	int dof = dofs%10;
	// Find the constraint equation for current (node,dof)
	MPC pDOF(MLGN[nit->second-1],dof);
	MPCIter mit = mpcs.find(&pDOF);
	if (mit == mpcs.end()) continue; // probably a deleted constraint

	// Find index to the control point value for this (node,dof) in dcrv
	RealArray::const_iterator cit = dcrv->coefs_begin();
	if (dcrv->dimension() > 1) // A vector field is specified
	  cit += (nit->first-1)*dcrv->dimension() + (dof-1);
	else // A scalar field is specified at this dof
	  cit += (nit->first-1);

	// Now update the prescribed value in the constraint equation
	(*mit)->setSlaveCoeff(*cit);
#if SP_DEBUG > 1
	std::cout <<"Updated constraint: "<< **mit;
#endif
      }
    delete dcrv;
  }

  // The parent class method takes care of the corner nodes with direct
  // evaluation of the Dirichlet functions (since they are interpolatory)
  return this->ASMbase::updateDirichlet(func,vfunc,time,g2l);
}


#define DERR -999.99

double ASMs2D::getParametricArea (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs2D::getParametricArea: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  int inod1 = MNPC[iel-1][surf->order_u()*surf->order_v()-1];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs2D::getParametricArea: Node index "<< inod1
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
    return DERR;
  }
#endif

  double du = surf->knotSpan(0,nodeInd[inod1].I);
  double dv = surf->knotSpan(1,nodeInd[inod1].J);

  return du*dv;
}


double ASMs2D::getParametricLength (int iel, int dir) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs2D::getParametricLength: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  int inod1 = MNPC[iel-1][surf->order_u()*surf->order_v()-1];
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nnod)
  {
    std::cerr <<" *** ASMs2D::getParametricLength: Node index "<< inod1
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
    return DERR;
  }
#endif

  switch (dir)
    {
    case 1: return surf->knotSpan(0,nodeInd[inod1].I);
    case 2: return surf->knotSpan(1,nodeInd[inod1].J);
    }

  std::cerr <<" *** ASMs2D::getParametricLength: Invalid edge direction "
	    << dir << std::endl;
  return DERR;
}


int ASMs2D::coeffInd (size_t inod) const
{
#ifdef INDEX_CHECK
  if (inod >= nnod)
  {
    std::cerr <<" *** ASMs2D::coeffInd: Node index "<< inod
	      <<" out of range [0,"<< nnod <<">."<< std::endl;
    return -1;
  }
#endif

  const int ni = nodeInd[inod].I;
  const int nj = nodeInd[inod].J;
  return nj*surf->numCoefs_u() + ni;
}


Vec3 ASMs2D::getCoord (size_t inod) const
{
  if (inod > nnod && inod <= MLGN.size())
  {
    // This is a node added due to constraints in local directions.
    // Find the corresponding original node (see constrainEdgeLocal)
    std::map<size_t,size_t>::const_iterator it = xnMap.find(inod);
    if (it != xnMap.end()) inod = it->second;
  }

  int ip = this->coeffInd(inod-1)*surf->dimension();
  if (ip < 0) return Vec3();

  return Vec3(&(*(surf->coefs_begin()+ip)),nsd);
}


bool ASMs2D::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs2D::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  X.resize(nsd,surf->order_u()*surf->order_v());

  RealArray::const_iterator cit = surf->coefs_begin();
  for (size_t n = 0; n < X.cols(); n++)
  {
    int ip = this->coeffInd(MNPC[iel-1][n])*surf->dimension();
    if (ip < 0) return false;

    for (size_t i = 0; i < nsd; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


void ASMs2D::getNodalCoordinates (Matrix& X) const
{
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  X.resize(nsd,n1*n2);

  RealArray::const_iterator cit = surf->coefs_begin();
  size_t inod = 1;
  for (int i2 = 0; i2 < n2; i2++)
    for (int i1 = 0; i1 < n1; i1++, inod++)
    {
      int ip = (i2*n1 + i1)*surf->dimension();
      for (size_t i = 0; i < nsd; i++)
	X(i+1,inod) = *(cit+(ip+i));
    }
}


bool ASMs2D::updateCoords (const Vector& displ)
{
  if (!surf) return true; // silently ignore empty patches
  if (shareFE) return true;

  if (displ.size() != nsd*MLGN.size())
  {
    std::cerr <<" *** ASMs2D::updateCoords: Invalid dimension "
	      << displ.size() <<" on displacement vector, should be "
	      << nsd*MLGN.size() << std::endl;
    return false;
  }

  surf->deform(displ,nsd);
  return true;
}


void ASMs2D::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                               int thick, bool local) const
{
  if (basis == 0)
    basis = 1;

  if (!this->getBasis(basis)) return; // silently ignore empty patches

#if SP_DEBUG > 1
  size_t last = nodes.size();
#endif

  int n1, n2, node = 1;
  for (char i = 1; i <= basis; i++)
    if (!this->getSize(n1,n2,i))
      return;
    else if (i < basis && !local)
      node += n1*n2;

  switch (lIndex)
    {
    case  2: // Right edge (positive I-direction)
      node += n1-thick;
    case 1: // Left edge (negative I-direction)
      for (int i2 = 1; i2 <= n2; i2++, node += n1)
        for (int t = 0; t < thick; ++t)
          nodes.push_back(local ? node+t : this->getNodeID(node+t));
      break;

    case  4: // Back edge (positive J-direction)
      node += n1*(n2-thick);
    case  3: // Front edge (negative J-direction)
      for (int i1 = 1; i1 <= n1; i1++, node++)
        for (int t = 0; t < thick; ++t)
          nodes.push_back(local ? node + t*n1 : this->getNodeID(node + t*n1));
      break;
    }

#if SP_DEBUG > 1
  std::cout <<"Boundary nodes in patch "<< idx+1 <<" edge "<< lIndex <<":";
  if (nodes.size() == last)
    std::cout <<" (none)";
  else for (size_t i = last; i < nodes.size(); i++)
    std::cout <<" "<< nodes[i];
  std::cout << std::endl;
#endif
}


bool ASMs2D::getOrder (int& p1, int& p2) const
{
  if (!surf) return false;

  p1 = surf->order_u();
  p2 = surf->order_v();
  return true;
}


bool ASMs2D::getOrder (int& p1, int& p2, int& p3) const
{
  p3 = 0;
  return this->getOrder(p1,p2);
}


bool ASMs2D::getSize (int& n1, int& n2, int& n3, int basis) const
{
  n3 = 0;
  return this->getSize(n1,n2,basis);
}


bool ASMs2D::getSize (int& n1, int& n2, int) const
{
  if (!surf) return false;

  n1 = surf->numCoefs_u();
  n2 = surf->numCoefs_v();
  return true;
}


size_t ASMs2D::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (ldim < 1 && lIndex > 0)
    return 1;

  switch (lIndex)
    {
    case 1:
    case 2:
      return surf->numCoefs_v() - surf->order_v() + 1;
    case 3:
    case 4:
      return surf->numCoefs_u() - surf->order_u() + 1;
    }

  return 0;
}


const Vector& ASMs2D::getGaussPointParameters (Matrix& uGP, int dir, int nGauss,
					       const double* xi) const
{
  int pm1 = (dir == 0 ? surf->order_u() : surf->order_v()) - 1;
  RealArray::const_iterator uit = surf->basis(dir).begin() + pm1;

  int nCol = (dir == 0 ? surf->numCoefs_u() : surf->numCoefs_v()) - pm1;
  uGP.resize(nGauss,nCol);

  double ucurr, uprev = *(uit++);
  for (int j = 1; j <= nCol; ++uit, j++)
  {
    ucurr = *uit;
    for (int i = 1; i <= nGauss; i++)
      uGP(i,j) = 0.5*((ucurr-uprev)*xi[i-1] + ucurr+uprev);
    uprev = ucurr;
  }

  return uGP;
}


void ASMs2D::getElementBorders (int i1, int i2, double* u, double* v) const
{
  RealArray::const_iterator uit = surf->basis(0).begin();
  RealArray::const_iterator vit = surf->basis(1).begin();
  for (int i = 0; i < 2; i++)
  {
    u[i] = uit[i1+i];
    v[i] = vit[i2+i];
  }
}


void ASMs2D::getElementCorners (int i1, int i2, Vec3Vec& XC) const
{
  // Fetch parameter values of the element (knot-span) corners
  RealArray u(2), v(2);
  this->getElementBorders(i1,i2,u.data(),v.data());

  // Evaluate the spline surface at the corners to find physical coordinates
  int dim = surf->dimension();
  RealArray XYZ(dim*4);
#pragma omp critical
  surf->gridEvaluator(XYZ,u,v);

  XC.clear();
  XC.reserve(4);
  const double* pt = &XYZ.front();
  for (int i = 0; i < 4; i++, pt += dim)
    XC.push_back(Vec3(pt,dim));
}


#if SP_DEBUG > 4
std::ostream& operator<< (std::ostream& os, const Go::BasisDerivsSf& bder)
{
  os <<" : (u,v) = "<< bder.param[0] <<" "<< bder.param[1]
     <<"  left_idx = "<< bder.left_idx[0] <<" "<< bder.left_idx[1] << std::endl;
  for (size_t i = 0; i < bder.basisValues.size(); i++)
    os << 1+i <<'\t'<< bder.basisValues[i] <<'\t'
       << bder.basisDerivs_u[i] <<'\t'<< bder.basisDerivs_v[i] << std::endl;
  return os;
}
#endif


bool ASMs2D::integrate (Integrand& integrand,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
  if (!surf) return true; // silently ignore empty patches

  PROFILE2("ASMs2D::integrate(I)");

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Get the reduced integration quadrature points, if needed
  const double* xr = nullptr;
  const double* wr = nullptr;
  int nRed = integrand.getReducedIntegration(nGauss);
  if (nRed > 0)
  {
    xr = GaussQuadrature::getCoord(nRed);
    wr = GaussQuadrature::getWeight(nRed);
    if (!xr && !wr) return false;
  }
  else if (nRed < 0)
    nRed = nGauss; // The integrand needs to know nGauss

  // Compute parameter values of the Gauss points over the whole patch
  std::array<Matrix,2> gpar, redpar;
  for (int d = 0; d < 2; d++)
  {
    this->getGaussPointParameters(gpar[d],d,nGauss,xg);
    if (xr)
      this->getGaussPointParameters(redpar[d],d,nRed,xr);
  }

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivsSf>  spline;
  std::vector<Go::BasisDerivsSf2> spline2;
  std::vector<Go::BasisDerivsSf>  splineRed;
  if (use2ndDer)
    surf->computeBasisGrid(gpar[0],gpar[1],spline2);
  else
    surf->computeBasisGrid(gpar[0],gpar[1],spline);
  if (xr)
    surf->computeBasisGrid(redpar[0],redpar[1],splineRed);

#if SP_DEBUG > 4
  for (size_t i = 0; i < spline.size(); i++)
    std::cout <<"\nBasis functions at integration point "<< 1+i << spline[i];
#endif

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int nel1 = n1 - p1 + 1;


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroups.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGroups[g].size(); t++)
    {
      FiniteElement fe(p1*p2);
      Matrix   dNdu, Xnod, Jac;
      Matrix3D d2Ndu2, Hess;
      double   dXidu[2];
      Vec4     X;
      for (size_t i = 0; i < threadGroups[g][t].size() && ok; i++)
      {
        int iel = threadGroups[g][t][i];
        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-area element

#ifdef SP_DEBUG
        if (dbgElm < 0 && 1+iel != -dbgElm)
          continue; // Skipping all elements, except for -dbgElm
#endif

        int i1 = p1 + iel % nel1;
        int i2 = p2 + iel / nel1;

        // Get element area in the parameter space
        double dA = 0.25*this->getParametricArea(++iel);
        if (dA < 0.0) // topology error (probably logic error)
        {
          ok = false;
          break;
        }

        // Set up control point (nodal) coordinates for current element
        if (!this->getElementCoordinates(Xnod,iel))
        {
          ok = false;
          break;
        }

        if (useElmVtx)
          this->getElementCorners(i1-1,i2-1,fe.XC);

        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        {
          // Element size in parametric space
          dXidu[0] = surf->knotSpan(0,i1-1);
          dXidu[1] = surf->knotSpan(1,i2-1);
        }

        if (integrand.getIntegrandType() & Integrand::AVERAGE)
        {
          // --- Compute average value of basis functions over the element -----

          fe.Navg.resize(p1*p2,true);
          double area = 0.0;
          int ip = ((i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
          for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
            for (int i = 0; i < nGauss; i++, ip++)
            {
              // Fetch basis function derivatives at current integration point
              SplineUtils::extractBasis(spline[ip],fe.N,dNdu);

              // Compute Jacobian determinant of coordinate mapping
              // and multiply by weight of current integration point
              double detJac = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu,false);
              double weight = dA*wg[i]*wg[j];

              // Numerical quadrature
              fe.Navg.add(fe.N,detJac*weight);
              area += detJac*weight;
            }

          // Divide by element area
          fe.Navg /= area;
        }

        if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
        {
          // Compute the element center
          double u0 = 0.5*(gpar[0](1,i1-p1+1) + gpar[0](nGauss,i1-p1+1));
          double v0 = 0.5*(gpar[1](1,i2-p2+1) + gpar[1](nGauss,i2-p2+1));
          SplineUtils::point(X,u0,v0,surf);
          if (!useElmVtx)
          {
            // When element corner coordinates are not needed, store coordinates
            // and parameters of the element center in XC, for material usage
            fe.XC.resize(2);
            fe.XC.front() = X;
            fe.XC.back() = Vec3(u0,v0,0.0);
          }
        }

        // Initialize element quantities
        LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
        if (!integrand.initElement(MNPC[iel-1],fe,X,nRed*nRed,*A))
        {
          A->destruct();
          ok = false;
          break;
        }

        if (xr)
        {
          // --- Selective reduced integration loop ----------------------------

          int ip = ((i2-p2)*nRed*nel1 + i1-p1)*nRed;
          for (int j = 0; j < nRed; j++, ip += nRed*(nel1-1))
            for (int i = 0; i < nRed; i++, ip++)
            {
              // Local element coordinates of current integration point
              fe.xi  = xr[i];
              fe.eta = xr[j];

              // Parameter values of current integration point
              fe.u = redpar[0](i+1,i1-p1+1);
              fe.v = redpar[1](j+1,i2-p2+1);

              // Fetch basis function derivatives at current point
              SplineUtils::extractBasis(splineRed[ip],fe.N,dNdu);

              // Compute Jacobian inverse and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

              // Cartesian coordinates of current integration point
              X = Xnod * fe.N;
              X.t = time.t;

              // Compute the reduced integration terms of the integrand
              fe.detJxW *= dA*wr[i]*wr[j];
              if (!integrand.reducedInt(*A,fe,X))
                ok = false;
            }
        }


        // --- Integration loop over all Gauss points in each direction --------

        int ip = ((i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
        int jp = ((i2-p2)*nel1 + i1-p1)*nGauss*nGauss;
        fe.iGP = firstIp + jp; // Global integration point counter

        for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
          for (int i = 0; i < nGauss; i++, ip++, fe.iGP++)
          {
            // Local element coordinates of current integration point
            fe.xi  = xg[i];
            fe.eta = xg[j];

            // Parameter values of current integration point
            fe.u = gpar[0](i+1,i1-p1+1);
            fe.v = gpar[1](j+1,i2-p2+1);

            // Fetch basis function derivatives at current integration point
            if (use2ndDer)
              SplineUtils::extractBasis(spline2[ip],fe.N,dNdu,d2Ndu2);
            else
              SplineUtils::extractBasis(spline[ip],fe.N,dNdu);

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
            {
              std::cout <<"\niel, ip = "<< iel <<" "<< ip
                        <<"\nN ="<< fe.N <<"dNdX ="<< fe.dNdX;
              if (!fe.d2NdX2.empty())
                std::cout <<"d2NdX2 ="<< fe.d2NdX2;
            }
#endif

            // Cartesian coordinates of current integration point
            X = Xnod * fe.N;
            X.t = time.t;

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= dA*wg[i]*wg[j];
#ifndef USE_OPENMP
            PROFILE3("Integrand::evalInt");
#endif
            if (!integrand.evalInt(*A,fe,time,X))
              ok = false;
          }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,time,firstIp+jp))
          ok = false;

        // Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();

#ifdef SP_DEBUG
	if (iel == -dbgElm) break; // Skipping all elements, except for -dbgElm
#endif
      }
    }
  }

  return ok;
}


bool ASMs2D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time,
                        const Real3DMat& itgPts)
{
  if (!surf) return true; // silently ignore empty patches

  if (integrand.getReducedIntegration(nGauss) != 0)
  {
    std::cerr <<" *** ASMs2D::integrate(Integrand&,GlobalIntegral&,"
              <<"const TimeDomain&,const Real3DMat&): Available for standard"
              <<" integrands only."<< std::endl;
    return false;
  }

  PROFILE2("ASMs2D::integrate(I)");

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  bool useElmVtx = integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS;

  // Evaluate basis function derivatives at all integration points
  size_t i, j, k;
  std::vector<size_t> MPitg(itgPts.size()+1,0);
  for (i = MPitg.front() = 0; i < itgPts.size(); i++)
    MPitg[i+1] = MPitg[i] + itgPts[i].size();
  size_t nPoints = MPitg.back();
  std::vector<Go::BasisDerivsSf>  spline(use2ndDer ? 0 : nPoints);
  std::vector<Go::BasisDerivsSf2> spline2(!use2ndDer ? 0 : nPoints);
  for (i = k = 0; i < itgPts.size(); i++)
    for (j = 0; j < itgPts[i].size(); j++, k++)
      if (use2ndDer)
        surf->computeBasis(itgPts[i][j][0],itgPts[i][j][1],spline2[k]);
      else
        surf->computeBasis(itgPts[i][j][0],itgPts[i][j][1],spline[k]);

#if SP_DEBUG > 4
  for (i = 0; i < spline.size(); i++)
    std::cout <<"\nBasis functions at integration point "<< 1+i << spline[i];
#endif

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int nel1 = surf->numCoefs_u() - p1 + 1;


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t g = 0; g < threadGroups.size() && ok; g++)
  {
#pragma omp parallel for schedule(static)
    for (size_t t = 0; t < threadGroups[g].size(); t++)
    {
      FiniteElement fe(p1*p2);
      Matrix   dNdu, Xnod, Jac;
      Matrix3D d2Ndu2, Hess;
      double   dXidu[2];
      Vec4     X;
      for (size_t e = 0; e < threadGroups[g][t].size() && ok; e++)
      {
        int iel = threadGroups[g][t][e];
        if (itgPts[iel].empty()) continue; // no points in this element

        fe.iel = MLGE[iel];
        if (fe.iel < 1) continue; // zero-area element

#ifdef SP_DEBUG
        if (dbgElm < 0 && iel != -dbgElm)
          continue; // Skipping all elements, except for -dbgElm
#endif

        int i1 = p1 + iel % nel1;
        int i2 = p2 + iel / nel1;

        // Get element area in the parameter space
        double dA = 0.25*this->getParametricArea(++iel);
        if (dA < 0.0) // topology error (probably logic error)
        {
          ok = false;
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
          this->getElementCorners(i1-1,i2-1,fe.XC);
          X = 0.25*(fe.XC[0]+fe.XC[1]+fe.XC[2]+fe.XC[3]);
        }
        else if (useElmVtx)
          this->getElementCorners(i1-1,i2-1,fe.XC);

        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        {
          // Element size in parametric space
          dXidu[0] = surf->knotSpan(0,i1-1);
          dXidu[1] = surf->knotSpan(1,i2-1);
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

        size_t jp = MPitg[iel-1]; // Patch-wise integration point counter
        fe.iGP = firstIp + jp;    // Global integration point counter

        const Real2DMat& elmPts = itgPts[iel-1]; // points for current element
        for (size_t ip = 0; ip < elmPts.size(); ip++, jp++, fe.iGP++)
        {
          // Parameter values of current integration point
          fe.u = elmPts[ip][0];
          fe.v = elmPts[ip][1];

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
            std::cout <<"\niel, jp = "<< iel <<" "<< jp
                      <<"\nN ="<< fe.N <<"dNdX ="<< fe.dNdX;
#endif

          // Cartesian coordinates of current integration point
          X = Xnod * fe.N;
          X.t = time.t;

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= dA*elmPts[ip][2];
#ifndef USE_OPENMP
          PROFILE3("Integrand::evalInt");
#endif
          if (!integrand.evalInt(*A,fe,time,X))
            ok = false;
        }

        // Finalize the element quantities
        if (ok && !integrand.finalizeElement(*A,time,firstIp+MPitg[iel]))
          ok = false;

        // Assembly of global system integral
        if (ok && !glInt.assemble(A->ref(),fe.iel))
          ok = false;

        A->destruct();

#ifdef SP_DEBUG
	if (iel == -dbgElm) break; // Skipping all elements, except for -dbgElm
#endif
      }
    }
  }

  return ok;
}


bool ASMs2D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time,
                        const InterfaceChecker& iChk)
{
  if (!surf) return true; // silently ignore empty patches
  if (!(integrand.getIntegrandType() & Integrand::INTERFACE_TERMS)) return true;

  PROFILE2("ASMs2D::integrate(J)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  FiniteElement fe(p1*p2);
  Matrix        dNdu, Xnod, Jac;
  Vector        dN;
  Vec4          X;
  Vec3          normal;
  double        u[2], v[2];
  bool          hasInterfaceElms = MLGE.size() > nel && MLGE.size() != 2*nel;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 0, jel = nel;
  for (int i2 = p2; i2 <= n2; i2++)
    for (int i1 = p1; i1 <= n1; i1++, iel++)
    {
      if (!hasInterfaceElms) jel = iel;

      fe.iel = abs(MLGE[jel]);
      if (fe.iel < 1) continue; // zero-area element

      short int status = iChk.hasContribution(i1,i2);
      if (!status) continue; // no interface contributions for this element

#if SP_DEBUG > 3
      std::cout <<"\n\nIntegrating interface terms for element "<< fe.iel
                << std::endl;
#endif

      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,1+iel)) return false;

      // Compute parameter values of the element edges
      this->getElementBorders(i1-1,i2-1,u,v);

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
        this->getElementCorners(i1-1,i2-1,fe.XC);

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(MNPC[jel].size(),fe.iel);
      bool ok = hasInterfaceElms ? true : integrand.initElement(MNPC[jel],*A);

      // Loop over the element edges with contributions
      for (int iedge = 1; iedge <= 4 && status > 0 && ok; iedge++, status /= 2)
        if (status%2 == 1)
        {
          // Find the parametric direction of the edge normal {-2,-1, 1, 2}
          const int edgeDir = (iedge+1)/(iedge%2 ? -2 : 2);
          const int t1 = abs(edgeDir);   // Tangent direction normal to the edge
          const int t2 = 3-abs(edgeDir); // Tangent direction along the edge

          // Get element edge length in the parameter space
          double dS = 0.5*this->getParametricLength(1+iel,t2);
          if (dS < 0.0) // topology error (probably logic error)
            ok = false;
          else if (hasInterfaceElms) // Initialize the interface element
            ok = integrand.initElement(MNPC[jel++],*A);

          // Find index of the neighboring element
          int kel = iel + (t1 == 1 ? 1 : n1-p1+1);


          // --- Integration loop over all Gauss points along the edge ---------

          for (int i = 0; i < nGauss && ok; i++)
          {
            // Local element coordinates and parameter values
            // of current integration point
            if (t1 == 1)
            {
              fe.xi = edgeDir;
              fe.eta = xg[i];
              fe.u = edgeDir > 0 ? u[1] : u[0];
              fe.v = 0.5*((v[1]-v[0])*xg[i] + v[1]+v[0]);
              fe.p = p1 - 1;
            }
            else
            {
              fe.xi = xg[i];
              fe.eta = edgeDir/2;
              fe.u = 0.5*((u[1]-u[0])*xg[i] + u[1]+u[0]);
              fe.v = edgeDir > 0 ? v[1] : v[0];
              fe.p = p2 - 1;
            }

            // Fetch basis function derivatives at current integration point
            this->extractBasis(fe.u,fe.v, fe.N, dNdu, edgeDir < 0);

            // Compute basis function derivatives and the edge normal
            fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
            if (fe.detJxW == 0.0) continue; // skip singular points

            if (edgeDir < 0) normal *= -1.0;

            // Cartesian coordinates of current integration point
            X = Xnod * fe.N;
            X.t = time.t;

            if (integrand.getIntegrandType() & Integrand::NORMAL_DERIVS)
            {
              // Compute the p'th order derivative in the normal direction
              if (fe.p == 1)
                fe.N = dNdu.getColumn(t1);
              else if (fe.p > 1)
                this->extractBasis(fe.u,fe.v,t1,fe.p, fe.N, edgeDir < 0);

              if (hasInterfaceElms)
              {
                // Compute derivative for the neighboring element
                this->extractBasis(fe.u,fe.v,t1,fe.p, dN, edgeDir > 0);
                utl::merge(fe.N,dN,MNPC[iel],MNPC[kel]);
              }

#if SP_DEBUG > 4
              std::cout <<"\niel, xi,eta = "<< fe.iel
                        <<" "<< fe.xi <<" "<< fe.eta
                        <<"\ndN ="<< fe.N <<"dNdX ="<< fe.dNdX;
#endif
            }

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= dS*wg[i];
            ok = integrand.evalInt(*A,fe,time,X,normal);
          }
        }

      // Assembly of global system integral
      if (ok && !glInt.assemble(A->ref(),fe.iel))
        ok = false;

      A->destruct();

      if (!ok) return false;
    }

  return true;
}


bool ASMs2D::integrate (Integrand& integrand, int lIndex,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
  if (!surf) return true; // silently ignore empty patches

  PROFILE2("ASMs2D::integrate(B)");

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = abs(edgeDir);   // Tangent direction normal to the patch edge
  const int t2 = 3-abs(edgeDir); // Tangent direction along the patch edge

  // Compute parameter values of the Gauss points along the whole patch edge
  std::array<Matrix,2> gpar;
  for (int d = 0; d < 2; d++)
    if (-1-d == edgeDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(d == 0 ? surf->startparam_u() : surf->startparam_v());
    }
    else if (1+d == edgeDir)
    {
      gpar[d].resize(1,1);
      gpar[d].fill(d == 0 ? surf->endparam_u() : surf->endparam_v());
    }
    else
      this->getGaussPointParameters(gpar[d],d,nGP,xg);

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivsSf> spline;
  surf->computeBasisGrid(gpar[0],gpar[1],spline);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  // Integrate the extraordinary elements?
  size_t doXelms = 0;
  if (integrand.getIntegrandType() & Integrand::XO_ELEMENTS)
    if ((doXelms = (n1-p1+1)*(n2-p2+1))*2 > MNPC.size())
    {
      std::cerr <<" *** ASMs2D::integrate: Too few XO-elements "
                << MNPC.size() - doXelms << std::endl;
      return false;
    }

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  FiniteElement fe(p1*p2);
  fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);

  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   normal;
  double dXidu[2];


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i2 = p2; i2 <= n2; i2++)
    for (int i1 = p1; i1 <= n1; i1++, iel++)
    {
      fe.iel = abs(MLGE[doXelms+iel-1]);
      if (fe.iel < 1) continue; // zero-area element

#ifdef SP_DEBUG
      if (dbgElm < 0 && iel != -dbgElm)
        continue; // Skipping all elements, except for -dbgElm
#endif

      // Skip elements that are not on current boundary edge
      bool skipMe = false;
      switch (edgeDir)
	{
	case -1: if (i1 > p1) skipMe = true; break;
	case  1: if (i1 < n1) skipMe = true; break;
	case -2: if (i2 > p2) skipMe = true; break;
	case  2: if (i2 < n2) skipMe = true; break;
	}
      if (skipMe) continue;

      // Get element edge length in the parameter space
      double dS = 0.5*this->getParametricLength(iel,t2);
      if (dS < 0.0) return false; // topology error (probably logic error)

      // Set up control point coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
        this->getElementCorners(i1-1,i2-1,fe.XC);

      if (integrand.getIntegrandType() & Integrand::G_MATRIX)
      {
        // Element size in parametric space
        dXidu[0] = surf->knotSpan(0,i1-1);
        dXidu[1] = surf->knotSpan(1,i2-1);
      }

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
      bool ok = integrand.initElementBou(MNPC[doXelms+iel-1],*A);


      // --- Integration loop over all Gauss points along the edge -------------

      int ip = (t1 == 1 ? i2-p2 : i1-p1)*nGP;
      fe.iGP = firstp + ip; // Global integration point counter

      for (int i = 0; i < nGP && ok; i++, ip++, fe.iGP++)
      {
	// Local element coordinates and parameter values
	// of current integration point
	if (gpar[0].size() > 1)
	{
	  fe.xi = xg[i];
	  fe.u = gpar[0](i+1,i1-p1+1);
	}
	if (gpar[1].size() > 1)
	{
	  fe.eta = xg[i];
	  fe.v = gpar[1](i+1,i2-p2+1);
	}

	// Fetch basis function derivatives at current integration point
	SplineUtils::extractBasis(spline[ip],fe.N,dNdu);

	// Compute basis function derivatives and the edge normal
	fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
	if (fe.detJxW == 0.0) continue; // skip singular points

	if (edgeDir < 0) normal *= -1.0;

	// Compute G-matrix
	if (integrand.getIntegrandType() & Integrand::G_MATRIX)
	  utl::getGmat(Jac,dXidu,fe.G);

	// Cartesian coordinates of current integration point
	X = Xnod * fe.N;
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= dS*wg[i];
	ok = integrand.evalBou(*A,fe,time,X,normal);
      }

      // Finalize the element quantities
      if (ok && !integrand.finalizeElementBou(*A,fe,time))
        ok = false;

      // Assembly of global system integral
      if (ok && !glInt.assemble(A->ref(),fe.iel))
	ok = false;

      A->destruct();

      if (!ok) return false;

#ifdef SP_DEBUG
      if (dbgElm < 0 && iel == -dbgElm)
        break; // Skipping all elements, except for -dbgElm
#endif
    }

  return true;
}


int ASMs2D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!surf) return -2;

  param[0] = (1.0-xi[0])*surf->startparam_u() + xi[0]*surf->endparam_u();
  param[1] = (1.0-xi[1])*surf->startparam_v() + xi[1]*surf->endparam_v();
  SplineUtils::point(X,param[0],param[1],surf);

  // Check if this point matches any of the control points (nodes)
  return this->searchCtrlPt(surf->coefs_begin(),surf->coefs_end(),
                            X,surf->dimension());
}


bool ASMs2D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
  if (!surf) return false;

  if (nSegPerSpan < 1)
  {
    std::cerr <<" *** ASMs2D::getGridParameters: Too few knot-span points "
	      << nSegPerSpan+1 <<" in direction "<< dir << std::endl;
    return false;
  }

  RealArray::const_iterator uit = surf->basis(dir).begin() + surf->basis(dir).order()-1;
  RealArray::const_iterator uend = surf->basis(dir).begin() + surf->basis(dir).numCoefs()+1;

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


bool ASMs2D::tesselate (ElementBlock& grid, const int* npe) const
{
  // Compute parameter values of the nodal points
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the spline surface at all points
  size_t nx = gpar[0].size();
  size_t ny = gpar[1].size();
  RealArray XYZ(surf->dimension()*nx*ny);
  surf->gridEvaluator(XYZ,gpar[0],gpar[1]);

  // Establish the block grid coordinates
  size_t i, j, l;
  grid.resize(nx,ny);
  for (i = l = 0; i < grid.getNoNodes(); i++, l += surf->dimension())
    for (j = 0; j < nsd; j++)
      grid.setCoor(i,j,XYZ[l+j]);

  // Establish the block grid topology
  int ie, nse1 = npe[0] - 1;
  int je, nse2 = npe[1] - 1;
  int nel1 = (nx-1)/nse1;
  int n[4], ip = 0;
  for (j = je = 1, n[1] = 0; j < ny; j++)
  {
    n[0] = n[1];
    n[1] = n[0] + 1;
    n[2] = n[1] + nx;
    n[3] = n[1] + nx-1;
    for (i = ie = 1; i < nx; i++)
    {
      for (l = 0; l < 4; l++)
	grid.setNode(ip++,n[l]++);
      grid.setElmId((j-1)*(nx-1)+i,(je-1)*nel1+ie);
      if (i%nse1 == 0) ie++;
    }
    if (j%nse2 == 0) je++;
  }

  return true;
}


void ASMs2D::scatterInd (int n1, int n2, int p1, int p2,
			 const int* start, IntVec& index)
{
  index.reserve(p1*p2);
  int ip = ((start[1]-p2+1))*n1 + (start[0]-p1+1);
  for (int i2 = 0; i2 < p2; i2++, ip += n1-p1)
    for (int i1 = 0; i1 < p1; i1++, ip++)
      index.push_back(ip);
}


bool ASMs2D::evalSolution (Matrix& sField, const Vector& locSol,
			   const int* npe) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  return this->evalSolution(sField,locSol,gpar.data());
}


bool ASMs2D::evalSolution (Matrix& sField, const Vector& locSol,
                           const RealArray* gpar, bool regular, int deriv) const
{
  // Evaluate the basis functions at all points
  size_t nPoints = gpar[0].size();
  std::vector<Go::BasisPtsSf>     spline0(regular || deriv != 0 ? 0 : nPoints);
  std::vector<Go::BasisDerivsSf>  spline1(regular || deriv != 1 ? 0 : nPoints);
  std::vector<Go::BasisDerivsSf2> spline2(regular || deriv != 2 ? 0 : nPoints);
  if (regular)
  {
    nPoints *= gpar[1].size();
    switch (deriv) {
    case 0:
      surf->computeBasisGrid(gpar[0],gpar[1],spline0);
      break;
    case 1:
      surf->computeBasisGrid(gpar[0],gpar[1],spline1);
      break;
    case 2:
      surf->computeBasisGrid(gpar[0],gpar[1],spline2);
      break;
    default:
      return false;
    }
  }
  else if (nPoints == gpar[1].size())
  {
    for (size_t i = 0; i < nPoints; i++)
      switch (deriv) {
      case 0:
        surf->computeBasis(gpar[0][i],gpar[1][i],spline0[i]);
        break;
      case 1:
        surf->computeBasis(gpar[0][i],gpar[1][i],spline1[i]);
        break;
      case 2:
        surf->computeBasis(gpar[0][i],gpar[1][i],spline2[i]);
        break;
      default:
        return false;
      }
  }
  else
    return false;

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  size_t nComp = locSol.size() / (n1*n2);

  Vector   ptSol;
  Matrix   dNdu, dNdX, Xnod, Xtmp, Jac, eSol, ptDer;
  Matrix3D d2Ndu2, d2NdX2, Hess, ptDer2;

  // Fetch nodal (control point) coordinates
  this->getNodalCoordinates(Xnod);

  // Evaluate the primary solution field at each point
  sField.resize(nComp*int(pow(nsd,deriv)),nPoints);
  sField.resize(nComp,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    IntVec ip;
    switch (deriv) {

    case 0: // Evaluate the solution
      scatterInd(n1,n2,p1,p2,spline0[i].left_idx,ip);
      utl::gather(ip,nComp,locSol,Xtmp);
      Xtmp.multiply(spline0[i].basisValues,ptSol);
      sField.fillColumn(1+i,ptSol);
      break;

    case 1: // Evaluate first derivatives of the solution
      scatterInd(n1,n2,p1,p2,spline1[i].left_idx,ip);
      SplineUtils::extractBasis(spline1[i],ptSol,dNdu);
      utl::gather(ip,nsd,Xnod,Xtmp);
      utl::Jacobian(Jac,dNdX,Xtmp,dNdu);
      utl::gather(ip,nComp,locSol,Xtmp);
      ptDer.multiply(Xtmp,dNdX);
      sField.fillColumn(1+i,ptDer);
      break;

    case 2: // Evaluate second derivatives of the solution
      scatterInd(n1,n2,p1,p2,spline2[i].left_idx,ip);
      SplineUtils::extractBasis(spline2[i],ptSol,dNdu,d2Ndu2);
      utl::gather(ip,nsd,Xnod,Xtmp);
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


bool ASMs2D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			   const int* npe, char project) const
{
  // Project the secondary solution onto the spline basis
  Go::SplineSurface* s = nullptr;
  if (project == 'S')
    s = this->scRecovery(integrand);
  else if (project == 'A')
    s = this->projectSolutionLocalApprox(integrand);
  else if (project == 'L')
    s = this->projectSolutionLocal(integrand);
  else if (project == 'W')
    s = this->projectSolutionLeastSquare(integrand);
  else if (project == 'D' || !npe)
    s = this->projectSolution(integrand);

  if (npe)
  {
    // Compute parameter values of the result sampling points
    std::array<RealArray,2> gpar;
    if (this->getGridParameters(gpar[0],0,npe[0]-1) &&
        this->getGridParameters(gpar[1],1,npe[1]-1))
    {
      if (!project)
        // Evaluate the secondary solution directly at all sampling points
        return this->evalSolution(sField,integrand,gpar.data());
      else if (s)
      {
        // Evaluate the projected field at the result sampling points
        const Vector& svec = sField; // using utl::matrix cast operator
        sField.resize(s->dimension(),gpar[0].size()*gpar[1].size());
        s->gridEvaluator(const_cast<Vector&>(svec),gpar[0],gpar[1]);
        delete s;
        return true;
      }
    }
    else if (s)
      delete s;
  }
  else if (s)
  {
    // Extract control point values from the spline object
    sField.resize(s->dimension(),s->numCoefs_u()*s->numCoefs_v());
    sField.fill(&(*s->coefs_begin()));
    delete s;
    return true;
  }

  std::cerr <<" *** ASMs2D::evalSolution: Failure!";
  if (project) std::cerr <<" project="<< project;
  std::cerr << std::endl;
  return false;
}


bool ASMs2D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			   const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  // Evaluate the basis functions and their derivatives at all points
  size_t nPoints = gpar[0].size();
  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  std::vector<Go::BasisDerivsSf>  spline1(regular ||  use2ndDer ? 0 : nPoints);
  std::vector<Go::BasisDerivsSf2> spline2(regular || !use2ndDer ? 0 : nPoints);
  if (regular)
  {
    nPoints *= gpar[1].size();
    if (use2ndDer)
      surf->computeBasisGrid(gpar[0],gpar[1],spline2);
    else
      surf->computeBasisGrid(gpar[0],gpar[1],spline1);
  }
  else if (gpar[0].size() == gpar[1].size())
  {
    for (size_t i = 0; i < nPoints; i++)
      if (use2ndDer)
        surf->computeBasis(gpar[0][i],gpar[1][i],spline2[i]);
      else
        surf->computeBasis(gpar[0][i],gpar[1][i],spline1[i]);
  }
  else
    return false;

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  FiniteElement fe(p1*p2,firstIp);
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
      scatterInd(n1,n2,p1,p2,spline2[i].left_idx,ip);
      fe.u = spline2[i].param[0];
      fe.v = spline2[i].param[1];
    }
    else
    {
      scatterInd(n1,n2,p1,p2,spline1[i].left_idx,ip);
      fe.u = spline1[i].param[0];
      fe.v = spline1[i].param[1];
    }

    // Fetch associated control point coordinates
    utl::gather(ip,nsd,Xnod,Xtmp);

    // Fetch basis function derivatives at current integration point
    if (use2ndDer)
      SplineUtils::extractBasis(spline2[i],fe.N,dNdu,d2Ndu2);
    else
      SplineUtils::extractBasis(spline1[i],fe.N,dNdu);

    // Compute the Jacobian inverse and derivatives
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xtmp,dNdu);

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (use2ndDer)
      if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xtmp,d2Ndu2,fe.dNdX))
        continue;

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,Xtmp*fe.N,ip))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMs2D::generateThreadGroups (const Integrand& integrand, bool silence,
                                   bool ignoreGlobalLM)
{
  const int p1 = surf->order_u() - 1;
  const int p2 = surf->order_v() - 1;

  generateThreadGroups(p1, p2, silence, ignoreGlobalLM);
}


void ASMs2D::generateThreadGroups (size_t strip1, size_t strip2,
                                   bool silence, bool ignoreGlobalLM)
{
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  const int p1 = surf->order_u() - 1;
  const int p2 = surf->order_v() - 1;

  std::vector<bool> el1, el2;
  el1.reserve(n1 - p1);
  el2.reserve(n2 - p2);

  int ii;
  for (ii = p1; ii < n1; ii++)
    el1.push_back(surf->knotSpan(0,ii) > 0.0);
  for (ii = p2; ii < n2; ii++)
    el2.push_back(surf->knotSpan(1,ii) > 0.0);

  threadGroups.calcGroups(el1,el2,strip1,strip2);
  if (silence || threadGroups.size() < 2) return;

  std::cout <<"\nMultiple threads are utilized during element assembly.";
  for (size_t i = 0; i < threadGroups.size(); i++)
  {
    std::vector< std::set<int> > nodes(threadGroups[i].size());

    std::cout <<"\n Thread group "<< i+1;
    for (size_t j = 0; j < threadGroups[i].size(); j++)
    {
      std::cout <<"\n\tthread "<< j+1
                << ": "<< threadGroups[i][j].size() <<" elements";
      size_t k, l, nzeroar = 0;
      for (k = 0; k < threadGroups[i][j].size(); k++)
      {
        int iel = threadGroups[i][j][k];
        if (MLGE[iel] > 0)
          for (l = 0; l < MNPC[iel].size(); l++)
            nodes[j].insert(MNPC[iel][l]);
        else
          nzeroar++;
      }
      if (nzeroar > 0)
        std::cout <<" ("<< threadGroups[i][j].size() - nzeroar <<" real)";

      // Verify that the nodes on this thread are not present on the others
      this->checkThreadGroups(nodes, j, ignoreGlobalLM);
    }
  }
  std::cout << std::endl;
}


bool ASMs2D::getNoStructElms (int& n1, int& n2, int& n3) const
{
  n1 = surf->numCoefs_u() - surf->order_u() + 1;
  n2 = surf->numCoefs_v() - surf->order_v() + 1;
  n3 = 0;

  return true;
}


void ASMs2D::extractBasis (double u, double v, Vector& N,
                           Matrix& dNdu, bool fromRight) const
{
  Go::BasisDerivsSf spline;
  surf->computeBasis(u,v,spline,fromRight);
  SplineUtils::extractBasis(spline,N,dNdu);
}


void ASMs2D::extractBasis (double u, double v, Vector& N,
                           Matrix& dNdu, Matrix3D& d2Ndu2, bool fromRight) const
{
  Go::BasisDerivsSf2 spline;
  surf->computeBasis(u,v,spline,fromRight);
  SplineUtils::extractBasis(spline,N,dNdu,d2Ndu2);
}


void ASMs2D::extractBasis (double u, double v, int dir, int p,
                           Vector& dN, bool fromRight) const
{
  Go::BasisDerivsSfU spline;
  surf->computeBasis(u,v,p,spline,fromRight);
  dN.resize(spline.values.size());
  dir += 2*p-2;
  for (size_t i = 0; i < spline.values.size(); i++)
    dN[i] = spline.values[i][dir];
}


short int ASMs2D::InterfaceChecker::hasContribution (int I, int J) const
{
  bool neighbor[4];
  neighbor[0] = I > myPatch.surf->order_u();    // West neighbor
  neighbor[1] = I < myPatch.surf->numCoefs_u(); // East neighbor
  neighbor[2] = J > myPatch.surf->order_v();    // South neighbor
  neighbor[3] = J < myPatch.surf->numCoefs_v(); // North neighbor

  // Check for existing neighbors
  short int status = 0, s = 1;
  for (short int i = 0; i < 4; i++, s *= 2)
    if (neighbor[i]) status += s;

  return status;
}


bool ASMs2D::evaluate (const RealFunc* func, RealArray& vec,
                       int basisNum, double time) const
{
  Go::SplineSurface* oldSurf = this->getBasis(basisNum);
  Go::SplineSurface* newSurf = SplineUtils::project(oldSurf,*func,time);
  if (!newSurf)
  {
    std::cerr <<" *** ASMs2D::evaluate: Projection failure."<< std::endl;
    return false;
  }

  vec.assign(newSurf->coefs_begin(),newSurf->coefs_end());
  delete newSurf;

  return true;
}


int ASMs2D::getCorner (int I, int J, int basis) const
{
  int n1, n2, node = 1;
  for (char i = 1; i <= basis; i++)
    if (!this->getSize(n1,n2,i))
      return -1;
    else if (i < basis)
      node += n1*n2;

  if (swapV) // Account for swapped parameter direction
    J = -J;

  if (I > 0) node += n1-1;
  if (J > 0) node += n1*(n2-1);

  return node;
}
