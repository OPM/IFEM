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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SurfaceInterpolator.h"

#include "ASMs2D.h"
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


ASMs2D::ASMs2D (const char* fileName, unsigned char n_s, unsigned char n_f)
  : ASMstruct(2,n_s,n_f), surf(0)
{
  std::cout <<"\nReading patch file "<< fileName << std::endl;
  std::ifstream is(fileName);
  if (!is.good())
    std::cerr <<" *** ASMs2D: Failure opening patch file"<< std::endl;
  else
    this->read(is);

  geo = surf;
}


ASMs2D::ASMs2D (std::istream& is, unsigned char n_s, unsigned char n_f)
  : ASMstruct(2,n_s,n_f), surf(0)
{
  this->read(is);

  geo = surf;
}


ASMs2D::ASMs2D (unsigned char n_s, unsigned char n_f)
  : ASMstruct(2,n_s,n_f), surf(0)
{
}


bool ASMs2D::read (std::istream& is)
{
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
    return false;
  }
  else if (surf->dimension() < nsd)
  {
    std::cout <<"  ** ASMs2D::read: The dimension of this surface patch "
	      << surf->dimension() <<" is less than nsd="<< nsd
	      <<".\n                   Resetting nsd to "<< surf->dimension()
	      <<" for this patch."<< std::endl;
    nsd = surf->dimension();
  }

  return true;
}


bool ASMs2D::write (std::ostream& os) const
{
  if (!surf) return false;

  os <<"200 1 0 0\n";
  os << *surf;

  return os.good();
}


void ASMs2D::clear ()
{
  // Erase spline data
  if (surf) delete surf;
  surf = 0;
  geo = 0;

  // Erase the FE data
  nodeInd.clear();
  ASMbase::clear();
}


bool ASMs2D::refine (int dir, const RealArray& xi)
{
  if (!surf || dir < 0 || dir > 1 || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;

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

  surf->raiseOrder(ru,rv);
  return true;
}


bool ASMs2D::generateFEMTopology ()
{
  if (!surf) return false;

  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  if (!nodeInd.empty())
  {
    if (nodeInd.size() == (size_t)n1*n2) return true;
    std::cerr <<" *** ASMs2D::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< nodeInd.size()
	      <<"\n     and the number of spline coefficients "<< n1*n2
	      <<" in the patch."<< std::endl;
    return false;
  }

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
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

  MLGE.resize((n1-p1+1)*(n2-p2+1),0);
  MLGN.resize(n1*n2);
  MNPC.resize(MLGE.size());
  nodeInd.resize(MLGN.size());

  int iel = 0;
  int inod = 0;
  for (int i2 = 1; i2 <= n2; i2++)
    for (int i1 = 1; i1 <= n1; i1++)
    {
      nodeInd[inod].I = i1-1;
      nodeInd[inod].J = i2-1;
      if (i1 >= p1 && i2 >= p2)
      {
	if (surf->knotSpan(0,i1-1) > 0.0)
	  if (surf->knotSpan(1,i2-1) > 0.0)
          {
	    MLGE[iel] = ++gEl; // global element number over all patches
	    MNPC[iel].resize(p1*p2,0);

	    int lnod = 0;
	    for (int j2 = p2-1; j2 >= 0; j2--)
	      for (int j1 = p1-1; j1 >= 0; j1--)
		MNPC[iel][lnod++] = inod - n1*j2 - j1;
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
	  MLGN[inod] = nodes.edges[2].next();
	else if (i == n1)
	  MLGN[inod] = nodes.edges[3].next();
	else
	  MLGN[inod] = nodes.next();
      }

#if SP_DEBUG > 1
  if (basis > 0) std::cout <<"\nBasis "<< basis <<":";
  for (int i = inod-n1*n2; i < inod; i++)
    std::cout <<"\nNode "<< i+1 <<"\t: "<< nodeInd[i].I <<" "<< nodeInd[i].J
	      <<"\tglobal no. "<< MLGN[i];
  std::cout << std::endl;
#endif
  return true;
}


bool ASMs2D::connectPatch (int edge, ASMs2D& neighbor, int nedge, bool rev)
{
  // Set up the slave node numbers for this surface patch

  int n1, n2;
  if (!this->getSize(n1,n2)) return false;
  int node = 1, i1 = 0;

  switch (edge)
    {
    case 2: // Positive I-direction
      node = n1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      break;

    case 4: // Positive J-direction
      node = n1*(n2-1)+1;
    case 3: // Negative J-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** ASMs2D::connectPatch: Invalid slave edge "
		<< edge << std::endl;
      return false;
    }

  int i;
  IntVec slaveNodes(n1,0);
  for (i = 0; i < n1; i++, node += i1)
    slaveNodes[i] = node;

  // Set up the master node numbers for the neighboring surface patch

  if (!neighbor.getSize(n1,n2)) return false;
  node = 1; i1 = 0;

  switch (nedge)
    {
    case 2: // Positive I-direction
      node = n1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      break;

    case 4: // Positive J-direction
      node = n1*(n2-1)+1;
    case 3: // Negative J-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** ASMs2D::connectPatch: Invalid master edge "
		<< nedge << std::endl;
      return false;
    }

  if (n1 != (int)slaveNodes.size())
  {
    std::cerr <<" *** ASMs2D::connectPatch: Non-matching edges, sizes "
	      << n1 <<" and "<< slaveNodes.size() << std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  for (i = 0; i < n1; i++, node += i1)
  {
    int k = rev ? n1-i-1 : i;
    if (!neighbor.getCoord(node).equal(this->getCoord(slaveNodes[k]),xtol))
    {
      std::cerr <<" *** ASMs2D::connectPatch: Non-matching nodes "
		<< node <<": "<< neighbor.getCoord(node)
		<<"\n                                          and "
		<< slaveNodes[k] <<": "<< this->getCoord(slaveNodes[k])
		<< std::endl;
      return false;
    }
    else
      ASMbase::collapseNodes(neighbor.MLGN[node-1],MLGN[slaveNodes[k]-1]);
  }

  return true;
}


void ASMs2D::closeEdges (int dir)
{
  int n1, n2, master = 1;
  if (!this->getSize(n1,n2)) return;

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


void ASMs2D::constrainEdge (int dir, int dof, int code)
{
  int n1, n2, node = 1;
  if (!this->getSize(n1,n2,1)) return;

  switch (dir)
    {
    case  1: // Right edge (positive I-direction)
      node += n1-1;
    case -1: // Left edge (negative I-direction)
      for (int i2 = 1; i2 <= n2; i2++, node += n1)
	this->prescribe(node,dof,code);
      break;

    case  2: // Back edge (positive J-direction)
      node += n1*(n2-1);
    case -2: // Front edge (negative J-direction)
      for (int i1 = 1; i1 <= n1; i1++, node++)
	this->prescribe(node,dof,code);
      break;
    }
}


void ASMs2D::constrainCorner (int I, int J, int dof, int code)
{
  int n1, n2;
  if (!this->getSize(n1,n2,1)) return;

  int node = 1;
  if (I > 0) node += n1-1;
  if (J > 0) node += n1*(n2-1);

  this->prescribe(node,dof,code);
}


void ASMs2D::constrainNode (double xi, double eta, int dof, int code)
{
  if (xi  < 0.0 || xi  > 1.0) return;
  if (eta < 0.0 || eta > 1.0) return;

  int n1, n2;
  if (!this->getSize(n1,n2,1)) return;

  int node = 1;
  if (xi  > 0.0) node += int(0.5+(n1-1)*xi);
  if (eta > 0.0) node += n1*int(0.5+(n2-1)*eta);

  this->prescribe(node,dof,code);
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

  int inod1 = MNPC[iel-1].back();
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nodeInd.size())
  {
    std::cerr <<" *** ASMs2D::getParametricArea: Node index "<< inod1
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
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

  int inod1 = MNPC[iel-1].back();
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= nodeInd.size())
  {
    std::cerr <<" *** ASMs2D::getParametricLength: Node index "<< inod1
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
    return DERR;
  }
#endif

  const int ni = nodeInd[inod1].I;
  const int nj = nodeInd[inod1].J;
  switch (dir)
    {
    case 1: return surf->knotSpan(1,nj);
    case 2: return surf->knotSpan(0,ni);
    }

  std::cerr <<" *** ASMs2D::getParametricLength: Invalid edge direction "
	    << dir << std::endl;
  return DERR;
}


int ASMs2D::coeffInd (size_t inod) const
{
#ifdef INDEX_CHECK
  if (inod >= nodeInd.size())
  {
    std::cerr <<" *** ASMs2D::coeffInd: Node index "<< inod
	      <<" out of range [0,"<< nodeInd.size() <<">."<< std::endl;
    return -1;
  }
#endif

  const int ni = nodeInd[inod].I;
  const int nj = nodeInd[inod].J;
  return nj*surf->numCoefs_u() + ni;
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

  const IntVec& mnpc = MNPC[iel-1];
  X.resize(nsd,mnpc.size());

  RealArray::const_iterator cit = surf->coefs_begin();
  for (size_t n = 0; n < mnpc.size(); n++)
  {
    int ip = this->coeffInd(mnpc[n])*surf->dimension();
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


Vec3 ASMs2D::getCoord (size_t inod) const
{
  Vec3 X;
  int ip = this->coeffInd(inod-1)*surf->dimension();
  if (ip < 0) return X;

  RealArray::const_iterator cit = surf->coefs_begin() + ip;
  for (unsigned char i = 0; i < nsd; i++, cit++)
    X[i] = *cit;

  return X;
}


bool ASMs2D::getSize (int& n1, int& n2, int) const
{
  if (!surf) return false;

  n1 = surf->numCoefs_u();
  n2 = surf->numCoefs_v();
  return true;
}


/*!
  \brief Computes the characteristic element length from nodal coordinates.
*/

static double getElmSize (int p1, int p2, const Matrix& X)
{
  int n = X.rows();
  int i, j, id1, id2;
  double value, v1, h = 1.0e12;

  // Y-direction
  for (i = 1; i <= p1; i++)
  {
    id1 = i;
    id2 = id1 + (p2-1)*p1;
    value = 0.0;
    for (j = 1; j <= n; j++)
    {
      v1 = X(j,id2) - X(j,id1);
      value += v1*v1;
    }
    if (value < h) h = value;
  }

  // X-direction
  for (j = 0; j < p2; j++)
  {
    id1 = j*p1 + 1;
    id2 = id1 + p1 - 1;
    value = 0.0;
    for (i = 1; i <= n; i++)
    {
      v1 = X(i,id2) - X(i,id1);
      value += v1*v1;
    }
    if (value < h) h = value;
  }

  return sqrt(h);
}


void ASMs2D::extractBasis (const Go::BasisDerivsSf& spline,
			   Vector& N, Matrix& dNdu)
{
  dNdu.resize(N.size(),2);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
     N  (n)   = spline.basisValues[jp];
    dNdu(n,1) = spline.basisDerivs_u[jp];
    dNdu(n,2) = spline.basisDerivs_v[jp];
  }
}


void ASMs2D::extractBasis (const Go::BasisDerivsSf2& spline,
			   Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2)
{
   dNdu .resize(N.size(),2);
  d2Ndu2.resize(N.size(),2,2);

  size_t jp, n = 1;
  for (jp = 0; jp < N.size(); jp++, n++)
  {
      N   (n)     = spline.basisValues[jp];
     dNdu (n,1)   = spline.basisDerivs_u[jp];
     dNdu (n,2)   = spline.basisDerivs_v[jp];
    d2Ndu2(n,1,1) = spline.basisDerivs_uu[jp];
    d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = spline.basisDerivs_uv[jp];
    d2Ndu2(n,2,2) = spline.basisDerivs_vv[jp];
  }
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
			const TimeDomain& time,
			const LintegralVec& locInt)
{
  if (!surf) return true; // silently ignore empty patches

  PROFILE2("ASMs2D::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  int dir;
  Matrix gpar[2];
  for (dir = 0; dir < 2; dir++)
  {
    int pm1 = (dir == 0 ? surf->order_u() : surf->order_v()) - 1;
    RealArray::const_iterator uit = surf->basis(dir).begin() + pm1;
    double ucurr, uprev = *(uit++);
    int nCol = (dir == 0 ? surf->numCoefs_u() : surf->numCoefs_v()) - pm1;
    gpar[dir].resize(nGauss,nCol);
    for (int j = 1; j <= nCol; uit++, j++)
    {
      ucurr = *uit;
      for (int i = 1; i <= nGauss; i++)
	gpar[dir](i,j) = 0.5*((ucurr-uprev)*xg[i-1] + ucurr+uprev);
      uprev = ucurr;
    }
  }

  // Evaluate basis function derivatives at all integration points
  std::vector<Go::BasisDerivsSf>  spline;
  std::vector<Go::BasisDerivsSf2> spline2;
  if (integrand.getIntegrandType() == 2)
    surf->computeBasisGrid(gpar[0],gpar[1],spline2);
  else
    surf->computeBasisGrid(gpar[0],gpar[1],spline);

#if SP_DEBUG > 4
  for (size_t i = 0; i < spline.size(); i++)
    std::cout <<"\nBasis functions at integration point "<< 1+i << spline[i];
#endif

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  const int nel1 = n1 - p1 + 1;

  FiniteElement fe(p1*p2);
  Matrix   dNdu, Xnod, Jac;
  Matrix3D d2Ndu2, Hess;
  Vec4     X;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 1;
  for (int i2 = p2; i2 <= n2; i2++)
    for (int i1 = p1; i1 <= n1; i1++, iel++)
    {
      if (MLGE[iel-1] < 1) continue; // zero-area element

      // Get element area in the parameter space
      double dA = this->getParametricArea(iel);
      if (dA < 0.0) return false; // topology error (probably logic error)

      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      // Compute characteristic element length, if needed
      if (integrand.getIntegrandType() == 2)
	fe.h = getElmSize(p1,p2,Xnod);

      else if (integrand.getIntegrandType() == 3)
      {
	// --- Compute average value of basis functions over the element -------

	fe.Navg.resize(p1*p2,true);
	double area = 0.0;
	int ip = ((i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
	for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
	  for (int i = 0; i < nGauss; i++, ip++)
	  {
	    // Fetch basis function derivatives at current integration point
	    extractBasis(spline[ip],fe.N,dNdu);

	    // Compute Jacobian determinant of coordinate mapping
	    // and multiply by weight of current integration point
	    double detJac = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu,false);
	    double weight = 0.25*dA*wg[i]*wg[j];

	    // Numerical quadrature
	    fe.Navg.add(fe.N,detJac*weight);
	    area += detJac*weight;
	  }

	// Divide by element area
	fe.Navg /= area;
      }

      else if (integrand.getIntegrandType() == 4)
      {
	// Compute the element center
	Go::Point X0;
	double u0 = 0.5*(gpar[0](1,i1-p1+1) + gpar[0](nGauss,i1-p1+1));
	double v0 = 0.5*(gpar[1](1,i2-p2+1) + gpar[1](nGauss,i2-p2+1));
	surf->point(X0,u0,v0);
	for (unsigned char i = 0; i < nsd; i++)
	  X[i] = X0[i];
      }

      // Initialize element quantities
      if (!integrand.initElement(MNPC[iel-1],X,nGauss*nGauss)) return false;

      // Caution: Unless locInt is empty, we assume it points to an array of
      // LocalIntegral pointers, of length at least the number of elements in
      // the model (as defined by the highest number in the MLGE array).
      // If the array is shorter than this, expect a segmentation fault.
      LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


      // --- Integration loop over all Gauss points in each direction ----------

      int ip = ((i2-p2)*nGauss*nel1 + i1-p1)*nGauss;
      for (int j = 0; j < nGauss; j++, ip += nGauss*(nel1-1))
	for (int i = 0; i < nGauss; i++, ip++)
	{
	  // Parameter values of current integration point
	  fe.u = gpar[0](i+1,i1-p1+1);
	  fe.v = gpar[1](j+1,i2-p2+1);

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
	  fe.detJxW *= 0.25*dA*wg[i]*wg[j];
	  if (!integrand.evalInt(elmInt,fe,time,X))
	    return false;
	}

      // Finalize the element quantities
      if (!integrand.finalizeElement(elmInt,time))
	return false;

      // Assembly of global system integral
      if (!glInt.assemble(elmInt,MLGE[iel-1]))
	return false;
    }

  return true;
}


bool ASMs2D::integrate (Integrand& integrand, int lIndex,
			GlobalIntegral& glInt,
			const TimeDomain& time,
			const LintegralVec& locInt)
{
  if (!surf) return true; // silently ignore empty patches

  PROFILE2("ASMs2D::integrate(B)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = abs(edgeDir);   // Tangent direction normal to the patch edge
  const int t2 = 3-abs(edgeDir); // Tangent direction along the patch edge

  // Compute parameter values of the Gauss points along the whole patch edge
  Matrix gpar[2];
  for (short int d = 0; d < 2; d++)
    if (-1-d == edgeDir)
    {
      gpar[d].resize(1,1);
      gpar[d](1,1) = d == 0 ? surf->startparam_u() : surf->startparam_v();
    }
    else if (1+d == edgeDir)
    {
      gpar[d].resize(1,1);
      gpar[d](1,1) = d == 0 ? surf->endparam_u() : surf->endparam_v();
    }
    else
    {
      int pm1 = (d == 0 ? surf->order_u() : surf->order_v()) - 1;
      RealArray::const_iterator uit = surf->basis(d).begin() + pm1;
      double ucurr, uprev = *(uit++);
      int nCol = (d == 0 ? surf->numCoefs_u() : surf->numCoefs_v()) - pm1;
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
  std::vector<Go::BasisDerivsSf> spline;
  surf->computeBasisGrid(gpar[0],gpar[1],spline);

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();

  FiniteElement fe(p1*p2);
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);

  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   normal;


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i2 = p2; i2 <= n2; i2++)
    for (int i1 = p1; i1 <= n1; i1++, iel++)
    {
      if (MLGE[iel-1] < 1) continue; // zero-area element

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
      double dS = this->getParametricLength(iel,t1);
      if (dS < 0.0) return false; // topology error (probably logic error)

      // Set up control point coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel)) return false;

      // Initialize element quantities
      if (!integrand.initElementBou(MNPC[iel-1])) return false;

      // Caution: Unless locInt is empty, we assume it points to an array of
      // LocalIntegral pointers, of length at least the number of elements in
      // the model (as defined by the highest number in the MLGE array).
      // If the array is shorter than this, expect a segmentation fault.
      LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


      // --- Integration loop over all Gauss points along the edge -------------

      int ip = (t1 == 1 ? i2-p2 : i1-p1)*nGauss;
      for (int i = 0; i < nGauss; i++, ip++)
      {
	// Parameter values of current integration point
	if (gpar[0].size() > 1) fe.u = gpar[0](i+1,i1-p1+1);
	if (gpar[1].size() > 1) fe.v = gpar[1](i+1,i2-p2+1);

	// Fetch basis function derivatives at current integration point
	extractBasis(spline[ip],fe.N,dNdu);

	// Compute basis function derivatives and the edge normal
	fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
	if (fe.detJxW == 0.0) continue; // skip singular points

	if (edgeDir < 0) normal *= -1.0;

	// Cartesian coordinates of current integration point
	X = Xnod * fe.N;
	X.t = time.t;

	// Evaluate the integrand and accumulate element contributions
	fe.detJxW *= 0.5*dS*wg[i];
	if (!integrand.evalBou(elmInt,fe,time,X,normal))
	  return false;
      }

      // Assembly of global system integral
      if (!glInt.assemble(elmInt,MLGE[iel-1]))
	return false;
    }

  return true;
}


int ASMs2D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!surf) return -2;

  param[0] = (1.0-xi[0])*surf->startparam_u() + xi[0]*surf->endparam_u();
  param[1] = (1.0-xi[1])*surf->startparam_v() + xi[1]*surf->endparam_v();

  Go::Point X0;
  surf->point(X0,param[0],param[1]);
  for (unsigned char d = 0; d < nsd; d++)
    X[d] = X0[d];

  // Check if this point matches any of the control points (nodes)
  Vec3 Xnod;
  size_t inod = 1;
  RealArray::const_iterator cit = surf->coefs_begin();
  for (int i = 0; cit != surf->coefs_end(); cit++, i++)
  {
    if (i < nsd) Xnod[i] = *cit;
    if (i+1 == surf->dimension())
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



bool ASMs2D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
  if (!surf) return false;

  if (nSegPerSpan < 1)
  {
    std::cerr <<" *** ASMs2D::getGridParameters: Too few knot-span points "
	      << nSegPerSpan+1 <<" in direction "<< dir << std::endl;
    return false;
  }

  RealArray::const_iterator uit = surf->basis(dir).begin();
  double ucurr = 0.0, uprev = *(uit++);
  while (uit != surf->basis(dir).end())
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


bool ASMs2D::getGrevilleParameters (RealArray& prm, int dir) const
{
  if (!surf) return false;

  const Go::BsplineBasis& basis = surf->basis(dir);

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}


bool ASMs2D::tesselate (ElementBlock& grid, const int* npe) const
{
  // Compute parameter values of the nodal points
  RealArray gpar[2];
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
  double X[3] = { 0.0, 0.0, 0.0 };
  grid.resize(nx,ny);
  for (i = j = 0; i < grid.getNoNodes(); i++, j += surf->dimension())
  {
    for (l = 0; l < nsd; l++)
      X[l] = XYZ[j+l];
    grid.setCoor(i,X[0],X[1],X[2]);
  }

  // Establish the block grid topology
  int n[4], ip = 0;
  for (j = 1, n[1] = 0; j < ny; j++)
  {
    n[0] = n[1];
    n[1] = n[0] + 1;
    n[2] = n[1] + nx;
    n[3] = n[1] + nx-1;
    for (i = 1; i < nx; i++)
      for (l = 0; l < 4; l++)
	grid.setNode(ip++,n[l]++);
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
  RealArray gpar[2];
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  return this->evalSolution(sField,locSol,gpar);
}


bool ASMs2D::evalSolution (Matrix& sField, const Vector& locSol,
			   const RealArray* gpar, bool regular) const
{
  // Evaluate the basis functions at all points
  std::vector<Go::BasisPtsSf> spline;
  if (regular)
    surf->computeBasisGrid(gpar[0],gpar[1],spline);
  else if (gpar[0].size() == gpar[1].size())
  {
    spline.resize(gpar[0].size());
    std::vector<Go::BasisPtsSf> tmpSpline(1);
    for (size_t i = 0; i < spline.size(); i++)
    {
      surf->computeBasisGrid(RealArray(1,gpar[0][i]),
			     RealArray(1,gpar[1][i]),
			     tmpSpline);
      spline[i] = tmpSpline.front();
    }
    // TODO: Request a GoTools method replacing the above:
    // void SplineSurface::computeBasisGrid(double param_u, double param_v,
    //                                      BasisPtsSf& result) const
    /*
    spline.resize(gpar[0].size());
    for (size_t i = 0; i < spline.size(); i++)
      surf->computeBasis(gpar[0][i],gpar[1][i],spline[i]);
    */
  }
  else
    return false;

  const int p1 = surf->order_u();
  const int p2 = surf->order_v();
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
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
    scatterInd(n1,n2,p1,p2,spline[i].left_idx,ip);

    utl::gather(ip,nComp,locSol,Xtmp);
    Xtmp.multiply(spline[i].basisValues,Ytmp);
    sField.fillColumn(1+i,Ytmp);
  }

  return true;
}


bool ASMs2D::evalSolution (Matrix& sField, const Integrand& integrand,
			   const int* npe) const
{
  bool retVal = true;

  if (npe)
  {
    // Compute parameter values of the result sampling points
    RealArray gpar[2];
    for (int dir = 0; dir < 2 && retVal; dir++)
      retVal = this->getGridParameters(gpar[dir],dir,npe[dir]-1);

    // Evaluate the secondary solution at all sampling points
    if (retVal)
      retVal = this->evalSolution(sField,integrand,gpar);
  }
  else
  {
    // Project the secondary solution onto the spline basis
    Go::SplineSurface* s = this->projectSolution(integrand);
    if (s)
    {
      // Extract control point values from the spline object
      sField.resize(s->dimension(),s->numCoefs_u()*s->numCoefs_v());
      sField.fill(&(*s->coefs_begin()));
      delete s;
    }
    else
      retVal = false;
  }

  return retVal;
}


Go::GeomObject* ASMs2D::evalSolution (const Integrand& integrand) const
{
  return this->projectSolution(integrand);
}


Go::SplineSurface* ASMs2D::projectSolution (const Integrand& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar[2];
  for (int dir = 0; dir < 2; dir++)
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
  if (surf->rational())
    surf->getWeights(weights);

  const Vector& vec = sValues;
  return Go::SurfaceInterpolator::regularInterpolation(surf->basis(0),
						       surf->basis(1),
						       gpar[0], gpar[1],
						       const_cast<Vector&>(vec),
						       sValues.rows(),
						       surf->rational(),
						       weights);
}


bool ASMs2D::evalSolution (Matrix& sField, const Integrand& integrand,
			   const RealArray* gpar, bool regular) const
{
  sField.resize(0,0);

  // Evaluate the basis functions and their derivatives at all points
  std::vector<Go::BasisDerivsSf> spline;
  if (regular)
    surf->computeBasisGrid(gpar[0],gpar[1],spline);
  else if (gpar[0].size() == gpar[1].size())
  {
    spline.resize(gpar[0].size());
    std::vector<Go::BasisDerivsSf> tmpSpline(1);
    for (size_t i = 0; i < spline.size(); i++)
    {
      surf->computeBasisGrid(RealArray(1,gpar[0][i]),
			     RealArray(1,gpar[1][i]),
			     tmpSpline);
      spline[i] = tmpSpline.front();
    }
    // TODO: Request a GoTools method replacing the above:
    // void SplineSurface::computeBasisGrid(double param_u, double param_v,
    //                                      BasisDerivsSf& result) const
    /*
    spline.resize(gpar[0].size());
    for (size_t i = 0; i < spline.size(); i++)
      surf->computeBasis(gpar[0][i],gpar[1][i],spline[i]);
    */
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

  Vector N(p1*p2), solPt;
  Matrix dNdu, dNdX, Jac;

  // Evaluate the secondary solution field at each point
  size_t nPoints = spline.size();
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntVec ip;
    scatterInd(n1,n2,p1,p2,spline[i].left_idx,ip);

    // Fetch associated control point coordinates
    utl::gather(ip,nsd,Xnod,Xtmp);

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
