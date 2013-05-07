// $Id$
//==============================================================================
//!
//! \file ASMs1D.C
//!
//! \date Apr 20 2010
//!
//! \author Einar Christensen / SINTEF
//!
//! \brief Driver for assembly of structured 1D spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveInterpolator.h"

#include "ASMs1D.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "SparseMatrix.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Vec3Oper.h"


ASMs1D::ASMs1D (unsigned char n_s, unsigned char n_f)
  : ASMstruct(1,n_s,n_f), curv(0)
{
}


ASMs1D::ASMs1D (const ASMs1D& patch, unsigned char n_f)
  : ASMstruct(patch,n_f), curv(patch.curv)
{
}


bool ASMs1D::read (std::istream& is)
{
  if (shareFE) return true;
  if (curv) delete curv;

  Go::ObjectHeader head;
  curv = new Go::SplineCurve;
  is >> head >> *curv;

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
    std::cerr <<" *** ASMs1D::read: Failure reading spline data"<< std::endl;
    delete curv;
    curv = 0;
    return false;
  }
  else if (curv->dimension() < 1)
  {
    std::cerr <<" *** ASMs1D::read: Invalid spline curve patch, dim="
	      << curv->dimension() << std::endl;
    delete curv;
    curv = 0;
    return false;
  }
  else if (curv->dimension() < nsd)
  {
    std::cout <<"  ** ASMs1D::read: The dimension of this curve patch "
	      << curv->dimension() <<" is less than nsd="<< nsd
	      <<".\n                   Resetting nsd to "<< curv->dimension()
	      <<" for this patch."<< std::endl;
    nsd = curv->dimension();
  }

  geo = curv;
  return true;
}


bool ASMs1D::write (std::ostream& os, int) const
{
  if (!curv) return false;

  os <<"100 1 0 0\n";
  os << *curv;

  return os.good();
}


void ASMs1D::clear (bool retainGeometry)
{
  if (!retainGeometry)
  {
    // Erase spline data
    if (curv && !shareFE) delete curv;
    curv = 0;
    geo = 0;
  }

  // Erase the FE data
  this->ASMbase::clear(retainGeometry);
}


bool ASMs1D::refine (const RealArray& xi)
{
  if (!curv || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = curv->basis().begin();
  double ucurr, uprev = *(uit++);
  while (uit != curv->basis().end())
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

  curv->insertKnot(extraKnots);
  return true;
}


bool ASMs1D::uniformRefine (int nInsert)
{
  if (!curv || nInsert < 1) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = curv->basis().begin();
  double ucurr, uprev = *(uit++);
  while (uit != curv->basis().end())
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

  curv->insertKnot(extraKnots);
  return true;
}


bool ASMs1D::raiseOrder (int ru)
{
  if (!curv) return false;
  if (shareFE) return true;

  curv->raiseOrder(ru);
  return true;
}


bool ASMs1D::generateFEMTopology ()
{
  if (!curv) return false;

  const int n1 = curv->numCoefs();
  if (!MLGN.empty())
  {
    if (MLGN.size() == (size_t)n1) return true;
    std::cerr <<" *** ASMs1D::generateFEMTopology: Inconsistency between the"
	      <<" number of FE nodes "<< MLGN.size()
	      <<"\n     and the number of spline coefficients "<< n1
	      <<" in the patch."<< std::endl;
    return false;
  }
  else if (shareFE)
    return true;

  const int p1 = curv->order();
#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1;
  std::cout <<"\norder: "<< p1;
  std::cout <<"\ndu:";
  for (int i = 0; i < n1; i++)
    std::cout <<' '<< this->getKnotSpan(i);
  std::cout << std::endl;
#endif
  // Consistency checks, just to be fool-proof
  if (n1 <  2) return false;
  if (p1 <  1) return false;
  if (p1 > n1) return false;

  myMLGE.resize(n1-p1+1,0);
  myMLGN.resize(n1);
  myMNPC.resize(myMLGE.size());

  int iel = 0;
  int inod = 0;
  for (int i1 = 1; i1 <= n1; i1++)
  {
    if (i1 >= p1)
    {
      if (this->getKnotSpan(i1-1) > 0.0)
      {
	myMLGE[iel] = ++gEl; // global element number over all patches
	myMNPC[iel].resize(p1,0);

	int lnod = 0;
	for (int j1 = p1-1; j1 >= 0; j1--)
	  myMNPC[iel][lnod++] = inod - j1;
      }

      iel++;
    }
    myMLGN[inod++] = ++gNod; // global node number over all patches
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< iel <<" NNOD = "<< inod << std::endl;
#endif
  return true;
}


bool ASMs1D::connectPatch (int vertex, ASMs1D& neighbor, int nvertex)
{
  if (!this->connectBasis(vertex,neighbor,nvertex))
    return false;

  this->addNeighbor(&neighbor);
  return true;
}


bool ASMs1D::connectBasis (int vertex, ASMs1D& neighbor, int nvertex,
			   int basis, int slave, int master)
{
  if (shareFE && neighbor.shareFE)
    return true;
  else if (shareFE || neighbor.shareFE)
  {
    std::cerr <<" *** ASMs1D::connectPatch: Logic error, cannot"
	      <<" connect a sharedFE patch with an unshared one"<< std::endl;
    return false;
  }

  // Set up the slave node number for this curve patch

  int n1 = this->getSize(basis);

  switch (vertex)
    {
    case 2: // Positive I-direction
      slave += n1;
    case 1: // Negative I-direction
      slave += 1;
      break;

    default:
      std::cerr <<" *** ASMs1D::connectPatch: Invalid slave vertex "
		<< vertex << std::endl;
      return false;
    }

  // Set up the master node number for the neighboring patch

  n1 = neighbor.getSize(basis);

  switch (nvertex)
    {
    case 2: // Positive I-direction
      master += n1;
    case 1: // Negative I-direction
      master += 1;
      break;

    default:
      std::cerr <<" *** ASMs1D::connectPatch: Invalid master vertex "
		<< nvertex << std::endl;
      return false;
    }

  const double xtol = 1.0e-4;
  if (!neighbor.getCoord(master).equal(this->getCoord(slave),xtol))
  {
    std::cerr <<" *** ASMs1D::connectPatch: Non-matching nodes "
	      << master <<": "<< neighbor.getCoord(master)
	      <<"\n                                          and "
	      << slave <<": "<< this->getCoord(slave) << std::endl;
    return false;
  }
  else
    ASMbase::collapseNodes(neighbor,master,*this,slave);

  return true;
}


void ASMs1D::closeEnds (int basis, int master)
{
  if (basis < 1) basis = 1;
  int n1 = this->getSize(basis < 1 ? 1 : basis);
  this->makePeriodic(1,master+n1-1);
}


void ASMs1D::constrainNode (double xi, int dof, int code)
{
  if (xi < 0.0 || xi > 1.0) return;

  int n1 = this->getSize();

  int node = 1;
  if (xi > 0.0) node += int(0.5+(n1-1)*xi);

  this->prescribe(node,dof,code);
}


#define DERR -999.99

double ASMs1D::getParametricLength (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs1D::getParametricLength: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  int inod1 = MNPC[iel-1].back();
#ifdef INDEX_CHECK
  if (inod1 < 0 || (size_t)inod1 >= MLGN.size())
  {
    std::cerr <<" *** ASMs1D::getParametricLength: Node index "<< inod1
	      <<" out of range [0,"<< MLGN.size() <<">."<< std::endl;
    return DERR;
  }
#endif

  return this->getKnotSpan(inod1);
}


double ASMs1D::getKnotSpan (int i) const
{
  if (i < 0 || i >= curv->numCoefs() + curv->order() - 1)
    return 0.0;

  RealArray::const_iterator uit = curv->basis().begin() + i;
  return *(uit+1) - *uit;
}


Vec3 ASMs1D::getCoord (size_t inod) const
{
  Vec3 X;
  int ip = (inod-1)*curv->dimension();
  if (ip < 0) return X;

  RealArray::const_iterator cit = curv->coefs_begin() + ip;
  for (unsigned char i = 0; i < nsd; i++, cit++)
    X[i] = *cit;

  return X;
}


bool ASMs1D::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMs1D::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  const IntVec& mnpc = MNPC[iel-1];
  X.resize(nsd,mnpc.size());

  RealArray::const_iterator cit = curv->coefs_begin();
  for (size_t n = 0; n < mnpc.size(); n++)
  {
    int ip = mnpc[n]*curv->dimension();
    if (ip < 0) return false;

    for (size_t i = 0; i < nsd; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


void ASMs1D::getNodalCoordinates (Matrix& X) const
{
  const int n1 = curv->numCoefs();

  X.resize(nsd,n1);

  RealArray::const_iterator cit = curv->coefs_begin();
  for (int inod = 0; inod < n1; inod++)
  {
    int ip = inod*curv->dimension();
    for (size_t i = 0; i < nsd; i++)
      X(i+1,inod+1) = *(cit+(ip+i));
  }
}


bool ASMs1D::updateCoords (const Vector& displ)
{
  if (!curv) return true; // silently ignore empty patches
  if (shareFE) return true;

  if (displ.size() != nsd*MLGN.size())
  {
    std::cerr <<" *** ASMs1D::updateCoords: Invalid dimension "
	      << displ.size() <<" on displ, should be "
	      << nsd*MLGN.size() << std::endl;
    return false;
  }

  curv->deform(displ,nsd);
  return true;
}


void ASMs1D::getBoundaryNodes (int lIndex, IntVec& glbNodes) const
{
  if (!curv) return; // silently ignore empty patches

  int iel = lIndex == 1 ? 0 : this->getNoElms()-1;
  if (MLGE[iel] > 0)
  {
    if (lIndex == 1)
      glbNodes.push_back(MLGN[MNPC[iel].back()]);
    else if (lIndex == 2)
      glbNodes.push_back(MLGN[MNPC[iel].front()]);
  }
}


bool ASMs1D::getOrder (int& p1, int& p2, int& p3) const
{
  p2 = p3 = 0;
  if (!curv) return false;

  p1 = curv->order();
  return true;
}


bool ASMs1D::getSize (int& n1, int& n2, int& n3, int basis) const
{
  n1 = this->getSize(basis);
  n2 = n3 = 0;
  return true;
}


int ASMs1D::getSize (int) const
{
  if (!curv) return 0;

  return curv->numCoefs();
}


const Vector& ASMs1D::getGaussPointParameters (Matrix& uGP, int nGauss,
					       const double* xi) const
{
  int pm1 = curv->order() - 1;
  RealArray::const_iterator uit = curv->basis().begin() + pm1;

  int nCol = curv->numCoefs() - pm1;
  uGP.resize(nGauss,nCol);

  double ucurr, uprev = *(uit++);
  for (int j = 1; j <= nCol; uit++, j++)
  {
    ucurr = *uit;
    for (int i = 1; i <= nGauss; i++)
      uGP(i,j) = 0.5*((ucurr-uprev)*xi[i-1] + ucurr+uprev);
    uprev = ucurr;
  }

  return uGP;
}


void ASMs1D::extractBasis (double u, Vector& N, Matrix& dNdu) const
{
  int p1 = curv->order();

  RealArray bas(p1*2);
  curv->basis().computeBasisValues(u,&bas.front(),1);

  N.resize(p1);
  dNdu.resize(p1,1);
  for (int i = 1; i <= p1; i++)
  {
     N(i)     = bas[2*i-2];
    dNdu(i,1) = bas[2*i-1];
  }
}


void ASMs1D::extractBasis (double u, Vector& N, Matrix& dNdu,
			   Matrix3D& d2Ndu2) const
{
  int p1 = curv->order();

  RealArray bas(p1*3);
  curv->basis().computeBasisValues(u,&bas.front(),2);

  N.resize(p1);
  dNdu.resize(p1,1);
  d2Ndu2.resize(p1,1,1);
  for (int i = 1; i <= p1; i++)
  {
      N(i)         = bas[3*i-3];
     dNdu(i,1)     = bas[3*i-2];
     d2Ndu2(i,1,1) = bas[3*i-1];
  }
}


bool ASMs1D::integrate (Integrand& integrand,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
  if (!curv) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Get the reduced integration quadrature points, if needed
  const double* xr = 0;
  int nRed = integrand.getReducedIntegration();
  if (nRed < 0)
    nRed = nGauss; // The integrand needs to know nGauss
  else if (nRed > 0 && !(xr = GaussQuadrature::getCoord(nRed)))
    return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gpar, redpar;
  this->getGaussPointParameters(gpar,nGauss,xg);
  if (xr)
    this->getGaussPointParameters(redpar,nRed,xr);

  const int p1 = curv->order();
  const int n1 = curv->numCoefs();

  FiniteElement fe(p1);
  Matrix   dNdu, Xnod, Jac;
  Matrix3D d2Ndu2, Hess;
  Vec4     X;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 1;
  for (int i1 = p1; i1 <= n1; i1++, iel++)
  {
    fe.iel = MLGE[iel-1];
    if (fe.iel < 1) continue; // zero-length element

    // Check that the current element has nonzero length
    double dL = this->getParametricLength(iel);
    if (dL < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel)) return false;

    // Initialize element matrices
    LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
    bool ok = integrand.initElement(MNPC[iel-1],X,nRed,*A);

    if (xr)
    {
      // --- Selective reduced integration loop --------------------------------

      for (int i = 0; i < nRed && ok; i++)
      {
	// Local element coordinates of current integration point
	fe.xi = xr[i];

	// Parameter values of current integration point
	fe.u = redpar(i+1,iel);

	// Fetch basis function derivatives at current point
	this->extractBasis(fe.u,fe.N,dNdu);

	// Compute Jacobian inverse and derivatives
	fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

	// Cartesian coordinates of current integration point
	X = Xnod * fe.N;
	X.t = time.t;

	// Compute the reduced integration terms of the integrand
	ok = integrand.reducedInt(*A,fe,X);
      }
    }


    // --- Integration loop over all Gauss points in current element -----------

    int jp = (iel-1)*nGauss;
    fe.iGP = firstIp + jp; // Global integration point counter

    for (int i = 0; i < nGauss && ok; i++, fe.iGP++)
    {
      // Local element coordinate of current integration point
      fe.xi = xg[i];

      // Parameter value of current integration point
      fe.u = gpar(i+1,iel);

      // Compute basis functions and derivatives
      if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
	this->extractBasis(fe.u,fe.N,dNdu,d2Ndu2);
      else
	this->extractBasis(fe.u,fe.N,dNdu);

      // Compute derivatives in terms of physical co-ordinates
      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
      if (fe.detJxW == 0.0) continue; // skip singular points

      // Compute Hessian of coordinate mapping and 2nd order derivatives
      if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
	if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,fe.dNdX))
	  ok = false;

      // Cartesian coordinates of current integration point
      X = Xnod * fe.N;
      X.t = time.t;

      // Evaluate the integrand and accumulate element contributions
      fe.detJxW *= 0.5*dL*wg[i];
      if (ok && !integrand.evalInt(*A,fe,time,X))
	ok = false;
    }

    // Finalize the element quantities
    if (ok && !integrand.finalizeElement(*A,time,firstIp+jp))
      ok = false;

    // Assembly of global system integral
    if (ok && !glInt.assemble(A->ref(),fe.iel))
      ok = false;

    A->destruct();

    if (!ok) return false;
  }

  return true;
}


bool ASMs1D::integrate (Integrand& integrand, int lIndex,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
  if (!curv) return true; // silently ignore empty patches

  // Integration of boundary point

  FiniteElement fe(curv->order());
  int iel;
  switch (lIndex)
    {
    case 1:
      fe.xi = -1.0;
      fe.u = curv->startparam();
      iel = 1;
      break;

    case 2:
      fe.xi = 1.0;
      fe.u = curv->endparam();
      iel = this->getNoElms();
      break;

    default:
      return false;
    }

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  fe.iGP = iit == firstBp.end() ? 0 : iit->second;
  fe.iel = MLGE[iel-1];
  if (fe.iel < 1) return true; // zero-length element

  // Set up control point coordinates for current element
  Matrix Xnod;
  if (!this->getElementCoordinates(Xnod,iel)) return false;

  // Initialize element matrices
  LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
  bool ok = integrand.initElementBou(MNPC[iel-1],*A);

  // Evaluate basis functions and corresponding derivatives
  Matrix dNdu;
  this->extractBasis(fe.u,fe.N,dNdu);

  // Cartesian coordinates of current integration point
  Vec4 X(Xnod*fe.N,time.t);

  // Compute basis function derivatives
  Matrix Jac, dNdX;
  utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

  // Set up the normal vector
  Vec3 normal;
  if (lIndex == 1)
    normal.x = -copysign(1.0,Jac(1,1));
  else
    normal.x = copysign(1.0,Jac(1,1));

  // Evaluate the integrand and accumulate element contributions
  if (ok && !integrand.evalBou(*A,fe,time,X,normal))
    ok = false;

  // Assembly of global system integral
  if (ok && !glInt.assemble(A->ref(),fe.iel))
    ok = false;

  A->destruct();
  return ok;
}


int ASMs1D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!curv) return -1;

  param[0] = (1.0-xi[0])*curv->startparam() + xi[0]*curv->endparam();

  Go::Point X0;
  curv->point(X0,param[0]);
  for (unsigned char d = 0; d < nsd; d++)
    X[d] = X0[d];

  // Check if this point matches any of the control points (nodes)
  Vec3 Xnod;
  size_t inod = 1;
  RealArray::const_iterator cit = curv->coefs_begin();
  for (int i = 0; cit != curv->coefs_end(); cit++, i++)
  {
    if (i < nsd) Xnod[i] = *cit;
    if (i+1 == curv->dimension())
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


bool ASMs1D::getGridParameters (RealArray& prm, int nSegPerSpan) const
{
  if (!curv) return false;

  if (nSegPerSpan < 1)
  {
    std::cerr <<" *** ASMs1D::getGridParameters: Too few knot-span points "
	      << nSegPerSpan+1 << std::endl;
    return false;
  }

  RealArray::const_iterator uit = curv->basis().begin();
  double ucurr = 0.0, uprev = *(uit++);
  while (uit != curv->basis().end())
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


bool ASMs1D::getGrevilleParameters (RealArray& prm) const
{
  if (!curv) return false;

  const Go::BsplineBasis& basis = curv->basis();

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}


bool ASMs1D::tesselate (ElementBlock& grid, const int* npe) const
{
  // Compute parameter values of the nodal points
  RealArray gpar;
  if (!this->getGridParameters(gpar,npe[0]-1))
    return false;

  // Evaluate the spline curve at all points
  size_t nx = gpar.size();
  RealArray XYZ(curv->dimension()*nx);
  curv->gridEvaluator(XYZ,gpar);

  // Establish the block grid coordinates
  size_t i, j, l;
  grid.resize(nx);
  for (i = l = 0; i < grid.getNoNodes(); i++, l += curv->dimension())
    for (j = 0; j < nsd; j++)
      grid.setCoor(i,j,XYZ[l+j]);

  // Establish the block grid topology
  int nse1 = npe[0] - 1;
  int n[2], ie = 1, ip = 0;
  n[0] = 0;
  n[1] = n[0] + 1;

  for (i = 1; i < nx; i++)
  {
    for (l = 0; l < 2; l++)
      grid.setNode(ip++,n[l]++);
    grid.setElmId(i,ie);
    if (i%nse1 == 0) ie++;
  }

  return true;
}


void ASMs1D::scatterInd (int p1, int start, IntVec& index)
{
  index.resize(p1);
  for (int i1 = 1; i1 <= p1; i1++)
    index[i1-1] = start-p1+i1;
}


bool ASMs1D::evalSolution (Matrix& sField, const Vector& locSol,
			   const int* npe) const
{
  // Compute parameter values of the result sampling points
  RealArray gpar;
  if (!this->getGridParameters(gpar,npe[0]-1))
    return false;

  // Evaluate the primary solution at all sampling points
  return this->evalSolution(sField,locSol,&gpar);
}


bool ASMs1D::evalSolution (Matrix& sField, const Vector& locSol,
			   const RealArray* gpar, bool) const
{
  const int p1 = curv->order();

//  size_t nComp = locSol.size() / this->getNoNodes();
//  if (nComp*this->getNoNodes() != locSol.size())
//    return false;
  size_t nComp = locSol.size() / curv->numCoefs();

  Matrix Xtmp;
  Vector Ytmp, basis(p1);

  // Evaluate the primary solution field at each point
  const RealArray& upar = *gpar;
  size_t nPoints = upar.size();
  sField.resize(nComp,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    curv->basis().computeBasisValues(upar[i],&basis.front());

    IntVec ip;
    scatterInd(p1,curv->basis().lastKnotInterval(),ip);

    utl::gather(ip,nComp,locSol,Xtmp);
    Xtmp.multiply(basis,Ytmp);
    sField.fillColumn(i+1,Ytmp);
  }

  return true;
}


bool ASMs1D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			   const int* npe, char project) const
{
  if (npe)
  {
    // Compute parameter values of the result sampling points
    RealArray gpar;
    if (this->getGridParameters(gpar,npe[0]-1))
      if (project)
      {
	// Project the secondary solution onto the spline basis
	Go::SplineCurve* c = this->projectSolution(integrand);
	if (c)
	{
	  // Evaluate the projected field at the result sampling points
	  const Vector& svec = sField; // using utl::matrix cast operator
	  sField.resize(c->dimension(),gpar.size());
	  c->gridEvaluator(const_cast<Vector&>(svec),gpar);
	  delete c;
	  return true;
	}
      }
      else
	// Evaluate the secondary solution directly at all sampling points
	return this->evalSolution(sField,integrand,&gpar);
  }
  else
  {
    // Project the secondary solution onto the spline basis
    Go::SplineCurve* c = this->projectSolution(integrand);
    if (c)
    {
      // Extract control point values from the spline object
      sField.resize(c->dimension(),c->numCoefs());
      sField.fill(&(*c->coefs_begin()));
      delete c;
      return true;
    }
  }

  std::cerr <<" *** ASMs1D::evalSolution: Failure!"<< std::endl;
  return false;
}


Go::GeomObject* ASMs1D::evalSolution (const IntegrandBase& integrand) const
{
  return this->projectSolution(integrand);
}


Go::SplineCurve* ASMs1D::projectSolution (const IntegrandBase& integrand) const
{
  // Compute parameter values of the result sampling points (Greville points)
  RealArray gpar;
  if (!this->getGrevilleParameters(gpar))
    return 0;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,&gpar))
    return 0;

  // Project the results onto the spline basis to find control point
  // values based on the result values evaluated at the Greville points.
  // Note that we here implicitly assume that the number of Greville points
  // equals the number of control points such that we don't have to resize
  // the result array. Think that is always the case, but beware if trying
  // other projection schemes later.

  RealArray weights;
  if (curv->rational())
    curv->getWeights(weights);

  const Vector& vec = sValues;
  return Go::CurveInterpolator::regularInterpolation(curv->basis(), gpar,
						     const_cast<Vector&>(vec),
						     sValues.rows(),
						     curv->rational(),
						     weights);
}


bool ASMs1D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
			   const RealArray* gpar, bool) const
{
  sField.resize(0,0);

  const int p1 = curv->order();

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  FiniteElement fe(p1);
  Vector        solPt;
  Matrix        dNdu, Jac;
  Matrix3D      d2Ndu2, Hess;

  // Evaluate the secondary solution field at each point
  const RealArray& upar = *gpar;
  size_t nPoints = upar.size();
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch indices of the non-zero basis functions at this point
    IntVec ip;
    scatterInd(p1,curv->basis().lastKnotInterval(),ip);
    fe.u = upar[i];
    fe.iGP = firstIp + i;

    // Fetch associated control point coordinates
    utl::gather(ip,nsd,Xnod,Xtmp);

    // Fetch basis function derivatives at current integration point
    if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
      this->extractBasis(upar[i],fe.N,dNdu,d2Ndu2);
    else
      this->extractBasis(upar[i],fe.N,dNdu);

    // Compute the Jacobian inverse and derivatives
    if (utl::Jacobian(Jac,fe.dNdX,Xtmp,dNdu) == 0.0) // Jac = (Xtmp * dNdu)^-1
      continue; // skip singular points

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
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


bool ASMs1D::globalL2projection (Matrix& sField,
				 const IntegrandBase& integrand) const
{
  if (!curv) return true; // silently ignore empty patches

  const int p1 = curv->order();
  const int n1 = curv->numCoefs();
  const int nel = n1 - p1 + 1;

  // Get Gaussian quadrature points
  const int ng1 = p1 - 1;
  const double* xg = GaussQuadrature::getCoord(ng1);
  if (!xg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gp;
  RealArray gpar = this->getGaussPointParameters(gp,ng1,xg);

  // Evaluate the secondary solution at all integration points
  if (!this->evalSolution(sField,integrand,&gpar))
    return false;

  // Set up the projection matrices
  const size_t nnod = this->getNoNodes();
  const size_t ncomp = sField.rows();
  SparseMatrix A(SparseMatrix::SUPERLU);
  StdVector B(nnod*ncomp);
  A.redim(nnod,nnod);

  RealArray phi(p1);


  // === Assembly loop over all elements in the patch ==========================

  int ip = 0;
  for (int iel = 1; iel <= nel; iel++)
  {
    if (MLGE[iel-1] < 1) continue; // zero-area element

    // --- Integration loop over all Gauss points in current element -----------

    for (int i = 1; i <= ng1; i++, ip++)
    {
      // Fetch basis function values at current integration point
      curv->basis().computeBasisValues(gp(i,iel),&phi.front());

      // Integrate the linear system A*x=B
      for (size_t ii = 0; ii < phi.size(); ii++)
      {
	int inod = MNPC[iel][ii]+1;
	for (size_t jj = 0; jj < phi.size(); jj++)
	{
	  int jnod = MNPC[iel][jj]+1;
	  A(inod,jnod) += phi[ii]*phi[jj];
	}
	for (size_t r = 1; r <= ncomp; r++)
	  B(inod+(r-1)*nnod) += phi[ii]*sField(r,ip+1);
      }
    }
  }

  // Solve the patch-global equation system
  if (!A.solve(B)) return false;

  // Store the control-point values of the projected field
  sField.resize(ncomp,nnod);
  for (size_t i = 1; i <= nnod; i++)
    for (size_t j = 1; j <= ncomp; j++)
      sField(j,i) = B(i+(j-1)*nnod);

  return true;
}
