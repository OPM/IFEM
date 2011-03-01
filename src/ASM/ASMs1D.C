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

#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"

#include "ASMs1D.h"
#include "TimeDomain.h"
#include "GlobalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Vec3Oper.h"
#include <ctype.h>
#include <fstream>

typedef std::vector<double>       DoubleVec;  //!< 1D double array
typedef DoubleVec::const_iterator DoubleIter; //!< Iterator over DoubleVec


ASMs1D::ASMs1D (const char* fileName, unsigned char n_s, unsigned char n_f)
  : ASMstruct(1,n_s,n_f)
{
  curv = 0;

  std::cout <<"\nReading patch file "<< fileName << std::endl;
  std::ifstream is(fileName);
  if (!is.good())
    std::cerr <<" *** ASMs1D: Failure opening patch file"<< std::endl;
  else
    this->read(is);

  geo = curv;
}


ASMs1D::ASMs1D (std::istream& is, unsigned char n_s, unsigned char n_f)
  : ASMstruct(1,n_s,n_f)
{
  curv = 0;

  this->read(is);

  geo = curv;
}


ASMs1D::ASMs1D (unsigned char n_s, unsigned char n_f)
  : ASMstruct(1,n_s,n_f)
{
  curv = 0;
}


bool ASMs1D::read (std::istream& is)
{
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

  return true;
}


bool ASMs1D::write (std::ostream& os) const
{
  if (!curv) return false;

  os <<"100 1 0 0\n";
  os << *curv;

  return os.good();
}


void ASMs1D::clear ()
{
  // Erase spline data
  if (curv) delete curv;
  curv = 0;
  geo = 0;

  // Erase the FE data
  ASMbase::clear();
}



bool ASMs1D::refine (const RealArray& xi)
{
  if (!curv || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;

  DoubleVec  extraKnots;
  DoubleIter uit = curv->basis().begin();
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

  DoubleVec  extraKnots;
  DoubleIter uit = curv->basis().begin();
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

  MLGE.resize(n1-p1+1,0);
  MLGN.resize(n1);
  MNPC.resize(MLGE.size());

  int iel = 0;
  int inod = 0;
  for (int i1 = 1; i1 <= n1; i1++)
  {
    if (i1 >= p1)
    {
      if (this->getKnotSpan(i1-1) > 0.0)
      {
	MLGE[iel] = ++gEl; // global element number over all patches
	MNPC[iel].resize(p1,0);

	int lnod = 0;
	for (int j1 = p1-1; j1 >= 0; j1--)
	  MNPC[iel][lnod++] = inod - j1;
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


bool ASMs1D::connectPatch (int vertex, ASMs1D& neighbor, int nvertex)
{
  // Set up the slave node number for this curve patch

  int n1 = this->getSize();
  int node = 1;

  switch (vertex)
    {
    case 2: // Positive I-direction
      node = n1;
    case 1: // Negative I-direction
      break;

    default:
      std::cerr <<" *** ASMs1D::connectPatch: Invalid slave vertex "
		<< vertex << std::endl;
      return false;
    }

  // Set up the master node number for the neighboring patch

  n1 = neighbor.getSize();
  int nnode = 1;

  switch (nvertex)
    {
    case 2: // Positive I-direction
      nnode = n1;
    case 1: // Negative I-direction
      break;

    default:
      std::cerr <<" *** ASMs1D::connectPatch: Invalid master vertex "
		<< nvertex << std::endl;
      return false;
    }

  const double xtol = 1.0e-4;
  if (!neighbor.getCoord(nnode).equal(this->getCoord(node),xtol))
  {
    std::cerr <<" *** ASMs1D::connectPatch: Non-matching nodes "
	      << nnode <<": "<< neighbor.getCoord(nnode)
	      <<"\n                                          and "
	      << node <<": "<< this->getCoord(node) << std::endl;
    return false;
  }
  else
    ASMbase::collapseNodes(neighbor.MLGN[nnode-1],MLGN[node-1]);

  return true;
}


void ASMs1D::closeEnds ()
{
  int n1 = this->getSize();
  this->makePeriodic(1,n1);
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

  DoubleIter uit = curv->basis().begin() + i;
  return *(uit+1) - *uit;
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

  DoubleIter cit = curv->coefs_begin();
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

  DoubleIter cit = curv->coefs_begin();
  for (int inod = 0; inod < n1; inod++)
  {
    int ip = inod*curv->dimension();
    for (size_t i = 0; i < nsd; i++)
      X(i+1,inod+1) = *(cit+(ip+i));
  }
}


Vec3 ASMs1D::getCoord (size_t inod) const
{
  Vec3 X;
  int ip = (inod-1)*curv->dimension();
  if (ip < 0) return X;

  DoubleIter cit = curv->coefs_begin() + ip;
  for (size_t i = 0; i < nsd; i++, cit++)
    X[i] = *cit;

  return X;
}


int ASMs1D::getSize () const
{
  if (!curv) return 0;

  return curv->numCoefs();
}


void ASMs1D::extractBasis (double u, Vector& N, Matrix& dNdu) const
{
  int p1 = curv->order();

  DoubleVec bas(p1*2);
  curv->basis().computeBasisValues(u,&bas.front(),1);

  N.resize(p1);
  dNdu.resize(p1,1);
  for (int i = 1; i <= p1; i++)
  {
     N(i)     = bas[2*i-2];
    dNdu(i,1) = bas[2*i-1];
  }
}


bool ASMs1D::integrate (Integrand& integrand,
			GlobalIntegral& glInt,
			const TimeDomain& time,
			const LintegralVec& locInt)
{
  if (!curv) return true; // silently ignore empty patches

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch

  const int p1 = curv->order();
  const int n1 = curv->numCoefs();

  Vector gpar;
  int pm1 = p1 - 1;
  DoubleIter uit = curv->basis().begin() + pm1;
  double ucurr, uprev = *(uit++);
  int nCol = n1 - pm1;
  gpar.reserve(nGauss*nCol);
  for (int j = 1; j <= nCol; uit++, j++)
  {
    ucurr = *uit;
    for (int i = 0; i < nGauss; i++)
      gpar.push_back(0.5*((ucurr-uprev)*xg[i] + ucurr+uprev));
    uprev = ucurr;
  }

  Vector N;
  Matrix dNdu, dNdX, Xnod, Jac;
  Vec4   X;


  // === Assembly loop over all elements in the patch ==========================

  int iel = 1;
  for (int i1 = p1; i1 <= n1; i1++, iel++)
  {
    if (MLGE[iel-1] < 1) continue; // zero-length element

    // Check that the current element has nonzero length
    double dL = this->getParametricLength(iel);
    if (dL < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel)) return false;

    // Initialize element matrices
    if (!integrand.initElement(MNPC[iel-1])) return false;

    // Caution: Unless locInt is empty, we assume it points to an array of
    // LocalIntegral pointers, of length at least the number of elements in
    // the model (as defined by the highest number in the MLGE array).
    // If the array is shorter than this, expect a segmentation fault.
    LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];


    // --- Integration loop over all Gauss points in current element -----------

    for (int i = 0; i < nGauss; i++)
    {
      // Weight of current integration point
      double weight = 0.5*dL*wg[i];

      // Compute basis functions and derivatives
      this->extractBasis(gpar((iel-1)*nGauss+i+1),N,dNdu);

      // Compute derivatives in terms of physical co-ordinates
      double detJ = utl::Jacobian(Jac,dNdX,Xnod,dNdu);

      // Cartesian coordinates of current integration point
      X = Xnod * N;
      X.t = time.t;

      // Evaluate the integrand and accumulate element contributions
      if (!integrand.evalInt(elmInt,time,detJ*weight,N,dNdX,X))
	return false;
    }

    // Assembly of global system integral
    if (!glInt.assemble(elmInt,MLGE[iel-1]))
      return false;
  }

  return true;
}


bool ASMs1D::integrate (Integrand& integrand, int lIndex,
			GlobalIntegral& glInt,
			const TimeDomain& time,
			const LintegralVec& locInt)
{
  if (!curv) return true; // silently ignore empty patches

  // Integration of boundary point

  int    iel;
  double param;
  switch (lIndex)
    {
    case 1:
      iel = 1;
      param = curv->startparam();
      break;

    case 2:
      iel = this->getNoElms();
      param = curv->endparam();
      break;

    default:
      return false;
    }

  if (MLGE[iel-1] < 1) return true; // zero-length element

  // Set up control point coordinates for current element
  Matrix Xnod;
  if (!this->getElementCoordinates(Xnod,iel)) return false;

  // Initialize element matrices
  if (!integrand.initElementBou(MNPC[iel-1])) return false;

  // Caution: Unless locInt is empty, we assume it points to an array of
  // LocalIntegral pointers, of length at least the number of elements in
  // the model (as defined by the highest number in the MLGE array).
  // If the array is shorter than this, expect a segmentation fault.
  LocalIntegral* elmInt = locInt.empty() ? 0 : locInt[MLGE[iel-1]-1];

  // Evaluate basis functions and corresponding derivatives
  Vector N;
  Matrix dNdu;
  this->extractBasis(param,N,dNdu);

  // Cartesian coordinates of current integration point
  Vec4 X(Xnod*N,time.t);

  // Compute basis function derivatives
  Matrix Jac, dNdX;
  utl::Jacobian(Jac,dNdX,Xnod,dNdu);

  // Set up the normal vector
  Vec3 normal;
  if (lIndex == 1)
    normal.x = -copysign(1.0,Jac(1,1));
  else
    normal.x = copysign(1.0,Jac(1,1));

  // Evaluate the integrand and accumulate element contributions
  if (!integrand.evalBou(elmInt,time,1.0,N,dNdX,X,normal))
    return false;

  // Assembly of global system integral
  return glInt.assemble(elmInt,MLGE[iel-1]);
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

  DoubleIter uit = curv->basis().begin();
  double ucurr, uprev = *(uit++);
  while (uit != curv->basis().end())
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

  prm.push_back(curv->basis().endparam());
  return true;
}


bool ASMs1D::tesselate (ElementBlock& grid, const int* npe) const
{
  // Compute parameter values of the nodal points
  DoubleVec gpar;
  if (!this->getGridParameters(gpar,npe[0]-1))
    return false;

  // Evaluate the spline curve at all points
  size_t nx = gpar.size();
  DoubleVec XYZ(3*nx,0.0);

  // Establish the block grid coordinates
  size_t i, j;
  for (i = j = 0; i < nx; i++, j += 3)
  {
    Go::Point pt;
    curv->point(pt,gpar[i]);
    for (int k = 0; k < pt.size(); k++)
      XYZ[j+k] = pt[k];
  }

  grid.resize(nx);
  for (i = j = 0; i < grid.getNoNodes(); i++, j += 3)
    grid.setCoor(i,XYZ[j],XYZ[j+1],XYZ[j+2]);

  // Establish the block grid topology
  int n[2], ip = 0;
  n[0] = 0;
  n[1] = n[0] + 1;

  for (i = 1; i < nx; i++)
    for (j = 0; j < 2; j++)
      grid.setNode(ip++,n[j]++);

  return true;
}


/*!
  \brief Auxilliary function for computation of GoTools basis function indices.
*/

static void scatterInd (int p1, int start, IntVec& index)
{
  index.resize(p1);
  for (int i1 = 1; i1 <= p1; i1++)
    index[i1-1] = start-p1+i1;
}


bool ASMs1D::evalSolution (Matrix& sField, const Vector& locSol,
			   const int* npe) const
{
  // Compute parameter values of the result sampling points
  DoubleVec gpar;
  if (!this->getGridParameters(gpar,npe[0]-1))
    return false;

  const int p1 = curv->order();
  size_t nComp = locSol.size() / this->getNoNodes();
  if (nComp*this->getNoNodes() != locSol.size())
    return false;

  Matrix Xtmp;
  Vector Ytmp, basis(p1);

  // Evaluate the primary solution field at each point
  size_t nPoints = gpar.size();
  sField.resize(nComp,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    curv->basis().computeBasisValues(gpar[i],&basis.front());

    IntVec ip;
    scatterInd(p1,curv->basis().lastKnotInterval(),ip);

    utl::gather(ip,nComp,locSol,Xtmp);
    Xtmp.multiply(basis,Ytmp);
    sField.fillColumn(i+1,Ytmp);
  }

  return true;
}


bool ASMs1D::evalSolution (Matrix& sField, const Integrand& integrand,
			   const int* npe) const
{
  sField.resize(0,0);

  // Compute parameter values of the result sampling points
  DoubleVec gpar;
  if (!this->getGridParameters(gpar,npe[0]-1))
    return false;

  const int p1 = curv->order();

  // Fetch nodal (control point) coordinates
  Matrix Xnod, Xtmp;
  this->getNodalCoordinates(Xnod);

  Vector N(p1), solPt;
  Matrix dNdu, dNdX, Jac;

  // Evaluate the secondary solution field at each point
  size_t nPoints = gpar.size();
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch basis function derivatives at current integration point
    this->extractBasis(gpar[i],N,dNdu);

    // Fetch indices of the non-zero basis functions at this point
    IntVec ip;
    scatterInd(p1,curv->basis().lastKnotInterval(),ip);

    // Fetch associated control point coordinates
    utl::gather(ip,nsd,Xnod,Xtmp);

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
