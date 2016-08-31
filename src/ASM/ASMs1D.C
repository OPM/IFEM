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
#include "SplineUtils.h"
#include "Utilities.h"
#include "Vec3Oper.h"


ASMs1D::ASMs1D (unsigned char n_s, unsigned char n_f)
  : ASMstruct(1,n_s,n_f), curv(nullptr), elmCS(myCS), nodalT(myT)
{
}


ASMs1D::ASMs1D (const ASMs1D& patch, unsigned char n_f)
  : ASMstruct(patch,n_f), curv(patch.curv), elmCS(patch.myCS), nodalT(patch.myT)
{
  // Need to set nnod here,
  // as hasXNodes might be invoked before the FE data is generated
  if (nnod == 0 && curv)
    nnod = curv->numCoefs();
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
    curv = nullptr;
    return false;
  }
  else if (curv->dimension() < 1)
  {
    std::cerr <<" *** ASMs1D::read: Invalid spline curve patch, dim="
	      << curv->dimension() << std::endl;
    delete curv;
    curv = nullptr;
    return false;
  }
  else if (curv->dimension() < nsd)
  {
    std::cout <<"  ** ASMs1D::read: The dimension of this curve patch "
	      << curv->dimension() <<" is less than nsd="<< (int)nsd
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
    geo = curv = nullptr;
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
  return this->generateOrientedFEModel(Vec3());
}


bool ASMs1D::generateOrientedFEModel (const Vec3& Zaxis)
{
  if (!curv) return false;

  const int n1 = curv->numCoefs();
  const int p1 = curv->order();

  if (!MLGN.empty())
  {
    nnod = n1;
    if (MLGN.size() != (size_t)nnod)
    {
      std::cerr <<" *** ASMs1D::generateFEMTopology: Inconsistency between the"
                <<" number of FE nodes "<< MLGN.size()
                <<"\n     and the number of spline coefficients "<< nnod
                <<" in the patch."<< std::endl;
      return false;
    }
    nel = n1-p1+1;
    return true;
  }
  else if (shareFE)
    return true;

#ifdef SP_DEBUG
  std::cout <<"numCoefs: "<< n1;
  std::cout <<"\norder: "<< p1;
  std::cout <<"\ndu:";
  for (int j = 0; j < n1; j++)
    std::cout <<' '<< this->getKnotSpan(j);
  std::cout << std::endl;
#endif
  // Consistency checks, just to be fool-proof
  if (n1 <  2) return false;
  if (p1 <  1) return false;
  if (p1 > n1) return false;

  myMLGE.resize(n1-p1+1,0);
  myMLGN.resize(n1);
  myMNPC.resize(myMLGE.size());
  if (nsd == 3 && nf == 6)
  {
    // This is a 3D beam problem, allocate the element/nodal rotation tensors.
    // The nodal rotations are updated during the simulation according to the
    // deformation state, whereas the element tensors are kept constant.
    myCS.resize(n1-p1+1,Tensor(3));
    myT.resize(n1,Tensor(3,true)); // Initialize nodal rotations to unity
  }

  nnod = nel = 0;
  for (int i1 = 1; i1 <= n1; i1++)
  {
    if (i1 >= p1)
    {
      if (this->getKnotSpan(i1-1) > 0.0)
      {
	myMLGE[nel] = ++gEl; // global element number over all patches
	myMNPC[nel].resize(p1,0);

	int lnod = 0;
	for (int j1 = p1-1; j1 >= 0; j1--)
	  myMNPC[nel][lnod++] = nnod - j1;
      }

      nel++;
    }
    myMLGN[nnod++] = ++gNod; // global node number over all patches
  }

#ifdef SP_DEBUG
  std::cout <<"NEL = "<< nel <<" NNOD = "<< nnod << std::endl;
#endif

  return myCS.empty() ? true : this->initLocalElementAxes(Zaxis);
}


bool ASMs1D::initLocalElementAxes (const Vec3& Zaxis)
{
  // Calculate local element axes for 3D beam elements
  for (size_t i = 0; i < myCS.size(); i++)
    if (MLGE[i] > 0)
    {
      Vec3 X1 = this->getCoord(1+MNPC[i].front());
      Vec3 X2 = this->getCoord(1+MNPC[i][curv->order()-1]);
      if (Zaxis.isZero())
        myCS[i] = Tensor(X2-X1,true);
      else
        myCS[i] = Tensor(X2-X1,Zaxis,false,true);
#ifdef SP_DEBUG
      std::cout <<"Local axes for beam element "<< MLGE[i]
                <<", from "<< X1 <<" to "<< X2 <<":\n"<< myCS[i];
#endif
    }

  return true;
}


bool ASMs1D::generateTwistedFEModel (const RealFunc& twist, const Vec3& Zaxis)
{
  if (!this->generateOrientedFEModel(Zaxis))
    return false;

  // Update the local element axes for 3D beam elements
  Tensor rotX(3);
  for (size_t i = 0; i < myCS.size(); i++)
    if (MLGE[i] > 0)
    {
      Vec3 X1 = this->getCoord(1+MNPC[i].front());
      Vec3 X2 = this->getCoord(1+MNPC[i][curv->order()-1]);
      double alpha = twist(0.5*(X1+X2)); // twist angle in the element mid-point
      myCS[i] *= Tensor(alpha*M_PI/180.0,1); // rotate about local X-axis
#ifdef SP_DEBUG
      std::cout <<"Twisted axes for beam element "<< MLGE[i]
                <<", from "<< X1 <<" to "<< X2 <<":\n"<< myCS[i];
#endif
    }

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
  int n1 = this->getSize(basis);
  this->makePeriodic(1,master+n1-1);
}


int ASMs1D::constrainNode (double xi, int dof, int code, char basis)
{
  if (xi < 0.0 || xi > 1.0) return 0;

  int n1 = this->getSize(basis);

  int node = 1;
  for (char i = 1; i < basis; i++)
    node += this->getSize(i);

  if (xi > 0.0) node += int(0.5+(n1-1)*xi);

  this->prescribe(node,dof,code);

  return node;
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

  int inod1 = MNPC[iel-1][curv->order()-1];
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
  int ip = (inod-1)*curv->dimension();
  if (ip < 0) return Vec3();

  return Vec3(&(*(curv->coefs_begin()+ip)),nsd);
}


Tensor ASMs1D::getRotation (size_t inod) const
{
  return inod < 1 || inod > nodalT.size() ? Tensor(nsd,true) : nodalT[inod-1];
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

  X.resize(nsd,curv->order());

  RealArray::const_iterator cit = curv->coefs_begin();
  for (size_t n = 0; n < X.cols(); n++)
  {
    int ip = MNPC[iel-1][n]*curv->dimension();
    if (ip < 0) return false;

    for (size_t i = 0; i < nsd; i++)
      X(i+1,n+1) = *(cit+(ip+i));
  }

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


bool ASMs1D::getElementNodalRotations (TensorVec& T, size_t iel) const
{
#ifdef INDEX_CHECK
  if (iel >= MNPC.size())
  {
    std::cerr <<" *** ASMs1D::getElementNodalRotations: Element index "<< iel
	      <<" out of range [0,"<< MNPC.size() <<">."<< std::endl;
    return false;
  }
#endif

  T.clear();
  if (nodalT.empty())
    return true;

  const IntVec& mnpc = MNPC[iel];
  Tensor Tgl(elmCS[iel],true);

  T.reserve(mnpc.size());
  for (size_t i = 0; i < mnpc.size(); i++)
    T.push_back(Tgl*nodalT[mnpc[i]]);

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


bool ASMs1D::updateRotations (const Vector& displ, bool reInit)
{
  if (shareFE || nf != 6) return true;

  if (displ.size() != 6*myT.size())
  {
    std::cerr <<" *** ASMs1D::updateRotations: Invalid dimension "
	      << displ.size() <<" on displ, should be "
	      << 6*myT.size() << std::endl;
    return false;
  }

  if (reInit)
  {
    if (prevT.empty())
    {
      for (size_t i = 0; i < myT.size(); i++)
        myT[i] = Tensor(displ[6*i+3],displ[6*i+4],displ[6*i+5]);
      return true;
    }
    else
      myT = prevT; // Restore rotation tensors from previous step
  }

  for (size_t i = 0; i < myT.size(); i++)
    myT[i].preMult(Tensor(displ[6*i+3],displ[6*i+4],displ[6*i+5]));

  return true;
}


void ASMs1D::getBoundaryNodes (int lIndex, IntVec& glbNodes,
                               int, bool local) const
{
  if (!curv) return; // silently ignore empty patches

  size_t iel = lIndex == 1 ? 0 : nel-1;
  if (MLGE[iel] > 0)
  {
    int node;
    if (lIndex == 1)
      node = MNPC[iel].front();
    else if (lIndex == 2)
      node = MNPC[iel][curv->order()-1];
    glbNodes.push_back(local ? node : MLGN[node]);
  }
}


std::pair<size_t,double> ASMs1D::findClosestNode (const Vec3& X) const
{
  if (!curv) return std::make_pair(0,-1.0); // silently ignore empty patches

  // Find the closest point on the spline curve
  double param, dist;
  Go::Point Xtarget(X.x,X.y,X.z), Xfound;
  curv->ParamCurve::closestPoint(Xtarget,param,Xfound,dist);

  // Check if point is inside parameter domain
  if (param <= curv->startparam())
    return std::make_pair(1,dist);
  else if (param >= curv->endparam())
    return std::make_pair(this->getNoNodes(),dist);

  // We are inside, now find which knot-span we are in and find closest node
  RealArray::iterator u0 = curv->basis().begin();
  RealArray::iterator u2 = std::lower_bound(u0,curv->basis().end(),param);
  RealArray::iterator u1 = u2-1;

  Go::Point X1, X2;
  curv->point(X1,*u1);
  curv->point(X2,*u2);
  double d1 = X1.dist2(Xfound);
  double d2 = X2.dist2(Xfound);
#ifdef SP_DEBUG
  std::cout <<"ASMs1D::findClosestNode("<< X
            <<"): Found "<< Xfound <<" at u="<< param
            <<" in ["<< *u1 <<","<< *u2
            <<"] d"<< u1-u0 <<"="<< d1 <<" d"<< u2-u0 <<"="<< d2 << std::endl;
#endif

  if (d1 < d2)
    return std::make_pair((u1-u0) - (curv->order()-2), sqrt(d1));
  else
    return std::make_pair((u2-u0) - (curv->order()-2), sqrt(d2));
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


void ASMs1D::getElementEnds (int i, Vec3Vec& XC) const
{
  RealArray::const_iterator uit = curv->basis().begin();

  // Fetch parameter values of the element ends (knots)
  RealArray u(2);
  u[0] = uit[i-1];
  u[1] = uit[i];

  // Evaluate the spline curve at the knots to find physical coordinates
  int dim = curv->dimension();
  RealArray XYZ(dim*2);
  curv->gridEvaluator(XYZ,u);

  XC.clear();
  XC.reserve(elmCS.empty() ? 2 : 3);
  const double* pt = &XYZ.front();
  for (int j = 0; j < 2; j++, pt += dim)
    XC.push_back(Vec3(pt,dim));

  if (elmCS.empty()) return;

  // Add the local Z-axis as the third vector
  int iel = i - curv->order();
  XC.push_back(elmCS[iel][2]);
}


void ASMs1D::extractBasis (double u, Vector& N) const
{
  N.resize(curv->order());
  RealArray basisDerivs;
  curv->computeBasis(u,N,basisDerivs);
}


void ASMs1D::extractBasis (double u, Vector& N, Matrix& dNdu) const
{
  int p1 = curv->order();

  N.resize(p1);
  dNdu.resize(p1,1);

  RealArray basisDerivs;
  curv->computeBasis(u,N,basisDerivs);
  dNdu.fillColumn(1,basisDerivs);
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
      N(i)        = bas[3*i-3];
     dNdu(i,1)    = bas[3*i-2];
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
  const double* xr = nullptr;
  const double* wr = nullptr;
  int nRed = integrand.getReducedIntegration(nGauss);
  if (nRed > 0)
  {
    xr = GaussQuadrature::getCoord(nRed);
    wr = GaussQuadrature::getWeight(nRed);
    if (!xr || !wr) return false;
  }
  else if (nRed < 0)
    nRed = nGauss; // The integrand needs to know nGauss

  if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
    if (curv->rational())
    {
      std::cerr <<" *** ASMs1D::integrate: Second-derivatives of NURBS "
                <<" is not implemented yet, sorry..."<< std::endl;
      return false;
    }

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gpar, redpar;
  this->getGaussPointParameters(gpar,nGauss,xg);
  if (xr)
    this->getGaussPointParameters(redpar,nRed,xr);

  const int p1 = curv->order();

  FiniteElement fe(p1);
  Matrix   dNdu, Jac;
  Matrix3D d2Ndu2, Hess;
  Vec4     X;

  if (nsd > 1 && (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES))
    fe.G.resize(nsd,2); // For storing d{X}/du and d2{X}/du2


  // === Assembly loop over all elements in the patch ==========================

  for (size_t iel = 0; iel < nel; iel++)
  {
    fe.iel = MLGE[iel];
    if (fe.iel < 1) continue; // zero-length element

    // Check that the current element has nonzero length
    double dL = this->getParametricLength(1+iel);
    if (dL < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(fe.Xn,1+iel)) return false;

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      this->getElementEnds(p1+iel,fe.XC);

    if (integrand.getIntegrandType() & Integrand::NODAL_ROTATIONS)
    {
      this->getElementNodalRotations(fe.Tn,iel);
      if (!elmCS.empty()) fe.Te = elmCS[iel];
    }

    // Initialize element matrices
    LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
    bool ok = integrand.initElement(MNPC[iel],fe,X,nRed,*A);

    if (xr)
    {
      // --- Selective reduced integration loop --------------------------------

      for (int i = 0; i < nRed && ok; i++)
      {
	// Local element coordinates of current integration point
	fe.xi = xr[i];

	// Parameter values of current integration point
	fe.u = redpar(1+i,1+iel);

        if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
          this->extractBasis(fe.u,fe.N);
        else
        {
          // Fetch basis function derivatives at current point
          this->extractBasis(fe.u,fe.N,dNdu);
          // Compute Jacobian inverse and derivatives
          dNdu.multiply(0.5*dL); // Derivatives w.r.t. xi=[-1,1]
          fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu)*wr[i];
        }

	// Cartesian coordinates of current integration point
	X = fe.Xn * fe.N;
	X.t = time.t;

	// Compute the reduced integration terms of the integrand
	ok = integrand.reducedInt(*A,fe,X);
      }
    }


    // --- Integration loop over all Gauss points in current element -----------

    int jp = iel*nGauss;
    fe.iGP = firstIp + jp; // Global integration point counter

    for (int i = 0; i < nGauss && ok; i++, fe.iGP++)
    {
      // Local element coordinate of current integration point
      fe.xi = xg[i];

      // Parameter value of current integration point
      fe.u = gpar(1+i,1+iel);

      // Compute basis functions and derivatives
      if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
        this->extractBasis(fe.u,fe.N);
      else if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
        this->extractBasis(fe.u,fe.N,dNdu,d2Ndu2);
      else
        this->extractBasis(fe.u,fe.N,dNdu);

      if (!dNdu.empty())
      {
        // Compute derivatives in terms of physical coordinates
        dNdu.multiply(0.5*dL); // Derivatives w.r.t. xi=[-1,1]
        fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu)*wg[i];
        if (fe.detJxW == 0.0) continue; // skip singular points

        // Compute Hessian of coordinate mapping and 2nd order derivatives
        if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
        {
          d2Ndu2.multiply(0.25*dL*dL); // 2nd derivatives w.r.t. xi=[-1,1]
          if (!utl::Hessian(Hess,fe.d2NdX2,Jac,fe.Xn,d2Ndu2,fe.dNdX))
            ok = false;
          else if (fe.G.cols() == 2)
          {
            // Store the first and second derivatives of {X} w.r.t.
            // the parametric coordinate (xi), in the G-matrix
            fe.G.fillColumn(1,Jac.ptr());
            fe.G.fillColumn(2,Hess.ptr());
          }
        }
      }

      // Cartesian coordinates of current integration point
      X = fe.Xn * fe.N;
      X.t = time.t;

      // Evaluate the integrand and accumulate element contributions
      if (ok && !integrand.evalInt(*A,fe,time,X))
        ok = false;
    }

    // Finalize the element quantities
    if (ok && !integrand.finalizeElement(*A,fe,time,firstIp+jp))
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
  size_t iel = 0;
  switch (lIndex)
    {
    case 1:
      fe.xi = -1.0;
      fe.u = curv->startparam();
      break;

    case 2:
      fe.xi = 1.0;
      fe.u = curv->endparam();
      iel = nel-1;
      break;

    default:
      return false;
    }

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  fe.iGP = iit == firstBp.end() ? 0 : iit->second;
  fe.iel = MLGE[iel];
  if (fe.iel < 1) return true; // zero-length element

  // Set up control point coordinates for current element
  if (!this->getElementCoordinates(fe.Xn,1+iel)) return false;

  if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
    this->getElementEnds(iel+curv->order(),fe.XC);

  if (integrand.getIntegrandType() & Integrand::NODAL_ROTATIONS)
  {
    this->getElementNodalRotations(fe.Tn,iel);
    if (!elmCS.empty()) fe.Te = elmCS[iel];
  }

  // Initialize element matrices
  LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
  bool ok = integrand.initElementBou(MNPC[iel],*A);

  Vec3 normal;

  // Evaluate basis functions and corresponding derivatives
  if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
    this->extractBasis(fe.u,fe.N);
  else
  {
    // Compute basis function derivatives
    Matrix dNdu, Jac;
    this->extractBasis(fe.u,fe.N,dNdu);
    utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu);

    // Set up the normal vector
    if (lIndex == 1)
      normal.x = -copysign(1.0,Jac(1,1));
    else
      normal.x = copysign(1.0,Jac(1,1));
  }

  // Cartesian coordinates of current integration point
  Vec4 X(fe.Xn*fe.N,time.t);

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
  SplineUtils::point(X,param[0],curv);

  // Check if this point matches any of the control points (nodes)
  return this->searchCtrlPt(curv->coefs_begin(),curv->coefs_end(),
                            X,curv->dimension());
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
                           const RealArray* gpar, bool, int deriv) const
{
  const int p1 = curv->order();
  size_t nComp = locSol.size() / curv->numCoefs();

  Vector   basis(p1), ptSol;
  Matrix   dNdu, dNdX, Xnod, Xtmp, Jac, eSol, ptDer;
  Matrix3D d2Ndu2, d2NdX2, Hess, ptDer2;

  // Fetch nodal (control point) coordinates
  this->getNodalCoordinates(Xnod);

  // Evaluate the primary solution field at each point
  const RealArray& upar = *gpar;
  size_t nPoints = upar.size();
  sField.resize(nComp*int(pow(nsd,deriv)),nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    IntVec ip;
    switch (deriv) {

    case 0: // Evaluate the solution
      this->extractBasis(upar[i],basis);
      scatterInd(p1,curv->basis().lastKnotInterval(),ip);
      utl::gather(ip,nComp,locSol,Xtmp);
      Xtmp.multiply(basis,ptSol);
      sField.fillColumn(1+i,ptSol);
      break;

    case 1: // Evaluate first derivatives of the solution
      this->extractBasis(upar[i],basis,dNdu);
      scatterInd(p1,curv->basis().lastKnotInterval(),ip);
      utl::gather(ip,nsd,Xnod,Xtmp);
      utl::Jacobian(Jac,dNdX,Xtmp,dNdu);
      utl::gather(ip,nComp,locSol,Xtmp);
      ptDer.multiply(Xtmp,dNdX);
      sField.fillColumn(1+i,ptDer);
      break;

    case 2: // Evaluate second derivatives of the solution
      this->extractBasis(upar[i],basis,dNdu,d2Ndu2);
      scatterInd(p1,curv->basis().lastKnotInterval(),ip);
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
    return nullptr;

  // Evaluate the secondary solution at all sampling points
  Matrix sValues;
  if (!this->evalSolution(sValues,integrand,&gpar))
    return nullptr;

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
  FiniteElement fe(p1,firstIp);
  this->getNodalCoordinates(fe.Xn);

  Vector   solPt;
  Matrix   dNdu, Jac, Xtmp;
  Matrix3D d2Ndu2, Hess;

  if (nsd > 1 && (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES))
    fe.G.resize(nsd,2); // For storing d{X}/du and d2{X}/du2

  // Evaluate the secondary solution field at each point
  const RealArray& upar = *gpar;
  size_t nPoints = upar.size();
  for (size_t i = 0; i < nPoints; i++, fe.iGP++)
  {
    fe.u = upar[i];

    // Fetch basis function derivatives at current integration point
    if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
      this->extractBasis(fe.u,fe.N);
    else if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
      this->extractBasis(fe.u,fe.N,dNdu,d2Ndu2);
    else
      this->extractBasis(fe.u,fe.N,dNdu);

    // Fetch indices of the non-zero basis functions at this point
    IntVec ip;
    scatterInd(p1,curv->basis().lastKnotInterval(),ip);

    // Fetch associated control point coordinates
    utl::gather(ip,nsd,fe.Xn,Xtmp);

    if (!dNdu.empty())
    {
      // Compute the Jacobian inverse and derivatives
      fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xtmp,dNdu);

      // Compute Hessian of coordinate mapping and 2nd order derivatives
      if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
      {
        if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xtmp,d2Ndu2,fe.dNdX))
          continue;
        else if (fe.G.cols() == 2)
        {
          // Store the first and second derivatives of {X} w.r.t.
          // the parametric coordinate (xi), in the G-matrix
          fe.G.fillColumn(1,Jac.ptr());
          fe.G.fillColumn(2,Hess.ptr());
        }
      }
    }

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
				 const IntegrandBase& integrand,
				 bool continuous) const
{
  if (!curv) return true; // silently ignore empty patches

  if (continuous)
  {
    std::cerr <<" *** ASMs1D::globalL2projection: Only discrete L2-projection"
              <<" is available in this method.\n                            "
              <<"     Use ASMbase::L2projection instead."<< std::endl;
  }

  const int p1 = curv->order();

  // Get Gaussian quadrature point coordinates
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

  Vector phi(p1);


  // === Assembly loop over all elements in the patch ==========================

  int ip = 0;
  for (size_t iel = 0; iel < nel; iel++)
  {
    if (MLGE[iel] < 1) continue; // zero-area element

    // --- Integration loop over all Gauss points in current element -----------

    for (int i = 1; i <= ng1; i++, ip++)
    {
      // Fetch basis function values at current integration point
      this->extractBasis(gp(i,1+iel),phi);

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


bool ASMs1D::getNoStructElms (int& n1, int& n2, int& n3) const
{
  n1 = nel;
  n2 = n3 = 0;

  return true;
}


bool ASMs1D::evaluate (const RealFunc* func, RealArray& values,
                       int, double time) const
{
  Go::SplineCurve* scrv = SplineUtils::project(curv,*func,time);
  if (!scrv)
  {
    std::cerr <<" *** ASMs1D::evaluate: Projection failure."<< std::endl;
    return false;
  }

  values.assign(scrv->coefs_begin(),scrv->coefs_end());
  delete scrv;

  return true;
}
