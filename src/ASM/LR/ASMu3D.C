// $Id$
//==============================================================================
//!
//! \file ASMu3D.C
//!
//! \date January 2013
//!
//! \author Kjetil A. Johannessen / NTNU
//!
//! \brief Driver for assembly of unstructured 3D spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariate/SplineVolume.h"

#include "LRSpline/LRSplineVolume.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

#include "ASMu3D.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "LagrangeInterpolator.h"
#include "ElementBlock.h"
#include "MPC.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Function.h"
#include "Vec3Oper.h"
#include <array>


ASMu3D::ASMu3D (unsigned char n_f)
  : ASMLRSpline(3,3,n_f), lrspline(nullptr), tensorspline(nullptr),
    myGeoBasis(1), bezierExtract(myBezierExtract)
{
  vMin = 0.0;
}


ASMu3D::ASMu3D (const ASMu3D& patch, unsigned char n_f)
  : ASMLRSpline(patch,n_f), lrspline(patch.lrspline), tensorspline(nullptr),
    myGeoBasis(1), bezierExtract(patch.myBezierExtract)
{
  vMin = 0.0;

  // Need to set nnod here,
  // as hasXNodes might be invoked before the FE data is generated
  if (nnod == 0 && lrspline)
    nnod = lrspline->nBasisFunctions();
}


bool ASMu3D::read (std::istream& is)
{
  if (shareFE) return true;

  // read inputfile as either an LRSpline file directly
  // or a tensor product B-spline and convert
  char firstline[256];
  is.getline(firstline, 256);
  if (strncmp(firstline, "# LRSPLINE", 10) == 0) {
    lrspline.reset(new LR::LRSplineVolume());
    is >> *lrspline;
  } else { // probably a SplineVolume, so we'll read that and convert
    tensorspline = new Go::SplineVolume();
    is >> *tensorspline;
    lrspline.reset(new LR::LRSplineVolume(tensorspline));
  }

  // Eat white-space characters to see if there is more data to read
  char c;
  while (is.get(c))
    if (!isspace(c)) {
      is.putback(c);
      break;
    }

  if (!is.good() && !is.eof())
  {
    std::cerr <<" *** ASMu3D::read: Failure reading spline data"<< std::endl;
    lrspline.reset();
    return false;
  }
  else if (lrspline->dimension() < 3)
  {
    std::cerr <<" *** ASMu3D::read: Invalid spline volume patch, dim="
              << lrspline->dimension() << std::endl;
    lrspline.reset();
    return false;
  }

  geo = lrspline.get();
  return true;
}


bool ASMu3D::write (std::ostream& os, int) const
{
  if (!lrspline) return false;

  os << *lrspline;

  return os.good();
}


void ASMu3D::clear (bool retainGeometry)
{
  if (!retainGeometry) {
    // Erase spline data
    if (!shareFE) {
      lrspline.reset();
      delete tensorspline;
    }
    geo = nullptr;
    tensorspline = nullptr;
  }

  // Erase the FE data
  this->ASMbase::clear(retainGeometry);
  this->dirich.clear();
}


bool ASMu3D::uniformRefine (int dir, int nInsert)
{
  if (!tensorspline || dir < 0 || dir > 2 || nInsert < 1) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = tensorspline->basis(dir).begin();
  double uprev = *(uit++);
  while (uit != tensorspline->basis(dir).end())
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

  tensorspline->insertKnot(dir,extraKnots);
  lrspline.reset(new LR::LRSplineVolume(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu3D::refine (int dir, const RealArray& xi)
{
  if (!tensorspline || dir < 0 || dir > 2 || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = tensorspline->basis(dir).begin();
  double uprev = *(uit++);
  while (uit != tensorspline->basis(dir).end())
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

  tensorspline->insertKnot(dir,extraKnots);
  lrspline.reset(new LR::LRSplineVolume(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu3D::raiseOrder (int ru, int rv, int rw)
{
  if (!tensorspline) return false;
  if (shareFE) return true;

  tensorspline->raiseOrder(ru,rv,rw);
  lrspline.reset(new LR::LRSplineVolume(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu3D::generateFEMTopology ()
{
  // At this point we are through with the tensor spline object,
  // so release it to avoid memory leakage
  delete tensorspline;
  tensorspline = nullptr;

  if (!lrspline) return false;

  nnod = lrspline->nBasisFunctions();
  nel  = lrspline->nElements();

  if (!MLGN.empty()) {
    if (MLGN.size() != nnod)
    {
      std::cerr <<" *** ASMu3D::generateFEMTopology: Inconsistency"
                <<" between the number of FE nodes "<< MLGN.size()
                <<"\n     and the number of basis functions "<< nnod
                <<" in the patch."<< std::endl;
      return false;
    }
    return true;
  }

  if (shareFE) return true;

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);
  const int p3 = lrspline->order(2);

  myMLGN.resize(nnod);
  myMLGE.resize(nel);
  myMNPC.resize(nel);

  myBezierExtract.resize(nel);
  lrspline->generateIDs();

  size_t iel = 0;
  RealArray extrMat;
  for (const LR::Element* elm : lrspline->getAllElements())
  {
    myMLGE[iel] = ++gEl; // global element number over all patches
    myMNPC[iel].resize(elm->nBasisFunctions());

    int lnod = 0;
    for (LR::Basisfunction* b : elm->support())
      myMNPC[iel][lnod++] = b->getId();

    // Get bezier extraction matrix
    PROFILE("Bezier extraction");
    lrspline->getBezierExtraction(iel,extrMat);
    myBezierExtract[iel].resize(elm->nBasisFunctions(),p1*p2*p3);
    myBezierExtract[iel++].fill(extrMat.data(),extrMat.size());
  }

  for (size_t inod = 0; inod < nnod; inod++)
    myMLGN[inod] = ++gNod;

  return true;
}


bool ASMu3D::connectPatch (int face, ASM3D& neighbor, int nface,
                           int norient, int, bool coordCheck, int thick)
{
  ASMu3D* neighU = dynamic_cast<ASMu3D*>(&neighbor);
  if (!neighU)
    return false;

  if (!this->connectBasis(face,*neighU,nface,norient,1,0,0,coordCheck,thick))
    return false;

  this->addNeighbor(neighU);
  return true;
}


bool ASMu3D::connectBasis (int face, ASMu3D& neighbor, int nface, int norient,
                           int basis, int slave, int master,
                           bool coordCheck, int /*thick*/)
{
  if (this->isShared() && neighbor.isShared())
    return true;
  else if (this->isShared() || neighbor.isShared())
  {
    std::cerr <<" *** ASMu3D::connectBasis: Logic error, cannot"
              <<" connect a shared patch with an unshared one"<< std::endl;
    return false;
  }

  // Set up the slave node numbers for this volume patch
  IntVec slaveNodes = this->getFaceNodes(face, basis, 0);
  for (int& it : slaveNodes)
    it += slave;

  // Set up the master node numbers for the neighboring volume patch
  IntVec masterNodes = neighbor.getFaceNodes(nface, basis, norient);
  for (int& it : masterNodes)
    it += master;

  if (masterNodes.empty() || masterNodes.size() != slaveNodes.size())
  {
    std::cerr <<" *** ASMu3D::connectBasis: Non-matching faces, sizes "
              << masterNodes.size() <<" and "<< slaveNodes.size() << std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  IntVec::const_iterator slvIt = slaveNodes.begin();
  for (int master : masterNodes)
    if (!coordCheck)
      ASMbase::collapseNodes(neighbor,master,*this,*(slvIt++));
    else if (neighbor.getCoord(master).equal(this->getCoord(*slvIt),xtol))
      ASMbase::collapseNodes(neighbor,master,*this,*(slvIt++));
    else
    {
      std::cerr <<" *** ASMu3D::connectBasis: Non-matching nodes "
                << master <<": "<< neighbor.getCoord(master)
                <<"\n                                          and "
                << *slvIt <<": "<< this->getCoord(*slvIt) << std::endl;
      return false;
    }

  return true;
}


/*!
  A negative \a code value implies direct evaluation of the Dirichlet condition
  function at the control point. Positive \a code implies projection onto the
  spline basis representing the boundary surface (needed for curved faces and/or
  non-constant functions).
*/

void ASMu3D::constrainFace (int dir, bool open, int dof, int code, char basis)
{
  if (basis < 1) basis = 1;

  // figure out function index offset (when using multiple basis)
  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  // figure out what edge we are at
  LR::parameterEdge face;
  switch (dir)
  {
  case -1: face = LR::WEST;   break;
  case  1: face = LR::EAST;   break;
  case -2: face = LR::SOUTH;  break;
  case  2: face = LR::NORTH;  break;
  case -3: face = LR::BOTTOM; break;
  case  3: face = LR::TOP;    break;
  default: return;
  }

  // fetch the right basis to consider
  LR::LRSplineVolume* lr = this->getBasis(basis);

  // get all elements and functions on this face
  std::vector<LR::Basisfunction*> faceFunctions;
  std::vector<LR::Element*>       faceElements;
  lr->getEdgeFunctions(faceFunctions,face);
  lr->getEdgeElements (faceElements ,face);

  // build up the local element/node correspondence needed by the projection
  // call on this face by ASMu3D::updateDirichlet()
  DirichletFace de(faceFunctions.size(), faceElements.size(), dof, code, basis);
  de.edg  = face;
  de.lr   = lr;

  // find the corners since these are not to be included in the L2-fitting
  // of the inhomogenuous dirichlet boundaries; corners are interpolatory.
  // Optimization note: loop over the "edge"-container to manually pick up
  // the end nodes. LRspline::getEdgeFunctions() does a global search.
  static const std::map<LR::parameterEdge,std::array<std::array<int,3>,4>> corners = {{
        { LR::WEST,{{{{-1, -1, -1}},
                     {{-1,  1, -1}},
                     {{-1, -1,  1}},
                     {{-1,  1,  1}}}}},
        { LR::EAST,{{{{ 1, -1, -1}},
                     {{ 1,  1, -1}},
                     {{ 1, -1,  1}},
                     {{ 1,  1,  1}}}}},
        {LR::SOUTH,{{{{-1, -1, -1}},
                     {{ 1, -1, -1}},
                     {{-1, -1,  1}},
                     {{ 1, -1,  1}}}}},
        {LR::NORTH,{{{{-1, -1, -1}},
                     {{ 1, -1, -1}},
                     {{-1, -1,  1}},
                     {{ 1, -1,  1}}}}},
       {LR::BOTTOM,{{{{-1, -1, -1}},
                     {{ 1, -1, -1}},
                     {{-1,  1, -1}},
                     {{ 1,  1, -1}}}}},
          {LR::TOP,{{{{-1, -1,  1}},
                     {{ 1, -1,  1}},
                     {{-1,  1,  1}},
                     {{ 1,  1,  1}}}}},
  }};

  int* c = de.corners;
  for (const std::array<int,3>& vx : corners.find(face)->second)
    *c++ = this->getCorner(vx[0], vx[1], vx[2], basis);

  int bcode = abs(code);
  for (LR::Basisfunction* b : faceFunctions)
  {
    de.MLGN.push_back(b->getId());
    int* cid = std::find(de.corners, de.corners+4, b->getId()+ofs);
    // skip corners for open boundaries
    if (open && cid == de.corners+4)
      continue;
    else
    {
      // corners are interpolated (positive 'code')
      if (cid != de.corners+4)
        this->prescribe(b->getId()+ofs, dof,  bcode);
      // inhomogenuous dirichlet conditions by function evaluation (negative 'code')
      else if (code > 0)
        this->prescribe(b->getId()+ofs, dof, -bcode);
      // (in)homogenuous constant dirichlet conditions
      else
        this->prescribe(b->getId()+ofs, dof,  bcode);
    }
  }

  // build MLGE and MNPC matrix
  for (size_t i=0; i < faceElements.size(); i++)
  {
    LR::Element* el = faceElements[i];

    // for mixed FEM models, let MLGE point to the *geometry* index
    if (de.lr != this->lrspline.get())
    {
      double umid = (el->umax() + el->umin())/2.0;
      double vmid = (el->vmax() + el->vmin())/2.0;
      double wmid = (el->wmax() + el->wmin())/2.0;
      de.MLGE[i] = lrspline->getElementContaining(umid, vmid, wmid);
    }
    else
    {
      de.MLGE[i] = el->getId();
    }
    for (LR::Basisfunction* b : el->support())
    {
      de.MNPC[i].push_back(-1);
      for (size_t j = 0; j < de.MLGN.size(); j++)
        if (b->getId() == de.MLGN[j])
          de.MNPC[i].back() = j;
    }
  }
  if (code > 0)
    dirich.push_back(de);
}

size_t ASMu3D::constrainFaceLocal(int dir, bool open, int dof, int code, bool project, char T1)
{
  std::cerr << "ASMu3D::constrainFaceLocal not implemented properly yet" << std::endl;
  exit(776654);
  return 0;
}


IntVec ASMu3D::getEdge (int lEdge, bool open, int basis, int orient) const
{
  // lEdge = 1-4, running index is u (vmin,wmin), (vmax,wmin), (vmin,wmax), (vmax,wmax)
  // lEdge = 5-8, running index is v (umin,wmin), (umax,wmin), (umin,wmax), (umax,wmax)
  // lEdge = 9-12, running index is w

  int edge = LR::NONE;
  if (lEdge == 1)
    edge = LR::BOTTOM | LR::SOUTH;
  else if (lEdge == 2)
    edge = LR::BOTTOM | LR::NORTH;
  else if (lEdge == 3)
    edge = LR::TOP    | LR::SOUTH;
  else if (lEdge == 4)
    edge = LR::TOP    | LR::NORTH;
  else if (lEdge == 5)
    edge = LR::BOTTOM | LR::WEST;
  else if (lEdge == 6)
    edge = LR::BOTTOM | LR::EAST;
  else if (lEdge == 7)
    edge = LR::TOP    | LR::WEST;
  else if (lEdge == 8)
    edge = LR::TOP    | LR::EAST;
  else if (lEdge == 9)
    edge = LR::SOUTH  | LR::WEST;
  else if (lEdge == 10)
    edge = LR::SOUTH  | LR::EAST;
  else if (lEdge == 11)
    edge = LR::NORTH  | LR::WEST;
  else if (lEdge == 12)
    edge = LR::NORTH  | LR::EAST;

  // figure out function index offset (when using multiple basis)
  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  // get all the boundary functions from the LRspline object
  std::vector<LR::Basisfunction*> thisEdge;
  this->getBasis(basis)->getEdgeFunctions(thisEdge, (LR::parameterEdge) edge, 1);
  if (orient > -1) {
    int dir = (edge-1)/4;
    int u = dir == 0;
    int v = 1 + (dir != 2);
    ASMLRSpline::Sort(u, v, orient, thisEdge);
  }

  IntVec result;
  for (LR::Basisfunction* b : thisEdge)
    result.push_back(b->getId()+ofs);

  return result;
}

void ASMu3D::constrainEdge (int lEdge, bool open, int dof, int code, char basis)
{
  if (open)
    std::cout <<"  ** ASMu3D::constrainEdge: Open boundary conditions are not"
              <<" supported for LR B-splines. Treated as closed."<< std::endl;

  if (basis < 1) basis = 1;

  // enforce the boundary conditions
  for (int node : this->getEdge(lEdge, open, basis, -1))
    this->prescribe(node,dof,code);
}


void ASMu3D::constrainLine (int, int, double, int, int, char)
{
  // We can't do this, since the LR meshes are unstructured by nature
  std::cout <<"  ** Constraining an interior line is not available"
            <<" for LR B-splines (ignored)."<< std::endl;
}


void ASMu3D::constrainCorner (int I, int J, int K, int dof,
                              int code, char basis)
{
  if (basis < 1) basis = 1;

  int corner = LR::NONE;
  corner |= (I < 0 ? LR::WEST   : LR::EAST );
  corner |= (J < 0 ? LR::SOUTH  : LR::NORTH);
  corner |= (K < 0 ? LR::BOTTOM : LR::TOP  );

  std::vector<LR::Basisfunction*> one_function;
  this->getBasis(basis)->getEdgeFunctions(one_function,
                                          LR::parameterEdge(corner));

  this->prescribe(one_function.front()->getId()+1, dof, code);
}


void ASMu3D::constrainNode (double xi, double eta, double zeta,
                            int dof, int code)
{
  // We can't do this, since the control point locations are unpredictable
  std::cout <<"  ** Constraining a nodal point is not available"
            <<" for LR B-splines (ignored)."<< std::endl;
}


#define DERR -999.99

double ASMu3D::getParametricArea (int iel, int dir) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMu3D::getParametricArea: Element index "<< iel
        <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  const LR::Element* el = lrspline->getElement(iel-1);
  double du = el->getParmax(0) - el->getParmin(0);
  double dv = el->getParmax(1) - el->getParmin(1);
  double dw = el->getParmax(2) - el->getParmin(2);
  switch (dir)
    {
    case 1: return dv * dw;
    case 2: return du * dw;
    case 3: return du * dv;
    }

  std::cerr <<" *** ASMu3D::getParametricArea: Invalid face direction "
            << dir << std::endl;
  return DERR;
}


double ASMu3D::getParametricVolume (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMu3D::getParametricArea: Element index "<< iel
              <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return DERR;
  }
#endif
  if (MNPC[iel-1].empty())
    return 0.0;

  const LR::Element* el = lrspline->getElement(iel-1);
  double du = el->getParmax(0) - el->getParmin(0);
  double dv = el->getParmax(1) - el->getParmin(1);
  double dw = el->getParmax(2) - el->getParmin(2);
  return du*dv*dw;
}


bool ASMu3D::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMu3D::getElementCoordinates: Element index "<< iel
              <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif

  const LR::Element* el = lrspline->getElement(iel-1);
  X.resize(3,el->nBasisFunctions());

  int n = 1;
  for (LR::Basisfunction* b : el->support())
    X.fillColumn(n++,&(*b->cp()));

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


void ASMu3D::getNodalCoordinates (Matrix& X) const
{
  X.resize(3,lrspline->nBasisFunctions());

  size_t inod = 1;
  for (LR::Basisfunction* b : lrspline->getAllBasisfunctions())
    X.fillColumn(inod++,&(*b->cp()));
}


Vec3 ASMu3D::getCoord (size_t inod) const
{
  LR::Basisfunction* basis = lrspline->getBasisfunction(inod-1);

  RealArray::const_iterator cit = basis->cp();
  return Vec3(*cit,*(cit+1),*(cit+2));
}


bool ASMu3D::updateCoords (const Vector& displ)
{
  std::cerr <<" *** ASMu3D::updateCoords: Not implemented!"<< std::endl;
  return false;
}


size_t ASMu3D::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (!lrspline) return 0;

  LR::parameterEdge edge;
  switch (lIndex)
  {
  case 1: edge = LR::WEST;   break;
  case 2: edge = LR::EAST;   break;
  case 3: edge = LR::SOUTH;  break;
  case 4: edge = LR::NORTH;  break;
  case 5: edge = LR::BOTTOM; break;
  case 6: edge = LR::TOP;    break;
  default:edge = LR::NONE;
  }

  std::vector<LR::Element*> edgeElms;
  lrspline->getEdgeElements(edgeElms, edge);

  return edgeElms.size();
}


void ASMu3D::getGaussPointParameters (RealArray& uGP, int dir, int nGauss,
                                      int iEl, const double* xi) const
{
  const LR::Element* el = lrspline->getElement(iEl-1);
  double start = el->getParmin(dir);
  double stop  = el->getParmax(dir);

  uGP.resize(nGauss);
  for (int i = 0; i < nGauss; i++)
    uGP[i] = 0.5*((stop-start)*xi[i] + stop+start);
}


double ASMu3D::getElementCorners (int iEl, Vec3Vec& XC) const
{
  const LR::Element* el = lrspline->getElement(iEl-1);
  double u[2] = { el->getParmin(0), el->getParmax(0) };
  double v[2] = { el->getParmin(1), el->getParmax(1) };
  double w[2] = { el->getParmin(2), el->getParmax(2) };

  XC.clear();
  XC.reserve(8);
  Go::Point pt;

  for (int k = 0; k < 2; k++)
    for (int j = 0; j < 2; j++)
      for (int i = 0; i < 2; i++)
      {
        lrspline->point(pt,u[i],v[j],w[k],iEl-1);
        XC.push_back(SplineUtils::toVec3(pt));
      }

  return getElementSize(XC);
}


void ASMu3D::evaluateBasis (int iel, int basis, double u, double v, double w,
                            Vector& N, Matrix& dNdu) const
{
  PROFILE3("ASMu3D::evalBasis(1)");

  std::vector<RealArray> result;
  this->getBasis(basis)->computeBasis(u, v, w, result, 1, iel);
  size_t n = 0, nBasis = result.size();

  N.resize(nBasis);
  dNdu.resize(nBasis,3);
  for (const RealArray& phi : result) {
    N   (++n) = phi[0];
    dNdu(n,1) = phi[1];
    dNdu(n,2) = phi[2];
    dNdu(n,3) = phi[3];
  }
}

void ASMu3D::evaluateBasis (int iel, FiniteElement& fe, Matrix& dNdu,
                            int basis) const
{
  this->evaluateBasis(iel, basis, fe.u, fe.v, fe.w, fe.basis(basis), dNdu);
}

void ASMu3D::evaluateBasis (FiniteElement& fe, Matrix& dNdu,
                            const Matrix& C, const Matrix& B, int basis) const
{
  PROFILE3("ASMu3D::evalBasis(BE)");

  Matrix N = C*B;
  dNdu.resize(N.rows(),3);
  fe.basis(basis) = N.getColumn(1);
  dNdu.fillColumn(1,N.getColumn(2));
  dNdu.fillColumn(2,N.getColumn(3));
  dNdu.fillColumn(3,N.getColumn(4));
}

void ASMu3D::evaluateBasis (int iel, FiniteElement& fe,
                            Matrix& dNdu, Matrix3D& d2Ndu2,
                            int basis) const
{
  PROFILE3("ASMu3D::evalBasis(2)");

  std::vector<RealArray> result;
  this->getBasis(basis)->computeBasis(fe.u, fe.v, fe.w, result, 2, iel);
  size_t n = 0, nBasis = result.size();
  Vector& N = fe.basis(basis);

  N.resize(nBasis);
  dNdu.resize(nBasis,3);
  d2Ndu2.resize(nBasis,3,3);
  for (const RealArray& phi : result) {
    N     (++n)   = phi[0];
    dNdu  (n,1)   = phi[1];
    dNdu  (n,2)   = phi[2];
    dNdu  (n,3)   = phi[3];
    d2Ndu2(n,1,1) = phi[4];
    d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = phi[5];
    d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = phi[6];
    d2Ndu2(n,2,2) = phi[7];
    d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = phi[8];
    d2Ndu2(n,3,3) = phi[9];
  }
}

void ASMu3D::evaluateBasis (int iel, FiniteElement& fe, int derivs,
                            int basis) const
{
  PROFILE3("ASMu2D::evalBasis");

  std::vector<RealArray> result;
  this->getBasis(basis)->computeBasis(fe.u, fe.v, fe.w, result, derivs, iel);
  size_t n = 0, nBasis = result.size();
  Vector& N = fe.basis(basis);
  Matrix& dNdu = fe.grad(basis);
  Matrix3D& d2Ndu2 = fe.hess(basis);

  N.resize(nBasis);
  if (derivs > 0)
    dNdu.resize(nBasis,3);
  if (derivs > 1)
    d2Ndu2.resize(nBasis,3,3);
  for (const RealArray& phi : result) {
    N(++n) = phi[0];
    if (derivs > 0) {
      dNdu(n,1) = phi[1];
      dNdu(n,2) = phi[2];
      dNdu(n,3) = phi[3];
    }
    if (derivs > 1) {
      d2Ndu2(n,1,1) = phi[4];
      d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = phi[5];
      d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = phi[6];
      d2Ndu2(n,2,2) = phi[7];
      d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = phi[8];
      d2Ndu2(n,3,3) = phi[9];
    }
  }
}


bool ASMu3D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu3D::integrate(I)");

  int p1 = lrspline->order(0);
  int p2 = lrspline->order(1);
  int p3 = lrspline->order(2);
  int pm = std::max(std::max(p1,p2),p3);

  // Get Gaussian quadrature points and weights
  int nGP = this->getNoGaussPt(pm);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Evaluate all gauss points on the bezier patch (-1, 1)
  double u[2*p1], v[2*p2], w[2*p3];
  Go::BsplineBasis basis1 = getBezierBasis(p1);
  Go::BsplineBasis basis2 = getBezierBasis(p2);
  Go::BsplineBasis basis3 = getBezierBasis(p3);
  Matrix BN(   p1*p2*p3, nGP*nGP*nGP), rBN;
  Matrix BdNdu(p1*p2*p3, nGP*nGP*nGP), rBdNdu;
  Matrix BdNdv(p1*p2*p3, nGP*nGP*nGP), rBdNdv;
  Matrix BdNdw(p1*p2*p3, nGP*nGP*nGP), rBdNdw;
  int ig = 1; // gauss point iterator
  for (int zeta = 0; zeta < nGP; zeta++)
    for (int eta = 0; eta < nGP; eta++)
      for (int xi = 0; xi < nGP; xi++, ig++) {
        basis1.computeBasisValues(xg[xi],   u, 1);
        basis2.computeBasisValues(xg[eta],  v, 1);
        basis3.computeBasisValues(xg[zeta], w, 1);
        int ib = 1; // basis function iterator
        for (int k = 0; k < p3; k++)
          for (int j = 0; j < p2; j++)
            for (int i = 0; i < p1; i++, ib++) {
              BN(ib,ig)    = u[2*i  ]*v[2*j  ]*w[2*k  ];
              BdNdu(ib,ig) = u[2*i+1]*v[2*j  ]*w[2*k  ];
              BdNdv(ib,ig) = u[2*i  ]*v[2*j+1]*w[2*k  ];
              BdNdw(ib,ig) = u[2*i  ]*v[2*j  ]*w[2*k+1];
            }
      }

  // Get the reduced integration quadrature points, if needed
  const double* xr = nullptr;
  const double* wr = nullptr;
  int nRed = integrand.getReducedIntegration(nGP);
  if (nRed > 0)
  {
    xr = GaussQuadrature::getCoord(nRed);
    wr = GaussQuadrature::getWeight(nRed);
    if (!xr || !wr) return false;

    rBN.resize(   p1*p2*p3, nRed*nRed*nRed);
    rBdNdu.resize(p1*p2*p3, nRed*nRed*nRed);
    rBdNdv.resize(p1*p2*p3, nRed*nRed*nRed);
    rBdNdw.resize(p1*p2*p3, nRed*nRed*nRed);
    int ig = 1; // gauss point iterator
    for (int zeta = 0; zeta < nRed; zeta++)
      for (int eta = 0; eta < nRed; eta++)
        for (int xi = 0; xi < nRed; xi++, ig++) {
          basis1.computeBasisValues(xr[xi],   u, 1);
          basis2.computeBasisValues(xr[eta],  v, 1);
          basis3.computeBasisValues(xr[zeta], w, 1);
          int ib = 1; // basis function iterator
          for (int k = 0; k < p3; k++)
            for (int j = 0; j < p2; j++)
              for (int i = 0; i < p1; i++, ib++) {
                rBN(ib,ig)    = u[2*i  ]*v[2*j  ]*w[2*k  ];
                rBdNdu(ib,ig) = u[2*i+1]*v[2*j  ]*w[2*k  ];
                rBdNdv(ib,ig) = u[2*i  ]*v[2*j+1]*w[2*k  ];
                rBdNdw(ib,ig) = u[2*i  ]*v[2*j  ]*w[2*k+1];
              }
        }
  }
  else if (nRed < 0)
    nRed = nGP; // The integrand needs to know nGauss

  ThreadGroups oneGroup;
  if (glInt.threadSafe()) oneGroup.oneGroup(nel);
  const IntMat& group = glInt.threadSafe() ? oneGroup[0] : threadGroups[0];


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t t = 0; t < group.size() && ok; t++)
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < group[t].size(); e++)
    {
      if (!ok)
        continue;

      int iel = group[t][e] + 1;
#if defined(SP_DEBUG) && !defined(USE_OPENMP)
      if (dbgElm < 0 && iel != -dbgElm)
        continue; // Skipping all elements, except for -dbgElm
#endif

      FiniteElement fe;
      fe.iel = MLGE[iel-1];
      fe.p   = p1 - 1;
      fe.q   = p2 - 1;
      fe.r   = p3 - 1;
      Matrix   B(p1*p2*p3,4); // Bezier evaluation points and derivatives
      const Matrix& C = bezierExtract[iel-1];
      Matrix   dNdu, Xnod, Jac;
      Matrix3D d2Ndu2, Hess;
      double   dXidu[3];
      double   param[3] = { 0.0, 0.0, 0.0 };
      Vec4     X(param);

      // Get element volume in the parameter space
      const LR::Element* el = lrspline->getElement(iel-1);
      double du = el->umax() - el->umin();
      double dv = el->vmax() - el->vmin();
      double dw = el->wmax() - el->wmin();
      double dV = el->volume();
      if (dV < 0.0)
      {
        ok = false;
        continue;
      }

      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel))
      {
        ok = false;
        continue;
      }

      // Compute parameter values of the Gauss points over this element
      std::array<RealArray,3> gpar, redpar;
      for (int d = 0; d < 3; d++)
      {
        this->getGaussPointParameters(gpar[d],d,nGP,iel,xg);
        if (xr)
          this->getGaussPointParameters(redpar[d],d,nRed,iel,xr);
      }

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
        fe.h = this->getElementCorners(iel,fe.XC);

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
      {
        // Compute the element center
        Go::Point X0;
        double u0 = 0.5*(el->getParmin(0) + el->getParmax(0));
        double v0 = 0.5*(el->getParmin(1) + el->getParmax(1));
        double w0 = 0.5*(el->getParmin(2) + el->getParmax(2));
#pragma omp critical
        lrspline->point(X0,u0,v0,w0);
        X = SplineUtils::toVec3(X0);
      }

      if (integrand.getIntegrandType() & Integrand::G_MATRIX)
      {
        // Element size in parametric space
        dXidu[0] = el->getParmin(0);
        dXidu[1] = el->getParmin(1);
        dXidu[2] = el->getParmin(2);
      }

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(MNPC[iel-1].size(),fe.iel);
      if (!integrand.initElement(MNPC[iel-1],fe,X,nRed*nRed*nRed,*A))
      {
        A->destruct();
        ok = false;
        continue;
      }

      if (xr)
      {
        // --- Selective reduced integration loop ------------------------------

        int ig = 1;
        for (int k = 0; k < nRed; k++)
          for (int j = 0; j < nRed; j++)
            for (int i = 0; i < nRed; i++, ig++)
            {
              // Local element coordinates of current integration point
              fe.xi   = xr[i];
              fe.eta  = xr[j];
              fe.zeta = xr[k];

              // Parameter values of current integration point
              fe.u = param[0] = redpar[0][i];
              fe.v = param[1] = redpar[1][j];
              fe.w = param[2] = redpar[2][k];

              // Extract bezier basis functions
              B.fillColumn(1, rBN.getColumn(ig));
              B.fillColumn(2, rBdNdu.getColumn(ig)*2.0/du);
              B.fillColumn(3, rBdNdv.getColumn(ig)*2.0/dv);
              B.fillColumn(4, rBdNdw.getColumn(ig)*2.0/dw);
              this->evaluateBasis(fe, dNdu, C, B);

              // Compute Jacobian inverse and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
              if (fe.detJxW == 0.0) continue; // skip singular points

              // Cartesian coordinates of current integration point
              X.assign(Xnod * fe.N);
              X.t = time.t;

              // Compute the reduced integration terms of the integrand
              fe.detJxW *= 0.125*dV*wr[i]*wr[j]*wr[k];
              if (!integrand.reducedInt(*A,fe,X))
                ok = false;
            }
      }


      // --- Integration loop over all Gauss points in each direction ----------

      int ig = 1;
      int jp = (iel-1)*nGP*nGP*nGP;
      fe.iGP = firstIp + jp; // Global integration point counter

      for (int k = 0; k < nGP; k++)
        for (int j = 0; j < nGP; j++)
          for (int i = 0; i < nGP; i++, fe.iGP++, ig++)
          {
            // Local element coordinates of current integration point
            fe.xi   = xg[i];
            fe.eta  = xg[j];
            fe.zeta = xg[k];

            // Parameter values of current integration point
            fe.u = param[0] = gpar[0][i];
            fe.v = param[1] = gpar[1][j];
            fe.w = param[2] = gpar[2][k];

            // Fetch basis function derivatives at current integration point
            if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
#pragma omp critical
              this->evaluateBasis(iel, fe, dNdu, d2Ndu2);
            else
            {
              // Extract bezier basis functions
              B.fillColumn(1, BN.getColumn(ig));
              B.fillColumn(2, BdNdu.getColumn(ig)*2.0/du);
              B.fillColumn(3, BdNdv.getColumn(ig)*2.0/dv);
              B.fillColumn(4, BdNdw.getColumn(ig)*2.0/dw);
              this->evaluateBasis(fe, dNdu, C, B);

#ifdef SP_DEBUG
              // Check for errors in the bezier extraction
              if (fabs(fe.N.sum()-1.0) > 1.0e-10) {
                std::cerr <<" *** N does not sum to one at integration point #"<< ig << std::endl;
                exit(123);
              }
              char u = 'u';
              for (size_t d = 1; d <= 3; d++, u++)
                if (fabs(dNdu.getColumn(d).sum()) > 1.0e-10) {
                  std::cerr <<" *** dNd"<< u <<" does not sum to zero at integration point #"<< ig << std::endl;
                  exit(123);
                }
              if (fabs(B.getColumn(1).sum()-1.0) > 1.0e-10) {
                std::cerr <<" *** Bezier basis does not sum to one at integration point #"<< ig << std::endl;
                exit(123);
              }
              if (fabs(B.getColumn(2).sum()) > 1.0e-10) {
                std::cerr <<" *** Bezier derivatives do not sum to zero at integration point #"<< ig << std::endl;
                exit(123);
              }
#endif
            }

            // Compute Jacobian inverse of coordinate mapping and derivatives
            fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
            if (fe.detJxW == 0.0) continue; // skip singular points

            // Compute Hessian of coordinate mapping and 2nd order derivatives
            if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
              if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
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
            X.t = time.t;

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= 0.125*dV*wg[i]*wg[j]*wg[k];
            PROFILE3("Integrand::evalInt");
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

#if defined(SP_DEBUG) && !defined(USE_OPENMP)
      if (iel == -dbgElm)
        break; // Skipping all elements, except for -dbgElm
#endif
    }

  return ok;
}


bool ASMu3D::integrate (Integrand& integrand, int lIndex,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu3D::integrate(B)");

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  // Get Gaussian quadrature points and weights
  // Use the largest polynomial order of the two tangent directions
  const int pmax = std::max(lrspline->order(t1-1),lrspline->order(t2-1));
  const int nG1 = this->getNoGaussPt(pmax,true);
  const int nGP = integrand.getBouIntegrationPoints(nG1);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  LR::parameterEdge edge;
  switch (lIndex%10)
  {
  case 1: edge = LR::WEST;   break;
  case 2: edge = LR::EAST;   break;
  case 3: edge = LR::SOUTH;  break;
  case 4: edge = LR::NORTH;  break;
  case 5: edge = LR::BOTTOM; break;
  case 6: edge = LR::TOP;    break;
  default:edge = LR::NONE;
  }

  // fetch all elements along the chosen edge
  std::vector<LR::Element*> edgeElms;
  lrspline->getEdgeElements(edgeElms,edge);


  // === Assembly loop over all elements on the patch face =====================

  bool ok = true;
  for (LR::Element* el : edgeElms)
  {
    int iEl = el->getId();
#ifdef SP_DEBUG
    if (dbgElm < 0 && iEl+1 != -dbgElm)
      continue; // Skipping all elements, except for -dbgElm
#endif
    if (!myElms.empty() && !glInt.threadSafe() &&
        std::find(myElms.begin(), myElms.end(), iEl) == myElms.end())
     continue;

    FiniteElement fe;
    fe.iel = MLGE[iEl];
    fe.p   = lrspline->order(0) - 1;
    fe.q   = lrspline->order(1) - 1;
    fe.r   = lrspline->order(2) - 1;

    // Compute parameter values of the Gauss points over the whole element
    std::array<Vector,3> gpar;
    for (int d = 0; d < 3; d++)
      if (-1-d == faceDir)
      {
        gpar[d].resize(1);
        gpar[d].fill(lrspline->startparam(d));
      }
      else if (1+d == faceDir)
      {
        gpar[d].resize(1);
        gpar[d].fill(lrspline->endparam(d));
      }
      else
        this->getGaussPointParameters(gpar[d],d,nGP,iEl+1,xg);

    fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
    fe.u  = gpar[0](1);
    fe.v  = gpar[1](1);
    fe.w  = gpar[2](1);

    Matrix dNdu, Xnod, Jac;
    double param[3] = { fe.u, fe.v, fe.w };
    Vec4   X(param);
    Vec3   normal;
    double dXidu[3];

    // Get element face area in the parameter space
    double dA = 0.25*this->getParametricArea(iEl+1,abs(faceDir));
    if (dA < 0.0) // topology error (probably logic error)
    {
      ok = false;
      break;
    }

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iEl+1))
    {
      ok = false;
      break;
    }

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      fe.h = this->getElementCorners(iEl+1,fe.XC);

    if (integrand.getIntegrandType() & Integrand::G_MATRIX)
    {
      // Element size in parametric space
      dXidu[0] = el->getParmax(0) - el->getParmin(0);
      dXidu[1] = el->getParmax(1) - el->getParmin(1);
      dXidu[2] = el->getParmax(2) - el->getParmin(2);
    }

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(MNPC[iEl].size(),fe.iel,true);
    if (!integrand.initElementBou(MNPC[iEl],*A))
    {
      A->destruct();
      ok = false;
      break;
    }

    // --- Integration loop over all Gauss points in each direction ------------

    fe.iGP = firstp; // Global integration point counter
    firstp += nGP*nGP;

    for (int j = 0; j < nGP; j++)
      for (int i = 0; i < nGP; i++, fe.iGP++)
      {
        // Local element coordinates and parameter values
        // of current integration point
        int k1, k2, k3;
        switch (abs(faceDir))
        {
          case 1: k2 = i; k3 = j; k1 = 0; break;
          case 2: k1 = i; k3 = j; k2 = 0; break;
          case 3: k1 = i; k2 = j; k3 = 0; break;
          default: k1 = k2 = k3 = 0;
        }
        if (gpar[0].size() > 1)
        {
          fe.xi = xg[k1];
          fe.u = param[0] = gpar[0](k1+1);
        }
        if (gpar[1].size() > 1)
        {
          fe.eta = xg[k2];
          fe.v = param[1] = gpar[1](k2+1);
        }
        if (gpar[2].size() > 1)
        {
          fe.zeta = xg[k3];
          fe.w = param[2] = gpar[2](k3+1);
        }

        // Fetch basis function derivatives at current integration point
        this->evaluateBasis(iEl, fe, dNdu);

        // Compute basis function derivatives and the face normal
        fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
        if (fe.detJxW == 0.0) continue; // skip singular points

        if (faceDir < 0) normal *= -1.0;

        // Compute G-matrix
        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
          utl::getGmat(Jac,dXidu,fe.G);

#if SP_DEBUG > 4
        if (iEl+1 == dbgElm || iEl+1 == -dbgElm || dbgElm == 0)
          std::cout <<"\n"<< fe;
#endif

        // Cartesian coordinates of current integration point
        X.assign(Xnod * fe.N);
        X.t = time.t;

        // Evaluate the integrand and accumulate element contributions
        fe.detJxW *= dA*wg[i]*wg[j];
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
    if (dbgElm < 0 && iEl+1 == -dbgElm)
      break; // Skipping all elements, except for -dbgElm
#endif
  }

  return ok;
}


bool ASMu3D::integrateEdge (Integrand& integrand, int lEdge,
                            GlobalIntegral& glInt,
                            const TimeDomain& time)
{
  std::cerr << "ASMu3D::integrateEdge(...) is not properly implemented yet :(" << std::endl;
  exit(776654);
  return false;
}


bool ASMu3D::diracPoint (Integrand& integrand, GlobalIntegral& glInt,
                         const double* param, const Vec3& pval)
{
  if (!lrspline) return false;

  int iel = lrspline->getElementContaining(param[0],param[1],param[2]);

  FiniteElement fe;
  fe.iel = MLGE[iel];
  fe.u   = param[0];
  fe.v   = param[1];
  fe.w   = param[2];
  this->evaluateBasis(iel,fe);

  LocalIntegral* A = integrand.getLocalIntegral(MNPC[iel].size(),fe.iel,true);
  bool ok = integrand.evalPoint(*A,fe,pval) && glInt.assemble(A,fe.iel);

  A->destruct();

  return ok;
}


int ASMu3D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!lrspline) return -3;

  int i;
  for (i = 0; i < 3; i++)
    param[i] = (1.0-xi[i])*lrspline->startparam(i) + xi[i]*lrspline->endparam(i);

  Go::Point X0;
  lrspline->point(X0,param[0],param[1],param[2]);
  for (i = 0; i < 3 && i < lrspline->dimension(); i++)
    X[i] = X0[i];

  return 0;
}


int ASMu3D::findElementContaining (const double* param) const
{
  return lrspline ? 1 + lrspline->getElementContaining(param[0],param[1],param[2]) : -2;
}


bool ASMu3D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
  if (!lrspline) return false;

  for (LR::Element* el : lrspline->getAllElements())
  {
    // evaluate element at element corner points
    double umin = el->umin();
    double umax = el->umax();
    double vmin = el->vmin();
    double vmax = el->vmax();
    double wmin = el->wmin();
    double wmax = el->wmax();
    for (int iw = 0; iw <= nSegPerSpan; iw++)
      for (int iv = 0; iv <= nSegPerSpan; iv++)
        for (int iu = 0; iu <= nSegPerSpan; iu++)
          if (dir == 0)
            prm.push_back(umin + (umax-umin)/nSegPerSpan*iu);
          else if (dir == 1)
            prm.push_back(vmin + (vmax-vmin)/nSegPerSpan*iv);
          else
            prm.push_back(wmin + (wmax-wmin)/nSegPerSpan*iw);
  }

  return true;
}


bool ASMu3D::tesselate (ElementBlock& grid, const int* npe) const
{
  if (!lrspline) return false;

  if (npe[0] != npe[1] || npe[0] != npe[2]) {
    std::cerr <<" *** ASMu3D::tesselate does not support different resolution in"
              <<" in u-, v- and w-direction.\n     nviz u = "<< npe[0]
              <<", nviz v = "<< npe[1] <<", nviz w = "<< npe[2] << std::endl;
    return false;
  }

  int nNodesPerElement =  npe[0]   * npe[1]   * npe[2];
  int nSubElPerElement = (npe[0]-1)*(npe[1]-1)*(npe[2]-1);
  int nElements        = lrspline->nElements();

  // output is written once for each element resulting in a lot of unnecessary storage
  // this is preferable to figuring out all element topology information
  grid.unStructResize(nElements * nSubElPerElement,
                      nElements * nNodesPerElement);

  int iel = 0, inod = 0;
  for (const LR::Element* el : lrspline->getAllElements())
  {
    // evaluate element at element corner points
    double umin = el->umin();
    double umax = el->umax();
    double vmin = el->vmin();
    double vmax = el->vmax();
    double wmin = el->wmin();
    double wmax = el->wmax();
    for (int iw = 0; iw < npe[2]; iw++)
      for (int iv = 0; iv < npe[1]; iv++)
        for (int iu = 0; iu < npe[0]; iu++, inod++) {
          double u = umin + (umax-umin)/(npe[0]-1)*iu;
          double v = vmin + (vmax-vmin)/(npe[1]-1)*iv;
          double w = wmin + (wmax-wmin)/(npe[2]-1)*iw;
          Go::Point pt;
          lrspline->point(pt, u,v,w, iel, iu!=npe[0]-1, iv!=npe[1]-1, iw!=npe[2]-1);
          grid.setParams(inod, u, v, w);
          grid.setCoor(inod, SplineUtils::toVec3(pt,nsd));
        }
    ++iel;
  }

  int iStart = iel = inod = 0;
  for (int i = 0; i < nElements; i++, iStart += nNodesPerElement)
    for (int iw = 0; iw < npe[2]-1; iw++)
      for (int iv = 0; iv < npe[1]-1; iv++)
        for (int iu = 0; iu < npe[0]-1; iu++) {
          // enumerate nodes counterclockwise around the hex
          grid.setNode(inod++, iStart + (iw  )*npe[0]*npe[1] + (iv  )*npe[0] + (iu  ) );
          grid.setNode(inod++, iStart + (iw  )*npe[0]*npe[1] + (iv  )*npe[0] + (iu+1) );
          grid.setNode(inod++, iStart + (iw  )*npe[0]*npe[1] + (iv+1)*npe[0] + (iu+1) );
          grid.setNode(inod++, iStart + (iw  )*npe[0]*npe[1] + (iv+1)*npe[0] + (iu  ) );
          grid.setNode(inod++, iStart + (iw+1)*npe[0]*npe[1] + (iv  )*npe[0] + (iu  ) );
          grid.setNode(inod++, iStart + (iw+1)*npe[0]*npe[1] + (iv  )*npe[0] + (iu+1) );
          grid.setNode(inod++, iStart + (iw+1)*npe[0]*npe[1] + (iv+1)*npe[0] + (iu+1) );
          grid.setNode(inod++, iStart + (iw+1)*npe[0]*npe[1] + (iv+1)*npe[0] + (iu  ) );
          grid.setElmId(++iel, i+1);
        }

  return true;
}


bool ASMu3D::evalSolution (Matrix& sField, const Vector& locSol,
                           const int* npe, int nf) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  return this->evalSolution(sField,locSol,gpar.data(),false,0,nf);
}


bool ASMu3D::evalSolution (Matrix& sField, const Vector& locSol,
                           const RealArray* gpar, bool, int deriv, int) const
{
  PROFILE2("ASMu3D::evalSol(P)");

  size_t nComp = locSol.size() / this->getNoNodes();
  if (nComp*this->getNoNodes() != locSol.size())
    return false;

  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size() || nPoints != gpar[2].size())
    return false;

  Vector   ptSol;
  Matrix   dNdu, dNdX, Jac, Xnod, eSol, ptDer;
  Matrix3D d2Ndu2, d2NdX2, Hess, ptDer2;

  Go::BasisPts     spline0;
  Go::BasisDerivs  spline1;
  Go::BasisDerivs2 spline2;
  int lel = -1;

  // Evaluate the primary solution field at each point
  sField.resize(nComp,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch element containing evaluation point.
    // Sadly, points are not always ordered in the same way as the elements.
    int iel = lrspline->getElementContaining(gpar[0][i],gpar[1][i],gpar[2][i]);
    if (iel < 0) {
      std::cerr <<" *** ASMu3D::evalSolution: Element at point ("<< gpar[0][i] <<", "
                << gpar[1][i] <<", "<< gpar[2][i] <<") not found."<< std::endl;
      return false;
    }
    utl::gather(MNPC[iel],nComp,locSol,eSol);

    if (iel != lel && deriv > 0)
    {
      lel = iel; // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel+1))
        return false;
    }

    // Evaluate basis function values/derivatives at current parametric point
    // and multiply with control point values to get the point-wise solution
    switch (deriv) {

    case 0: // Evaluate the solution
      lrspline->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline0,iel);
      eSol.multiply(spline0.basisValues,ptSol);
      sField.fillColumn(1+i,ptSol);
      break;

    case 1: // Evaluate first derivatives of the solution
      lrspline->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline1,iel);
      SplineUtils::extractBasis(spline1,ptSol,dNdu);
      utl::Jacobian(Jac,dNdX,Xnod,dNdu);
      ptDer.multiply(eSol,dNdX);
      sField.fillColumn(1+i,ptDer);
      break;

    case 2: // Evaluate second derivatives of the solution
      lrspline->computeBasis(gpar[0][i],gpar[1][i],gpar[2][i],spline2,iel);
      SplineUtils::extractBasis(spline2,ptSol,dNdu,d2Ndu2);
      utl::Jacobian(Jac,dNdX,Xnod,dNdu);
      utl::Hessian(Hess,d2NdX2,Jac,Xnod,d2Ndu2,dNdu);
      ptDer2.multiply(eSol,d2NdX2);
      sField.fillColumn(1+i,ptDer2);
      break;

    default:
      return false;
    }
  }

  return true;
}


bool ASMu3D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                           const int* npe, char project) const
{
  if (npe && (npe[0] != npe[1] || npe[0] != npe[2] || npe[1] != npe[2]))
  {
    std::cerr <<" *** ASMu3D::evalSolution: LR B-splines require the"
              <<" same number of evaluation points in u-, v- and w-direction."
              << std::endl;
    return false;
  }

  // Project the secondary solution onto the spline basis
  LR::LRSplineVolume* v = nullptr;
  if (project == 'S')
    v = this->scRecovery(integrand);
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
        Go::Point p;
        sField.resize(v->dimension(),gpar[0].size()*gpar[1].size()*gpar[2].size());

        int iel = 0; // evaluation points are always structured in element order
        for (size_t i = 0; i < gpar[0].size(); i++)
        {
          if ((i+1)%npe[0] == 0) iel++;
          v->point(p,gpar[0][i],gpar[1][i],gpar[2][i],iel);
          sField.fillColumn(i+1,p.begin());
        }
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
    sField.resize(v->dimension(),v->nBasisFunctions());
    for (LR::Basisfunction* bf : v->getAllBasisfunctions())
      sField.fillColumn(bf->getId()+1,&(*bf->cp()));
    delete v;
    return true;
  }

  std::cerr <<" *** ASMu3D::evalSolution: Failure!"<< std::endl;
  return false;
}


bool ASMu3D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                           const RealArray* gpar, bool) const
{
  sField.resize(0,0);
  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size() || nPoints != gpar[2].size())
    return false;

  PROFILE2("ASMu3D::evalSol(S)");

  // TODO: investigate the possibility of doing "regular" refinement by
  //       uniform tesselation grid and ignoring LR mesh lines

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;

  Vector   solPt;
  Matrix   dNdu, Jac, Xnod;
  Matrix3D d2Ndu2, Hess;

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);
  const int p3 = lrspline->order(2);
  Matrix B(p1*p2*p3, 4); // Bezier evaluation points and derivatives

  // Evaluate the secondary solution field at each point
  int lel = -1;
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch element containing evaluation point
    // sadly, points are not always ordered in the same way as the elements
    int iel = lrspline->getElementContaining(gpar[0][i],gpar[1][i],gpar[2][i]);
    if (iel < 0) {
      std::cerr <<" *** ASMu3D::evalSolution: Element at point ("
                << gpar[0][i] <<", "<< gpar[1][i] <<", "<< gpar[2][i]
                <<") not found."<< std::endl;
      return false;
    }
    Go::BsplineBasis basis1 = getBezierBasis(p1, lrspline->getElement(iel)->umin(), lrspline->getElement(iel)->umax());
    Go::BsplineBasis basis2 = getBezierBasis(p2, lrspline->getElement(iel)->vmin(), lrspline->getElement(iel)->vmax());
    Go::BsplineBasis basis3 = getBezierBasis(p3, lrspline->getElement(iel)->wmin(), lrspline->getElement(iel)->wmax());

    // Evaluate the basis functions at current parametric point
    FiniteElement fe(lrspline->getElement(iel)->nBasisFunctions());
    fe.u   = gpar[0][i];
    fe.v   = gpar[1][i];
    fe.w   = gpar[2][i];

    double u[2*p1];
    double v[2*p2];
    double w[2*p3];
    basis1.computeBasisValues(gpar[0][i], u, 1);
    basis2.computeBasisValues(gpar[1][i], v, 1);
    basis3.computeBasisValues(gpar[2][i], w, 1);

    size_t ib=1; // basis function iterator
    for(int kk=0; kk<p3; kk++) {
      for(int jj=0; jj<p2; jj++) {
        for(int ii=0; ii<p1; ii++, ib++) {
          B(ib,1) = u[2*ii  ]*v[2*jj  ]*w[2*kk  ];
          B(ib,2) = u[2*ii+1]*v[2*jj  ]*w[2*kk  ];
          B(ib,3) = u[2*ii  ]*v[2*jj+1]*w[2*kk  ];
          B(ib,4) = u[2*ii  ]*v[2*jj  ]*w[2*kk+1];
        }
      }
    }

    if (use2ndDer)
      this->evaluateBasis(iel, fe, dNdu, d2Ndu2);
    else
      this->evaluateBasis(fe, dNdu, bezierExtract[iel], B);

    if (iel != lel)
    {
      lel = iel; // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel+1))
        return false;
    }

    // Compute the Jacobian inverse
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
    if (fe.detJxW == 0.0) continue;

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (use2ndDer)
      if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
        continue;

#if SP_DEBUG > 4
    if (1+iel == dbgElm || dbgElm == 0)
      std::cout <<"\n"<< fe;
#endif

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,Xnod*fe.N,MNPC[iel]))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


IntVec ASMu3D::getFaceNodes (int face, int basis, int orient) const
{
  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  LR::parameterEdge edge;
  switch (face) {
  case 1: edge = LR::WEST; break;
  case 2: edge = LR::EAST; break;
  case 3: edge = LR::SOUTH; break;
  case 4: edge = LR::NORTH; break;
  case 5: edge = LR::BOTTOM; break;
  case 6: edge = LR::TOP; break;
  default: return IntVec();
  }

  std::vector<LR::Basisfunction*> edgeFunctions;
  this->getBasis(basis)->getEdgeFunctions(edgeFunctions, edge);
  if (orient > -1) {
    int dir = (face-1)/2;
    int u = dir == 0;
    int v = 1 + (dir != 2);
    ASMLRSpline::Sort(u, v, orient, edgeFunctions);
  }

  IntVec result(edgeFunctions.size());
  std::transform(edgeFunctions.begin(), edgeFunctions.end(), result.begin(),
                 [ofs](LR::Basisfunction* a) { return a->getId()+ofs; });

  return result;
}


void ASMu3D::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                               int, int orient, bool local) const
{
  if (basis == 0)
    basis = 1;

  if (!this->getBasis(basis)) return; // silently ignore empty patches

  nodes = this->getFaceNodes(lIndex, basis, orient);

#if SP_DEBUG > 1
  std::cout <<"Boundary nodes in patch "<< idx+1 <<" edge "<< lIndex <<":";
  for (int n : nodes) std::cout <<" "<< n;
  std::cout << std::endl;
#endif

  if (!local)
    for (int& node : nodes)
      node = this->getNodeID(node);
}


bool ASMu3D::getOrder (int& p1, int& p2, int& p3) const
{
  p1 = geo->order(0);
  p2 = geo->order(1);
  p3 = geo->order(2);

  return true;
}


int ASMu3D::getCorner(int I, int J, int K, int basis) const
{
  int edge = (I > 0 ? LR::EAST  : LR::WEST  ) |
             (J > 0 ? LR::NORTH : LR::SOUTH ) |
             (K > 0 ? LR::TOP   : LR::BOTTOM);

  const LR::LRSplineVolume* vol = this->getBasis(basis);

  std::vector<LR::Basisfunction*> corner; // vector of one function for corner-input
  vol->getEdgeFunctions(corner, (LR::parameterEdge) edge);

  if ( corner.empty() )
    return -1;

  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  return corner.front()->getId()+ofs;
}


void ASMu3D::generateThreadGroups (const Integrand& integrand, bool silence,
                                   bool ignoreGlobalLM)
{
  LR::generateThreadGroups(threadGroups, this->getBasis(1));
  if (silence || threadGroups[0].size() < 2) return;

  std::cout <<"\nMultiple threads are utilized during element assembly.";
  for (size_t i = 0; i < threadGroups[0].size(); i++)
    std::cout <<"\n Color "<< i+1 <<": "
              << threadGroups[0][i].size() <<" elements";
}


bool ASMu3D::updateDirichlet (const std::map<int,RealFunc*>& func,
                              const std::map<int,VecFunc*>& vfunc, double time,
                              const std::map<int,int>* g2l)
{

  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;

  for (size_t i = 0; i < dirich.size(); i++)
  {
    // figure out function index offset (when using multiple basis)
    size_t ofs = 1;
    for (int j = 1; j < dirich[i].basis; j++)
      ofs += this->getNoNodes(j);

    Real2DMat edgeControlmatrix;
    if ((fit = func.find(dirich[i].code)) != func.end())
      this->faceL2projection(dirich[i], *fit->second, edgeControlmatrix, time);
    else if ((vfit = vfunc.find(dirich[i].code)) != vfunc.end())
      this->faceL2projection(dirich[i], *vfit->second, edgeControlmatrix, time);
    else
    {
      std::cerr <<" *** ASMu3D::updateDirichlet: Code "<< dirich[i].code
                <<" is not associated with any function."<< std::endl;
      return false;
    }
    if (edgeControlmatrix.empty())
    {
      std::cerr <<" *** ASMu3D::updateDirichlet: Projection failure."
                << std::endl;
      return false;
    }

    // Loop over the nodes of this boundary curve
    int j = 0;
    for (int node : dirich[i].MLGN)
    {
      // skip corner nodes, since these are special cased (interpolatory)
      auto it = std::find(dirich[i].corners, dirich[i].corners+4, node+ofs);
      if (it != dirich[i].corners+4) {
        ++j;
        continue;
      }

      for (int dofs = dirich[i].dof; dofs > 0; dofs /= 10)
      {
        int dof = dofs%10;
        // Find the constraint equation for current (node,dof)
        MPC pDOF(MLGN[node+ofs-1],dof);
        MPCIter mit = mpcs.find(&pDOF);
        if (mit == mpcs.end()) continue; // probably a deleted constraint

        // Now update the prescribed value in the constraint equation
        (*mit)->setSlaveCoeff(edgeControlmatrix[dirich[i].dof > 10 ? dof-1 : 0][j]);
#if SP_DEBUG > 1
        std::cout <<"Updated constraint: "<< **mit;
#endif
      }
      ++j;
    }
  }

  // The parent class method takes care of the corner nodes with direct
  // evaluation of the Dirichlet functions; since they are interpolatory
  return this->ASMbase::updateDirichlet(func,vfunc,time,g2l);
}


void ASMu3D::remapErrors (RealArray& errors,
                          const RealArray& origErr, bool elemErrors) const
{
  const LR::LRSplineVolume* basis = this->getBasis(1);

  if (elemErrors) {
    errors = origErr;
    return;
  }

  for (const LR::Element* elm : basis->getAllElements())
    for (const LR::Basisfunction* b : elm->support())
      errors[b->getId()] += origErr[elm->getId()];
}


bool ASMu3D::evaluate (const FunctionBase* func, RealArray& vec,
                       int, double time) const
{
  Matrix ctrlPvals;
  ASMu3D* pch = const_cast<ASMu3D*>(this);
  bool ok = pch->L2projection(ctrlPvals,const_cast<FunctionBase*>(func),time);
  vec = ctrlPvals;
  return ok;
}


bool ASMu3D::transferGaussPtVars (const LR::LRSpline* old_basis,
                                  const RealArray& oldVars, RealArray& newVars,
                                  int nGauss) const
{
  const LR::LRSplineVolume* newBasis = this->getBasis();
  const LR::LRSplineVolume* oldBasis = static_cast<const LR::LRSplineVolume*>(old_basis);

  size_t nGp = nGauss*nGauss*nGauss;
  newVars.resize(newBasis->nElements()*nGp);

  const double* xi = GaussQuadrature::getCoord(nGauss);
  LagrangeInterpolator interp(Vector(xi,nGauss));

  int iEl = 0;
  for (const LR::Element* newEl : newBasis->getAllElements())
  {
    double u_center = 0.5*(newEl->umin() + newEl->umax());
    double v_center = 0.5*(newEl->vmin() + newEl->vmax());
    double w_center = 0.5*(newEl->wmin() + newEl->wmax());
    int iOld = oldBasis->getElementContaining(u_center,v_center,w_center);
    if (iOld < 0)
    {
      std::cerr <<" *** ASMu3D: Failed to locate element "<< newEl->getId()
                <<" of the new mesh in the old mesh."<< std::endl;
      return false;
    }
    const LR::Element* oldEl = oldBasis->getElement(iOld);

    std::array<Matrix,3> I;
    for (int i = 0; i < 3; i++)
    {
      RealArray UGP;
      LR::getGaussPointParameters(newBasis, UGP, i, nGauss, iEl+1, xi);
      double pmin = oldEl->getParmin(i);
      double pmax = oldEl->getParmax(i);
      for (size_t j = 0; j < UGP.size(); j++)
        UGP[j] = -1.0 + 2.0*(UGP[j]-pmin)/(pmax-pmin);

      // lagrangian interpolation
      I[i] = interp.get(UGP);
    }

    Matrix data(nGauss,nGauss*nGauss);
    data.fill(&oldVars[nGp*iOld],nGp);

    Matrix I1 = I[0]*data;
    for (int i = 0; i < nGauss; ++i)
      for (int j = 0; j < nGauss; ++j)
        for (int k = 0; k < nGauss; ++k)
          data(j + 1, i + k*(nGauss) + 1) = I1(i+1, 1 + j + k*(nGauss));
    I1 = I[1]*data;
    for (int i = 0; i < nGauss; ++i)
      for (int j = 0; j < nGauss; ++j)
        for (int k = 0; k < nGauss; ++k)
          data(k + 1, i + j*(nGauss) + 1) = I1(j+1, 1 + i + k*(nGauss));

    I1 = I[2]*data;
    for (int i = 0; i < nGauss; ++i)
      for (int j = 0; j < nGauss; ++j)
        for (int k = 0; k < nGauss; ++k)
          newVars[iEl*nGp+i+(j + k*nGauss)*nGauss] = I1(k+1, 1+i + j*nGauss);

    ++iEl;
  }

  return true;
}


bool ASMu3D::transferGaussPtVarsN (const LR::LRSpline* old_basis,
                                   const RealArray& oldVars, RealArray& newVars,
                                   int nGauss) const
{
  const LR::LRSplineVolume* newBasis = this->getBasis();
  const LR::LRSplineVolume* oldBasis = static_cast<const LR::LRSplineVolume*>(old_basis);

  size_t nGP = nGauss*nGauss*nGauss;
  newVars.clear();
  newVars.reserve(nGP*newBasis->nElements());

  struct Param{ double u,v,w; };
  std::vector<Param> oGP(nGP);

  const double* xi = GaussQuadrature::getCoord(nGauss);
  for (const LR::Element* newEl : newBasis->getAllElements())
  {
    double u_center = 0.5*(newEl->umin() + newEl->umax());
    double v_center = 0.5*(newEl->vmin() + newEl->vmax());
    double w_center = 0.5*(newEl->wmin() + newEl->wmax());
    int iOld = oldBasis->getElementContaining(u_center,v_center,w_center);
    if (iOld < 0)
    {
      std::cerr <<" *** ASMu3D: Failed to locate element "<< newEl->getId()
                <<" of the new mesh in the old mesh."<< std::endl;
      return false;
    }
    const LR::Element* oldEl = oldBasis->getElement(iOld);

    // find parameters of old gauss points
    double umin = oldEl->umin();
    double vmin = oldEl->vmin();
    double wmin = oldEl->wmin();
    double du = 0.5*(oldEl->umax() - umin);
    double dv = 0.5*(oldEl->vmax() - vmin);
    double dw = 0.5*(oldEl->wmax() - wmin);
    size_t l = 0;
    for (int k = 0; k < nGauss; ++k)
      for (int j = 0; j < nGauss; ++j)
        for (int i = 0; i < nGauss; ++i, ++l) {
          oGP[l].u = umin + du * (xi[i] + 1.0);
          oGP[l].v = vmin + dv * (xi[j] + 1.0);
          oGP[l].w = wmin + dw * (xi[k] + 1.0);
        }

    // parameters of new gauss points
    umin = newEl->umin();
    vmin = newEl->vmin();
    wmin = newEl->wmin();
    du = 0.5*(newEl->umax() - umin);
    dv = 0.5*(newEl->vmax() - vmin);
    dw = 0.5*(newEl->wmax() - wmin);
    for (int k = 0; k < nGauss; ++k)
      for (int j = 0; j < nGauss; ++j)
        for (int i = 0; i < nGauss; ++i) {
          double u = umin + du * (xi[i] + 1.0);
          double v = vmin + dv * (xi[j] + 1.0);
          double w = wmin + dw * (xi[k] + 1.0);
          double dist = 1.0e16;
          size_t near = 0;
          for (size_t l = 0; l < nGP; ++l) {
            Vec3 d(oGP[l].u-u, oGP[l].v-v, oGP[l].w-w);
            double nd = d.length();
            if (nd < dist) {
              near = l;
              dist = nd;
            }
          }
          newVars.push_back(oldVars[iOld*nGP+near]);
      }
  }

  return true;
}


bool ASMu3D::transferCntrlPtVars (const LR::LRSpline* old_basis,
                                  RealArray& newVars, int nGauss) const
{
  const LR::LRSplineVolume* newBasis = this->getBasis();
  const LR::LRSplineVolume* oldBasis = static_cast<const LR::LRSplineVolume*>(old_basis);

  newVars.clear();
  newVars.reserve(newBasis->nElements()*nGauss*nGauss*nGauss*oldBasis->dimension());
  const double* xi = GaussQuadrature::getCoord(nGauss);

  for (int iel = 0; iel < newBasis->nElements(); iel++)
  {
    RealArray U, V, W, ptVar;
    LR::getGaussPointParameters(newBasis, U, 0, nGauss, iel, xi);
    LR::getGaussPointParameters(newBasis, V, 1, nGauss, iel, xi);
    LR::getGaussPointParameters(newBasis, W, 2, nGauss, iel, xi);
    for (int k = 0; k < nGauss; k++)
      for (int j = 0; j < nGauss; j++)
        for (int i = 0; i < nGauss; i++)
        {
          oldBasis->point(ptVar,U[i],V[j],W[k]);
          for (size_t l = 0; l < ptVar.size(); l++)
            newVars.push_back(ptVar[l]);
        }
  }

  return true;
}


/*!
  Refines all elements for which refC(X0) < refTol,
  where X0 is the element center.
*/

bool ASMu3D::refine (const RealFunc& refC, double refTol)
{
  Go::Point X0;
  int iel = 0;
  IntVec elements;
  for (const LR::Element* elm : lrspline->getAllElements())
  {
    double u0 = 0.5*(elm->umin() + elm->umax());
    double v0 = 0.5*(elm->vmin() + elm->vmax());
    double w0 = 0.5*(elm->wmin() + elm->wmax());
    lrspline->point(X0,u0,v0,w0);
    if (refC(SplineUtils::toVec3(X0,nsd)) < refTol)
      elements.push_back(iel);
    ++iel;
  }

  Vectors dummySol;
  LR::RefineData prm(true);
  prm.options = { 10, 1, 2 };
  prm.elements = this->getFunctionsForElements(elements);
  return this->refine(prm,dummySol);
}


double ASMu3D::getMinimumSize (int nrefinements) const
{
  if (vMin > 0.0 || lrspline->nElements() <= 0)
    return vMin;

  double redMax = pow(2.0,nrefinements);
  vMin = lrspline->getElement(0)->volume()/(redMax*redMax*redMax);
  return vMin;
}


bool ASMu3D::checkElementSize (int elmId, bool globalNum) const
{
  if (globalNum)
  {
    IntVec::const_iterator it = std::find(MLGE.begin(),MLGE.end(),1+elmId);
    if (it == MLGE.end()) return false;

    elmId = it - MLGE.begin();
  }

  if (elmId < lrspline->nElements())
    return lrspline->getElement(elmId)->volume() > vMin+1.0e-12;
  else
    return false;
}


void ASMu3D::extendRefinementDomain (IntSet& refineIndices,
                                     const IntSet& neighborIndices) const
{
  const int nedge = 12;
  const int nface =  6;

  IntVec bndry0;
  for (int K = -1; K < 2; K += 2)
    for (int J = -1; J < 2; J += 2)
      for (int I = -1; I < 2; I += 2)
        bndry0.push_back(this->getCorner(I,J,K,1));

  std::vector<IntVec> bndry1;
  for (int j = 1; j <= nedge; j++)
    bndry1.push_back(this->getEdge(j, true, 1, 0));

  std::vector<IntVec> bndry2(nface);
  for (int j = 1; j <= nface; j++)
    this->getBoundaryNodes(j, bndry2[j-1], 1, 1, 0, true);

  // Add refinement from neighbors
  for (int j : neighborIndices)
  {
    bool done_with_this_node = false;

    // Check if node is a corner node,
    // compute large extended domain (all directions)
    for (int edgeNode : bndry0)
      if (edgeNode-1 == j)
      {
        IntVec secondary = this->getOverlappingNodes(j);
        refineIndices.insert(secondary.begin(),secondary.end());
        done_with_this_node = true;
        break;
      }

    // Check if node is an edge node,
    // compute moderate extended domain (2 directions)
    int allowedDir;
    for (int edge = 0; edge < nedge && !done_with_this_node; edge++)
      for (int edgeNode : bndry1[edge])
        if (edgeNode-1 == j)
        {
          if (edge < 4)
            allowedDir = 6; // bin(110), allowed to grow in v- and w-direction
          else if (edge < 8)
            allowedDir = 5; // bin(101), allowed to grow in u- and w-direction
          else
            allowedDir = 3; // bin(011), allowed to grow in u- and v-direction
          IntVec secondary = this->getOverlappingNodes(j,allowedDir);
          refineIndices.insert(secondary.begin(),secondary.end());
          done_with_this_node = true;
          break;
        }

    // Check if node is a face node,
    // compute small extended domain (1 direction)
    for (int face = 0; face < nface && !done_with_this_node; face++)
      for (int edgeNode : bndry2[face])
        if (edgeNode-1 == j)
        {
          allowedDir = 1 << face/2;
          IntVec secondary = this->getOverlappingNodes(j,allowedDir);
          refineIndices.insert(secondary.begin(),secondary.end());
          done_with_this_node = true;
          break;
        }
  }
}


void ASMu3D::getElmConnectivities (IntMat& neigh) const
{
  const LR::LRSplineVolume* lr = this->getBasis(1);
  for (const LR::Element* m : lr->getAllElements()) {
    int gEl = MLGE[m->getId()]-1;
    for (auto edge : {LR::WEST, LR::EAST, LR::SOUTH, LR::NORTH, LR::BOTTOM, LR::TOP}) {
      std::set<int> elms = lr->getElementNeighbours(m->getId(), edge);
      for (int elm : elms)
        neigh[gEl].push_back(MLGE[elm]-1);
    }
  }
}


void ASMu3D::getBoundaryElms (int lIndex, int orient, IntVec& elms) const
{
  LR::parameterEdge edge;
  switch (lIndex)
  {
  case 1: edge = LR::WEST;   break;
  case 2: edge = LR::EAST;   break;
  case 3: edge = LR::SOUTH;  break;
  case 4: edge = LR::NORTH;  break;
  case 5: edge = LR::BOTTOM; break;
  case 6: edge = LR::TOP;    break;
  default:edge = LR::NONE;
  }

  std::vector<LR::Element*> elements;
  this->getBasis(1)->getEdgeElements(elements, edge);

  std::sort(elements.begin(), elements.end(),
            [lIndex,orient](const LR::Element* a, const LR::Element* b)
            {
              int dir = (lIndex - 1) / 2;
              int u = dir == 0 ? 1 : 0;
              int v = 1 + (dir != 2 ? 1 : 0);
              int idx = (orient & 4) ? v : u;
              auto A = a->midpoint();
              auto B = b->midpoint();
              if (A[idx] != B[idx])
                return (orient & 2) ? A[idx] > B[idx] : A[idx] < B[idx];

              idx = (orient & 4) ? u : v;
              if (A[idx] != B[idx])
                return (orient & 1) ? A[idx] > B[idx] : A[idx] < B[idx];

              return false;
            });


  for (const LR::Element* elem : elements)
    elms.push_back(MLGE[elem->getId()]-1);
}


void ASMu3D::generateThreadGroupsFromElms(const std::vector<int>& elms)
{
  myElms.clear();
  for (int elm : elms)
    if (this->getElmIndex(elm+1) > 0)
      myElms.push_back(this->getElmIndex(elm+1)-1);

  threadGroups = threadGroups.filter(myElms);
}
