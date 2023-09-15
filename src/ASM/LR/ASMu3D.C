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
#include "LRSplineField3D.h"
#include "LRSplineFields3D.h"
#include "ElementBlock.h"
#include "MPC.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "Point.h"
#include "IFEM.h"
#include <array>
#include <utility>


ASMu3D::ASMu3D (unsigned char n_f)
  : ASMLRSpline(3,3,n_f), lrspline(nullptr),
    bezierExtract(myBezierExtract)
{
  vMin = 0.0;
  tensorspline = tensorPrjBas = nullptr;
}


ASMu3D::ASMu3D (const ASMu3D& patch, unsigned char n_f)
  : ASMLRSpline(patch,n_f), lrspline(patch.lrspline),
    bezierExtract(patch.myBezierExtract)
{
  vMin = 0.0;
  tensorspline = tensorPrjBas = nullptr;

  // Need to set nnod here,
  // as hasXNodes might be invoked before the FE data is generated
  if (nnod == 0 && lrspline)
    nnod = lrspline->nBasisFunctions();
}


const LR::LRSplineVolume* ASMu3D::getBasis (int basis) const
{
  switch (basis) {
    case ASM::GEOMETRY_BASIS:
      return static_cast<const LR::LRSplineVolume*>(geomB.get());
    case ASM::PROJECTION_BASIS:
      return static_cast<const LR::LRSplineVolume*>(projB.get());
    case ASM::PROJECTION_BASIS_2:
      return static_cast<const LR::LRSplineVolume*>(projB2.get());
    case ASM::REFINEMENT_BASIS:
      return static_cast<const LR::LRSplineVolume*>(refB.get());
    default:
      return lrspline.get();
  }
}


LR::LRSplineVolume* ASMu3D::getBasis (int basis)
{
  if (tensorspline)
    this->createLRfromTensor();

  return const_cast<LR::LRSplineVolume*>(std::as_const(*this).getBasis(basis));
}


bool ASMu3D::read (std::istream& is)
{
  if (shareFE) return true;

  // Read the input file as either an LRSpline file directly,
  // or a tensor product B-spline and convert
  char firstline[256];
  is.getline(firstline, 256);
  if (strncmp(firstline, "# LRSPLINE", 10) == 0) {
    lrspline.reset(new LR::LRSplineVolume());
    is >> *lrspline;
  }
  else {
    // Probably a SplineVolume, so we'll read that and convert
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

  geomB = lrspline;
  return true;
}


bool ASMu3D::write (std::ostream& os, int basis) const
{
  if (!lrspline) return false;
  if (basis > static_cast<int>(this->getNoBasis())) return false;
  const LR::LRSplineVolume* spline = this->getBasis(basis);
  if (!spline) return false;

  os << *spline;

  return os.good();
}


void ASMu3D::clear (bool retainGeometry)
{
  if (!retainGeometry) {
    // Erase spline data
    if (!shareFE) {
      lrspline.reset();
      delete tensorspline;
      delete tensorPrjBas;
    }
    geomB = nullptr;
    tensorspline = tensorPrjBas = nullptr;
  }

  // Erase the FE data
  this->ASMbase::clear(retainGeometry);
  this->dirich.clear();
  projThreadGroups = ThreadGroups();

  myCache.clear();
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
  return true;
}


bool ASMu3D::raiseOrder (int ru, int rv, int rw, bool setOrder)
{
  if (!tensorspline) return false;
  if (shareFE) return true;

  if (setOrder)
  {
    ru -= tensorspline->order(0);
    rv -= tensorspline->order(1);
    rw -= tensorspline->order(2);
  }
  tensorspline->raiseOrder(ru,rv,rw);
  return true;
}


/*!
  This method is supposed to be invoked twice during the model generation.
  In the first call, with \a init = \e true, the spline volume object
  is cloned and the two pointers are then swapped, such that the subsequent
  refine and raiseOrder operations will apply to the projection basis
  and not on the geometry basis.
  In the second call, the pointers are swapped back.

  The method can also be invoked twice with \a init = \e false in case the
  projection basis is to be read from a file.
*/

bool ASMu3D::createProjectionBasis (bool init)
{
  if (!tensorspline)
    return false;
  else if (init && !tensorPrjBas)
    tensorPrjBas = tensorspline->clone();

  std::swap(tensorspline,tensorPrjBas);
  std::swap(geomB,projB);
  lrspline = std::static_pointer_cast<LR::LRSplineVolume>(geomB);
  return true;
}


std::shared_ptr<LR::LRSplineVolume> ASMu3D::createLRfromTensor ()
{
  if (tensorspline)
  {
    lrspline.reset(new LR::LRSplineVolume(tensorspline));
    delete tensorspline;
    tensorspline = nullptr;
  }

  return lrspline;
}


bool ASMu3D::generateFEMTopology ()
{
  refB = geomB = this->createLRfromTensor();

  if (tensorPrjBas)
  {
    projB.reset(new LR::LRSplineVolume(tensorPrjBas));
    delete tensorPrjBas;
    tensorPrjBas = nullptr;
  }
  else if (!projB)
    projB = lrspline;

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
  // force cache creation
  lrspline->getElementContaining(lrspline->getElement(0)->midpoint());

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
                           bool coordCheck, int thick)
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
  IntVec slaveNodes;
  this->getBoundaryNodes(face, slaveNodes, basis, thick, 0, true);
  for (int& it : slaveNodes)
    it += slave;

  // Set up the master node numbers for the neighboring volume patch
  IntVec masterNodes;
  neighbor.getBoundaryNodes(nface, masterNodes, basis, thick, norient, true);
  for (int& it : masterNodes)
    it += master;

  if (masterNodes.empty() || masterNodes.size() != slaveNodes.size())
  {
    std::cerr <<" *** ASMu3D::connectBasis: Non-matching faces, sizes "
              << masterNodes.size() <<" and "<< slaveNodes.size() << std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  for (size_t i = 0; i < masterNodes.size(); i++)
  {
    int node = masterNodes[i];
    int slvn = slaveNodes[i];
    if (!coordCheck)
      ASMbase::collapseNodes(neighbor,node,*this,slvn);
    else if (neighbor.getCoord(node).equal(this->getCoord(slvn),xtol))
      ASMbase::collapseNodes(neighbor,node,*this,slvn);
    else
    {
      std::cerr <<" *** ASMu3D::connectBasis: Non-matching nodes "
                << node <<": "<< neighbor.getCoord(node)
                <<"\n                                          and "
                << slvn <<": "<< this->getCoord(slvn) << std::endl;
      return false;
    }
  }

  return true;
}


ASMu3D::DirichletFace::DirichletFace (LR::LRSplineVolume* sv,
                                      int dir, int d, int c, int offset)
  : lr(sv), edg(LR::NONE), dof(d), code(c), corners{0,0,0,0}
{
  // Figure out what face we are at
  switch (dir) {
  case -1: edg = LR::WEST;   break;
  case  1: edg = LR::EAST;   break;
  case -2: edg = LR::SOUTH;  break;
  case  2: edg = LR::NORTH;  break;
  case -3: edg = LR::BOTTOM; break;
  case  3: edg = LR::TOP;    break;
  default: return;
  }

  // Find the corners since these are not to be included in the L2-fitting
  // of the inhomogenuous dirichlet boundaries; corners are interpolatory.
  // Optimization note: loop over the "edge"-container to manually pick up
  // the end nodes. LRspline::getEdgeFunctions() does a global search.
  typedef std::array<int,3> Vrtx;
  typedef std::vector<Vrtx> FaceVrtx;
  static const std::map<LR::parameterEdge,FaceVrtx> faces = {{
      {LR::WEST,  {{{{-1, -1, -1}},
                    {{-1,  1, -1}},
                    {{-1, -1,  1}},
                    {{-1,  1,  1}}}}},
      {LR::EAST,  {{{{ 1, -1, -1}},
                    {{ 1,  1, -1}},
                    {{ 1, -1,  1}},
                    {{ 1,  1,  1}}}}},
      {LR::SOUTH, {{{{-1, -1, -1}},
                    {{ 1, -1, -1}},
                    {{-1, -1,  1}},
                    {{ 1, -1,  1}}}}},
      {LR::NORTH, {{{{-1,  1, -1}},
                    {{ 1,  1, -1}},
                    {{-1,  1,  1}},
                    {{ 1,  1,  1}}}}},
      {LR::BOTTOM,{{{{-1, -1, -1}},
                    {{ 1, -1, -1}},
                    {{-1,  1, -1}},
                    {{ 1,  1, -1}}}}},
      {LR::TOP,   {{{{-1, -1,  1}},
                    {{ 1, -1,  1}},
                    {{-1,  1,  1}},
                    {{ 1,  1,  1}}}}}
    }};

  int i = 0;
  for (const Vrtx& vx : faces.find(edg)->second)
  {
    std::vector<LR::Basisfunction*> corner;

    lr->getEdgeFunctions(corner,
                         (vx[0] > 0 ? LR::EAST  : LR::WEST  ) |
                         (vx[1] > 0 ? LR::NORTH : LR::SOUTH ) |
                         (vx[2] > 0 ? LR::TOP   : LR::BOTTOM));

    corners[i++] = corner.empty() ? 0 : offset + corner.front()->getId();
  }
}


bool ASMu3D::DirichletFace::isCorner (int b) const
{
  return std::find(corners,corners+4,b) != corners+4;
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

  // Figure out function index offset (when using multiple basis)
  int offset = 1;
  for (int i = 1; i < basis; i++)
    offset += this->getNoNodes(i);

  // Figure out what face we are at
  DirichletFace de(this->getBasis(basis), dir, dof, code, offset);

  // Get all basis functions on this face
  std::vector<LR::Basisfunction*> faceFunctions;
  de.lr->getEdgeFunctions(faceFunctions,de.edg);

  // Add constraints for all basis functions on the face
  for (LR::Basisfunction* b : faceFunctions)
    if (!de.isCorner(b->getId()+offset))
      this->prescribe(b->getId()+offset, dof, -code);
    else if (!open) // skip corners for open boundaries
      // corners are always interpolated (positive 'code')
      this->prescribe(b->getId()+offset, dof, abs(code));

  if (code <= 0) return; // If no projection, we're done

  // Build up the local face-node correspondence for this face
  for (LR::Basisfunction* b : faceFunctions)
    de.MLGN.push_back(b->getId()+offset);

  // Get all elements connected to this face
  std::vector<LR::Element*> faceElements;
  de.lr->getEdgeElements(faceElements,de.edg);

  // Build the MLGE and MNPC arrays
  de.MLGE.reserve(faceElements.size());
  de.MNPC.reserve(faceElements.size());
  for (LR::Element* el : faceElements)
  {
    // for mixed FEM models, let MLGE point to the integration basis
    if (de.lr != this->lrspline.get())
      de.MLGE.push_back(lrspline->getElementContaining(el->midpoint()));
    else
      de.MLGE.push_back(el->getId());

    IntVec mnpc; mnpc.reserve(el->support().size());
    for (LR::Basisfunction* b : el->support())
      mnpc.push_back(utl::findIndex(de.MLGN,b->getId()+offset));
    de.MNPC.push_back(mnpc);
  }

  dirich.push_back(de);
}


size_t ASMu3D::constrainFaceLocal (int dir, bool open, int dof, int code,
                                   bool project, char T1)
{
  return 0; // TODO...
}


void ASMu3D::getBoundary1Nodes (int lEdge, IntVec& nodes,
                                int basis, int orient, bool local, bool) const
{
  if (basis < 1) basis = 1;

  // lEdge = 1-4, running index is u (vmin,wmin), (vmax,wmin), (vmin,wmax), (vmax,wmax)
  // lEdge = 5-8, running index is v (umin,wmin), (umax,wmin), (umin,wmax), (umax,wmax)
  // lEdge = 9-12, running index is w

  auto&& getEdge = [](int in)
  {
    switch (in) {
      case  1: return LR::BOTTOM | LR::SOUTH;
      case  2: return LR::BOTTOM | LR::NORTH;
      case  3: return LR::TOP    | LR::SOUTH;
      case  4: return LR::TOP    | LR::NORTH;
      case  5: return LR::BOTTOM | LR::WEST;
      case  6: return LR::BOTTOM | LR::EAST;
      case  7: return LR::TOP    | LR::WEST;
      case  8: return LR::TOP    | LR::EAST;
      case  9: return LR::SOUTH  | LR::WEST;
      case 10: return LR::SOUTH  | LR::EAST;
      case 11: return LR::NORTH  | LR::WEST;
      case 12: return LR::NORTH  | LR::EAST;
      default: return LR::NONE;
    }
  };

  const LR::parameterEdge edge = getEdge(lEdge);

  // figure out function index offset (when using multiple basis)
  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  // get all the boundary functions from the LRspline object
  std::vector<LR::Basisfunction*> thisEdge;
  this->getBasis(basis)->getEdgeFunctions(thisEdge, edge, 1);
  if (orient >= 0) {
    int u = lEdge <= 4 ? 1 : 0;
    int v = lEdge >= 9 ? 1 : 2;
    ASMLRSpline::Sort(u, v, orient, thisEdge);
  }

  for (LR::Basisfunction* b : thisEdge)
    nodes.push_back(local ? ofs+b->getId() : this->getNodeID(ofs+b->getId()));
}


void ASMu3D::constrainEdge (int lEdge, bool open, int dof, int code, char basis)
{
  if (open)
    std::cout <<"  ** ASMu3D::constrainEdge: Open boundary conditions are not"
              <<" supported for LR B-splines. Treated as closed."<< std::endl;

  IntVec edgeNodes;
  this->getBoundary1Nodes(lEdge,edgeNodes,basis,-1,true);
  for (int node : edgeNodes)
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


bool ASMu3D::getElementCoordinates (Matrix& X, int iel, bool forceItg) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > MNPC.size())
  {
    std::cerr <<" *** ASMu3D::getElementCoordinates: Element index "<< iel
              <<" out of range [1,"<< MNPC.size() <<"]."<< std::endl;
    return false;
  }
#endif
  const LR::LRSplineVolume* spline = this->getBasis(forceItg ? ASM::INTEGRATION_BASIS
                                                             : ASM::GEOMETRY_BASIS);
  if (spline != lrspline.get())
    iel = spline->getElementContaining(lrspline->getElement(iel-1)->midpoint()) + 1;

  const LR::Element* el = spline->getElement(iel-1);
  X.resize(3,el->nBasisFunctions());

  int n = 1;
  for (LR::Basisfunction* b : el->support())
    X.fillColumn(n++,&(*b->cp()));

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


void ASMu3D::getNodalCoordinates (Matrix& X, bool) const
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


/*!
  \brief Static helper mapping a face index into a \a parameterEdge enum value.
*/

static LR::parameterEdge getFaceEnum (int faceIndex)
{
  switch (faceIndex) {
  case 1: return LR::WEST;
  case 2: return LR::EAST;
  case 3: return LR::SOUTH;
  case 4: return LR::NORTH;
  case 5: return LR::BOTTOM;
  case 6: return LR::TOP;
  }

  return LR::NONE;
}


size_t ASMu3D::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (!lrspline || ldim < 2)
    return 0;

  std::vector<LR::Element*> edgeElms;
  if (lIndex > 0 && lIndex <= 6)
    lrspline->getEdgeElements(edgeElms,getFaceEnum(lIndex));

  return edgeElms.size();
}


void ASMu3D::getGaussPointParameters (RealArray& uGP, int dir, int nGauss,
                                      int iEl, const double* xi,
                                      const LR::LRSplineVolume* spline) const
{
  if (!spline)
    spline = lrspline.get();

  const LR::Element* el = spline->getElement(iEl-1);
  double start = el->getParmin(dir);
  double stop  = el->getParmax(dir);

  uGP.resize(nGauss);
  for (int i = 0; i < nGauss; i++)
    uGP[i] = 0.5*((stop-start)*xi[i] + stop+start);
}


double ASMu3D::getElementCorners (int iEl, Vec3Vec& XC, RealArray* uC) const
{
  const LR::Element* el = lrspline->getElement(iEl-1);
  double u[2] = { el->getParmin(0), el->getParmax(0) };
  double v[2] = { el->getParmin(1), el->getParmax(1) };
  double w[2] = { el->getParmin(2), el->getParmax(2) };

  XC.clear();
  XC.reserve(8);
  if (uC) uC->reserve(24);
  Go::Point pt;

  for (int k = 0; k < 2; k++)
    for (int j = 0; j < 2; j++)
      for (int i = 0; i < 2; i++)
      {
        lrspline->point(pt,u[i],v[j],w[k],iEl-1);
        XC.push_back(SplineUtils::toVec3(pt));
        if (uC)
        {
          uC->push_back(u[i]);
          uC->push_back(v[j]);
          uC->push_back(w[k]);
        }
      }

  return getElementSize(XC);
}


void ASMu3D::getCornerPoints (int iel, PointVec& XC) const
{
  RealArray uC;
  Vec3Vec  XYZ;
  this->getElementCorners(iel,XYZ,&uC);

  XC.clear();
  XC.reserve(8);
  for (int i = 0; i < 8; i++)
    XC.push_back(utl::Point(XYZ[i], { uC[3*i], uC[3*i+1], uC[3*i+2] }));
}


void ASMu3D::evaluateBasis (int iel, double u, double v, double w,
                            Vector& N, Matrix& dNdu, int basis) const
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
  this->evaluateBasis(iel, fe.u, fe.v, fe.w, fe.basis(basis), dNdu, basis);
}

void ASMu3D::evaluateBasis (Vector& N, Matrix& dNdu,
                            const Matrix& C, const Matrix& B) const
{
  PROFILE3("ASMu3D::evalBasis(BE)");

  Matrix CB = C*B;
  dNdu.resize(CB.rows(),3);
  N = CB.getColumn(1);
  dNdu.fillColumn(1,CB.getColumn(2));
  dNdu.fillColumn(2,CB.getColumn(3));
  dNdu.fillColumn(3,CB.getColumn(4));
}

void ASMu3D::evaluateBasis (int iel, double u, double v, double w,
                            Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2,
                            int basis) const
{
  PROFILE3("ASMu3D::evalBasis(2)");

  std::vector<RealArray> result;
  this->getBasis(basis)->computeBasis(u, v, w, result, 2, iel);
  size_t n = 0, nBasis = result.size();

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


bool ASMu3D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu3D::integrate(I)");


  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;

  if (myCache.empty())
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this, cachePolicy, 1));

  BasisFunctionCache& cache = *myCache.front();
  cache.setIntegrand(&integrand);
  if (!cache.init(use2ndDer ? 2 : 1))
    return false;

  const std::array<int,3>& ng = cache.nGauss();
  const std::array<const double*,3>& xg = cache.coord();
  const std::array<const double*,3>& wg = cache.weight();

  // Get the reduced integration quadrature points, if needed
  const double* xr = cache.coord(true)[0];
  const double* wr = cache.weight(true)[0];

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);
  const int p3 = lrspline->order(2);

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
      Matrix   Xnod, Jac;
      Matrix3D Hess;
      double   dXidu[3];
      double   param[3] = { 0.0, 0.0, 0.0 };
      Vec4     X(param,time.t);

      // Get element volume in the parameter space
      const LR::Element* el = lrspline->getElement(iel-1);
      double dV = 0.125*el->volume();

      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel))
      {
        ok = false;
        continue;
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
        lrspline->point(X0,u0,v0,w0,iel-1);
        X.assign(SplineUtils::toVec3(X0));
      }

      if (integrand.getIntegrandType() & Integrand::G_MATRIX)
        // Element size in parametric space
        for (int i = 0; i < 3; i++)
          dXidu[i] = el->getParmax(i) - el->getParmin(i);

      size_t nen = el->support().size();
      if (integrand.getIntegrandType() & Integrand::AVERAGE)
      {
        // --- Compute average value of basis functions over the element -----

        int ip = (iel-1)*ng[0]*ng[1]*ng[2] + firstIp;
        fe.Navg.resize(nen,true);
        double vol = 0.0;
        int ig = 0;
        for (int k = 0; k < ng[2]; k++)
          for (int j = 0; j < ng[1]; j++)
            for (int i = 0; i < ng[0]; i++, ++ip, ++ig)
            {
              const BasisFunctionVals& bfs = cache.getVals(iel-1,ig);

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

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(nen,fe.iel);
      int nRed = cache.nGauss(true)[0];
      if (!integrand.initElement(MNPC[iel-1],fe,X,nRed*nRed*nRed,*A))
      {
        A->destruct();
        ok = false;
        continue;
      }

      if (xr)
      {
        // --- Selective reduced integration loop ------------------------------

        int ig = 0;
        for (int k = 0; k < nRed; k++)
          for (int j = 0; j < nRed; j++)
            for (int i = 0; i < nRed; i++, ig++)
            {
              // Local element coordinates of current integration point
              fe.xi   = xr[i];
              fe.eta  = xr[j];
              fe.zeta = xr[k];

              // Parameter values of current integration point
              fe.u = param[0] = cache.getParam(0,iel-1,i,true);
              fe.v = param[1] = cache.getParam(1,iel-1,j,true);
              fe.w = param[2] = cache.getParam(2,iel-1,k,true);

              const BasisFunctionVals& bfs = cache.getVals(iel-1,ig,true);
              fe.N = bfs.N;

              // Compute Jacobian inverse and derivatives
              fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,bfs.dNdu);
              if (fe.detJxW == 0.0) continue; // skip singular points

              // Cartesian coordinates of current integration point
              X.assign(Xnod * fe.N);

              // Compute the reduced integration terms of the integrand
              fe.detJxW *= dV*wr[i]*wr[j]*wr[k];
              if (!integrand.reducedInt(*A,fe,X))
                ok = false;
            }
      }


      // --- Integration loop over all Gauss points in each direction ----------

      int ig = 0;
      int jp = (iel-1)*ng[0]*ng[1]*ng[2];
      fe.iGP = firstIp + jp; // Global integration point counter

      for (int k = 0; k < ng[2]; k++)
        for (int j = 0; j < ng[1]; j++)
          for (int i = 0; i < ng[0]; i++, fe.iGP++, ig++)
          {
            // Local element coordinates of current integration point
            fe.xi   = xg[0][i];
            fe.eta  = xg[1][j];
            fe.zeta = xg[2][k];

            // Parameter values of current integration point
            fe.u = param[0] = cache.getParam(0,iel-1,i);
            fe.v = param[1] = cache.getParam(1,iel-1,j);
            fe.w = param[2] = cache.getParam(2,iel-1,k);

            const BasisFunctionVals& bfs = cache.getVals(iel-1,ig);
            fe.N = bfs.N;

            // Compute Jacobian inverse of coordinate mapping and derivatives
            fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,bfs.dNdu);
            if (fe.detJxW == 0.0) continue; // skip singular points

            // Compute Hessian of coordinate mapping and 2nd order derivatives
            if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
              if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,bfs.d2Ndu2,bfs.dNdu))
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
            PROFILE3("Integrand::evalInt");
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

#if defined(SP_DEBUG) && !defined(USE_OPENMP)
      if (iel == -dbgElm)
        break; // Skipping all elements, except for -dbgElm
#endif
    }

  cache.finalizeAssembly();
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

  // Fetch all elements on the chosen face
  std::vector<LR::Element*> edgeElms;
  lrspline->getEdgeElements(edgeElms,getFaceEnum(lIndex%10));


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
    Vec4   X(param,time.t);
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
      // Element size in parametric space
      for (int i = 0; i < 3; i++)
        dXidu[i] = el->getParmax(i) - el->getParmin(i);

    // Initialize element quantities
    size_t nen = el->support().size();
    LocalIntegral* A = integrand.getLocalIntegral(nen,fe.iel,true);
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
  std::cerr <<" *** ASMu3D::integrateEdge is not implemented yet :("<< std::endl;
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
  Go::BasisPts pt;
  this->getBasis()->computeBasis(fe.u, fe.v, fe.w, pt, iel);
  fe.N = pt.basisValues;

  LocalIntegral* A = integrand.getLocalIntegral(MNPC[iel].size(),fe.iel,true);
  bool ok = integrand.evalPoint(*A,fe,pval) && glInt.assemble(A,fe.iel);

  A->destruct();

  return ok;
}


int ASMu3D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  const LR::LRSplineVolume* geo = this->getBasis(ASM::GEOMETRY_BASIS);
  if (!geo) return -3;

  for (int i = 0; i < 3; i++)
    param[i] = (1.0-xi[i])*geo->startparam(i) + xi[i]*geo->endparam(i);

  int iel = 0;
  return this->evalPoint(iel,param,X);
}


int ASMu3D::evalPoint (int, const double* param, Vec3& X) const
{
  const LR::LRSplineVolume* geo = this->getBasis(ASM::GEOMETRY_BASIS);
  Go::Point X0;
  geo->point(X0,param[0],param[1],param[2]);
  for (int i = 0; i < 3 && i < geo->dimension(); i++)
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

  double eps = ElementBlock::eps;

  for (const LR::Element* el : lrspline->getAllElements())
  {
    // Get parametric element evaluation points, optionally with some shrinkage
    std::pair<double,double> u;
    if (dir == 0)
      u = { el->umin(), el->umax() };
    else if (dir == 1)
      u = { el->vmin(), el->vmax() };
    else
      u = { el->wmin(), el->wmax() };

    double du = (u.second - u.first)*(1.0-2.0*eps)/nSegPerSpan;
    u.first = u.first*(1.0-eps) + u.second*eps;
    for (int k = 0; k <= nSegPerSpan; k++)
      for (int j = 0; j <= nSegPerSpan; j++)
        for (int i = 0; i <= nSegPerSpan; i++)
          prm.push_back(u.first + du*double(dir == 0 ? i : (dir == 1 ? j : k)));
  }

  return true;
}


/*!
  Each nodal point in \a grid is generated once for each element using it.
  This results in a lot of unnecessary duplicates (if ElementBlock::eps is 0.0),
  but is preferable instead of figuring out all element topology information.

  Setting ElementBlock::eps > 0.0 will handle the case of internal
  C<sup>-1</sup> continuities automatically, in that result quantities along
  the discontinuity will not be unique.
*/

bool ASMu3D::tesselate (ElementBlock& grid, const int* npe) const
{
  if (!lrspline) return false;

  int nNodesPerElement =  npe[0]   * npe[1]   * npe[2];
  int nSubElPerElement = (npe[0]-1)*(npe[1]-1)*(npe[2]-1);
  int nElements        = lrspline->nElements();
  grid.unStructResize(nElements * nSubElPerElement,
                      nElements * nNodesPerElement);

  double eps = ElementBlock::eps;

  int iel = 0, inod = 0;
  for (const LR::Element* el : lrspline->getAllElements())
  {
    // Evaluate element at corner points, optionally with some shrinkage (eps)
    double umin = el->umin()*(1.0-eps) + el->umax()*eps;
    double vmin = el->vmin()*(1.0-eps) + el->vmax()*eps;
    double wmin = el->wmin()*(1.0-eps) + el->wmax()*eps;
    double du = (el->umax() - el->umin())*(1.0-2.0*eps)/(npe[0]-1);
    double dv = (el->vmax() - el->vmin())*(1.0-2.0*eps)/(npe[1]-1);
    double dw = (el->wmax() - el->wmin())*(1.0-2.0*eps)/(npe[2]-1);
    for (int iw = 0; iw < npe[2]; iw++)
      for (int iv = 0; iv < npe[1]; iv++)
        for (int iu = 0; iu < npe[0]; iu++, inod++) {
          double u = umin + du*iu;
          double v = vmin + dv*iv;
          double w = wmin + dw*iw;
          Go::Point pt;
          lrspline->point(pt, u,v,w, iel, iu!=npe[0]-1, iv!=npe[1]-1, iw!=npe[2]-1);
          grid.setCoor(inod, SplineUtils::toVec3(pt,nsd));
          grid.setParams(inod, u, v, w);
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
                           const int* npe, int nf, bool piola) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  if (piola)
    return this->evalSolutionPiola(sField,locSol,gpar.data(),false);
  else
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


bool ASMu3D::evalProjSolution (Matrix& sField, const Vector& locSol,
                               const int* npe, int nf) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the projected solution at all sampling points
  if (!this->separateProjectionBasis())
    return this->evalSolution(sField,locSol,gpar.data(),false,0,nf);

  // The projection uses a separate basis, need to interpolate
  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size() || nPoints != gpar[2].size())
    return false;

  Fields* f = this->getProjectedFields(locSol);
  if (!f) return false;

  // Evaluate the projected solution field at each point
  Vector vals;
  sField.resize(f->getNoFields(),nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    f->valueFE(ItgPoint(gpar[0][i],gpar[1][i],gpar[2][i]),vals);
    sField.fillColumn(1+i,vals);
  }

  delete f;
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
      this->evaluateBasis(iel, fe.u, fe.v, fe.w, fe.basis(1), dNdu, d2Ndu2);
    else
      this->evaluateBasis(fe.basis(1), dNdu, bezierExtract[iel], B);

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

    // Cartesian coordinates of current integration point
    utl::Point X4(Xnod*fe.N, {fe.u, fe.v, fe.w});

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,X4,MNPC[iel]))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMu3D::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                               int, int orient, bool local) const
{
  if (basis == 0)
    basis = 1;

  const LR::LRSplineVolume* vol = this->getBasis(basis);
  if (!vol) return; // silently ignore empty patches

  std::vector<LR::Basisfunction*> edgeFunctions;
  if (lIndex > 0 && lIndex <= 6) {
    this->getBasis(basis)->getEdgeFunctions(edgeFunctions, getFaceEnum(lIndex));
    if (orient >= 0) {
      int u = lIndex <= 2 ? 1 : 0;
      int v = lIndex >= 5 ? 1 : 2;
      ASMLRSpline::Sort(u, v, orient, edgeFunctions);
    }
  }

  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  for (LR::Basisfunction* b : edgeFunctions)
    nodes.push_back(local ? b->getId()+ofs : this->getNodeID(b->getId()+ofs));

#if SP_DEBUG > 1
  std::cout <<"Boundary nodes in patch "<< idx+1 <<" face "<< lIndex <<":";
  for (int n : nodes) std::cout <<" "<< n;
  std::cout << std::endl;
#endif
}


bool ASMu3D::getOrder (int& p1, int& p2, int& p3) const
{
  p1 = lrspline->order(0);
  p2 = lrspline->order(1);
  p3 = lrspline->order(2);

  return true;
}


int ASMu3D::getCorner(int I, int J, int K, int basis) const
{
  std::vector<LR::Basisfunction*> corner; // vector of one function for corner-input
  this->getBasis(basis)->getEdgeFunctions(corner,
                                          (I > 0 ? LR::EAST  : LR::WEST  ) |
                                          (J > 0 ? LR::NORTH : LR::SOUTH ) |
                                          (K > 0 ? LR::TOP   : LR::BOTTOM));
  if (corner.empty())
    return -1;

  int nodeId = corner.front()->getId() + 1;
  for (int i = 1; i < basis; i++)
    nodeId += this->getNoNodes(i);
  return nodeId;
}


void ASMu3D::getNoBouPoints (size_t& nPt, char ldim, char lindx)
{
  size_t nGp = 1;
  if (nGauss > 0 && nGauss <= 10)
    for (char d = 0; d < ldim; d++)
      nGp *= nGauss;
  else if (ldim == 2)
  {
    // Use polynomial order to define number of quadrature points
    int p[3] = { 0, 0, 0 };
    this->getOrder(p[0],p[1],p[2]);
    p[(lindx-1)/2] = 1;
    int nG = std::max(std::max(p[0],p[1]),p[2]);
    nGp = nG*nG;
  }
  else
    nGp = 0;

  firstBp[lindx] = nPt;

  nPt += this->getNoBoundaryElms(lindx,ldim)*nGp;
}


void ASMu3D::generateThreadGroups (const Integrand& integrand, bool silence,
                                   bool ignoreGlobalLM)
{
  LR::generateThreadGroups(threadGroups, this->getBasis(1));
  if (this->separateProjectionBasis())
    LR::generateThreadGroups(projThreadGroups, projB.get());
  if (silence || threadGroups[0].size() < 2) return;

  IFEM::cout <<"\nMultiple threads are utilized during element assembly.";
#if SP_DEBUG
  for (size_t i = 0; i < threadGroups[0].size(); i++)
    IFEM::cout <<"\n Color "<< i+1 <<": "
               << threadGroups[0][i].size() <<" elements";
  IFEM::cout << std::endl;
#else
  this->analyzeThreadGroups(threadGroups[0]);
#endif
}


bool ASMu3D::updateDirichlet (const std::map<int,RealFunc*>& func,
                              const std::map<int,VecFunc*>& vfunc, double time,
                              const std::map<int,int>* g2l)
{
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;

  for (const DirichletFace& dfac : dirich)
  {
    Real2DMat controlPts;
    if ((fit = func.find(dfac.code)) != func.end())
      this->faceL2projection(dfac, *fit->second, controlPts, time);
    else if ((vfit = vfunc.find(dfac.code)) != vfunc.end())
      this->faceL2projection(dfac, *vfit->second, controlPts, time);
    else
    {
      std::cerr <<" *** ASMu3D::updateDirichlet: Code "<< dfac.code
                <<" is not associated with any function."<< std::endl;
      return false;
    }
    if (controlPts.empty())
    {
      std::cerr <<" *** ASMu3D::updateDirichlet: Projection failure."
                << std::endl;
      return false;
    }

    // Loop over the (non-corner) nodes of this boundary face
    for (size_t j = 0; j < dfac.MLGN.size(); j++)
      if (!dfac.isCorner(dfac.MLGN[j]))
        for (int dofs = dfac.dof; dofs > 0; dofs /= 10)
        {
          int dof = dofs%10;
          // Find the constraint equation for current (node,dof)
          MPC pDOF(MLGN[dfac.MLGN[j]-1],dof);
          MPCIter mit = mpcs.find(&pDOF);
          if (mit != mpcs.end())
          {
            // Now update the prescribed value in the constraint equation
            if (dfac.dof < 10) dof = 1; // scalar condition
            (*mit)->setSlaveCoeff(controlPts[dof-1][j]);
#if SP_DEBUG > 1
            std::cout <<"Updated constraint: "<< **mit;
#endif
          }
        }
  }

  // The parent class method takes care of the corner nodes with direct
  // evaluation of the Dirichlet functions; since they are interpolatory
  return this->ASMbase::updateDirichlet(func,vfunc,time,g2l);
}


size_t ASMu3D::getNoNodes (int basis) const
{
  if (basis == 0)
    return this->ASMbase::getNoNodes(0);
  else
    return lrspline->nBasisFunctions();
}


size_t ASMu3D::getNoProjectionNodes () const
{
  return projB->nBasisFunctions();
}


bool ASMu3D::separateProjectionBasis () const
{
  return this->getBasis(ASM::PROJECTION_BASIS) != this->getBasis(1);
}


Field* ASMu3D::getProjectedField (const Vector& coefs) const
{
  if (coefs.size() == this->getNoProjectionNodes())
    return new LRSplineField3D(this->getBasis(ASM::PROJECTION_BASIS),coefs);

  std::cerr <<" *** ASMu3D::getProjectedFields: Non-matching coefficent array,"
            <<" size="<< coefs.size() <<" nnod="<< this->getNoProjectionNodes()
            << std::endl;
  return nullptr;
}


Fields* ASMu3D::getProjectedFields (const Vector& coefs, size_t) const
{
  if (!this->separateProjectionBasis())
    return nullptr;

  size_t ncmp = coefs.size() / this->getNoProjectionNodes();
  if (ncmp*this->getNoProjectionNodes() == coefs.size())
    return new LRSplineFields3D(this->getBasis(ASM::PROJECTION_BASIS),coefs,ncmp);

  std::cerr <<" *** ASMu3D::getProjectedFields: Non-matching coefficent array,"
            <<" size="<< coefs.size() <<" nnod="<< this->getNoProjectionNodes()
            << std::endl;
  return nullptr;
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

  for (int iel = 1; iel <= newBasis->nElements(); iel++)
  {
    RealArray U, V, W, ptVar;
    LR::getGaussPointParameters(newBasis, U, 0, nGauss, iel, xi);
    LR::getGaussPointParameters(newBasis, V, 1, nGauss, iel, xi);
    LR::getGaussPointParameters(newBasis, W, 2, nGauss, iel, xi);
    for (int k = 0; k < nGauss; k++)
      for (int j = 0; j < nGauss; j++)
        for (int i = 0; i < nGauss; i++)
        {
          oldBasis->point(ptVar,U[i],V[j],W[k],iel-1);
          newVars.insert(newVars.end(),ptVar.begin(),ptVar.end());
        }
  }

  return true;
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

  std::vector<IntVec> bndry1(nedge);
  for (int j = 1; j <= nedge; j++)
    this->getBoundary1Nodes(j , bndry1[j-1], 1, 0, true);

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


bool ASMu3D::refine (const LR::RefineData& prm, Vectors& sol)
{
  if (!this->ASMLRSpline::refine(prm,sol))
    return false;

  // check if refinement was actually done
  if (prm.elements.size() + prm.errors.size() == 0)
    return true;

  if (!this->separateProjectionBasis())
    return true;

  LR::LRSplineVolume* proj = this->getBasis(ASM::PROJECTION_BASIS);
  for (const LR::MeshRectangle* rect : lrspline->getAllMeshRectangles())
    proj->insert_line(rect->copy());

  proj->generateIDs();

  IFEM::cout <<"Refined projection basis: "<< proj->nElements()
             <<" elements "<< proj->nBasisFunctions() <<" nodes."
             << std::endl;

  return true;
}


void ASMu3D::getElmConnectivities (IntMat& neigh) const
{
  const LR::LRSplineVolume* lr = this->getBasis(1);
  for (const LR::Element* m : lr->getAllElements()) {
    int iel = m->getId();
    IntVec& neighbor = neigh[MLGE[iel]-1];
    for (int face = 1; face <= 6; face++) {
      std::set<int> elms = lr->getElementNeighbours(iel,getFaceEnum(face));
      for (int elm : elms)
        neighbor.push_back(MLGE[elm]-1);
    }
  }
}


void ASMu3D::findBoundaryElms (IntVec& elms, int lIndex, int orient) const
{
  elms.clear();
  if (lIndex < 1 || lIndex > 6)
    return;

  std::vector<LR::Element*> elements;
  this->getBasis(1)->getEdgeElements(elements,getFaceEnum(lIndex));

  if (orient >= 0)
    std::sort(elements.begin(), elements.end(),
              [lIndex,orient](const LR::Element* a, const LR::Element* b)
              {
                int u = lIndex <= 2 ? 1 : 0;
                int v = lIndex >= 5 ? 1 : 2;
                int idx = (orient & 4) ? v : u;
                std::vector<double> A = a->midpoint();
                std::vector<double> B = b->midpoint();
                if (A[idx] != B[idx])
                  return (orient & 2) ? A[idx] > B[idx] : A[idx] < B[idx];

                idx = (orient & 4) ? u : v;
                if (A[idx] != B[idx])
                  return (orient & 1) ? A[idx] > B[idx] : A[idx] < B[idx];

                return false;
              });

  for (const LR::Element* elem : elements)
    elms.push_back(elem->getId());
}


void ASMu3D::generateThreadGroupsFromElms (const IntVec& elms)
{
  myElms.clear();
  for (int elm : elms)
    if (this->getElmIndex(elm+1) > 0)
      myElms.push_back(this->getElmIndex(elm+1)-1);

  // We need myElms to be non-empty to flag that partitioning is enabled
  if (myElms.empty())
    myElms.push_back(-1);

  if (projThreadGroups.size() == 0 || projThreadGroups[0].empty())
    projThreadGroups = threadGroups;

  threadGroups = threadGroups.filter(myElms);
}


ASMu3D::BasisFunctionCache::BasisFunctionCache (const ASMu3D& pch,
                                                ASM::CachePolicy plcy,
                                                int b, bool useBezier) :
  ::BasisFunctionCache<3>(plcy),
  bezierEnabled(useBezier),
  patch(pch),
  basis(b)
{
}


ASMu3D::BasisFunctionCache::BasisFunctionCache (const BasisFunctionCache& cache,
                                                int b) :
  ::BasisFunctionCache<3>(cache),
  bezierEnabled(cache.bezierEnabled),
  patch(cache.patch),
  basis(b)
{
}


bool ASMu3D::BasisFunctionCache::internalInit ()
{
  if (!mainQ->xg[0])
    this->setupQuadrature();

  std::array<int,3> order = { patch.getBasis(basis)->order(0),
                              patch.getBasis(basis)->order(1),
                              patch.getBasis(basis)->order(2) };
  // Evaluate all gauss points on the bezier patch (-1, 1)
  auto&& extractBezier = [order](const Quadrature& q,
                                 BezierExtract& b)
  {
    double u[2*order[0]], v[2*order[1]], w[2*order[2]];
    Go::BsplineBasis basis1 = getBezierBasis(order[0]);
    Go::BsplineBasis basis2 = getBezierBasis(order[1]);
    Go::BsplineBasis basis3 = getBezierBasis(order[2]);
    int P = order[0]*order[1]*order[2];
    int N = q.ng[0]*q.ng[1]*q.ng[2];
    b.N.resize(P,N);
    b.dNdu.resize(P,N);
    b.dNdv.resize(P,N);
    b.dNdw.resize(P,N);
    int ig = 1; // gauss point iterator
    for (int zeta = 0; zeta < q.ng[2]; zeta++)
      for (int eta = 0; eta < q.ng[1]; eta++)
        for (int xi = 0; xi < q.ng[0]; xi++, ig++) {
          basis1.computeBasisValues(q.xg[0][xi],   u, 1);
          basis2.computeBasisValues(q.xg[1][eta],  v, 1);
          basis3.computeBasisValues(q.xg[2][zeta], w, 1);
          int ib = 1; // basis function iterator
          for (int k = 0; k < order[2]; k++)
            for (int j = 0; j < order[1]; j++)
              for (int i = 0; i < order[0]; i++, ib++) {
                b.N(ib,ig)    = u[2*i  ]*v[2*j  ]*w[2*k  ];
                b.dNdu(ib,ig) = u[2*i+1]*v[2*j  ]*w[2*k  ];
                b.dNdv(ib,ig) = u[2*i  ]*v[2*j+1]*w[2*k  ];
                b.dNdw(ib,ig) = u[2*i  ]*v[2*j  ]*w[2*k+1];
              }
        }
  };

  if (bezierEnabled) {
    extractBezier(*mainQ, mainB);
    if (reducedQ->xg[0])
      extractBezier(*reducedQ, reducedB);
  }

  nTotal = patch.nel*mainQ->ng[0]*mainQ->ng[1]*mainQ->ng[2];
  if (reducedQ->xg[0])
    nTotalRed = patch.nel*reducedQ->ng[0]*reducedQ->ng[1]*reducedQ->ng[2];

  return true;
}


void ASMu3D::BasisFunctionCache::internalCleanup ()
{
  if (basis == 1) {
    mainQ->reset();
    reducedQ->reset();
  }
}


bool ASMu3D::BasisFunctionCache::setupQuadrature ()
{
  std::array<int,3> order = { patch.getBasis(basis)->order(0),
                              patch.getBasis(basis)->order(1),
                              patch.getBasis(basis)->order(2) };

  // Get Gaussian quadrature points and weights
  for (int d = 0; d < 3; d++)
  {
    mainQ->ng[d] = patch.getNoGaussPt(order[d]);
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
  } else
    nRed = mainQ->ng[0];

  reducedQ->ng[0] = reducedQ->ng[1] = reducedQ->ng[2] = nRed;

  // Compute parameter values of the Gauss points over the whole patch
  mainQ->gpar[0].resize(mainQ->ng[0],patch.nel);
  mainQ->gpar[1].resize(mainQ->ng[1],patch.nel);
  mainQ->gpar[2].resize(mainQ->ng[2],patch.nel);
  if (reducedQ->xg[0]) {
    reducedQ->gpar[0].resize(reducedQ->ng[0],patch.nel);
    reducedQ->gpar[1].resize(reducedQ->ng[1],patch.nel);
    reducedQ->gpar[2].resize(reducedQ->ng[2],patch.nel);
  }
  for (size_t iel = 1; iel <= patch.nel; ++iel)
  {
    RealArray u, v, w;
    patch.getGaussPointParameters(u,0,mainQ->ng[0],iel,mainQ->xg[0]);
    patch.getGaussPointParameters(v,1,mainQ->ng[1],iel,mainQ->xg[1]);
    patch.getGaussPointParameters(w,2,mainQ->ng[2],iel,mainQ->xg[2]);
    mainQ->gpar[0].fillColumn(iel,u.data());
    mainQ->gpar[1].fillColumn(iel,v.data());
    mainQ->gpar[2].fillColumn(iel,w.data());

    if (reducedQ->xg[0])
    {
      patch.getGaussPointParameters(u,0,reducedQ->ng[0],iel,reducedQ->xg[0]);
      patch.getGaussPointParameters(v,1,reducedQ->ng[1],iel,reducedQ->xg[0]);
      patch.getGaussPointParameters(w,2,reducedQ->ng[2],iel,reducedQ->xg[0]);
      reducedQ->gpar[0].fillColumn(iel,u.data());
      reducedQ->gpar[1].fillColumn(iel,v.data());
      reducedQ->gpar[2].fillColumn(iel,w.data());
    }
  }
  return true;
}


BasisFunctionVals ASMu3D::BasisFunctionCache::calculatePt (size_t el,
                                                           size_t gp,
                                                           bool reduced) const
{
  PROFILE2("Spline evaluation");
  const std::array<size_t,3> gpIdx = this->gpIndex(gp,reduced);
  FiniteElement fe;
  fe.u = this->getParam(0,el,gpIdx[0],reduced);
  fe.v = this->getParam(1,el,gpIdx[1],reduced);
  fe.w = this->getParam(2,el,gpIdx[2],reduced);

  const LR::Element* elm = patch.lrspline->getElement(el);
  std::array<double,3> du;
  du[0] = 0.5*(elm->umax() - elm->umin());
  du[1] = 0.5*(elm->vmax() - elm->vmin());
  du[2] = 0.5*(elm->wmax() - elm->wmin());

  if (patch.lrspline.get() != patch.getBasis(basis))
    el = patch.getBasis(basis)->getElementContaining(elm->midpoint());

  return this->calculatePrm(fe,du,el,gp,reduced);
}


BasisFunctionVals ASMu3D::BasisFunctionCache::
calculatePrm (FiniteElement& fe,
              const std::array<double,3>& du,
              size_t el, size_t gp, bool reduced) const
{
  BasisFunctionVals result;
  if (nderiv == 1 || reduced) {
    if (bezierEnabled) {
      const BezierExtract& b = reduced ? reducedB : mainB;
      RealArray extrMat;
      patch.getBasis(basis)->getBezierExtraction(el,extrMat);
      Matrix C;
      const LR::Element* elm = patch.getBasis(basis)->getElement(el);
      C.resize(elm->nBasisFunctions(), b.N.rows());
      C.fill(extrMat.data(),extrMat.size());
      Matrix B(b.N.rows(), 4);
      B.fillColumn(1, b.N.getColumn(gp+1));
      B.fillColumn(2, b.dNdu.getColumn(gp+1) / du[0]);
      B.fillColumn(3, b.dNdv.getColumn(gp+1) / du[1]);
      B.fillColumn(4, b.dNdw.getColumn(gp+1) / du[2]);

      patch.evaluateBasis(result.N, result.dNdu, C, B);
    } else
      patch.evaluateBasis(el, fe.u, fe.v, fe.w, result.N, result.dNdu, basis);
  } else if (nderiv == 2)
    patch.evaluateBasis(el, fe.u, fe.v, fe.w,
                        result.N, result.dNdu, result.d2Ndu2, basis);

  return result;
}


void ASMu3D::BasisFunctionCache::calculateAll ()
{
  PROFILE2("Spline evaluation");
  // Evaluate basis function values and derivatives at all integration points.
  // We do this before the integration point loop to exploit multi-threading
  // in the integrand evaluations, which may be the computational bottleneck.
  size_t iel, jp, rp;
  for (iel = jp = rp = 0; iel < patch.nel; iel++)
  {
    for (int k = 0; k < mainQ->ng[2]; k++)
      for (int j = 0; j < mainQ->ng[1]; j++)
        for (int i = 0; i < mainQ->ng[0]; i++, jp++)
          values[jp] = this->calculatePt(iel,(k*mainQ->ng[1]+j)*mainQ->ng[0]+i,false);

    if (reducedQ->xg[0])
      for (int k = 0; k < reducedQ->ng[2]; k++)
        for (int j = 0; j < reducedQ->ng[1]; j++)
          for (int i = 0; i < reducedQ->ng[0]; i++, rp++)
            valuesRed[rp] = this->calculatePt(iel,
                                              (k*reducedQ->ng[1]+j)*reducedQ->ng[0]+i,true);
  }
}
