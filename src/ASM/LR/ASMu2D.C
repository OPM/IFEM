// $Id$
//==============================================================================
//!
//! \file ASMu2D.C
//!
//! \date September 2011
//!
//! \author Kjetil Andre Johannessen / SINTEF
//!
//! \brief Driver for assembly of unstructured 2D spline FE models.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

#include "ASMu2D.h"
#include "TimeDomain.h"
#include "FiniteElement.h"
#include "GlobalIntegral.h"
#include "GlobalNodes.h"
#include "LocalIntegral.h"
#include "IntegrandBase.h"
#include "CoordinateMapping.h"
#include "GaussQuadrature.h"
#include "LagrangeInterpolator.h"
#include "LRSplineField2D.h"
#include "LRSplineFields2D.h"
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
#include <fstream>
#include <utility>


ASMu2D::ASMu2D (unsigned char n_s, unsigned char n_f)
  : ASMLRSpline(2,n_s,n_f), lrspline(nullptr),
    bezierExtract(myBezierExtract)
{
  aMin = 0.0;
  tensorspline = tensorPrjBas = nullptr;
}


ASMu2D::ASMu2D (const ASMu2D& patch, unsigned char n_f)
  : ASMLRSpline(patch,n_f), lrspline(patch.lrspline),
    bezierExtract(patch.myBezierExtract)
{
  aMin = 0.0;
  tensorspline = tensorPrjBas = nullptr;
  is_rational = patch.is_rational;

  // Need to set nnod here,
  // as hasXNodes might be invoked before the FE data is generated
  if (nnod == 0 && lrspline)
    nnod = lrspline->nBasisFunctions();
}


const LR::LRSplineSurface* ASMu2D::getBasis (int basis) const
{
  switch (basis) {
    case ASM::GEOMETRY_BASIS:
      return static_cast<const LR::LRSplineSurface*>(geomB.get());
    case ASM::PROJECTION_BASIS:
      return static_cast<const LR::LRSplineSurface*>(projB.get());
    case ASM::PROJECTION_BASIS_2:
      return static_cast<const LR::LRSplineSurface*>(projB2.get());
    case ASM::REFINEMENT_BASIS:
      return static_cast<const LR::LRSplineSurface*>(refB.get());
    default:
      return lrspline.get();
  }
}


LR::LRSplineSurface* ASMu2D::getBasis (int basis)
{
  if (tensorspline)
    this->createLRfromTensor();

  return const_cast<LR::LRSplineSurface*>(std::as_const(*this).getBasis(basis));
}


bool ASMu2D::read (std::istream& is)
{
  if (shareFE) return true;

  // Read the input file as either an LRSpline file directly,
  // or a tensor product B-spline and convert
  char firstline[256];
  is.getline(firstline, 256);
  if (strncmp(firstline, "# LRSPLINE", 10) == 0) {
    lrspline.reset(new LR::LRSplineSurface());
    is >> *lrspline;
  }
  else {
    // Probably a SplineSurface, so we'll read that and convert
    tensorspline = new Go::SplineSurface();
    is >> *tensorspline;
    lrspline.reset(new LR::LRSplineSurface(tensorspline));
  }

  // Eat white-space characters to see if there is more data to read
  char c;
  while (is.get(c))
    if (!isspace(c)) {
      is.putback(c);
      break;
    }

  int readDim = lrspline->dimension();
  if (!is.good() && !is.eof())
  {
    std::cerr <<" *** ASMu2D::read: Failure reading spline data"<< std::endl;
    lrspline.reset();
    return false;
  }
  else if (readDim < 2)
  {
    std::cerr <<" *** ASMu2D::read: Invalid lrspline patch, dim="
              << readDim << std::endl;
    lrspline.reset();
    return false;
  }
  else if (readDim < nsd)
  {
    std::cout <<"  ** ASMu2D::read: The dimension of this lrspline patch "
              << readDim <<" is less than nsd="<< nsd
              <<".\n                   Resetting nsd to "<< readDim
              <<" for this patch."<< std::endl;
    nsd = readDim;
  }

  geomB = lrspline;
  return true;
}


bool ASMu2D::write (std::ostream& os, int basis) const
{
  if (!lrspline) return false;
  if (basis > static_cast<int>(this->getNoBasis())) return false;
  const LR::LRSplineSurface* spline = this->getBasis(basis);
  if (!spline) return false;

  os << *spline;

  return os.good();
}


void ASMu2D::clear (bool retainGeometry)
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


bool ASMu2D::cornerRefine (int minBasisfunctions)
{
  if (!lrspline) return false;
  if (shareFE) return true;

  double h = 1.0;
  int nBasis = lrspline->nBasisFunctions();
  double unif_step_h = 1.0 / ((minBasisfunctions - nBasis) / 3.0 + 1.0);
  while(lrspline->nBasisFunctions() < minBasisfunctions) {
    lrspline->insert_const_u_edge(h-unif_step_h, 0, h);
    lrspline->insert_const_v_edge(h-unif_step_h, 0, h);
    h -= unif_step_h;
  }

  std::ofstream paramMeshFile("mesh_param.eps");
  std::ofstream physicalMeshFile("mesh_physical.eps");
  lrspline->writePostscriptMesh(paramMeshFile);
  lrspline->writePostscriptMesh(physicalMeshFile);
  return true;
}

bool ASMu2D::diagonalRefine (int minBasisfunctions)
{
  if (!lrspline) return false;
  if (shareFE) return true;

  double end1 = lrspline->endparam(0);
  double end2 = lrspline->endparam(1);
  double h = 1.0;
  int iter = 0;
  double u = h/2.0;
  double v = h/2.0;
  while(lrspline->nBasisFunctions() < minBasisfunctions) {
    lrspline->insert_const_u_edge(u, (iter-1<0) ? 0 : (iter-1)*h, ((iter+2)*h>end2) ? end2 : (iter+2)*h);
    lrspline->insert_const_v_edge(v, (iter-1<0) ? 0 : (iter-1)*h, ((iter+2)*h>end1) ? end1 : (iter+2)*h);
    u += h;
    v += h;
    iter++;
    if ( u>end1 ) {
      h /= 2.0;
      iter = 0;
      u = h/2.0;
      v = h/2.0;
    }
  }

  std::ofstream meshFile("mesh.eps");
  lrspline->writePostscriptMesh(meshFile);
  return true;
}

bool ASMu2D::uniformRefine (int minBasisfunctions)
{
  if (!lrspline) return false;
  if (shareFE) return true;

  double end1 = lrspline->endparam(0);
  double end2 = lrspline->endparam(1);
  double h = 1.0;
  bool step_u = true;
  double u = h/2.0;
  double v = h/2.0;
  while(lrspline->nBasisFunctions() < minBasisfunctions) {
    if (step_u) {
      lrspline->insert_const_u_edge(u, 0, end2);
      u += h;
      if (u > end1) {
        step_u = !step_u;
        u = h/4.0;
      }
    } else {
      lrspline->insert_const_v_edge(v, 0, end1);
      v += h;
      if (v > end2) {
        step_u = !step_u;
        v = h/4.0;
        h /= 2.0;
      }
    }
  }

  std::ofstream meshFile("mesh.eps");
  lrspline->writePostscriptMesh(meshFile);
  return true;
}


bool ASMu2D::uniformRefine (int dir, int nInsert)
{
  if (!tensorspline || dir < 0 || dir > 1 || nInsert < 1) return false;
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

  if (dir == 0)
    tensorspline->insertKnot_u(extraKnots);
  else
    tensorspline->insertKnot_v(extraKnots);

  return true;
}


bool ASMu2D::refine (int dir, const RealArray& xi, double scale)
{
  if (!tensorspline || dir < 0 || dir > 1 || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > scale || scale < 1.0) return false;
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
          extraKnots.push_back(ucurr*xi[i]/scale + uprev*(1.0-xi[i]/scale));

    uprev = ucurr;
  }

  if (dir == 0)
    tensorspline->insertKnot_u(extraKnots);
  else
    tensorspline->insertKnot_v(extraKnots);

  return true;
}


double ASMu2D::getMinimumSize (int nrefinements) const
{
  if (aMin > 0.0 || lrspline->nElements() <= 0)
    return aMin;

  double redMax = pow(2.0,nrefinements);
  aMin = lrspline->getElement(0)->area()/(redMax*redMax);
  return aMin;
}


bool ASMu2D::checkElementSize (int elmId, bool globalNum) const
{
  if (globalNum)
    elmId = utl::findIndex(MLGE,1+elmId);

  if (elmId >= 0 && elmId < lrspline->nElements())
    return lrspline->getElement(elmId)->area() > aMin+1.0e-12;
  else
    return false;
}


void ASMu2D::extendRefinementDomain (IntSet& refineIndices,
                                     const IntSet& neighborIndices) const
{
  // OPTIMIZATION NOTE: If we by some clever data structures already knew
  // which edge each node in conformingIndices was on, then we don't have
  // to brute-force search for it like we do here.
  // getBoundaryNodes() seems to compute this, but it is only for
  // the sending patch boundary index, not the recieving patch boundary index.

  IntVec              bndry0;
  std::vector<IntVec> bndry1(4);
  for (int i = 1; i <= 4; i++)
  {
    bndry0.push_back(GlobalNodes::getBoundaryNodes(*this->getBasis(ASM::REFINEMENT_BASIS),
                                                   0, i, 0).front());
    bndry1[i-1] = GlobalNodes::getBoundaryNodes(*this->getBasis(ASM::REFINEMENT_BASIS),
                                                  1, i, 0);
  }

  // Add refinement from neighbors
  for (int j : neighborIndices)
  {
    bool done_with_this_node = false;

    // Check if node is a corner node,
    // compute large extended domain (all directions)
    for (int edgeNode : bndry0)
      if (edgeNode == j)
      {
        IntVec secondary = this->getOverlappingNodes(j);
        refineIndices.insert(secondary.begin(),secondary.end());
        done_with_this_node = true;
        break;
      }

    // Check if node is an edge node,
    // compute small extended domain (one direction)
    for (int edge = 0; edge < 4 && !done_with_this_node; edge++)
      for (int edgeNode : bndry1[edge])
        if (edgeNode == j)
        {
          IntVec secondary = this->getOverlappingNodes(j, edge/2+1);
          refineIndices.insert(secondary.begin(),secondary.end());
          done_with_this_node = true;
          break;
        }
  }
}


bool ASMu2D::raiseOrder (int ru, int rv)
{
  if (!tensorspline) return false;
  if (shareFE) return true;

  tensorspline->raiseOrder(ru,rv);
  return true;
}


/*!
  This method is supposed to be invoked twice during the model generation.
  In the first call, with \a init = \e true, the spline surface object
  is cloned and the two pointers are then swapped, such that the subsequent
  refine and raiseOrder operations will apply to the projection basis
  and not on the geometry basis.
  In the second call, the pointers are swapped back.

  The method can also be invoked twice with \a init = \e false in case the
  projection basis is to be read from a file.
*/

bool ASMu2D::createProjectionBasis (bool init)
{
  if (!tensorspline)
    return false;
  else if (init && !tensorPrjBas)
    tensorPrjBas = tensorspline->clone();

  std::swap(tensorspline,tensorPrjBas);
  std::swap(geomB,projB);
  lrspline = std::static_pointer_cast<LR::LRSplineSurface>(geomB);
  return true;
}


bool ASMu2D::evaluateBasis (int iel, FiniteElement& fe, int derivs) const
{
  if (is_rational)
    return this->evaluateBasisNurbs(iel, fe, derivs);

  PROFILE3("ASMu2D::evalBasis");
#ifdef INDEX_CHECK
  if (iel < 0 || iel >= lrspline->nElements())
  {
    std::cerr <<" *** ASMu2D::evaluateBasis: Element index "<< 1+iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return false;
  }
#endif

  const LR::Element* el = lrspline->getElement(iel);
  fe.xi  = 2.0*(fe.u - el->umin()) / (el->umax() - el->umin()) - 1.0;
  fe.eta = 2.0*(fe.v - el->vmin()) / (el->vmax() - el->vmin()) - 1.0;
  RealArray Nu, Nv;
#pragma omp critical
  {
    Nu = bezier_u.computeBasisValues(fe.xi, derivs);
    Nv = bezier_v.computeBasisValues(fe.eta,derivs);
  }

  Vector B(lrspline->order(0)*lrspline->order(1)); // Bezier basis functions
  const Matrix& C = bezierExtract[iel];

  ++derivs;
  size_t i, j, k;
  for (j = k = 0; j < Nv.size(); j += derivs)
    for (i = 0; i < Nu.size(); i += derivs, k++)
      B[k] = Nu[i]*Nv[j];

  fe.N = C*B;

#ifdef SP_DEBUG
  if (fabs(fe.N.sum()-1.0) > 1.0e-10) {
    std::cerr <<"fe.N do not sum to one at integration point #"
              << fe.iGP << std::endl;
    return false;
  }
  else if (fabs(B.sum()-1.0) > 1.0e-10) {
    std::cerr <<"Bezier basis do not sum to one at integration point #"
              << fe.iGP << std::endl;
    return false;
  }
#endif

  if (derivs <= 1)
    return true;

  Vector Bu(B.size()), Bv(B.size()); // Bezier basis function derivatives
  for (j = k = 0; j < Nv.size(); j += derivs)
    for (i = 0; i < Nu.size(); i += derivs, k++) {
      Bu[k] = Nu[i+1]*Nv[j  ];
      Bv[k] = Nu[i  ]*Nv[j+1];
    }

  Matrix dNdu(fe.N.size(),2);
  dNdu.fillColumn(1,C*Bu);
  dNdu.fillColumn(2,C*Bv);

#ifdef SP_DEBUG
  if (fabs(dNdu.sum(-1)) > 1.0e-10) {
    std::cerr <<"dNdu do not sum to zero at integration point #"
              << fe.iGP << std::endl;
    return false;
  }
  else if (fabs(dNdu.sum(-2)) > 1.0e-10) {
    std::cerr <<"dNdv do not sums to zero at integration point #"
              << fe.iGP << std::endl;
    return false;
  }
  else if (fabs(Bu.sum()) > 1.0e-10 || fabs(Bv.sum()) > 1.0e-10) {
    std::cerr <<"Bezier derivatives do not sum to zero at integration point #"
              << fe.iGP << std::endl;
    return false;
  }
#endif

  Matrix Xnod, Jac;
  this->getElementCoordinates(Xnod,1+iel);
  fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

  return true;
}


std::shared_ptr<LR::LRSplineSurface>
ASMu2D::createLRNurbs (const Go::SplineSurface& srf)
{
  return std::make_shared<LR::LRSplineSurface>(srf.numCoefs_u(), srf.numCoefs_v(),
                                               srf.order_u(), srf.order_v(),
                                               srf.basis_u().begin(),
                                               srf.basis_v().begin(),
                                               srf.rcoefs_begin(),
                                               srf.dimension()+1);
}


std::shared_ptr<LR::LRSplineSurface> ASMu2D::createLRfromTensor ()
{
  if (tensorspline)
  {
    if (tensorspline->rational())
    {
      lrspline = createLRNurbs(*tensorspline);
      is_rational = true;
    }
    else if (tensorspline->dimension() > nsd)
    {
      // Remove superfluous components from the tensor spline coefficients.
      // We need to do this because getElementCoordinates uses the dimension
      // of the lrspline object to determine whether it is rational or not.
      RealArray::iterator it = tensorspline->coefs_begin();
      RealArray coefs;
      coefs.reserve(tensorspline->numCoefs_u()*tensorspline->numCoefs_v()*nsd);
      for (; it != tensorspline->coefs_end(); it += tensorspline->dimension())
        coefs.insert(coefs.end(),it,it+nsd);
      lrspline.reset(new LR::LRSplineSurface(tensorspline->numCoefs_u(),
                                             tensorspline->numCoefs_v(),
                                             tensorspline->order_u(),
                                             tensorspline->order_v(),
                                             tensorspline->basis_u().begin(),
                                             tensorspline->basis_v().begin(),
                                             coefs.begin(), nsd));
    }
    else
      lrspline.reset(new LR::LRSplineSurface(tensorspline));

    delete tensorspline;
    tensorspline = nullptr;
  }

  return lrspline;
}


bool ASMu2D::generateFEMTopology ()
{
  refB = geomB = this->createLRfromTensor();

  if (tensorPrjBas)
  {
    projB = tensorPrjBas->rational() ? createLRNurbs(*tensorPrjBas)
                                     : std::make_shared<LR::LRSplineSurface>(tensorPrjBas);
    projB->generateIDs();
    delete tensorPrjBas;
    tensorPrjBas = nullptr;
  }
  else if (!projB)
    projB = lrspline;

  if (!lrspline) return false;

  nnod = lrspline->nBasisFunctions();
  nel  = lrspline->nElements();

  this->generateBezierBasis();

  if (!MLGN.empty()) {
    if (MLGN.size() != nnod)
    {
      std::cerr <<" *** ASMu2D::generateFEMTopology: Inconsistency"
                <<" between the number of FE nodes "<< MLGN.size()
                <<"\n     and the number of basis functions "<< nnod
                <<" in the patch."<< std::endl;
      return false;
    }
    return true;
  }

  if (shareFE) return true;

  myMLGN.resize(nnod);
  myMLGE.resize(nel);
  myMNPC.resize(nel);

  lrspline->generateIDs();

  size_t iel = 0;
  for (const LR::Element* elm : lrspline->getAllElements())
  {
    myMLGE[iel] = ++gEl; // global element number over all patches
    myMNPC[iel].resize(elm->nBasisFunctions());

    int lnod = 0;
    for (LR::Basisfunction* b : elm->support())
      myMNPC[iel][lnod++] = b->getId();
    ++iel;
  }

  for (size_t inod = 0; inod < nnod; inod++)
    myMLGN[inod] = ++gNod;

  this->generateBezierExtraction();

  return true;
}


bool ASMu2D::connectPatch (int edge, ASM2D& neighbor, int nedge,
                           bool revers, int, bool coordCheck, int thick)
{
  ASMu2D* neighU = dynamic_cast<ASMu2D*>(&neighbor);
  if (!neighU)
    return false;

  if (!this->connectBasis(edge,*neighU,nedge,revers,1,0,0,coordCheck,thick))
    return false;

  this->addNeighbor(neighU);
  return true;
}


bool ASMu2D::connectBasis (int edge, ASMu2D& neighbor, int nedge, bool revers,
                           int basis, int slave, int master,
                           bool coordCheck, int thick)
{
  if (this->isShared() && neighbor.isShared())
    return true;
  else if (this->isShared() || neighbor.isShared())
  {
    std::cerr <<" *** ASMu2D::connectBasis: Logic error, cannot"
              <<" connect a shared patch with an unshared one"<< std::endl;
    return false;
  }

  // Set up the slave node numbers for this surface patch
  IntVec slaveNodes;
  this->getBoundaryNodes(edge, slaveNodes, basis, thick, 0, true);
  for (int& it : slaveNodes)
    it += slave;

  // Set up the master node numbers for the neighboring surface patch
  IntVec masterNodes;
  neighbor.getBoundaryNodes(nedge, masterNodes, basis, thick, 0, true);
  for (int& it : masterNodes)
    it += master;

  if (masterNodes.empty() || masterNodes.size() != slaveNodes.size())
  {
    std::cerr <<" *** ASMu2D::connectBasis: Non-matching edges, sizes "
              << masterNodes.size() <<" and "<< slaveNodes.size() << std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  for (size_t i = 0; i < masterNodes.size(); ++i)
  {
    int node = masterNodes[i];
    int slvn = slaveNodes[revers ? slaveNodes.size()-i-1 : i];
    if (!coordCheck)
      ASMbase::collapseNodes(neighbor,node,*this,slvn);
    else if (neighbor.getCoord(node).equal(this->getCoord(slvn),xtol))
      ASMbase::collapseNodes(neighbor,node,*this,slvn);
    else
    {
      std::cerr <<" *** ASMu2D::connectBasis: Non-matching nodes "
                << node <<": "<< neighbor.getCoord(node)
                <<"\n                                          and "
                << slvn <<": "<< this->getCoord(slvn) << std::endl;
      return false;
    }
  }

  return true;
}


ASMu2D::DirichletEdge::DirichletEdge (LR::LRSplineSurface* sf,
                                      int dir, int d, int c, int offset)
  : lr(sf), edg(LR::NONE), dof(d), code(c)
{
  // Figure out what edge we are at, and
  // find the corners since these are not to be included in the L2-fitting
  // of the inhomogenuous dirichlet boundaries; corners are interpolatory.
  // Optimization note: loop over the "edge"-container to manually pick up
  // the end nodes. LRspline::getEdgeFunctions() does a global search.
  std::vector<LR::Basisfunction*> c1, c2;
  switch (dir)
  {
  case -2:
    edg = LR::SOUTH;
    lr->getEdgeFunctions(c1, LR::SOUTH_WEST);
    lr->getEdgeFunctions(c2, LR::SOUTH_EAST);
    break;
  case -1:
    edg = LR::WEST;
    lr->getEdgeFunctions(c1, LR::SOUTH_WEST);
    lr->getEdgeFunctions(c2, LR::NORTH_WEST);
    break;
  case 1:
    edg = LR::EAST;
    lr->getEdgeFunctions(c1, LR::NORTH_EAST);
    lr->getEdgeFunctions(c2, LR::SOUTH_EAST);
    break;
  case 2:
    edg = LR::NORTH;
    lr->getEdgeFunctions(c1, LR::NORTH_WEST);
    lr->getEdgeFunctions(c2, LR::NORTH_EAST);
    break;
  default:
    corners[0] = corners[1] = 0;
    return;
  }

  corners[0] = c1.front()->getId() + offset;
  corners[1] = c2.front()->getId() + offset;
}


/*!
  A negative \a code value implies direct evaluation of the Dirichlet condition
  function at the control point. Positive \a code implies projection onto the
  spline basis representing the boundary curve (needed for curved edges and/or
  non-constant functions).
*/

void ASMu2D::constrainEdge (int dir, bool open, int dof, int code, char basis)
{
  if (basis < 1) basis = 1;

  // Figure out function index offset (when using multiple basis)
  int offset = 1;
  for (int i = 1; i < basis; i++)
    offset += this->getNoNodes(i);

  // Figure out what edge we are at
  DirichletEdge de(this->getBasis(basis), dir, dof, code, offset);

  // Get all basis functions on this edge
  std::vector<LR::Basisfunction*> edgeFunctions;
  de.lr->getEdgeFunctions(edgeFunctions,de.edg);

  // Add constraints for all basis functions on the edge
  for (LR::Basisfunction* b : edgeFunctions)
    if (!de.isCorner(b->getId()+offset))
      this->prescribe(b->getId()+offset, dof, -code);
    else if (!open) // skip corners for open boundaries
      // corners are always interpolated (positive 'code')
      this->prescribe(b->getId()+offset, dof, abs(code));

  if (code <= 0) return; // If no projection, we're done

  // Build up the local edge-node correspondence for this edge
  for (LR::Basisfunction* b : edgeFunctions)
    de.MLGN.push_back(b->getId()+offset);

  // Get all elements connected to this edge
  std::vector<LR::Element*> edgeElements;
  de.lr->getEdgeElements(edgeElements,de.edg);

  // Build the MLGE and MNPC arrays
  de.MLGE.reserve(edgeElements.size());
  de.MNPC.reserve(edgeElements.size());
  for (LR::Element* el : edgeElements)
  {
    // for mixed FEM models, let MLGE point to the integration basis index
    if (de.lr != this->lrspline.get())
      de.MLGE.push_back(lrspline->getElementContaining(el->midpoint()));
    else
      de.MLGE.push_back(el->getId());

    IntVec mnpc; mnpc.reserve(el->support().size());
    for (LR::Basisfunction* b : el->support())
      mnpc.push_back(utl::findIndex(de.MLGN,b->getId()+offset));
    de.MNPC.emplace_back(std::move(mnpc));
  }

  dirich.push_back(de);
}


size_t ASMu2D::constrainEdgeLocal (int dir, bool open, int dof, int code,
                                   bool project)
{
  return 0; // TODO...
}


int ASMu2D::getCorner (int I, int J, int basis) const
{
  std::vector<LR::Basisfunction*> edgeFunctions;

  const LR::LRSplineSurface* srf = this->getBasis(basis);

  // Note: Corners are identified by "coordinates" {-1,-1} {-1,1} {1,-1} {1,1}.
  int dir = (I > 0 ? LR::EAST : LR::WEST) | (J > 0 ? LR::NORTH : LR::SOUTH);
  srf->getEdgeFunctions(edgeFunctions, static_cast<LR::parameterEdge>(dir));

  if (edgeFunctions.empty()) {
    std::cerr <<" *** ASMu2D::constrainCorner: Invalid corner I,J="
              << I <<","<< J << std::endl;
    return 0;
  }

  if (edgeFunctions.size() > 1)
    std::cerr <<"  ** ASMu2D::constrainCorner: "<< edgeFunctions.size()
              <<" corners returned from LRSplineSurface::getEdgeFunctions()"
              << std::endl;

  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  return edgeFunctions.front()->getId()+ofs;
}


void ASMu2D::constrainCorner (int I, int J, int dof, int code, char basis)
{
  if (basis < 1) basis = 1;

  int corner = this->getCorner(I,J,basis);
  if (corner > 0)
    this->prescribe(corner,dof,code);
}


void ASMu2D::constrainNode (double xi, double eta, int dof, int code)
{
  double xip[2] = { xi, eta };
  double param[2] = { 0.0, 0.0 };
  Vec3 X;

  // Check if the point has a matching node
  if (this->evalPoint(xip,param,X) >= 0)
    for (size_t inod = 1; inod <= nnod; inod++)
      if (this->getCoord(inod).equal(X))
      {
        IFEM::cout <<"\tPrescribing node "<< inod
                   <<" at X = "<< this->getCoord(inod) << std::endl;
        this->prescribe(inod,dof,code);
        return;
      }

  std::cerr <<"  ** ASMu2D::constrainNode: No matching node at u,v = "
            << param[0] <<","<< param[1] << std::endl;
}


#define DERR -999.99

double ASMu2D::getParametricLength (int iel, int dir) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > lrspline->nElements())
  {
    std::cerr <<" *** ASMu2D::getParametricLength: Element index "<< iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return DERR;
  }
#endif

  const LR::Element* el = lrspline->getElement(iel-1);
  switch (dir)
    {
    case 1: return el->umax() - el->umin();
    case 2: return el->vmax() - el->vmin();
    }

  std::cerr <<" *** ASMu2D::getParametricLength: Invalid edge direction "
            << dir << std::endl;
  return DERR;
}


/*!
  In case of NURBS, this is not strictly correct as we then return the rational
  spline coefficients rather than nodal (control point) coordinates.
  But since these are used as coefficients, and not coordinates at caller sites,
  we do this for simplicity.

  Consider introducing ASMbase::getCoefficients() to lessen confusion.
*/

bool ASMu2D::getCoordinates (Matrix& X, unsigned char nsd,
                             const LR::LRSplineSurface& spline,
                             int iel)
{
  const LR::Element* elm = iel > 0 ? spline.getElement(iel-1) : nullptr;
  X.resize(nsd, iel > 0 ? elm->nBasisFunctions() : spline.nBasisFunctions());

  // Lambda-function inserting coordinates for a given basis function
  // into the array X, accounting for the weights in case of NURBS
  auto&& insertCoords = [&X,iel](int n, size_t inod, LR::Basisfunction* b)
  {
    X.fillColumn(inod,&(*b->cp()));
    double weight = b->dim() == n+1 ? b->cp(n) : 1.0;
    if (weight <= 0.0)
    {
      std::cerr <<" *** ASMu2D::getCoordinates: Zero weight for node "<< inod;
      if (iel > 0) std::cerr <<" in element "<< iel;
      std::cerr << std::endl;
      return false;
    }
    else if (weight != 1.0)
      for (int j = 1; j <= n; j++)
        X(j,inod) /= weight;

    return true;
  };

  bool status = true;
  size_t inod = 0;
  if (iel > 0)
    for (LR::Basisfunction* b : elm->support())
      status &= insertCoords(nsd,++inod,b);
  else
    for (LR::Basisfunction* b : spline.getAllBasisfunctions())
      status &= insertCoords(nsd,++inod,b);

  return status;
}


bool ASMu2D::getElementCoordinates (Matrix& X, int iel, bool forceItg) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > lrspline->nElements())
  {
    std::cerr <<" *** ASMu2D::getElementCoordinates: Element index "<< iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return false;
  }
#endif

  bool stat;
  int basis = forceItg ? ASM::INTEGRATION_BASIS : ASM::GEOMETRY_BASIS;
  const LR::LRSplineSurface* spline = this->getBasis(basis);
  if (spline != lrspline.get()) {
    const LR::Element* el = lrspline->getElement(iel-1);
    stat = this->getCoordinates(X,nsd,*spline,
                                spline->getElementContaining(el->midpoint())+1);
  }
  else
    stat = this->getCoordinates(X,nsd,*spline,iel);

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return stat;
}


Vec3 ASMu2D::getCoord (size_t inod) const
{
  LR::Basisfunction* basis = lrspline->getBasisfunction(inod-1);
  if (!basis) {
    std::cerr << "Asked to get coordinate for node " << inod
              << ", but only have " << lrspline->nBasisFunctions() << std::endl;
    return Vec3();
  }
  return Vec3(&(*basis->cp()),nsd);
}


bool ASMu2D::updateCoords (const Vector& displ)
{
  std::cerr <<" *** ASMu2D::updateCoords: Not implemented!"<< std::endl;
  return false;
}



/*!
  \brief Static helper mapping an edge index into a \a parameterEdge enum value.
*/

static LR::parameterEdge getEdgeEnum (int edgeIndex)
{
  switch (edgeIndex) {
  case 1: return LR::WEST;
  case 2: return LR::EAST;
  case 3: return LR::SOUTH;
  case 4: return LR::NORTH;
  }

  return LR::NONE;
}


size_t ASMu2D::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (!lrspline || lIndex < 1 || lIndex > 4)
    return 0;
  else if (ldim < 1)
    return 1;

  std::vector<LR::Element*> edgeElms;
  lrspline->getEdgeElements(edgeElms, getEdgeEnum(lIndex));

  return edgeElms.size();
}


void ASMu2D::getGaussPointParameters (RealArray& uGP, int dir, int nGauss,
                                      int iel, const double* xi,
                                      const LR::LRSplineSurface* spline) const
{
  if (!spline)
    spline = lrspline.get();
  LR::getGaussPointParameters(spline, uGP, dir, nGauss, iel, xi);
}


double ASMu2D::getElementCorners (int iel, Vec3Vec& XC, RealArray* uC) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > lrspline->nElements())
  {
    std::cerr <<" *** ASMu2D::getElementCorners: Element index "<< iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return 0.0;
  }
#endif

  const LR::Element* el = lrspline->getElement(iel-1);
  double u[4] = { el->umin(), el->umax(), el->umin(), el->umax() };
  double v[4] = { el->vmin(), el->vmin(), el->vmax(), el->vmax() };

  XC.resize(4);
  if (uC) uC->resize(8,0.0);

  for (int i = 0; i < 4; i++)
  {
    double xi[2] = { u[i], v[i] };
    if (this->evalPoint(iel-1,xi,XC[i]) < 0)
      return 0.0;
    else if (uC)
    {
      uC->at(2*i)   = u[i];
      uC->at(2*i+1) = v[i];
    }
  }

  return getElementSize(XC);
}


void ASMu2D::getCornerPoints (int iel, PointVec& XC) const
{
  RealArray uC;
  Vec3Vec  XYZ;
  this->getElementCorners(iel,XYZ,&uC);

  XC.clear();
  XC.reserve(4);
  for (int i = 0; i < 4; i++)
    XC.push_back(utl::Point(XYZ[i], { uC[2*i], uC[2*i+1] }));
}


bool ASMu2D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu2D::integrate(I)");

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  bool use3rdDer = integrand.getIntegrandType() & Integrand::THIRD_DERIVATIVES;

  if (myCache.empty())
    myCache.emplace_back(std::make_unique<BasisFunctionCache>(*this));

  BasisFunctionCache& cache = *myCache.front();
  cache.setIntegrand(&integrand);
  if (!cache.init(use3rdDer ? 3 : (use2ndDer ? 2 : 1)))
    return false;

  const std::array<int,2>& ng = cache.nGauss();
  const std::array<const double*,2>& xg = cache.coord();
  const std::array<const double*,2>& wg = cache.weight();

  // Get the reduced integration quadrature points, if needed
  const double* xr = cache.coord(true)[0];
  const double* wr = cache.weight(true)[0];

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);

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
      Matrix   Xnod, Jac;
      Matrix3D Hess;
      double   dXidu[2];
      double   param[3] = { 0.0, 0.0, 0.0 };
      Vec4     X(param,time.t);

      const LR::Element* el = lrspline->getElement(iel-1);

      // Get element area in the parameter space
      double dA = 0.25*el->area();
      if (dA < 0.0)
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

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
        fe.h = this->getElementCorners(iel,fe.XC);

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
      {
        // Compute the element center
        param[0] = 0.5*(cache.getParam(0,iel-1,0)+cache.getParam(0,iel-1,ng[0]-1));
        param[1] = 0.5*(cache.getParam(1,iel-1,0)+cache.getParam(1,iel-1,ng[1]-1));
        if (this->evalPoint(iel-1,param,X) < 0)
          ok = false;
      }

      if (integrand.getIntegrandType() & Integrand::G_MATRIX)
      {
        // Element size in parametric space
        dXidu[0] = el->umax() - el->umin();
        dXidu[1] = el->vmax() - el->vmin();
      }

      size_t nen = el->support().size();
      if (integrand.getIntegrandType() & Integrand::AVERAGE)
      {
        // --- Compute average value of basis functions over the element -----

        fe.Navg.resize(nen,true);
        double area = 0.0;
        size_t ip = 0;
        for (int j = 0; j < ng[1]; j++)
          for (int i = 0; i < ng[0]; i++, ip++)
          {
            const BasisFunctionVals& bfs = cache.getVals(iel-1,ip);

            // Compute Jacobian determinant of coordinate mapping
            // and multiply by weight of current integration point
            double detJac = utl::Jacobian(Jac,fe.dNdX,Xnod,bfs.dNdu,false);
            double weight = dA*wg[0][i]*wg[1][j];

            // Numerical quadrature
            fe.Navg.add(bfs.N,detJac*weight);
            area += detJac*weight;
          }

        // Divide by element area
        fe.Navg /= area;
      }

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(nen,fe.iel);
      int nRed = cache.nGauss(true)[0];
      if (!integrand.initElement(MNPC[iel-1],fe,X,nRed*nRed,*A))
      {
        A->destruct();
        ok = false;
        continue;
      }

      if (integrand.getIntegrandType() & Integrand::UPDATED_NODES)
        if (!time.first || time.it > 0)
          if (!this->deformedConfig(Xnod,A->vec))
          {
            A->destruct();
            ok = false;
            continue;
          }

      if (xr)
      {
        // --- Selective reduced integration loop ------------------------------

        size_t ip = 0;
        for (int j = 0; j < nRed; j++)
          for (int i = 0; i < nRed; i++, ++ip)
          {
            // Local element coordinates of current integration point
            fe.xi  = xr[i];
            fe.eta = xr[j];

            // Parameter values of current integration point
            fe.u = param[0] = cache.getParam(0,iel-1,i,true);
            fe.v = param[1] = cache.getParam(1,iel-1,j,true);

            const BasisFunctionVals& bfs = cache.getVals(iel-1,ip,true);
            fe.N = bfs.N;

            // Compute Jacobian inverse and derivatives
            fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,bfs.dNdu);
            if (fe.detJxW == 0.0) continue; // skip singular points

            // Store tangent vectors in fe.G for shells
            if (nsd > 2) fe.G = Jac;

            // Cartesian coordinates of current integration point
            X.assign(Xnod * fe.N);

            // Compute the reduced integration terms of the integrand
            fe.detJxW *= dA*wr[i]*wr[j];
            if (!integrand.reducedInt(*A,fe,X))
              ok = false;
          }
      }


      // --- Integration loop over all Gauss points in each direction ----------

      size_t ip = 0;
      fe.iGP = firstIp + (iel-1)*ng[0]*ng[1]; // Global integration point counter

      for (int j = 0; j < ng[1]; j++)
        for (int i = 0; i < ng[0]; i++, fe.iGP++, ++ip)
        {
          // Local element coordinates of current integration point
          fe.xi  = xg[0][i];
          fe.eta = xg[1][j];

          // Parameter values of current integration point
          fe.u = param[0] = cache.getParam(0,iel-1,i);
          fe.v = param[1] = cache.getParam(1,iel-1,j);

          const BasisFunctionVals& bfs = cache.getVals(iel-1,ip);
          fe.N = bfs.N;
#if SP_DEBUG > 4
          if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
          {
            std::cout <<"\nBasis functions at a integration point "
                      <<" : (u,v) = "<< fe.u <<" "<< fe.v;
            for (size_t n = 1; n <= bfs.N.size(); n++)
              std::cout <<'\n'<< n <<'\t'<< bfs.N(n) <<'\t'
                        << bfs.dNdu(n,1) <<'\t'<< bfs.dNdu(n,2);
            std::cout << std::endl;
          }
#endif

          // Compute Jacobian inverse of coordinate mapping and derivatives
          fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,bfs.dNdu);
          if (fe.detJxW == 0.0) continue; // skip singular points

          // Compute Hessian of coordinate mapping and 2nd order derivatives
          if (use2ndDer)
          {
            if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,bfs.d2Ndu2,bfs.dNdu))
              ok = false;
            else if (nsd > 2)
              utl::Hessian(Hess,fe.H);
          }

          // Compute 3rd order derivatives
          if (use3rdDer)
            ok &= utl::Hessian2(fe.d3NdX3,Jac,bfs.d3Ndu3);

          // Compute G-matrix
          if (integrand.getIntegrandType() & Integrand::G_MATRIX)
            utl::getGmat(Jac,dXidu,fe.G);
          else if (nsd > 2)
            fe.G = Jac; // Store tangent vectors in fe.G for shells

#if SP_DEBUG > 4
          if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
            std::cout <<"\n"<< fe;
#endif

          // Cartesian coordinates of current integration point
          X.assign(Xnod * fe.N);

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= dA*wg[0][i]*wg[1][j];
#ifndef USE_OPENMP
          PROFILE3("Integrand::evalInt");
#endif
          if (!integrand.evalInt(*A,fe,time,X))
            ok = false;
        }

      // Finalize the element quantities
      if (ok && !integrand.finalizeElement(*A,fe,time,firstIp+(iel-1)*ng[0]*ng[1]))
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


bool ASMu2D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time,
                        const Real3DMat& itgPts)
{
  if (!lrspline) return true; // silently ignore empty patches

  if (integrand.getReducedIntegration(2) != 0)
  {
    std::cerr <<" *** ASMu2D::integrate(Integrand&,GlobalIntegral&,"
              <<"const TimeDomain&,const Real3DMat&): Available for standard"
              <<" integrands only."<< std::endl;
    return false;
  }

  PROFILE2("ASMu2D::integrate(I)");

  std::vector<size_t> MPitg(itgPts.size()+1,0);
  for (size_t i = MPitg.front() = 0; i < itgPts.size(); i++)
    MPitg[i+1] = MPitg[i] + itgPts[i].size();

  // Evaluate basis function values and derivatives at all integration points.
  // We do this before the integration point loop to exploit multi-threading
  // in the integrand evaluations, which may be the computational bottleneck.

  std::vector<Go::BasisDerivsSf>  spline1;
  std::vector<Go::BasisDerivsSf2> spline2;
  if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
    spline2.resize(MPitg.back());
  else
    spline1.resize(MPitg.back());

  size_t iel, ip, jp = 0;
  for (iel = 0; iel < itgPts.size(); iel++)
    for (ip = 0; ip < itgPts[iel].size(); ip++, jp++)
    {
      double u = itgPts[iel][ip][0];
      double v = itgPts[iel][ip][1];
      if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
        this->computeBasis(u,v,spline2[jp],iel);
      else
        this->computeBasis(u,v,spline1[jp],iel);
    }

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
      fe.p   = lrspline->order(0) - 1;
      fe.q   = lrspline->order(1) - 1;
      Matrix   dNdu, Xnod, Jac;
      Matrix3D d2Ndu2, Hess;
      Vec4     X(nullptr,time.t);

      const LR::Element* el = lrspline->getElement(iel-1);

      // Get element area in the parameter space
      double dA = 0.25*el->area();
      if (dA < 0.0)
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

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
        fe.h = this->getElementCorners(iel,fe.XC);

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
      {
        // Compute the element center
        fe.h = this->getElementCorners(iel,fe.XC);
        X.assign(0.25*(fe.XC[0]+fe.XC[1]+fe.XC[2]+fe.XC[3]));
      }

      // Initialize element quantities
      size_t nen = lrspline->getElement(iel-1)->support().size();
      LocalIntegral* A = integrand.getLocalIntegral(nen,fe.iel);
      if (!integrand.initElement(MNPC[iel-1],fe,X,0,*A))
      {
        A->destruct();
        ok = false;
        continue;
      }

      if (integrand.getIntegrandType() & Integrand::UPDATED_NODES)
        if (!time.first || time.it > 0)
          if (!this->deformedConfig(Xnod,A->vec))
          {
            A->destruct();
            ok = false;
            continue;
          }

      // --- Integration loop over all quadrature points in this element -------

      size_t jp = MPitg[iel-1]; // Patch-wise integration point counter
      fe.iGP = firstIp + jp;    // Global integration point counter

      const Real2DMat& elmPts = itgPts[iel-1]; // points for current element
      for (size_t ip = 0; ip < elmPts.size(); ip++, jp++, fe.iGP++)
      {
        // Parameter values of current integration point
        fe.u = elmPts[ip][0];
        fe.v = elmPts[ip][1];

        // Compute basis function derivatives at current integration point
        if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
          SplineUtils::extractBasis(spline2[jp],fe.N,dNdu,d2Ndu2);
        else
          SplineUtils::extractBasis(spline1[jp],fe.N,dNdu);

        // Compute Jacobian inverse of coordinate mapping and derivatives
        fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
        if (fe.detJxW == 0.0) continue; // skip singular points

        // Compute Hessian of coordinate mapping and 2nd order derivatives
        if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
        {
          if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
            ok = false;
          else if (nsd > 2)
            utl::Hessian(Hess,fe.H);
        }

        // Store tangent vectors in fe.G for shells
        if (nsd > 2) fe.G = Jac;

#if SP_DEBUG > 4
        if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
          std::cout <<"\n"<< fe;
#endif

        // Cartesian coordinates of current integration point
        X.assign(Xnod * fe.N);
        X.u = elmPts[ip].data();

        // Evaluate the integrand and accumulate element contributions
        fe.detJxW *= dA*elmPts[ip][2];
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

#if defined(SP_DEBUG) && !defined(USE_OPENMP)
      if (iel == -dbgElm)
        break; // Skipping all elements, except for -dbgElm
#endif
    }

  return ok;
}


bool ASMu2D::integrate (Integrand& integrand, int lIndex,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu2D::integrate(B)");

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);

  // Get Gaussian quadrature points and weights
  int nG1 = this->getNoGaussPt(lIndex%10 < 3 ? p2 : p1, true);
  int nGP = integrand.getBouIntegrationPoints(nG1);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex%10+1)/((lIndex%2) ? -2 : 2);

  const int t1 = abs(edgeDir);   // Tangent direction normal to the patch edge
  const int t2 = 3-abs(edgeDir); // Tangent direction along the patch edge

  std::array<Vector,2> gpar;
  for (int d = 0; d < 2; d++)
    if (-1-d == edgeDir)
    {
      gpar[d].resize(nGP);
      gpar[d].fill(d == 0 ? lrspline->startparam(0) : lrspline->startparam(1));
    }
    else if (1+d == edgeDir)
    {
      gpar[d].resize(nGP);
      gpar[d].fill(d == 0 ? lrspline->endparam(0) : lrspline->endparam(1));
    }

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  FiniteElement fe;
  fe.p  = p1 - 1;
  fe.q  = p2 - 1;
  fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
  double param[3] = { 0.0, 0.0, 0.0 };

  Matrix dNdu, Xnod, Jac;
  Vec4   X(param,time.t);
  Vec3   normal;


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 0;
  for (const LR::Element* el : lrspline->getAllElements())
  {
    if (!myElms.empty() && !glInt.threadSafe() &&
        std::find(myElms.begin(), myElms.end(), iel) == myElms.end()) {
        ++iel;
      continue;
    }

    fe.iel = MLGE[iel++];
#ifdef SP_DEBUG
    if (dbgElm < 0 && iel != -dbgElm)
      continue; // Skipping all elements, except for -dbgElm
#endif

    // Skip elements that are not on current boundary edge
    bool skipMe = false;
    switch (edgeDir)
    {
    case -1: if (el->umin() != lrspline->startparam(0)) skipMe = true; break;
    case  1: if (el->umax() != lrspline->endparam(0)  ) skipMe = true; break;
    case -2: if (el->vmin() != lrspline->startparam(1)) skipMe = true; break;
    case  2: if (el->vmax() != lrspline->endparam(1)  ) skipMe = true; break;
    }
    if (skipMe) continue;

    // Get element edge length in the parameter space
    double dS = 0.5*this->getParametricLength(iel,t2);
    if (dS < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel)) return false;

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      fe.h = this->getElementCorners(iel,fe.XC);

    // Initialize element quantities
    size_t nen = el->support().size();
    LocalIntegral* A = integrand.getLocalIntegral(nen,fe.iel,true);
    bool ok = integrand.initElementBou(MNPC[iel-1],*A);

    // Get integration gauss points over this element
    this->getGaussPointParameters(gpar[t2-1],t2-1,nGP,iel,xg);


    // --- Integration loop over all Gauss points along the edge ---------------

    fe.iGP = firstp; // Global integration point counter
    firstp += nGP;

    for (int i = 0; i < nGP && ok; i++, fe.iGP++)
    {
      // Local element coordinates and parameter values
      // of current integration point
      if (t1 == 2)
        fe.xi = xg[i];
      else
        fe.eta = xg[i];
      fe.u = param[0] = gpar[0][i];
      fe.v = param[1] = gpar[1][i];

      // Evaluate basis function derivatives at current integration points
      Go::BasisDerivsSf spline;
      this->computeBasis(fe.u, fe.v, spline, iel-1);

      // Fetch basis function derivatives at current integration point
      SplineUtils::extractBasis(spline,fe.N,dNdu);

      // Compute basis function derivatives and the edge normal
      fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
      if (fe.detJxW == 0.0) continue; // skip singular points

      if (edgeDir < 0) normal *= -1.0;

      // Store tangent vectors in fe.G for shells
      if (nsd > 2) fe.G = Jac;

#if SP_DEBUG > 4
      if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
        std::cout <<"\n"<< fe;
#endif

      // Cartesian coordinates of current integration point
      X.assign(Xnod * fe.N);

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


bool ASMu2D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time,
                        const ASM::InterfaceChecker& iChkgen)
{
  if (!geomB) return true; // silently ignore empty patches
  if (!(integrand.getIntegrandType() & Integrand::INTERFACE_TERMS)) return true;

  PROFILE2("ASMu2D::integrate(J)");

  const InterfaceChecker& iChk = static_cast<const InterfaceChecker&>(iChkgen);

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  FiniteElement fe;
  Matrix Xnod, Xnod2, dNdu, dN2du, Jac;
  double param[3] = { 0.0, 0.0, 0.0 };
  Vec4   X(param,time.t);
  Vec3   normal;

  int iel = 0;
  for (const LR::Element* elm : lrspline->getAllElements())
  {
    fe.iel = abs(MLGE[iel]);
    short int status = iChk.hasContribution(++iel);
    if (!status) continue; // no interface contributions for this element

    status &= iChk.elmBorderMask(elm->umin(),elm->umax(),
                                 elm->vmin(),elm->vmax());
    if (!status) continue; // no interface contributions for this element

#if SP_DEBUG > 3
    std::cout <<"\n\nIntegrating interface terms for element "<< fe.iel
              << std::endl;
#endif

    // Set up control point (nodal) coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel))
      return false;

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(MNPC[iel-1].size(),iel);
    bool ok = integrand.initElement(MNPC[iel-1],*A);

    // Loop over the element edges with contributions
    int bit = 8;
    for (int iedge = 4; iedge > 0 && ok; iedge--, bit /= 2)
      if (status & bit)
      {
        // Find the parametric direction of the edge normal {-2,-1, 1, 2}
        const int edgeDir = (iedge+1)/((iedge%2) ? -2 : 2);
        const int t1 = abs(edgeDir);   // Tangent direction normal to the edge
        const int t2 = 3-abs(edgeDir); // Tangent direction along the edge

        // Set up parameters
        double u1 = iedge != 2 ? elm->umin() : elm->umax();
        double v1 = iedge < 4  ? elm->vmin() : elm->vmax();
        double u2 = u1;
        double v2 = v1;

        for (double uv : iChk.getIntersections(iel,iedge))
        {
          if (iedge <= 2)
          {
            v1 = v2;
            v2 = uv;
          }
          else
          {
            u1 = u2;
            u2 = uv;
          }

          // Get element edge length in the parameter space
          double dS = 0.5*(iedge <= 2 ? v2 - v1 : u2 - u1);


          // --- Integration loop over all Gauss points along the edge ---------

          for (int g = 0; g < nGP && ok; g++)
          {
            // Local element coordinates and parameter values
            // of current integration point
            fe.xi  = t1 == 1 ? edgeDir : xg[g];
            fe.eta = t1 == 1 ? xg[g] : edgeDir/2;
            fe.u = param[0] = iedge <= 2 ? u1 : 0.5*((u2-u1)*xg[g] + u2 + u1);
            fe.v = param[1] = iedge >= 3 ? v1 : 0.5*((v2-v1)*xg[g] + v2 + v1);

            // Evaluate basis function derivatives at current integration points
            Go::BasisDerivsSf spline;
            this->computeBasis(fe.u, fe.v, spline, iel-1);
            SplineUtils::extractBasis(spline, fe.N, dNdu);

            // Compute Jacobian inverse of the coordinate mapping and
            // basis function derivatives w.r.t. Cartesian coordinates
            fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
            if (fe.detJxW == 0.0) continue; // skip singular points

            if (edgeDir < 0) normal *= -1.0;

            // Store tangent vectors in fe.G for shells
            if (nsd > 2) fe.G = std::move(Jac);

            // Cartesian coordinates of current integration point
            X.assign(Xnod * fe.N);

#if SP_DEBUG > 4
            std::cout <<"\n"<< fe;
#endif

            // Evaluate the integrand and accumulate element contributions
            fe.detJxW *= dS*wg[g];
            ok = integrand.evalInt(*A,fe,time,X,normal);
          }
        }
      }

    // Finalize the element quantities
    if (ok && !integrand.finalizeElement(*A,fe,time))
      ok = false;

    // Assembly of global system integral
    if (ok && !glInt.assemble(A,iel))
      ok = false;

    A->destruct();

    if (!ok) return false;
  }

  return true;
}


bool ASMu2D::diracPoint (Integrand& integrand, GlobalIntegral& glInt,
                         const double* param, const Vec3& pval)
{
  if (!lrspline) return false;

  int iel = lrspline->getElementContaining(param[0],param[1]);

  FiniteElement fe;
  fe.iel = MLGE[iel];
  fe.u   = param[0];
  fe.v   = param[1];
  if (!this->evaluateBasis(iel,fe))
    return false;

  LocalIntegral* A = integrand.getLocalIntegral(MNPC[iel].size(),fe.iel,true);
  bool ok = integrand.evalPoint(*A,fe,pval) && glInt.assemble(A,fe.iel);

  A->destruct();

  return ok;
}


/*!
  If \a param is NULL, then \a xi is taken to be in the parameter domain of
  the LR-spline object instead of the dimensionless [0,1] domain.
*/

int ASMu2D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  const LR::LRSplineSurface* geo = this->getBasis(ASM::GEOMETRY_BASIS);
  if (!geo) return -2;

  double u[2];
  if (param)
    for (int i = 0; i < 2; i++)
      u[i] = param[i] = (1.0-xi[i])*geo->startparam(i) + xi[i]*geo->endparam(i);
  else
  {
    u[0] = xi[0];
    u[1] = xi[1];
  }

  int iel;
#pragma omp critical
  iel = lrspline->getElementContaining(u[0],u[1]);

  return this->evalPoint(iel,u,X);
}


int ASMu2D::evalPoint (int iel, const double* param, Vec3& X) const
{
  const LR::LRSplineSurface* geo = this->getBasis(ASM::GEOMETRY_BASIS);

  int ielG = geo->getElementContaining(lrspline->getElement(iel)->midpoint());
  Go::BasisPtsSf bas;
  this->computeBasis(param[0],param[1],bas,ielG,geo);

  Matrix Xnod;
  if (!this->getElementCoordinates(Xnod,1+iel))
    return -1;

  X = Xnod * bas.basisValues;

  return 0;
}


int ASMu2D::findElementContaining (const double* param) const
{
  return lrspline ? 1 + lrspline->getElementContaining(param[0],param[1]) : -2;
}


double ASMu2D::findPoint (Vec3&, double*) const
{
  std::cerr <<" *** ASMu2D::findPoint: Not implemented yet"<< std::endl;
  return -1.0;
}


bool ASMu2D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
  if (!lrspline) return false;

  if (outputMaster)
    return outputMaster->getGridParameters(prm, dir, nSegPerSpan);

  double eps = ElementBlock::eps;

  for (const LR::Element* el : lrspline->getAllElements())
  {
    // Get parametric element evaluation points, optionally with some shrinkage
    std::pair<double,double> u;
    if (dir == 0)
      u = { el->umin(), el->umax() };
    else
      u = { el->vmin(), el->vmax() };

    double du = (u.second - u.first)*(1.0-2.0*eps)/nSegPerSpan;
    u.first = u.first*(1.0-eps) + u.second*eps;
    for (int j = 0; j <= nSegPerSpan; j++)
      for (int i = 0; i <= nSegPerSpan; i++)
        prm.push_back(u.first + du*double(dir == 0 ? i : j));
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

bool ASMu2D::tesselate (ElementBlock& grid, const int* npe) const
{
  if (!lrspline) return false;

  int nNodesPerElement =  npe[0]   * npe[1];
  int nSubElPerElement = (npe[0]-1)*(npe[1]-1);
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
    double du = (el->umax() - el->umin())*(1.0-2.0*eps)/(npe[0]-1);
    double dv = (el->vmax() - el->vmin())*(1.0-2.0*eps)/(npe[1]-1);
    for (int iv = 0; iv < npe[1]; iv++)
      for (int iu = 0; iu < npe[0]; iu++, inod++) {
        Vec3 Xpt;
        double U[2] = { umin + du*iu, vmin + dv*iv };
        this->evalPoint(iel, U, Xpt);
        grid.setCoor(inod, Xpt);
        grid.setParams(inod, U[0], U[1]);
      }
    ++iel;
  }

  int iStart = iel = inod = 0;
  for (int i = 0; i < nElements; i++, iStart += nNodesPerElement)
    for (int iv = 0; iv < npe[1]-1; iv++)
      for (int iu = 0; iu < npe[0]-1; iu++) {
        // enumerate nodes counterclockwise around the quad
        grid.setNode(inod++, iStart + (iv  )*npe[0] + (iu  ) );
        grid.setNode(inod++, iStart + (iv  )*npe[0] + (iu+1) );
        grid.setNode(inod++, iStart + (iv+1)*npe[0] + (iu+1) );
        grid.setNode(inod++, iStart + (iv+1)*npe[0] + (iu  ) );
        grid.setElmId(++iel, i+1);
      }

  return true;
}


bool ASMu2D::evalSolution (Matrix& sField, const Vector& locSol,
                           const int* npe, int n_f, bool piola) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  if (piola)
    return this->evalSolutionPiola(sField,locSol,gpar.data(),false);
  else
    return this->evalSolution(sField,locSol,gpar.data(),false,0,n_f);
}


bool ASMu2D::evalSolution (Matrix& sField, const Vector& locSol,
                           const RealArray* gpar, bool, int deriv, int) const
{
  PROFILE2("ASMu2D::evalSol(P)");

  size_t nComp = locSol.size() / this->getNoNodes(-1);
  if (nComp*this->getNoNodes(-1) > locSol.size())
    return false;

  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size())
    return false;

  FiniteElement fe;
  fe.p = lrspline->order(0) - 1;
  fe.q = lrspline->order(1) - 1;
  Vector   ptSol;
  Matrix   dNdu, dNdX, Jac, Xnod, eSol, ptDer;
  Matrix3D d2Ndu2, d2NdX2, Hess, ptDer2;

  Go::BasisDerivsSf2 spline2;
  int lel = -1;

  // Evaluate the primary solution field at each point
  sField.resize(nComp,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch element containing evaluation point.
    // Sadly, points are not always ordered in the same way as the elements.
    fe.u = gpar[0][i];
    fe.v = gpar[1][i];
    int iel = lrspline->getElementContaining(fe.u,fe.v);

    if (iel != lel && deriv == 2)
    {
      lel = iel; // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel+1))
        return false;
    }

    // Evaluate basis function values/derivatives at current parametric point
    // and multiply with control point values to get the point-wise solution
    int nskip = MNPC[iel].size() - lrspline->getElement(iel)->nBasisFunctions();
    utl::gather(MNPC[iel],nComp,locSol,eSol,0,nskip);
    switch (deriv)
    {
    case 0: // Evaluate the solution
      if (!this->evaluateBasis(iel,fe,deriv))
        return false;
      sField.fillColumn(1+i, eSol * fe.N);
      break;

    case 1: // Evaluate first derivatives of the solution
      if (!this->evaluateBasis(iel,fe,deriv))
        return false;
      ptDer.multiply(eSol,fe.dNdX);
      sField.fillColumn(1+i,ptDer);
      break;

    case 2: // Evaluate second derivatives of the solution
      this->computeBasis(fe.u,fe.v,spline2,iel);
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


bool ASMu2D::evalProjSolution (Matrix& sField, const Vector& locSol,
                               const int* npe, int n_f) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the projected solution at all sampling points
  if (!this->separateProjectionBasis())
    return this->evalSolution(sField,locSol,gpar.data(),false,0,n_f);

  // The projection uses a separate basis, need to interpolate
  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size())
    return false;

  Fields* f = this->getProjectedFields(locSol);
  if (!f) return false;

  // Evaluate the projected solution field at each point
  Vector vals;
  sField.resize(f->getNoFields(),nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    f->valueFE(ItgPoint(gpar[0][i],gpar[1][i]),vals);
    sField.fillColumn(1+i,vals);
  }

  delete f;
  return true;
}


bool ASMu2D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                           const int* npe, char project) const
{
  if (npe != nullptr && npe[0] != npe[1])
  {
    std::cerr <<" *** ASMu2D::evalSolution: LR B-splines require the"
              <<" same number of evaluation points in u- and v-direction."
              << std::endl;
    return false;
  }

  // Project the secondary solution onto the spline basis
  LR::LRSplineSurface* s = nullptr;
  if (project == 'S')
    s = this->scRecovery(integrand);
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
        Go::Point p;
        sField.resize(s->dimension(),gpar[0].size()*gpar[1].size());

        int iel = 0; // evaluation points are always structured in element order
        for (size_t i = 0; i < gpar[0].size(); i++)
        {
          if ((i+1)%npe[0] == 0) iel++;
          s->point(p,gpar[0][i],gpar[1][i],iel);
          sField.fillColumn(i+1,p.begin());
        }
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
    sField.resize(s->dimension(),s->nBasisFunctions());
    for (LR::Basisfunction* bf : s->getAllBasisfunctions())
      sField.fillColumn(bf->getId()+1,&(*bf->cp()));
    delete s;
    return true;
  }

  std::cerr <<" *** ASMu2D::evalSolution: Failure!"<< std::endl;
  return false;
}


bool ASMu2D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                           const RealArray* gpar, bool) const
{
  PROFILE2("ASMu2D::evalSol(S)");

  sField.resize(0,0);
  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size())
    return false;

  // TODO: investigate the possibility of doing "regular" refinement by
  //       uniform tesselation grid and ignoring LR mesh lines

  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  bool use3rdDer = integrand.getIntegrandType() & Integrand::THIRD_DERIVATIVES;

  FiniteElement fe(0,firstIp);
  fe.p = lrspline->order(0) - 1;
  fe.q = lrspline->order(1) - 1;
  Vector   solPt;
  Matrix   dNdu, Jac, Xnod;
  Matrix3D d2Ndu2, Hess;
  Matrix4D d3Ndu3;

  if (integrand.getIntegrandType() & Integrand::UPDATED_NODES)
  {
    // Calculate updated control point coordinates for the entire patch,
    // stored as the second primary solution vector in the integrand object
    // while the first vector being the current total displacement vector
    this->getNodalCoordinates(Xnod);
    Vectors& eV = const_cast<IntegrandBase&>(integrand).getSolutions();
    if (!this->deformedConfig(Xnod,eV,true))
      return false;
  }

  // Evaluate the secondary solution field at each point
  int lel = -1;
  for (size_t i = 0; i < nPoints; i++, fe.iGP++)
  {
    // Fetch element containing evaluation point
    // sadly, points are not always ordered in the same way as the elements
    fe.u = gpar[0][i];
    fe.v = gpar[1][i];
    int iel = lrspline->getElementContaining(fe.u,fe.v);

    // Evaluate the basis functions at current parametric point
    if (use3rdDer)
    {
      Go::BasisDerivsSf3 spline;
      this->computeBasis(fe.u,fe.v,spline,iel);
      SplineUtils::extractBasis(spline,fe.N,dNdu,d2Ndu2,d3Ndu3);
    }
    else if (use2ndDer)
    {
      Go::BasisDerivsSf2 spline;
      this->computeBasis(fe.u,fe.v,spline,iel);
      SplineUtils::extractBasis(spline,fe.N,dNdu,d2Ndu2);
    }
    else
    {
      Go::BasisDerivsSf spline;
      this->computeBasis(fe.u,fe.v,spline,iel);
      SplineUtils::extractBasis(spline,fe.N,dNdu);
    }

    if (iel != lel)
    {
      lel = iel; // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,iel+1))
        return false;
    }

    // Compute the Jacobian inverse
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (use2ndDer)
    {
      if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
        continue;
      else if (nsd > 2)
        utl::Hessian(Hess,fe.H);
    }

    // Compute 3rd order derivatives
    if (use3rdDer)
      utl::Hessian2(fe.d3NdX3,Jac,d3Ndu3);

    // Store tangent vectors in fe.G for shells
    if (nsd > 2) fe.G = std::move(Jac);

#if SP_DEBUG > 4
    std::cout <<"\n"<< fe;
#endif

    // Cartesian coordinates of current integration point
    utl::Point X4(Xnod*fe.N,{fe.u,fe.v});

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,X4,MNPC[iel]))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMu2D::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                               int, int orient, bool local) const
{
  if (basis == 0)
    basis = 1;

  const LR::LRSplineSurface* srf = this->getBasis(basis);
  if (!srf) return; // silently ignore empty patches

  std::vector<LR::Basisfunction*> edgeFunctions;
  if (lIndex > 0 && lIndex <= 4) {
    srf->getEdgeFunctions(edgeFunctions, getEdgeEnum(lIndex));
    if (orient >= 0) {
      int u = lIndex <= 2 ? 1 : 0;
      ASMLRSpline::Sort(u, 1-u, orient, edgeFunctions);
    }
  }

  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  for (LR::Basisfunction* b : edgeFunctions)
    nodes.push_back(local ? b->getId()+ofs : this->getNodeID(b->getId()+ofs));

#if SP_DEBUG > 1
  std::cout <<"Boundary nodes in patch "<< idx+1 <<" edge "<< lIndex <<":";
  for (int n : nodes) std::cout <<" "<< n;
  std::cout << std::endl;
#endif
}


bool ASMu2D::getOrder (int& p1, int& p2, int& p3) const
{
  p1 = lrspline->order(0);
  p2 = lrspline->order(1);
  p3 = 0;

  return true;
}


bool ASMu2D::updateDirichlet (const std::map<int,RealFunc*>& func,
                              const std::map<int,VecFunc*>& vfunc, double time,
                              const std::map<int,int>* g2l)
{
  std::map<int,RealFunc*>::const_iterator fit;
  std::map<int,VecFunc*>::const_iterator vfit;

  for (const DirichletEdge& dedg : dirich)
  {
    Real2DMat controlPts;
    if ((fit = func.find(dedg.code)) != func.end())
      this->edgeL2projection(dedg, *fit->second, controlPts, time);
    else if ((vfit = vfunc.find(dedg.code)) != vfunc.end())
      this->edgeL2projection(dedg, *vfit->second, controlPts, time);
    else
    {
      std::cerr <<" *** ASMu2D::updateDirichlet: Code "<< dedg.code
                <<" is not associated with any function."<< std::endl;
      return false;
    }
    if (controlPts.empty())
    {
      std::cerr <<" *** ASMu2D::updateDirichlet: Projection failure."
                << std::endl;
      return false;
    }

    // Loop over the (non-corner) nodes of this boundary curve
    for (size_t j = 0; j < dedg.MLGN.size(); j++)
      if (!dedg.isCorner(dedg.MLGN[j]))
        for (int dofs = dedg.dof; dofs > 0; dofs /= 10)
        {
          int dof = dofs%10;
          // Find the constraint equation for current (node, local DOF)
          MPC pDOF(MLGN[dedg.MLGN[j]-1],dof);
          MPCIter mit = mpcs.find(&pDOF);
          if (mit != mpcs.end())
          {
            // Now update the prescribed value in the constraint equation
            if (fit != func.end()) dof = 1; // scalar condition
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


size_t ASMu2D::getNoNodes (int basis) const
{
  if (basis == 0)
    return this->ASMbase::getNoNodes(basis);
  else
    return lrspline->nBasisFunctions();
}


size_t ASMu2D::getNoProjectionNodes () const
{
  return projB->nBasisFunctions();
}


bool ASMu2D::separateProjectionBasis () const
{
  return this->getBasis(ASM::PROJECTION_BASIS) != this->getBasis(1);
}


Field* ASMu2D::getProjectedField (const Vector& coefs) const
{
  if (coefs.size() == this->getNoProjectionNodes())
    return new LRSplineField2D(this->getBasis(ASM::PROJECTION_BASIS),
                               coefs,is_rational);

  std::cerr <<" *** ASMu2D::getProjectedFields: Non-matching coefficent array,"
            <<" size="<< coefs.size() <<" nnod="<< this->getNoProjectionNodes()
            << std::endl;
  return nullptr;
}


Fields* ASMu2D::getProjectedFields (const Vector& coefs, size_t) const
{
  size_t ncmp = coefs.size() / this->getNoProjectionNodes();
  if (ncmp*this->getNoProjectionNodes() == coefs.size())
    return new LRSplineFields2D(this->getBasis(ASM::PROJECTION_BASIS),
                                coefs,ncmp,is_rational);

  std::cerr <<" *** ASMu2D::getProjectedFields: Non-matching coefficent array,"
            <<" size="<< coefs.size() <<" nnod="<< this->getNoProjectionNodes()
            << std::endl;
  return nullptr;
}


bool ASMu2D::transferGaussPtVars (const LR::LRSpline* old_basis,
                                  const RealArray& oldVars, RealArray& newVars,
                                  int nGauss) const
{
  const LR::LRSplineSurface* newBasis = this->getBasis();
  const LR::LRSplineSurface* oldBasis = static_cast<const LR::LRSplineSurface*>(old_basis);

  size_t nGp = nGauss*nGauss;
  newVars.resize(newBasis->nElements()*nGp);

  const double* xi = GaussQuadrature::getCoord(nGauss);
  LagrangeInterpolator interp(Vector(xi,nGauss));

  int iEl = 0;
  for (const LR::Element* newEl : newBasis->getAllElements())
  {
    double u_center = 0.5*(newEl->umin() + newEl->umax());
    double v_center = 0.5*(newEl->vmin() + newEl->vmax());
    int iOld = oldBasis->getElementContaining(u_center,v_center);
    if (iOld < 0)
    {
      std::cerr <<" *** ASMu2D: Failed to locate element "<< newEl->getId()
                <<" of the new mesh in the old mesh."<< std::endl;
      return false;
    }
    const LR::Element* oldEl = oldBasis->getElement(iOld);

    std::array<Matrix,2> I;
    for (int i = 0; i < 2; i++)
    {
      RealArray UGP;
      LR::getGaussPointParameters(newBasis, UGP, i, nGauss, iEl+1, xi);
      double pmin = i == 0 ? oldEl->umin() : oldEl->vmin();
      double pmax = i == 0 ? oldEl->umax() : oldEl->vmax();
      for (double& p : UGP)
        p = -1.0 + 2.0*(p-pmin)/(pmax-pmin);

      // lagrangian interpolation
      I[i] = interp.get(UGP);
    }

    Matrix data(nGauss,nGauss);
    data.fill(&oldVars[nGp*iOld],nGp);

    Matrix newdata;
    newdata.multiply(I[0]*data,I[1],false,true);
    std::copy(newdata.ptr(), newdata.ptr()+nGp, newVars.begin()+iEl*nGp);

    ++iEl;
  }

  return true;
}


bool ASMu2D::transferGaussPtVarsN (const LR::LRSpline* old_basis,
                                   const RealArray& oldVars, RealArray& newVars,
                                   int nGauss) const
{
  const LR::LRSplineSurface* newBasis = this->getBasis();
  const LR::LRSplineSurface* oldBasis = static_cast<const LR::LRSplineSurface*>(old_basis);

  size_t nGP = nGauss*nGauss;
  newVars.clear();
  newVars.reserve(nGP*newBasis->nElements());

  struct Param{ double u,v; };
  std::vector<Param> oGP(nGP);

  const double* xi = GaussQuadrature::getCoord(nGauss);
  for (const LR::Element* newEl : newBasis->getAllElements())
  {
    double u_center = 0.5*(newEl->umin() + newEl->umax());
    double v_center = 0.5*(newEl->vmin() + newEl->vmax());
    int iOld = oldBasis->getElementContaining(u_center,v_center);
    if (iOld < 0)
    {
      std::cerr <<" *** ASMu2D: Failed to locate element "<< newEl->getId()
                <<" of the new mesh in the old mesh."<< std::endl;
      return false;
    }
    const LR::Element* oldEl = oldBasis->getElement(iOld);

    // find parameters of old gauss points
    double umin = oldEl->umin();
    double vmin = oldEl->vmin();
    double du = 0.5*(oldEl->umax() - umin);
    double dv = 0.5*(oldEl->vmax() - vmin);
    size_t k = 0;
    for (int j = 0; j < nGauss; ++j)
      for (int i = 0; i < nGauss; ++i, ++k) {
        oGP[k].u = umin + du * (xi[i] + 1.0);
        oGP[k].v = vmin + dv * (xi[j] + 1.0);
      }

    // parameters of new gauss points
    umin = newEl->umin();
    vmin = newEl->vmin();
    du = 0.5*(newEl->umax() - umin);
    dv = 0.5*(newEl->vmax() - vmin);
    for (int j = 0; j < nGauss; ++j)
      for (int i = 0; i < nGauss; ++i) {
        double u = umin + du * (xi[i] + 1.0);
        double v = vmin + dv * (xi[j] + 1.0);
        double dist = 1.0e16;
        size_t near = 0;
        for (size_t k = 0; k < nGP; ++k) {
          double nd = hypot(oGP[k].u-u,oGP[k].v-v);
          if (nd < dist) {
            near = k;
            dist = nd;
          }
        }
        newVars.push_back(oldVars[iOld*nGP+near]);
      }
  }

  return true;
}


bool ASMu2D::transferCntrlPtVars (const LR::LRSpline* old_basis,
                                  RealArray& newVars, int nGauss) const
{
  const LR::LRSplineSurface* newBasis = this->getBasis();
  const LR::LRSplineSurface* oldBasis = static_cast<const LR::LRSplineSurface*>(old_basis);

  newVars.clear();
  newVars.reserve(newBasis->nElements()*nGauss*nGauss*oldBasis->dimension());
  const double* xi = GaussQuadrature::getCoord(nGauss);

  for (int iel = 1; iel <= newBasis->nElements(); iel++)
  {
    RealArray U, V, ptVar;
    LR::getGaussPointParameters(newBasis, U, 0, nGauss, iel, xi);
    LR::getGaussPointParameters(newBasis, V, 1, nGauss, iel, xi);
    for (int j = 0; j < nGauss; j++)
      for (int i = 0; i < nGauss; i++)
      {
        oldBasis->point(ptVar,U[i],V[j],iel-1);
        newVars.insert(newVars.end(),ptVar.begin(),ptVar.end());
      }
  }

  return true;
}


void ASMu2D::generateThreadGroups (const Integrand& integrand, bool silence,
                                   bool ignoreGlobalLM)
{
  LR::generateThreadGroups(threadGroups, this->getBasis(1));
  if (this->separateProjectionBasis())
    LR::generateThreadGroups(projThreadGroups, projB.get());
  if (silence || threadGroups[0].size() < 2) return;

  IFEM::cout <<"\nMultiple threads are utilized during element assembly.";
#ifdef SP_DEBUG
  for (size_t i = 0; i < threadGroups[0].size(); i++)
    IFEM::cout <<"\n Color "<< i+1 <<": "
               << threadGroups[0][i].size() <<" elements";
  IFEM::cout << std::endl;
#else
  this->analyzeThreadGroups(threadGroups[0]);
#endif
}


void ASMu2D::changeNumThreads ()
{
  for (std::unique_ptr<BasisFunctionCache>& c : myCache)
    c->resizeThreadBuffers();
}


void ASMu2D::remapErrors (RealArray& errors,
                          const RealArray& origErr, bool elemErrors) const
{
  const LR::LRSplineSurface* basis = this->getBasis(1);

  if (elemErrors) {
    errors = origErr;
    return;
  }

  for (const LR::Element* elm : basis->getAllElements())
    for (const LR::Basisfunction* b : elm->support())
      errors[b->getId()] += origErr[elm->getId()];
}


bool ASMu2D::evaluate (const FunctionBase* func, RealArray& vec,
                       int, double time) const
{
  Matrix ctrlPvals;
  ASMu2D* pch = const_cast<ASMu2D*>(this);
  bool ok = pch->L2projection(ctrlPvals,const_cast<FunctionBase*>(func),time);
  vec = ctrlPvals;
  return ok;
}


ASMu2D::InterfaceChecker::InterfaceChecker (const ASMu2D& pch) : myPatch(pch)
{
  const double epsilon = 1.0e-6;
  const LR::LRSplineSurface* lr = myPatch.getBasis(1);

  for (LR::Meshline* m : lr->getAllMeshlines()) {
    RealArray isectpts(1,0.0);
    for (LR::Meshline* m2 : lr->getAllMeshlines())
      if (m->intersects(m2,&isectpts.back()))
        isectpts.push_back(0);

    isectpts.pop_back();
    std::sort(isectpts.begin(), isectpts.end());
    auto end = std::unique(isectpts.begin(), isectpts.end());
    isectpts.erase(end,isectpts.end());

    // find elements where this intersection lives
    RealArray parval_left(2), parval_right(2);
    for (size_t i = 0; i < isectpts.size()-1; i++) {
      if (m->is_spanning_u()) {
#if SP_DEBUG > 2
        std::cout << "Line piece from ("<< isectpts[i] << ", " << m->const_par_ << ") to ("
                  << isectpts[i+1] << ", " << m->const_par_ << ")" << std::endl;
#endif
        parval_left[0]  = (isectpts[i]+isectpts[i+1])/2.0;
        parval_right[0] = (isectpts[i]+isectpts[i+1])/2.0;
        parval_left[1]  = m->const_par_ - epsilon;
        parval_right[1] = m->const_par_ + epsilon;
      } else {
#if SP_DEBUG > 2
        std::cout << "Line piece from ("<< m->const_par_ << ", " << isectpts[i] << ") to ("
                  << m->const_par_ << ", " << isectpts[i+1] << ")" << std::endl;
#endif
        parval_left[0]  = m->const_par_ - epsilon;
        parval_right[0] = m->const_par_ + epsilon;
        parval_left[1]  = (isectpts[i]+isectpts[i+1])/2.0;
        parval_right[1] = (isectpts[i]+isectpts[i+1])/2.0;
      }
      int el1 = lr->getElementContaining(parval_left);
      int el2 = lr->getElementContaining(parval_right);
#if SP_DEBUG > 2
      std::cout << "\t elem1 " << el1 << " elem2 " << el2 << std::endl;
#endif
      if (m->is_spanning_u()) {
        if (el1 > -1 && el2 > -1) {
          intersections[el2*16 + 3].continuity =
          intersections[el1*16 + 4].continuity = lr->order(1)-m->multiplicity_-1;
          intersections[el2*16 + 3].pts.push_back(isectpts[i+1]);
          intersections[el1*16 + 4].pts.push_back(isectpts[i+1]);
        }
      } else {
        if (el1 > -1 && el2 > -1) {
          intersections[el2*16 + 1].continuity =
          intersections[el1*16 + 2].continuity = lr->order(0)-m->multiplicity_-1;
          intersections[el2*16 + 1].pts.push_back(isectpts[i+1]);
          intersections[el1*16 + 2].pts.push_back(isectpts[i+1]);
        }
      }
    }
  }

  for (std::pair<const int,Intersection>& it : intersections) {
    RealArray& points = it.second.pts;
    std::sort(points.begin(),points.end());
    auto end = std::unique(points.begin(),points.end());
    points.erase(end,points.end());
  }
}


short int ASMu2D::InterfaceChecker::hasContribution (int e, int, int, int) const
{
  const LR::Element* elm = myPatch.lrspline->getElement(e-1);

  bool neighbor[4];
  neighbor[0] = elm->getParmin(0) != myPatch.lrspline->startparam(0); // West
  neighbor[1] = elm->getParmax(0) != myPatch.lrspline->endparam(0);   // East
  neighbor[2] = elm->getParmin(1) != myPatch.lrspline->startparam(1); // South
  neighbor[3] = elm->getParmax(1) != myPatch.lrspline->endparam(1);   // North

  // Check for existing neighbors
  short int status = 0, s = 1;
  for (short int i = 0; i < 4; i++, s *= 2)
    if (neighbor[i]) status += s;

  return status;
}


const RealArray& ASMu2D::InterfaceChecker::getIntersections (int iel, int edge,
                                                             int* cont) const
{
  auto it = intersections.find((iel-1)*16 + edge);
  if (it == intersections.end())
  {
    static RealArray empty;
    return empty;
  }

  if (cont)
    *cont = it->second.continuity;

  return it->second.pts;
}


bool ASMu2D::refine (const LR::RefineData& prm, Vectors& sol)
{
  if (!this->ASMLRSpline::refine(prm,sol))
    return false;

  // check if refinement was actually done
  if (prm.elements.size() + prm.errors.size() == 0)
    return true;

  if (!this->separateProjectionBasis())
    return true;

  LR::LRSplineSurface* proj = this->getBasis(ASM::PROJECTION_BASIS);
  for (const LR::Meshline* line : lrspline->getAllMeshlines())
    if (line->span_u_line_)
      proj->insert_const_v_edge(line->const_par_,
                                line->start_, line->stop_,
                                line->multiplicity());
    else
      proj->insert_const_u_edge(line->const_par_,
                                line->start_, line->stop_,
                                line->multiplicity());

  proj->generateIDs();

  IFEM::cout <<"Refined projection basis: "<< proj->nElements()
             <<" elements "<< proj->nBasisFunctions() <<" nodes."
             << std::endl;

  return true;
}


void ASMu2D::generateBezierBasis ()
{
  bezier_u = this->getBezierBasis(lrspline->order(0));
  bezier_v = this->getBezierBasis(lrspline->order(1));
}


void ASMu2D::generateBezierExtraction ()
{
  PROFILE2("Bezier extraction");

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);

  myBezierExtract.resize(lrspline->nElements());
  RealArray extrMat;
  int iel = 0;
  for (const LR::Element* elm : lrspline->getAllElements())
  {
    // Get bezier extraction matrix
    lrspline->getBezierExtraction(iel,extrMat);
    myBezierExtract[iel].resize(elm->nBasisFunctions(),p1*p2);
    myBezierExtract[iel++].fill(extrMat.data(),extrMat.size());
  }
}


void ASMu2D::computeBasis (double u, double v, Go::BasisPtsSf& bas,
                           int iel, const LR::LRSplineSurface* spline) const
{
  if (!spline)
    spline = lrspline.get();

  if (is_rational) {
    this->computeBasisNurbs(u, v, bas, iel, *spline);
    return;
  }

  PROFILE3("ASMu2D::compBasis(0)");
  spline->computeBasis(u,v,bas,iel);
}


void ASMu2D::computeBasis (double u, double v, Go::BasisDerivsSf& bas,
                           int iel, const LR::LRSplineSurface* spline) const
{
  if (!spline)
    spline = lrspline.get();

  if (is_rational) {
    this->computeBasisNurbs(u, v, bas, iel, *spline);
    return;
  }

  PROFILE3("ASMu2D::compBasis(1)");
  spline->computeBasis(u,v,bas,iel);
}


void ASMu2D::computeBasis (double u, double v, Go::BasisDerivsSf2& bas,
                           int iel, const LR::LRSplineSurface* spline) const
{
  if (!spline)
    spline = lrspline.get();

  if (is_rational) {
    this->computeBasisNurbs(u, v, bas, iel, *spline);
    return;
  }

  PROFILE3("ASMu2D::compBasis(2)");
  spline->computeBasis(u,v,bas,iel);
}


void ASMu2D::computeBasis (double u, double v, Go::BasisDerivsSf3& bas,
                           int iel, const LR::LRSplineSurface* spline) const
{
  if (!spline)
    spline = lrspline.get();

  if (is_rational) {
    this->computeBasisNurbs(u, v, bas, iel, *spline);
    return;
  }

  PROFILE3("ASMu2D::compBasis(3)");
  spline->computeBasis(u,v,bas,iel);
}


void ASMu2D::getElmConnectivities (IntMat& neigh) const
{
  const double epsilon = 1.0e-6;
  const LR::LRSplineSurface* lr = this->getBasis(1);

  for (LR::Meshline* m : lr->getAllMeshlines()) {
    RealArray isectpts(1,0.0);
    for (LR::Meshline* m2 : lr->getAllMeshlines())
      if (m->intersects(m2,&isectpts.back()))
        isectpts.push_back(0);

    isectpts.pop_back();
    std::sort(isectpts.begin(), isectpts.end());
    auto end = std::unique(isectpts.begin(), isectpts.end());
    isectpts.erase(end,isectpts.end());

    // find elements where this intersection lives
    RealArray parval_left(2), parval_right(2);
    for (size_t i = 0; i < isectpts.size()-1; i++) {
      if (m->is_spanning_u()) {
        parval_left[0]  = (isectpts[i]+isectpts[i+1])/2.0;
        parval_right[0] = (isectpts[i]+isectpts[i+1])/2.0;
        parval_left[1]  = m->const_par_ - epsilon;
        parval_right[1] = m->const_par_ + epsilon;
      } else {
        parval_left[0]  = m->const_par_ - epsilon;
        parval_right[0] = m->const_par_ + epsilon;
        parval_left[1]  = (isectpts[i]+isectpts[i+1])/2.0;
        parval_right[1] = (isectpts[i]+isectpts[i+1])/2.0;
      }
      int el1 = lr->getElementContaining(parval_left);
      int el2 = lr->getElementContaining(parval_right);
      if (el1 > -1 && el2 > -1) {
        neigh[MLGE[el1]-1].push_back(MLGE[el2]-1);
        neigh[MLGE[el2]-1].push_back(MLGE[el1]-1);
      }
    }
  }
}


void ASMu2D::findBoundaryElms (IntVec& elms, int lIndex, int orient) const
{
  elms.clear();
  if (lIndex < 1 || lIndex > 4)
    return;

  std::vector<LR::Element*> elements;
  this->getBasis(1)->getEdgeElements(elements, getEdgeEnum(lIndex));

  if (orient >= 0)
    std::sort(elements.begin(), elements.end(),
              [orient,lIndex](LR::Element* a, LR::Element* b)
              {
                int index = lIndex < 3 ? 1 : 0;
                double am = a->midpoint()[index];
                double bm = b->midpoint()[index];
                return orient == 1 ? bm < am : am < bm;
              });

  for (const LR::Element* elem : elements)
    elms.push_back(elem->getId());
}


void ASMu2D::generateThreadGroupsFromElms (const IntVec& elms)
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


void ASMu2D::storeMesh (const std::string& fName, int fType) const
{
  if (fType%2)
  {
    std::ofstream meshFile("param_"+fName+".eps");
    lrspline->writePostscriptMesh(meshFile);
  }

  if ((fType/=2)%2)
  {
    std::ofstream meshFile("physical_"+fName+".eps");
    if (is_rational)
      const_cast<ASMu2D*>(this)->writePostscriptElementsNurbs(lrspline, meshFile);
    else
      lrspline->writePostscriptElements(meshFile);
  }

  if ((fType/=2)%2)
  {
    std::ofstream meshFile("param_dot_"+fName+".eps");
    lrspline->writePostscriptFunctionSpace(meshFile);
  }

  if ((fType / 2)%2)
  {
    std::ofstream meshFile("physical_dot_"+fName+".eps");
    if (is_rational)
      const_cast<ASMu2D*>(this)->writePostscriptMeshWithControlPointsNurbs(lrspline, meshFile);
    else
      lrspline->writePostscriptMeshWithControlPoints(meshFile);
  }
}


void ASMu2D::extractElmRes (const Matrix& globRes, Matrix& elmRes,
                            size_t internalFirst) const
{
  if (outputMaster)
    outputMaster->extractElmRes(globRes, elmRes, internalFirst);
  else
    this->ASMbase::extractElmRes(globRes, elmRes, internalFirst);
}


ASMu2D::BasisFunctionCache::BasisFunctionCache (const ASMu2D& pch) :
  ::BasisFunctionCache<2>(), patch(pch)
{
}


ASMu2D::BasisFunctionCache::BasisFunctionCache (const BasisFunctionCache& cache,
                                                int b) :
  ::BasisFunctionCache<2>(cache), patch(cache.patch)
{
  basis = b;
}


bool ASMu2D::BasisFunctionCache::internalInit ()
{
  if (!mainQ->xg[0])
    this->setupQuadrature();

  nTotal = patch.nel*mainQ->ng[0]*mainQ->ng[1];
  if (reducedQ->xg[0])
    nTotalRed = patch.nel*reducedQ->ng[0]*reducedQ->ng[1];

  return true;
}


void ASMu2D::BasisFunctionCache::internalCleanup ()
{
  if (basis == 1) {
    mainQ->reset();
    reducedQ->reset();
  }
}


bool ASMu2D::BasisFunctionCache::setupQuadrature ()
{
  const int p1 = patch.lrspline->order(0);
  const int p2 = patch.lrspline->order(1);

  // Get Gaussian quadrature points and weights
  for (int d = 0; d < 2; d++)
  {
    mainQ->ng[d] = patch.getNoGaussPt(d == 0 ? p1 : p2);
    mainQ->xg[d] = GaussQuadrature::getCoord(mainQ->ng[d]);
    mainQ->wg[d] = GaussQuadrature::getWeight(mainQ->ng[d]);
    if (!mainQ->xg[d] || !mainQ->wg[d]) return false;
  }

  // Get the reduced integration quadrature points, if needed
  int nRed = integrand ? integrand->getReducedIntegration(mainQ->ng[0]) : 0;
  if (nRed > 0)
  {
    reducedQ->xg[0] = reducedQ->xg[1] = GaussQuadrature::getCoord(nRed);
    reducedQ->wg[0] = reducedQ->wg[1] = GaussQuadrature::getWeight(nRed);
    if (!reducedQ->xg[0] || !reducedQ->wg[0]) return false;
  } else if (nRed < 0)
    nRed = mainQ->ng[0]; // The integrand needs to know nGauss

  reducedQ->ng[0] = reducedQ->ng[1] = nRed;

  // Compute parameter values of the Gauss points over the whole patch
  mainQ->gpar[0].resize(mainQ->ng[0],patch.nel);
  mainQ->gpar[1].resize(mainQ->ng[1],patch.nel);
  if (reducedQ->xg[0]) {
    reducedQ->gpar[0].resize(reducedQ->ng[0],patch.nel);
    reducedQ->gpar[1].resize(reducedQ->ng[1],patch.nel);
  }
  for (size_t iel = 1; iel <= patch.nel; ++iel)
  {
    RealArray u, v;
    patch.getGaussPointParameters(u,0,mainQ->ng[0],iel,mainQ->xg[0]);
    patch.getGaussPointParameters(v,1,mainQ->ng[1],iel,mainQ->xg[1]);
    mainQ->gpar[0].fillColumn(iel,u.data());
    mainQ->gpar[1].fillColumn(iel,v.data());

    if (reducedQ->xg[0])
    {
      patch.getGaussPointParameters(u,0,reducedQ->ng[0],iel,reducedQ->xg[0]);
      patch.getGaussPointParameters(v,1,reducedQ->ng[1],iel,reducedQ->xg[0]);
      reducedQ->gpar[0].fillColumn(iel,u.data());
      reducedQ->gpar[1].fillColumn(iel,v.data());
    }
  }
  return true;
}


BasisFunctionVals ASMu2D::BasisFunctionCache::calculatePt (size_t el,
                                                           size_t gp,
                                                           bool reduced) const
{
  PROFILE2("Spline evaluation");
  const std::array<size_t,2> gpIdx = this->gpIndex(gp,reduced);
  double u = this->getParam(0,el,gpIdx[0],reduced);
  double v = this->getParam(1,el,gpIdx[1],reduced);

  if (patch.lrspline.get() != patch.getBasis(basis)) {
    const LR::Element* el1 = patch.lrspline->getElement(el);
    el = patch.getBasis(basis)->getElementContaining(el1->midpoint());
  }

  BasisFunctionVals result;
  if (nderiv == 1 || reduced) {
    Go::BasisDerivsSf spline;
    patch.computeBasis(u,v,spline,el,patch.getBasis(basis));
    SplineUtils::extractBasis(spline,result.N,result.dNdu);
  } else if (nderiv == 2) {
    Go::BasisDerivsSf2 spline;
    patch.computeBasis(u,v,spline,el,patch.getBasis(basis));
    SplineUtils::extractBasis(spline,result.N,result.dNdu,result.d2Ndu2);
  } else if (nderiv == 3) {
    Go::BasisDerivsSf3 spline;
    patch.computeBasis(u,v,spline,el,patch.getBasis(basis));
    SplineUtils::extractBasis(spline,result.N,result.dNdu,result.d2Ndu2,result.d3Ndu3);
  }

  return result;
}


void ASMu2D::BasisFunctionCache::calculateAll ()
{
  PROFILE2("Spline evaluation");
  // Evaluate basis function values and derivatives at all integration points.
  // We do this before the integration point loop to exploit multi-threading
  // in the integrand evaluations, which may be the computational bottleneck.
  size_t iel, jp, rp;
  for (iel = jp = rp = 0; iel < patch.nel; iel++)
  {
    for (int j = 0; j < mainQ->ng[1]; j++)
      for (int i = 0; i < mainQ->ng[0]; i++, jp++)
        values[jp] = this->calculatePt(iel,j*mainQ->ng[0]+i,false);

    if (reducedQ->xg[0])
      for (int j = 0; j < reducedQ->ng[1]; j++)
        for (int i = 0; i < reducedQ->ng[0]; i++, rp++)
          valuesRed[rp] = this->calculatePt(iel,j*reducedQ->ng[0]+i,true);
  }
}
