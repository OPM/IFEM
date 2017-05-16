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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"
#include "LRSpline/Basisfunction.h"

#include "ASMu2D.h"
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
#include "Function.h"
#include "Vec3Oper.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include <array>
#include <fstream>


ASMu2D::ASMu2D (unsigned char n_s, unsigned char n_f)
  : ASMunstruct(2,n_s,n_f), lrspline(nullptr), tensorspline(nullptr),
    bezierExtract(myBezierExtract)
{
}


ASMu2D::ASMu2D (const ASMu2D& patch, unsigned char n_f)
  : ASMunstruct(patch,n_f), lrspline(patch.lrspline), tensorspline(nullptr),
    bezierExtract(patch.myBezierExtract)
{
  // Need to set nnod here,
  // as hasXNodes might be invoked before the FE data is generated
  if (nnod == 0 && lrspline)
    nnod = lrspline->nBasisFunctions();
}


bool ASMu2D::read (std::istream& is)
{
  if (shareFE) return false;

  // read inputfile as either an LRSpline file directly
  // or a tensor product B-spline and convert
  char firstline[256];
  is.getline(firstline, 256);
  if (strncmp(firstline, "# LRSPLINE", 10) == 0) {
    lrspline.reset(new LR::LRSplineSurface());
    is >> *lrspline;
  } else { // probably a SplineSurface, so we'll read that and convert
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

  if (!is.good() && !is.eof())
  {
    std::cerr <<" *** ASMu2D::read: Failure reading spline data"<< std::endl;
    lrspline.reset();
    return false;
  }
  else if (lrspline->dimension() < 2)
  {
    std::cerr <<" *** ASMu2D::read: Invalid spline lrsplineace patch, dim="
        << lrspline->dimension() << std::endl;
    lrspline.reset();
    return false;
  }
  else if (lrspline->dimension() < nsd)
  {
    std::cout <<"  ** ASMu2D::read: The dimension of this lrsplineace patch "
        << lrspline->dimension() <<" is less than nsd="<< nsd
        <<".\n                   Resetting nsd to "<< lrspline->dimension()
        <<" for this patch."<< std::endl;
    nsd = lrspline->dimension();
  }

  geo = lrspline.get();
  return true;
}


bool ASMu2D::write (std::ostream& os, int) const
{
  if (!lrspline) return false;

  os << *lrspline;

  return os.good();
}


void ASMu2D::clear (bool retainGeometry)
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
  double ucurr, uprev = *(uit++);
  while (uit != tensorspline->basis(dir).end())
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
    tensorspline->insertKnot_u(extraKnots);
  else
    tensorspline->insertKnot_v(extraKnots);

  lrspline.reset(new LR::LRSplineSurface(tensorspline));
  geo = lrspline.get();

  return true;
}

bool ASMu2D::refine (int dir, const RealArray& xi)
{
  if (!tensorspline || dir < 0 || dir > 1 || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = tensorspline->basis(dir).begin();
  double ucurr, uprev = *(uit++);
  while (uit != tensorspline->basis(dir).end())
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
    tensorspline->insertKnot_u(extraKnots);
  else
    tensorspline->insertKnot_v(extraKnots);

  lrspline.reset(new LR::LRSplineSurface(tensorspline));
  geo = lrspline.get();

  return true;
}


/*!
  Refines all elements for which refC(X0) < refTol,
  where X0 is the element center.
*/

bool ASMu2D::refine (const RealFunc& refC, double refTol)
{
  Go::Point X0;
  std::vector<int> elements;
  std::vector<LR::Element*>::const_iterator eit = lrspline->elementBegin();
  for (int iel = 0; eit != lrspline->elementEnd(); iel++, ++eit)
  {
    double u0 = 0.5*((*eit)->umin() + (*eit)->umax());
    double v0 = 0.5*((*eit)->vmin() + (*eit)->vmax());
    lrspline->point(X0,u0,v0);
    if (refC(SplineUtils::toVec3(X0,nsd)) < refTol)
      elements.push_back(iel);
  }

  Vectors dummySol;
  LR::RefineData prm(true);
  prm.options = { 10, 1, 2 };
  prm.elements = this->getFunctionsForElements(elements);
  return this->refine(prm,dummySol);
}


bool ASMu2D::raiseOrder (int ru, int rv)
{
  if (!tensorspline) return false;
  if (shareFE) return true;

  tensorspline->raiseOrder(ru,rv);
  lrspline.reset(new LR::LRSplineSurface(tensorspline));
  geo = lrspline.get();
  return true;
}


bool ASMu2D::evaluateBasis (FiniteElement& fe, int derivs) const
{
  LR::Element* el = lrspline->getElement(fe.iel-1);
  if (!el) return false;

  fe.xi  = 2.0*(fe.u - el->umin()) / (el->umax() - el->umin()) - 1.0;
  fe.eta = 2.0*(fe.v - el->vmin()) / (el->vmax() - el->vmin()) - 1.0;
  RealArray Nu = bezier_u.computeBasisValues(fe.xi, derivs);
  RealArray Nv = bezier_v.computeBasisValues(fe.eta,derivs);
  const Matrix& C = bezierExtract[fe.iel-1];

  if (derivs < 1) {
    Matrix B;
    B.outer_product(Nu,Nv);
    fe.N = C*static_cast<const Vector&>(B);

#if SP_DEBUG > 2
    if (fabs(fe.N.sum()-1.0) > 1.0e-10)
      std::cerr <<"fe.N do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(static_cast<const Vector&>(B).sum()-1.0) > 1.0e-10)
      std::cerr <<"Bezier basis do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else
      return true; // The basis is OK

    return false;
#endif
  }
  else {
    int p = lrspline->order(0)*lrspline->order(1);

    Vector B(p);
    Vector Bu(p); // Bezier basis functions differentiated wrt u
    Vector Bv(p); // Bezier basis functions differentiated wrt v

    size_t i, j, k = 0;
    for (j = 0; j < Nv.size(); j+=(derivs+1))
      for (i = 0; i < Nu.size(); i+=(derivs+1), k++) {
        B[k]  = Nu[i  ]*Nv[j  ];
        Bu[k] = Nu[i+1]*Nv[j  ];
        Bv[k] = Nu[i  ]*Nv[j+1];
      }

    fe.N = C*B;

    Matrix dNdu(el->nBasisFunctions(),2);
    dNdu.fillColumn(1,C*Bu);
    dNdu.fillColumn(2,C*Bv);
    Matrix Xnod, Jac;
    this->getElementCoordinates(Xnod,fe.iel);
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

#if SP_DEBUG > 2
    if (fabs(fe.N.sum()-1.0) > 1.0e-10)
      std::cerr <<"fe.N do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(B.sum()-1.0) > 1.0e-10)
      std::cerr <<"Bezier basis do not sum to one at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(dNdu.getColumn(1).sum()) > 1.0e-10)
      std::cerr <<"dNdu not sums to zero at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(dNdu.getColumn(1).sum()) > 1.0e-10)
      std::cerr <<"dNdv not sums to zero at integration point #"
                << fe.iGP << std::endl;
    else if (fabs(Bu.sum()) > 1.0e-10 || fabs(Bv.sum()) > 1.0e-10)
      std::cerr <<"Bezier derivatives do not sum to zero at integration point #"
                << fe.iGP << std::endl;
    else
      return true; // The basis is OK

    return false;
#endif
  }
  return true;
}


bool ASMu2D::generateFEMTopology ()
{
  // At this point we are through with the tensor spline object,
  // so release it to avoid memory leakage
  delete tensorspline;
  tensorspline = nullptr;

  if (!lrspline) return false;

  nnod = lrspline->nBasisFunctions();
  nel  = lrspline->nElements();

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);

  bezier_u = getBezierBasis(p1);
  bezier_v = getBezierBasis(p2);

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

  myBezierExtract.resize(nel);
  lrspline->generateIDs();

  RealArray extrMat;
  std::vector<LR::Element*>::const_iterator eit = lrspline->elementBegin();
  for (size_t iel = 0; iel < nel; iel++, ++eit)
  {
    myMLGE[iel] = ++gEl; // global element number over all patches
    myMNPC[iel].resize((*eit)->nBasisFunctions());

    int lnod = 0;
    for (LR::Basisfunction *b : (*eit)->support())
      myMNPC[iel][lnod++] = b->getId();

    {
      PROFILE("Bezier extraction");

      // Get bezier extraction matrix
      lrspline->getBezierExtraction(iel,extrMat);
      myBezierExtract[iel].resize((*eit)->nBasisFunctions(),p1*p2);
      myBezierExtract[iel].fill(extrMat.data(),extrMat.size());
    }
  }

  for (size_t inod = 0; inod < nnod; inod++)
    myMLGN[inod] = ++gNod;

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
  this->getBoundaryNodes(edge, slaveNodes, basis, thick, true);
  for (int& it : slaveNodes)
    it += slave;

  // Set up the master node numbers for the neighboring surface patch
  IntVec masterNodes;
  neighbor.getBoundaryNodes(nedge, masterNodes, basis, thick, true);
  for (int& it : masterNodes)
    it += master;

  if (masterNodes.size() != slaveNodes.size())
  {
    std::cerr <<" *** ASMu2D::connectBasis: Non-matching edges, sizes "
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
      std::cerr <<" *** ASMu2D::connectBasis: Non-matching nodes "
                << node <<": "<< neighbor.getCoord(node)
                <<"\n                                          and "
                << slave <<": "<< this->getCoord(slave) << std::endl;
      return false;
    }
  }

  return true;
}


/*
void ASMu2D::closeEdges (int dir, int basis, int master)
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
*/


std::vector<int> ASMu2D::getEdgeNodes (int edge, int basis) const
{
  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  std::vector<LR::Basisfunction*> edgeFunctions;
  this->getBasis(basis)->getEdgeFunctions(edgeFunctions,
                                          static_cast<LR::parameterEdge>(edge));

  ASMunstruct::Sort(edgeFunctions);
  std::vector<int> result(edgeFunctions.size());
  std::transform(edgeFunctions.begin(), edgeFunctions.end(), result.begin(),
                 [ofs](LR::Basisfunction* a) { return a->getId()+ofs; });

  return result;
}


void ASMu2D::constrainEdge (int dir, bool open, int dof, int code, char basis)
{
  // figure out function index offset (when using multiple basis)
  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  // figure out what edge we are at
  LR::parameterEdge edge;
  switch (dir) {
  case -2: edge = LR::SOUTH; break;
  case -1: edge = LR::WEST;  break;
  case  1: edge = LR::EAST;  break;
  case  2: edge = LR::NORTH; break;
  default: return;
  }

  // fetch the right basis to consider
  LR::LRSplineSurface* lr = this->getBasis(basis);

  // get all elements and functions on this edge
  std::vector<LR::Basisfunction*> edgeFunctions;
  std::vector<LR::Element*>       edgeElements;
  lr->getEdgeFunctions(edgeFunctions,edge);
  lr->getEdgeElements (edgeElements ,edge);

  // find the corners since these are not to be included in the L2-fitting
  // of the inhomogenuous dirichlet boundaries; corners are interpolatory.
  // Optimization note: loop over the "edge"-container to manually pick up
  // the end nodes. LRspine::getEdgeFunctions() does a global search.
  std::vector<LR::Basisfunction*> c1, c2;
  switch (edge)
  {
  case LR::SOUTH:
    lr->getEdgeFunctions(c1, LR::SOUTH_WEST);
    lr->getEdgeFunctions(c2, LR::SOUTH_EAST);
    break;
  case LR::WEST:
    lr->getEdgeFunctions(c1, LR::SOUTH_WEST);
    lr->getEdgeFunctions(c2, LR::NORTH_WEST);
    break;
  case LR::EAST:
    lr->getEdgeFunctions(c1, LR::NORTH_EAST);
    lr->getEdgeFunctions(c2, LR::SOUTH_EAST);
    break;
  case LR::NORTH:
    lr->getEdgeFunctions(c1, LR::NORTH_WEST);
    lr->getEdgeFunctions(c2, LR::NORTH_EAST);
    break;
  default: return;
  }

  // build up the local element/node correspondence needed by the projection
  // call on this edge by ASMu2D::updateDirichlet()
  DirichletEdge de(edgeFunctions.size(), edgeElements.size(), dof, code, basis);
  de.corners[0] = c1[0]->getId();
  de.corners[1] = c2[0]->getId();
  de.edg  = edge;
  de.lr   = lr;
  int bcode = abs(code);

  int j = 0;
  for (auto b : edgeFunctions)
  {
    de.MLGN[j++] = b->getId();
    // skip corners for open boundaries
    if (open && (b->getId() == de.corners[0] || b->getId() == de.corners[1]))
      continue;
    else
    {
      // corners are interpolated (positive 'code')
      if (b->getId() == de.corners[0] || b->getId() == de.corners[1])
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
  for (size_t i=0; i<edgeElements.size(); i++)
  {
    LR::Element* el = edgeElements[i];

    // for mixed FEM models, let MLGE point to the *geometry* index
    if (de.lr != this->lrspline.get())
    {
      double umid = (el->umax() + el->umin())/2.0;
      double vmid = (el->vmax() + el->vmin())/2.0;
      de.MLGE[i] = lrspline->getElementContaining(umid, vmid);
    }
    else
    {
      de.MLGE[i] = el->getId();
    }
    for (auto b : el->support())
    {
      de.MNPC[i].push_back(-1);
      for (auto l2g : de.MLGN)
        // can't do map::find, since this works on first
        if (b->getId() == l2g.second)
          de.MNPC[i].back() = l2g.first;
    }
  }
  if (code > 0)
    dirich.push_back(de);
}


size_t ASMu2D::constrainEdgeLocal (int dir, bool open, int dof, int code,
                                   bool project)
{
  return 0; // TODO...
}


int ASMu2D::getCorner(int I, int J, int basis) const
{
  std::vector<LR::Basisfunction*> edgeFunctions;

  const LR::LRSplineSurface* srf = this->getBasis(basis);

  // Note: Corners are identified by "coordinates" {-1,-1} {-1,1} {1,-1} {1,1}.
  if (I < 0) {
    if (J < 0)
      srf->getEdgeFunctions(edgeFunctions, LR::SOUTH_WEST);
    else if (J > 0)
      srf->getEdgeFunctions(edgeFunctions, LR::NORTH_WEST);
  }
  else if (I > 0) {
    if (J < 0)
      srf->getEdgeFunctions(edgeFunctions, LR::SOUTH_EAST);
    else if (J > 0)
      srf->getEdgeFunctions(edgeFunctions, LR::NORTH_EAST);
  }

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
  int corner = getCorner(I, J, basis);
  if (corner == 0)
    return;
  else
    this->prescribe(corner,dof,code);
}


// Hopefully we don't have to constrain non-corner single nodes inside patches.
// KMO: Actually, we would like to have this, to prescribe mid-edge points, etc.
// Can it be done, Kjetil?
void ASMu2D::constrainNode (double xi, double eta, int dof, int code)
{
  std::cerr <<" *** ASMu2D::constrainNode: Not implemented yet!"<< std::endl;
  /*
  if (xi  < 0.0 || xi  > 1.0) return;
  if (eta < 0.0 || eta > 1.0) return;

  int n1, n2;
  if (!this->getSize(n1,n2,1)) return;

  int node = 1;
  if (xi  > 0.0) node += int(0.5+(n1-1)*xi);
  if (eta > 0.0) node += n1*int(0.5+(n2-1)*eta);

  this->prescribe(node,dof,code);
  */
}


#define DERR -999.99

double ASMu2D::getParametricArea (int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > lrspline->nElements())
  {
    std::cerr <<" *** ASMu2D::getParametricArea: Element index "<< iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return DERR;
  }
#endif

  return lrspline->getElement(iel-1)->area();
}


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

  LR::Element* el = lrspline->getElement(iel-1);
  switch (dir)
  {
  case 1: return el->vmax() - el->vmin();
  case 2: return el->umax() - el->umin();
  }

  std::cerr <<" *** ASMu2D::getParametricLength: Invalid edge direction "
            << dir << std::endl;
  return DERR;
}


bool ASMu2D::getElementCoordinates (Matrix& X, int iel) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > lrspline->nElements())
  {
    std::cerr <<" *** ASMu2D::getElementCoordinates: Element index "<< iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return false;
  }
#endif

  LR::Element* el = lrspline->getElement(iel-1);
  X.resize(nsd,el->nBasisFunctions());

  int n = 1;
  for (LR::Basisfunction* b : el->support())
    X.fillColumn(n++,&(*b->cp()));

#if SP_DEBUG > 2
  std::cout <<"\nCoordinates for element "<< iel << X << std::endl;
#endif
  return true;
}


void ASMu2D::getNodalCoordinates (Matrix& X) const
{
  X.resize(nsd,lrspline->nBasisFunctions());

  int inod = 1;
  for (LR::Basisfunction* b : lrspline->getAllBasisfunctions())
    X.fillColumn(inod++,&(*b->cp()));
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


size_t ASMu2D::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (!lrspline)
    return 0;
  else if (ldim < 1 && lIndex > 0)
    return 1;

  LR::parameterEdge edge;
  switch(lIndex)
  {
  case 1: edge = LR::WEST;  break;
  case 2: edge = LR::EAST;  break;
  case 3: edge = LR::SOUTH; break;
  case 4: edge = LR::NORTH; break;
  default:edge = LR::NONE;
  }

  std::vector<LR::Element*> edgeElms;
  lrspline->getEdgeElements(edgeElms, edge);

  return edgeElms.size();
}


void ASMu2D::getGaussPointParameters (RealArray& uGP, int dir, int nGauss,
                                      int iel, const double* xi) const
{
  LR::getGaussPointParameters(lrspline.get(), uGP, dir, nGauss, iel, xi);
}


double ASMu2D::getElementCorners (int iel, Vec3Vec& XC) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || iel > lrspline->nElements())
  {
    std::cerr <<" *** ASMu2D::getElementCorners: Element index "<< iel
              <<" out of range [1,"<< lrspline->nElements() <<"]."<< std::endl;
    return 0.0;
  }
#endif

  LR::Element* el = lrspline->getElement(iel-1);
  double u[4] = { el->umin(), el->umax(), el->umin(), el->umax() };
  double v[4] = { el->vmin(), el->vmin(), el->vmax(), el->vmax() };

  XC.clear();
  XC.reserve(4);
  Go::Point point;

  for (int i = 0; i < 4; i++)
  {
    lrspline->point(point,u[i],v[i],iel-1);
    XC.push_back(SplineUtils::toVec3(point,nsd));
  }

  return getElementSize(XC);
}


bool ASMu2D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu2D::integrate(I)");

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

  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t t = 0; t < threadGroups[0].size() && ok; ++t)
  {
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < threadGroups[0][t].size(); ++e)
    {
      if (!ok)
        continue;

      int iel = threadGroups[0][t][e] + 1;

#ifdef SP_DEBUG
      if (dbgElm < 0 && iel != -dbgElm)
        continue; // Skipping all elements, except for -dbgElm
#endif
      FiniteElement fe(MNPC[iel-1].size());
      fe.iel = MLGE[iel-1];
      Matrix   dNdu, Xnod, Jac;
      Matrix3D d2Ndu2, Hess;
      Vec4     X;

      // Get element area in the parameter space
      double dA = this->getParametricArea(iel);
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

      // Compute parameter values of the Gauss points over this element
      std::array<RealArray,2> gpar, redpar;
      for (int d = 0; d < 2; d++)
      {
        this->getGaussPointParameters(gpar[d],d,nGauss,iel,xg);
        if (xr)
          this->getGaussPointParameters(redpar[d],d,nRed,iel,xr);
      }

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
        fe.h = this->getElementCorners(iel,fe.XC);

      if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
      {
        // Compute the element center
        Go::Point X0;
        double u0 = 0.5*(gpar[0].front() + gpar[0].back());
        double v0 = 0.5*(gpar[1].front() + gpar[1].back());
        lrspline->point(X0,u0,v0);
        for (unsigned char i = 0; i < nsd; i++)
          X[i] = X0[i];
      }

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
      if (!integrand.initElement(MNPC[iel-1],fe,X,nRed*nRed,*A))
      {
        ok = false;
        continue;
      }

      if (xr)
      {
        // --- Selective reduced integration loop ------------------------------

        for (int j = 0; j < nRed; j++)
          for (int i = 0; i < nRed; i++)
          {
            // Local element coordinates of current integration point
            fe.xi  = xr[i];
            fe.eta = xr[j];

            // Parameter values of current integration point
            fe.u = redpar[0][i];
            fe.v = redpar[1][j];

            // Compute basis function derivatives at current point
            Go::BasisDerivsSf spline;
            lrspline->computeBasis(fe.u,fe.v,spline,iel-1);
            SplineUtils::extractBasis(spline,fe.N,dNdu);

            // Compute Jacobian inverse and derivatives
            fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

            // Cartesian coordinates of current integration point
            X = Xnod * fe.N;
            X.t = time.t;

            // Compute the reduced integration terms of the integrand
            fe.detJxW *= 0.25*dA*wr[i]*wr[j];
            if (!integrand.reducedInt(*A,fe,X))
            {
              ok = false;
              continue;
            }
          }
      }

      // --- Integration loop over all Gauss points in each direction ----------

      int jp = (iel-1)*nGauss*nGauss;
      fe.iGP = firstIp + jp; // Global integration point counter

      for (int j = 0; j < nGauss; j++)
        for (int i = 0; i < nGauss; i++, fe.iGP++)
        {
          // Local element coordinates of current integration point
          fe.xi  = xg[i];
          fe.eta = xg[j];

          // Parameter values of current integration point
          fe.u = gpar[0][i];
          fe.v = gpar[1][j];

          // Compute basis function derivatives at current integration point
          if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES) {
            Go::BasisDerivsSf2 spline;
            lrspline->computeBasis(fe.u,fe.v,spline,iel-1);
            SplineUtils::extractBasis(spline,fe.N,dNdu,d2Ndu2);
          }
          else {
            Go::BasisDerivsSf spline;
            lrspline->computeBasis(fe.u,fe.v,spline, iel-1);
            SplineUtils::extractBasis(spline,fe.N,dNdu);
#if SP_DEBUG > 4
            if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
            {
              std::cout <<"\nBasis functions at a integration point "
                        <<" : (u,v) = "<< spline.param[0] <<" "<< spline.param[1]
                        <<"  left_idx = "<< spline.left_idx[0]
                        <<" "<< spline.left_idx[1];
              for (size_t ii = 0; ii < spline.basisValues.size(); ii++)
                std::cout <<'\n'<< 1+ii <<'\t' << spline.basisValues[ii] <<'\t'
                          << spline.basisDerivs_u[ii] <<'\t'
                          << spline.basisDerivs_v[ii];
              std::cout << std::endl;
            }
#endif
          }

          // Compute Jacobian inverse of coordinate mapping and derivatives
          fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
          if (fe.detJxW == 0.0) continue; // skip singular points

          // Compute Hessian of coordinate mapping and 2nd order derivatives
          if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
            if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
            {
              ok = false;
              continue;
            }

#if SP_DEBUG > 4
          if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
            std::cout <<"\nN ="<< fe.N <<"dNdX ="<< fe.dNdX;
#endif

          // Cartesian coordinates of current integration point
          X = Xnod * fe.N;
          X.t = time.t;

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= 0.25*dA*wg[i]*wg[j];
          PROFILE3("Integrand::evalInt");
          if (!integrand.evalInt(*A,fe,time,X))
          {
            ok = false;
            continue;
          }
        }

      // Finalize the element quantities
      if (!integrand.finalizeElement(*A,time,firstIp+jp))
      {
        ok = false;
        continue;
      }

      // Assembly of global system integral
      if (!glInt.assemble(A->ref(),fe.iel))
      {
        ok = false;
        continue;
      }

      A->destruct();

#ifdef SP_DEBUG
      if (iel == -dbgElm)
        continue; // Skipping all elements, except for -dbgElm
#endif
    }
  }

  return ok;
}


bool ASMu2D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time,
                        const Real3DMat& itgPts)
{
  if (!lrspline) return true; // silently ignore empty patches

  if (integrand.getReducedIntegration(nGauss) != 0)
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

  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t t = 0; t < threadGroups[0].size() && ok; ++t)
  {
#pragma omp parallel for schedule(static)
    for (size_t e = 0; e < threadGroups[0][t].size(); ++e)
    {
      if (!ok)
        continue;

      int iel = threadGroups[0][t][e] + 1;
#ifdef SP_DEBUG
      if (dbgElm < 0 && iel != -dbgElm)
        continue; // Skipping all elements, except for -dbgElm
#endif
      FiniteElement fe(MNPC[iel-1].size());
      fe.iel = MLGE[iel-1];
      Matrix   dNdu, Xnod, Jac;
      Matrix3D d2Ndu2, Hess;
      Vec4     X;

      // Get element area in the parameter space
      double dA = this->getParametricArea(iel);
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
        X = 0.25*(fe.XC[0]+fe.XC[1]+fe.XC[2]+fe.XC[3]);
      }

      // Initialize element quantities
      LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
      if (!integrand.initElement(MNPC[iel-1],fe,X,0,*A))
      {
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
        if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES) {
          Go::BasisDerivsSf2 spline;
          lrspline->computeBasis(fe.u,fe.v,spline,iel-1);
          SplineUtils::extractBasis(spline,fe.N,dNdu,d2Ndu2);
        }
        else {
          Go::BasisDerivsSf spline;
          lrspline->computeBasis(fe.u,fe.v,spline,iel-1);
          SplineUtils::extractBasis(spline,fe.N,dNdu);
        }

        // Compute Jacobian inverse of coordinate mapping and derivatives
        fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);
        if (fe.detJxW == 0.0) continue; // skip singular points

        // Compute Hessian of coordinate mapping and 2nd order derivatives
        if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
          if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
          {
            ok = false;
            continue;
          }

#if SP_DEBUG > 4
        if (iel == dbgElm || iel == -dbgElm || dbgElm == 0)
          std::cout <<"\niel, ip = "<< iel <<" "<< ip
                    <<"\nN ="<< fe.N <<"dNdX ="<< fe.dNdX;
#endif

        // Cartesian coordinates of current integration point
        X = Xnod * fe.N;
        X.t = time.t;

        // Evaluate the integrand and accumulate element contributions
        fe.detJxW *= 0.25*dA*elmPts[ip][2];
        PROFILE3("Integrand::evalInt");
        if (!integrand.evalInt(*A,fe,time,X))
        {
          ok = false;
          continue;
        }
      }

      // Finalize the element quantities
      if (!integrand.finalizeElement(*A,time,firstIp+MPitg[iel]))
      {
        ok = false;
        continue;
      }

      // Assembly of global system integral
      if (!glInt.assemble(A->ref(),fe.iel))
      {
        ok = false;
        continue;
      }

      A->destruct();
    }
  }

  return ok;
}


bool ASMu2D::integrate (Integrand& integrand, int lIndex,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu2D::integrate(B)");

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the edge normal {-2,-1, 1, 2}
  const int edgeDir = (lIndex%10+1)/(lIndex%2 ? -2 : 2);

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

  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   normal;


  // === Assembly loop over all elements on the patch edge =====================

  std::vector<LR::Element*>::iterator el = lrspline->elementBegin();
  for (int iel = 1; el != lrspline->elementEnd(); el++, iel++)
  {
#ifdef SP_DEBUG
    if (dbgElm < 0 && iel != -dbgElm)
      continue; // Skipping all elements, except for -dbgElm
#endif

    // Skip elements that are not on current boundary edge
    bool skipMe = false;
    switch (edgeDir)
    {
    case -1: if ((*el)->umin() != lrspline->startparam(0)) skipMe = true; break;
    case  1: if ((*el)->umax() != lrspline->endparam(0)  ) skipMe = true; break;
    case -2: if ((*el)->vmin() != lrspline->startparam(1)) skipMe = true; break;
    case  2: if ((*el)->vmax() != lrspline->endparam(1)  ) skipMe = true; break;
    }
    if (skipMe) continue;

    // Get element edge length in the parameter space
    double dS = this->getParametricLength(iel,t1);
    if (dS < 0.0) return false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel)) return false;

    // Initialize element quantities
    FiniteElement fe((**el).nBasisFunctions());
    fe.iel = MLGE[iel-1];
    fe.xi = fe.eta = edgeDir < 0 ? -1.0 : 1.0;
    LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
    if (!integrand.initElementBou(MNPC[iel-1],*A)) return false;

    // Get integration gauss points over this element
    this->getGaussPointParameters(gpar[t2-1],t2-1,nGP,iel,xg);

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      fe.h = this->getElementCorners(iel,fe.XC);

    // --- Integration loop over all Gauss points along the edge ---------------

    fe.iGP = firstp; // Global integration point counter
    firstp += nGP;

    for (int i = 0; i < nGP; i++, fe.iGP++)
    {
      // Local element coordinates and parameter values
      // of current integration point
      fe.xi = xg[i];
      fe.eta = xg[i];
      fe.u = gpar[0][i];
      fe.v = gpar[1][i];

      // Evaluate basis function derivatives at current integration points
      Go::BasisDerivsSf spline;
      lrspline->computeBasis(fe.u, fe.v, spline, iel-1);

      // Fetch basis function derivatives at current integration point
      SplineUtils::extractBasis(spline,fe.N,dNdu);

      // Compute basis function derivatives and the edge normal
      fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
      if (fe.detJxW == 0.0) continue; // skip singular points

      if (edgeDir < 0) normal *= -1.0;

      // Cartesian coordinates of current integration point
      X = Xnod * fe.N;
      X.t = time.t;

      // Evaluate the integrand and accumulate element contributions
      fe.detJxW *= 0.5*dS*wg[i];
      if (!integrand.evalBou(*A,fe,time,X,normal))
        return false;
    }

    // Finalize the element quantities
    if (!integrand.finalizeElementBou(*A,fe,time))
      return false;

    // Assembly of global system integral
    if (!glInt.assemble(A,fe.iel))
      return false;

    A->destruct();

#ifdef SP_DEBUG
    if (dbgElm < 0 && iel == -dbgElm)
      break; // Skipping all elements, except for -dbgElm
#endif
  }

  return true;
}


int ASMu2D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  if (!lrspline) return -2;

  param[0] = (1.0-xi[0])*lrspline->startparam(0) + xi[0]*lrspline->endparam(0);
  param[1] = (1.0-xi[1])*lrspline->startparam(1) + xi[1]*lrspline->endparam(1);

  FiniteElement fe;
  fe.iel = 1 + lrspline->getElementContaining(param[0],param[1]);
  fe.u   = param[0];
  fe.v   = param[1];

  Matrix Xnod;
  if (!this->getElementCoordinates(Xnod,fe.iel))
    return -1;

  if (!this->evaluateBasis(fe))
    return -1;

  X = Xnod * fe.N;

  return 0;
}


bool ASMu2D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
#ifdef SP_DEBUG
  std::cout << "ASMu2D::getGridParameters(  )\n";
#endif

  // output is written once for each element resulting in a lot of unnecessary storage
  // this is preferable to figuring out all element topology information

  std::vector<LR::Element*>::iterator el;
  for(el=lrspline->elementBegin(); el<lrspline->elementEnd(); el++) {
    // evaluate element at element corner points
    double umin = (**el).umin();
    double umax = (**el).umax();
    double vmin = (**el).vmin();
    double vmax = (**el).vmax();
    for(int iv=0; iv<=nSegPerSpan; iv++) {
      for(int iu=0; iu<=nSegPerSpan; iu++) {
        double u = umin + (umax-umin)/nSegPerSpan*iu;
        double v = vmin + (vmax-vmin)/nSegPerSpan*iv;
        if (dir==0)
          prm.push_back(u);
        else
          prm.push_back(v);
      }
    }
  }
  return true;
}

#if 0
  if (!lrspline) return false;

  if (nSegPerSpan < 1)
  {
    std::cerr <<" *** ASMu2D::getGridParameters: Too few knot-span points "
              << nSegPerSpan+1 <<" in direction "<< dir << std::endl;
    return false;
  }

  RealArray::const_iterator uit = lrspline->basis(dir).begin();
  double ucurr = 0.0, uprev = *(uit++);
  while (uit != lrspline->basis(dir).end())
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


bool ASMu2D::getGrevilleParameters (RealArray& prm, int dir) const
{
  if (!lrspline) return false;

  const Go::BsplineBasis& basis = lrspline->basis(dir);

  prm.resize(basis.numCoefs());
  for (size_t i = 0; i < prm.size(); i++)
    prm[i] = basis.grevilleParameter(i);

  return true;
}
#endif


bool ASMu2D::tesselate (ElementBlock& grid, const int* npe) const
{
#ifdef SP_DEBUG
  std::cout << "ASMu2D::tesselate(  )\n";
#endif
  if (!lrspline) return false;

  if (npe[0] != npe[1]) {
    std::cerr << "ASMu2D::tesselate does not support different tesselation resolution in "
              << "u- and v-direction. nviz u = " << npe[0] << ", nviz v = " << npe[1] << "\n";
    return false;
  }

  int nNodesPerElement =  npe[0]   * npe[1];
  int nSubElPerElement = (npe[0]-1)*(npe[1]-1);
  int nElements        = lrspline->nElements();

  // output is written once for each element resulting in a lot of unnecessary storage
  // this is preferable to figuring out all element topology information
  grid.unStructResize(nElements * nSubElPerElement,
                      nElements * nNodesPerElement);

  std::vector<LR::Element*>::iterator el;
  int inod = 0;
  int iel = 0;
  for(el=lrspline->elementBegin(); el<lrspline->elementEnd(); el++, iel++) {
    // evaluate element at element corner points
    double umin = (**el).umin();
    double umax = (**el).umax();
    double vmin = (**el).vmin();
    double vmax = (**el).vmax();
    for(int iv=0; iv<npe[1]; iv++) {
      for(int iu=0; iu<npe[0]; iu++) {
        double u = umin + (umax-umin)/(npe[0]-1)*iu;
        double v = vmin + (vmax-vmin)/(npe[1]-1)*iv;
        Go::Point pt;
        lrspline->point(pt, u,v, iel, iu!=npe[0]-1, iv!=npe[1]-1);
        for(int dim=0; dim<nsd; dim++)
          grid.setCoor(inod, dim, pt[dim]);
        inod++;
      }
    }
  }

  int ip = 0;
  iel = 0;
  for(int i=0; i<lrspline->nElements(); i++) {
    int iStart = i*nNodesPerElement;
    for(int iv=0; iv<npe[1]-1; iv++) {
      for(int iu=0; iu<npe[0]-1; iu++, iel++) {
        // enumerate nodes counterclockwise around the quad
        grid.setNode(ip++, iStart + (iv  )*npe[0] + (iu  ) );
        grid.setNode(ip++, iStart + (iv  )*npe[0] + (iu+1) );
        grid.setNode(ip++, iStart + (iv+1)*npe[0] + (iu+1) );
        grid.setNode(ip++, iStart + (iv+1)*npe[0] + (iu  ) );
        grid.setElmId(iel+1, i+1);
      }
    }
  }

  return true;
}


bool ASMu2D::evalSolution (Matrix& sField, const Vector& locSol,
                           const int* npe) const
{
#ifdef SP_DEBUG
  std::cout <<"ASMu2D::evalSolution(Matrix&,const Vector&,const int*)\n";
#endif
  // Compute parameter values of the result sampling points
  std::array<RealArray,2> gpar;
  for (int dir = 0; dir < 2; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  return this->evalSolution(sField,locSol,gpar.data());
}


bool ASMu2D::evalSolution (Matrix& sField, const Vector& locSol,
                           const RealArray* gpar, bool, int deriv) const
{
#ifdef SP_DEBUG
  std::cout <<"ASMu2D::evalSolution(Matrix&,const Vector&,const RealArray*,"
            <<"bool,int)"<< std::endl;
#endif
  size_t nComp = locSol.size() / this->getNoNodes();
  if (nComp*this->getNoNodes() != locSol.size())
    return false;

  size_t nPoints = gpar[0].size();
  if (nPoints != gpar[1].size())
    return false;

  Vector   ptSol;
  Matrix   dNdu, dNdX, Jac, Xnod, eSol, ptDer;
  Matrix3D d2Ndu2, d2NdX2, Hess, ptDer2;

  Go::BasisPtsSf     spline0;
  Go::BasisDerivsSf  spline1;
  Go::BasisDerivsSf2 spline2;

  // Evaluate the primary solution field at each point
  sField.resize(nComp,nPoints);
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch element containing evaluation point.
    // Sadly, points are not always ordered in the same way as the elements.
    int iel = lrspline->getElementContaining(gpar[0][i],gpar[1][i]);
    FiniteElement fe(lrspline->getElement(iel)->nBasisFunctions());
    fe.iel = iel + 1;
    fe.u   = gpar[0][i];
    fe.v   = gpar[1][i];
    utl::gather(MNPC[iel],nComp,locSol,eSol);

    // Set up control point (nodal) coordinates for current element
    if (deriv > 0 && !this->getElementCoordinates(Xnod,iel+1))
      return false;

    // Evaluate basis function values/derivatives at current parametric point
    // and multiply with control point values to get the point-wise solution
    switch (deriv)
    {
    case 0: // Evaluate the solution
      if (!this->evaluateBasis(fe,deriv))
        return false;
      sField.fillColumn(1+i, eSol * fe.N);
      break;

    case 1: // Evaluate first derivatives of the solution
      if (!this->evaluateBasis(fe,deriv))
        return false;
      ptDer.multiply(eSol,fe.dNdX);
      sField.fillColumn(1+i,ptDer);
      break;

    case 2: // Evaluate second derivatives of the solution
      lrspline->computeBasis(gpar[0][i],gpar[1][i],spline2,iel);
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
    for (int i = 0; i < s->nBasisFunctions(); i++)
      sField.fillColumn(i+1,&(*s->getBasisfunction(i)->cp()));
    delete s;
    return true;
  }

  std::cerr <<" *** ASMu2D::evalSolution: Failure!"<< std::endl;
  return false;
}


bool ASMu2D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                           const RealArray* gpar, bool) const
{
#ifdef SP_DEBUG
  std::cout <<"ASMu2D::evalSolution(Matrix&,const IntegrandBase&,const RealArray*,bool)\n";
#endif

  sField.resize(0,0);

  // TODO: investigate the possibility of doing "regular" refinement by
  //       uniform tesselation grid and ignoring LR mesh lines

  size_t nPoints = gpar[0].size();
  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  if (nPoints != gpar[1].size())
    return false;

  Vector   solPt;
  Matrix   dNdu, Jac, Xnod;
  Matrix3D d2Ndu2, Hess;

  // Evaluate the secondary solution field at each point
  for (size_t i = 0; i < nPoints; i++)
  {
    // Fetch element containing evaluation point
    // sadly, points are not always ordered in the same way as the elements
    int iel = lrspline->getElementContaining(gpar[0][i],gpar[1][i]);

    // Evaluate the basis functions at current parametric point
    FiniteElement fe(lrspline->getElement(iel)->nBasisFunctions(),firstIp+i);
    if (use2ndDer)
    {
      Go::BasisDerivsSf2 spline;
      lrspline->computeBasis(gpar[0][i],gpar[1][i],spline,iel);
      SplineUtils::extractBasis(spline,fe.N,dNdu,d2Ndu2);
    }
    else
    {
      Go::BasisDerivsSf spline;
      lrspline->computeBasis(gpar[0][i],gpar[1][i],spline,iel);
      SplineUtils::extractBasis(spline,fe.N,dNdu);
    }

    // Set up control point (nodal) coordinates for current element
    if (!this->getElementCoordinates(Xnod,iel+1)) return false;

    // Compute the Jacobian inverse
    fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

    // Compute Hessian of coordinate mapping and 2nd order derivatives
    if (use2ndDer)
      if (!utl::Hessian(Hess,fe.d2NdX2,Jac,Xnod,d2Ndu2,dNdu))
        continue;

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,Xnod*fe.N,MNPC[iel]))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


void ASMu2D::getBoundaryNodes (int lIndex, IntVec& nodes, int basis,
                               int, bool local) const
{
  if (basis == 0)
    basis = 1;

  if (!this->getBasis(basis)) return; // silently ignore empty patches

  LR::parameterEdge edge;
  switch (lIndex) {
  case 1: edge = LR::WEST; break;
  case 2: edge = LR::EAST; break;
  case 3: edge = LR::SOUTH; break;
  case 4: edge = LR::NORTH; break;
  default: return;
  }

  nodes = this->getEdgeNodes(edge, basis);

#if SP_DEBUG > 1
  std::cout <<"Boundary nodes in patch "<< idx+1 <<" edge "<< lIndex <<":";
  for (size_t i = 0; i < nodes.size(); i++)
    std::cout <<" "<< nodes[i];
  std::cout << std::endl;
#endif
}


bool ASMu2D::getOrder (int& p1, int& p2, int& p3) const
{
  p1 = geo->order(0);
  p2 = geo->order(1);
  p3 = 0;

  return true;
}


bool ASMu2D::updateDirichlet (const std::map<int,RealFunc*>& func,
                              const std::map<int,VecFunc*>& vfunc, double time,
                              const std::map<int,int>* g2l)
{

  std::map<int,RealFunc*>::const_iterator    fit;
  std::map<int,VecFunc*>::const_iterator     vfit;
  std::map<int,int>::const_iterator          nit;

  for (size_t i = 0; i < dirich.size(); i++)
  {
    // figure out function index offset (when using multiple basis)
    size_t ofs = 0;
    for (int j = 1; j < dirich[i].basis; j++)
      ofs += this->getNoNodes(j);

    RealArray edgeControlpoints;
    Real2DMat edgeControlmatrix;
    if ((fit = func.find(dirich[i].code)) != func.end())
      edgeL2projection(dirich[i], *fit->second, edgeControlpoints, time);
    else if ((vfit = vfunc.find(dirich[i].code)) != vfunc.end())
      edgeL2projection(dirich[i], *vfit->second, edgeControlmatrix, time);
    else
    {
      std::cerr <<" *** ASMu2D::updateDirichlet: Code "<< dirich[i].code
                <<" is not associated with any function."<< std::endl;
      return false;
    }
    if (edgeControlpoints.empty() && edgeControlmatrix.empty())
    {
      std::cerr <<" *** ASMu2D::updateDirichlet: Projection failure."
                << std::endl;
      return false;
    }

    // Loop over the nodes of this boundary curve
    for (nit = dirich[i].MLGN.begin(); nit != dirich[i].MLGN.end(); nit++)
    {
      // skip corner nodes, since these are special cased (interpolatory)
      if (nit->second == dirich[i].corners[0] ||
          nit->second == dirich[i].corners[1])
        continue;
      for (int dofs = dirich[i].dof; dofs > 0; dofs /= 10)
      {
        int dof = dofs%10;
        // Find the constraint equation for current (node,dof)
        MPC pDOF(MLGN[nit->second]+ofs,dof);
        MPCIter mit = mpcs.find(&pDOF);
        if (mit == mpcs.end()) continue; // probably a deleted constraint

        // Now update the prescribed value in the constraint equation
        if (edgeControlpoints.empty()) // vector conditions
          (*mit)->setSlaveCoeff(edgeControlmatrix[dof-1][nit->first]);
        else                           //scalar condition
          (*mit)->setSlaveCoeff(edgeControlpoints[nit->first]);
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


size_t ASMu2D::getNoNodes (int) const
{
  return lrspline->nBasisFunctions();
}


bool ASMu2D::transferGaussPtVars (const LR::LRSplineSurface* oldBasis,
                                  const RealArray& oldVars, RealArray& newVars,
                                  int nGauss) const
{
  const LR::LRSplineSurface* newBasis = this->getBasis();

  size_t nGp = nGauss*nGauss;
  newVars.resize(newBasis->nElements()*nGp);

  const double* xi = GaussQuadrature::getCoord(nGauss);
  LagrangeInterpolator interp(Vector(xi,nGauss));

  for (int iEl = 0; iEl < newBasis->nElements(); iEl++)
  {
    const LR::Element* newEl = newBasis->getElement(iEl);
    double u_center = 0.5*(newEl->umin() + newEl->umax());
    double v_center = 0.5*(newEl->vmin() + newEl->vmax());
    int iOld = oldBasis->getElementContaining(u_center,v_center);
    if (iOld < 0)
    {
      std::cerr <<" *** ASMu2D: Failed to locate element "<< iEl
                <<" of the old mesh in the new mesh."<< std::endl;
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
      for (size_t j = 0; j < UGP.size(); j++)
        UGP[j] = -1.0 + 2.0*(UGP[j]-pmin)/(pmax-pmin);

      // lagrangian interpolation
      I[i] = interp.get(UGP);
    }

    Matrix data(nGauss,nGauss);
    data.fill(&oldVars[nGp*iOld],nGp);

    Matrix newdata;
    newdata.multiply(I[0]*data,I[1],false,true);
    std::copy(newdata.ptr(), newdata.ptr()+nGp, newVars.begin()+iEl*nGp);
  }

  return true;
}


bool ASMu2D::transferCntrlPtVars (LR::LRSplineSurface* oldBasis,
                                  const RealArray& oldVars, RealArray& newVars,
                                  int nGauss) const
{
  oldBasis->rebuildDimension(1);
  oldBasis->setControlPoints(const_cast<RealArray&>(oldVars));
  return this->transferCntrlPtVars(oldBasis,newVars,nGauss);
}


bool ASMu2D::transferCntrlPtVars (const LR::LRSplineSurface* oldBasis,
                                  RealArray& newVars, int nGauss) const
{
  const LR::LRSplineSurface* newBasis = this->getBasis();

  newVars.clear();
  newVars.reserve(newBasis->nElements()*nGauss*nGauss*oldBasis->dimension());
  const double* xi = GaussQuadrature::getCoord(nGauss);

  for (int iEl = 0; iEl < newBasis->nElements(); iEl++)
  {
    RealArray U, V, ptVar;
    LR::getGaussPointParameters(newBasis, U, 0, nGauss, iEl+1, xi);
    LR::getGaussPointParameters(newBasis, V, 1, nGauss, iEl+1, xi);
    for (int j = 0; j < nGauss; j++)
      for (int i = 0; i < nGauss; i++)
      {
        oldBasis->point(ptVar,U[i],V[j]);
        for (size_t k = 0; k < ptVar.size(); k++)
          newVars.push_back(ptVar[k]);
      }
  }

  return true;
}


void ASMu2D::generateThreadGroups (const Integrand& integrand, bool silence,
                                   bool ignoreGlobalLM)
{
  LR::generateThreadGroups(threadGroups, this->getBasis(1));
  if (silence || threadGroups[0].size() < 2) return;

  std::cout <<"\nMultiple threads are utilized during element assembly.";
  for (size_t i = 0; i < threadGroups[0].size(); i++)
  {
    std::cout <<"\n Color "<< i+1;
    std::cout << ": "<< threadGroups[0][i].size() <<" elements";
  }
}
