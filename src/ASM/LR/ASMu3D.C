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
#include "ElementBlock.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Profiler.h"
#include "Vec3Oper.h"
#include <array>


ASMu3D::ASMu3D (unsigned char n_f)
  : ASMunstruct(3,3,n_f), lrspline(nullptr), tensorspline(nullptr),
    bezierExtract(myBezierExtract)
{
  ASMunstruct::resetNumbering(); // Replace this when going multi-patch...
}


ASMu3D::ASMu3D (const ASMu3D& patch, unsigned char n_f)
  : ASMunstruct(patch,n_f), lrspline(patch.lrspline), tensorspline(nullptr),
    bezierExtract(patch.myBezierExtract)
{
  // Need to set nnod here,
  // as hasXNodes might be invoked before the FE data is generated
  if (nnod == 0 && lrspline)
    nnod = lrspline->nBasisFunctions();
}


bool ASMu3D::read (std::istream& is)
{
  if (shareFE) return true;
        lrspline.reset();

  // read inputfile as either an LRSpline file directly or a tensor product B-spline and convert
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
}


bool ASMu3D::checkRightHandSystem ()
{
  std::cout <<"ASMu3D::checkRightHandSystem(): Not available for LR-splines (ignored)."<< std::endl;
  return false;
}


bool ASMu3D::refine (int dir, const RealArray& xi)
{
  if (!tensorspline || dir < 0 || dir > 2 || xi.empty()) return false;
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

  tensorspline->insertKnot(dir,extraKnots);
        lrspline.reset(new LR::LRSplineVolume(tensorspline));
        geo = lrspline.get();
  return true;
}

bool ASMu3D::uniformRefine (int dir, int nInsert)
{
  if (!tensorspline || dir < 0 || dir > 2 || nInsert < 1) return false;
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

  if (!MLGN.empty() && MLGN.size() != nnod)
  {
    std::cerr <<" *** ASMu3D::generateFEMTopology: Inconsistency"
              <<" between the number of FE nodes "<< MLGN.size()
              <<"\n     and the number of basis functions "<< nnod
              <<" in the patch."<< std::endl;
    return false;
  }
  else if (shareFE)
    return true;

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);
  const int p3 = lrspline->order(2);

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
      myBezierExtract[iel].resize((*eit)->nBasisFunctions(),p1*p2*p3);
      myBezierExtract[iel].fill(extrMat.data(),extrMat.size());
    }
  }

  for (size_t inod = 0; inod < nnod; inod++)
    myMLGN[inod] = ++gNod;

  return true;
}


/*
// this is connecting multiple patches and handling deformed geometries.
// We'll deal with at a later time, for now we only allow single patch models

bool ASMu3D::connectPatch (int face, ASMu3D& neighbor, int nface, int norient)
{
  return this->connectBasis(face,neighbor,nface,norient);
}


bool ASMu3D::connectBasis (int face, ASMu3D& neighbor, int nface, int norient,
                           int basis, int slave, int master)
{
  if (shareFE && neighbor.shareFE)
    return true;
  else if (shareFE || neighbor.shareFE)
  {
    std::cerr <<" *** ASMu3D::connectPatch: Logic error, cannot"
        <<" connect a sharedFE patch with an unshared one"<< std::endl;
    return false;
  }

  // Set up the slave node numbers for this volume patch

  int n1, n2, n3;
  if (!this->getSize(n1,n2,n3,basis)) return false;
  int node = slave+1, i1 = 0, i2 = 0;

  switch (face)
    {
    case 2: // Positive I-direction
      node += n1-1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      n2 = n3;
      break;

    case 4: // Positive J-direction
      node += n1*(n2-1);
    case 3: // Negative J-direction
      i2 = n1*(n2-1);
      i1 = 1;
      n2 = n3;
      break;

    case 6: // Positive K-direction
      node += n1*n2*(n3-1);
    case 5: // Negative K-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** ASMu3D::connectPatch: Invalid slave face "
    << face << std::endl;
      return false;
    }

  int i, j;
  IntMat slaveNodes(n1,IntVec(n2,0));
  for (j = 0; j < n2; j++, node += i2)
    for (i = 0; i < n1; i++, node += i1)
      slaveNodes[i][j] = node;

  // Set up the master node numbers for the neighboring volume patch

  if (!neighbor.getSize(n1,n2,n3,basis)) return false;
  node = master+1; i1 = i2 = 0;

  switch (nface)
    {
    case 2: // Positive I-direction
      node += n1-1;
    case 1: // Negative I-direction
      i1 = n1;
      n1 = n2;
      n2 = n3;
      break;

    case 4: // Positive J-direction
      node += n1*(n2-1);
    case 3: // Negative J-direction
      i2 = n1*(n2-1);
      i1 = 1;
      n2 = n3;
      break;

    case 6: // Positive K-direction
      node += n1*n2*(n3-1);
    case 5: // Negative K-direction
      i1 = 1;
      break;

    default:
      std::cerr <<" *** ASMu3D::connectPatch: Invalid master face "
    << nface << std::endl;
      return false;
    }

  if (norient < 0 || norient > 7)
  {
    std::cerr <<" *** ASMu3D::connectPatch: Orientation flag "
        << norient <<" is out of range [0,7]"<< std::endl;
    return false;
  }

  int m1 = slaveNodes.size();
  int m2 = slaveNodes.front().size();
  if (norient < 4 ? (n1 != m1 || n2 != m2) : (n2 != m1 || n1 != m2))
  {
    std::cerr <<" *** ASMu3D::connectPatch: Non-matching faces, sizes "
        << n1 <<","<< n2 <<" and "<< m1 <<","<< m2 << std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  for (j = 0; j < n2; j++, node += i2)
    for (i = 0; i < n1; i++, node += i1)
    {
      int k = i, l = j;
      switch (norient)
  {
  case  1: k =    i  ; l = n2-j-1; break;
  case  2: k = n1-i-1; l =    j  ; break;
  case  3: k = n1-i-1; l = n2-j-1; break;
  case  4: k =    j  ; l =    i  ; break;
  case  5: k =    j  ; l = n1-i-1; break;
  case  6: k = n2-j-1; l =    i  ; break;
  case  7: k = n2-j-1; l = n1-i-1; break;
  default: k =    i  ; l = j     ;
  }

      int slave = slaveNodes[k][l];
      if (!neighbor.getCoord(node).equal(this->getCoord(slave),xtol))
      {
  std::cerr <<" *** ASMu3D::connectPatch: Non-matching nodes "
      << node <<": "<< neighbor.getCoord(node)
      <<"\n                                          and "
      << slave <<": "<< this->getCoord(slave) << std::endl;
  return false;
      }
      else
  ASMbase::collapseNodes(neighbor,node,*this,slave);
    }

  return true;
}


void ASMu3D::closeFaces (int dir, int basis, int master)
{
  int n1, n2, n3;
  if (basis < 1) basis = 1;
  if (!this->getSize(n1,n2,n3,basis)) return;

  switch (dir)
    {
    case 1: // Faces are closed in I-direction
      for (int i3 = 1; i3 <= n3; i3++)
  for (int i2 = 1; i2 <= n2; i2++, master += n1)
    this->makePeriodic(master,master+n1-1);
      break;

    case 2: // Faces are closed in J-direction
      for (int i3 = 1; i3 <= n3; i3++, master += n1*(n2-1))
  for (int i1 = 1; i1 <= n1; i1++, master++)
    this->makePeriodic(master,master+n1*(n2-1));
      break;

    case 3: // Faces are closed in K-direction
      for (int i2 = 1; i2 <= n2; i2++)
  for (int i1 = 1; i1 <= n1; i1++, master++)
    this->makePeriodic(master,master+n1*n2*(n3-1));
      break;
    }
}
*/


/*!
  A negative \a code value implies direct evaluation of the Dirichlet condition
  function at the control point. Positive \a code implies projection onto the
  spline basis representing the boundary surface (needed for curved faces and/or
  non-constant functions).
*/

void ASMu3D::constrainFace (int dir, bool open, int dof, int code, char)
{
  if (open)
    std::cerr << "\nWARNING: ASMu3D::constrainFace, open boundary conditions not supported yet. Treating it as closed" << std::endl;

  int bcode = code;
  if (code > 0) {// Dirichlet projection will be performed
    std::cerr << "\nWARNING: Projective (nonhomogenuous) dirichlet boundary conditions not implemented. ";
    std::cerr << "\n         Performing variational diminishing approximation instead" << std::endl;
    // dirich.push_back(DirichletFace(this->getBoundary(dir),dof,code));
  }
  else if (code < 0)
    bcode = -code;

  // get all the boundary functions from the LRspline object
  std::vector<LR::Basisfunction*> thisEdge;
  if (dir == -1)
    lrspline->getEdgeFunctions(thisEdge, LR::WEST, 1);
  else if (dir == 1)
    lrspline->getEdgeFunctions(thisEdge, LR::EAST, 1);
  else if (dir == -2)
    lrspline->getEdgeFunctions(thisEdge, LR::SOUTH, 1);
  else if (dir == 2)
    lrspline->getEdgeFunctions(thisEdge, LR::NORTH, 1);
  else if (dir == -3)
    lrspline->getEdgeFunctions(thisEdge, LR::BOTTOM, 1);
  else if (dir == 3)
    lrspline->getEdgeFunctions(thisEdge, LR::TOP, 1);

  std::cout << "\nNumber of constraints: " << thisEdge.size() << std::endl;

  // enforce the boundary conditions
  for(LR::Basisfunction* b  : thisEdge)
    this->prescribe(b->getId()+1,dof,bcode);
}

size_t ASMu3D::constrainFaceLocal(int dir, bool open, int dof, int code, bool project, char T1)
{
  std::cerr << "ASMu3D::constrainFaceLocal not implemented properly yet" << std::endl;
  exit(776654);
  return 0;
}


std::vector<int> ASMu3D::getEdge(int lEdge, bool open, int) const
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

  // get all the boundary functions from the LRspline object
  std::vector<LR::Basisfunction*> thisEdge;
  lrspline->getEdgeFunctions(thisEdge, (LR::parameterEdge) edge, 1);
  std::vector<int> result;
  for (LR::Basisfunction* b  : thisEdge)
    result.push_back(b->getId()+1);

  return result;
}

void ASMu3D::constrainEdge (int lEdge, bool open, int dof, int code, char)
{
  if (open)
    std::cerr << "\nWARNING: ASMu3D::constrainEdge, open boundary conditions not supported yet. Treating it as closed" << std::endl;

  // enforce the boundary conditions
  std::vector<int> nodes = this->getEdge(lEdge, open, 1);
  for (const int& node : nodes)
    this->prescribe(node,dof,code);
}


void ASMu3D::constrainLine (int fdir, int ldir, double xi, int dof,
                            int code, char)
{
  std::cerr << "ASMu3D::constrainLine not implemented properly yet" << std::endl;
  exit(776654);
#if 0
  if (xi < 0.0 || xi > 1.0) return;

  int n1, n2, n3, node = 1;
  if (!this->getSize(n1,n2,n3,1)) return;

  switch (fdir)
    {
    case  1: // Right face (positive I-direction)
      node += n1-1;
    case -1: // Left face (negative I-direction)
      if (ldir == 2)
      {
  // Line goes in J-direction
  node += n1*n2*int(0.5+(n3-1)*xi);
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
  node += n1*n2*int(0.5+(n3-1)*xi);
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
#endif
}


void ASMu3D::constrainCorner (int I, int J, int K, int dof, int code, char)
{
  int corner = LR::NONE;
  if (  I < 0) corner |= LR::WEST;
  else        corner |= LR::EAST;
  if (  J < 0) corner |= LR::SOUTH;
  else        corner |= LR::NORTH;
  if (  K < 0) corner |= LR::BOTTOM;
  else        corner |= LR::TOP;

  std::vector<LR::Basisfunction*> one_function;
  lrspline->getEdgeFunctions(one_function, LR::parameterEdge(corner));
  int node = one_function[0]->getId() + 1;

  this->prescribe(node, dof, code);
}


void ASMu3D::constrainNode (double xi, double eta, double zeta,
                            int dof, int code, char)
{
  std::cerr << "ASMu3D::constrainNode not implemented properly yet" << std::endl;
  exit(776654);
  if (xi   < 0.0 || xi   > 1.0) return;
  if (eta  < 0.0 || eta  > 1.0) return;
  if (zeta < 0.0 || zeta > 1.0) return;

#if 0
  int n1, n2, n3;
  if (!this->getSize(n1,n2,n3,1)) return;

  int node = 1;
  if (xi   > 0.0) node += int(0.5+(n1-1)*xi);
  if (eta  > 0.0) node += n1*int(0.5+(n2-1)*eta);
  if (zeta > 0.0) node += n1*n2*int(0.5+(n3-1)*zeta);

  this->prescribe(node,dof,code);
#endif
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

  LR::Element *el = lrspline->getElement(iel-1);
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

  LR::Element *el = lrspline->getElement(iel-1);
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

  LR::Element* el = lrspline->getElement(iel-1);
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
  std::cerr <<" *** ASMu3D::updateCoords not implemented!"<< std::endl;
  return false;
}


size_t ASMu3D::getNoBoundaryElms (char lIndex, char ldim) const
{
  if (!lrspline) return 0;

  LR::parameterEdge edge;
  switch(lIndex)
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
  LR::Element* el = lrspline->getElement(iEl);
  double start = el->getParmin(dir);
  double stop  = el->getParmax(dir);

  uGP.resize(nGauss);
  for (int i = 0; i < nGauss; i++)
    uGP[i] = 0.5*((stop-start)*xi[i] + stop+start);
}


void ASMu3D::getElementCorners (int iEl, Vec3Vec& XC) const
{
  LR::Element* el = lrspline->getElement(iEl);
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
        lrspline->point(pt,u[i],v[j],w[k],iEl);
        XC.push_back(SplineUtils::toVec3(pt));
      }
}


void ASMu3D::evaluateBasis (FiniteElement &el, Matrix &dNdu) const
{
  PROFILE2("Spline evaluation");
  size_t nBasis = lrspline->getElement(el.iel-1)->nBasisFunctions();

  std::vector<RealArray> result;
  lrspline->computeBasis(el.u, el.v, el.w, result, 1, el.iel-1);

  el.N.resize(nBasis);
  dNdu.resize(nBasis,3);
  size_t jp, n = 1;
  for (jp = 0; jp < nBasis; jp++, n++) {
    el.N(n)   = result[jp][0];
    dNdu(n,1) = result[jp][1];
    dNdu(n,2) = result[jp][2];
    dNdu(n,3) = result[jp][3];
  }
}

void ASMu3D::evaluateBasis (FiniteElement &el, Matrix &dNdu, Matrix &C, Matrix &B) const
{
  PROFILE2("BeSpline evaluation");
  Matrix N = C*B;
  size_t nBasis = lrspline->getElement(el.iel-1)->nBasisFunctions();
  el.N = N.getColumn(1);
  dNdu.resize(nBasis,3);
  dNdu.fillColumn(1,N.getColumn(2));
  dNdu.fillColumn(2,N.getColumn(3));
  dNdu.fillColumn(3,N.getColumn(4));
}

void ASMu3D::evaluateBasis (FiniteElement &el, Matrix &dNdu, Matrix3D &d2Ndu2) const
{
  PROFILE2("Spline evaluation");
  size_t nBasis = lrspline->getElement(el.iel-1)->nBasisFunctions();

  std::vector<RealArray> result;
  lrspline->computeBasis(el.u, el.v, el.w, result, 2, el.iel-1 );

  el.N.resize(nBasis);
  dNdu.resize(nBasis,3);
  d2Ndu2.resize(nBasis,3,3);
  size_t jp, n = 1;
  for (jp = 0; jp < nBasis; jp++, n++) {
    el.N   (n)    = result[jp][0];
    dNdu (n,1)    = result[jp][1];
    dNdu (n,2)    = result[jp][2];
    dNdu (n,3)    = result[jp][3];
    d2Ndu2(n,1,1) = result[jp][4];
    d2Ndu2(n,1,2) = d2Ndu2(n,2,1) = result[jp][5];
    d2Ndu2(n,1,3) = d2Ndu2(n,3,1) = result[jp][6];
    d2Ndu2(n,2,2) = result[jp][7];
    d2Ndu2(n,2,3) = d2Ndu2(n,3,2) = result[jp][8];
    d2Ndu2(n,3,3) = result[jp][9];
  }
}

void ASMu3D::evaluateBasis (FiniteElement &el, int derivs) const
{
  PROFILE2("Spline evaluation");
  size_t nBasis = lrspline->getElement(el.iel-1)->nBasisFunctions();

  std::vector<RealArray> result;
  lrspline->computeBasis(el.u, el.v, el.w, result, derivs, el.iel-1 );

  el.N.resize(nBasis);
  el.dNdX.resize(nBasis,3);
  el.d2NdX2.resize(nBasis,3,3);
  size_t jp, n = 1;
  for (jp = 0; jp < nBasis; jp++, n++) {
    el.N   (n)      = result[jp][0];
    if (derivs > 0) {
      el.dNdX (n,1)    = result[jp][1];
      el.dNdX (n,2)    = result[jp][2];
      el.dNdX (n,3)    = result[jp][3];
    }
    if (derivs > 1) {
      el.d2NdX2(n,1,1) = result[jp][4];
      el.d2NdX2(n,1,2) = el.d2NdX2(n,2,1) = result[jp][5];
      el.d2NdX2(n,1,3) = el.d2NdX2(n,3,1) = result[jp][6];
      el.d2NdX2(n,2,2) = result[jp][7];
      el.d2NdX2(n,2,3) = el.d2NdX2(n,3,2) = result[jp][8];
      el.d2NdX2(n,3,3) = result[jp][9];
    }
  }
}

bool ASMu3D::integrate (Integrand& integrand,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu3D::integrate(I)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // evaluate all gauss points on the bezier patch (-1, 1)
  int p1 = lrspline->order(0);
  int p2 = lrspline->order(1);
  int p3 = lrspline->order(2);
  Go::BsplineBasis basis1 = getBezierBasis(p1);
  Go::BsplineBasis basis2 = getBezierBasis(p2);
  Go::BsplineBasis basis3 = getBezierBasis(p3);

  Matrix BN(   p1*p2*p3, nGauss*nGauss*nGauss);
  Matrix BdNdu(p1*p2*p3, nGauss*nGauss*nGauss);
  Matrix BdNdv(p1*p2*p3, nGauss*nGauss*nGauss);
  Matrix BdNdw(p1*p2*p3, nGauss*nGauss*nGauss);
  int ig=1; // gauss point iterator
  for(int zeta=0; zeta<nGauss; zeta++) {
    for(int eta=0; eta<nGauss; eta++) {
      for(int xi=0; xi<nGauss; xi++, ig++) {
        double u[2*p1];
        double v[2*p2];
        double w[2*p3];
        basis1.computeBasisValues(xg[xi],   u, 1);
        basis2.computeBasisValues(xg[eta],  v, 1);
        basis3.computeBasisValues(xg[zeta], w, 1);
        int ib=1; // basis function iterator
        double sum = 0;
        for(int k=0; k<p3; k++) {
          for(int j=0; j<p2; j++) {
            for(int i=0; i<p1; i++, ib++) {
              BN(ib,ig)    = u[2*i  ]*v[2*j  ]*w[2*k  ];
              BdNdu(ib,ig) = u[2*i+1]*v[2*j  ]*w[2*k  ];
              BdNdv(ib,ig) = u[2*i  ]*v[2*j+1]*w[2*k  ];
              BdNdw(ib,ig) = u[2*i  ]*v[2*j  ]*w[2*k+1];
              sum += BN(ib,ig);
            }
          }
        }
      }
    }
  }

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

  bool ok=true;
  for(LR::Element *el : lrspline->getAllElements()) {
    int iEl = el->getId();
    int nBasis = el->nBasisFunctions();
    FiniteElement fe(nBasis);
    fe.iel = iEl+1;

    Matrix   C = bezierExtract[iEl];
    Matrix   dNdu, Xnod, Jac;
    Matrix3D d2Ndu2, Hess;
    double   dXidu[3];
    Vec4     X;
    // Get element volume in the parameter space
    double du = el->umax() - el->umin();
    double dv = el->vmax() - el->vmin();
    double dw = el->wmax() - el->wmin();
    double vol = el->volume();
    if (vol < 0.0)
    {
      ok = false; // topology error (probably logic error)
      break;
    }

    // Set up control point (nodal) coordinates for current element
    if (!this->getElementCoordinates(Xnod,iEl+1))
    {
      ok = false;
      break;
    }

    // Compute parameter values of the Gauss points over the whole element
    std::array<RealArray,3> gpar, redpar;
    for (int d = 0; d < 3; d++)
    {
      this->getGaussPointParameters(gpar[d],d,nGauss,iEl,xg);
      if (xr)
        this->getGaussPointParameters(redpar[d],d,nRed,iEl,xr);
    }


    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      this->getElementCorners(iEl, fe.XC);

    if (integrand.getIntegrandType() & Integrand::G_MATRIX)
    {
      // Element size in parametric space
      dXidu[0] = el->getParmin(0);
      dXidu[1] = el->getParmin(1);
      dXidu[2] = el->getParmin(2);
    }
    else if (integrand.getIntegrandType() & Integrand::AVERAGE)
    {
      // --- Compute average value of basis functions over the element -----

      fe.Navg.resize(nBasis,true);
      double vol = 0.0;
      for (int k = 0; k < nGauss; k++)
        for (int j = 0; j < nGauss; j++)
          for (int i = 0; i < nGauss; i++)
          {
            // Fetch basis function derivatives at current integration point
            evaluateBasis(fe, dNdu);

            // Compute Jacobian determinant of coordinate mapping
            // and multiply by weight of current integration point
            double detJac = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu,false);
            double weight = 0.125*vol*wg[i]*wg[j]*wg[k];

            // Numerical quadrature
            fe.Navg.add(fe.N,detJac*weight);
            vol += detJac*weight;
      }

      // Divide by element volume
      fe.Navg /= vol;
    }

    else if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
    {
      // Compute the element center
      Go::Point X0;
      double u0 = 0.5*(el->getParmin(0) + el->getParmax(0));
      double v0 = 0.5*(el->getParmin(1) + el->getParmax(1));
      double w0 = 0.5*(el->getParmin(2) + el->getParmax(2));
      lrspline->point(X0,u0,v0,w0);
      X = SplineUtils::toVec3(X0);
    }

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
    if (!integrand.initElement(MNPC[iEl],fe,X,nRed*nRed*nRed,*A))
    {
      A->destruct();
      ok = false;
      break;
    }

    if (xr)
    {
      std::cerr << "Haven't really figured out what this part does yet\n";
      exit(42142);
#if 0
      // --- Selective reduced integration loop ----------------------------

      int ip = (((i3-p3)*nRed*nel2 + i2-p2)*nRed*nel1 + i1-p1)*nRed;
      for (int k = 0; k < nRed; k++, ip += nRed*(nel2-1)*nRed*nel1)
        for (int j = 0; j < nRed; j++, ip += nRed*(nel1-1))
          for (int i = 0; i < nRed; i++, ip++)
          {
            // Local element coordinates of current integration point
            fe.xi   = xr[i];
            fe.eta  = xr[j];
            fe.zeta = xr[k];

            // Parameter values of current integration point
            fe.u = redpar[0](i+1,i1-p1+1);
            fe.v = redpar[1](j+1,i2-p2+1);
            fe.w = redpar[2](k+1,i3-p3+1);

            // Fetch basis function derivatives at current point
            evaluateBasis(fe, 1);

            // Compute Jacobian inverse and derivatives
            fe.detJxW = utl::Jacobian(Jac,fe.dNdX,Xnod,dNdu);

            // Cartesian coordinates of current integration point
            X = Xnod * fe.N;
            X.t = time.t;

            // Compute the reduced integration terms of the integrand
            fe.detJxW *= 0.125*dV*wr[i]*wr[j]*wr[k];
            if (!integrand.reducedInt(*A,fe,X))
              ok = false;
      }
#endif
    }


    // --- Integration loop over all Gauss points in each direction --------

    fe.iGP = iEl*nGauss*nGauss; // Global integration point counter

    Matrix B(p1*p2*p3, 4); // Bezier evaluation points and derivatives
    ig = 1;
    for (int k = 0; k < nGauss; k++)
      for (int j = 0; j < nGauss; j++)
        for (int i = 0; i < nGauss; i++, fe.iGP++, ig++)
        {
          // Local element coordinates of current integration point
          fe.xi   = xg[i];
          fe.eta  = xg[j];
          fe.zeta = xg[k];

          // Parameter values of current integration point
          fe.u = gpar[0][i];
          fe.v = gpar[1][j];
          fe.w = gpar[2][k];

          // Extract bezier basis functions
          B.fillColumn(1, BN.getColumn(ig));
          B.fillColumn(2, BdNdu.getColumn(ig)*2.0/du);
          B.fillColumn(3, BdNdv.getColumn(ig)*2.0/dv);
          B.fillColumn(4, BdNdw.getColumn(ig)*2.0/dw);

          // Fetch basis function derivatives at current integration point
          if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
            evaluateBasis(fe, dNdu, d2Ndu2);
          else
            evaluateBasis(fe, dNdu, C, B) ;

          // look for errors in bezier extraction
          /*
          int N    = nBasis;
          int allP = p1*p2*p3;
          double sum = 0;
          for(int qq=1; qq<=N; qq++) sum+= fe.N(qq);
          if (fabs(sum-1) > 1e-10) {
            std::cerr << "fe.N not sums to one at integration point #" << ig << std::endl;
            exit(123);
          }
          sum = 0;
          for(int qq=1; qq<=N; qq++) sum+= dNdu(qq,1);
          if (fabs(sum) > 1e-10) {
            std::cerr << "dNdu not sums to zero at integration point #" << ig << std::endl;
            exit(123);
          }
          sum = 0;
          for(int qq=1; qq<=N; qq++) sum+= dNdu(qq,2);
          if (fabs(sum) > 1e-10) {
            std::cerr << "dNdv not sums to zero at integration point #" << ig << std::endl;
            exit(123);
          }
          sum = 0;
          for(int qq=1; qq<=N; qq++) sum+= dNdu(qq,3);
          if (fabs(sum) > 1e-10) {
            std::cerr << "dNdw not sums to zero at integration point #" << ig << std::endl;
            exit(123);
          }
          sum = 0;
          for(int qq=1; qq<=allP; qq++) sum+= B(qq,1);
          if (fabs(sum-1) > 1e-10) {
            std::cerr << "Bezier basis not sums to one at integration point #" << ig << std::endl;
            exit(123);
          }
          sum = 0;
          for(int qq=1; qq<=allP; qq++) sum+= B(qq,2);
          if (fabs(sum) > 1e-10) {
            std::cerr << "Bezier derivatives not sums to zero at integration point #" << ig << std::endl;
            exit(123);
          }
          */

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

          // Cartesian coordinates of current integration point
          X   = Xnod * fe.N;
          X.t = time.t;

          // Evaluate the integrand and accumulate element contributions
          fe.detJxW *= 0.125*vol*wg[i]*wg[j]*wg[k];
          if (!integrand.evalInt(*A,fe,time,X))
            ok = false;

    } // end gauss integrand

    // Finalize the element quantities
    if (ok && !integrand.finalizeElement(*A,time,0))
      ok = false;

    // Assembly of global system integral
    if (ok && !glInt.assemble(A->ref(),fe.iel))
      ok = false;

    A->destruct();
  }

  return ok;
}


bool ASMu3D::integrate (Integrand& integrand, int lIndex,
                        GlobalIntegral& glInt,
                        const TimeDomain& time)
{
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu3D::integrate(B)");

  // Get Gaussian quadrature points and weights
  int nGP = integrand.getBouIntegrationPoints(nGauss);
  const double* xg = GaussQuadrature::getCoord(nGP);
  const double* wg = GaussQuadrature::getWeight(nGP);
  if (!xg || !wg) return false;

  // Find the parametric direction of the face normal {-3,-2,-1, 1, 2, 3}
  const int faceDir = (lIndex+1)/(lIndex%2 ? -2 : 2);

  const int t1 = 1 + abs(faceDir)%3; // first tangent direction
  const int t2 = 1 + t1%3;           // second tangent direction

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  LR::parameterEdge edge;
  switch(lIndex)
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
  lrspline->getEdgeElements(edgeElms, (LR::parameterEdge) edge);

  // iterate over all edge elements
  bool ok = true;
  for(LR::Element *el : edgeElms) {
    int iEl = el->getId();
    int nBasis = el->nBasisFunctions();
    FiniteElement fe(nBasis);
    fe.iel = iEl+1;

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
        this->getGaussPointParameters(gpar[d],d,nGP,iEl,xg);

    fe.xi = fe.eta = fe.zeta = faceDir < 0 ? -1.0 : 1.0;
    fe.u = gpar[0](1);
    fe.v = gpar[1](1);
    fe.w = gpar[2](1);

    Matrix dNdu, Xnod, Jac;
    Vec4   X;
    Vec3   normal;
    double dXidu[3];

    // Get element face area in the parameter space
    double dA = this->getParametricArea(iEl+1,abs(faceDir));
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
      this->getElementCorners(iEl,fe.XC);

    if (integrand.getIntegrandType() & Integrand::G_MATRIX)
    {
      // Element size in parametric space
      dXidu[0] = el->getParmax(0) - el->getParmin(0);
      dXidu[1] = el->getParmax(1) - el->getParmin(1);
      dXidu[2] = el->getParmax(2) - el->getParmin(2);
    }

    // Initialize element quantities
    LocalIntegral* A = integrand.getLocalIntegral(nBasis,fe.iel,true);
    if (!integrand.initElementBou(MNPC[iEl],*A))
    {
      A->destruct();
      ok = false;
      break;
    }

    // --- Integration loop over all Gauss points in each direction --------

    fe.iGP = firstp; // Global integration point counter
    int k1,k2,k3;
    for (int j = 0; j < nGP; j++)
      for (int i = 0; i < nGP; i++, fe.iGP++)
      {
        // Local element coordinates and parameter values
        // of current integration point
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
          fe.u = gpar[0](k1+1);
        }
        if (gpar[1].size() > 1)
        {
          fe.eta = xg[k2];
          fe.v = gpar[1](k2+1);
        }
        if (gpar[2].size() > 1)
        {
          fe.zeta = xg[k3];
          fe.w = gpar[2](k3+1);
        }

        // Fetch basis function derivatives at current integration point
        evaluateBasis(fe, dNdu);

        // Compute basis function derivatives and the face normal
        fe.detJxW = utl::Jacobian(Jac,normal,fe.dNdX,Xnod,dNdu,t1,t2);
        if (fe.detJxW == 0.0) continue; // skip singular points

        if (faceDir < 0) normal *= -1.0;

        // Compute G-matrix
        if (integrand.getIntegrandType() & Integrand::G_MATRIX)
          utl::getGmat(Jac,dXidu,fe.G);

        // Cartesian coordinates of current integration point
        X = Xnod * fe.N;
        X.t = time.t;

        // Evaluate the integrand and accumulate element contributions
        fe.detJxW *= 0.25*dA*wg[i]*wg[j];
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

    firstp += nGP*nGP;
  }

  return ok;
}


bool ASMu3D::integrateEdge (Integrand& integrand, int lEdge,
                            GlobalIntegral& glInt,
                            const TimeDomain& time)
{
  std::cerr << "ASMu3D::integrateEdge(...) is not properly implemented yet :(" << std::endl;
  exit(776654);
#if 0
  if (!lrspline) return true; // silently ignore empty patches

  PROFILE2("ASMu3D::integrate(E)");

  // Get Gaussian quadrature points and weights
  const double* xg = GaussQuadrature::getCoord(nGauss);
  const double* wg = GaussQuadrature::getWeight(nGauss);
  if (!xg || !wg) return false;

  // Compute parameter values of the Gauss points along the whole patch edge
  std::array<Matrix,3> gpar;
  for (int d = 0; d < 3; d++)
    if (lEdge < d*4+1 || lEdge >= d*4+5)
    {
      gpar[d].resize(1,1);
      if (lEdge%4 == 1)
  gpar[d].fill(lrspline->startparam(d));
      else if (lEdge%4 == 0)
  gpar[d].fill(lrspline->endparam(d));
      else if (lEdge == 6 || lEdge == 10)
  gpar[d].fill(d == 0 ? lrspline->endparam(d) : lrspline->startparam(d));
      else if (lEdge == 2 || lEdge == 11)
  gpar[d].fill(d == 1 ? lrspline->endparam(d) : lrspline->startparam(d));
      else if (lEdge == 3 || lEdge == 7)
  gpar[d].fill(d == 2 ? lrspline->endparam(d) : lrspline->startparam(d));
    }
    else
    {
      int pm1 = lrspline->order(d) - 1;
      RealArray::const_iterator uit = lrspline->basis(d).begin() + pm1;
      double ucurr, uprev = *(uit++);
      int nCol = lrspline->numCoefs(d) - pm1;
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
  std::vector<Go::BasisDerivs> spline;
  {
    PROFILE2("Spline evaluation");
    lrspline->computeBasisGrid(gpar[0],gpar[1],gpar[2],spline);
  }

  const int n1 = lrspline->numCoefs(0);
  const int n2 = lrspline->numCoefs(1);
  const int n3 = lrspline->numCoefs(2);

  const int p1 = lrspline->order(0);
  const int p2 = lrspline->order(1);
  const int p3 = lrspline->order(2);

  std::map<char,size_t>::const_iterator iit = firstBp.find(lEdge);
  size_t firstp = iit == firstBp.end() ? 0 : iit->second;

  FiniteElement fe(p1*p2*p3);
  fe.u = gpar[0](1,1);
  fe.v = gpar[1](1,1);
  fe.w = gpar[2](1,1);
  if (gpar[0].size() == 1) fe.xi = fe.u == lrspline->startparam(0) ? -1.0 : 1.0;
  if (gpar[1].size() == 1) fe.eta = fe.v == lrspline->startparam(1) ? -1.0 : 1.0;
  if (gpar[2].size() == 1) fe.zeta = fe.w == lrspline->startparam(2) ? -1.0 : 1.0;

  Matrix dNdu, Xnod, Jac;
  Vec4   X;
  Vec3   tang;


  // === Assembly loop over all elements on the patch edge =====================

  int iel = 1;
  for (int i3 = p3; i3 <= n3; i3++)
    for (int i2 = p2; i2 <= n2; i2++)
      for (int i1 = p1; i1 <= n1; i1++, iel++)
      {
  fe.iel = MLGE[iel-1];
  if (fe.iel < 1) continue; // zero-volume element

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
  int ip = MNPC[iel-1].back();
#ifdef INDEX_CHECK
  if (ip < 0 || (size_t)ip >= nodeInd.size()) return false;
#endif
  if (lEdge < 5)
  {
    dS = lrspline->knotSpan(0,nodeInd[ip].I);
    ip = (i1-p1)*nGauss;
  }
  else if (lEdge < 9)
  {
    dS = lrspline->knotSpan(1,nodeInd[ip].J);
    ip = (i2-p2)*nGauss;
  }
  else if (lEdge < 13)
  {
    dS = lrspline->knotSpan(2,nodeInd[ip].K);
    ip = (i3-p3)*nGauss;
  }

  // Set up control point coordinates for current element
  if (!this->getElementCoordinates(Xnod,iel)) return false;

  // Initialize element quantities
  LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
  bool ok = integrand.initElementBou(MNPC[iel-1],*A);


  // --- Integration loop over all Gauss points along the edge -----------

  fe.iGP = firstp + ip; // Global integration point counter

  for (int i = 0; i < nGauss && ok; i++, ip++, fe.iGP++)
  {
    // Parameter values of current integration point
    if (gpar[0].size() > 1) fe.u = gpar[0](i+1,i1-p1+1);
    if (gpar[1].size() > 1) fe.v = gpar[1](i+1,i2-p2+1);
    if (gpar[2].size() > 1) fe.w = gpar[2](i+1,i3-p3+1);

    // Fetch basis function derivatives at current integration point
    evaluateBasis(fe, dNdu);

    // Compute basis function derivatives and the edge tang
    fe.detJxW = utl::Jacobian(Jac,tang,fe.dNdX,Xnod,dNdu,1+(lEdge-1)/4);
    if (fe.detJxW == 0.0) continue; // skip singular points

    // Cartesian coordinates of current integration point
    X = Xnod * fe.N;
    X.t = time.t;

    // Evaluate the integrand and accumulate element contributions
    fe.detJxW *= 0.5*dS*wg[i];
          ok = integrand.evalBou(*A,fe,time,X,tang);
  }

  // Assembly of global system integral
  if (ok && !glInt.assemble(A->ref(),fe.iel))
    ok = false;

  A->destruct();

  if (!ok) return false;
  }

  return true;
#endif
  return false;
}


int ASMu3D::evalPoint (const double* xi, double* param, Vec3& X) const
{
  std::cerr << "ASMu3D::evalPoint(...) is not properly implemented yet :(" << std::endl;
  exit(776654);
#if 0
  if (!lrspline) return -3;

  int i;
  for (i = 0; i < 3; i++)
    param[i] = (1.0-xi[i])*lrspline->startparam(i) + xi[i]*lrspline->endparam(i);

  Go::Point X0;
  lrspline->point(X0,param[0],param[1],param[2]);
  for (i = 0; i < 3 && i < lrspline->dimension(); i++)
    X[i] = X0[i];

  // Check if this point matches any of the control points (nodes)
  Vec3 Xnod;
  size_t inod = 1;
  RealArray::const_iterator cit = lrspline->coefs_begin();
  for (i = 0; cit != lrspline->coefs_end(); cit++, i++)
  {
    if (i < 3) Xnod[i] = *cit;
    if (i+1 == lrspline->dimension())
      if (X.equal(Xnod,0.001))
  return inod;
      else
      {
  inod++;
  i = -1;
      }
  }

  return 0;
#endif
  return -1;
}


bool ASMu3D::getGridParameters (RealArray& prm, int dir, int nSegPerSpan) const
{
  if (!lrspline) return false;

  // output is written once for each element resulting in a lot of unnecessary storage
  // this is preferable to figuring out all element topology information

  for(LR::Element *el : lrspline->getAllElements() ) {
    // evaluate element at element corner points
    double umin = el->umin();
    double umax = el->umax();
    double vmin = el->vmin();
    double vmax = el->vmax();
    double wmin = el->wmin();
    double wmax = el->wmax();
    for(int iw=0; iw<=nSegPerSpan; iw++) {
      for(int iv=0; iv<=nSegPerSpan; iv++) {
        for(int iu=0; iu<=nSegPerSpan; iu++) {
          double u = umin + (umax-umin)/nSegPerSpan*iu;
          double v = vmin + (vmax-vmin)/nSegPerSpan*iv;
          double w = wmin + (wmax-wmin)/nSegPerSpan*iw;
          if (dir==0)
            prm.push_back(u);
          else if (dir==1)
            prm.push_back(v);
          else
            prm.push_back(w);
        }
      }
    }
  }
  return true;
}

bool ASMu3D::tesselate (ElementBlock& grid, const int* npe) const
{
  if (!lrspline) return false;

  if (npe[0] != npe[1] || npe[0] != npe[2]) {
    std::cerr << "ASMu2D::tesselate does not support different tesselation resolution in "
              << "u- and v-direction. nviz u = " << npe[0] << ", nviz v = " << npe[1]
              << ", nviz w = " << npe[2] << std::endl;
    return false;
  }

  int nNodesPerElement =  npe[0]   * npe[1]   * npe[2];
  int nSubElPerElement = (npe[0]-1)*(npe[1]-1)*(npe[2]-1);
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
    double wmin = (**el).wmin();
    double wmax = (**el).wmax();
    for(int iw=0; iw<npe[2]; iw++) {
      for(int iv=0; iv<npe[1]; iv++) {
        for(int iu=0; iu<npe[0]; iu++) {
          double u = umin + (umax-umin)/(npe[0]-1)*iu;
          double v = vmin + (vmax-vmin)/(npe[1]-1)*iv;
          double w = wmin + (wmax-wmin)/(npe[2]-1)*iw;
          Go::Point pt;
          lrspline->point(pt, u,v,w, iel, iu!=npe[0]-1, iv!=npe[1]-1, iw!=npe[2]-1);
          for(int dim=0; dim<nsd; dim++)
            grid.setCoor(inod, dim, pt[dim]);
          inod++;
        }
      }
    }
  }

  int ip = 0;
  iel = 0;
  for(int i=0; i<lrspline->nElements(); i++) {
    int iStart = i*nNodesPerElement;
    for(int iw=0; iw<npe[2]-1; iw++) {
      for(int iv=0; iv<npe[1]-1; iv++) {
        for(int iu=0; iu<npe[0]-1; iu++, iel++) {
          // enumerate nodes counterclockwise around the hex
          grid.setNode(ip++, iStart + (iw  )*npe[0]*npe[1] + (iv  )*npe[0] + (iu  ) );
          grid.setNode(ip++, iStart + (iw  )*npe[0]*npe[1] + (iv  )*npe[0] + (iu+1) );
          grid.setNode(ip++, iStart + (iw  )*npe[0]*npe[1] + (iv+1)*npe[0] + (iu+1) );
          grid.setNode(ip++, iStart + (iw  )*npe[0]*npe[1] + (iv+1)*npe[0] + (iu  ) );
          grid.setNode(ip++, iStart + (iw+1)*npe[0]*npe[1] + (iv  )*npe[0] + (iu  ) );
          grid.setNode(ip++, iStart + (iw+1)*npe[0]*npe[1] + (iv  )*npe[0] + (iu+1) );
          grid.setNode(ip++, iStart + (iw+1)*npe[0]*npe[1] + (iv+1)*npe[0] + (iu+1) );
          grid.setNode(ip++, iStart + (iw+1)*npe[0]*npe[1] + (iv+1)*npe[0] + (iu  ) );
          grid.setElmId(iel+1, i+1);
        }
      }
    }
  }

  return true;
}


bool ASMu3D::evalSolution (Matrix& sField, const Vector& locSol,
                           const int* npe) const
{
  // Compute parameter values of the result sampling points
  std::array<RealArray,3> gpar;
  for (int dir = 0; dir < 3; dir++)
    if (!this->getGridParameters(gpar[dir],dir,npe[dir]-1))
      return false;

  // Evaluate the primary solution at all sampling points
  return this->evalSolution(sField,locSol,gpar.data());
}


bool ASMu3D::evalSolution (Matrix& sField, const Vector& locSol,
                           const RealArray* gpar, bool, int deriv) const
{
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

    // Set up control point (nodal) coordinates for current element
    if (deriv > 0 && !this->getElementCoordinates(Xnod,iel+1))
      return false;

    // Evaluate basis function values/derivatives at current parametric point
    // and mulitiply with control point values to get the point-wise solution
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
  if (npe == nullptr || npe[0] != npe[1] || npe[0] != npe[2] || npe[1] != npe[2])
  {
    std::cerr <<" *** ASMu2D::evalSolution: LR B-splines require the"
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
    for (int i = 0; i < v->nBasisFunctions(); i++)
      sField.fillColumn(i+1,&(*v->getBasisfunction(i)->cp()));
    delete v;
    return true;
  }

  std::cerr <<" *** ASMu2D::evalSolution: Failure!"<< std::endl;
  return false;
}


bool ASMu3D::evalSolution (Matrix& sField, const IntegrandBase& integrand,
                           const RealArray* gpar, bool) const
{
#ifdef SP_DEBUG
  std::cout <<"ASMu3D::evalSolution(Matrix&,const IntegrandBase&,const RealArray*,bool)\n";
#endif

  sField.resize(0,0);

  // TODO: investigate the possibility of doing "regular" refinement by
  //       uniform tesselation grid and ignoring LR mesh lines

  size_t nPoints = gpar[0].size();
  bool use2ndDer = integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES;
  if (nPoints != gpar[1].size() || nPoints != gpar[2].size())
    return false;

  Vector   solPt;
  Matrix   dNdu, Jac, Xnod;
  Matrix3D d2Ndu2, Hess;

  // Evaluate the secondary solution field at each point
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

    // Evaluate the basis functions at current parametric point
    FiniteElement fe(lrspline->getElement(iel)->nBasisFunctions());
    fe.u   = gpar[0][i];
    fe.v   = gpar[1][i];
    fe.w   = gpar[2][i];
    fe.iel = iel+1;
    if (use2ndDer)
      evaluateBasis(fe, dNdu, d2Ndu2);
    else
      evaluateBasis(fe, dNdu);

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


std::vector<int> ASMu3D::getFaceNodes (int face, int basis) const
{
  size_t ofs = 1;
  for (int i = 1; i < basis; i++)
    ofs += this->getNoNodes(i);

  std::vector<LR::Basisfunction*> edgeFunctions;
  this->getBasis(basis)->getEdgeFunctions(edgeFunctions,
                                          static_cast<LR::parameterEdge>(face));

  std::vector<int> result(edgeFunctions.size());
  std::transform(edgeFunctions.begin(), edgeFunctions.end(), result.begin(),
                 [ofs](LR::Basisfunction* a) { return a->getId()+ofs; });

  return result;
}


void ASMu3D::getBoundaryNodes (int lIndex, IntVec& nodes, int basis, int, bool) const
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
  case 5: edge = LR::BOTTOM; break;
  case 6: edge = LR::TOP; break;
  default: return;
  }

  nodes = this->getFaceNodes(edge, basis);

#if SP_DEBUG > 1
  std::cout <<"Boundary nodes in patch "<< idx+1 <<" edge "<< lIndex <<":";
  for (size_t i = 0; i < nodes.size(); i++)
    std::cout <<" "<< nodes[i];
  std::cout << std::endl;
#endif
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
  int edge = LR::NONE;
  if (I > 0) edge |= LR::EAST;
  else       edge |= LR::WEST;
  if (J > 0) edge |= LR::NORTH;
  else       edge |= LR::SOUTH;
  if (K > 0) edge |= LR::TOP;
  else       edge |= LR::BOTTOM;

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
