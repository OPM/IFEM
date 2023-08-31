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
#include "GlbL2projector.h"
#include "SparseMatrix.h"
#include "SplineFields1D.h"
#include "ElementBlock.h"
#include "SplineUtils.h"
#include "Utilities.h"
#include "Function.h"
#include "Vec3Oper.h"
#include "Tensor.h"
#include "IFEM.h"
#include <numeric>


ASMs1D::ASMs1D (unsigned char n_s, unsigned char n_f)
  : ASMstruct(1,n_s,n_f), elmCS(myCS), nodalT(myT)
{
  curv = proj = nullptr;
  updatedT = false;
}


ASMs1D::ASMs1D (const ASMs1D& patch, unsigned char n_f)
  : ASMstruct(patch,n_f), elmCS(patch.myCS), nodalT(patch.myT)
{
  curv = patch.curv;
  proj = patch.proj;
  updatedT = false;

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
    std::cerr <<"  ** ASMs1D::read: The dimension of this curve patch "
	      << curv->dimension() <<" is less than nsd="<< (int)nsd
	      <<".\n                   Resetting nsd to "<< curv->dimension()
	      <<" for this patch."<< std::endl;
    nsd = curv->dimension();
  }

  geomB = curv;
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
    if (proj && proj != curv) delete proj;
    if (curv && !shareFE) delete curv;
    geomB = projB = curv = proj = nullptr;
  }

  // Erase the FE data
  this->ASMbase::clear(retainGeometry);
}


size_t ASMs1D::getNoNodes (int basis) const
{
  size_t n = this->ASMbase::getNoNodes(basis);
  if (n > 0 || basis < 1 || !curv) return n;

  // We request the number of nodes before the FE topology has been generated
  return curv->numCoefs();
}


bool ASMs1D::refine (const LR::RefineData& prm, Vectors&)
{
  if (!curv) return false;

  if (shareFE && !prm.refShare)
  {
    nnod = curv->numCoefs();
    return true;
  }
  else if (prm.elements.empty())
    return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = curv->basis().begin();
  for (size_t i = curv->order(); i <= nnod; i++)
    if (uit[i-1] < uit[i])
      if (std::find(prm.elements.begin(),prm.elements.end(),
                    MLGE[i-curv->order()]-1) != prm.elements.end())
      {
        extraKnots.push_back(0.5*(uit[i-1] + uit[i]));
        if (prm.elements.size() < 100)
          IFEM::cout <<"Refining element "<< MLGE[i-curv->order()]
                     <<"\tu="<< extraKnots.back()
                     <<"\th="<< uit[i] - uit[i-1] << std::endl;
      }

  curv->insertKnot(extraKnots);
  if (proj != curv)
    proj->insertKnot(extraKnots);

  return true;
}


bool ASMs1D::refine (const RealArray& xi)
{
  if (!curv || xi.empty()) return false;
  if (xi.front() < 0.0 || xi.back() > 1.0) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = curv->basis().begin();
  double uprev = *(uit++);
  while (uit != curv->basis().end())
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

  curv->insertKnot(extraKnots);
  return true;
}


bool ASMs1D::uniformRefine (int nInsert)
{
  if (!curv || nInsert < 1) return false;
  if (shareFE) return true;

  RealArray extraKnots;
  RealArray::const_iterator uit = curv->basis().begin();
  double uprev = *(uit++);
  while (uit != curv->basis().end())
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


/*!
  This method is supposed to be invoked twice during the model generation.
  In the first call, with \a init = \e true, the spline curve object \a *curv
  is cloned into \a *proj and the two pointers are then swapped, such that
  the subsequent refine and raiseOrder operations will apply to the projection
  basis and not on the geometry basis.
  In the second call, the pointers are swapped back.

  The method can also be invoked twice with \a init = \e false in case the
  projection basis is to be read from a file.
*/

bool ASMs1D::createProjectionBasis (bool init)
{
  if (!curv)
    return false;
  else if (init && !proj)
    projB = proj = curv->clone();

  std::swap(geomB,projB);
  std::swap(curv,proj);
  return true;
}


bool ASMs1D::generateFEMTopology ()
{
  return this->generateOrientedFEModel(Vec3());
}


bool ASMs1D::generateOrientedFEModel (const Vec3& Zaxis)
{
  if (!curv) return false;

  if (!proj)
    proj = curv;
  else
  {
    RealArray simple1, simple2;
    curv->basis().knotsSimple(simple1);
    proj->basis().knotsSimple(simple2);
    if (simple1 != simple2)
    {
      std::cerr <<" *** FE basis and projection basis must have the same knots."
                << std::endl;
      return false;
    }
  }

  int n1 = curv->numCoefs();
  int p1 = curv->order();

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

  int pgEl = gEl;
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

  if (proj != curv)
  {
    n1 = proj->numCoefs();
    p1 = proj->order();
    projMLGE.resize(n1-p1+1,0);
    projMNPC.resize(projMLGE.size());
    int nnod_p = 0, nel_p = 0;
    for (int i1 = 1; i1 <= n1; i1++)
    {
      if (i1 >= p1)
      {
        if (*(proj->basis().begin()+i1) > *(proj->basis().begin()+i1-1))
        {
          projMLGE[nel_p] = ++pgEl;
          projMNPC[nel_p].resize(p1,0);

          int lnod = 0;
          for (int j1 = p1-1; j1 >= 0; j1--)
            projMNPC[nel_p][lnod++] = nnod_p - j1;
        }

        ++nel_p;
      }
      ++nnod_p;
    }
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
      Vec3 X2 = this->getCoord(1+MNPC[i].back());
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
      Vec3 X2 = this->getCoord(1+MNPC[i].back());
      double alpha = twist(0.5*(X1+X2)); // twist angle in the element mid-point
      myCS[i] *= Tensor(alpha*M_PI/180.0,1); // rotate about local X-axis
#ifdef SP_DEBUG
      std::cout <<"Twisted axes for beam element "<< MLGE[i]
                <<", from "<< X1 <<" to "<< X2 <<":\n"<< myCS[i];
#endif
    }

  return true;
}


bool ASMs1D::connectPatch (int vertex, ASM1D& neighbor, int nvertex, int thick)
{
  ASMs1D* neighS = dynamic_cast<ASMs1D*>(&neighbor);
  if (!neighS)
    return false;

  int slave  = vertex  == 1 ? 1 : 1+this->getSize()-thick;
  int master = nvertex == 1 ? 1 : 1+neighS->getSize()-thick;
  if (!this->connectBasis(*neighS,slave,master,thick))
    return false;

  this->addNeighbor(neighS);
  return true;
}


bool ASMs1D::connectBasis (ASMs1D& neighbor, int slave, int master, int thick)
{
  if (shareFE && neighbor.shareFE)
    return true;
  else if (shareFE || neighbor.shareFE)
  {
    std::cerr <<" *** ASMs1D::connectBasis: Logic error, cannot"
	      <<" connect a sharedFE patch with an unshared one"<< std::endl;
    return false;
  }

  const double xtol = 1.0e-4;
  for (int i = 0; i < thick; i++, slave++, master++)
    if (!neighbor.getCoord(master).equal(this->getCoord(slave),xtol))
    {
      std::cerr <<" *** ASMs1D::connectBasis: Non-matching nodes "
                << master <<": "<< neighbor.getCoord(master)
                <<"\n                                          and "
                << slave <<": "<< this->getCoord(slave) << std::endl;
      return false;
    }
    else
      ASMbase::collapseNodes(neighbor,master,*this,slave);

  return true;
}


void ASMs1D::closeBoundaries (int, int basis, int master)
{
  if (basis < 1) basis = 1;
  int n1 = this->getSize(basis);
  this->makePeriodic(1,master+n1-1);
}


int ASMs1D::constrainNode (double xi, int dof, int code)
{
  if (xi < 0.0 || xi > 1.0)
    return 0;

  int n1 = this->getSize();
  int node = xi > 0.0 ? 1+int(0.5+(n1-1)*xi) : 1;

  this->prescribe(node,dof,code);

  return node;
}


size_t ASMs1D::constrainEndLocal (int dir, int dof, int code)
{
  if (nf < 2 || this->allDofs(dof))
  {
    // If all DOFs are constrained, local axis directions are irrelevant
    this->constrainNode(dir > 0 ? 1.0 : 0.0, dof, code);
    return 0;
  }
  else if (shareFE == 'F')
  {
    std::cerr <<"\n *** ASMs1D::constrainEndLocal: Logic error, can not have"
              <<" constraints in local axes for shared patches."<< std::endl;
    return 0;
  }

  // We need an extra node representing the local (master) DOFs at this point
  int iMnod = myMLGN.size();
  int iSnod = dir > 0 ? this->getSize()-1 : 0;

  // Create an extra node for the local DOFs. The new node, for which
  // the Dirichlet boundary conditions will be defined, then inherits
  // the global node number of the original node. The original node, which
  // do not enter the equation system, receives a new global node number.
  std::map<int,int>::const_iterator xit = xNode.find(MLGN[iSnod]);
  if (xit != xNode.end())
  {
    // This node has already been processed by another patch
    myMLGN.push_back(xit->second);
    return 0;
  }

  // Store the original-to-extra node number mapping in ASMstruct::xNode
  myMLGN.push_back(++gNod);
  xNode[MLGN[iSnod]] = gNod;
  std::swap(myMLGN[iMnod],myMLGN[iSnod]);

  xnMap[1+iMnod] = 1+iSnod; // Store nodal connection needed by getCoord
  nxMap[1+iSnod] = 1+iMnod; // Store nodal connection needed by getNodeID

  // Add Dirichlet condition on the local DOF(s) of the added node
  this->prescribe(1+iMnod,dof,code);

  // Find the local-to-global transformation matrix at this end point,
  // where the X-axis corresponds to the tangent direction of the curve,
  // and establish constraint equations relating the global and local DOFs
  double uEnd = dir > 0 ? curv->endparam() : curv->startparam();
  this->addLocal2GlobalCpl(iSnod,MLGN[iMnod],this->getLocal2Global(uEnd));

  return 1;
}


Tensor ASMs1D::getLocal2Global (double u) const
{
  if (nsd < 2)
    return Tensor(1,true);

  // Find the local-to-global transformation matrix at this point,
  // where the X-axis corresponds to the tangent direction of the curve
  Tensor Tlg(nsd);
  if (nsd == 2)
  {
    std::vector<Go::Point> pts(2);
    curv->point(pts,u,1);
    double dXlen = pts[1].length();
    Tlg(1,1) = pts[1][0] / dXlen;
    Tlg(2,1) = pts[1][1] / dXlen;
    Tlg(1,2) = -Tlg(2,1);
    Tlg(2,2) =  Tlg(1,1);
  }
  else
  {
    std::vector<Go::Point> pts(3);
    curv->point(pts,u,2);
    Vec3 dX(SplineUtils::toVec3(pts[1]));
    if (pts[2].length2() < 1.0e-12) // zero curvature ==> straight line
      Tlg = Tensor(dX,true);
    else
    {
      // Let dX define the local X-axis
      // and the bi-normal vector define a point in the local XZ-plane
      Vec3 biNormal(dX,SplineUtils::toVec3(pts[2]));
      Tlg = Tensor(dX,biNormal,false,true);
    }
  }

#if SP_DEBUG > 2
  std::cout <<"Local-to-global transformation at u="<< u <<":\n"<< Tlg;
#endif
  return Tlg;
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
  if (inod == 0) return Vec3();

  std::map<size_t,XYZ>::const_iterator it = myRmaster.find(inod);
  if (it != myRmaster.end()) return Vec3(it->second.data());

  int ip = (inod-1)*curv->dimension();
  if (ip < 0) return Vec3();

  return Vec3(&(*(curv->coefs_begin()+ip)),nsd);
}


Tensor ASMs1D::getRotation (size_t inod) const
{
  return inod < 1 || inod > nodalT.size() ? Tensor(nsd,true) : nodalT[inod-1];
}


bool ASMs1D::getElementCoordinates (Matrix& X, int iel, bool) const
{
  return this->getElementCoordinates(X,iel,MNPC,curv);
}


bool ASMs1D::getElementCoordinates (Matrix& X, int iel, const IntMat& mnpc,
                                    const Go::SplineCurve* crv) const
{
#ifdef INDEX_CHECK
  if (iel < 1 || (size_t)iel > mnpc.size())
  {
    std::cerr <<" *** ASMs1D::getElementCoordinates: Element index "<< iel
	      <<" out of range [1,"<< mnpc.size() <<"]."<< std::endl;
    return false;
  }
#endif

  X.resize(nsd,crv->order());

  RealArray::const_iterator cit = crv->coefs_begin();
  for (size_t n = 0; n < X.cols(); n++)
  {
    int ip = mnpc[iel-1][n]*crv->dimension();
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

  size_t nno = curv->numCoefs();
  if (displ.size() != nsd*nno && displ.size() != nsd*MLGN.size())
  {
    std::cerr <<" *** ASMs1D::updateCoords: Invalid dimension "
              << displ.size() <<" on displacement vector, should be ";
    if (nno != MLGN.size())
      std::cerr <<"either "<< nsd*MLGN.size() <<" or ";
    std::cerr << nsd*nno << std::endl;
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

  updatedT = true;
  return true;
}


Vec3 ASMs1D::getElementCenter (int iel) const
{
  Vec3 X0;
  if (curv)
  {
    double u[2];
    this->getElementBorders(iel,u);
    SplineUtils::point(X0,0.5*(u[0]+u[1]),curv);
  }
  return X0;
}


void ASMs1D::getBoundaryNodes (int lIndex, IntVec& nodes,
                               int, int thick, int, bool local) const
{
  if (!curv) return; // silently ignore empty patches

  size_t iel = lIndex == 1 ? 0 : nel-1;
  if (MLGE[iel] > 0 && (lIndex == 1 || lIndex == 2))
  {
    int offset = lIndex == 1 ? 0 : curv->order() - thick;
    for (int i = 0; i < thick; i++)
    {
      int node = MNPC[iel][offset+i];
      nodes.push_back(local ? node+1 : MLGN[node]);
    }
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
    return std::make_pair(curv->numCoefs(),dist);

  // We are inside, now find which knot-span we are in and find closest node
  RealArray::iterator it = std::lower_bound(curv->basis().begin(),
                                            curv->basis().end(),param);
  size_t mnod = it - curv->basis().begin();
  size_t jnod = 0;
  double dmin = 0.0;
  for (size_t inod = mnod-curv->order(); inod < mnod; inod++)
  {
    RealArray::const_iterator p = curv->coefs_begin() + inod*curv->dimension();
    double d2 = Go::Point(p,p+curv->dimension()).dist2(Xfound);
    if (d2 < dmin || jnod == 0)
    {
      jnod = inod+1;
      dmin = d2;
    }
  }

#ifdef SP_DEBUG
  std::cout <<"ASMs1D::findClosestNode("<< X
            <<"): Found "<< Xfound <<" at u="<< param
            <<" inod="<< jnod <<" distance="<< sqrt(dmin) << std::endl;
#endif
  return std::make_pair(jnod,sqrt(dmin));
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


size_t ASMs1D::getNoProjectionNodes () const
{
  if (!proj) return 0;

  return proj->numCoefs();
}


bool ASMs1D::getParameterDomain (Real2DMat& u, IntVec* corners) const
{
  u.resize(1,RealArray(2));
  u.front().front() = curv->basis().startparam();
  u.front().back() = curv->basis().endparam();

  if (corners)
  {
    corners->resize(2);
    corners->front() = 1;
    corners->back() = curv->numCoefs();
  }

  return true;
}


const Vector& ASMs1D::getGaussPointParameters (Matrix& uGP, int nGauss,
                                               const double* xi,
                                               const Go::SplineCurve* crv,
                                               bool skipNullSpans) const
{
  if (!crv) crv = curv;

  int pm1 = crv->order() - 1;
  RealArray::const_iterator uit = crv->basis().begin() + pm1;

  int nCol = crv->numCoefs() - pm1;
  uGP.resize(nGauss,nCol);

  int iel = 0;
  double uprev = *(uit++);
  for (int j = 1; j <= nCol; ++uit, j++)
  {
    double ucurr = *uit;
    if (!skipNullSpans || ucurr > uprev)
    {
      ++iel;
      for (int i = 1; i <= nGauss; i++)
        uGP(i,iel) = 0.5*((ucurr-uprev)*xi[i-1] + ucurr+uprev);
    }
    uprev = ucurr;
  }

  if (iel < nCol)
    uGP.resize(nGauss,iel);

  return uGP;
}


void ASMs1D::getElementBorders (int iel, double* ub) const
{
  RealArray::const_iterator uit = curv->basis().begin();

  // Fetch parameter values of the element ends (knots)
  int i = iel-1 + curv->order();
  ub[0] = uit[i-1];
  ub[1] = uit[i];
}


double ASMs1D::getElementEnds (int i, Vec3Vec& XC) const
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
    XC.push_back(Vec3(pt,nsd));

  // Calculate the characteristic element size
  double h = getElementSize(XC);
  if (elmCS.empty()) return h;

  // Add the local Z-axis as the third vector
  int iel = i - curv->order();
  XC.push_back(elmCS[iel][2]);
  return h;
}


void ASMs1D::evaluateBasis (double u, double, double, Vector& N) const
{
  this->extractBasis(u,N);
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

  N.resize(p1);
  dNdu.resize(p1,1);
  d2Ndu2.resize(p1,1,1);

  RealArray basisDerivs, basisDerivs2;
  curv->computeBasis(u,N,basisDerivs,basisDerivs2);
  dNdu.fillColumn(1,basisDerivs);
  d2Ndu2.fillColumn(1,1,basisDerivs2);
}


void ASMs1D::extractBasis (double u, Vector& N, Matrix& dNdu,
                           Matrix3D& d2Ndu2, Matrix4D& d3Ndu3) const
{
  int p1 = curv->order();

  N.resize(p1);
  dNdu.resize(p1,1);
  d2Ndu2.resize(p1,1,1);
  d3Ndu3.resize(p1,1,1,1);

  RealArray basisDerivs, basisDerivs2, basisDerivs3;
  curv->computeBasis(u,N,basisDerivs,basisDerivs2,basisDerivs3);
  dNdu.fillColumn(1,basisDerivs);
  d2Ndu2.fillColumn(1,1,basisDerivs2);
  d3Ndu3.fillColumn(1,1,1,basisDerivs3);
}


bool ASMs1D::integrate (Integrand& integrand,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
  if (!curv) return true; // silently ignore empty patches

  const int p1 = curv->order();

  // Get Gaussian quadrature points and weights
  const int     ng = this->getNoGaussPt(p1);
  const double* xg = GaussQuadrature::getCoord(ng);
  const double* wg = GaussQuadrature::getWeight(ng);
  if (!xg || !wg) return false;

  // Get the reduced integration quadrature points, if needed
  const double* xr = nullptr;
  const double* wr = nullptr;
  int nRed = integrand.getReducedIntegration(ng);
  if (nRed > 0)
  {
    xr = GaussQuadrature::getCoord(nRed);
    wr = GaussQuadrature::getWeight(nRed);
    if (!xr || !wr) return false;
  }
  else if (nRed < 0)
    nRed = ng; // The integrand needs to know nGauss

  // Compute parameter values of the Gauss points over the whole patch
  Matrix gpar, redpar;
  this->getGaussPointParameters(gpar,ng,xg);
  if (xr)
    this->getGaussPointParameters(redpar,nRed,xr);

  FiniteElement fe(p1);
  Matrix   dNdu, Jac;
  Matrix3D d2Ndu2, Hess;
  Matrix4D d3Ndu3;
  double   param[3] = { 0.0, 0.0, 0.0 };
  Vec4     X(param);

  if (nsd > 1 && (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES))
    fe.G.resize(nsd,2); // For storing d{X}/du and d2{X}/du2


  // === Assembly loop over all elements in the patch ==========================

  bool ok = true;
  for (size_t iel = 0; iel < nel && ok; iel++)
  {
    fe.iel = MLGE[iel];
    if (fe.iel < 1) continue; // zero-length element

#ifdef SP_DEBUG
    int ielm = 1+iel;
    if (dbgElm < 0 && ielm != -dbgElm)
      continue; // Skipping all elements, except for -dbgElm
#endif

    LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel);
    if (!A) continue; // no integrand contributions for this element

    // Check that the current element has nonzero length
    double dL = 0.5*this->getParametricLength(1+iel);
    if (dL < 0.0) ok = false; // topology error (probably logic error)

    // Set up control point coordinates for current element
    ok &= this->getElementCoordinates(fe.Xn,1+iel);

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
      fe.h = this->getElementEnds(p1+iel,fe.XC);

    if (integrand.getIntegrandType() & Integrand::NODAL_ROTATIONS)
    {
      this->getElementNodalRotations(fe.Tn,iel);
      if (!elmCS.empty()) fe.Te = elmCS[iel];
    }

    if (integrand.getIntegrandType() & Integrand::ELEMENT_CENTER)
    {
      // Compute the element center
      param[0] = 0.5*(gpar(1,1+iel) + gpar(ng,1+iel));
      SplineUtils::point(X,param[0],curv);
    }

    // Initialize element matrices
    ok &= integrand.initElement(MNPC[iel],fe,X,nRed,*A);

    if (xr)
    {
      // --- Selective reduced integration loop --------------------------------

      for (int i = 0; i < nRed && ok; i++)
      {
	// Local element coordinates of current integration point
	fe.xi = xr[i];

        // Parameter values of current integration point
        fe.u = param[0] = redpar(1+i,1+iel);

        if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
          this->extractBasis(fe.u,fe.N);
        else
        {
          // Fetch basis function derivatives at current point
          this->extractBasis(fe.u,fe.N,dNdu);
          // Compute Jacobian inverse and derivatives
          dNdu.multiply(dL); // Derivatives w.r.t. xi=[-1,1]
          fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu)*wr[i];
        }

	// Cartesian coordinates of current integration point
        X.assign(fe.Xn * fe.N);
	X.t = time.t;

	// Compute the reduced integration terms of the integrand
	ok = integrand.reducedInt(*A,fe,X);
      }
    }


    // --- Integration loop over all Gauss points in current element -----------

    int jp = iel*ng;
    fe.iGP = firstIp + jp; // Global integration point counter

    for (int i = 0; i < ng && ok; i++, fe.iGP++)
    {
      // Local element coordinate of current integration point
      fe.xi = xg[i];

      // Parameter value of current integration point
      fe.u = param[0] = gpar(1+i,1+iel);

      // Compute basis functions and derivatives
      if (integrand.getIntegrandType() & Integrand::NO_DERIVATIVES)
        this->extractBasis(fe.u,fe.N);
      else if (integrand.getIntegrandType() & Integrand::THIRD_DERIVATIVES)
        this->extractBasis(fe.u,fe.N,dNdu,d2Ndu2,d3Ndu3);
      else if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
        this->extractBasis(fe.u,fe.N,dNdu,d2Ndu2);
      else
        this->extractBasis(fe.u,fe.N,dNdu);

      if (!dNdu.empty())
      {
        // Compute derivatives in terms of physical coordinates
        dNdu.multiply(dL); // Derivatives w.r.t. xi=[-1,1]
        fe.detJxW = utl::Jacobian(Jac,fe.dNdX,fe.Xn,dNdu)*wg[i];
        if (fe.detJxW == 0.0) continue; // skip singular points

        // Compute Hessian of coordinate mapping and 2nd order derivatives
        if (integrand.getIntegrandType() & Integrand::SECOND_DERIVATIVES)
        {
          d2Ndu2.multiply(dL*dL); // 2nd derivatives w.r.t. xi=[-1,1]
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

        if (integrand.getIntegrandType() & Integrand::THIRD_DERIVATIVES)
        {
          d3Ndu3.multiply(dL*dL*dL); // 3rd derivatives w.r.t. xi=[-1,1]
          ok &= utl::Hessian2(fe.d3NdX3,Jac,d3Ndu3);
        }
      }

#if SP_DEBUG > 4
      if (ielm == dbgElm || ielm == -dbgElm || dbgElm == 0)
        std::cout <<"\n"<< fe;
#endif

      // Cartesian coordinates of current integration point
      X.assign(fe.Xn * fe.N);
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

#ifdef SP_DEBUG
    if (ielm == -dbgElm)
      break; // Skipping all elements, except for -dbgElm
#endif
  }

  return ok;
}


bool ASMs1D::integrate (Integrand& integrand, int lIndex,
			GlobalIntegral& glInt,
			const TimeDomain& time)
{
  if (!curv) return true; // silently ignore empty patches

  // Extract the Neumann order flag (1 or higher) for the integrand
  integrand.setNeumannOrder(1 + lIndex/10);

  // Integration of boundary point

  FiniteElement fe(curv->order());
  size_t iel = 0;
  switch (lIndex%10)
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

#ifdef SP_DEBUG
  int ielm = 1+iel;
  if (dbgElm < 0 && ielm != -dbgElm)
    return true; // Skipping all elements, except for -dbgElm
#endif

  fe.iel = MLGE[iel];
  if (fe.iel < 1) return true; // zero-length element

  LocalIntegral* A = integrand.getLocalIntegral(fe.N.size(),fe.iel,true);
  if (!A) return true; // no integrand contributions for this element

  std::map<char,size_t>::const_iterator iit = firstBp.find(lIndex%10);
  fe.iGP = iit == firstBp.end() ? 0 : iit->second;

  // Set up control point coordinates for current element
  bool ok = this->getElementCoordinates(fe.Xn,1+iel);

  if (integrand.getIntegrandType() & Integrand::ELEMENT_CORNERS)
    fe.h = this->getElementEnds(iel+curv->order(),fe.XC);

  if (integrand.getIntegrandType() & Integrand::NODAL_ROTATIONS)
  {
    this->getElementNodalRotations(fe.Tn,iel);
    if (!elmCS.empty()) fe.Te = elmCS[iel];
  }

  // Initialize element matrices
  ok &= integrand.initElementBou(MNPC[iel],*A);

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
    if (lIndex%10 == 1)
      normal.x = -copysign(1.0,Jac(1,1));
    else
      normal.x = copysign(1.0,Jac(1,1));
  }

#if SP_DEBUG > 4
  if (ielm == dbgElm || ielm == -dbgElm || dbgElm == 0)
    std::cout <<"\n"<< fe;
#endif

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


int ASMs1D::findElementContaining (const double* param) const
{
  return curv ? 2 + curv->basis().knotInterval(param[0]) - curv->order() : -1;
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

  RealArray::const_iterator uit  = curv->basis().begin() + curv->order()-1;
  RealArray::const_iterator uend = curv->basis().begin() + curv->numCoefs()+1;
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
  std::iota(index.begin(),index.end(),start-p1+1);
}


bool ASMs1D::getSolution (Matrix& sField, const Vector& locSol,
                          const IntVec& nodes) const
{
  if (!this->ASMbase::getSolution(sField,locSol,nodes))
    return false;
  else if (nf < 6 || !updatedT)
    return true;

  // Extract the total angular rotations as components 4-6
  for (size_t i = 0; i < nodes.size(); i++)
    if (nodes[i] > 0 && (size_t)nodes[i] <= nodalT.size())
    {
      Vec3 rot = nodalT[nodes[i]-1].rotVec();
      sField(4,1+i) = rot.x;
      sField(5,1+i) = rot.y;
      sField(6,1+i) = rot.z;
    }

  return true;
}


bool ASMs1D::evalSolution (Matrix& sField, const Vector& locSol,
                           const int* npe, int) const
{
  // Compute parameter values of the result sampling points
  RealArray gpar;
  if (!this->getGridParameters(gpar,npe[0]-1))
    return false;

  // Evaluate the primary solution at all sampling points
  return this->evalSolution(sField,locSol,&gpar);
}


bool ASMs1D::evalSolution (Matrix& sField, const Vector& locSol,
                           const RealArray* gpar, bool, int deriv, int) const
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


bool ASMs1D::evalProjSolution (Matrix& sField, const Vector& locSol,
                               const int* npe, int) const
{
  // Compute parameter values of the result sampling points
  RealArray gpar;
  if (!this->getGridParameters(gpar,npe[0]-1))
    return false;

  // Evaluate the projected solution at all sampling points
  if (!this->separateProjectionBasis())
    return this->evalSolution(sField,locSol,&gpar,true);

  Fields* f = this->getProjectedFields(locSol);
  if (!f) return false;

  // Evaluate the projected solution field at each point
  Vector vals;
  sField.resize(f->getNoFields(),gpar.size());

  size_t ipt = 0;
  for (double u : gpar)
  {
    f->valueFE(ItgPoint(u),vals);
    sField.fillColumn(++ipt,vals);
  }

  delete f;
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
  Matrix4D d3Ndu3;

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
    else if (integrand.getIntegrandType() & Integrand::THIRD_DERIVATIVES)
      this->extractBasis(fe.u,fe.N,dNdu,d2Ndu2,d3Ndu3);
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
      if (fe.detJxW == 0.0) continue; // skip singular points

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

      if (integrand.getIntegrandType() & Integrand::THIRD_DERIVATIVES)
        utl::Hessian2(fe.d3NdX3,Jac,d3Ndu3);
    }

#if SP_DEBUG > 4
    std::cout <<"\n"<< fe;
#endif

    // Now evaluate the solution field
    if (!integrand.evalSol(solPt,fe,Xtmp*fe.N,ip))
      return false;
    else if (sField.empty())
      sField.resize(solPt.size(),nPoints,true);

    sField.fillColumn(1+i,solPt);
  }

  return true;
}


bool ASMs1D::assembleL2matrices (SparseMatrix& A, StdVector& B,
                                 const L2Integrand& integrand,
                                 bool continuous) const
{
  const size_t nnod = this->getNoProjectionNodes();
  const int p1 = proj->order();

  // Get Gaussian quadrature point coordinates (and weights if continuous)
  const int     ng = continuous ? this->getNoGaussPt(p1,true) : p1 - 1;
  const double* xg = GaussQuadrature::getCoord(ng);
  const double* wg = continuous ? GaussQuadrature::getWeight(ng) : nullptr;
  if (!xg) return false;
  if (continuous && !wg) return false;

  // Compute parameter values of the Gauss points over the whole patch
  // and evaluate the secondary solution at all integration points
  Matrix gp, sField;
  RealArray gpar = this->getGaussPointParameters(gp,ng,xg,proj,true);
  if (!integrand.evaluate(sField,&gpar))
  {
    std::cerr <<" *** ASMs1D::assembleL2matrices: Failed for patch "<< idx+1
              <<" nPoints="<< gpar.size() << std::endl;
    return false;
  }

  double dL = 0.5;
  Vector phi(p1);
  Matrix dNdu, Xnod, J;

  const IntVec& mlge = projMLGE.empty() ? MLGE : projMLGE;
  const IntMat& mnpc = projMNPC.empty() ? MNPC : projMNPC;


  // === Assembly loop over all elements in the patch ==========================

  int ip = 0;
  for (size_t iel = 0; iel < mlge.size(); iel++)
  {
    if (mlge[iel] < 1) continue; // zero-length element

    if (continuous)
    {
      // Set up control point (nodal) coordinates for current element
      if (!this->getElementCoordinates(Xnod,1+iel,mnpc,proj))
        return false;

      int inod1 = mnpc[iel][proj->order()-1];
      dL = 0.5*(*(proj->basis().begin()+inod1+1) -
                *(proj->basis().begin()+inod1));
      if (dL < 0.0)
        return false; // topology error (probably logic error)
    }

    // --- Integration loop over all Gauss points in current element -----------

    for (int i = 0; i < ng; i++, ip++)
    {
      double dJw = 1.0;

      // Fetch basis function values at current integration point
      RealArray basisDerivs;
      proj->computeBasis(gpar[ip],phi,basisDerivs);
      if (continuous)
      {
        dNdu.resize(phi.size(),1);
        dNdu.fillColumn(1,basisDerivs);

        // Compute the Jacobian inverse and derivatives
        dJw = dL*wg[i]*utl::Jacobian(J,dNdu,Xnod,dNdu,false);
        if (dJw == 0.0) continue; // skip singular points
      }

      // Integrate the linear system A*x=B
      for (size_t ii = 0; ii < phi.size(); ii++)
      {
        int inod = mnpc[iel][ii]+1;
        for (size_t jj = 0; jj < phi.size(); jj++)
        {
          int jnod = mnpc[iel][jj]+1;
          A(inod,jnod) += phi[ii]*phi[jj]*dJw;
        }
        for (size_t r = 1; r <= sField.rows(); r++)
          B(inod+(r-1)*nnod) += phi[ii]*sField(r,ip+1)*dJw;
      }
    }
  }

  return true;
}


bool ASMs1D::getNoStructElms (int& n1, int& n2, int& n3) const
{
  n1 = nel;
  n2 = n3 = 0;

  return true;
}


bool ASMs1D::evaluate (const FunctionBase* func, RealArray& values,
                       int, double time) const
{
  Go::SplineCurve* scrv = SplineUtils::project(curv,*func,func->dim(),time);
  if (!scrv)
  {
    std::cerr <<" *** ASMs1D::evaluate: Projection failure."<< std::endl;
    return false;
  }

  values.assign(scrv->coefs_begin(),scrv->coefs_end());
  delete scrv;

  return true;
}


Fields* ASMs1D::getProjectedFields (const Vector& coefs, size_t) const
{
  if (proj == curv || this->getNoProjectionNodes() == 0)
    return nullptr;

  size_t ncmp = coefs.size() / this->getNoProjectionNodes();
  if (ncmp*this->getNoProjectionNodes() == coefs.size())
    return new SplineFields1D(proj,coefs,ncmp);

  std::cerr <<" *** ASMs1D::getProjectedFields: Non-matching coefficent array,"
            <<" size="<< coefs.size() <<" nnod="<< this->getNoProjectionNodes()
            << std::endl;
  return nullptr;
}


void ASMs1D::getElmConnectivities (IntMat& neigh) const
{
  neigh[MLGE[    0]-1] = { -1, MLGE[1]-1 };
  neigh[MLGE[nel-1]-1] = { MLGE[nel-2]-1, -1 };
  for (size_t i = 1; i+1 < nel; i++)
    neigh[MLGE[i]-1] = { MLGE[i-1]-1, MLGE[i+1]-1 };
}


void ASMs1D::getBoundaryElms (int lIndex, int, IntVec& elms) const
{
  if (lIndex == 1)
    elms = {0};
  else
    elms = {MLGE[nel-1]-1};
}
