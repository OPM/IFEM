// $Id$
//==============================================================================
//!
//! \file ASMu2DIB.C
//!
//! \date Dec 18 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of unstructured 2D spline FE models with immersed boundary.
//!
//==============================================================================

#include "LRSpline/LRSplineSurface.h"
#include "LRSpline/Element.h"

#include "ASMu2DIB.h"
#include "IBGeometries.h"
#include "ElementBlock.h"
#include "Point.h"


ASMu2DIB::ASMu2DIB (unsigned char n_s, unsigned char n_f, int max_depth)
  : ASMu2D(n_s,n_f)
{
  maxDepth = max_depth;
  myGeometry = nullptr;
}


ASMu2DIB::ASMu2DIB (const ASMu2DIB& patch, unsigned char n_f)
  : ASMu2D(patch,n_f), quadPoints(patch.quadPoints)
{
  maxDepth = patch.maxDepth;
  myGeometry = nullptr;
}


ASMu2DIB::~ASMu2DIB ()
{
  delete myGeometry;
}


void ASMu2DIB::addHole (double R, double Xc, double Yc)
{
  std::cout <<"\tHole Xc={"<< Xc <<","<< Yc <<"} R="<< R << std::endl;

  if (!myGeometry)
    myGeometry = new Hole2D(R,Xc,Yc);
  else
  {
    PerforatedPlate2D* plate = dynamic_cast<PerforatedPlate2D*>(myGeometry);
    if (!plate)
    {
      Hole2D* hole = dynamic_cast<Hole2D*>(myGeometry);
      if (hole)
        myGeometry = plate = new PerforatedPlate2D(hole);
      else
        return;
    }
    plate->addHole(R,Xc,Yc);
  }
}


void ASMu2DIB::addHole (double R, double X1, double Y1, double X2, double Y2)
{
  std::cout <<"\tOval X1={"<< X1 <<","<< Y1
            <<"} X2={"<< X2 <<","<< Y2 <<"} R="<< R << std::endl;

  if (!myGeometry)
    myGeometry = new Oval2D(R,X1,Y1,X2,Y2);
  else
  {
    PerforatedPlate2D* plate = dynamic_cast<PerforatedPlate2D*>(myGeometry);
    if (!plate)
    {
      Hole2D* hole = dynamic_cast<Hole2D*>(myGeometry);
      if (hole)
        myGeometry = plate = new PerforatedPlate2D(hole);
      else
        return;
    }
    plate->addHole(R,X1,Y1,X2,Y2);
  }
}


bool ASMu2DIB::setGeometry (RealFunc* f, double power, double threshold)
{
  if (myGeometry) return false;

  myGeometry = new GeoFunc2D(f,power,threshold);

  return true;
}


ElementBlock* ASMu2DIB::immersedGeometry () const
{
  return myGeometry ? myGeometry->tesselate() : nullptr;
}


void ASMu2DIB::getNoIntPoints (size_t& nPt, size_t& nIPt)
{
  firstIp = nPt;
  this->ASMbase::getNoIntPoints(nPt,nIPt);

  if (myGeometry)
  {
    nPt = firstIp;
    for (const Real2DMat& qp : quadPoints)
      nPt += qp.size();
  }
}


bool ASMu2DIB::generateFEMTopology ()
{
  if (!this->ASMu2D::generateFEMTopology())
    return false;
  else if (!myGeometry)
    return true;

  size_t i, e, n;
  const int nBasis = lrspline->nBasisFunctions();
  const int nElements = lrspline->nElements();

  // Find the corner points of each element
  std::vector<PointVec> elmCorners(nElements);
  for (int iel = 1; iel <= nElements; iel++)
    this->getCornerPoints(iel,elmCorners[iel-1]);

  // Calculate coordinates and weights of the integration points
  bool ok = Immersed::getQuadraturePoints(*myGeometry,elmCorners,
                                          maxDepth,nGauss,quadPoints);

  // Map Gauss point coordinates from bi-unit square to parameter domain (u,v)
  e = 0;
  for (const LR::Element* el : lrspline->getAllElements())
  {
    double u0 = el->umin();
    double v0 = el->vmin();
    double u1 = el->umax();
    double v1 = el->vmax();
#if SP_DEBUG > 1
    std::cout <<"\n Element "<< MLGE[e] <<":\n";
#endif
    for (i = 0; i < quadPoints[e].size(); i++)
    {
      double& xi  = quadPoints[e][i][0];
      double& eta = quadPoints[e][i][1];
#if SP_DEBUG > 1
      std::cout <<"\tItg.point "<< i+1 <<": xi,eta = "<< xi <<" "<< eta;
#endif
      xi  = 0.5*((u1-u0)*xi  + u1+u0);
      eta = 0.5*((v1-v0)*eta + v1+v0);
#if SP_DEBUG > 1
      std::cout <<" --> u,v = "<< xi <<" "<< eta << std::endl;
#endif
    }
    ++e;
  }

  // Find all nodes with contributions
  std::set<int> activeNodes;
  for (e = 0; e < quadPoints.size(); e++)
    if (quadPoints[e].empty())
      std::cout <<"\n Element "<< MLGE[e] <<" is completely outside the domain";
    else for (n = 0; n < MNPC[e].size(); n++)
      activeNodes.insert(MNPC[e][n]+1);
  std::cout <<"\n"<< std::endl;

  // Automatically fix all the non-active nodes
  for (int inod = 1; inod <= nBasis; inod++)
    if (activeNodes.find(inod) == activeNodes.end())
      this->fix(inod);

  return ok;
}


bool ASMu2DIB::integrate (Integrand& integrand,
                          GlobalIntegral& glbInt, const TimeDomain& time)
{
  if (!myGeometry)
    return this->ASMu2D::integrate(integrand,glbInt,time);

  return this->integrate(integrand,glbInt,time,quadPoints);
}


void ASMu2DIB::filterResults (Matrix& field, const ElementBlock* grid) const
{
  if (!myGeometry) return;

  for (size_t n = 0; n < field.cols() && n < grid->getNoNodes(); n++)
    if (myGeometry->Alpha(Vec4(grid->getCoord(n),0.0,grid->getParam(n))) < 1.0)
      for (size_t r = 1; r <= field.rows(); r++)
        field(r,n+1) = 0.0;
}
