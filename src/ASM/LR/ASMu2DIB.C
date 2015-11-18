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
#include "Vec3.h"


ASMu2DIB::ASMu2DIB (unsigned char n_s, unsigned char n_f, int max_depth)
  : ASMu2D(n_s,n_f)
{
  maxDepth = max_depth;
  myGeometry = nullptr;
}


ASMu2DIB::ASMu2DIB (const ASMu2DIB& patch, unsigned char n_f)
  : ASMu2D(patch,n_f)
{
  quadPoints = patch.quadPoints;
  myGeometry = nullptr; // because we don't allow multiple patches (yet)
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


ElementBlock* ASMu2DIB::immersedGeometry () const
{
  return myGeometry ? myGeometry->tesselate() : nullptr;
}


void ASMu2DIB::getNoIntPoints (size_t& nPt, size_t&)
{
  nPt = 0;
  for (size_t e = 0; e < quadPoints.size(); e++)
    nPt += quadPoints[e].size();
}


bool ASMu2DIB::generateFEMTopology ()
{
  if (!myGeometry)
  {
    std::cerr <<" *** ASMu2DIB::generateFEMTopology: No geometry description."
              << std::endl;
    return false;
  }

  if (!this->ASMu2D::generateFEMTopology())
    return false;

  size_t i, e, n;
  const int nBasis = lrspline->nBasisFunctions();
  const int nElements = lrspline->nElements();

  // Find the corner point coordinates of each element
  Real3DMat elmCorners;
  for (int iel = 1; iel <= nElements; iel++)
  {
    Real2DMat XC;
    Vec3Vec XCmat;
    this->getElementCorners(iel,XCmat);
    for (i = 0; i < XCmat.size(); i++)
      XC.push_back(RealArray(XCmat[i].ptr(),XCmat[i].ptr()+nsd));

    elmCorners.push_back(XC);
  }

  // Calculate coordinates and weights of the integration points
  bool ok = Immersed::getQuadraturePoints(*myGeometry,elmCorners,
                                          maxDepth,nGauss,quadPoints);

  // Map Gauss point coordinates from bi-unit square to parameter domain (u,v)
  std::vector<LR::Element*>::iterator eit = lrspline->elementBegin();
  for (e = 0; eit != lrspline->elementEnd(); eit++, e++)
  {
    double u0 = (*eit)->umin();
    double v0 = (*eit)->vmin();
    double u1 = (*eit)->umax();
    double v1 = (*eit)->vmax();
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
  return this->ASMu2D::integrate(integrand,glbInt,time,quadPoints);
}


void ASMu2DIB::filterResults (Matrix& field, const ElementBlock* grid) const
{
  if (!myGeometry) return;

  Vec3Vec::const_iterator it = grid->begin_XYZ();
  for (size_t c = 1; c <= field.cols() && it != grid->end_XYZ(); c++, ++it)
    if (myGeometry->Alpha(it->x,it->y) < 1.0)
      for (size_t r = 1; r <= field.rows(); r++)
        field(r,c) = 0.0;
}
