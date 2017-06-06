// $Id$
//==============================================================================
//!
//! \file ASMs2DIB.C
//!
//! \date Dec 18 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Assembly of structured 2D spline FE models with immersed boundaries.
//!
//==============================================================================

#include "GoTools/geometry/SplineSurface.h"

#include "ASMs2DIB.h"
#include "IBGeometries.h"
#include "ElementBlock.h"
#include "Vec3.h"


ASMs2DIB::ASMs2DIB (unsigned char n_s, unsigned char n_f, int max_depth)
  : ASMs2D(n_s,n_f)
{
  maxDepth = max_depth;
  myGeometry = nullptr;
  myLines = nullptr;
}


ASMs2DIB::ASMs2DIB (const ASMs2DIB& patch, unsigned char n_f)
  : ASMs2D(patch,n_f), quadPoints(patch.quadPoints)
{
  myGeometry = nullptr;
  myLines = nullptr;
}


ASMs2DIB::~ASMs2DIB ()
{
  delete myGeometry;
  delete myLines;
}


void ASMs2DIB::addHole (double R, double Xc, double Yc)
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


void ASMs2DIB::addHole (double R, double X1, double Y1, double X2, double Y2)
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


ElementBlock* ASMs2DIB::immersedGeometry () const
{
  if (!myGeometry) return nullptr;

  ElementBlock* geo = myGeometry->tesselate();
  if (myLines)
  {
    std::vector<int> nodes;
    geo->merge(myLines,nodes);
  }

  return geo;
}


void ASMs2DIB::getNoIntPoints (size_t& nPt, size_t& nIPt)
{
  firstIp = nPt;
  this->ASMbase::getNoIntPoints(nPt,nIPt);

  if (myGeometry)
  {
    nPt = firstIp;
    for (size_t e = 0; e < quadPoints.size(); e++)
      nPt += quadPoints[e].size();
  }

#ifdef SP_DEBUG
  std::cout <<"Number of quadrature points in patch "<< idx+1
            <<": "<< nPt-firstIp << std::endl;
#endif
}


bool ASMs2DIB::isIntersected (int iel, bool checkIfInDomainOnly) const
{
  if (iel < 0 || (size_t)iel >= quadPoints.size())
    return false;
  else if (quadPoints[iel].empty())
    return false;
  else if (checkIfInDomainOnly)
    return true;

  return quadPoints[iel].size() > (size_t)nGauss*nGauss;
}


bool ASMs2DIB::generateFEMTopology ()
{
  if (!this->ASMs2D::generateFEMTopology())
    return false;
  else if (!myGeometry)
    return true;

  size_t i, e, n;
  int i1, i2, inod;
  const int n1 = surf->numCoefs_u();
  const int n2 = surf->numCoefs_v();
  const int p1 = surf->order_u();
  const int p2 = surf->order_v();

  // Find the corner point coordinates of each element
  Real3DMat elmCorners;
  for (i2 = p2-1; i2 < n2; i2++)
    for (i1 = p1-1; i1 < n1; i1++)
    {
      Real2DMat XC;
      if (surf->knotSpan(0,i1) > 0.0)
        if (surf->knotSpan(1,i2) > 0.0)
        {
          std::vector<Vec3> XCmat;
          this->getElementCorners(i1,i2,XCmat);
          for (i = 0; i < XCmat.size(); i++)
            XC.push_back(RealArray(XCmat[i].ptr(),XCmat[i].ptr()+nsd));
        }
      elmCorners.push_back(XC);
    }

  // Calculate coordinates and weights of the integration points
  if (Immersed::plotCells)
    myLines = new ElementBlock(2);
  bool ok = Immersed::getQuadraturePoints(*myGeometry,elmCorners,
                                          maxDepth,nGauss,quadPoints,myLines);

  // Map Gauss point coordinates from bi-unit square to parameter domain (u,v)
  RealArray::const_iterator vit = surf->basis(1).begin() + p2-1;
  double vcurr, vprev = *(vit++);
  for (e = 0, i2 = p2-1; i2 < n2; i2++, ++vit)
  {
    vcurr = *vit;
    RealArray::const_iterator uit = surf->basis(0).begin() + p1-1;
    double ucurr, uprev = *(uit++);
    for (i1 = p1-1; i1 < n1; i1++, ++uit, e++)
    {
#if SP_DEBUG > 1
      std::cout <<"\n Element "<< MLGE[e] <<":\n";
#endif
      ucurr = *uit;
      for (i = 0; i < quadPoints[e].size(); i++)
      {
        double& xi  = quadPoints[e][i][0];
        double& eta = quadPoints[e][i][1];
#if SP_DEBUG > 1
        std::cout <<"\tItg.point "<< i+1 <<": xi,eta = "<< xi <<" "<< eta;
#endif
        xi  = 0.5*((ucurr-uprev)*xi  + ucurr+uprev);
        eta = 0.5*((vcurr-vprev)*eta + vcurr+vprev);
#if SP_DEBUG > 1
        std::cout <<" --> u,v = "<< xi <<" "<< eta << std::endl;
#endif
      }
      uprev = ucurr;
    }
    vprev = vcurr;
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
  for (inod = 1; inod <= n1*n2; inod++)
    if (activeNodes.find(inod) == activeNodes.end())
      this->fix(inod,12);

  if (Immersed::stabilization == Immersed::ALL_INTERFACES)
    ok &= this->addInterfaceElms(Intersected(*this,true));
  else if (Immersed::stabilization == Immersed::SUBDIV_INTERFACES)
    ok &= this->addInterfaceElms(Intersected(*this,false));

  return ok;
}


short int ASMs2DIB::Intersected::hasContribution (int I, int J) const
{
  const int n1 = myPatch.surf->numCoefs_u();
  const int n2 = myPatch.surf->numCoefs_v();
  const int p1 = myPatch.surf->order_u();
  const int p2 = myPatch.surf->order_v();

  int iel = I-p1 + (n1-p1+1)*(J-p2); // Zero-based index of this element
  if (myPatch.quadPoints[iel].empty())
    return 0; // This element is completely outside the domain

  int jel[4];
  jel[0] = I > p1 ? iel - 1         : 0; // West neighbor
  jel[1] = I < n1 ? iel + 1         : 0; // East neighbor
  jel[2] = J > p2 ? iel - (n1-p1+1) : 0; // South neighbor
  jel[3] = J < n2 ? iel + (n1-p1+1) : 0; // North neighbor

  // Check if the element itself is intersected
  bool isCut = myAll || myPatch.isIntersected(iel);

  // Check if any of the neighboring elements are intersected, or if
  // this element is intersected, only check if the neighbor is in the domain
  short int status = 0, s = 1;
  for (short int i = 0; i < 4; i++, s *= 2)
    if (alsoSW || i%2 == 1)
      if (jel[i] && myPatch.isIntersected(jel[i],isCut))
        status += s;

  return status;
}


bool ASMs2DIB::integrate (Integrand& integrand,
                          GlobalIntegral& glbInt, const TimeDomain& time)
{
  if (!myGeometry)
    return this->ASMs2D::integrate(integrand,glbInt,time);

  if (!this->integrate(integrand,glbInt,time,quadPoints))
    return false;

  switch (Immersed::stabilization) {
  case Immersed::ALL_INTERFACES:
    return this->integrate(integrand,glbInt,time,Intersected(*this,true));
  case Immersed::SUBDIV_INTERFACES:
    return this->integrate(integrand,glbInt,time,Intersected(*this,false));
  default:
    return true;
  }
}


void ASMs2DIB::filterResults (Matrix& field, const ElementBlock* grid) const
{
  if (!myGeometry) return;

  Vec3Vec::const_iterator it = grid->begin_XYZ();
  for (size_t c = 1; c <= field.cols() && it != grid->end_XYZ(); c++, ++it)
    if (myGeometry->Alpha(it->x,it->y) < 1.0)
      for (size_t r = 1; r <= field.rows(); r++)
        field(r,c) = 0.0;
}
