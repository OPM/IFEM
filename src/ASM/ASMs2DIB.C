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
#include "ImmersedBoundaries.h"
#include "Vec3.h"


ASMs2DIB::ASMs2DIB (unsigned char n_s, unsigned char n_f, int max_depth)
  : ASMs2D(n_s,n_f)
{
  maxDepth = max_depth;
}


ASMs2DIB::ASMs2DIB (const ASMs2DIB& patch, unsigned char n_f)
  : ASMs2D(patch,n_f)
{
  quadPoints = patch.quadPoints;
}


bool ASMs2DIB::generateFEMTopology ()
{
  if (!this->ASMs2D::generateFEMTopology())
    return false;

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
  if (!Immersed::getQuadraturePoints(elmCorners,maxDepth,nGauss,quadPoints))
    return false;

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

  // Find all nodes within with contributions
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
      this->fix(inod);

  return true;
}


bool ASMs2DIB::integrate (Integrand& integrand,
                          GlobalIntegral& glbInt, const TimeDomain& time)
{
  return this->ASMs2D::integrate(integrand,glbInt,time,quadPoints);
}
