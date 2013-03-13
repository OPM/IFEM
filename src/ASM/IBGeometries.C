
// $Id$
//==============================================================================
//!
//! \file IBGeometries.C
//!
//! \date Jan 24 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Physical geometries for immersed boundary simulations.
//!
//==============================================================================

#include "IBGeometries.h"
#include "ElementBlock.h"


double Hole2D::Alpha (double X, double Y, double) const
{
  // Determine radius of the point to be checked
  double r2 = (X-Xc)*(X-Xc) + (Y-Yc)*(Y-Yc);

  // Determine if point is located within the hole or not
  double alpha = r2 > R*R ? 1.0 : 0.0;

  // Alpha is the penalization parameter in the sense of equation (29) in the
  // immersed boundary paper. It can be used either to indicate where a point
  // is located, or can be used as a penalization value for the local stiffness
  // contribution at the current integration point.
  // In the latter case, it needs to have data type double.
  return alpha;
}


double PerforatedPlate2D::Alpha (double X, double Y, double) const
{
  // Determine if point is located within any of the holes or not
  double alpha = 1.0;
  for (size_t i = 0; i < holes.size() && alpha > 0.0; i++)
    alpha = holes[i].Alpha(X,Y);

  return alpha;
}


ElementBlock* Hole2D::tesselate () const
{
  size_t i, nseg = 360;
  double theta = 0.0, dt = 2.0*M_PI/(double)nseg;

  ElementBlock* grid = new ElementBlock(2);
  grid->unStructResize(nseg,nseg);

  for (i = 0; i < nseg; i++, theta += dt)
    grid->setCoor(i,Xc+R*cos(theta),Yc+R*sin(theta),0.01);

  int n[2] = { 0, 1 };
  int l, ip = 0;
  for (i = 1; i < nseg; i++)
    for (l = 0; l < 2; l++)
      grid->setNode(ip++,n[l]++);

  grid->setNode(ip++,n[0]);
  grid->setNode(ip++,0);

  return grid;
}


ElementBlock* PerforatedPlate2D::tesselate () const
{
  if (holes.empty())
    return new ElementBlock(2);

  ElementBlock* grid = holes.front().tesselate();
  for (size_t i = 1; i < holes.size(); i++)
  {
    ElementBlock* g2 = holes[i].tesselate();
    std::vector<int> nodes;
    grid->merge(g2,nodes);
    delete g2;
  }

  return grid;
}
