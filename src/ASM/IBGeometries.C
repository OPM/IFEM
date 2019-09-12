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
#include "Function.h"


Oval2D::Oval2D (double r, double x0, double y0, double x1, double y1)
  : Hole2D(r,x0,y0), X1(x1), Y1(y1)
{
  double L = sqrt((X1-Xc)*(X1-Xc) + (Y1-Yc)*(Y1-Yc));
  D2 = L*L + R*R;
  LR = L*R;
}


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


double Oval2D::Alpha (double X, double Y, double) const
{
  // Check if point is located within the first circle
  double r1 = (X-Xc)*(X-Xc) + (Y-Yc)*(Y-Yc);
  if (r1 <= R*R) return 0.0;

  // Check if point is located within the second circle
  double r2 = (X-X1)*(X-X1) + (Y-Y1)*(Y-Y1);
  if (r2 <= R*R) return 0.0;

  // Check if point is outside the extended rectangular domain
  double d = (X1-Xc)*(Y1-Y) - (X1-X)*(Y1-Yc);
  if (d > LR || d < -LR) return 1.0;

  // Finally, check of point is within the two end points
  return r1 <= D2 && r2 <= D2 ? 0.0 : 1.0;
}


double PerforatedPlate2D::Alpha (double X, double Y, double) const
{
  // Determine if point is located within any of the holes or not
  double alpha = 1.0;
  for (size_t i = 0; i < holes.size() && alpha > 0.0; i++)
    alpha = holes[i]->Alpha(X,Y);

  return alpha;
}


ElementBlock* Hole2D::tesselate () const
{
  size_t i, nseg = 360;
  double theta = 0.0, dt = 2.0*M_PI/(double)nseg;

  ElementBlock* grid = new ElementBlock(2);
  grid->unStructResize(nseg,nseg);

  for (i = 0; i < nseg; i++, theta += dt)
    grid->setCoor(i,Xc+R*cos(theta),Yc+R*sin(theta),0.001);

  int n[2] = { 0, 1 };
  int l, ip = 0;
  for (i = 1; i < nseg; i++)
    for (l = 0; l < 2; l++)
      grid->setNode(ip++,n[l]++);

  grid->setNode(ip++,n[0]);
  grid->setNode(ip++,0);

  return grid;
}


ElementBlock* Oval2D::tesselate () const
{
  size_t i, nseg = 360;
  double dt = 2.0*M_PI/(double)nseg;
  double theta = atan2(Y1-Yc,X1-Xc) + 0.5*M_PI;

  ElementBlock* grid = new ElementBlock(2);
  grid->unStructResize(nseg+2,nseg+2);

  for (i = 0; i < nseg/2; i++, theta += dt)
    grid->setCoor(i,Xc+R*cos(theta),Yc+R*sin(theta),0.001);
  grid->setCoor(i,Xc+R*cos(theta),Yc+R*sin(theta),0.001);
  for (++i; i < nseg+2; i++, theta += dt)
    grid->setCoor(i,X1+R*cos(theta),Y1+R*sin(theta),0.001);

  int n[2] = { 0, 1 };
  int l, ip = 0;
  for (i = 1; i < nseg+2; i++)
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

  ElementBlock* grid = holes.front()->tesselate();
  for (size_t i = 1; i < holes.size(); i++)
  {
    ElementBlock* g2 = holes[i]->tesselate();
    std::vector<int> nodes;
    grid->merge(g2,nodes);
    delete g2;
  }

  return grid;
}


PerforatedPlate2D::~PerforatedPlate2D ()
{
  for (Hole2D* h : holes) delete h;
}


void PerforatedPlate2D::addHole (double r, double x, double y)
{
  holes.push_back(new Hole2D(r,x,y));
}


void PerforatedPlate2D::addHole (double r,
                                 double x0, double y0,
                                 double x1, double y1)
{
  holes.push_back(new Oval2D(r,x0,y0,x1,y1));
}


GeoFunc2D::GeoFunc2D (RealFunc* f, double p, double eps)
{
  myAlpha = f;
  myExponent = p;
  threshold = eps;
}


GeoFunc2D::~GeoFunc2D ()
{
  delete myAlpha;
}


double GeoFunc2D::Alpha (double X, double Y, double Z) const
{
  double value = myAlpha ? pow((*myAlpha)(Vec3(X,Y,Z)),myExponent) : 0.0;
  return value < threshold ? 0.0 : 1.0;
}
