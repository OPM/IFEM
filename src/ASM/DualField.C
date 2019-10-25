// $Id$
//==============================================================================
//!
//! \file DualField.C
//!
//! \date Oct 22 2017
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of dual fields for goal-oriented error estimation.
//!
//==============================================================================

#include "DualField.h"
#include "ASMbase.h"
#include "Vec3Oper.h"


DualRealFunc::DualRealFunc (const Vec3& o, const Vec3& n, const Vec3& XZp,
                            double d, double w, ASMbase* p) : X0(o), normal(n)
{
  tangent.cross(XZp-X0,normal);
  normal.normalize();
  tangent.normalize();
  depth = d;
  width = w;
  patch = p;
}


DualRealFunc::DualRealFunc (const Vec3& o, const Vec3& n,
                            double d, double w, ASMbase* p) : X0(o)
{
  normal.x = n.x;
  normal.y = n.y;
  normal.normalize();
  tangent.x = -normal.y;
  tangent.y =  normal.x;
  depth = d;
  width = w;
  patch = p;
}


DualRealFunc::DualRealFunc (const utl::Point& o, const Vec3Pair& d, ASMbase* p,
                            double eps) : X0(o), Xll(d.first), Xur(d.second)
{
  normal.x = tangent.y = 1.0;
  depth = -eps;
  width = 0.0;
  patch = p;
}


bool DualRealFunc::initPatch (size_t idx)
{
  if (!patch)
    return true;
  else if (idx != patch->idx)
    return false;
  else if (depth > 0.0 || (Xur-Xll).length2() > 0.0 || !X0.u)
    return true;

  // Find the element(s) containing or close to the point X0
  Delem.clear();
  double uv[3] = { X0.u[0], X0.u[1], X0.u[2] };
  Delem.insert(patch->findElementContaining(uv));

  // If depth is negative, it is interpreted as a small epsilon in the
  // parameter domain and used to make sure we find all elements in case
  // the point is on the interface between two more more elements
  if (patch->getNoParamDim() == 1 && depth < 0.0)
  {
    uv[0] = X0.u[0] - depth;
    Delem.insert(patch->findElementContaining(uv));
    uv[0] = X0.u[0] + depth;
    Delem.insert(patch->findElementContaining(uv));
  }
  else if (patch->getNoParamDim() == 2 && depth < 0.0)
    for (size_t i = 0; i < 8; i++)
    {
      double theta = 0.25*M_PI*i;
      uv[0] = X0.u[0] - depth*cos(theta);
      uv[1] = X0.u[1] - depth*sin(theta);
      Delem.insert(patch->findElementContaining(uv));
    }

#ifdef SP_DEBUG
  std::cout <<"DualRealFunc::initPatch: Element(s) in domain:";
  for (int e : Delem) if (e > 0) std::cout <<" "<< e;
  std::cout << std::endl;
#endif
  return true;
}


bool DualRealFunc::inDomain (const Vec3& X) const
{
  if (depth > 0.0)
  {
    double x = (X-X0)*normal;
    if (x > 0.0 || -x >= depth)
      return false;
    else if (width <= 0.0)
      return true;

    double y = (X-X0)*tangent;
    return y <= 0.5*width && y >= -0.5*width;
  }

  if (patch && !Delem.empty())
  {
    const Vec4* Xp = dynamic_cast<const Vec4*>(&X);
#ifdef SP_DEBUG
    if (!Xp || !Xp->u)
    {
      std::cerr <<" *** DualRealFunc::inDomain: The point "<< X
                <<" lacks parameters."<< std::endl;
      return false;
    }
#endif
    int iel = patch->findElementContaining(Xp->u);
    if (iel < 1 || Delem.find(iel) == Delem.end()) return false;
#ifdef SP_DEBUG
    std::cout <<"DualRealFunc::inDomain("<< X <<") --> iel="<< iel << std::endl;
#endif
    return true;
  }

  // Rectangular point domain
  return (X.x >= Xll.x && X.x <= Xur.x) &&
         (X.y >= Xll.y && X.y <= Xur.y) &&
         (X.z >= Xll.z && X.z <= Xur.z);
}


double DualRealFunc::value (const Vec3& X, bool ignoreDomain) const
{
  if (ignoreDomain || this->inDomain(X))
    return depth > 0.0 ? 1.0 + (X-X0)*normal/depth : 1.0;

  return 0.0; // outside the function domain
}


DualVecFunc::DualVecFunc (int c, const Vec3& o, const Vec3& n, const Vec3& XZp,
                          double d, double w, ASMbase* p)
  : VecFunc(3), W(o,n,XZp,d,w,p)
{
  comp = c;
}


DualVecFunc::DualVecFunc (int c, const Vec3& o, const Vec3& n,
                          double d, double w, ASMbase* p)
  : VecFunc(2), W(o,n,d,w,p)
{
  comp = c > 2 ? 6 : c;
}


DualVecFunc::DualVecFunc (int c, const utl::Point& o, const Vec3Pair& d,
                          ASMbase* p, double eps)
  : VecFunc(2), W(o,d,p,eps)
{
  comp = c;
}


/*!
  The field does not evaluate to zero outside the defined function domain, since
  it is used to find the control point values of a spline representation of it.
*/

Vec3 DualVecFunc::evaluate (const Vec3& X) const
{
  double w = W.value(X,true);
  if (fabs(w) < 1.0e-12)
    return Vec3();

  if (W.isPointExtraction())
    switch (comp) {
    case 1: return Vec3(W.ecc(X,1),0.0,0.0);
    case 2: return Vec3(0.0,W.ecc(X,2),0.0);
    case 3: return Vec3(0.0,0.0,W.ecc(X,3));
    case 4: return Vec3(0.0,0.5*W.ecc(X,3),0.5*W.ecc(X,2));
    case 5: return Vec3(0.5*W.ecc(X,3),0.0,0.5*W.ecc(X,1));
    case 6: return Vec3(0.5*W.ecc(X,2),0.5*W.ecc(X,1),0.0);
    }
  else
    switch (comp) {
    case 1: return  w*W.x();
    case 2: return  w*W.y();
    case 3: return  w*W.z();
    case 4: return -w*W.ecc(X,3)*W.y() + w*W.ecc(X,2)*W.z();
    case 5: return  w*W.ecc(X,3)*W.x() - w*W.ecc(X,1)*W.z();
    case 6: return -w*W.ecc(X,2)*W.x() + w*W.ecc(X,1)*W.y();
    }

  return Vec3();
}
