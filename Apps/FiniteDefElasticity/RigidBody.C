// $Id$
//==============================================================================
//!
//! \file RigidBody.C
//!
//! \date Jun 05 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of rigid bodies in contact analysis.
//!
//==============================================================================

#include "RigidBody.h"
#include "ElementBlock.h"
#include "Utilities.h"
#include "Tensor.h"
#include "Vec3Oper.h"

#ifndef epsZ
//! \brief Zero tolerance for point coalescence.
#define epsZ 1.0e-16
#endif


RigidBody::RigidBody (unsigned char n, size_t np) :  nsd(n)
{
  X0.resize(np);
  Xn.resize(np);
  MLGN.resize(np,0);

  eps = 1.0;
  code = -1;
}


void RigidBody::initPoints (const std::vector<Vec3>& p)
{
  for (size_t i = 0; i < X0.size() && i < p.size(); i++)
    X0[i] = p[i];
}


void RigidBody::initNodes (const std::vector<int>& nodes)
{
  for (size_t i = 0; i < MLGN.size() && i < nodes.size(); i++)
    MLGN[i] = nodes[i];
}


void RigidBody::renumberNodes (const std::map<int,int>& old2new)
{
  for (size_t i = 0; i < MLGN.size(); i++)
    utl::renumber(MLGN[i],old2new);
}


void RigidBody::print (std::ostream& os) const
{
  for (size_t i = 0; i < MLGN.size(); i++)
    os <<"\n\tP"<< i+1 <<" (node "<< MLGN[i] <<"): "<< X0[i];
  os <<"\n\tProperty code: "<< code
     <<"\n\tPenalty parameter: "<< eps << std::endl;
}


bool RigidBody::update (const RealArray& displ)
{
  for (size_t i = 0; i < MLGN.size(); i++)
    if (displ.empty())
      Xn[i] = X0[i];
    else
    {
      int n = MLGN[i];
#ifdef INDEX_CHECK
      if (n < 1 || nsd*(size_t)n > displ.size())
      {
	std::cerr <<" *** RigidBody::update: Global DOF "<< nsd*n
		  <<" is out of range [1,"<< displ.size() <<"]"<< std::endl;
	return false;
      }
#endif
      for (unsigned char j = 0; j < nsd; j++)
	Xn[i][j] = X0[i][j] + displ[nsd*(n-1)+j];
    }

  return true;
}


Vec3 RigidBody::getPosition () const
{
  return Xn.front() - X0.front();
}


void RigidSphere::print (std::ostream& os) const
{
  const_cast<RigidSphere*>(this)->update(RealArray());

  os <<"\tRigid Sphere: R = "<< R;
  this->RigidBody::print(os);
}


ElementBlock* RigidSphere::tesselate () const
{
  const size_t ntheta = 180; // Number of elements around equator
  const size_t nphi   =  60; // Number of elements from pole to pole

  ElementBlock* g = new ElementBlock(4);
  g->unStructResize(nphi*ntheta,2+(nphi-1)*ntheta);

  size_t m, n, ip = 2;
  Vec3 Xo(X0.front());

  g->setCoor(0,Xo.x,Xo.y,Xo.z+R);
  g->setCoor(1,Xo.x,Xo.y,Xo.z-R);
  for (n = 0; n < ntheta; n++)
  {
    double theta = M_PI*n*2.0/(double)ntheta;
    double Rct = R*cos(theta);
    double Rst = R*sin(theta);
    for (m = 1; m < nphi; m++, ip++)
    {
      double phi = M_PI*m/(double)nphi;
      g->setCoor(ip,Xo.x+Rct*sin(phi),Xo.y+Rst*sin(phi),Xo.z+R*cos(phi));
    }
  }

  for (n = ip = 0; n+1 < ntheta; n++)
  {
    g->setNode(ip++,0);
    g->setNode(ip++,nphi* n   -n+2);
    g->setNode(ip++,nphi*(n+1)-n+1);
    g->setNode(ip++,0);
    for (m = 1; m+1 < nphi; m++)
    {
      g->setNode(ip++,nphi* n   -n+m+1);
      g->setNode(ip++,nphi* n   -n+m+2);
      g->setNode(ip++,nphi*(n+1)-n+m+1);
      g->setNode(ip++,nphi*(n+1)-n+m);
    }
    g->setNode(ip++,nphi*(n+1)-n);
    g->setNode(ip++,1);
    g->setNode(ip++,1);
    g->setNode(ip++,nphi*(n+2)-n-1);
  }

  g->setNode(ip++,0);
  g->setNode(ip++,(nphi-1)*(ntheta-1)+2);
  g->setNode(ip++,2);
  g->setNode(ip++,0);
  for (m = 1; m+1 < nphi; m++)
  {
    g->setNode(ip++,(nphi-1)*(ntheta-1)+m+1);
    g->setNode(ip++,(nphi-1)*(ntheta-1)+m+2);
    g->setNode(ip++,m+2);
    g->setNode(ip++,m+1);
  }
  g->setNode(ip++,(nphi-1)*ntheta+1);
  g->setNode(ip++,1);
  g->setNode(ip++,1);
  g->setNode(ip++,nphi);

  return g;
}


double RigidSphere::evalGap (const Vec3& X, Vec3& normal, RealArray& N) const
{
  normal = X - Xn[0];

  N.resize(1);
  N[0] = 1.0;

  return normal.normalize() - R;
}


void RigidSphere::geometricStiffness (const Vec3& X, Vec3& normal,
				      Matrix& Tg, RealArray& N) const
{
  Tg.resize(nsd,nsd,true);

  Vec3 dX = X - Xn[0];
  double d = dX.normalize(epsZ);

  if (d > epsZ)
    normal = dX;
  else
  {
    normal *= -1.0;
    d = 1.0;
  }

  for (unsigned char j = 0; j < nsd; j++)
    for (unsigned char i = 0; i < nsd; i++)
      Tg(i+1,j+1) = ((i == j ? 1.0 : 0.0) - dX[i]*dX[j]) / d;

  N.resize(1);
  N[0] = 1.0;
}


RigidCylinder::RigidCylinder (double r, unsigned char n)
  : RigidBody(n,n-1), R(r), L(1.0)
{
  if (nsd == 2) dirC.z = 1.0;
}


void RigidCylinder::print (std::ostream& os) const
{
  const_cast<RigidCylinder*>(this)->update(RealArray());

  os <<"\tRigid Cylinder: R = "<< R <<", L = "<< L <<", dirC = "<< dirC;
  this->RigidBody::print(os);
}


ElementBlock* RigidCylinder::tesselate () const
{
  const size_t ntheta = 180; // Number of element in circular direction

  ElementBlock* g = new ElementBlock(4);
  g->unStructResize(3*ntheta,2+2*ntheta);

  size_t n, ip;
  Tensor Tlg(dirC);
  Vec3 X1 = X0.size() > 1 ? X0[1] : X0[0] + L*dirC;

  g->setCoor(0,X0[0]);
  g->setCoor(1+ntheta,X1);
  for (n = 0; n < ntheta; n++)
  {
    double theta = M_PI*n*2.0/(double)ntheta;
    Vec3 X(R*cos(theta),R*sin(theta),0.0);
    g->setCoor(1+n,X0[0]+Tlg*X);
    g->setCoor(2+ntheta+n,X1+Tlg*X);
  }

  for (n = ip = 0; n+1 < ntheta; n++)
  {
    g->setNode(ip++,0);
    g->setNode(ip++,2+n);
    g->setNode(ip++,1+n);
    g->setNode(ip++,0);

    g->setNode(ip++,2+n);
    g->setNode(ip++,1+n);
    g->setNode(ip++,2+ntheta+n);
    g->setNode(ip++,3+ntheta+n);

    g->setNode(ip++,1+ntheta);
    g->setNode(ip++,2+ntheta+n);
    g->setNode(ip++,3+ntheta+n);
    g->setNode(ip++,1+ntheta);
  }

  g->setNode(ip++,0);
  g->setNode(ip++,1);
  g->setNode(ip++,ntheta);
  g->setNode(ip++,0);

  g->setNode(ip++,1);
  g->setNode(ip++,ntheta);
  g->setNode(ip++,1+ntheta+ntheta);
  g->setNode(ip++,2+ntheta);

  g->setNode(ip++,1+ntheta);
  g->setNode(ip++,1+ntheta+ntheta);
  g->setNode(ip++,2+ntheta);
  g->setNode(ip++,1+ntheta);

  return g;
}


bool RigidCylinder::update (const RealArray& displ)
{
  if (!this->RigidBody::update(displ))
    return false;

  if (nsd == 3)
  {
    // Update the cylinder axis direction vector
    dirC = Xn[1] - Xn[0];
    if ((L = dirC.normalize(epsZ)) <= epsZ)
    {
      std::cerr <<" *** RigidCylinder::update: Degenerated cylinder.\n";
      return false;
    }
  }

  return true;
}


double RigidCylinder::evalGap (const Vec3& X, Vec3& normal, RealArray& N) const
{
  normal = X - Xn.front();
  double z = normal*dirC;
  normal -= z*dirC;

  N.resize(nsd-1);
  if (nsd == 2)
    N[0] = 1.0;
  else
  {
    N[0] = 1.0 - z/L;
    N[1] = z/L;
  }

  return normal.normalize() - R;
}


void RigidCylinder::geometricStiffness (const Vec3& X, Vec3& normal,
					Matrix& Tg, RealArray& N) const
{
  Tg.resize(nsd,nsd,true);

  Vec3 dX = X - Xn.front();
  double z = dX*dirC;
  dX -= z*dirC;
  double d = dX.normalize(epsZ);

  if (d > epsZ)
    normal = dX;
  else
  {
    normal *= -1.0;
    d = 1.0;
  }

  for (unsigned char j = 0; j < nsd; j++)
    for (unsigned char i = 0; i < nsd; i++)
      Tg(i+1,j+1) = ((i == j ? 1.0 : 0.0) - dirC[i]*dirC[j] - dX[i]*dX[j]) / d;

  N.resize(nsd-1);
  if (nsd == 2)
    N[0] = 1.0;
  else
  {
    N[0] = 1.0 - z/L;
    N[1] = z/L;
  }
}


void RigidPlane::print (std::ostream& os) const
{
  const_cast<RigidPlane*>(this)->update(RealArray());

  os <<"\tRigid Plane: Triangle area = "<< Area <<", normal = "<< nVec;
  this->RigidBody::print(os);
}


ElementBlock* RigidPlane::tesselate () const
{
  // One single rectangular element is sufficient
  ElementBlock* g = new ElementBlock(4);
  g->resize(2,2);

  g->setCoor(0,X0[0]);
  g->setCoor(1,X0[1]);
  if (X0.size() > 2)
  {
    g->setCoor(3,X0[2]);
    g->setCoor(2,X0[1]+X0[2]-X0[0]);
  }
  else
  {
    Vec3 X2(X0.front());
    X2.z += X0[1].x - X0[0].x;
    g->setCoor(3,X2);
    g->setCoor(2,X0[1]+X2-X0[0]);
  }
  for (int i = 0; i < 4; i++)
    g->setNode(i,i);

  return g;
}


bool RigidPlane::update (const RealArray& displ)
{
  if (!this->RigidBody::update(displ))
    return false;

  // Calculate the plane normal
  Vec3 t1 = Xn[1] - Xn.front();
  if (nsd == 2)
  {
    if ((Area = t1.normalize(epsZ)) <= epsZ)
    {
      std::cerr <<" *** RigidPlane::update: Degenerated plane.\n";
      return false;
    }
    nVec.x = -t1.y;
    nVec.y =  t1.x;
  }
  else
  {
    Vec3 t2 = Xn[2] - Xn.front();
    if ((Area = nVec.cross(t1,t2).normalize(epsZ)) <= epsZ)
    {
      std::cerr <<" *** RigidPlane::update: Degenerated plane.\n";
      return false;
    }
  }

  return true;
}


double RigidPlane::evalGap (const Vec3& X, Vec3& normal, RealArray& N) const
{
  normal = nVec;

  Vec3 dX = X - Xn.front();
  double gN = dX*nVec;
  dX -= gN*nVec;

  N.resize(nsd);
  if (nsd == 2)
  {
    double xi = dX.length()/Area;
    N[0] = 1.0 - xi;
    N[1] = xi;
  }
  else
  {
    // Area coordinates
    double A2 = Vec3(dX, Xn[1] - Xn.front()).length()/Area;
    double A3 = Vec3(Xn[2] - Xn.front(), dX).length()/Area;

    N[0] = 1.0 - A2 - A3;
    N[1] = A2;
    N[2] = A3;
  }

  return gN;
}


void RigidPlane::geometricStiffness (const Vec3& X, Vec3& normal,
				     Matrix& Tg, RealArray& N) const
{
  // No geometric stiffness, but invoke evalGap to evaluate N
  Tg.resize(0,0,true);
  this->evalGap(X,normal,N);
}
