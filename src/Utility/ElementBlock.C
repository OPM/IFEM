// $Id$
//==============================================================================
//!
//! \file ElementBlock.C
//!
//! \date Jun 03 2009
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Storage of a standard FE grid block of uniform element type.
//!
//==============================================================================

#include "ElementBlock.h"
#include "Vec3Oper.h"
#include "Tensor.h"
#include <algorithm>
#include <numeric>


double ElementBlock::eps = 0.0;

//! \brief List of supported number of element nodes.
static const std::array<size_t,7> legalNENs = { 2, 3, 4, 6, 8, 9, 27 };


ElementBlock::ElementBlock (size_t nenod)
{
  if (std::find(legalNENs.begin(),legalNENs.end(),nenod) == legalNENs.end())
  {
    std::cout <<"ElementBlock: Invalid number of element nodes "<< nenod
              <<" reset to 8"<< std::endl;
    nenod = 8;
  }
  nen = nenod;
}


void ElementBlock::resize (size_t nI, size_t nJ, size_t nK)
{
  coord.resize(nI*nJ*nK);
  param.resize(nI*nJ*nK);
  if (nen == 2 && nJ < 2 && nK < 2)
    MMNPC.resize(2*(nI-1));
  else if (nen == 3 && nK < 2)
    MMNPC.resize(6*(nI-1)*(nJ-1));
  else if (nen == 4 && nK < 2)
    MMNPC.resize(4*(nI-1)*(nJ-1));
  else if (nen == 6 && nK < 2)
    MMNPC.resize(12*(nI-1)*(nJ-1));
  else if (nen == 9 && nK < 2)
    MMNPC.resize(9*(nI-1)*(nJ-1)/4);
  else if (nen == 8)
    MMNPC.resize(8*(nI-1)*(nJ-1)*(nK-1));
  else if (nen == 27)
    MMNPC.resize(27*(nI-1)*(nJ-1)*(nK-1)/8);

  MINEX.resize(MMNPC.size()/nen,0);
  std::iota(MINEX.begin(),MINEX.end(),1);
}


void ElementBlock::unStructResize (size_t nEL, size_t nPts, size_t nMNPC)
{
  if (nMNPC > 0)
  {
    nen = nMNPC/nEL;
    if (nen*nEL != nMNPC)
      nen = 0; // Unstructured grid with varying element types
  }
  else
    nMNPC = nen*nEL;

  coord.resize(nPts);
  param.resize(nPts);
  MMNPC.resize(nen > 0 ? nMNPC : nMNPC+nEL, 0);
  MINEX.resize(nEL,0);
  std::iota(MINEX.begin(),MINEX.end(),1);
}


bool ElementBlock::setCoor (size_t i, size_t j, Real x)
{
  if (i >= coord.size() || j >= 3) return false;

  coord[i][j] = x;
  return true;
}


bool ElementBlock::setCoor (size_t i, const Vec3& X)
{
  if (i >= coord.size()) return false;

  coord[i] = X;
  return true;
}


bool ElementBlock::setParams (size_t i, Real u, Real v, Real w)
{
  if (i >= param.size()) return false;

  param[i][0] = u;
  param[i][1] = v;
  param[i][2] = w;
  return true;
}


bool ElementBlock::setNode (size_t i, int nodeNumb)
{
  if (i >= MMNPC.size()) return false;

  MMNPC[i] = nodeNumb;
  return true;
}


bool ElementBlock::endOfElm (size_t& i)
{
  if (nen == 0)
  {
    if (i >= MMNPC.size()) return false;

    MMNPC[i++] = -1;

    if (i == MMNPC.size())
    {
      // Compute the internal element mapping due to that mixed element types
      // are stored blockwise. Same logic as in the getElements() method.
      elmIdx.resize(MINEX.size(),0);
      size_t npc, iEx, iEl = 0;
      for (size_t nenod : legalNENs)
      {
        npc = iEx = 0;
        std::vector<int>::const_iterator it, itb = MMNPC.begin();
        for (it = MMNPC.begin(); it != MMNPC.end(); ++it)
          if (*it < 0)
          {
            if (npc == nenod)
              elmIdx[iEx] = iEl++;
            itb = it + 1;
            npc = 0;
            ++iEx;
          }
          else
            ++npc;
        if (iEl >= MINEX.size())
          break;
      }
    }
  }
  return true;
}


bool ElementBlock::addLine (Real x1, Real y1, Real z1,
                            Real x2, Real y2, Real z2)
{
  if (nen != 2) return false;

  coord.push_back(Vec3(x1,y1,z1));
  coord.push_back(Vec3(x2,y2,z2));
  MMNPC.push_back(coord.size()-2);
  MMNPC.push_back(coord.size()-1);
  MINEX.push_back(MINEX.size()+1);
  return true;
}


size_t ElementBlock::addLine (size_t i1, const Vec3& X2, int elmId)
{
  if (nen != 2 || i1 >= coord.size()) return 0;

  coord.push_back(X2);
  MMNPC.push_back(i1);
  MMNPC.push_back(coord.size()-1);
  MINEX.push_back(elmId < 0 ? MINEX.size()+1 : elmId);
  return coord.size()-1;
}


void ElementBlock::merge (const ElementBlock* other,
                          std::vector<int>& nodeNums, bool uniqNodes)
{
  nodeNums.clear();
  nodeNums.reserve(other->coord.size());

  std::vector<Vec3>::const_iterator cit;
  for (const Vec3& X : other->coord)
    if (!uniqNodes || (cit = find(coord.begin(),coord.end(),X)) == coord.end())
    {
      nodeNums.push_back(coord.size());
      coord.push_back(X);
    }
    else
      nodeNums.push_back(cit - coord.begin());

  for (int inod : other->MMNPC)
    MMNPC.push_back(inod < 0 ? inod : nodeNums[inod]);

  MINEX.insert(MINEX.end(),other->MINEX.begin(),other->MINEX.end());
}


void ElementBlock::merge (const ElementBlock& other, bool uniqNodes)
{
  std::vector<int> nodes;
  this->merge(&other,nodes,uniqNodes);
}


void ElementBlock::transform (const Matrix& Tlg)
{
  if (!Tlg.empty())
    for (Vec3& X : coord)
      X = Tlg * X;
}


const Real* ElementBlock::getParam (size_t i) const
{
  return i < param.size() ? param[i].data() : nullptr;
}


bool ElementBlock::getElements (std::vector<int>& mnpc, size_t nenod) const
{
  mnpc.clear();
  if (nen > 0)
  {
    if (nen != nenod)
      return false;

    mnpc = MMNPC;
    return true;
  }

  size_t npc = 0;
  std::vector<int>::const_iterator it, itb = MMNPC.begin();
  for (it = MMNPC.begin(); it != MMNPC.end(); ++it)
    if (*it < 0)
    {
      if (npc == nenod)
        mnpc.insert(mnpc.end(),itb,itb+npc);
      itb = it + 1;
      npc = 0;
    }
    else
      ++npc;

  return !mnpc.empty();
}


utl::Point ElementBlock::getCenter (size_t i) const
{
  if (i < 1 || i > MINEX.size())
    return utl::Point();

  utl::Point XC;
  if (nen > 0)
  {
    for (size_t j = nen*i-nen; j < nen*i; j++)
      if (MMNPC[j] < (int)param.size())
      {
        const Prm3& uu = param[MMNPC[j]];
        XC += utl::Point(coord[MMNPC[j]],{uu[0],uu[1],uu[2]});
      }
      else // No spline parameters
        XC += utl::Point(coord[MMNPC[j]]);
    XC /= nen;
  }
  else
  {
    size_t elm = 0, nenod = 0;
    for (int inod : MMNPC)
      if (inod < 0)
      {
        if (++elm == i)
          break;
      }
      else if (elm == i-1)
      {
        if (inod < (int)param.size())
        {
          const Prm3& uu = param[inod];
          XC += utl::Point(coord[inod],{uu[0],uu[1],uu[2]});
        }
        else // No spline parameters
          XC += utl::Point(coord[inod]);
        ++nenod;
      }
    XC /= nenod;
  }

  return XC;
}


void ElementBlock::removeElement (size_t i)
{
  if (i < 1 || i > MINEX.size())
    return;

  std::vector<int>::iterator it1 = MMNPC.end(), it2 = MMNPC.end();
  if (nen > 0)
  {
    it2 = MMNPC.begin()+nen*i;
    it1 = it2 - nen;
  }
  else
  {
    size_t elm = 0;
    for (size_t ip = 0; ip < MMNPC.size() && it2 == MMNPC.end(); ip++)
      if (MMNPC[ip] < 0 && ++elm == i)
        it2 = MMNPC.begin() + ip+1;
      else if (elm == i-1 && it1 == MMNPC.end())
        it1 = MMNPC.begin() + ip;
  }

  MMNPC.erase(it1,it2);
  MINEX.erase(MINEX.begin()+i-1);
}


PlaneBlock::PlaneBlock (const Vec3& X0, const Vec3& X1,
                        const Vec3& X2) : ElementBlock(4)
{
  this->resize(2,2);
  this->setCoor(0,X0);
  this->setCoor(1,X1);
  this->setCoor(2,X1+X2-X0);
  this->setCoor(3,X2);
  for (int i = 0; i < 4; i++)
    this->setNode(i,i);
  this->setElmId(1,1);
}


CubeBlock::CubeBlock (const Vec3& X0, double dX) : ElementBlock(8)
{
  dX *= 0.5;
  this->resize(2,2,2);
  this->setCoor(0, X0.x-dX, X0.y-dX, X0.z-dX);
  this->setCoor(1, X0.x+dX, X0.y-dX, X0.z-dX);
  this->setCoor(2, X0.x+dX, X0.y+dX, X0.z-dX);
  this->setCoor(3, X0.x-dX, X0.y+dX, X0.z-dX);
  this->setCoor(4, X0.x-dX, X0.y-dX, X0.z+dX);
  this->setCoor(5, X0.x+dX, X0.y-dX, X0.z+dX);
  this->setCoor(6, X0.x+dX, X0.y+dX, X0.z+dX);
  this->setCoor(7, X0.x-dX, X0.y+dX, X0.z+dX);
  for (int i = 0; i < 8; i++)
    this->setNode(i,i);
  this->setElmId(1,1);
}


SphereBlock::SphereBlock (const Vec3& X0, double R,
                          size_t nTheta, size_t nPhi) : ElementBlock(4)
{
  this->unStructResize(nPhi*nTheta,2+(nPhi-1)*nTheta);

  const double dTheta = M_PI*2.0/static_cast<double>(nTheta);
  const double dPhi   = M_PI/static_cast<double>(nPhi);

  double theta, phi;
  size_t m, n, ip = 2;

  this->setCoor(0,X0.x,X0.y,X0.z+R);
  this->setCoor(1,X0.x,X0.y,X0.z-R);
  for (theta = 0.0, n = 0; n < nTheta; n++, theta += dTheta)
  {
    double Rct = R*cos(theta);
    double Rst = R*sin(theta);
    for (phi = dPhi, m = 1; m < nPhi; m++, ip++, phi += dPhi)
      this->setCoor(ip,X0.x+Rct*sin(phi),X0.y+Rst*sin(phi),X0.z+R*cos(phi));
  }

  for (n = ip = 0; n+1 < nTheta; n++)
  {
    this->setNode(ip++,0);
    this->setNode(ip++,nPhi* n   -n+2);
    this->setNode(ip++,nPhi*(n+1)-n+1);
    this->setNode(ip++,0);
    for (m = 1; m+1 < nPhi; m++)
    {
      this->setNode(ip++,nPhi* n   -n+m+1);
      this->setNode(ip++,nPhi* n   -n+m+2);
      this->setNode(ip++,nPhi*(n+1)-n+m+1);
      this->setNode(ip++,nPhi*(n+1)-n+m);
    }
    this->setNode(ip++,nPhi*(n+1)-n);
    this->setNode(ip++,1);
    this->setNode(ip++,1);
    this->setNode(ip++,nPhi*(n+2)-n-1);
  }

  this->setNode(ip++,0);
  this->setNode(ip++,(nPhi-1)*(nTheta-1)+2);
  this->setNode(ip++,2);
  this->setNode(ip++,0);
  for (m = 1; m+1 < nPhi; m++)
  {
    this->setNode(ip++,(nPhi-1)*(nTheta-1)+m+1);
    this->setNode(ip++,(nPhi-1)*(nTheta-1)+m+2);
    this->setNode(ip++,m+2);
    this->setNode(ip++,m+1);
  }
  this->setNode(ip++,(nPhi-1)*nTheta+1);
  this->setNode(ip++,1);
  this->setNode(ip++,1);
  this->setNode(ip++,nPhi);
}


CylinderBlock::CylinderBlock (const Vec3& X0, const Vec3& X1,
                              double R, size_t nTheta) : ElementBlock(4)
{
  this->unStructResize(3*nTheta,2+2*nTheta);

  const double dTheta = M_PI*2.0/static_cast<double>(nTheta);

  double theta;
  size_t n, ip;
  Tensor Tlg(X1-X0);

  this->setCoor(0,X0);
  this->setCoor(1+nTheta,X1);
  for (theta = 0.0, n = 0; n < nTheta; n++, theta += dTheta)
  {
    Vec3 X(R*cos(theta),R*sin(theta),0.0);
    this->setCoor(1+n,X0+Tlg*X);
    this->setCoor(2+nTheta+n,X1+Tlg*X);
  }

  for (n = ip = 0; n+1 < nTheta; n++)
  {
    this->setNode(ip++,0);
    this->setNode(ip++,2+n);
    this->setNode(ip++,1+n);
    this->setNode(ip++,0);

    this->setNode(ip++,2+n);
    this->setNode(ip++,1+n);
    this->setNode(ip++,2+nTheta+n);
    this->setNode(ip++,3+nTheta+n);

    this->setNode(ip++,1+nTheta);
    this->setNode(ip++,2+nTheta+n);
    this->setNode(ip++,3+nTheta+n);
    this->setNode(ip++,1+nTheta);
  }

  this->setNode(ip++,0);
  this->setNode(ip++,1);
  this->setNode(ip++,nTheta);
  this->setNode(ip++,0);

  this->setNode(ip++,1);
  this->setNode(ip++,nTheta);
  this->setNode(ip++,1+nTheta+nTheta);
  this->setNode(ip++,2+nTheta);

  this->setNode(ip++,1+nTheta);
  this->setNode(ip++,1+nTheta+nTheta);
  this->setNode(ip++,2+nTheta);
  this->setNode(ip++,1+nTheta);
}
