// $Id$
//==============================================================================
//!
//! \file ASMLagBase.C
//!
//! \date May 1 2024
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Common base class for %Lagrange FE models.
//!
//==============================================================================

#include "ASMLagBase.h"


ASMLagBase::ASMLagBase (const ASMLagBase& patch, bool cpCoord) : coord(myCoord)
{
  if (cpCoord) myCoord = patch.coord;
}


bool ASMLagBase::nodalField (Matrix& field, const Vector& sol, size_t nno) const
{
  size_t nPnts = coord.size();
  size_t nComp = sol.size() / nno;
  if (nno < nPnts || nComp*nno != sol.size())
    return false;

  field.resize(nComp,nPnts);
  const double* u = sol.ptr();
  for (size_t iPt = 1; iPt <= nPnts; iPt++, u += nComp)
    field.fillColumn(iPt,u);

  return true;
}


Vec3 ASMLagBase::getGeometricCenter (const std::vector<int>& MNPC) const
{
  Vec3 X0;

  for (int inod : MNPC)
    X0 += coord[inod];
  X0 *= 1.0 / static_cast<double>(MNPC.size());

  return X0;
}


bool ASMLagBase::updateCoords (const Vector& displ, unsigned char nsd)
{
  if (displ.size() != nsd*coord.size())
  {
    std::cerr <<" *** ASMLagBase::updateCoords: Invalid dimension "
              << displ.size() <<" on displ, should be "
              << nsd*coord.size() << std::endl;
    return false;
  }

  const double* dpt = displ.ptr();
  for (Vec3& XYZ : myCoord)
  {
    XYZ += RealArray(dpt,dpt+nsd);
    dpt += nsd;
  }

  return true;
}


void ASMLagBase::updateOrigin (const Vec3& origin)
{
  for (Vec3& XYZ : myCoord)
    XYZ += origin;
}
