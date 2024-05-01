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
