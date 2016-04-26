// $Id$
//==============================================================================
//!
//! \file SIMgeneric.C
//!
//! \date Aug 28 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Generic SIM class with some added functionalities.
//!
//==============================================================================

#include "SIMgeneric.h"
#include "ASMbase.h"


void SIMgeneric::createDefaultModel ()
{
  if (!myModel.empty()) return;

  nGlPatches = 1;
  myModel.resize(1,this->createDefaultGeometry(nullptr));
}


Vector SIMgeneric::getSolution (const Vector& psol, const double* par,
                                int deriv, int patch) const
{
  if (psol.empty() || !par || opt.discretization < ASM::Spline)
    return Vector();

  ASMbase* pch = this->getPatch(this->getLocalPatchIndex(patch));
  if (!pch) return Vector();

  size_t ndim = pch->getNoParamDim();
  std::vector<RealArray> params(ndim);
  for (size_t i = 0; i < ndim; i++)
    params[i].resize(1,par[i]);

  Matrix tmpVal;
  Vector localVec;
  pch->extractNodeVec(psol,localVec);
  if (!pch->evalSolution(tmpVal,localVec,&params.front(),false,deriv))
    return Vector();

  return tmpVal.getColumn(1);
}


int SIMgeneric::evalPoint (const double* xi, Vec3& X, double* param,
                           int patch) const
{
  ASMbase* pch = this->getPatch(this->getLocalPatchIndex(patch));
  if (!pch) return -1;

  double dummy[3];
  return pch->evalPoint(xi, param ? param : dummy, X);
}
