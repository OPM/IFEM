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
#include "ModelGenerator.h"
#include "ASMbase.h"


ASMbase* SIMgeneric::createDefaultModel ()
{
  if (!myModel.empty()) return nullptr;

  ModelGenerator* gen = this->getModelGenerator(nullptr);
  if (!gen) return nullptr;

  bool okGen = gen->createGeometry(*this);
  nGlPatches = myModel.size();
  delete gen;

  return okGen && nGlPatches > 0 ? myModel.front() : nullptr;
}


Vector SIMgeneric::getSolution (const Vector& psol, const double* par,
                                int deriv, int patch) const
{
  if (psol.empty() || !par || opt.discretization < ASM::Spline)
    return Vector();

  ASMbase* pch = this->getPatch(patch,true);
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
                           int patch, bool global) const
{
  ASMbase* pch = this->getPatch(patch,true);
  if (!pch) return -1;

  double dummy[3] = {};
  int inod = pch->evalPoint(xi, param ? param : dummy, X);
  return inod > 0 && global ? pch->getNodeID(inod) : inod;
}


int SIMgeneric::findElementContaining (const double* param,
                                       int patch, bool global) const
{
  ASMbase* pch = this->getPatch(patch,true);
  if (!pch) return -1;

  int iel = pch->findElementContaining(param);
  return iel > 0 && global ? pch->getElmID(iel) : iel;
}
