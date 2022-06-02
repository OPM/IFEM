// $Id$
//==============================================================================
//!
//! \file BasisFunctionCache.h
//!
//! \date Jun 1 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basis function cache.
//!
//==============================================================================

#ifndef _BASIS_FUNCTION_CACHE_H
#define _BASIS_FUNCTION_CACHE_H

#include "ASMenums.h"
#include "matrixnd.h"
#include "MatVec.h"


struct BasisFunctionVals {
  Vector N;
  Matrix dNdu;
  Matrix3D d2Ndu2;
  Matrix4D d3Ndu3;
};


class BasisFunctionCache
{
public:
  BasisFunctionCache(ASM::CachePolicy plcy) : policy(plcy) {}
  void finalizeAssembly();

  bool init(int nd);

  const BasisFunctionVals& getVals(size_t el, size_t gp, bool reduced = false);

  virtual bool hasReduced() const  = 0;

  virtual size_t index (size_t el, size_t gp, bool reduced) const = 0;

protected:
  virtual bool internalInit() = 0;
  virtual BasisFunctionVals calculatePt(size_t el, size_t gp,
                                        bool reduced = false) const = 0;
  virtual void calculateAll() = 0;

  ASM::CachePolicy policy;
  std::vector<BasisFunctionVals> values;
  std::vector<BasisFunctionVals> valuesRed;
  int nderiv = 0;
  size_t nTotal = 0;
  size_t nTotalRed = 0;
};

#endif
