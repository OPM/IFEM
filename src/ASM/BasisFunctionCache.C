// $Id$
//==============================================================================
//!
//! \file BasisFunctionCache.C
//!
//! \date Jun 1 2022
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Basis function cache.
//!
//==============================================================================

#include "BasisFunctionCache.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif


bool BasisFunctionCache::init (int nd)
{
  if (!values.empty() && nd <= nderiv)
    return true;

  values.clear();
  valuesRed.clear();
  nderiv = nd;
  if (!this->internalInit())
    return false;

  if (policy == ASM::NO_CACHE) {
#ifdef USE_OPENMP
    size_t size = omp_get_max_threads();
#else
    size_t size = 1;
#endif
    values.resize(size);
    if (this->hasReduced())
      valuesRed.resize(size);
    return true;
  }
  values.resize(nTotal);
  valuesRed.resize(nTotalRed);
  if (policy != ASM::ON_THE_FLY)
    this->calculateAll();

  return true;
}


void BasisFunctionCache::finalizeAssembly ()
{
  if (policy == ASM::PRE_CACHE)
    values.clear();
}


const BasisFunctionVals& BasisFunctionCache::getVals (size_t el, size_t gp, bool reduced)
{
  std::vector<BasisFunctionVals>& vals = reduced ? valuesRed : values;
  if (policy == ASM::NO_CACHE) {
#ifdef USE_OPENMP
    size_t idx = omp_get_thread_num();
#else
    size_t idx = 0;
#endif
    vals[idx] = this->calculatePt(el,gp,reduced);
    return vals[idx];
  }

  size_t idx = this->index(el, gp, reduced);
  if (policy == ASM::ON_THE_FLY && vals[idx].N.empty())
    vals[idx] = this->calculatePt(el,gp,reduced);

  return vals[idx];
}
