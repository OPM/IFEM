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


template<size_t Dim>
BasisFunctionCache<Dim>::BasisFunctionCache ()
  : mainQ(std::make_shared<Quadrature>()),
    reducedQ(std::make_shared<Quadrature>())
{
}

template<size_t Dim>
BasisFunctionCache<Dim>::BasisFunctionCache (const BasisFunctionCache<Dim>& rhs)
  : basis(rhs.basis), mainQ(rhs.mainQ), reducedQ(rhs.reducedQ)
{
}


template<size_t Dim>
bool BasisFunctionCache<Dim>::init (int nd)
{
  if (!values.empty() && nd <= nderiv)
    return true;

  values.clear();
  valuesRed.clear();
  this->internalCleanup();
  nderiv = nd;
  if (!this->internalInit())
    return false;

  if (ASM::cachePolicy == ASM::NO_CACHE) {
    this->resizeThreadBuffers();
    return true;
  }

  values.resize(nTotal);
  if (this->hasReduced())
    valuesRed.resize(nTotalRed);
  if (ASM::cachePolicy != ASM::ON_THE_FLY)
    this->calculateAll();

  return true;
}


template<size_t Dim>
void BasisFunctionCache<Dim>::finalizeAssembly ()
{
  if (ASM::cachePolicy == ASM::PRE_CACHE) {
    values.clear();
    valuesRed.clear();
    internalCleanup();
  }
}


template<size_t Dim>
const BasisFunctionVals& BasisFunctionCache<Dim>::
getVals (size_t el, size_t gp, bool reduced)
{
  std::vector<BasisFunctionVals>& vals = reduced ? valuesRed : values;
  if (ASM::cachePolicy == ASM::NO_CACHE) {
#ifdef USE_OPENMP
    size_t idx = omp_get_thread_num();
#else
    size_t idx = 0;
#endif
    vals[idx] = this->calculatePt(el,gp,reduced);
    return vals[idx];
  }

  size_t idx = this->index(el, gp, reduced);
  if (ASM::cachePolicy == ASM::ON_THE_FLY && vals[idx].N.empty())
    vals[idx] = this->calculatePt(el,gp,reduced);

  return vals[idx];
}


template<size_t Dim>
double BasisFunctionCache<Dim>::getParam (int dir, size_t el,
                                          size_t gp, bool reduced) const
{
  const Quadrature& q = reduced ? *reducedQ : *mainQ;

  return q.gpar[dir][gp+q.ng[dir]*el];
}


template<size_t Dim>
size_t BasisFunctionCache<Dim>::index (size_t el, size_t gp, bool reduced) const
{
  const Quadrature& q = reduced ? *reducedQ : *mainQ;

  if constexpr (Dim == 2)
    return el*q.ng[0]*q.ng[1] + gp;
  else if constexpr (Dim == 3)
    return el*q.ng[0]*q.ng[1]*q.ng[2] + gp;
}


template<size_t Dim>
std::array<size_t,Dim>
BasisFunctionCache<Dim>::gpIndex (size_t gp, bool reduced) const
{
  const Quadrature& q = reduced ? *reducedQ : *mainQ;
  if constexpr (Dim == 1)
    return {gp};
  else if constexpr (Dim == 2)
    return { gp % q.ng[0], gp / q.ng[0] };
  else if constexpr (Dim == 3)
    return { gp % q.ng[0], (gp / q.ng[0]) % q.ng[1],  gp / (q.ng[0]*q.ng[1]) };
}


template<size_t Dim>
void BasisFunctionCache<Dim>::resizeThreadBuffers ()
{
  if (ASM::cachePolicy == ASM::NO_CACHE) {
#ifdef USE_OPENMP
    size_t size = omp_get_max_threads();
#else
    size_t size = 1;
#endif
    values.resize(size);
    if (this->hasReduced())
      valuesRed.resize(size);
  }
}


template class BasisFunctionCache<2>;
template class BasisFunctionCache<3>;
