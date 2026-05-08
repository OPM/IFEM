//==============================================================================
//!
//! \file MeshUtils.h
//!
//! \date Feb 16 2015
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Various helpers for mesh calculations
//!
//==============================================================================

#ifndef MESHUTILS_H_
#define MESHUTILS_H_

#include "MatVec.h"

class SIMbase;

namespace MeshUtils
{
  //! \brief Compute element aspect ratios for a mesh
  //! \param[out] elmAspects The element aspect ratios
  //! \param[in] model The model holding the mesh
  //! \param[in] displacement A displacement to apply to mesh coordinates
  bool computeAspectRatios(std::vector<double>& elmAspects,
                           const SIMbase& model,
                           const Vector& displacement=Vector());

  //! \brief Compute element skewness for a mesh
  //! \param[out] elmSkewness The element skewness values
  //! \param[in] model The model holding the mesh
  //! \param[in] displacement A displacement to apply to mesh coordinates
  bool computeMeshSkewness(std::vector<double>& elmSkewness,
                           const SIMbase& model,
                           const Vector& displacement=Vector());
}

#endif
