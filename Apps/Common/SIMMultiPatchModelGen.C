// $Id$
//==============================================================================
//!
//! \file SIMMultiPatchModelGen.C
//!
//! \date Sep 5 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for simulators equipped with multi-patch model generators.
//!
//==============================================================================

#include "SIMMultiPatchModelGen.h"
#include "MultiPatchModelGenerator.h"
#include "IFEM.h"
#include "SIM2D.h"
#include "SIM3D.h"


template<>
ModelGenerator* SIMMultiPatchModelGen<SIM2D>::getModelGenerator (const TiXmlElement* geo) const
{
  IFEM::cout <<"  Using multi-patch model generator" << std::endl;
  return new MultiPatchModelGenerator2D(geo);
}


template<>
ModelGenerator* SIMMultiPatchModelGen<SIM3D>::getModelGenerator (const TiXmlElement* geo) const
{
  IFEM::cout <<"  Using multi-patch model generator" << std::endl;
  return new MultiPatchModelGenerator3D(geo);
}
