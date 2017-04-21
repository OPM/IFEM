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

#include "ModelGenerator.h"
#include "SIMMultiPatchModelGen.h"
#include "MultiPatchModelGenerator.h"
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "IFEM.h"


template<>
ModelGenerator* SIMMultiPatchModelGen<SIM1D>::getModelGenerator(const TiXmlElement* geo) const
{
  IFEM::cout <<"  Using 1D multi-patch model generator."<< std::endl;
  return new MultiPatchModelGenerator1D(geo);
}


template<>
ModelGenerator* SIMMultiPatchModelGen<SIM2D>::getModelGenerator(const TiXmlElement* geo) const
{
  IFEM::cout <<"  Using 2D multi-patch model generator."<< std::endl;
  return new MultiPatchModelGenerator2D(geo);
}


template<>
ModelGenerator* SIMMultiPatchModelGen<SIM3D>::getModelGenerator(const TiXmlElement* geo) const
{
  IFEM::cout <<"  Using 3D multi-patch model generator."<< std::endl;
  return new MultiPatchModelGenerator3D(geo);
}
