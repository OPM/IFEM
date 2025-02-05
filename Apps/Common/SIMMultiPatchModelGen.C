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
#include "SIM1D.h"
#include "SIM2D.h"
#include "SIM3D.h"
#include "IFEM.h"


//! \brief Template specialization for 1D.
template<>
ModelGenerator* SIMMultiPatchModelGen<SIM1D>::getModelGenerator(const tinyxml2::XMLElement* geo) const
{
  if (!useMultiPatchGen(geo))
    return this->SIM1D::getModelGenerator(geo);

  IFEM::cout <<"  Using 1D multi-patch model generator."<< std::endl;
  return new MultiPatchModelGenerator1D(geo);
}


//! \brief Template specialization for 2D.
template<>
ModelGenerator* SIMMultiPatchModelGen<SIM2D>::getModelGenerator(const tinyxml2::XMLElement* geo) const
{
  if (!useMultiPatchGen(geo))
    return this->SIM2D::getModelGenerator(geo);

  IFEM::cout <<"  Using 2D multi-patch model generator."<< std::endl;
  return new MultiPatchModelGenerator2D(geo);
}


//! \brief Template specialization for 3D.
template<>
ModelGenerator* SIMMultiPatchModelGen<SIM3D>::getModelGenerator(const tinyxml2::XMLElement* geo) const
{
  if (!useMultiPatchGen(geo))
    return this->SIM3D::getModelGenerator(geo);

  IFEM::cout <<"  Using 3D multi-patch model generator."<< std::endl;
  return new MultiPatchModelGenerator3D(geo);
}
