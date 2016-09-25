// $Id$
//==============================================================================
//!
//! \file SIMMultiPatchModelGen.h
//!
//! \date Sep 5 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for simulators equipped with multi-patch model generators.
//!
//==============================================================================

#ifndef _SIM_MULTI_PATCH_MODEL_GEN_H_
#define _SIM_MULTI_PATCH_MODEL_GEN_H_

#include "SIMbase.h"

class ModelGenerator;


/*!
  \brief Inherit this class to equip your SIM with multi-patch model generators.
*/

template<class Dim>
class SIMMultiPatchModelGen : public Dim
{
public:
  //! \brief Constructor for standard problems.
  //! \param[in] n1 Number of fields
  SIMMultiPatchModelGen(int n1) : Dim(n1) {}

  //! \brief Constructor for mixed problems.
  //! \param[in] unf Number of fields on bases
  SIMMultiPatchModelGen(const SIMbase::CharVec& unf) : Dim(unf) {}

  //! \brief Empty destructor.
  virtual ~SIMMultiPatchModelGen() {}

protected:
  //! \brief Instantiates a generator for the finite element model.
  //! \param[in] geo XML element containing geometry defintion
  virtual ModelGenerator* getModelGenerator(const TiXmlElement* geo) const;
};

#endif
