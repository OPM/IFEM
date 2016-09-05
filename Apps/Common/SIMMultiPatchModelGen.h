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

#ifndef _SIM_MULTPATCH_MODEL_GEN_H_
#define _SIM_MULTIPATCH_MODEL_GEN_H_

#include "SIMbase.h"

class TiXmlElement;


//! \brief Inherit this class to equip your SIM multi-patch model generators.
template<class Dim>
class SIMMultiPatchModelGen : public Dim {
public:
  //! \brief Default constructor.
  //! \param[in] n1 Number of fields
  //! \param[in] check If \e true, ensure the model is in a right-hand system
  SIMMultiPatchModelGen(int n1,
                        bool checkRHS=false) :
    Dim(n1, checkRHS) {}

  //! \brief Default constructor.
  //! \param[in] unf Number of fields on bases
  //! \param[in] check If \e true, ensure the model is in a right-hand system
  SIMMultiPatchModelGen(const SIMbase::CharVec& unf,
                        bool checkRHS=false) :
    Dim(unf,checkRHS) {}

  //! \brief Empty destructor
  virtual ~SIMMultiPatchModelGen() {}

protected:
  //! \brief Instantiate a generator for the finite element model.
  //! \param[in] geo XML element containing geometry defintion
  ModelGenerator* createModelGenerator(const TiXmlElement* geo) const override;
};

#endif
