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

#include <vector>

class ModelGenerator;
class TiXmlElement;


/*!
  \brief Inherit this class to equip your SIM with multi-patch model generators.
*/

template<class Dim>
class SIMMultiPatchModelGen : public Dim
{
public:
  //! \brief Constructor for standard problems.
  //! \param[in] n1 Dimension of the primary solution field
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  explicit SIMMultiPatchModelGen(int n1, bool checkRHS = false) : Dim(n1,checkRHS) {}

  //! \brief Constructor for mixed problems.
  //! \param[in] unf Dimension of the primary solution field
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  explicit SIMMultiPatchModelGen(const std::vector<unsigned char>& unf,
                                 bool checkRHS = false) : Dim(unf,checkRHS) {}

  //! \brief Empty destructor.
  virtual ~SIMMultiPatchModelGen() {}

protected:
  //! \brief Instantiates a FEM model generator for a multi-patch model.
  //! \param[in] geo XML element containing geometry definition
  virtual ModelGenerator* getModelGenerator(const TiXmlElement* geo) const;
};

#endif
