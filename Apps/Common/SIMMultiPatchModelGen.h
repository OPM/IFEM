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
#include "tinyxml2.h"

class ModelGenerator;


/*!
  \brief Inherit this class to equip your SIM with multi-patch model generators.
*/

template<class Dim>
class SIMMultiPatchModelGen : public Dim
{
protected:
  //! \brief Constructor for standard problems.
  //! \param[in] n1 Dimension of the primary solution field
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  explicit SIMMultiPatchModelGen(int n1, bool checkRHS = false) : Dim(n1,checkRHS) {}

  //! \brief Constructor for mixed problems.
  //! \param[in] unf Dimension of the primary solution field
  //! \param[in] checkRHS If \e true, ensure the model is in a right-hand system
  explicit SIMMultiPatchModelGen(const std::vector<unsigned char>& unf,
                                 bool checkRHS = false) : Dim(unf,checkRHS) {}

  //! \brief Instantiates a FEM model generator for a multi-patch model.
  //! \param[in] geo XML element containing geometry definition
  virtual ModelGenerator* getModelGenerator(const tinyxml2::XMLElement* geo) const;

private:
  //! \brief Private helper defining if the multi-patch generator must be used.
  //! \param[in] geo XML element containing geometry definition
  static bool useMultiPatchGen(const tinyxml2::XMLElement* geo)
  {
    if (!geo) return false;

    // If any of these attribute keywords are present on the geometry tag,
    // the multi-patch model generator must be used.
    // Otherwise, use the standard single-patch generators.
    static const char* keywords[9] = {
      "nx", "ny", "nz",
      "periodic_x", "periodic_x", "periodic_x",
      "subdivision", "sets", nullptr
    };

    for (const char** q = keywords; *q; ++q)
      if (geo->Attribute(*q)) return true;

    return geo->FirstChildElement("subdivision");
  }
};

#endif
