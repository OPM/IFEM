// $Id$
//==============================================================================
//!
//! \file MultiPatchModelGenerator.h
//!
//! \date Sep 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Multi-patch model generators for NURBS-based FEM simulators.
//!
//==============================================================================

#ifndef _MULTIPATCH_MODEL_GENERATOR_H
#define _MULTIPATCH_MODEL_GENERATOR_H

#include "ModelGenerator.h"
#include <string>


/*!
  \brief 2D multi-patch model generator for FEM simulators.
  \details Generate a rectangle split in a given number of blocks.
*/

class MultiPatchModelGenerator2D : public ModelGenerator
{
public:
  //! \brief Constructor initializes common members.
  //!\ param[in] elem XML element to parse
  MultiPatchModelGenerator2D(const TiXmlElement* elem);

  //! \brief Empty destructor.
  virtual ~MultiPatchModelGenerator2D() {}

  //! \brief Creates a geometry.
  //! \param[in] sim SIM with patch read function to use
  SIMdependency::PatchVec createGeometry(const SIMbase& sim) const override;

  //! \brief Creates topology for geometry.
  //! \param sim Simulator to apply topology to
  bool createTopology(SIMbase& sim) const override;

  //! \brief Creates topology sets for geometry.
  TopologySet createTopologySets(const SIMbase& sim) const override;

protected:
  //! \brief Generates the G2 description of the geometry.
  //! \param nsd Number of spatial dimension
  std::string createG2 (int nsd = 2) const;

  int nx; //!< Number of blocks in x
  int ny; //!< Number of blocks in y
  int periodic_x; //!< If non-zero, make model periodic in x for given bases
  int periodic_y; //!< If non-zero, make model periodic in y for given bases
};


/*!
 \brief 3D multi-patch model generator for FEM simulators.
 \details Generates a hexahedra split in a given number of blocks.
 */

class MultiPatchModelGenerator3D : public ModelGenerator
{
public:
  //! \brief Constructor initializes common members.
  //! \param[in] elem XML element to parse
  MultiPatchModelGenerator3D(const TiXmlElement* geo);

  //! \brief Empty destructor.
  virtual ~MultiPatchModelGenerator3D() {}

  //! \brief Creates a geometry.
  //! \param[in] sim SIM with patch read function to use
  SIMdependency::PatchVec createGeometry(const SIMbase& sim) const override;

  //! \brief Creates topology for geometry.
  //! \param[in] geo XML element containing geometry defintion
  //! \param sim Simulator to apply topology to
  bool createTopology(SIMbase& sim) const override;

  //! \brief Creates topology sets for geometry.
  //! \param[in] SIM Simulator with patch ownerships
  virtual TopologySet createTopologySets(const SIMbase& sim) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  //! \param nsd Number of spatial dimension
  std::string createG2 (int nsd = 3) const;

  int nx; //!< Number of blocks in x
  int ny; //!< Number of blocks in y
  int nz; //!< Number of blocks in z
  int periodic_x; //!< If non-zero, make model periodic in x for given bases
  int periodic_y; //!< If non-zero, make model periodic in y for given bases
  int periodic_z; //!< If non-zero, make model periodic in z for given bases
};

#endif
