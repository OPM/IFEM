// $Id$
//==============================================================================
//!
//! \file ModelGenerator.h
//!
//! \date Sep 2 2016
//!
//! \author Arne Morten Kvarving / SINTEF
//!
//! \brief Base class for model generators for NURBS-based FEM simulators.
//!
//==============================================================================

#ifndef _MODEL_GENERATOR_H
#define _MODEL_GENERATOR_H

#include "SIMdependency.h"
#include "TopologySet.h"
#include <string>

class SIMbase;
class TiXmlElement;

/*!
  \brief Base class for model generators for FEM simulators.
*/

class ModelGenerator
{
public:
  //! \brief Constructor initializes common members
  //!\ param elem XML element to parse
  ModelGenerator(const TiXmlElement* elem);

  //! \brief Empty destructor.
  virtual ~ModelGenerator() {}

  //! \brief Creates a geometry.
  //! \param[in] sim SIM with patch read function to use
  virtual SIMdependency::PatchVec createGeometry(const SIMbase& sim) const = 0;

  //! \brief Creates topology for geometry.
  //! \param[in] geo XML element containing geometry defintion
  //! \param sim Simulator to apply topology to
  virtual bool createTopology(SIMbase& sim) const = 0;

  //! \brief Creates topology sets for geometry.
  //! \param[in] sim Simulator with patch ownerships
  virtual TopologySet createTopologySets(const SIMbase& sim) const = 0;

protected:
  bool sets; //!< Whether to generate topologysets or not
  const TiXmlElement* geo; //!< Pointer to xml element describing geometry
};


/*!
  \brief Default model generator for 1D FEM simulators.
  \details Generates a line.
*/

class DefaultGeometry1D : public ModelGenerator {
public:
  //! \brief The constructor forwards to the base class.
  //! \param[in] geo XML element containing geometry defintion
  DefaultGeometry1D(const TiXmlElement* geo) : ModelGenerator(geo) {}

  //! \brief Creates a 1D single-patch geometry.
  //! \param[in] sim SIM with patch read function to use
  SIMdependency::PatchVec createGeometry(const SIMbase& sim) const override;

  //! \brief Creates the topology
  //! \details No topology information for single patch models
  bool createTopology(SIMbase&) const override
  { return true; }

  //! \brief Creates topology sets for geometry.
  TopologySet createTopologySets(const SIMbase&) const override;

protected:
  //! \brief Generates the G2 description of the geometry.
  //! \param nsd Number of spatial dimension
  std::string createG2 (int nsd = 2) const;
};


/*!
  \brief Default model generator for 2D FEM simulators.
  \details Generates a rectangle.
*/

class DefaultGeometry2D : public ModelGenerator {
public:
  //! \brief The constructor forwards to the base class.
  //! \param[in] geo XML element containing geometry defintion
  DefaultGeometry2D(const TiXmlElement* geo) : ModelGenerator(geo) {}

  //! \brief Creates a 2D rectangular single-patch geometry.
  //! \param[in] sim SIM with patch read function to use
  SIMdependency::PatchVec createGeometry(const SIMbase& sim) const override;

  //! \brief Creates the topology
  //! \param sim Simulator to apply topology to
  //! \details No topology information for single patch models
  bool createTopology(SIMbase&) const override
  { return true; }

  //! \brief Creates topology sets for geometry.
  TopologySet createTopologySets(const SIMbase&) const override;

protected:
  //! \brief Generates the G2 description of the geometry.
  //! \param nsd Number of spatial dimension
  std::string createG2 (int nsd = 3) const;
};


/*!
  \brief Default model generator for 3D FEM simulators.
  \details Generates a hexahedra.
*/

class DefaultGeometry3D : public ModelGenerator {
public:
  //! \brief The constructor forwards to the base class.
  //! \param[in] geo XML element containing geometry defintion
  DefaultGeometry3D(const TiXmlElement* geo) : ModelGenerator(geo) {}

  //! \brief Creates a 3D hexahedral single-patch geometry.
  //! \param[in] sim Simulator with patch read function to use
  SIMdependency::PatchVec createGeometry(const SIMbase& sim) const override;

  //! \brief Creates the topology
  //! \param sim Simulator to apply topology to
  //! \details No topology information for single patch models
  bool createTopology(SIMbase&) const override
  { return true; }

  //! \brief Creates topology sets for geometry.
  TopologySet createTopologySets(const SIMbase&) const override;

protected:
  //! \brief Generates the G2 description of the geometry.
  std::string createG2 (int = 3) const;
};

#endif
