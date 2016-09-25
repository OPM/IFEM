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

class SIMinput;
class TiXmlElement;


/*!
  \brief Base class for model generators for FEM simulators.
*/

class ModelGenerator
{
public:
  //! \brief The constructor initializes the common members.
  //!\ param elem XML element to parse
  ModelGenerator(const TiXmlElement* elem) : geo(elem) {}
  //! \brief Empty destructor.
  virtual ~ModelGenerator() {}

  //! \brief Creates a geometry.
  //! \param[in] m Simulator object with patch read function to use
  virtual SIMdependency::PatchVec createGeometry(const SIMinput& m) const;

  //! \brief Creates topology for multi-patch geometries.
  virtual bool createTopology(SIMinput&) const { return true; }

  //! \brief Creates topology sets for geometry.
  //! \param[in] sim Simulator object with patch ownerships
  virtual TopologySet createTopologySets(const SIMinput& sim) const = 0;

protected:
  //! \brief Generates the G2 description of the geometry.
  //! \param nsd Number of spatial dimension
  virtual std::string createG2(int nsd) const { return ""; }

  //! \brief Returns \e true if topology sets is to be generated.
  bool topologySets() const;

protected:
  const TiXmlElement* geo; //!< Pointer to xml element describing geometry
};


/*!
  \brief Default model generator for 1D FEM simulators.
  \details Generates a line domain.
*/

class DefaultGeometry1D : public ModelGenerator
{
public:
  //! \brief The constructor forwards to the base class.
  //! \param[in] geo XML element containing geometry defintion
  DefaultGeometry1D(const TiXmlElement* geo) : ModelGenerator(geo) {}
  //! \brief Empty destructor.
  virtual ~DefaultGeometry1D() {}

  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput&) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  //! \param nsd Number of spatial dimension
  virtual std::string createG2(int nsd) const;
};


/*!
  \brief Default model generator for 2D FEM simulators.
  \details Generates a quadrilateral domain.
*/

class DefaultGeometry2D : public ModelGenerator
{
public:
  //! \brief The constructor forwards to the base class.
  //! \param[in] geo XML element containing geometry defintion
  DefaultGeometry2D(const TiXmlElement* geo) : ModelGenerator(geo) {}
  //! \brief Empty destructor.
  virtual ~DefaultGeometry2D() {}

  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput&) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  //! \param nsd Number of spatial dimension
  virtual std::string createG2(int nsd) const;
};


/*!
  \brief Default model generator for 3D FEM simulators.
  \details Generates a hexahedral domain.
*/

class DefaultGeometry3D : public ModelGenerator
{
public:
  //! \brief The constructor forwards to the base class.
  //! \param[in] geo XML element containing geometry defintion
  DefaultGeometry3D(const TiXmlElement* geo) : ModelGenerator(geo) {}
  //! \brief Empty destructor.
  virtual ~DefaultGeometry3D() {}

  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput&) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int) const;
};

#endif
