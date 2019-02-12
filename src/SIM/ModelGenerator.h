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

#include "TopologySet.h"
#include <vector>
#include <string>

class SIMinput;
class ASMbase;
class TiXmlElement;


/*!
  \brief Base class for model generators for FEM simulators.
*/

class ModelGenerator
{
public:
  //! \brief The constructor initializes the common members.
  //! \param[in] elem XML element containing geometry definition
  explicit ModelGenerator(const TiXmlElement* elem) : geo(elem) {}
  //! \brief Empty destructor.
  virtual ~ModelGenerator() {}

  //! \brief Creates a geometry.
  //! \param[in] m Simulator object with patch read function to use
  virtual std::vector<ASMbase*> createGeometry(const SIMinput& m) const;

  //! \brief Creates topology for multi-patch geometries.
  virtual bool createTopology(SIMinput&) const { return true; }

  //! \brief Creates topology sets for geometry.
  //! \param[in] sim Simulator object with patch ownerships
  virtual TopologySet createTopologySets(const SIMinput& sim) const = 0;

protected:
  //! \brief Generates the G2 description of the geometry.
  //! \param[in] nsd Number of spatial dimension
  //! \param[in] rational If \e true, create a NURBS geometry basis
  virtual std::string createG2(int nsd, bool rational) const { return ""; }

  //! \brief Returns \e true if topology sets is to be generated.
  bool topologySets() const;

protected:
  const TiXmlElement* geo; //!< Pointer to XML element describing geometry
};


/*!
  \brief Default model generator for 1D FEM simulators.
  \details Generates a line domain.
*/

class DefaultGeometry1D : public ModelGenerator
{
public:
  //! \brief The constructor forwards to the base class.
  explicit DefaultGeometry1D(const TiXmlElement* geo) : ModelGenerator(geo) {}
  //! \brief Empty destructor.
  virtual ~DefaultGeometry1D() {}

  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput&) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int nsd, bool rational = false) const;
};


/*!
  \brief Default model generator for 2D FEM simulators.
  \details Generates a quadrilateral domain.
*/

class DefaultGeometry2D : public ModelGenerator
{
public:
  //! \brief The constructor forwards to the base class.
  explicit DefaultGeometry2D(const TiXmlElement* geo) : ModelGenerator(geo) {}
  //! \brief Empty destructor.
  virtual ~DefaultGeometry2D() {}

  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput&) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int nsd, bool rational = false) const;
};


/*!
  \brief Default model generator for 3D FEM simulators.
  \details Generates a hexahedral domain.
*/

class DefaultGeometry3D : public ModelGenerator
{
public:
  //! \brief The constructor forwards to the base class.
  explicit DefaultGeometry3D(const TiXmlElement* geo) : ModelGenerator(geo) {}
  //! \brief Empty destructor.
  virtual ~DefaultGeometry3D() {}

  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput&) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int, bool rational = false) const;
};

#endif
