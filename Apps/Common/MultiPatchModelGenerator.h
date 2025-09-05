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

#ifndef _MULTI_PATCH_MODEL_GENERATOR_H
#define _MULTI_PATCH_MODEL_GENERATOR_H

#include "ModelGenerator.h"
#include <GoTools/geometry/SplineCurve.h>
#include <GoTools/geometry/SplineSurface.h>
#include <GoTools/trivariate/SplineVolume.h>


/*!
  \brief 1D multi-patch model generator for FEM simulators.
  \details Generate a line split in a given number of blocks.
*/

class MultiPatchModelGenerator1D : public ModelGenerator
{
public:
  //! \brief The constructor initializes common members.
  //!\ param[in] elem XML element to parse
  explicit MultiPatchModelGenerator1D(const tinyxml2::XMLElement* elem);
  //! \brief Empty destructor.
  virtual ~MultiPatchModelGenerator1D() {}

  //! \brief Creates a geometry.
  virtual bool createGeometry(SIMinput& sim) const;
  //! \brief Creates topology for geometry.
  virtual bool createTopology(SIMinput& sim) const;
  //! \brief Creates topology sets for geometry.
  virtual bool createTopologySets(SIMinput& sim) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int nsd, bool rational = false) const;

private:
  size_t nx; //!< Number of blocks in x
  int periodic_x; //!< If non-zero, make model periodic in x for given bases
};


/*!
  \brief 2D multi-patch model generator for FEM simulators.
  \details Generate a rectangle split in a given number of blocks.
*/

class MultiPatchModelGenerator2D : public ModelGenerator
{
public:
  //! \brief The constructor initializes common members.
  //!\ param[in] elem XML element to parse
  explicit MultiPatchModelGenerator2D(const tinyxml2::XMLElement* elem);
  //! \brief Empty destructor.
  virtual ~MultiPatchModelGenerator2D() {}

  //! \brief Creates a geometry.
  virtual bool createGeometry(SIMinput& sim) const;
  //! \brief Creates topology for geometry.
  virtual bool createTopology(SIMinput& sim) const;
  //! \brief Creates topology sets for geometry.
  virtual bool createTopologySets(SIMinput& sim) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int nsd, bool rational = false) const;

private:
  size_t nx; //!< Number of blocks in x
  size_t ny; //!< Number of blocks in y
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
  //! \brief The constructor initializes common members.
  //! \param[in] geo XML element to parse
  explicit MultiPatchModelGenerator3D(const tinyxml2::XMLElement* geo);
  //! \brief Empty destructor.
  virtual ~MultiPatchModelGenerator3D() {}

  //! \brief Creates a geometry.
  virtual bool createGeometry(SIMinput& sim) const;
  //! \brief Creates topology for geometry.
  virtual bool createTopology(SIMinput& sim) const;
  //! \brief Creates topology sets for geometry.
  virtual bool createTopologySets(SIMinput& sim) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int, bool rational = false) const;

private:
  size_t nx; //!< Number of blocks in x
  size_t ny; //!< Number of blocks in y
  size_t nz; //!< Number of blocks in z
  int periodic_x; //!< If non-zero, make model periodic in x for given bases
  int periodic_y; //!< If non-zero, make model periodic in y for given bases
  int periodic_z; //!< If non-zero, make model periodic in z for given bases
};

#endif
