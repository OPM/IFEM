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


//! \brief Extend basis (raise order and continuity) tuned for subdivided patches.
//! \param[in] old Spline basis to extend
//! \return New spline basis
Go::BsplineBasis extendedBasis(const Go::BsplineBasis& old);


/*!
  \brief 1D multi-patch model generator for FEM simulators.
  \details Generate a line split in a given number of blocks.
*/

class MultiPatchModelGenerator1D : public ModelGenerator
{
public:
  //! \brief The constructor initializes common members.
  //!\ param[in] elem XML element to parse
  MultiPatchModelGenerator1D(const TiXmlElement* elem);
  //! \brief Empty destructor.
  virtual ~MultiPatchModelGenerator1D() {}

  //! \brief Creates a geometry.
  //! \param[in] m Simulator object with patch read function to use
  SIMdependency::PatchVec createGeometry(const SIMinput& m) const;
  //! \brief Creates topology for geometry.
  virtual bool createTopology(SIMinput& sim) const;
  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput& sim) const;

  //! \brief Generates knot vectors for subdivision.
  //! \param[in] cur Univariate patch to extract subpatch from
  //! \param[in] startu Index in source knot vector
  //! \param[in] numcoefsu Number of DOFs (include overlap DOFs)
  //! \param[in] orderu Polynomial order ("p+1")
  static Go::SplineCurve getSubPatch(const Go::SplineCurve* cur,
      const size_t startu, const size_t numcoefsu, const int orderu);

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int nsd) const;

private:
  size_t nx; //!< Number of blocks in x
  int periodic_x; //!< If non-zero, make model periodic in x for given bases
  bool subdivision; //!< Use patch-subdivision and not multiple blocks.
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
  MultiPatchModelGenerator2D(const TiXmlElement* elem);
  //! \brief Empty destructor.
  virtual ~MultiPatchModelGenerator2D() {}

  //! \brief Creates a geometry.
  //! \param[in] m Simulator object with patch read function to use
  SIMdependency::PatchVec createGeometry(const SIMinput& m) const;
  //! \brief Creates topology for geometry.
  virtual bool createTopology(SIMinput& sim) const;
  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput& sim) const;

  //! \brief Generates knot vectors for subdivision.
  //! \param[in] srf Bivariate patch to extract subpatch from
  //! \param[in] startu Index in source u-knot vector
  //! \param[in] numcoefsu Number of DOFs in u-direction (include overlap DOFs)
  //! \param[in] orderu Polynomial order in u-direction ("p+1")
  //! \param[in] startv Index in source v-knot vector
  //! \param[in] numcoefsv Number of DOFs in v-direction
  //! \param[in] orderv Polynomial order in v-direction
  static Go::SplineSurface getSubPatch(const Go::SplineSurface* srf,
      const size_t startu, const size_t numcoefsu, const int orderu,
      const size_t startv, const size_t numcoefsv, const int orderv);
  //! \brief Establishes mixed subdivision bases.
  //! \param sim Simulator object holding MixedType information
  bool establishSubdivisionBases (SIMinput& sim) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int nsd) const;

private:
  size_t nx; //!< Number of blocks in x
  size_t ny; //!< Number of blocks in y
  int periodic_x; //!< If non-zero, make model periodic in x for given bases
  int periodic_y; //!< If non-zero, make model periodic in y for given bases
  bool subdivision; //!< Use patch-subdivision and not multiple blocks.
};


/*!
 \brief 3D multi-patch model generator for FEM simulators.
 \details Generates a hexahedra split in a given number of blocks.
*/

class MultiPatchModelGenerator3D : public ModelGenerator
{
public:
  //! \brief The constructor initializes common members.
  //! \param[in] elem XML element to parse
  MultiPatchModelGenerator3D(const TiXmlElement* geo);
  //! \brief Empty destructor.
  virtual ~MultiPatchModelGenerator3D() {}

  //! \brief Creates a geometry.
  //! \param[in] m Simulator object with patch read function to use
  SIMdependency::PatchVec createGeometry(const SIMinput& m) const;
  //! \brief Creates topology for geometry.
  virtual bool createTopology(SIMinput& sim) const;
  //! \brief Creates topology sets for geometry.
  virtual TopologySet createTopologySets(const SIMinput& sim) const;

  //! \brief Generates knot vectors for subdivision.
  //! \param[in] vol Trivariate patch to extract subpatch from
  //! \param[in] startu Index in source u-knot vector
  //! \param[in] numcoefsu Number of DOFs in u-direction (include overlap DOFs)
  //! \param[in] orderu Polynomial order in u-direction ("p+1")
  //! \param[in] startv Index in source v-knot vector
  //! \param[in] numcoefsv Number of DOFs in v-direction
  //! \param[in] orderv Polynomial order in v-direction
  //! \param[in] startw Index in source w-knot vector
  //! \param[in] numcoefsw Number of DOFs in w-direction
  //! \param[in] orderw Polynomial order in w-direction
  static Go::SplineVolume getSubPatch(const Go::SplineVolume* vol,
      const size_t startu, const size_t numcoefsu, const int orderu,
      const size_t startv, const size_t numcoefsv, const int orderv,
      const size_t startw, const size_t numcoefsw, const int orderw);
  //! \brief Establishes mixed subdivision bases.
  //! \param sim Simulator object holding MixedType information
  bool establishSubdivisionBases (SIMinput& sim) const;

protected:
  //! \brief Generates the G2 description of the geometry.
  virtual std::string createG2(int) const;

private:
  size_t nx; //!< Number of blocks in x
  size_t ny; //!< Number of blocks in y
  size_t nz; //!< Number of blocks in z
  int periodic_x; //!< If non-zero, make model periodic in x for given bases
  int periodic_y; //!< If non-zero, make model periodic in y for given bases
  int periodic_z; //!< If non-zero, make model periodic in z for given bases
  bool subdivision; //!< Use patch-subdivision and not multiple blocks.
};

#endif
