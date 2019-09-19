// $Id$
//==============================================================================
//!
//! \file IBGeometries.h
//!
//! \date Jan 24 2013
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Physical geometries for immersed boundary simulations.
//!
//==============================================================================

#ifndef _IB_GEOMETRIES_H
#define _IB_GEOMETRIES_H

#include "ImmersedBoundaries.h"

class RealFunc;


/*!
  \brief Class representing the perforated plate benchmark.

  \details See the example Fig. 38 in the recent immersed boundary paper
  coming out of Tom's group).

  The following global information is needed:
  (a) The center of the hole; (b) the radius of the hole denoted by variable R.
  The function Alpha receives the global coordinates X and Y of the point
  under consideration.
  (may be a vertex of an integration element during the set up of the adaptive
  integration structure,
  or an intergration point during the integration of the stiffness matrix).
*/

class Hole2D : public Immersed::Geometry
{
public:
  //! \brief Default constructor initializing the radius and center of the hole.
  Hole2D(double r = 1.0, double x = 0.0, double y = 0.0) : R(r), Xc(x), Yc(y) {}
  //! \brief Empty destructor.
  virtual ~Hole2D() {}

  //! \brief Performs the inside-outside test for the perforated plate object.
  //! \details Alpha is used as an indicator here:
  //! Alpha = 0.0 if the point is lying outside the physical domain
  //! Alpha = 0.0 if the point is lying directly on the boundary
  //! Alpha = 1.0 if the point is lying inside the physical domain
  virtual double Alpha(double X, double Y, double = 0.0) const;

  //! \brief Creates a finite element model of the geometry for visualization.
  virtual ElementBlock* tesselate() const;

protected:
  double R;  //!< Hole radius
  double Xc; //!< X-coordinate of hole center
  double Yc; //!< Y-coordinate of hole center
};


/*!
  \brief Class representing a plate perforated by an oval hole.
*/

class Oval2D : public Hole2D
{
public:
  //! \brief Default constructor initializing the radius and center of the hole.
  Oval2D(double r, double x0, double y0, double x1, double y1);
  //! \brief Empty destructor.
  virtual ~Oval2D() {}

  //! \brief Performs the inside-outside test for the perforated plate object.
  virtual double Alpha(double X, double Y, double) const;

  //! \brief Creates a finite element model of the geometry for visualization.
  virtual ElementBlock* tesselate() const;

private:
  double X1; //!< X-coordinate of second circle center
  double Y1; //!< Y-coordinate of second circle center
  double D2; //!< Auxilliary parameter used by method Alpha
  double LR; //!< Auxilliary parameter used by method Alpha
};


/*!
  \brief Class representing a plate perforated by multiple holes.
*/

class PerforatedPlate2D : public Immersed::Geometry
{
public:
  //! \brief Default constructor.
  PerforatedPlate2D() {}
  //! \brief Constructor creating a single hole.
  explicit PerforatedPlate2D(Hole2D* hole) : holes({hole}) { }
  //! \brief The destructor deletes the holes.
  virtual ~PerforatedPlate2D();

  //! \brief Adds a hole to the perforated plate.
  void addHole(double r, double x, double y);
  //! \brief Adds a hole to the perforated plate.
  void addHole(double r, double x0, double y0, double x1, double y1);

  //! \brief Performs the inside-outside test for the perforated plate object.
  virtual double Alpha(double X, double Y, double) const;

  //! \brief Creates a finite element model of the geometry for visualization.
  virtual ElementBlock* tesselate() const;

private:
  std::vector<Hole2D*> holes; //!< The holes that perforate the plate
};


/*!
  \brief Class describing the immersed boundary as a scalar function.
*/

class GeoFunc2D : public Immersed::Geometry
{
public:
  //! \brief Default constructor.
  explicit GeoFunc2D(RealFunc* f = nullptr, double p = 1.0, double eps = 0.5);
  //! \brief The destructor deletes the function.
  virtual ~GeoFunc2D();

  //! \brief Performs the inside-outside test for the geometric object.
  virtual double Alpha(const Vec3& X) const;

private:
  RealFunc* myAlpha; //!< Function describing the inside-outside state
  double myExponent; //!< The exponent to apply on the \a myAlpha function
  double  threshold; //!< Inside/outside threshold
};

#endif
