// $Id$
//==============================================================================
//!
//! \file ImmersedBoundaries.h
//!
//! \date Dec 18 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Utilities for immersed boundary calculations.
//!
//==============================================================================

#ifndef _IMMERSED_BOUNDARIES_H
#define _IMMERSED_BOUNDARIES_H

#ifndef Real
#define Real double //!< The floating point type to use
#endif

#include <vector>

class Vec3;

namespace utl {
  class Point;
}

class ElementBlock;

//! A real-valued array without algebraic operations
typedef std::vector<Real>      RealArray;
//! A real-valued two-dimensional array without algebraic operations
typedef std::vector<RealArray> Real2DMat;
//! A real-valued three-dimensional array without algebraic operations
typedef std::vector<Real2DMat> Real3DMat;
//! An array of point vectors
typedef std::vector<utl::Point> PointVec;


namespace Immersed //! Utilities for immersed boundary calculations
{
  //! \brief Interface class representing a geometric object.
  class Geometry
  {
  protected:
    //! \brief The default constructor is protected to allow sub-classes only.
    Geometry() {}
  public:
    //! \brief Empty destructor.
    virtual ~Geometry() {}

    //! \brief Performs the inside-outside test for the geometric object.
    //! \details Alpha is used as an indicator here:
    //! Alpha = 0.0 if the point is lying outside the physical domain
    //! Alpha = 0.0 if the point is lying directly on the boundary
    //! Alpha = 1.0 if the point is lying inside the physical domain
    virtual double Alpha(double X, double Y, double Z) const { return 0.0; }
    //! \brief Performs the inside-outside test for the geometric object.
    virtual double Alpha(const Vec3& X) const;

    //! \brief Creates a finite element model of the geometry for visualization.
    virtual ElementBlock* tesselate() const { return nullptr; }
  };

  //! \brief Returns the coordinates and weights for the quadrature points.
  //! \param[in] geo Object describing the boundary of the physical geometry.
  //! The objects returns the inside/outside status of a given spatial points
  //! through its virtual member function \a Alpha.
  //! \param[in] elmCorner Cartesian coordinates of the element corners;
  //! first index is the element counter (0 to number of elements minus 1),
  //! second index is the corner point counter (0 to 3 in 2D, 0 to 7 in 3D),
  //! third index is the coordinate counter (0 to 2)
  //! \param[in] max_depth Maximum depth up to which you want to refine
  //! \param[in] p Order of the Gauss integration
  //! \param[out] quadPoints the quadrature point coordinates and weights;
  //! first index is the element counter (0 to number of elements minus 1),
  //! second index is the quadrature point counter for each element (0 to number
  //! of quadrature points minus 1 for element identified by the first index),
  //! third index is the coordinate/weight index (0=xi, 1=eta, 2=weight in 2D,
  //! 0=xi, 1=eta, 2=zeta, 3=weight in 3D)
  //! \param grid Points to an \a ElementBlock plotting the added grid lines
  //!
  //! \details The element corner points are ordered according to a standard
  //! tensor-product definition of the element, i.e., the index runs fastest
  //! in the first parameter direction, then in the second direction, and
  //! finally (in 3D) the third direction.
  //! The coordinates returned are assumed to be referring to the bi-unit square
  //! (tri-unit cube in 3D) of each element, and the weights are standard Gauss
  //! quadrature weights, which summs to 2 in the power of number of dimensions.
  bool getQuadraturePoints(const Geometry& geo,
                           const std::vector<PointVec>& elmCorner,
                           int max_depth, int p,
                           Real3DMat& quadPoints, ElementBlock* grid = nullptr);

  //! \brief Returns the quadrature points for a 2D element.
  bool getQuadraturePoints(const Geometry& geo,
                           const utl::Point& X1, const utl::Point& X2,
                           const utl::Point& X3, const utl::Point& X4,
                           int max_depth, int nGauss,
                           RealArray& GP1, RealArray& GP2,
                           RealArray& GPw, ElementBlock* grid = nullptr);

  //! \brief Returns the quadrature points for a 3D element.
  bool getQuadraturePoints(const Geometry& geo,
                           const utl::Point& X1, const utl::Point& X2,
                           const utl::Point& X3, const utl::Point& X4,
                           const utl::Point& X5, const utl::Point& X6,
                           const utl::Point& X7, const utl::Point& X8,
                           int max_depth, int nGauss,
                           RealArray& GP1, RealArray& GP2, RealArray& GP3,
                           RealArray& GPw);

  //! \brief Enum defining different stabilizations.
  enum Stab
  {
    NO_STAB = 0,
    ALL_INTERFACES = 1,
    SUBDIV_INTERFACES = 2
  };

  extern int stabilization; //!< Stabilization option

  extern bool plotCells; //!< Flags whether subcells should be plotted or not
}

#endif
