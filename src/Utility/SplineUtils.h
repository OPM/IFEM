// $Id$
//==============================================================================
//!
//! \file SplineUtils.h
//!
//! \date Mar 29 2012
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Various utility functions on spline objects - GoTools extensions.
//!
//==============================================================================

#ifndef _SPLINE_UTILS_H
#define _SPLINE_UTILS_H

#include "MatVec.h"

class FunctionBase;
class Vec4;
class Vec3;

namespace Go {
  class Point;
  struct BasisDerivsSf;
  struct BasisDerivsSf2;
  struct BasisDerivsSf3;
#if GoTools_VERSION_MAJOR > 4
  struct BasisDerivsVol;
  struct BasisDerivsVol2;
  struct BasisPtsVol;
#else
  struct BasisDerivs;
  struct BasisDerivs2;
  struct BasisPts;
#endif
  class BsplineBasis;
  class SplineCurve;
  class SplineSurface;
  class SplineVolume;
}

#if GoTools_VERSION_MAJOR > 4
using GoBasisPtsVol = Go::BasisPtsVol; //!< Forward for new version
using GoBasisDerivsVol = Go::BasisDerivsVol; //!< Forward for new version
using GoBasisDerivsVol2 = Go::BasisDerivsVol2; //!< Forward for new version
#else
using GoBasisPtsVol = Go::BasisPts; //!< Type alias for old version
using GoBasisDerivsVol = Go::BasisDerivs; //!< Type alias for old version
using GoBasisDerivsVol2 = Go::BasisDerivs2; //!< Type alias for old version
#endif


namespace SplineUtils //! Various utility functions on spline objects.
{
  //! \brief Helper method for casting a \a Go::Point object to Vec3.
  Vec3 toVec3(const Go::Point& X, int nsd = 3);
  //! \brief Helper method for casting a \a Go::Point and time object to Vec4.
  Vec4 toVec4(const Go::Point& X, Real time = Real(0),
              const double* u = nullptr);

  //! \brief Evaluates given spline curve at a parametric point.
  void point(Vec3& X, double u, const Go::SplineCurve* curve);
  //! \brief Evaluates given spline surface at a parametric point.
  void point(Vec3& X, double u, double v, const Go::SplineSurface* surf);
  //! \brief Evaluates given spline colume at a parametric point.
  void point(Vec3& X, double u, double v, double w, const Go::SplineVolume* vol);

  //! \brief Establishes matrices with basis functions and 1st derivatives.
  void extractBasis(const Go::BasisDerivsSf& spline,
                    Vector& N, Matrix& dNdu);
  //! \brief Establishes matrices with basis functions, 1st and 2nd derivatives.
  void extractBasis(const Go::BasisDerivsSf2& spline,
                    Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2);
  //! \brief Establishes matrices with basis functions, 1st, 2nd
  //! and 3rd derivatives.
  void extractBasis(const Go::BasisDerivsSf3& spline,
                    Vector& N, Matrix& dNdu,
                    Matrix3D& d2Ndu2, Matrix4D& d3Ndu3);

  //! \brief Establishes matrices with basis functions and 1st derivatives.
  void extractBasis(const GoBasisDerivsVol& spline,
                    Vector& N, Matrix& dNdu);
  //! \brief Establishes matrices with basis functions, 1st and 2nd derivatives.
  void extractBasis(const GoBasisDerivsVol2& spline,
                    Vector& N, Matrix& dNdu, Matrix3D& d2Ndu2);

  //! \brief Extracts parameter values of the Gauss points for a spline basis.
  //! \param[out] uGP Parameter values for all points
  //! \param[in] nGP Number of Gauss points along a knot-span
  //! \param[in] xi Dimensionless Gauss point coordinates [-1,1]
  //! \param[in] basis Spline basis to calculate Gauss point parameters for
  //! \param[in] skipNullSpans If \e true, consider non-zero knot spans only
  void getGaussParameters(RealArray& uGP, int nGP, const double* xi,
                          const Go::BsplineBasis& basis,
                          bool skipNullSpans = false);

  //! \brief Projects a spatial function onto a spline curve.
  Go::SplineCurve* project(const Go::SplineCurve* curve,
                           const FunctionBase& f,
                           int nComp = 1, Real time = Real(0));

  //! \brief Projects a spatial function onto a spline surface.
  Go::SplineSurface* project(const Go::SplineSurface* surface,
                             const FunctionBase& f,
                             int nComp = 1, Real time = Real(0));

  //! \brief Projects a spatial function onto a spline volume.
  Go::SplineVolume* project(const Go::SplineVolume* volume,
                            const FunctionBase& f,
                            int nComp = 1, Real time = Real(0));

  //! \brief Enumeration of basis adjustment operations.
  enum class AdjustOp {
    Lower,    //!< Lower order by one
    Original, //!< Keep original
    Raise,    //!< Raise order by one
  };

  //! \brief Returns a basis adjusted according to the given operation.
  Go::BsplineBasis adjustBasis(const Go::BsplineBasis& basis,
                               AdjustOp adjust);

  //! \brief Builds a knot vector from a given polynomial order, knots and continuities.
  std::vector<double> buildKnotVector(int p, const std::vector<double>& simple_knots,
                                      const std::vector<int>& continuities);
}

#endif
