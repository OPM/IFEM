// $Id$
//==============================================================================
//!
//! \file SplineInterpolator.h
//!
//! \date Jan 2012
//!
//! \author Anette Stahl
//!
//! \brief Implementation of interpolation/projection schemes for B-splines.
//!
//==============================================================================

#ifndef _SPLINE_INTERPOLATOR_H
#define _SPLINE_INTERPOLATOR_H

#include <vector>

namespace Go {
  class BsplineBasis;
};


namespace SplineInterpolator //! Spline interpolation/projection schemes.
{
  //! \brief Global spline interpolation method.
  //! \param[in] params Parameters for the data points.
  //! Its size is equal to the total number of data points.
  //! \param[in] points The data points.
  //! Its size is equal to the total number of data points multiplied with
  //! the spatial dimension. (NB: The algorithm deduces the spatial dimension to
  //! use by dividing the size of \a points with the size of \a params).
  //! \param[in] tangent_points Tangents associated with those data points that
  //! are referred to in \a params.
  //! Its size is equal to the total number of datapoints that has tangents,
  //! multiplied by the dimension of the space.
  //! \param[in] basis The basis of the B-spline
  //! \param[out] coefs Control point values of the generated spline curve.
  //! (Its basis can be obtained by calling the basis() method).
  void interpolate(const std::vector<double>& params,
		   const std::vector<double>& points,
		   const std::vector<double>& tangent_points,
		   const Go::BsplineBasis& basis,
		   std::vector<double>& coefs);

  //! \brief Local spline approximation method (Quasi-Interpolation).
  //! \param[in] params Parameters for the data points.
  //! Its size is equal to the local number of data points.
  //! \param[in] points The local data points
  //! \param[in] tangent_points Tangents associated with those data points that
  //! are referred to in \a params
  //! \param[in] index Index of control point to extract in the sliding window
  //! \param[in] basis The basis of the B-spline
  //! \param[out] coefs Local control point values of the generated spline curve
  void quasiinterpolate(const std::vector<double>& params,
			const std::vector<double>& points,
			const std::vector<double>& tangent_points,
			const Go::BsplineBasis& basis, int index,
			std::vector<double>& coefs);

  //! \brief Global spline approximation method (Least-Square Fit).
  //! \param[in] params Parameters for the data points
  //! \param[in] paramsweights Parameter Gauss weights
  //! \param[in] points Data points evaluated at the Gauss points
  //! \param[in] tangent_points Tangents associated with those data points that
  //! are referred to in \a params
  //! \param[in] basis The basis of the B-spline
  //! \param[out] coefs Control point values of the generated spline curve
  void leastsquare_approximation(const std::vector<double>& params,
				 const std::vector<double>& paramsweights,
				 const std::vector<double>& points,
				 const std::vector<double>& tangent_points,
				 const Go::BsplineBasis& basis,
				 std::vector<double>& coefs);
};

#endif
