#pragma once
//==============================================================================
//!
//! \file SplineInterpolator.h
//!
//! \date Jan 2012
//!
//! \author Anette Stahl
//!
//! \brief Implementation of interpolation/projection schemes for B-splines
//!
//==============================================================================

#include <vector>

namespace Go {
  class BsplineBasis;
};


//! \brief Helper class for different spline interpolation/projection schemes
class SplineInterpolator {
  public:
    //! \brief Global spline interpolation Method
    //! \param[in] params Vector containing the parameters for the data points. Its size is equal to the total number of datapoints.
    //! \param[in] points Vector containing the data points. Its size is equal to the total number of datapoints multiplied with the spatial dimension. (NB: This is how the algorithm deduces the spatial dimension to use: divide the size of 'points' with the size of 'params').
    //! \param[in] tangent_points A vector containing the tangents associated with those data points that are referred to in 'param_index'. Its size is equal to the total number of datapoints that has tangents, multiplied by the dimension of the space.
    //! \param[in] basis The basis of the B-spline
    //! \param[out] coefs Upon function completion, this vector will hold the coordinates of the control points of the generated spline curve. (Its basis can be obtained by calling the basis() function). 
    static void interpolate(const std::vector<double>& params,
                            const std::vector<double>& points,
                            const std::vector<double>& tangent_points,
                            const Go::BsplineBasis& basis,
                            std::vector<double>& coefs);

    //! \brief Local spline approximation method (Quasi-Interpolation). 
    //! \param[in] params Vector containing the parameters for the data points. Its size is equal to the local number of datapoints.
    //! \param[in] points Vector containing the local data points.
    //! \param[in] tangent_points A vector containing the tangents associated with those data points that are referred to in 'param_index'. 
    //! \param[in] index The index of the control point to extract in the sliding window
    //! \param[in] basis The basis of the B-spline
    //! \param[out] coefs Upon function completion, this vector will hold the coordinates of the local control points of the generated spline curve. 
    static void quasiinterpolate(const std::vector<double>& params,
                                 const std::vector<double>& points,
                                 const std::vector<double>& tangent_points,
                                 const Go::BsplineBasis& basis,
                                 int index, std::vector<double>& coefs);

    //! \brief Global spline approximation method (Least-Square Fit).
    //! \param[in] params Vector containing the parameters for the data points.
    //! \param[in] paramsweights Vector containing the parameter Gauss weights, respectively.
    //! \param[in] points Vector containing the data points evaluated at the Gauss points.
    //! \param[in] tangent_points A vector containing the tangents associated with those data points that are referred to in 'param_index'. 
    //! \param[in] basis The basis of the B-spline
    //! \param[out] coefs Upon function completion, this vector will hold the coordinates of the control points of the generated spline curve. 
    static void leastsquare_approximation(const std::vector<double>& params,
                                          const std::vector<double>& paramsweights,
                                          const std::vector<double>& points,
                                          const std::vector<double>& tangent_points,
                                          const Go::BsplineBasis& basis,
                                          std::vector<double>& coefs);
};
