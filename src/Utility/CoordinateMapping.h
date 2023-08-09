// $Id$
//==============================================================================
//!
//! \file CoordinateMapping.h
//!
//! \date May 18 2010
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Utilities for coordinate mapping transformations.
//!
//==============================================================================

#ifndef _COORDINATE_MAPPING_H
#define _COORDINATE_MAPPING_H

#include "matrixnd.h"

class Vec3;


namespace utl
{
  //! \brief Set up the Jacobian matrix of the coordinate mapping.
  //! \param[out] J The inverse of the Jacobian matrix
  //! \param[out] dNdX First order derivatives of basis functions, w.r.t. X
  //! \param[in] X Matrix of element nodal coordinates
  //! \param[in] dNdu First order derivatives of basis functions
  //! \param[in] computeGradient If \e false, skip calculation of \a dNdX
  //! \return The Jacobian determinant
  Real Jacobian(matrix<Real>& J, matrix<Real>& dNdX,
                const matrix<Real>& X, const matrix<Real>& dNdu,
                bool computeGradient = true);

  //! \brief Set up the Jacobian matrix of the coordinate mapping along an edge.
  //! \param[out] J The inverse of the Jacobian matrix
  //! \param[out] t Unit tangent vector along the edge
  //! \param[out] dNdX 1st order derivatives of basis functions, w.r.t. X
  //! \param[in] X Matrix of element nodal coordinates
  //! \param[in] dNdu First order derivatives of basis functions
  //! \param[in] tangent Parametric tangent direction along the edge
  //! \return The curve dilation of the edge
  Real Jacobian(matrix<Real>& J, Vec3& t, matrix<Real>& dNdX,
                const matrix<Real>& X, const matrix<Real>& dNdu,
                size_t tangent);

  //! \brief Set up the Jacobian matrix of the coordinate mapping on a boundary.
  //! \param[out] J The inverse of the Jacobian matrix
  //! \param[out] n Outward-directed unit normal vector on the boundary
  //! \param[out] dNdX 1st order derivatives of basis functions, w.r.t. X
  //! \param[in] X Matrix of element nodal coordinates
  //! \param[in] dNdu First order derivatives of basis functions
  //! \param[in] t1 First parametric tangent direction of the boundary
  //! \param[in] t2 Second parametric tangent direction of the boundary
  //! \return The surface/curve dilation of the boundary
  Real Jacobian(matrix<Real>& J, Vec3& n, matrix<Real>& dNdX,
                const matrix<Real>& X, const matrix<Real>& dNdu,
                size_t t1, size_t t2);

  //! \brief Set up the Hessian matrix of the coordinate mapping.
  //! \param[out] H The Hessian matrix
  //! \param[out] d2NdX2 Second order derivatives of basis functions, w.r.t. X
  //! \param[in] Ji The inverse of the Jacobian matrix
  //! \param[in] X Matrix of element nodal coordinates
  //! \param[in] d2Ndu2 Second order derivatives of basis functions
  //! \param[in] dNdX First order derivatives of basis functions
  //! \param[in] geoMapping If \e true, calculate geometry mapping
  //! \return \e false if matrix dimensions are incompatible, otherwise \e true
  //!
  //! \details If geoMapping is \e true, H is output,
  //! else H is input and assumed to be already calculated in a previous call.
  bool Hessian(matrix3d<Real>& H, matrix3d<Real>& d2NdX2,
               const matrix<Real>& Ji, const matrix<Real>& X,
               const matrix3d<Real>& d2Ndu2, const matrix<Real>& dNdX,
               bool geoMapping = true);
  //! \brief Convert a Hessian from a matrix3d to a matrix assuming symmetry.
  void Hessian(const matrix3d<Real>& Hess, matrix<Real>& H);

  //! \brief Compute the stabilization matrix \b G from the Jacobian inverse.
  //! \param[in] Ji The inverse of the Jacobian matrix
  //! \param[in] du Element lengths in each parametric direction
  //! \param[out] G The stabilization matrix (used in CFD simulators)
  void getGmat(const matrix<Real>& Ji, const Real* du, matrix<Real>& G);

  //! \brief Set up third-order derivatives of the coordinate mapping.
  //! \param[out] d3NdX3 Third order derivatives of basis functions, w.r.t. X
  //! \param[in] Ji The inverse of the Jacobian matrix
  //! \param[in] d3Ndu3 Third order derivatives of basis functions
  //! \return \e false if matrix dimensions are incompatible, otherwise \e true
  bool Hessian2(matrix4d<Real>& d3NdX3,
                const matrix<Real>& Ji, const matrix4d<Real>& d3Ndu3);

  //! \brief Calculates the derivatives of the Jacobian of the coordinate mapping.
  //! \param[in] dudX Derivatives of the geometry basis
  //! \param[in] d2Xdu2 Second order derivatives of the geometry basis
  //! \param[out] dJdX Derivatives of the Jacobian wrt physical coordinates
  void JacobianGradient(const matrix<Real>& dudX,
                        const matrix3d<Real>& d2Xdu2,
                        std::vector<matrix<Real>>& dJdX);

  //! \brief Calculates the derivatives of determinant of the Jacobian.
  //! \param[in] J Jacobian of the geometry mapping
  //! \param[in] Ji Inverse jacobian of the geometry mapping
  //! \param[in] H Hessian of the geometry mapping wrt parameters
  //! \param[out] ddet Derivatives of the determinant
  void detJacGradient(const matrix<Real>& J, const matrix<Real>& Ji,
                      const matrix3d<Real>& H, std::vector<Real>& ddet);
}

#endif
