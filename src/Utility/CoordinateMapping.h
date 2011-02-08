// $Id: CoordinateMapping.h,v 1.2 2010-10-14 20:01:55 kmo Exp $
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

#include "matrix.h"

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
  real Jacobian(matrix<real>& J, matrix<real>& dNdX,
		const matrix<real>& X, const matrix<real>& dNdu,
		bool computeGradient = true);

  //! \brief Set up the Jacobian matrix of the coordinate mapping along an edge.
  //! \param[out] J The inverse of the Jacobian matrix
  //! \param[out] t Unit tangent vector along the edge
  //! \param[out] dNdX 1st order derivatives of basis functions, w.r.t. X
  //! \param[in] X Matrix of element nodal coordinates
  //! \param[in] dNdu First order derivatives of basis functions
  //! \param[in] tangent Parametric tangent direction along the edge
  //! \return The curve dilation of the edge
  real Jacobian(matrix<real>& J, Vec3& t, matrix<real>& dNdX,
		const matrix<real>& X, const matrix<real>& dNdu,
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
  real Jacobian(matrix<real>& J, Vec3& n, matrix<real>& dNdX,
		const matrix<real>& X, const matrix<real>& dNdu,
		size_t t1, size_t t2);

  //! \brief Set up the Hessian matrix of the coordinate mapping.
  //! \param[out] H The Hessian matrix
  //! \param[out] d2NdX2 Second order derivatives of basis functions, w.r.t. X
  //! \param[in] Ji The inverse of the Jacobian matrix
  //! \param[in] X Matrix of element nodal coordinates
  //! \param[in] d2Ndu2 Second order derivatives of basis functions
  //! \param[in] dNdu First order derivatives of basis functions
  //! \param[in] computeGradient If \e false, skip calculation of \a d2NdX2
  //! \return \e false if matrix dimensions are incompatible, otherwise \e true
  bool Hessian(matrix3d<real>& H, matrix3d<real>& d2NdX2,
	       const matrix<real>& Ji, const matrix<real>& X,
	       const matrix3d<real>& d2Ndu2, const matrix<real>& dNdu,
	       bool computeGradient = true);
};

#endif
