// $Id$
//==============================================================================
//!
//! \file Vec3Oper.h
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Global algebraic operators involving Vec3 objects.
//!
//==============================================================================

#ifndef _VEC3_OPER_H
#define _VEC3_OPER_H

#include "matrix.h"

class Vec3;


//! \brief Multiplication of a matrix and a vector.
Vec3 operator*(const utl::matrix<real>& A, const std::vector<real>& x);

//! \brief Multiplication of a vector and a scalar.
Vec3 operator*(const Vec3& a, real value);

//! \brief Multiplication of a scalar and a vector.
Vec3 operator*(real value, const Vec3& a);

//! \brief Division of a vector by a scalar,
Vec3 operator/(const Vec3& a, real value);

//! \brief Dot product between two vectors.
real operator*(const Vec3& a, const Vec3& b);

//! \brief Summation of two vectors.
Vec3 operator+(const Vec3& a, const Vec3& b);

//! \brief Subraction of two vectors.
Vec3 operator-(const Vec3& a, const Vec3& b);

//! \brief Equality operator.
bool operator==(const Vec3& a, const Vec3& b);

//! \brief Unequality operator.
bool operator!=(const Vec3& a, const Vec3& b);

//! \brief Less-than operator.
bool operator<(const Vec3& a, const Vec3& b);

//! \brief Output stream operator.
std::ostream& operator<<(std::ostream& os, const Vec3& a);

#endif
