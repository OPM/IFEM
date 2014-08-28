// $Id$
//==============================================================================
//!
//! \file Tensor.C
//!
//! \date Dec 17 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of second-order tensors with some basic operations.
//!
//==============================================================================

#include "Tensor.h"
#include "Vec3.h"
#include <algorithm>

#ifndef epsZ
//! \brief Zero tolerance for the incremental rotations.
#define epsZ 1.0e-16
#endif


std::ostream& Tensor::print (std::ostream& os) const
{
  switch (n) {
  case 1:
    return os << v[0] << std::endl;
  case 2:
    return os << v[0] <<' '<< v[2] <<'\n'
	      << v[1] <<' '<< v[3] << std::endl;
  case 3:
    return os << v[0] <<' '<< v[3] <<' '<< v[6] <<'\n'
	      << v[1] <<' '<< v[4] <<' '<< v[7] <<'\n'
	      << v[2] <<' '<< v[5] <<' '<< v[8] << std::endl;
  }

  return os;
}


Tensor::Tensor (const t_ind nsd, bool identity) : n(nsd)
{
  v.resize(n*n,Real(0));

  if (identity)
    for (t_ind i = 0; i < n*n; i += n+1)
      v[i] = Real(1);
}


/*!
  The provided vector \a vn is taken as the local Z-axis.
  The local X- and Y-axes are then defined by projecting either the global X-
  or Y-axis onto the normal plane, depending on whether \a vn points mostly in
  the global Y- or X-direction, respectively.

  If \a vnIsX is \e true, \a vn is taken as the local X-axis instead and the
  other two axes are shifted accordingly.
*/

Tensor::Tensor (const Vec3& vn, bool vnIsX) : n(3)
{
  Vec3 v1, v2, v3(vn);
  v3.normalize();

  if (vnIsX)
  {
    if (fabs(v3.y) < fabs(v3.z))
    {
      // Define the Y-axis by projecting the global Y-axis
      // onto the normal plane of v3
      v2.x = -v3.y*v3.x;
      v2.y =  v3.x*v3.x + v3.z*v3.z;
      v2.z = -v3.y*v3.z;
      // Define the Z-axis as the cross product of X-axis and Y-axis
      v1.cross(v3,v2);
    }
    else
    {
      // Define the Z-axis by projecting the global Z-axis
      // onto the normal plane of v3
      v1.x = -v3.z*v3.x;
      v1.y = -v3.z*v3.y;
      v1.z =  v3.x*v3.x + v3.y*v3.y;
      // Define the Y-axis as the cross product of Z-axis and X-axis
      v2.cross(v1,v3);
    }
    v1.normalize();
    v2.normalize();
    this->define3Dtransform(v3,v2,v1);
  }
  else
  {
    if (fabs(v3.y) < fabs(v3.x))
    {
      // Define the Y-axis by projecting global Y-axis
      // onto the normal plane of v3
      v2.x = -v3.y*v3.x;
      v2.y =  v3.x*v3.x + v3.z*v3.z;
      v2.z = -v3.y*v3.z;
      // Define the X-axis as the cross product of Y-axis and Z-axis
      v1.cross(v2,v3);
    }
    else
    {
      // Define the X-axis by projecting global X-axis
      // onto the normal plane of v3
      v1.x =  v3.y*v3.y + v3.z*v3.z;
      v1.y = -v3.x*v3.y;
      v1.z = -v3.x*v3.z;
      // Define the Y-axis as the cross product of Z-axis and X-axis
      v2.cross(v3,v1);
    }
    v1.normalize();
    v2.normalize();
    this->define3Dtransform(v1,v2,v3);
  }
}


/*!
  The first tangent vector \a t1 is taken as the local X-axis (or Z-axis, if
  \a t1isZ is \e true). The local Z-axis (X-axis) is then defined as the cross
  product between \a t1 and \a t2. If \a t2isXZ is \e true, the local Y-axis is
  instead defined as the cross product between \a t2 and \a t1.
*/

Tensor::Tensor (const Vec3& t1, const Vec3& t2, bool t1isZ, bool t2isXZ) : n(3)
{
  Vec3 v1(t1), v2, v3;

  v1.normalize();
  if (t2isXZ)
  {
    v2.cross(t2,t1).normalize();
    v3.cross(v1,v2);
  }
  else
  {
    v3.cross(t1,t2).normalize();
    v2.cross(v3,v1);
  }

  if (t1isZ)
    this->define3Dtransform(v2,v3,v1);
  else
    this->define3Dtransform(v1,v2,v3);
}


/*!
  This constructor assumes the three vectors provided form an orthonormal basis.
*/

Tensor::Tensor (const Vec3& v1, const Vec3& v2, const Vec3& v3) : n(3)
{
  this->define3Dtransform(v1,v2,v3);
}


/*!
  This constructor computes an incremental rotation tensor from the given
  rotation angles via a quaternion representation of the rotation.
  The angles provided are those related to a Rodrigues parameterization.
*/

Tensor::Tensor (Real a1, Real a2, Real a3) : n(3)
{
  Vec3 q(a1,a2,a3);
  double theta = q.length();
  q  *= (theta < epsZ ? 0.5 : sin(0.5*theta)/theta);
  Real q0 = cos(0.5*theta);
  Real ql = sqrt(q0*q0 + q.x*q.x + q.y*q.y + q.z*q.z);
  q0 /= ql;
  q  /= ql;

  v.resize(9);
  v[0] = 2.0*(q.x*q.x + q0*q0) - 1.0;
  v[4] = 2.0*(q.y*q.y + q0*q0) - 1.0;
  v[8] = 2.0*(q.z*q.z + q0*q0) - 1.0;

  v[3] = 2.0*(q.x*q.y - q.z*q0);
  v[6] = 2.0*(q.x*q.z + q.y*q0);
  v[7] = 2.0*(q.y*q.z - q.x*q0);

  v[1] = 2.0*(q.y*q.x + q.z*q0);
  v[2] = 2.0*(q.z*q.x - q.y*q0);
  v[5] = 2.0*(q.z*q.y + q.x*q0);

#if SP_DEBUG > 2
  std::cout <<"Tensor("<< a1 <<","<< a2 <<","<< a3 <<"):\n" << *this;
#endif
}


Tensor::Tensor (Real alpha, t_ind axis) : n(3)
{
  v.resize(9,Real(0));
  for (t_ind i = 0; i < 9; i += 4)
    v[i] = Real(1);

  double ca = cos(alpha);
  double sa = sin(alpha);

  switch (axis) {
  case 1:
    v[4] = v[8] = ca;
    v[5] = sa;
    v[7] = -sa;
    break;
  case 2:
    v[0] = v[8] = ca;
    v[2] = -sa;
    v[6] = -sa;
    break;
  case 3:
    v[0] = v[4] = ca;
    v[1] = sa;
    v[3] = -sa;
    break;
  }
}


Tensor::Tensor (const Tensor& T, bool transpose) : n(T.n)
{
  v.resize(n*n);
  std::copy(T.v.begin(),T.v.end(),v.begin());

  if (transpose)
    this->transpose();
}


Tensor::Tensor (const std::vector<Real>& a, bool transpose) : n(sqrt(a.size()))
{
  if (size_t(n*n) < a.size() || n > 3)
    std::cerr <<" *** Tensor(const std::vector<Real>&): Invalid value size "
              << a.size() << std::endl;

  v.resize(n*n);
  this->operator=(a);

  if (transpose)
    this->transpose();
}


void Tensor::define3Dtransform (const Vec3& v1, const Vec3& v2, const Vec3& v3)
{
  v.resize(9);
  for (t_ind i = 0; i < 3; i++)
  {
    v[  i] = v1[i];
    v[3+i] = v2[i];
    v[6+i] = v3[i];
  }
}


Vec3 Tensor::operator[] (t_ind i) const
{
  if (i >= n)
    return Vec3();

  return Vec3(&v[0]+n*i,n);
}


Tensor& Tensor::operator= (const Tensor& T)
{
  if (v.size() == T.v.size())
    std::copy(T.v.begin(),T.v.end(),v.begin());
  else
  {
    // Handle different dimensions and/or symmetry
    t_ind ndim = T.n;
    if (n > ndim)
      this->zero();
    else
      ndim = n;
    for (t_ind i = 1; i <= ndim; i++)
      for (t_ind j = (this->symmetric() ? i : 1); j <= ndim; j++)
	v[this->index(i,j)] = T(i,j);

    if (this->symmetric() && v.size() == 4 && T.v.size() >= 4)
      v[2] = T(3,3);
  }

  return *this;
}


Tensor& Tensor::operator= (const std::vector<Real>& val)
{
  if (val.size() == v.size())
    std::copy(val.begin(),val.end(),v.begin());
  else if (val.size() > v.size())
    std::copy(val.begin(),val.begin()+v.size(),v.begin());
  else
  {
    std::copy(val.begin(),val.end(),v.begin());
    std::fill(v.begin()+val.size(),v.end(),Real(0));
  }

  return *this;
}


Tensor& Tensor::operator= (Real val)
{
  this->zero();

  t_ind i, j;
  for (i = j = 0; i < n; i++, j += n+1)
    v[j] = val;

  return *this;
}


Tensor& Tensor::operator+= (const Tensor& T)
{
  if (T.n > 0)
    if (v.size() == T.v.size())
      for (t_ind i = 0; i < v.size(); i++)
	v[i] += T.v[i];
    else
      std::cerr <<"Tensor::operator+=(const Tensor&): "
		<<"Not implemented for tensors of different size."<< std::endl;

  return *this;
}


Tensor& Tensor::operator+= (Real val)
{
  t_ind i, j, inc = this->symmetric() ? 1 : n+1;
  for (i = j = 0; i < n; i++, j += inc)
    v[j] += val;

  if (inc == 1 && v.size() == 4)
    v[2] += val;

  return *this;
}


Tensor& Tensor::operator-= (const Tensor& T)
{
  if (T.n > 0)
    if (v.size() == T.v.size())
      for (t_ind i = 0; i < v.size(); i++)
	v[i] -= T.v[i];
    else
      std::cerr <<"Tensor::operator-=(const Tensor&): "
		<<"Not implemented for tensors of different size."<< std::endl;

  return *this;
}


Tensor& Tensor::operator-= (Real val)
{
  t_ind i, j, inc = this->symmetric() ? 1 : n+1;
  for (i = j = 0; i < n; i++, j += inc)
    v[j] -= val;

  if (inc == 1 && v.size() == 4)
    v[2] -= val;

  return *this;
}


Tensor& Tensor::operator*= (Real val)
{
  for (t_ind i = 0; i < v.size(); i++)
    v[i] *= val;

  return *this;
}


Tensor& Tensor::postMult (const Tensor& B)
{
  switch (std::min(n,B.n)) {
  case 1:
    v[0] *= B.v[0];
    break;
  case 2:
    {
      Tensor A(*this);
      for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
          v[this->index(i,j)] = A(i,1)*B(1,j) + A(i,2)*B(2,j);
      break;
    }
  case 3:
    {
      Tensor A(*this);
      for (int i = 1; i <= 3; i++)
        for (int j = 1; j <= 3; j++)
          v[this->index(i,j)] = A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j);
      break;
    }
  }

  return *this;
}


Tensor& Tensor::preMult (const Tensor& A)
{
  switch (std::min(n,A.n)) {
  case 1:
    v[0] *= A.v[0];
    break;
  case 2:
    {
      Tensor B(*this);
      for (int i = 1; i <= 2; i++)
        for (int j = 1; j <= 2; j++)
          v[this->index(i,j)] = A(i,1)*B(1,j) + A(i,2)*B(2,j);
      break;
    }
  case 3:
    {
      Tensor B(*this);
      for (int i = 1; i <= 3; i++)
        for (int j = 1; j <= 3; j++)
          v[this->index(i,j)] = A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j);
      break;
    }
  }

  return *this;
}


Real Tensor::innerProd (const Tensor& T) const
{
  Real value = Real(0);
  if (v.size() == T.v.size())
    for (t_ind i = 0; i < v.size(); i++)
      value += v[i]*T.v[i];
  else if (this->symmetric() && T.symmetric())
  {
    // Handle symmetric tensors with different dimensions
    t_ind ndim = n < T.n ? n : T.n;
    for (t_ind i = 1; i <= ndim; i++)
      for (t_ind j = i; j <= ndim; j++)
	value += v[this->index(i,j)]*T(i,j);
  }
  else if (n > 0 && T.n > 0)
    std::cerr <<"Tensor::innerProd(const Tensor&) const: "
	      <<"Not implemented for tensors of different size."<< std::endl;

  return value;
}


bool Tensor::isZero (Real tol) const
{
  for (t_ind i = 0; i < v.size(); i++)
    if (v[i] > tol || -v[i] > tol)
      return false;

  return true;
}


Tensor& Tensor::transpose ()
{
  switch (n) {
  case 2:
    std::swap(v[1],v[2]);
    break;
  case 3:
    std::swap(v[1],v[3]);
    std::swap(v[2],v[6]);
    std::swap(v[5],v[7]);
  }

  return *this;
}


Tensor& Tensor::symmetrize ()
{
  switch (n) {
  case 2:
    v[1] = v[2] = 0.5*(v[1]+v[2]);
    break;
  case 3:
    v[1] = v[3] = 0.5*(v[1]+v[3]);
    v[2] = v[6] = 0.5*(v[2]+v[6]);
    v[5] = v[7] = 0.5*(v[5]+v[7]);
  }

  return *this;
}


Tensor& Tensor::shift (short int idx)
{
  if (this->symmetric() || idx <= -n || (n+idx)%n == 0)
    return *this;

  if (n == 2)
  {
    std::swap(v[0],v[2]);
    std::swap(v[1],v[3]);
  }
  else if (n == 3)
  {
    t_ind j = (3+idx)%3;
    for (t_ind r = 0; r < 3; r++)
    {
      std::swap(v[r],v[r+3*j]);
      std::swap(v[r],v[r-3*j+9]);
    }
  }

  return *this;
}


Real Tensor::trace () const
{
  if (n == 3)
    return v[0] + v[4] + v[8];
  else if (n == 2)
    return v[0] + v[3];
  else if (n == 1)
    return v[0];

  return Real(0);
}


Real Tensor::det () const
{
  if (n == 3)
    return v[0]*(v[4]*v[8] - v[5]*v[7])
      -    v[3]*(v[1]*v[8] - v[2]*v[7])
      +    v[6]*(v[1]*v[5] - v[2]*v[4]);
  else if (n == 2)
    return v[0]*v[3] - v[1]*v[2];
  else if (n == 1)
    return v[0];

  return Real(0);
}


Real Tensor::inverse (Real tol)
{
  Real det = this->det();
  if (det <= tol && det >= -tol)
  {
    std::cerr <<"Tensor::inverse: Singular tensor |T|="<< det << std::endl;
    return Real(0);
  }

  if (n == 3)
  {
    Real T11 = v[0]; Real T12 = v[3]; Real T13 = v[6];
    Real T21 = v[1]; Real T22 = v[4]; Real T23 = v[7];
    Real T31 = v[2]; Real T32 = v[5]; Real T33 = v[8];
    v[0] =  (T22*T33 - T32*T23) / det;
    v[1] = -(T21*T33 - T31*T23) / det;
    v[2] =  (T21*T32 - T31*T22) / det;
    v[3] = -(T12*T33 - T32*T13) / det;
    v[4] =  (T11*T33 - T31*T13) / det;
    v[5] = -(T11*T32 - T31*T12) / det;
    v[6] =  (T12*T23 - T22*T13) / det;
    v[7] = -(T11*T23 - T21*T13) / det;
    v[8] =  (T11*T22 - T21*T12) / det;
  }
  else if (n == 2)
  {
    Real T11 = v[0]; Real T12 = v[2];
    Real T21 = v[1]; Real T22 = v[3];
    v[0] =  T22 / det;
    v[1] = -T21 / det;
    v[2] = -T12 / det;
    v[3] =  T11 / det;
  }
  else if (n == 1)
    v[0] = Real(1) / det;

  return det;
}


/*!
  \brief Multiplication between a Tensor and a point vector.
*/

Vec3 operator* (const Tensor& T, const Vec3& v)
{
  switch (T.n) {
  case 1:
    return Vec3(T.v[0]*v.x, v.y, v.z);
  case 2:
    return Vec3(T(1,1)*v.x + T(1,2)*v.y,
		T(2,1)*v.x + T(2,2)*v.y, v.z);
  case 3:
    return Vec3(T(1,1)*v.x + T(1,2)*v.y + T(1,3)*v.z,
		T(2,1)*v.x + T(2,2)*v.y + T(2,3)*v.z,
		T(3,1)*v.x + T(3,2)*v.y + T(3,3)*v.z);
  }
  return v;
}


/*!
  \brief Multiplication between a point vector and transpose of a Tensor.
*/

Vec3 operator* (const Vec3& v, const Tensor& T)
{
  switch (T.n) {
  case 1:
    return Vec3(T.v[0]*v.x, v.y, v.z);
  case 2:
    return Vec3(T(1,1)*v.x + T(2,1)*v.y,
		T(1,2)*v.x + T(2,2)*v.y, v.z);
  case 3:
    return Vec3(T(1,1)*v.x + T(2,1)*v.y + T(3,1)*v.z,
		T(1,2)*v.x + T(2,2)*v.y + T(3,2)*v.z,
		T(1,3)*v.x + T(2,3)*v.y + T(3,3)*v.z);
  }
  return v;
}


/*!
  \brief Multiplication between two Tensors.
*/

Tensor operator* (const Tensor& A, const Tensor& B)
{
  Tensor C(std::min(A.n,B.n));

  switch (C.n) {
  case 1:
    C.v[0] = A.v[0]*B.v[0];
    break;
  case 2:
    for (int i = 1; i <= 2; i++)
      for (int j = 1; j <= 2; j++)
        C(i,j) = A(i,1)*B(1,j) + A(i,2)*B(2,j);
    break;
  case 3:
    for (int i = 1; i <= 3; i++)
      for (int j = 1; j <= 3; j++)
        C(i,j) = A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j);
    break;
  }

  return C;
}


std::ostream& SymmTensor::print (std::ostream& os) const
{
  switch (n) {
  case 1:
    return os << v.front() << std::endl;
  case 2:
    os << v.front() <<'\n'<< v.back() <<' '<< v[1];
    if (v.size() == 4) os <<"\n0 0 "<< v[2];
    return os << std::endl;
  case 3:
    return os << v[0] <<'\n'
	      << v[3] <<' '<< v[1] <<'\n'
	      << v[5] <<' '<< v[4] <<' '<< v[2] << std::endl;
  }

  return os;
}


SymmTensor::SymmTensor (const t_ind nsd, bool with33) : Tensor(0)
{
  this->redim(nsd,with33);
}


bool SymmTensor::redim (const t_ind nsd, bool with33)
{
  if (n == nsd) return false;

  const_cast<t_ind&>(n) = nsd;
  v.resize(nsd == 2 && with33 ? 4 : n*(n+1)/2, Real(0));
  return true;
}


SymmTensor::SymmTensor (const std::vector<Real>& vec) : Tensor(0)
{
  if (vec.size() > 5)
    this->redim(3);
  else if (vec.size() > 2)
    this->redim(2, vec.size() == 4);
  else if (vec.size() > 0)
    this->redim(1);
  else
    this->redim(0);

  std::copy(vec.begin(),vec.begin()+v.size(),v.begin());
}



SymmTensor& SymmTensor::operator= (Real val)
{
  for (t_ind i = 0; i < v.size(); i++)
    v[i] = i < n || (i == 2 && v.size() == 4) ? val : 0.0;

  return *this;
}


void SymmTensor::copy (const SymmTensor& T)
{
  this->redim(T.n, T.v.size() == 4);
  std::copy(T.v.begin(),T.v.end(),v.begin());
}


/*!
  Perform the triple matrix product \f[{\bf A} = {\bf T A T}^T\f]
  where \b A = \a *this
*/

SymmTensor& SymmTensor::transform (const Tensor& T)
{
  if (T.symmetric()) return *this;

  Real S11, S12, S13, S21, S22, S23, S31, S32, S33;
  switch (n) {
  case 2:
    if (T.dim() > 1)
    {
      S11 = T(1,1)*v.front() + T(1,2)*v.back();
      S12 = T(1,1)*v.back()  + T(1,2)*v[1];
      S21 = T(2,1)*v.front() + T(2,2)*v.back();
      S22 = T(2,1)*v.back()  + T(2,2)*v[1];

      v.front() = S11*T(1,1) + S12*T(1,2);
      v[1]      = S21*T(2,1) + S22*T(2,2);
      v.back()  = S11*T(2,1) + S12*T(2,2);
    }
    break;

  case 3:
    if (T.dim() > 2)
    {
      S11 = T(1,1)*v[0] + T(1,2)*v[3] + T(1,3)*v[5];
      S12 = T(1,1)*v[3] + T(1,2)*v[1] + T(1,3)*v[4];
      S13 = T(1,1)*v[5] + T(1,2)*v[4] + T(1,3)*v[2];
      S21 = T(2,1)*v[0] + T(2,2)*v[3] + T(2,3)*v[5];
      S22 = T(2,1)*v[3] + T(2,2)*v[1] + T(2,3)*v[4];
      S23 = T(2,1)*v[5] + T(2,2)*v[4] + T(2,3)*v[2];
      S31 = T(3,1)*v[0] + T(3,2)*v[3] + T(3,3)*v[5];
      S32 = T(3,1)*v[3] + T(3,2)*v[1] + T(3,3)*v[4];
      S33 = T(3,1)*v[5] + T(3,2)*v[4] + T(3,3)*v[2];

      v[0] = S11*T(1,1) + S12*T(1,2) + S13*T(1,3);
      v[1] = S21*T(2,1) + S22*T(2,2) + S23*T(2,3);
      v[2] = S31*T(3,1) + S32*T(3,2) + S33*T(3,3);

      v[3] = S11*T(2,1) + S12*T(2,2) + S13*T(2,3);
      v[4] = S21*T(3,1) + S22*T(3,2) + S23*T(3,3);
      v[5] = S31*T(1,1) + S32*T(1,2) + S33*T(1,3);
    }
    else if (T.dim() == 2)
    {
      S11 = T(1,1)*v[0] + T(1,2)*v[3];
      S12 = T(1,1)*v[3] + T(1,2)*v[1];
      S13 = T(1,1)*v[5] + T(1,2)*v[4];
      S21 = T(2,1)*v[0] + T(2,2)*v[3];
      S22 = T(2,1)*v[3] + T(2,2)*v[1];
      S23 = T(2,1)*v[5] + T(2,2)*v[4];
      S31 = v[5];
      S32 = v[4];
      S33 = v[2];

      v[0] = S11*T(1,1) + S12*T(1,2);
      v[1] = S21*T(2,1) + S22*T(2,2);
      v[3] = S11*T(2,1) + S12*T(2,2);
      v[4] = S23;
      v[5] = S31*T(1,1) + S32*T(1,2);
    }
  }

  return *this;
}


Real SymmTensor::trace () const
{
  Real t = Real(0);

  if (n == 3 || v.size() == 4)
    t = v[0] + v[1] + v[2];
  else if (n == 2)
    t = v[0] + v[1];
  else if (n == 1)
    t = v[0];

  return t;
}


Real SymmTensor::det () const
{
  Real d = Real(0);

  if (n == 3)
    d = v[0]*(v[1]*v[2] - v[4]*v[4])
      - v[3]*(v[3]*v[2] - v[5]*v[4])
      + v[5]*(v[3]*v[4] - v[5]*v[1]);
  else if (n == 2)
    d = (v.front()-v.back())*v[1];
  else if (n == 1)
    d = v.front();

  if (v.size() == 4)
    d *= v[2];

  return d;
}


Real SymmTensor::inverse (Real tol)
{
  Real det = this->det();
  if (det <= tol && det >= -tol)
  {
    std::cerr <<"SymmTensor::inverse: Singular tensor |T|="<< det << std::endl;
    return Real(0);
  }

  if (n == 3)
  {
    Real T11 = v[0];
    Real T21 = v[3]; Real T22 = v[1];
    Real T31 = v[5]; Real T32 = v[4]; Real T33 = v[2];
    v[0] =  (T22*T33 - T32*T32) / det;
    v[1] =  (T11*T33 - T31*T31) / det;
    v[2] =  (T11*T22 - T21*T21) / det;
    v[3] = -(T21*T33 - T31*T32) / det;
    v[4] = -(T11*T32 - T31*T21) / det;
    v[5] =  (T21*T32 - T31*T22) / det;
  }
  else if (n == 2)
  {
    Real T11 = v.front();
    Real T21 = v.back(); Real T22 = v[1];
    v.front()=  T22 / det;
    v[1]     =  T11 / det;
    v.back() = -T21 / det;
    if (v.size() == 4)
    {
      v[0] *= v[2];
      v[1] *= v[2];
      v[3] *= v[2];
      v[2] = (T11*T22 - T21*T21) / det;
    }
  }
  else if (n == 1)
    v.front() = Real(1) / det;

  return det;
}


SymmTensor& SymmTensor::rightCauchyGreen (const Tensor& F)
{
  this->redim(F.dim());

  switch (n) {
  case 1:
    v[0] = F(1,1)*F(1,1);
    break;

  case 2:
    v[0] = F(1,1)*F(1,1) + F(2,1)*F(2,1);
    v[1] = F(1,2)*F(1,2) + F(2,2)*F(2,2);
    v[2] = F(1,1)*F(1,2) + F(2,1)*F(2,2);
    break;

  case 3:
    v[0] = F(1,1)*F(1,1) + F(2,1)*F(2,1) + F(3,1)*F(3,1);
    v[1] = F(1,2)*F(1,2) + F(2,2)*F(2,2) + F(3,2)*F(3,2);
    v[2] = F(1,3)*F(1,3) + F(2,3)*F(2,3) + F(3,3)*F(3,3);
    v[3] = F(1,1)*F(1,2) + F(2,1)*F(2,2) + F(3,1)*F(3,2);
    v[4] = F(1,2)*F(1,3) + F(2,2)*F(2,3) + F(3,2)*F(3,3);
    v[5] = F(1,1)*F(1,3) + F(2,1)*F(2,3) + F(3,1)*F(3,3);
  }

  return *this;
}


Real SymmTensor::L2norm (bool doSqrt) const
{
  double l2n = 0.0;
  for (t_ind i = 0; i < v.size(); i++)
    l2n += (i < n || (i == 2 && v.size() == 4) ? 1.0 : 2.0)*v[i]*v[i];

  return doSqrt ? sqrt(l2n) : l2n;
}


/*!
  The von Mises value of a symmetric 3D (stress) tensor is defined as follows:
  \f[ s_{\rm vm} = \sqrt{
  s_{11}(s_{11}-s_{22}) +
  s_{22}(s_{22}-s_{33}) +
  s_{33}(s_{33}-s_{11}) + 3(s_{12}^2 + s_{23}^2 + s_{31}^2)}
  \f]
*/

Real SymmTensor::vonMises (bool doSqrt) const
{
  double vms = 0.0;
  if (n == 3)
  {
    double s11 = v[0];
    double s22 = v[1];
    double s33 = v[2];
    double s12 = v[3];
    double s23 = v[4];
    double s31 = v[5];
    vms = s11*(s11-s22) + s22*(s22-s33) + s33*(s33-s11) +
          3.0*(s12*s12 + s23*s23 + s31*s31);
  }
  else if (n == 2)
  {
    double s11 = v[0];
    double s22 = v[1];
    double s33 = v.size() == 4 ? v[2] : 0.0;
    double s12 = v.size() == 4 ? v[3] : v[2];
    vms = s11*(s11-s22) + s22*(s22-s33) + s33*(s33-s11) + 3.0*s12*s12;
  }
  else if (n == 1)
    return doSqrt ? v.front() : v.front()*v.front();

  return doSqrt ? sqrt(vms) : vms;
}


void SymmTensor::principal (Vec3& p) const
{
  if (n < 2) return;

  // Compute mean and deviatoric (upper triangular part) tensors
  const Real tol(1.0e-12);

  Real b1 = this->trace() / Real(3);
  Real s1 = v[0] - b1;
  Real s2 = v[1] - b1;
  Real s3 = v.size() > 3 ? v[2] - b1 : Real(0);
  Real s4 = n > 2 ? v[3] : v.back();
  Real s5 = n > 2 ? v[4] : Real(0);
  Real s6 = n > 2 ? v[5] : Real(0);

  // Compute 2nd and 3rd invariants of deviator J_2 and J_3

  Real c1 = s4*s4;
  Real c2 = s5*s5;
  Real c3 = s6*s6;
  Real b2 = (s1*s1 + s2*s2 + s3*s3)/Real(2) + c1 + c2 + c3;
  if (b2 <= tol*b1*b1)
  {
    p = s1;
    return;
  }

  Real b3 = s1*s2*s3 + (s4+s4)*s5*s6 + s1*(c1-c2) + s2*(c1-c3);

  // Set constants

  c1 = Real(2)*sqrt(b2/Real(3));
  c2 = Real(4)*b3;
  c3 = c1*c1*c1;
  Real al = atan2(sqrt(fabs(c3*c3-c2*c2)),c2)/Real(3);
  Real pi23 = M_PI*Real(2)/Real(3);

  // Set principal values

  p.x = b1 + c1*cos(al);
  p.y = b1 + c1*cos(al-pi23);
  p.z = b1 + c1*cos(al+pi23);

  // Ensure that p.x >= p.y >= p.z

  if (p.x < p.y) std::swap(p.x,p.y);
  if (p.y < p.z) std::swap(p.y,p.z);
  if (p.x < p.y) std::swap(p.x,p.y);
}


/*!
  \brief Adding a scaled unit tensor to a symmetric tensor.
*/

SymmTensor operator+ (const SymmTensor& T, Real a)
{
  SymmTensor S(T);

  for (unsigned short int i = 0; i < S.n; i++)
    S.v[i] += a;

  if (S.v.size() == 4)
    S.v[2] += a;

  return S;
}


/*!
  \brief Subtracting a scaled unit tensor from a symmetric tensor.
*/

SymmTensor operator- (const SymmTensor& T, Real a)
{
  SymmTensor S(T);

  for (unsigned short int i = 0; i < S.n; i++)
    S.v[i] -= a;

  if (S.v.size() == 4)
    S.v[2] -= a;

  return S;
}


/*!
  \brief Multiplication between a scalar and a symmetric tensor.
*/

SymmTensor operator* (Real a, const SymmTensor& T)
{
  SymmTensor S(T.dim(), T.v.size() == 4);

  for (size_t i = 0; i < T.v.size(); i++)
    S.v[i] = a*T.v[i];

  return S;
}


SymmTensor4::SymmTensor4 (const std::vector<Real>& x, t_ind nsd)
  : n(nsd), m(0), v(x)
{
  if (n == 3)
    m = 6;
  else if (n == 2)
    m = 3;
  else
    std::cerr <<" *** Invalid fourth-order tensor, dim="<< n << std::endl;

  if (v.size() < (size_t)m*m)
    std::cerr <<" *** Invalid fourth-order tensor,"
	      <<" matrix represention too small, size="<< v.size() << std::endl;

  ptr = (Real*)&v.front();
}


const Real& SymmTensor4::operator() (t_ind i, t_ind j, t_ind k, t_ind l) const
{
  return v[this->index(i,j)*m+this->index(k,l)];
}


Real& SymmTensor4::operator() (t_ind i, t_ind j, t_ind k, t_ind l)
{
  return ptr[this->index(i,j)*m+this->index(k,l)];
}


int LocalSystem::patch = 0;
