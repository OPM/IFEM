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


Tensor::Tensor (const std::vector<real>& t1, const std::vector<real>& t2) : n(3)
{
  Vec3 v1(t1), v2(t2), v3;

  v1.normalize();
  v3.cross(v1,v2).normalize();
  v2.cross(v3,v1);

  v.resize(9);
  for (t_ind i = 0; i < 3; i++)
  {
    v[  i] = v1[i];
    v[3+i] = v2[i];
    v[6+i] = v3[i];
  }
}


Tensor::Tensor (const Tensor& T) : n(T.n)
{
  v.resize(n*n);
  this->operator=(T);
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


Tensor& Tensor::operator= (const std::vector<real>& val)
{
  if (val.size() == v.size())
    std::copy(val.begin(),val.end(),v.begin());
  else if (val.size() > v.size())
    std::copy(val.begin(),val.begin()+v.size(),v.begin());
  else
  {
    std::copy(val.begin(),val.end(),v.begin());
    std::fill(v.begin()+val.size(),v.end(),real(0));
  }

  return *this;
}


Tensor& Tensor::operator= (real val)
{
  this->zero();

  t_ind i, j, inc = this->symmetric() ? 1 : n+1;
  for (i = j = 0; i < n; i++, j += inc)
    v[j] = val;

  if (inc == 1 && v.size() == 4)
    v[2] = val;

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


Tensor& Tensor::operator+= (real val)
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


Tensor& Tensor::operator-= (real val)
{
  t_ind i, j, inc = this->symmetric() ? 1 : n+1;
  for (i = j = 0; i < n; i++, j += inc)
    v[j] -= val;

  if (inc == 1 && v.size() == 4)
    v[2] -= val;

  return *this;
}


Tensor& Tensor::operator*= (real val)
{
  for (t_ind i = 0; i < v.size(); i++)
    v[i] *= val;

  return *this;
}


real Tensor::innerProd (const Tensor& T) const
{
  real value = real(0);
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


bool Tensor::isZero (real tol) const
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


real Tensor::trace () const
{
  if (n == 3)
    return v[0] + v[4] + v[8];
  else if (n == 2)
    return v[0] + v[3];
  else if (n == 1)
    return v[0];

  return real(0);
}


real Tensor::det () const
{
  if (n == 3)
    return v[0]*(v[4]*v[8] - v[5]*v[7])
      -    v[3]*(v[1]*v[8] - v[2]*v[7])
      +    v[6]*(v[1]*v[5] - v[2]*v[4]);
  else if (n == 2)
    return v[0]*v[3] - v[1]*v[2];
  else if (n == 1)
    return v[0];

  return real(0);
}


real Tensor::inverse (real tol)
{
  real det = this->det();
  if (det <= tol && det >= -tol)
  {
    std::cerr <<"Tensor::inverse: Singular tensor |T|="<< det << std::endl;
    return real(0);
  }

  if (n == 3)
  {
    real T11 = v[0]; real T12 = v[3]; real T13 = v[6];
    real T21 = v[1]; real T22 = v[4]; real T23 = v[7];
    real T31 = v[2]; real T32 = v[5]; real T33 = v[8];
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
    real T11 = v[0]; real T12 = v[2];
    real T21 = v[1]; real T22 = v[3];
    v[0] =  T22 / det;
    v[1] = -T21 / det;
    v[2] = -T12 / det;
    v[3] =  T11 / det;
  }
  else if (n == 1)
    v[0] = real(1) / det;

  return det;
}


/*!
  \brief Multiplication between a Tensor and a point vector.
*/

Vec3 operator* (const Tensor& T, const Vec3& v)
{
  switch (T.n) {
  case 1:
    return Vec3(T(1,1)*v.x, v.y, v.z);
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
    return Vec3(T(1,1)*v.x, v.y, v.z);
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


std::ostream& SymmTensor::print (std::ostream& os) const
{
  switch (n) {
  case 1:
    return os << v.front() << std::endl;
  case 2:
    return os << v.front() <<'\n'
	      << v.back()  <<' '<< v[1] << std::endl;
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
  v.resize(nsd == 2 && with33 ? 4 : n*(n+1)/2, real(0));
  return true;
}


SymmTensor::SymmTensor (const std::vector<real>& vec) : Tensor(0)
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

  real S11, S12, S13, S21, S22, S23, S31, S32, S33;
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


real SymmTensor::trace () const
{
  real t = real(0);

  if (n == 3 || v.size() == 4)
    t = v[0] + v[1] + v[2];
  else if (n == 2)
    t = v[0] + v[1];
  else if (n == 1)
    t = v[0];

  return t;
}


real SymmTensor::det () const
{
  real d = real(0);

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


real SymmTensor::inverse (real tol)
{
  real det = this->det();
  if (det <= tol && det >= -tol)
  {
    std::cerr <<"SymmTensor::inverse: Singular tensor |T|="<< det << std::endl;
    return real(0);
  }

  if (n == 3)
  {
    real T11 = v[0];
    real T21 = v[3]; real T22 = v[1];
    real T31 = v[5]; real T32 = v[4]; real T33 = v[2];
    v[0] =  (T22*T33 - T32*T32) / det;
    v[1] =  (T11*T33 - T31*T31) / det;
    v[2] =  (T11*T22 - T21*T21) / det;
    v[3] = -(T21*T33 - T31*T32) / det;
    v[4] = -(T11*T32 - T31*T21) / det;
    v[5] =  (T21*T32 - T31*T22) / det;
  }
  else if (n == 2)
  {
    real T11 = v.front();
    real T21 = v.back(); real T22 = v[1];
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
    v.front() = real(1) / det;

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


/*!
  The von Mises value of a symmetric 3D (stress) tensor is defined as follows:
  \f[ s_{\rm vm} = \sqrt{
  s_{11}(s_{11}-s_{22}) +
  s_{22}(s_{22}-s_{33}) +
  s_{33}(s_{33}-s_{11}) + 3(s_{12}^2 + s_{23}^2 + s_{31}^2)}
  \f]
*/

real SymmTensor::vonMises (bool doSqrt) const
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


real SymmTensor::L2norm (bool doSqrt) const
{
  double l2n = 0.0;
  for (t_ind i = 0; i < v.size(); i++)
    l2n += (i < n || (i == 2 && v.size() == 4) ? 1.0 : 2.0)*v[i]*v[i];

  return doSqrt ? sqrt(l2n) : l2n;
}


/*!
  \brief Adding a scaled unit tensor to a symmetric tensor.
*/

SymmTensor operator+ (const SymmTensor& T, real a)
{
  SymmTensor S(T);

  for (Tensor::t_ind i = 0; i < S.n; i++)
    S.v[i] += a;

  if (S.v.size() == 4)
    S.v[2] += a;

  return S;
}


/*!
  \brief Subtracting a scaled unit tensor from a symmetric tensor.
*/

SymmTensor operator- (const SymmTensor& T, real a)
{
  SymmTensor S(T);

  for (Tensor::t_ind i = 0; i < S.n; i++)
    S.v[i] -= a;

  if (S.v.size() == 4)
    S.v[2] -= a;

  return S;
}


/*!
  \brief Multiplication between a scalar and a symmetric tensor.
*/

SymmTensor operator* (real a, const SymmTensor& T)
{
  SymmTensor S(T.dim(), T.v.size() == 4);

  for (Tensor::t_ind i = 0; i < T.v.size(); i++)
    S.v[i] = a*T.v[i];

  return S;
}


SymmTensor4::SymmTensor4 (const std::vector<real>& x, t_ind nsd)
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

  ptr = (real*)&v.front();
}


const real& SymmTensor4::operator() (t_ind i, t_ind j, t_ind k, t_ind l) const
{
  return v[index(i,j)*m+index(k,l)];
}


real& SymmTensor4::operator() (t_ind i, t_ind j, t_ind k, t_ind l)
{
  return ptr[index(i,j)*m+index(k,l)];
}


int LocalSystem::patch = 0;
