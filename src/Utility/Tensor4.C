// $Id$
//==============================================================================
//!
//! \file Tensor4.C
//!
//! \date Oct 27 2015
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of fourth-order tensors with some basic operations.
//!
//==============================================================================

#include "Tensor4.h"
#include "matrix.h"


Tensor4::Tensor4 (t_ind nd, Real scale, bool makeJ)
{
  this->redim(nd);

  v.resize(m*m,Real(0));
  if (scale != Real(0))
  {
    if (makeJ) // Define the "identity" tensor, J
      for (t_ind i = 1; i <= n; i++)
        for (t_ind j = 1; j <= n; j++)
          v[this->index(i,i,j,j)] = scale;

    else // Define the identity tensor, I
      for (t_ind i = 1; i <= n; i++)
        v[this->index(i,i,i,i)] = scale;
  }
}


Tensor4::Tensor4 (const std::vector<Real>& x, t_ind nd) : v(x)
{
  this->redim(nd);

  if (v.size() < (size_t)m*m)
    std::cerr <<" *** Invalid fourth-order tensor,"
              <<" matrix represention too small, size="<< v.size() << std::endl;
  else if (v.size() > (size_t)m*m)
    v.resize(m*m);
}


std::ostream& Tensor4::print (std::ostream& os) const
{
  utl::matrix<Real> A(m,m);
  A.fill(&v.front());
  return os << A;;
}


const Real& Tensor4::operator() (t_ind i, t_ind j, t_ind k, t_ind l) const
{
  return v[this->index(i,j,k,l)];
}


Real& Tensor4::operator() (t_ind i, t_ind j, t_ind k, t_ind l)
{
  return v[this->index(i,j,k,l)];
}


Tensor4& Tensor4::operator= (const Tensor4& T)
{
  this->redim(T.n);
  v.resize(m*m,Real(0));
  std::copy(T.v.begin(),T.v.end(),v.begin());

  return *this;
}


Tensor4& Tensor4::operator= (Real val)
{
  this->zero();

  for (t_ind i = 1; i <= n; i++)
    v[this->index(i,i,i,i)] = val;

  return *this;
}


Tensor4& Tensor4::operator+= (Real val)
{
  for (t_ind i = 1; i <= n; i++)
    v[this->index(i,i,i,i)] += val;

  return *this;
}


SymmTensor4::SymmTensor4 (t_ind nd, bool makeJ) : Tensor4(nd,Real(0))
{
  if (makeJ) // Define the "identity" tensor, J
    for (t_ind i = 0; i < n; i++)
      for (t_ind j = 0; j < n; j++)
        v[i*m+j] = Real(1);

  else // Define the identity tensor, I
    for (t_ind i = 0; i < n; i++)
      v[i*m+i] = Real(1);
}


SymmTensor4::SymmTensor4 (const std::vector<Real>& x, t_ind nd) : Tensor4(x,nd)
{
}


std::ostream& SymmTensor4::print (std::ostream& os) const
{
  utl::matrix<Real> A(m,m);
  for (t_ind l = 1; l <= n; l++)
    for (t_ind k = 1; k <= n; k++)
      for (t_ind j = 1; j <= n; j++)
        for (t_ind i = 1; i <= n; i++)
          A(i+n*(j-1),k+l*(n-1)) = v[this->index(i,j)*m+this->index(k,l)];

  return os << A;;
}


void SymmTensor4::redim (t_ind ndim)
{
  switch (n = ndim) {
  case 1: m = 1; break;
  case 2: m = 3; break;
  case 3: m = 6; break;
  default: m = 0;
    std::cerr <<" *** Invalid fourth-order tensor, dim="<< n << std::endl;
  }
}


const Real& SymmTensor4::operator() (t_ind i, t_ind j, t_ind k, t_ind l) const
{
  return v[this->index(i,j)*m+this->index(k,l)];
}


Real& SymmTensor4::operator() (t_ind i, t_ind j, t_ind k, t_ind l)
{
  return v[this->index(i,j)*m+this->index(k,l)];
}
