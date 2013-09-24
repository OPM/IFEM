#include "matrix.h"

utl::vector<Real> utl::operator*(const utl::vector<Real>& X, Real c)
{
  utl::vector<Real> result(X);
  return result *= c;
}

utl::vector<Real> utl::operator*(Real c, const utl::vector<Real>& X) 
{
  utl::vector<Real> result(X);
  return result *= c;
}

Real utl::operator*(const utl::vector<Real>& X, const utl::vector<Real>& Y) 
{
  return X.dot(Y);
}

utl::vector<Real> utl::operator/(const utl::vector<Real>& X, Real d)
{
  utl::vector<Real> result(X);
  return result /= d;
}

utl::vector<Real> utl::operator+(const utl::vector<Real>& X, const utl::vector<Real>& Y)
{
  utl::vector<Real> result(X);
  return result += Y;
}

utl::vector<Real> utl::operator-(const utl::vector<Real>& X, const utl::vector<Real>& Y)
{
  utl::vector<Real> result(X);
  return result -= Y;
}

utl::matrix<Real> utl::operator*(const utl::matrix<Real>& A, Real c)
{
  utl::matrix<Real> B(A);
  return B.multiply(c);
}

utl::matrix<Real> utl::operator*(Real c, const utl::matrix<Real>& A)
{
  utl::matrix<Real> B(A);
  return B.multiply(c);
}

utl::vector<Real> utl::operator*(const utl::matrix<Real>& A, const utl::vector<Real>& X)
{
  utl::vector<Real> Y;
  A.multiply(X,Y);
  return Y;
}

utl::matrix<Real> utl::operator*(const utl::matrix<Real>& A, const utl::matrix<Real>& B)
{
  utl::matrix<Real> C;
  C.multiply(A, B);
  return C;
}
