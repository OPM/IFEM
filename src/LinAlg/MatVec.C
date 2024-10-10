// $Id$
//==============================================================================
//!
//! \file MatVec.C
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Global algebraic operations on index 1-based matrices and vectors.
//!
//! The two transformation functions defined in this file are rewritten versions
//! of the Fortran subroutines MATTRA and VECTRA from the SAM library.
//!
//==============================================================================

#include "MatVec.h"
#include "LAPack.h"


// Set default values on some global variables.

int utl::nval_per_line = 6;
double utl::zero_print_tol = 1.0e-6;


// Global algebraic matrix-vector operators.

Vector utl::operator* (const Vector& X, Real c)
{
  Vector result(X);
  return result *= c;
}


Vector utl::operator+ (const Vector& X, const Vector& Y)
{
  Vector result(X);
  return result.add(Y);
}


Vector utl::operator- (const Vector& X, const Vector& Y)
{
  Vector result(X);
  return result.add(Y,Real(-1));
}


Matrix utl::operator* (const Matrix& A, Real c)
{
  Matrix B(A);
  return B.multiply(c);
}


Vector utl::operator* (const Matrix& A, const Vector& X)
{
  Vector Y;
  A.multiply(X,Y);
  return Y;
}


Vector utl::operator* (const Vector& X, const Matrix& A)
{
  Vector Y;
  A.multiply(X,Y,true);
  return Y;
}


Matrix utl::operator* (const Matrix& A, const Matrix& B)
{
  Matrix C;
  C.multiply(A,B);
  return C;
}


/*!
  The following matrix multiplication is performed by this function:
  \f[ {\bf A} = {\bf T}{\bf A}{\bf T}^T \f]
  where \b A is a full, symmetric matrix, and \b T is an identity matrix
  with the nodal sub-matrix \b Tn inserted on the diagonal.
*/

bool utl::transform (Matrix& A, const Matrix& Tn, size_t K)
{
  size_t M = A.rows();
  size_t N = Tn.rows();
  if (M < A.cols() || N > Tn.cols()) return false;
  if (K < 1 || K+N-1 > M || N > 3) return false;

  Real* WA = (Real*)alloca(N*sizeof(Real));
  size_t i, ii, j, jj, l, KN = K+N-1;
  for (jj = K; jj <= M; jj++)
  {
    for (i = 1; i <= N; i++)
    {
      WA[i-1] = Real(0);
      for (l = 1, ii = K; l <= N; l++, ii++)
        WA[i-1] += Tn(i,l)*A(ii,jj);
    }
    ii = K;
    for (i = 1; i <= N; i++, ii++)
    {
      A(ii,jj) = WA[i-1];
      if (jj > KN) A(jj,ii) = WA[i-1];
    }
  }

  for (ii = 1; ii <= KN; ii++)
  {
    size_t JS = ii > K ? ii-K+1 : 1;
    for (j = JS; j <= N; j++)
    {
      WA[j-1] = Real(0);
      for (l = 1, jj = K; l <= N; l++, jj++)
        WA[j-1] += A(ii,jj)*Tn(j,l);
    }
    jj = ii > K ? ii : K;
    for (j = JS; j <= N; j++, jj++)
    {
      A(ii,jj) = WA[j-1];
      if (ii < K) A(jj,ii) = WA[j-1];
    }
  }

  for (ii = K; ii <= KN; ii++)
    for (jj = ii; jj <= KN; jj++)
      A(jj,ii) = A(ii,jj);

  return true;
}


/*!
  The vector \b V is pre-multiplied with the transformation matrix \b T which is
  the identity matrix with the nodal sub-matrix \b Tn inserted on the diagonal.
*/

bool utl::transform (Vector& V, const Matrix& Tn, size_t K, bool transpose)
{
  size_t M = V.size();
  size_t N = Tn.rows();
  if (N > Tn.cols()) return false;
  if (K < 1 || K+N-1 > M || N > 3) return false;

  Real* WA = (Real*)alloca(N*sizeof(Real));
  size_t i, ii, j;
  if (transpose)
    for (i = 1; i <= N; i++)
    {
      WA[i-1] = Real(0);
      for (j = 1, ii = K; j <= N; j++, ii++)
        WA[i-1] += Tn(j,i)*V(ii);
    }

  else
    for (i = 1; i <= N; i++)
    {
      WA[i-1] = Real(0);
      for (j = 1, ii = K; j <= N; j++, ii++)
        WA[i-1] += Tn(i,j)*V(ii);
    }

  for (i = 0, ii = K; i < N; i++, ii++)
    V(ii) = WA[i];

  return true;
}


bool utl::invert (Matrix& A)
{
  if (A.rows() <= 3)
    return A.inverse() > 0.0;

#ifdef HAS_BLAS
  if (sizeof(Real) != sizeof(double))
  {
    std::cerr <<" *** utl::invert: Available for double precision matrices"
              <<" only."<< std::endl;
    return false;
  }

  // Use LAPack/BLAS for larger matrices
  int INFO, N = A.rows() > A.cols() ? A.rows() : A.cols();
  int* IPIV = (int*)alloca(N*sizeof(int));
  dgetrf_(N,N,A.ptr(),A.rows(),IPIV,INFO);
  if (INFO != 0)
  {
    std::cerr <<" *** utl::invert:DGETRF: INFO = "<< INFO << std::endl;
    return false;
  }
  double NWORK;
  dgetri_(N,A.ptr(),A.rows(),IPIV,&NWORK,-1,INFO);
  double* WORK = new double[int(NWORK)];
  dgetri_(N,A.ptr(),A.rows(),IPIV,WORK,int(NWORK),INFO);
  delete[] WORK;
  if (INFO == 0) return true;
  std::cerr <<" *** utl::invert:DGETRI: INFO = "<< INFO << std::endl;
#else
  std::cerr <<"DGETRI not available - linked without LAPack/BLAS"<< std::endl;
#endif
  return false;
}
