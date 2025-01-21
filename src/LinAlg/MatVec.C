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


RealArray utl::operator* (const Matrix& A, const Vector& X)
{
  RealArray Y;
  A.multiply(X,Y);
  return Y;
}


RealArray utl::operator* (const Vector& X, const Matrix& A)
{
  RealArray Y;
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
  where \b A is a full, symmetric matrix, and the transformation matrix \b T
  has the nodal sub-matrix \b Tn repeated on the diagonal and otherwise zero.
*/

bool utl::transform (Matrix& A, const Matrix& Tn)
{
  size_t M = A.rows();
  size_t N = Tn.rows();
  if (M < A.cols() || Tn.cols() < N || M < N)
    return false;

  Matrix B(N,N);
  size_t i, j, k;
  for (size_t K = 0; K < M; K += N)
    for (size_t L = 0; L < M; L += N)
    {
      B.fill(Real(0));

      for (i = 1; i <= N; i++)
        for (j = 1; j <= N; j++)
          for (k = 1; k <= N; k++)
            B(j,i) += Tn(j,k)*A(L+k,K+i);

      for (i = 1; i <= N; i++)
        for (k = 1; k <= N; k++)
          A(L+k,K+i) = Real(0);

      for (i = 1; i <= N; i++)
        for (k = 1; k <= N; k++)
          for (j = 1; j <= N; j++)
            A(L+k,K+i) += B(k,j)*Tn(i,j);
    }

  return true;
}


/*!
  The vector \b V is pre-multiplied with the transformation matrix \b T which
  has the nodal sub-matrix \b Tn (or its transpose) repeated on the diagonal
  and otherwise zero.
*/

bool utl::transform (Vector& V, const Matrix& Tn, bool transpose)
{
  size_t M = V.size();
  size_t N = Tn.rows();
  if (Tn.cols() < N || M < N)
    return false;

  Vector WA(N);
  for (size_t K = 0; K < M; K += N)
  {
    WA.fill(Real(0));

    for (size_t i = 1; i <= N; i++)
      for (size_t j = 1; j <= N; j++)
        WA(i) += (transpose ? Tn(j,i) : Tn(i,j)) * V(K+j);

    for (size_t j = 1; j <= N; j++)
      V(K+j) = WA(j);
  }

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
