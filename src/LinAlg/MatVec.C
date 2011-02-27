// $Id$
//==============================================================================
//!
//! \file MatVec.C
//!
//! \date Oct 1 2007
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Index 1-based matrices and vectors for algebraic operations.
//!
//! The two transformation functions defined in this file are rewritten versions
//! of the Fortran subroutines MATTRA and VECTRA from the SAM library.
//!
//==============================================================================

#include "MatVec.h"

#if defined(_WIN32)
#define dgetrf_ DGETRF
#define dgetri_ DGITRI
#elif defined(_AIX)
#define dgetrf_ dgetrf
#define dgetri_ dgetri
#endif

extern "C" {
//! \brief Computes an LU factorization of a general M-by-N matrix A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetrf_(const int& M, const int& N, double* A, const int& LDA,
             int* IPIV, int& INFO);
//! \brief Computes the inverse of a matrix based on its LU factorization.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetri_(const int& N, double* A, const int& LDA,
             int* IPIV, double* WORK, const int& LWORK, int& INFO);
}


// Set default values on some global variables.

int utl::nval_per_line = 6;
double utl::zero_print_tol = 1.0e-6;


bool utl::transform (Matrix& A, const Matrix& T, size_t K)
{
  size_t M = A.rows();
  size_t N = T.rows();
  if (M < A.cols() || N < T.cols()) return false;
  if (K < 1 || K+N-1 > M || N > 3) return false;

  real WA[3];
  size_t i, ii, j, jj, l, KN = K+N-1;
  for (jj = K; jj <= M; jj++)
  {
    for (i = 1; i <= N; i++)
    {
      WA[i-1] = real(0);
      for (l = 1, ii = K; l <= N; l++, ii++)
	WA[i-1] += T(l,i)*A(ii,jj);
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
      WA[j-1] = real(0);
      for (l = 1, jj = K; l <= N; l++, jj++)
	WA[j-1] += A(ii,jj)*T(l,j);
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


bool utl::transform (Vector& V, const Matrix& T, size_t K, bool transpose)
{
  size_t M = V.size();
  size_t N = T.rows();
  if (N < T.cols()) return false;
  if (K < 1 || K+N-1 > M || N > 3) return false;

  real WA[3];
  size_t i, ii, j;
  if (transpose)
    for (i = 1; i <= N; i++)
    {
      WA[i-1] = real(0);
      for (j = 1, ii = K; j <= N; j++, ii++)
	WA[i-1] += T(j,i)*V(ii);
    }

  else
    for (i = 1; i <= N; i++)
    {
      WA[i-1] = real(0);
      for (j = 1, ii = K; j <= N; j++, ii++)
	WA[i-1] += T(i,j)*V(ii);
    }

  for (i = 0, ii = K; i < N; i++, ii++)
    V(ii) = WA[i];

  return true;
}


bool utl::invert (Matrix& A)
{
  if (A.rows() <= 3)
    return A.inverse() > 0.0;

#ifdef USE_CBLAS
  // Use LAPack/BLAS for larger matrices
  int INFO, N = A.rows() > A.cols() ? A.rows() : A.cols();
  int* IPIV = new int[N];
  dgetrf_(N,N,A.ptr(),A.rows(),IPIV,INFO);
  if (INFO != 0)
  {
    delete[] IPIV;
    std::cerr <<" *** utl::invert:DGETRF: INFO = "<< INFO << std::endl;
    return false;
  }
  double NWORK;
  dgetri_(N,A.ptr(),A.rows(),IPIV,&NWORK,-1,INFO);
  double* WORK = new double[int(NWORK)];
  dgetri_(N,A.ptr(),A.rows(),IPIV,WORK,int(NWORK),INFO);
  delete[] IPIV;
  delete[] WORK;
  if (INFO == 0) return true;
  std::cerr <<" *** utl::invert:DGETRI: INFO = "<< INFO << std::endl;
#else
  std::cerr <<"DGETRI not available - linked without LAPack/BLAS"<< std::endl;
#endif
  return false;
}
