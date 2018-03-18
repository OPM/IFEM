// $Id$
//==============================================================================
//!
//! \file LAPack.h
//!
//! \date Jan 10 2016
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief LAPack support for various platforms.
//!
//==============================================================================

#ifndef UTL_LAPACK_H
#define UTL_LAPACK_H

#include "BLAS.h"

#if HAS_BLAS == 3
#define Subroutine static inline void

//! \brief Estimates the reciprocal of the condition number of the matrix \b A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dgecon (char norm, int n, const Real* A, int lda,
                   Real anorm, Real* rcond, Real* work, int* iwork, int& info)
{ dgecon_(&norm,&n,const_cast<Real*>(A),&lda,&anorm,rcond,work,iwork,&info); }

//! \brief Solves the linear equation system \a A*x=b.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dgesv (int n, int nrhs,
                  Real* A, int lda, int* ipiv,
                  Real* B, int ldb, int& info)
{ dgesv_(&n,&nrhs,A,&lda,ipiv,B,&ldb,&info); }

//! \brief Computes an LU factorization of a general M-by-N matrix A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dgetrf (int m, int n, Real* A, int lda,
                   int* ipiv, int& info)
{ dgetrf_(&m,&n,A,&lda,ipiv,&info); }

//! \brief Computes the inverse of a matrix based on its LU factorization.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dgetri (int n, Real* A, int lda,
                   int* ipiv, Real* work, int lwork, int& info)
{ dgetri_(&n,A,&lda,ipiv,work,&lwork,&info); }

//! \brief Solves the equation system \a A*x=b when \a A is already factorized.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dgetrs (char trans, int n, int nrhs,
                   Real* A, int lda, int* ipiv,
                   Real* B, int ldb, int& info)
{ dgetrs_(&trans,&n,&nrhs,A,&lda,ipiv,B,&ldb,&info); }

//! \brief Returns a norm of the matrix \b A.
//! \details This is a FORTRAN-77 function in the LAPack library.
//! \sa LAPack library documentation.
static inline double dlange (char norm, int m, int n, const Real* A,
                             int lda, Real* work)
{ return dlange_(&norm,&m,&n,const_cast<Real*>(A),&lda,work); }

//! \brief Solves the symmetric linear equation system \a A*x=b.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dposv (char uplo, int n, int nrhs,
                  Real* A, int lda, Real* B, int ldb, int& info)
{ dposv_(&uplo,&n,&nrhs,A,&lda,B,&ldb,&info); }

//! \brief Solves the symmetric equation system \a A*x=b for prefactored \b A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dpotrs (char uplo, int n, int nrhs,
                   Real* A, int lda, Real* B, int ldb, int& info)
{ dpotrs_(&uplo,&n,&nrhs,A,&lda,B,&ldb,&info); }

//! \brief Solves the standard eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dsyev (char jobz, char uplo,
                  int n, double* A, int lda, double* w,
                  double* work, int lwork, int& info)
{ dsyev_(&jobz,&uplo,&n,A,&lda,w,work,&lwork,&info); }

//! \brief Solves the standard eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dsyevx (char jobz, char range, char uplo,
                   int n, Real* A, int lda,
                   Real vl, Real vu, int il, int iu,
                   Real abstol, int& m, Real* w, Real* z, int ldz,
                   Real* work, int lwork, int* iwork, int* ifail, int& info)
{ dsyevx_(&jobz,&range,&uplo,&n,A,&lda,&vl,&vu,&il,&iu,&abstol,
          &m,w,z,&ldz,work,&lwork,iwork,ifail,&info); }

//! \brief Solves the generalized eigenproblem \a A*x=(lambda)*B*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dsygvx (int itype, char jobz, char range,
                   char uplo, int n, Real* A, int lda,
                   Real* B, int ldb, Real vl, Real vu,
                   int il, int iu, Real abstol,
                   int& m, Real* w, Real* z, int ldz,
                   Real* work, int lwork, int* iwork, int* ifail, int& info)
{ dsygvx_(&itype,&jobz,&range,&uplo,&n,A,&lda,B,&ldb,&vl,&vu,&il,&iu,&abstol,
          &m,w,z,&ldz,work,&lwork,iwork,ifail,&info); }

//! \brief Solves the non-symmetric eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
Subroutine dgeev (char jobvl, char jobvr,
                  int n, Real* A, int lda,
                  Real* wr, Real* wi, Real* vl, int ldvl,
                  Real* vr, int ldvr,
                  Real* work, int lwork, int& info)
{ dgeev_(&jobvl,&jobvr,&n,A,&lda,wr,wi,vl,&ldvl,vr,&ldvr,work,&lwork,&info); }

#else

#if defined(_WIN32) && !defined(__MINGW32__) && !defined(__MINGW64__)
#define dgecon_ DGECON
#define dgesv_  DGESV
#define dgetrf_ DGETRF
#define dgetri_ DGETRI
#define dgetrs_ DGETRS
#define dlange_ DLANGE
#define dposv_  DPOSV
#define dpotrs_ DPOTRS
#define dsyev_  DSYEV
#define dsyevx_ DSYEVX
#define dsygvx_ DSYGVX
#define dgeev_  DGEEV
#elif defined(_AIX)
#define dgecon_ dgecon
#define dgesv_  dgesv
#define dgetrf_ dgetrf
#define dgetri_ dgetri
#define dgetrs_ dgetrs
#define dlange_ dlange
#define dposv_  dposv
#define dpotrs_ dpotrs
#define dsyev_  dsyev
#define dsyevx_ dsyevx
#define dsygvx_ dsygvx
#define dgeev_  dgeev
#endif

extern "C" {
//! \brief Estimates the reciprocal of the condition number of the matrix \b A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgecon_(const char& norm, const int& n, const Real* A, const int& lda,
             const Real& anorm, Real* rcond, Real* work, int* iwork, int& info);

//! \brief Solves the linear equation system \a A*x=b.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgesv_(const int& n, const int& nrhs,
            Real* A, const int& lda, int* ipiv,
            Real* B, const int& ldb, int& info);

//! \brief Computes an LU factorization of a general M-by-N matrix A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetrf_(const int& m, const int& n, Real* A, const int& lda,
             int* ipiv, int& info);

//! \brief Computes the inverse of a matrix based on its LU factorization.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetri_(const int& n, Real* A, const int& lda,
             int* ipiv, Real* work, const int& lwork, int& info);

//! \brief Solves the equation system \a A*x=b when \a A is already factorized.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetrs_(const char& trans, const int& n, const int& nrhs,
             Real* A, const int& lda, int* ipiv,
             Real* B, const int& ldb, int& info);

//! \brief Returns a norm of the matrix \b A.
//! \details This is a FORTRAN-77 function in the LAPack library.
//! \sa LAPack library documentation.
double dlange_(const char& norm, const int& m, const int& n, const Real* A,
               const int& lda, Real* work);

//! \brief Solves the symmetric linear equation system \a A*x=b.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dposv_(const char& uplo, const int& n, const int& nrhs,
            Real* A, const int& lda, Real* B, const int& ldb, int& info);

//! \brief Solves the symmetric equation system \a A*x=b for prefactored \b A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dpotrs_(const char& uplo, const int& n, const int& nrhs,
             Real* A, const int& lda, Real* B, const int& ldb, int& info);

//! \brief Solves the standard eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dsyev_(const char& jobz, const char& uplo,
            const int& n, double* A, const int& lda, double* w,
            double* work, const int& lwork, int& info);

//! \brief Solves the standard eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dsyevx_(const char& jobz, const char& range, const char& uplo,
             const int& n, Real* a, const int& lda,
             const Real& vl, const Real& vu, const int& il, const int& iu,
             const Real& abstol, int& m, Real* w, Real* z, const int& ldz,
             Real* work, const int& lwork, int* iwork, int* ifail, int& info);

//! \brief Solves the generalized eigenproblem \a A*x=(lambda)*B*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dsygvx_(const int& itype, const char& jobz, const char& range,
             const char& uplo, const int& n, Real* a, const int& lda,
             Real* b, const int& ldb, const Real& vl, const Real& vu,
             const int& il, const int& iu, const Real& abstol,
             int& m, Real* w, Real* z, const int& ldz,
             Real* work, const int& lwork, int* iwork, int* ifail, int& info);

//! \brief Solves the non-symmetric eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgeev_(const char& jobvl, const char& jobvr,
            const int& n, Real* a, const int& lda,
            Real* wr, Real* wi, Real* vl, const int& ldvl,
            Real* vr, const int& ldvr,
            Real* work, const int& lwork, int& info);
}
#endif

#endif
