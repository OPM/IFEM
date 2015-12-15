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
#define dgecon DGECON
#define dgesv  DGESV
#define dgetrf DGETRF
#define dgetri DGETRI
#define dgetrs DGETRS
#define dlange DLANGE
#define dposv  DPOSV
#define dpotrs DPOTRS
#define dsyev  DSYEV
#define dsyevx DSYEVX
#define dsygvx DSYGVX
#define dgeev  DGEEV
#elif !defined(_AIX)
#define dgecon dgecon_
#define dgesv  dgesv_
#define dgetrf dgetrf_
#define dgetri dgetri_
#define dgetrs dgetrs_
#define dlange dlange_
#define dposv  dposv_
#define dpotrs dpotrs_
#define dsyev  dsyev_
#define dsyevx dsyevx_
#define dsygvx dsygvx_
#define dgeev  dgeev_
#endif

extern "C" {
//! \brief Estimates the reciprocal of the condition number of the matrix \b A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgecon (const char& norm, const int& n, const Real* A, const int& lda,
             const Real& anorm, Real* rcond, Real* work, int* iwork, int& info);

//! \brief Solves the linear equation system \a A*x=b.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgesv (const int& n, const int& nrhs,
            Real* A, const int& lda, int* ipiv,
            Real* B, const int& ldb, int& info);

//! \brief Computes an LU factorization of a general M-by-N matrix A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetrf (const int& m, const int& n, Real* A, const int& lda,
             int* ipiv, int& info);

//! \brief Computes the inverse of a matrix based on its LU factorization.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetri (const int& n, Real* A, const int& lda,
             int* ipiv, Real* work, const int& lwork, int& info);

//! \brief Solves the equation system \a A*x=b when \a A is already factorized.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetrs (const char& trans, const int& n, const int& nrhs,
             Real* A, const int& lda, int* ipiv,
             Real* B, const int& ldb, int& info);

//! \brief Returns a norm of the matrix \b A.
//! \details This is a FORTRAN-77 function in the LAPack library.
//! \sa LAPack library documentation.
double dlange (const char& norm, const int& m, const int& n, const Real* A,
               const int& lda, Real* work);

//! \brief Solves the symmetric linear equation system \a A*x=b.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dposv (const char& uplo, const int& n, const int& nrhs,
            Real* A, const int& lda, Real* B, const int& ldb, int& info);

//! \brief Solves the symmetric equation system \a A*x=b for prefactored \b A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dpotrs (const char& uplo, const int& n, const int& nrhs,
             Real* A, const int& lda, Real* B, const int& ldb, int& info);

//! \brief Solves the standard eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dsyev (const char& jobz, const char& uplo,
            const int& n, double* A, const int& lda, double* w,
            double* work, const int& lwork, int& info);

//! \brief Solves the standard eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dsyevx (const char& jobz, const char& range, const char& uplo,
             const int& n, Real* a, const int& lda,
             const Real& vl, const Real& vu, const int& il, const int& iu,
             const Real& abstol, int& m, Real* w, Real* z, const int& ldz,
             Real* work, const int& lwork, int* iwork, int* ifail, int& info);

//! \brief Solves the generalized eigenproblem \a A*x=(lambda)*B*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dsygvx (const int& itype, const char& jobz, const char& range,
             const char& uplo, const int& n, Real* a, const int& lda,
             Real* b, const int& ldb, const Real& vl, const Real& vu,
             const int& il, const int& iu, const Real& abstol,
             int& m, Real* w, Real* z, const int& ldz,
             Real* work, const int& lwork, int* iwork, int* ifail, int& info);

//! \brief Solves the non-symmetric eigenproblem \a A*x=(lambda)*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgeev (const char& jobvl, const char& jobvr,
            const int& n, Real* a, const int& lda,
            Real* wr, Real* wi, Real* vl, const int& ldvl,
            Real* vr, const int& ldvr,
            Real* work, const int& lwork, int& info);
}
#endif

#endif
