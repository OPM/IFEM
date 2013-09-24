// $Id$
//==============================================================================
//!
//! \file DenseMatrix.C
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Dense system matrix representation.
//!
//==============================================================================

#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "SAM.h"

#if defined(_WIN32) && !defined(__MINGW32__) && !defined(__MINGW64__)
#define addem2_ ADDEM2
#define dgecon_ DGECON
#define dgesv_  DGESV
#define dgetrs_ DGETRS
#define dlange_ DLANGE
#define dposv_  DPOSV
#define dpotrs_ DPOTRS
#define dsyevx_ DSYEVX
#define dsygvx_ DSYGVX
#define dgeev_  DGEEV
#elif defined(_AIX)
#define addem2_ addem2
#define dgecon_ dgecon
#define dgesv_  dgesv
#define dgetrs_ dgetrs
#define dlange_ dlange
#define dposv_  dposv
#define dpotrs_ dpotrs
#define dsyevx_ dsyevx
#define dsygvx_ dsygvx
#define dgeev_  dgeev
#endif

extern "C" {
//! \brief Adds an element matrix \a eK into the system matrix \a sysK.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void addem2_(const Real* eK, const Real* ttcc, const int* mpar,
             const int* madof, const int* meqn, const int* mpmnpc,
             const int* mmnpc, const int* mpmceq, const int* mmceq,
             const int& iel, const int& nedof, const int& neq,
             const int& lpu, const int& nrhs, Real* sysK, Real* sysRHS,
             int* work, int& ierr);

//! \brief Estimates the reciprocal of the condition number of the matrix \b A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgecon_(const char& norm, const int& n, const Real* A, const int& lda,
             const Real& anorm, Real& rcond, Real* work, int* iwork, int& info);

//! \brief Solves the linear equation system \a A*x=b.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgesv_(const int& n, const int& nrhs,
	    Real* A, const int& lda, int* ipiv,
	    Real* B, const int& ldb, int& info);

//! \brief Solves the equation system \a A*x=b when \a A is already factorized.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dgetrs_(const char& trans, const int& n, const int& nrhs,
	     Real* A, const int& lda, int* ipiv,
	     Real* B, const int& ldb, int& info);

//! \brief Returns a norm of the matrix \b A.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
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
void dsyevx_(const char& jobz, const char& range, const char& uplo,
	     const int& n, Real* a, const int& lda,
	     const Real& vl, const Real& vu, const int& il, const int& iu,
	     const Real& abstol, const int& m, Real* w, Real* z, const int& ldz,
	     Real* work, const int& lwork, int* iwork, int* ifail, int& info);

//! \brief Solves the generalized eigenproblem \a A*x=(lambda)*B*x.
//! \details This is a FORTRAN-77 subroutine in the LAPack library.
//! \sa LAPack library documentation.
void dsygvx_(const int& itype, const char& jobz, const char& range,
	     const char& uplo, const int& n, Real* a, const int& lda,
	     Real* b, const int& ldb, const Real& vl, const Real& vu,
	     const int& il, const int& iu, const Real& abstol,
	     const int& m, Real* w, Real* z, const int& ldz,
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


DenseMatrix::DenseMatrix (size_t m, size_t n, bool s) : myMat(m,n), ipiv(NULL)
{
  symm = s && m == n;
}


DenseMatrix::DenseMatrix (const DenseMatrix& A) : ipiv(NULL)
{
  myMat = A.myMat;
  symm = A.symm;
  if (A.ipiv)
    std::cerr <<"DenseMatrix constructor: Copying factored matrix"<< std::endl;
}


DenseMatrix::DenseMatrix (const RealArray& data, size_t nrows) : ipiv(NULL)
{
  size_t ndata = data.size();
  if (nrows == 0) nrows = (size_t)sqrt((double)ndata);
  size_t ncols = nrows ? ndata/nrows : 0;

  myMat.resize(nrows,ncols);
  memcpy(myMat.ptr(),&data.front(),nrows*ncols*sizeof(Real));
  symm = false;
}


size_t DenseMatrix::dim (int idim) const
{
  if (idim == 1)
    return myMat.rows();
  else if (idim == 2)
    return myMat.cols();
  else
    return myMat.size();
}


void DenseMatrix::initAssembly (const SAM& sam, bool)
{
  myMat.resize(sam.neq,sam.neq,true);
}


void DenseMatrix::init ()
{
  myMat.fill(Real(0));

  // Delete pivotation vector of old factorization, if any
  delete[] ipiv;
  ipiv = NULL;
}


void DenseMatrix::dump (std::ostream& os, char format, const char* label)
{
  switch (format)
    {
    case 'M':
    case 'm':
      utl::writeMatlab(label,myMat,os);
      break;

    default:
      if (label) os << label <<" =";
      os << myMat;
    }
}


/*!
  \brief This is a C++ version of the F77 subroutine ADDEM2 (SAM library).
  \details It performs exactly the same tasks, except that \a NRHS always is 1.
*/
#ifndef USE_F77SAM
static void assemDense (const Matrix& eM, Matrix& SM, Vector& SV,
			const std::vector<int>& meen, const int* meqn,
			const int* mpmceq, const int* mmceq, const Real* ttcc)
{
  // Add elements corresponding to free dofs in eM into SM
  int i, j, ip, nedof = meen.size();
  for (j = 1; j <= nedof; j++)
  {
    int jeq = meen[j-1];
    if (jeq < 1) continue;

    SM(jeq,jeq) += eM(j,j);

    for (i = 1; i < j; i++)
    {
      int ieq = meen[i-1];
      if (ieq < 1) continue;

      SM(ieq,jeq) += eM(i,j);
      SM(jeq,ieq) += eM(j,i);
    }
  }

  // Add (appropriately weighted) elements corresponding to constrained
  // (dependent and prescribed) dofs in eM into SM and/or SV
  for (j = 1; j <= nedof; j++)
  {
    int jceq = -meen[j-1];
    if (jceq < 1) continue;

    int jp = mpmceq[jceq-1];
    Real c0 = ttcc[jp-1];

    // Add contributions to SV (right-hand-side)
    if (!SV.empty())
      for (i = 1; i <= nedof; i++)
      {
	int ieq = meen[i-1];
	int iceq = -ieq;
	if (ieq > 0)
	  SV(ieq) -= c0*eM(i,j);
	else if (iceq > 0)
	  for (ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
	    if (mmceq[ip] > 0)
	    {
	      ieq = meqn[mmceq[ip]-1];
	      SV(ieq) -= c0*ttcc[ip]*eM(i,j);
	    }
      }

    // Add contributions to SM
    for (jp = mpmceq[jceq-1]; jp < mpmceq[jceq]-1; jp++)
      if (mmceq[jp] > 0)
      {
	int jeq = meqn[mmceq[jp]-1];
	for (i = 1; i <= nedof; i++)
	{
	  int ieq = meen[i-1];
	  int iceq = -ieq;
	  if (ieq > 0)
	  {
	    SM(ieq,jeq) += ttcc[jp]*eM(i,j);
	    SM(jeq,ieq) += ttcc[jp]*eM(j,i);
	  }
	  else if (iceq > 0)
	    for (ip = mpmceq[iceq-1]; ip < mpmceq[iceq]-1; ip++)
	      if (mmceq[ip] > 0)
	      {
		ieq = meqn[mmceq[ip]-1];
		SM(ieq,jeq) += ttcc[ip]*ttcc[jp]*eM(i,j);
	      }
	}
      }
  }
}
#endif

bool DenseMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  if (myMat.rows() != (size_t)sam.neq || myMat.cols() < (size_t)sam.neq)
    return false;

  int ierr = 0;
#ifdef USE_F77SAM
  Real dummyRHS;
  int* work = new int[eM.rows()];
  addem2_(eM.ptr(), sam.ttcc, sam.mpar,
	  sam.madof, sam.meqn, sam.mpmnpc, sam.mmnpc, sam.mpmceq, sam.mmceq,
	  e, eM.rows(), sam.neq, 6, 0, myMat.ptr(), &dummyRHS, work, ierr);
  delete[] work;
#else
  Vector dummyB;
  std::vector<int> meen;
  if (sam.getElmEqns(meen,e,eM.rows()))
    assemDense(eM,myMat,dummyB,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);
  else
    ierr = 1;
#endif
  return ierr == 0;
}


bool DenseMatrix::assemble (const Matrix& eM, const SAM& sam,
			    SystemVector& B, int e)
{
  if (myMat.rows() != (size_t)sam.neq || myMat.cols() < (size_t)sam.neq)
    return false;

  StdVector* Bptr = dynamic_cast<StdVector*>(&B);
  if (!Bptr) return false;

  int ierr = 0;
#ifdef USE_F77SAM
  int* work = new int[eM.rows()];
  addem2_(eM.ptr(), sam.ttcc, sam.mpar,
	  sam.madof, sam.meqn, sam.mpmnpc, sam.mmnpc, sam.mpmceq, sam.mmceq,
	  e, eM.rows(), sam.neq, 6, 1, myMat.ptr(), Bptr->ptr(), work, ierr);
  delete[] work;
#else
  std::vector<int> meen;
  if (sam.getElmEqns(meen,e,eM.rows()))
    assemDense(eM,myMat,*Bptr,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);
  else
    ierr = 1;
#endif
  return ierr == 0;
}


bool DenseMatrix::assemble (const Matrix& eM, const SAM& sam,
			    SystemVector& B, const std::vector<int>& meen)
{
#ifndef USE_F77SAM
  if (myMat.rows() != (size_t)sam.neq || myMat.cols() < (size_t)sam.neq)
    return false;

  if (eM.rows() < meen.size() || eM.cols() < meen.size())
    return false;

  StdVector* Bptr = dynamic_cast<StdVector*>(&B);
  if (!Bptr) return false;

  assemDense(eM,myMat,*Bptr,meen,sam.meqn,sam.mpmceq,sam.mmceq,sam.ttcc);
  return true;
#else
  return this->SystemMatrix::assemble(eM,sam,B,meen); // for error message
#endif
}


bool DenseMatrix::augment (const SystemMatrix& B, size_t r0, size_t c0)
{
  const SparseMatrix* sB = dynamic_cast<const SparseMatrix*>(&B);
  if (sB) return this->augment(*sB,r0,c0);

  const DenseMatrix* dB = dynamic_cast<const DenseMatrix*>(&B);
  if (dB) return this->augment(dB->myMat,r0,c0);

  return false;
}


bool DenseMatrix::augment (const Matrix& B, size_t r0, size_t c0)
{
  size_t newRow = r0 + B.rows();
  size_t newCol = c0 + B.cols();
  size_t newDim = newRow > newCol ? newRow : newCol;

  this->redim(newDim,newDim);
  for (size_t c = 1; c <= B.cols(); c++)
    for (size_t r = 1; r <= B.rows(); r++)
    {
      myMat(r0+r,c0+c) += B(r,c);
      myMat(c0+c,r0+r) += B(r,c);
  }

  return true;
}


bool DenseMatrix::augment (const SparseMatrix& B, size_t r0, size_t c0)
{
  const ValueMap& elem = B.getValues();
  if (elem.empty()) return false;

  size_t newRow = r0 + B.rows();
  size_t newCol = c0 + B.cols();
  size_t newDim = newRow > newCol ? newRow : newCol;

  this->redim(newDim,newDim);
  for (ValueIter it = elem.begin(); it != elem.end(); it++)
  {
    myMat(r0+it->first.first,c0+it->first.second) += it->second;
    myMat(c0+it->first.second,r0+it->first.first) += it->second;
  }

  return true;
}


bool DenseMatrix::redim (size_t r, size_t c)
{
  if (r == myMat.rows())
    if (c == myMat.cols())
      return false;
    else if (c > myMat.cols())
    {
      myMat.resize(r,c);
      return true;
    }

  Matrix tmp(myMat);
  myMat.resize(r,c);
  for (size_t i = 1; i <= c && i <= tmp.cols(); i++)
    myMat.fillColumn(i,tmp.getColumn(i));

  if (r != c) symm = false;

  return true;
}


bool DenseMatrix::add (const SystemMatrix& B, Real alpha)
{
  const DenseMatrix* Bptr = dynamic_cast<const DenseMatrix*>(&B);
  if (!Bptr) return false;

  if (myMat.rows() != Bptr->myMat.rows()) return false;
  if (myMat.cols() <  Bptr->myMat.cols()) return false;

  myMat.add(Bptr->myMat,alpha);
  return true;
}


bool DenseMatrix::add (Real sigma)
{
  Real* v = myMat.ptr();
  size_t inc = myMat.rows()+1;
  for (size_t i = 0; i < myMat.size(); i += inc)
    v[i] += sigma;
  return true;
}


bool DenseMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  const StdVector* Bptr = dynamic_cast<const StdVector*>(&B);
  if (!Bptr) return false;
  StdVector*       Cptr = dynamic_cast<StdVector*>(&C);
  if (!Cptr) return false;

  return myMat.multiply(*Bptr,*Cptr);
}


bool DenseMatrix::solve (SystemVector& B, bool, Real* rc)
{
  size_t nrhs = myMat.rows() > 0 ? B.dim()/myMat.rows() : 1;
  return this->solve(B.getPtr(),nrhs,rc);
}


bool DenseMatrix::solve (Matrix& B)
{
  return this->solve(B.ptr(),B.cols());
}


bool DenseMatrix::solve (Real* B, size_t nrhs, Real* rcond)
{
  const size_t n = myMat.rows();
  if (n < 1 || nrhs < 1) return true; // Nothing to solve
  if (n > myMat.cols()) return false; // More equations than unknowns

  const char* dsolv = symm ? "DGESV" : "DPOSV";
#ifdef USE_CBLAS
  int info = 0;
  if (symm)
  {
    // Matrix is marked as symmetric - use Cholesky instead of full LU
    if (ipiv)
      dpotrs_('U',n,nrhs,myMat.ptr(),n,B,n,info);
    else
    {
      ipiv = new int[1]; // dummy allocation to flag factorization reuse
      dposv_('U',n,nrhs,myMat.ptr(),n,B,n,info);
    }
  }
  else if (ipiv)
    dgetrs_('N',n,nrhs,myMat.ptr(),n,ipiv,B,n,info);
  else
  {
    Real anorm;
    if (rcond) // Evaluate the 1-norm of the original LHS-matrix
      anorm = dlange_('1',n,n,myMat.ptr(),n,rcond);
    ipiv = new int[n];
    dgesv_(n,nrhs,myMat.ptr(),n,ipiv,B,n,info);
    if (rcond && info == 0)
    {
      // Estimate the condition number
      Real* work = new Real[4*n];
      int* iwork = new int[n];
      dgecon_('1',n,myMat.ptr(),n,anorm,*rcond,work,iwork,info);
      delete[] work;
      delete[] iwork;
    }
  }
  if (info == 0) return true;

  std::cerr <<"LAPACK::"<< dsolv <<": ";
  if (info < 0)
    std::cerr <<"Invalid argument #"<< -info << std::endl;
  else
    std::cerr <<"Singular stiffness matrix, pivot "<< info
	      <<" (of total "<< n <<") is zero."<< std::endl;
#else
  std::cerr << dsolv <<" not available - built without LAPack/BLAS"<< std::endl;
#endif
  return false;
}


bool DenseMatrix::solveEig (RealArray& val, Matrix& vec, int nv)
{
  const size_t n = myMat.rows();
  if (n < 1 || nv < 1) return true; // No equations to solve
  if (n > myMat.cols()) return false;

#ifdef USE_CBLAS
  std::cout <<"  Solving dense eigenproblem using LAPACK::DSYEVX"<< std::endl;
  int m, info = 0;
  Real dummy = 0.0;
  Real abstol = 0.0;
  // Invoke with Lwork = -1 to estimate work space size
  dsyevx_('V','I','U',n,myMat.ptr(),n,dummy,dummy,1,nv,
          abstol,m,&val.front(),vec.ptr(),n,&dummy,-1,0,0,info);

  if (info == 0)
  {
    // Allocate work space
    int  Lwork = int(dummy);
    Real* work = new Real[Lwork];
    int* Iwork = new int[6*n];
    val.resize(n);
    vec.resize(n,nv);
    // Solve the eigenproblem
    dsyevx_('V','I','U',n,myMat.ptr(),n,dummy,dummy,1,nv,
	    abstol,m,&val.front(),vec.ptr(),n,work,Lwork,Iwork+n,Iwork,info);
    delete[] work;
    delete[] Iwork;
    val.resize(nv);
    if (info == 0) return true;
  }

  std::cerr <<"LAPACK::DSYEVX: ";
  if (info < 0)
    std::cerr <<"Invalid argument #"<< -info << std::endl;
  else
    std::cerr << info <<" eigenvectors failed to converge."<< std::endl;
#else
  std::cerr <<"DSYEVX not available - built without LAPack/BLAS"<< std::endl;
#endif
  return false;
}


bool DenseMatrix::solveEig (DenseMatrix& B, RealArray& val, Matrix& vec, int nv,
			    Real)
{
  const size_t n = myMat.rows();
  if (n < 1 || nv < 1) return true; // No equations to solve
  if (n > myMat.cols()) return false;

#ifdef USE_CBLAS
  std::cout <<"  Solving dense eigenproblem using LAPACK::DSYGVX"<< std::endl;
  int m, info = 0;
  Real dummy = 0.0;
  Real abstol = 0.0;
  // Invoke with Lwork = -1 to estimate work space size
  dsygvx_(1,'V','I','U',n,myMat.ptr(),n,B.myMat.ptr(),n,
          dummy,dummy,1,nv,abstol,m,&val.front(),vec.ptr(),n,
          &dummy,-1,0,0,info);

  if (info == 0)
  {
    // Allocate work space
    int  Lwork = int(dummy);
    Real* work = new Real[Lwork];
    int* Iwork = new int[6*n];
    val.resize(n);
    vec.resize(n,nv);
    // Solve the eigenproblem
    dsygvx_(1,'V','I','U',n,myMat.ptr(),n,B.myMat.ptr(),n,
	    dummy,dummy,1,nv,abstol,m,&val.front(),vec.ptr(),n,
	    work,Lwork,Iwork+n,Iwork,info);
    delete[] work;
    delete[] Iwork;
    val.resize(nv);
    if (info == 0) return true;
  }

  std::cerr <<"LAPACK::DSYGVX: ";
  if (info < 0)
    std::cerr <<"Invalid argument #"<< -info << std::endl;
  else if ((size_t)info <= n)
    std::cerr << info <<" eigenvectors failed to converge."<< std::endl;
  else if ((size_t)info <= 2*n)
    std::cerr <<"The leading minor of order "<< info-n
	      <<" of matrix B is not positive definite."<< std::endl;
#else
  std::cerr <<"DSYGVX not available - built without LAPack/BLAS"<< std::endl;
#endif
  return false;
}


bool DenseMatrix::solveEigNon (RealArray& r_val, RealArray& c_val)
{
  const size_t n = myMat.rows();
  if (n < 1) return true; // No equations to solve

#ifdef USE_CBLAS
  int  info  = 0;
  Real dummy = 0.0;
  // Invoke with Lwork = -1 to estimate work space size
  dgeev_('N','N',n,myMat.ptr(),n,&r_val.front(),&c_val.front(),
	 &dummy,1,&dummy,1,&dummy,-1,info);

  if (info == 0)
  {
    // Allocate work space
    int  Lwork = int(dummy);
    Real* work = new Real[Lwork];
    r_val.resize(n);
    c_val.resize(n);
    // Solve the eigenproblem
    dgeev_('N','N',n,myMat.ptr(),n,&r_val.front(),&c_val.front(),
	   &dummy,1,&dummy,1,work,Lwork,info);
    delete[] work;
    if (info == 0) return true;
  }

#else
  std::cerr <<"DGEEV not available - built without LAPack/BLAS"<< std::endl;
#endif
  return false;
}

DenseMatrix operator*(Real alpha, const DenseMatrix& A)
{
  DenseMatrix B(A);
  B.getMat() *= alpha;
  return B;
}

DenseMatrix operator*(const DenseMatrix& A, Real alpha)
{
  DenseMatrix B(A);
  B.getMat() *= alpha;
  return B;
}
