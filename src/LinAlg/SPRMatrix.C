// $Id$
//==============================================================================
//!
//! \file SPRMatrix.C
//!
//! \date Jan 4 2008
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Representation of the system matrix on the SPR-format.
//!
//==============================================================================

#include "SPRMatrix.h"
#include "SAM.h"
#include <numeric>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#if defined(_WIN32) && !defined(__MINGW32__) && !defined(__MINGW64__)
#define sprprp_ SPRPRP
#define sprsas_ SPRSAS
#define sprrnm_ SPRRNM
#define sprtrs_ SPRTRS
#define sprsmb_ SPRSMB
#define sprpmp_ SPRPMP
#define spradm_ SPRADM
#define sprdad_ SPRDAD
#define sprprm_ SPRPRM
#define sprsol_ SPRSOL
#define sprlax_ SPRLAX
#elif defined(_AIX)
#define sprprp_ sprprp
#define sprsas_ sprsas
#define sprrnm_ sprrnm
#define sprtrs_ sprtrs
#define sprsmb_ sprsmb
#define sprpmp_ sprpmp
#define spradm_ spradm
#define sprdad_ sprdad
#define sprprm_ sprprm
#define sprsol_ sprsol
#define sprlax_ sprlax
#endif

extern "C" {
//! \brief Prepares the control information for the sparse assembly process.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprprp_(const int* madof, const int* minex,  const int* mpmnpc,
             const int* mmnpc, const int* mpmceq, const int* mmceq,
             const int* msc,   const int& nspar,  const int& lpu,
             const int* mpar,  int* mspar,
             const int* meqn,  int* iWork, int& ierr);
//! \brief Computes the SPR-version of the connectivity array \a mmnpc.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsas_(const int* mpar,  const int* mpmnpc, const int* mmnpc,
             const int* madof, const int* msc,    const int* mpmceq,
             const int* mmceq, const int* meqn,
             int* mspar, int* msica, int* iWork,
             const int& nspar, const int& lpu, int& ierr);
//! \brief Reorders the equations by means of the SPR-node partition.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprrnm_(int* mspar, int* msica, int* iWork, const int& lpu, int& ierr);
//! \brief Analyses the sparsity pattern of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprtrs_(int* mspar, int* msica, int* mtrees, const int* meqn,
             const int& nspar, const int& nmsica, const int& ntrees,
             const int& ndof, const int& niwork, int* iWork, Real* rinfo,
             const int& lpu, int& ierr);
//! \brief Performs the symbolic factorization of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsmb_(int* mspar, int* msica, int* mtrees, int* msifa,
             const int* meqn, int* iWork, const int& lpu, int& ierr);
//! \brief Finalizes the SPR control arrays \a msica, \a mtrees and \a msifa.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprpmp_(const int* msica, const int* mtrees, const int* meqn,
             const int* mpar, int* mspar, int* mvarnc);
//! \brief Assembles an element matrix \a EM into the system matrix \a SM.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void spradm_(const Real* eK, const Real* ttcc, const int* mpar,
             const int* mspar, const int* madof, const int* meqn,
             const int* mpmnpc, const int* mmnpc, const int* mpmceq,
             const int* mmceq, const int* msica, const int* mtrees,
             const int* msifa, const int* mvarnc, Real* values, Real* sysRHS,
             int* work, const int& iel, const int& nedof, const int& lpu,
             const int& nrhs, int& ierr);
//! \brief Adds a scalar value into the diagonal of the system matrix.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprdad_(const int* mpar, const int* mtrees, const int* msifa,
             Real* values, const Real& sigma, const int& lpu, int& ierr);
//! \brief Performs a matrix-matrix multiplication.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprprm_(const Real* A, const Real* B, Real* C, Real* rWork,
             const int* mspar, const int* mtrees, const int* msifa,
             const int& m, const int& n, const int& ksa, const int& iflag,
             const int& lpu, int& ierr);
//! \brief Solves the linear equation system \a Ax=b.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprsol_(const int& iop, const int* mspar, const int* mtrees,
             const int* msifa, Real* value, Real* B,
             const int& ldB, const int& nrhs, Real* tol,
             int* iWork, Real* rWork, const int& lpu, int& ierr);
//! \brief Solves the generalized eigenproblem \a Ax=&lambda;Bx.
//! \details This is a FORTRAN-77 subroutine in the SAM library.
//! \sa SAM library documentation.
void sprlax_(Real* A, Real* B, const Real* tol,
             const int* mparA, const int* mtreeA, const int* msifA,
             const int* mparB, const int* mtreeB, const int* msifB,
             const int& iop, Real* val, Real* vec, Real* rWork, int* iWork,
             const Real& shift, const int& n, const int& nv, const int& maxlan,
             const int& lpu, const int& ipsw, int& ierr);
}


SPRMatrix::SPRMatrix () : rWork(1)
{
  ierr = 0;
  memset(mpar,0,NS*sizeof(int));
  msica = msifa = mtrees = mvarnc = nullptr;
  values = nullptr;

#ifdef USE_OPENMP
  jWork = new std::vector<int>[omp_get_max_threads()];
#else
  jWork = &iWork;
#endif
}


SPRMatrix::SPRMatrix (const SPRMatrix& A) : SystemMatrix(A), rWork(1)
{
  ierr   = 0;
  msica  = new int[A.mpar[1]];
  msifa  = new int[A.mpar[2]];
  mtrees = new int[A.mpar[35]];
  mvarnc = new int[2*A.mpar[7]];
  values = new Real[A.mpar[7] + A.mpar[15]];

  memcpy(mpar  ,A.mpar  ,NS*sizeof(int));
  memcpy(msica ,A.msica ,A.mpar[1]*sizeof(int));
  memcpy(msifa ,A.msifa ,A.mpar[2]*sizeof(int));
  memcpy(mtrees,A.mtrees,A.mpar[35]*sizeof(int));
  memcpy(mvarnc,A.mvarnc,2*A.mpar[7]*sizeof(int));
  memcpy(values,A.values,(A.mpar[7]+A.mpar[15])*sizeof(Real));

#ifdef USE_OPENMP
  jWork = new std::vector<int>[omp_get_max_threads()];
#else
  jWork = &iWork;
#endif
}


SPRMatrix::~SPRMatrix ()
{
  delete[] msica;
  delete[] msifa;
  delete[] mtrees;
  delete[] mvarnc;
  delete[] values;
#ifdef USE_OPENMP
  delete[] jWork;
#endif
}


size_t SPRMatrix::dim (int idim) const
{
  switch (idim) {
  case 1:
  case 2:
    return mpar[7];
  case 3:
    return mpar[7] * mpar[7];
  default:
    return mpar[7] + mpar[15];
  }
}


/*!
  Must be called once before the element assembly loop.
  The SPR data structures are initialized and the all symbolic operations
  that are need before the actual assembly can start are performed.
*/

void SPRMatrix::initAssembly (const SAM& sam, bool)
{
  memset(mpar,0,NS*sizeof(int));
  msica = msifa = mtrees = mvarnc = nullptr;
  values = nullptr;

#ifdef HAS_SPR
  if (!sam.madof || !sam.mpmnpc)
  {
    std::cerr <<"SPRMatrix: SAM object is not properly initialized"<< std::endl;
    return;
  }

  iWork.resize(2*sam.nnod + 2*sam.ndof);
  int* minex = sam.minex;
  if (!minex)
  {
    // External node numbers are not provided - generate an identity array
    minex = iWork.data();
    std::iota(minex,minex+sam.nnod,1);
  }

  int* msc = sam.msc;
  for (int idof = 0; idof < sam.ndof; idof++)
    if (msc[idof] < 0)
    {
      // Make a temporary copy of msc with non-negative values only
      if (msc == sam.msc)
      {
        msc = new int[sam.ndof];
        memcpy(msc,sam.msc,sam.ndof*sizeof(int));
      }
      msc[idof] = 0;
    }

  sprprp_(sam.madof, minex, sam.mpmnpc, sam.mmnpc,
          sam.mpmceq, sam.mmceq, msc, NS, 6, sam.mpar,
          mpar, sam.meqn, iWork.data()+sam.nnod, ierr);
  if (msc != sam.msc) delete[] msc;
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRPRM: Failure "<< ierr << std::endl;
    return;
  }
  else if (sam.neq != mpar[7])
  {
    std::cerr <<"SAM::SPRPRM: Internal error, NEQ = "<< sam.neq
              <<" != "<< mpar[7] << std::endl;
    return;
  }
  else if (mpar[7] < 1)
    return; // No equations to solve

  // Perform symbolic assembly
  std::vector<int> itemp(mpar[34]);
  sprsas_(sam.mpar, sam.mpmnpc, sam.mmnpc, sam.madof, sam.msc, sam.mpmceq,
          sam.mmceq, sam.meqn, mpar, itemp.data(), iWork.data(), NS, 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRSAS: Failure "<< ierr << std::endl;
    return;
  }

  // Reallocation of msica with the correct size
  msica = new int[mpar[1]];
  memcpy(msica,itemp.data(),mpar[1]*sizeof(int));

  // Perform nodal reordering
  iWork.resize(mpar[36]);
  sprrnm_(mpar, msica, iWork.data(), 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRRNM: Failure "<< ierr << std::endl;
    return;
  }

  // Analyze the sparsity pattern
  mtrees = new int[mpar[35]];
  iWork.resize(mpar[37]);
  Real rinfo[4];
  sprtrs_(mpar, msica, mtrees, sam.meqn, NS, mpar[1], mpar[35], sam.ndof,
          mpar[37], iWork.data(), rinfo, 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRTRS: Failure "<< ierr << std::endl;
    return;
  }

  // Perform symbolic factorization
  itemp.resize(mpar[2]);
  iWork.resize(3*mpar[5] + 2*mpar[7] + 1);
  sprsmb_(mpar, msica, mtrees, itemp.data(), sam.meqn,
          iWork.data(), 6, ierr);
  if (ierr < 0)
  {
    std::cerr <<"SAM::SPRSMB: Failure "<< ierr << std::endl;
    return;
  }

  // Reallocation of msifa with the correct size
  msifa = new int[mpar[2]];
  memcpy(msifa,itemp.data(),mpar[2]*sizeof(int));

  // Finalize the SPR datastructure
  mvarnc = new int[2*mpar[7]];
  sprpmp_(msica, mtrees, sam.meqn, sam.mpar, mpar, mvarnc);

  // Allocate space for the matrix itself
  values = new Real[mpar[7] + mpar[15]];
  memset(values,0,(mpar[7] + mpar[15])*sizeof(Real));
#else
  std::cerr <<"SPRMatrix: SPR solver not available"<< std::endl;
#endif
}


void SPRMatrix::init ()
{
  if (mpar[0] < 4) return;

  mpar[0] = 4;
  memset(values,0,(mpar[7] + mpar[15])*sizeof(Real));
}


bool SPRMatrix::assemble (const Matrix& eM, const SAM& sam, int e)
{
  return this->assemble(e,eM,sam);
}


bool SPRMatrix::assemble (const Matrix& eM, const SAM& sam,
                          SystemVector& B, int e)
{
  if (B.getType() != LinAlg::DENSE) return false;

  return this->assemble(e,eM,sam,B.getPtr());
}


bool SPRMatrix::assemble (const Matrix&, const SAM&,
                          SystemVector&, const std::vector<int>&)
{
  std::cerr <<"SPRMatrix::assemble(const Matrix&,const SAM&,"
            <<"SystemVector&,const std::vector<int>&): Not implemented."
            << std::endl;
  return false;
}


bool SPRMatrix::assemble (int e, const Matrix& eM, const SAM& sam, Real* B)
{
#ifdef HAS_SPR
#ifdef USE_OPENMP
  std::vector<int>& IWORK = jWork[omp_in_parallel() ? omp_get_thread_num() : 0];
#else
  std::vector<int>& IWORK = *jWork;
#endif
  IWORK.resize(mpar[38]);
  spradm_(eM.ptr(), sam.ttcc, sam.mpar, mpar,
          sam.madof, sam.meqn, sam.mpmnpc, sam.mmnpc, sam.mpmceq, sam.mmceq,
          msica, mtrees, msifa, mvarnc, values, B ? B : rWork.data(),
          IWORK.data(), e, eM.rows(), 6, B ? 1 : 0, ierr);
  if (ierr == 0) return haveContributions = true;

  std::cerr <<"SAM::SPRADM: Failure "<< ierr << std::endl;
#endif
  return false;
}


void SPRMatrix::mult (Real alpha)
{
  cblas_dscal(mpar[7]+mpar[15], alpha, values, 1);
}


bool SPRMatrix::add (const SystemMatrix& B, Real alpha)
{
  const SPRMatrix* Bptr = dynamic_cast<const SPRMatrix*>(&B);
  if (!Bptr) return false;

  if (mpar[7] != Bptr->mpar[7] || mpar[15] != Bptr->mpar[15]) return false;

  if (B.isZero()) return true;

  cblas_daxpy(mpar[7]+mpar[15], alpha, Bptr->values, 1, values, 1);

  return haveContributions = true;
}


bool SPRMatrix::add (Real sigma)
{
#ifdef HAS_SPR
  sprdad_(mpar, mtrees, msifa, values, sigma, 6, ierr);
  if (ierr == 0) return haveContributions = true;

  std::cerr <<"SAM::SPRDAD: Failure "<< ierr << std::endl;
#endif
  return false;
}


bool SPRMatrix::multiply (size_t n, const Real* b, Real* c)
{
#ifdef HAS_SPR
  if (rWork.size() < n) rWork.resize(n);
  sprprm_(values+n, b, c, rWork.data(), mpar, mtrees, msifa,
          n, 1, 1, 1, 6, ierr);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRPRM: Failure "<< ierr << std::endl;
#endif
  return false;
}


bool SPRMatrix::multiply (const SystemVector& B, SystemVector& C) const
{
  if (mpar[7] <= 0) return true;

  if (B.getType() != LinAlg::DENSE) return false;
  if (C.getType() != LinAlg::DENSE) return false;

  C.redim(mpar[7]);
  return const_cast<SPRMatrix*>(this)->multiply(C.dim(),B.getRef(),C.getPtr());
}


bool SPRMatrix::solve (SystemVector& B, Real*)
{
  if (mpar[7] < 1) return true; // No equations to solve

  if (B.getType() != LinAlg::DENSE) return false;

#ifdef HAS_SPR
  Real tol[3] = { Real(1.0e-12), Real(0), Real(0) };
  iWork.resize(std::max(mpar[12],mpar[13]+1));
  rWork.resize(std::max(mpar[16],mpar[13]));
  int iop = mpar[0] < 5 ? 3 : 4;
  sprsol_(iop, mpar, mtrees, msifa, values, B.getPtr(),
          B.dim(), 1, tol, iWork.data(), rWork.data(), 6, ierr);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRSOL: Failure "<< ierr << std::endl;
#endif
  return false;
}


/*!
  The eigenproblem is assumed to be on the form
  \b A \b x = &lambda; \b B \b x where \b A ( = \a *this ) and \b B
  both are assumed to be symmetric and \b B also to be positive definite.
  The eigenproblem is solved by the SAM library subroutine \a SPRLAN.
  \sa SAM library documentation.
*/

bool SPRMatrix::solveEig (SPRMatrix& B, RealArray& val, Matrix& vec, int nev,
                          Real shift, int iop)
{
  const int n = mpar[7];
  if (n < 1 || nev == 0) return true; // No equations to solve
  if (n != B.mpar[7]) return false; // Incompatible matrices

#ifdef USE_F77SAM
  const Real tol[3] = { Real(1.0e-8), Real(1.0e-12), Real(1.0e-8) };
  int maxlan = 3*nev+12;
  if (maxlan > n) maxlan = n;
  int nrWork = std::max(maxlan*(maxlan+7), std::max(mpar[13],mpar[16]));
  val.resize(n);
  vec.resize(n,maxlan);
  iWork.resize(std::max(2*maxlan+mpar[12],n));
  rWork.resize(maxlan + 2*n + nrWork);
  sprlax_(values, B.values, tol, mpar, mtrees, msifa, B.mpar, B.mtrees, B.msifa,
          iop, val.data(), vec.ptr(), rWork.data(), iWork.data(),
          shift, n, nev, maxlan, 6, 0, ierr);
  val.resize(nev);
  vec.resize(n,nev);
  if (!ierr) return true;

  std::cerr <<"SAM::SPRLAX: Failure "<< ierr << std::endl;
#else
  std::cerr <<"SAM::SPRLAX eigensolver not available"<< std::endl;
#endif
  return false;
}


Real SPRMatrix::Linfnorm () const
{
  RealArray B(mpar[7],Real(1)), C(mpar[7],Real(0));
  if (const_cast<SPRMatrix*>(this)->multiply(B.size(),B.data(),C.data()))
    return *std::max_element(C.begin(),C.end());

  return Real(0);
}
